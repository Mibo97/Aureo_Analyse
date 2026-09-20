"""
bud_size.py
===========
Groessenkriterium der Mutter/Bud-Heuristik: wie gross ist eine neu
auftauchende Zelle beim ERSTEN Auftreten, relativ zu der Zelle, der sie
zugeordnet wuerde?

WARUM
-----
lineage.classify_mother_bud() macht jede Zelle, die nach dem ersten Frame
einer Kammer neu auftaucht, zum Knospen-Kandidaten und ordnet sie der
naechsten etablierten Zelle im adaptiven Radius zu. In den Kammern werden
aber laufend Blastokonidien angespuelt, die aus anderen Kammern stammen. Sie
tauchen ebenfalls "neu" auf, oft direkt neben einer sitzenden Zelle, und
bestehen die raeumliche Zuordnung wie eine Knospe. Das ist kein zufaelliger
Fehler: die Anspuelrate haengt von Stroemung und Zelldichte ab, also vom Chip
und von der Bedingung - validate_lineage.py (lv_03) zeigt fuer die
Glc-Oszillationen eine mit dem Batch konfundierte Erkennungsrate.

Eine echte Knospe ist beim ersten Auftreten deutlich KLEINER als ihre Mutter;
eine angespuelte Zelle ist in der Regel etwa so gross wie die Zelle, neben der
sie landet. Das Verhaeltnis bud_area / mother_area beim ersten Auftreten
trennt die beiden Faelle.

DIE SCHWELLE KOMMT AUS DEN DATEN
--------------------------------
Ueber alle Kandidaten ist das Verhaeltnis zweigipflig: kleine echte Knospen,
etwa gleich grosse angespuelte Zellen. Die Schwelle ist der Antimodus (das
Dichteminimum) zwischen den beiden Gipfeln, geschaetzt per Kerndichte auf
log10(Verhaeltnis) bei der KLEINSTEN Bandbreite, bei der die Dichte genau zwei
Gipfel hat (Silvermans kritische Bandbreite, Start bei Scotts Faktor). Drei
Gruende fuehren zum Rueckfall: zu wenige Kandidaten, keine erkennbare
Zweigipfligkeit (nur ein Gipfel, oder das Tal zwischen den Gipfeln ist zu
flach) oder ein Antimodus ausserhalb des plausiblen Bereichs
(config.BUD_SIZE_PLAUSIBLE_RANGE). Rueckfall heisst per Default: KEIN
Groessenfilter (config.BUD_MAX_AREA_FRACTION_FALLBACK = None, Schwelle inf).
Ein fester Rueckfallwert waere ein willkuerlicher Schnitt - auf den echten
Daten liegt die einzige Mode bei 0.4-0.5, ein Wert von 0.5 halbierte dort
die Events und machte die Erkennungsrate chip-abhaengiger, nicht weniger.
Welcher Fall eintrat, steht in 20_bud_size_threshold.csv (Spalte 'source').
In die Dichte gehen
nur Verhaeltnisse zwischen 0.05 und 1.5 ein (kde_window): groessere sind
sicher keine Knospen (eine grosse Zelle, die neben einer kleinen etablierten
Zelle landet) und wuerden nur einen dritten Gipfel rechts erzeugen, der die
Zwei-Gipfel-Logik stoert; kleinere sind Segmentierungsreste. Beide zaehlen
in n_outside_window und bleiben in der Abbildung sichtbar.

EINE Schwelle fuer alle Zweige (Oszillation, statisch, PKO, no_qc): das
Groessenverhaeltnis beim Abknospen ist eine Eigenschaft des Organismus und der
Optik, nicht der Bedingung. Die Tabelle fuehrt die Antimoden je Stamm/osc_type
zusaetzlich als Kontrolle auf - weicht einer stark ab, ist das ein Hinweis auf
eine andere Optik oder eine andere Population, nicht auf eine andere Biologie
des Knospens.

AUSGABEN (run_bud_size_threshold, in OUTPUT_DIR)
------------------------------------------------
    20_bud_size_at_appearance.csv   ein Kandidat pro Zeile: bud_area,
                                    mother_area, bud_area_fraction (+ Metadaten)
    20_bud_size_threshold.csv       Schwelle, Quelle, Gipfel, Antimodus -
                                    global (angewandt) und je Gruppe (Kontrolle)
    20_bud_size_at_appearance.pdf   Verteilung mit Dichte, Gipfeln und Schwelle
"""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from lineage import LineageParams, FluxChannelConfig, classify_mother_bud

logger = logging.getLogger(__name__)

CANDIDATE_COLS = ["exp_id", "mother_cell_uid", "bud_cell_uid", "budding_frame",
                  "mother_area", "bud_area", "bud_area_fraction",
                  # Diagnose (lineage.classify_mother_bud): Flaechenbilanz der
                  # Mutter, Wachstum des Kandidaten, Kontakt, Tracklaengen.
                  "mother_area_prev", "mother_area_next", "mother_area_drop",
                  "mother_area_drop_over_bud", "mother_age_frames",
                  "bud_area_plus1", "bud_area_plus3", "contact_ratio",
                  "distance_px", "adaptive_radius_px", "bud_final_track_length",
                  "bud_was_washed_out", "mother_eccentricity", "bud_eccentricity"]
META_COLS = ["biosensor", "osc_type", "osc_freq", "condition", "chip", "chip_family", "medium"]
GROUP_COLS = ["biosensor", "osc_type"]

# Auswertungsgitter der Kerndichte in log10(Verhaeltnis): 0.02 bis 5.
_GRID = np.linspace(np.log10(0.02), np.log10(5.0), 600)
_TICKS = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]
# Abbildungstexte englisch wie in den uebrigen Abbildungen der Pipeline;
# die 'source'-Codes ebenfalls, weil sie im Titel der Abbildung stehen.
XLABEL = "area at first appearance / area of the assigned mother"
YLABEL = "density (log10 scale)"


@dataclass
class BudSizeThreshold:
    """Ergebnis von resolve_bud_size_threshold() - eine Zeile in 20_bud_size_threshold.csv."""
    threshold: float
    source: str
    n_candidates: int
    antimode: float = np.nan
    peak_low: float = np.nan
    peak_high: float = np.nan
    valley_depth: float = np.nan
    bw_factor: float = np.nan
    fallback: float = np.nan
    plausible_low: float = np.nan
    plausible_high: float = np.nan
    window_low: float = np.nan
    window_high: float = np.nan
    n_outside_window: int = 0

    def as_row(self) -> dict:
        return asdict(self)


def _log_fractions(fractions) -> np.ndarray:
    x = pd.Series(fractions, dtype=float).to_numpy()
    x = x[np.isfinite(x) & (x > 0)]
    return np.log10(x)


def _density(logx: np.ndarray, bw_factor: float) -> np.ndarray:
    return stats.gaussian_kde(logx, bw_method=bw_factor)(_GRID)


def _modes(dens: np.ndarray, min_rel_height: float = 0.05) -> list[int]:
    """Indizes lokaler Maxima; kleine Randhoecker (< 5 % des Hauptgipfels) zaehlen nicht."""
    inner = np.arange(1, len(dens) - 1)
    is_peak = (dens[inner] > dens[inner - 1]) & (dens[inner] >= dens[inner + 1])
    return [int(i) for i in inner[is_peak] if dens[i] >= min_rel_height * dens.max()]


def resolve_bud_size_threshold(
    fractions,
    fallback: Optional[float] = None,
    min_candidates: int = 50,
    plausible: Sequence[float] = (0.15, 0.9),
    min_valley_depth: float = 0.2,
    max_bandwidth_steps: int = 40,
    kde_window: Sequence[float] = (0.05, 1.5),
) -> BudSizeThreshold:
    """Antimodus der zweigipfligen Verteilung von bud_area / mother_area.

    fractions : Verhaeltnisse (ein Wert pro Kandidat); NaN, 0 und negative
                Werte werden ignoriert.
    fallback  : Schwelle, wenn kein Antimodus bestimmbar ist. None (Default)
                = unendlich = KEIN Groessenfilter.
    min_candidates : darunter keine Schaetzung, sondern fallback.
    plausible : Bereich, in dem der Antimodus liegen muss (sonst fallback).
    min_valley_depth : das Tal muss mindestens diesen Anteil unter dem
                niedrigeren der beiden Gipfel liegen (0.2 = 20 %); sonst gilt
                die Verteilung als nicht zweigipflig.
    kde_window : nur Verhaeltnisse in diesem Bereich gehen in die Dichte ein
                (siehe Modul-Docstring); die uebrigen zaehlen in
                n_outside_window.
    """
    fallback = float("inf") if fallback is None else float(fallback)
    logx_all = _log_fractions(fractions)
    lo_w, hi_w = np.log10(kde_window[0]), np.log10(kde_window[1])
    logx = logx_all[(logx_all >= lo_w) & (logx_all <= hi_w)]
    n = int(len(logx))
    base = dict(n_candidates=int(len(logx_all)), n_outside_window=int(len(logx_all) - n),
                fallback=fallback,
                plausible_low=float(plausible[0]), plausible_high=float(plausible[1]),
                window_low=float(kde_window[0]), window_high=float(kde_window[1]))
    if n < min_candidates:
        return BudSizeThreshold(float(fallback),
                                f"fallback: only {n} candidates inside the window (< {min_candidates})", **base)
    if np.ptp(logx) < 1e-9:
        return BudSizeThreshold(float(fallback), "fallback: all ratios identical", **base)

    # Kritische Bandbreite: bei Scotts Faktor starten und so lange glaetten,
    # bis hoechstens zwei Gipfel uebrig sind.
    factor = float(stats.gaussian_kde(logx).factor)
    dens = _density(logx, factor)
    peaks = _modes(dens)
    steps = 0
    while len(peaks) > 2 and steps < max_bandwidth_steps:
        factor *= 1.1
        dens = _density(logx, factor)
        peaks = _modes(dens)
        steps += 1
    if len(peaks) != 2:
        note = "one mode" if len(peaks) < 2 else "more than two modes even after smoothing"
        peak = float(10 ** _GRID[peaks[0]]) if peaks else np.nan
        return BudSizeThreshold(float(fallback), f"fallback: distribution not bimodal ({note})",
                                peak_low=peak, bw_factor=factor, **base)

    lo, hi = peaks
    valley = lo + int(np.argmin(dens[lo:hi + 1]))
    depth = float(1.0 - dens[valley] / min(dens[lo], dens[hi]))
    antimode = float(10 ** _GRID[valley])
    detail = dict(antimode=antimode, peak_low=float(10 ** _GRID[lo]), peak_high=float(10 ** _GRID[hi]),
                  valley_depth=depth, bw_factor=factor)
    if depth < min_valley_depth:
        return BudSizeThreshold(
            float(fallback),
            f"fallback: valley between the modes too shallow ({depth:.2f} < {min_valley_depth})",
            **detail, **base,
        )
    if not (plausible[0] <= antimode <= plausible[1]):
        return BudSizeThreshold(
            float(fallback),
            f"fallback: antimode {antimode:.2f} outside [{plausible[0]}, {plausible[1]}]",
            **detail, **base,
        )
    return BudSizeThreshold(antimode, "antimode", **detail, **base)


def bud_size_candidates(
    cells: pd.DataFrame,
    params: Optional[LineageParams] = None,
    flux_config: Optional[FluxChannelConfig] = None,
) -> pd.DataFrame:
    """Alle Kandidaten, die die raeumliche Zuordnung bestehen, OHNE
    Groessenfilter - genau die Population, aus der die Schwelle kommt."""
    events = classify_mother_bud(cells, params, flux_config=flux_config, bud_size_threshold=np.inf)
    if events.empty:
        return pd.DataFrame(columns=CANDIDATE_COLS)
    out = events[[c for c in CANDIDATE_COLS if c in events.columns]].copy()
    meta = [c for c in META_COLS if c in cells.columns]
    if meta:
        out = out.merge(cells[["exp_id"] + meta].drop_duplicates("exp_id"), on="exp_id", how="left")
    return out


def bud_size_threshold_table(
    candidates: pd.DataFrame,
    global_result: BudSizeThreshold,
    group_cols: Sequence[str] = GROUP_COLS,
    **resolve_kwargs,
) -> pd.DataFrame:
    """Eine Zeile 'global' (die angewandte Schwelle) plus eine Kontrollzeile je
    Gruppe: haette dieselbe Regel pro Stamm/osc_type dasselbe ergeben?"""
    rows = [{"scope": "global (applied)", "group": "all", **global_result.as_row()}]
    cols = [c for c in group_cols if c in candidates.columns]
    if cols and not candidates.empty:
        for keys, grp in candidates.groupby(cols, dropna=False):
            keys = keys if isinstance(keys, tuple) else (keys,)
            res = resolve_bud_size_threshold(grp["bud_area_fraction"], **resolve_kwargs)
            rows.append({"scope": "/".join(cols) + " (check only)",
                         "group": "/".join(map(str, keys)), **res.as_row()})
    return pd.DataFrame(rows)


def _windowed(logx: np.ndarray, result: BudSizeThreshold) -> np.ndarray:
    if not (np.isfinite(result.window_low) and np.isfinite(result.window_high)):
        return logx
    return logx[(logx >= np.log10(result.window_low)) & (logx <= np.log10(result.window_high))]


def _format_axis(ax: plt.Axes, result: BudSizeThreshold) -> None:
    ax.set_xlim(_GRID[0], _GRID[-1])
    for edge in (result.window_low, result.window_high):
        if np.isfinite(edge):
            ax.axvline(np.log10(edge), color="0.6", ls=":", lw=1)
    ax.set_xticks(np.log10(_TICKS))
    ax.set_xticklabels([f"{t:g}" for t in _TICKS])
    ax.set_xlabel(XLABEL)
    ax.set_ylabel(YLABEL)
    if np.isfinite(result.plausible_low) and np.isfinite(result.plausible_high):
        ax.axvspan(np.log10(result.plausible_low), np.log10(result.plausible_high),
                   color="C3", alpha=0.05, lw=0)


def plot_bud_size_distribution(
    candidates: pd.DataFrame,
    result: BudSizeThreshold,
    out_path: Path,
    group_cols: Sequence[str] = GROUP_COLS,
) -> None:
    """Links: gepoolte Verteilung mit Kerndichte, Gipfeln und Schwelle.
    Rechts: Kerndichte je Gruppe, gleiche Schwelle - sitzt sie ueberall im Tal?"""
    logx = _log_fractions(candidates["bud_area_fraction"]) if not candidates.empty else np.array([])
    has_threshold = np.isfinite(result.threshold) and result.threshold > 0
    thr_log = float(np.log10(result.threshold)) if has_threshold else np.nan
    thr_label = (f"threshold {result.threshold:.2f} (antimode)" if result.source == "antimode"
                 else f"threshold {result.threshold:.2f} (fallback)")

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.3))
    if logx.size:
        bins = np.linspace(_GRID[0], _GRID[-1], 71)
        ax.hist(logx, bins=bins, density=True, color="0.8", edgecolor="white",
                label=f"candidates (n = {logx.size})")
        logx_in = _windowed(logx, result)
        if np.isfinite(result.bw_factor) and logx_in.size >= 2 and np.ptp(logx_in) > 1e-9:
            # Dichte auf dem Fenster, auf den Anteil der Kandidaten im Fenster
            # skaliert, damit sie zum Histogramm ueber ALLE Kandidaten passt.
            dens = _density(logx_in, result.bw_factor) * (logx_in.size / logx.size)
            ax.plot(_GRID, dens, color="C0", lw=2,
                    label="kernel density (critical bandwidth, dotted = window)")
            for peak, name in ((result.peak_low, "buds"), (result.peak_high, "washed in")):
                if np.isfinite(peak):
                    i = int(np.argmin(np.abs(_GRID - np.log10(peak))))
                    ax.annotate(name, (_GRID[i], dens[i]), xytext=(0, 6), textcoords="offset points",
                                ha="center", fontsize=8, color="C0")
    if has_threshold:
        ax.axvline(thr_log, color="C3", ls="--", lw=1.5, label=thr_label)
    else:
        ax.plot([], [], " ", label="no size filter applied (not bimodal)")
    ax.set_title("all branches pooled", fontsize=10)
    _format_axis(ax, result)
    ax.legend(fontsize=8, loc="upper left")

    cols = [c for c in group_cols if c in candidates.columns]
    n_lines = 0
    if cols and not candidates.empty:
        for keys, grp in candidates.groupby(cols, dropna=False):
            lx = _windowed(_log_fractions(grp["bud_area_fraction"]), result)
            if lx.size < 20 or np.ptp(lx) < 1e-9:
                continue
            factor = result.bw_factor if np.isfinite(result.bw_factor) else None
            d = stats.gaussian_kde(lx, bw_method=factor)(_GRID)
            keys = keys if isinstance(keys, tuple) else (keys,)
            ax2.plot(_GRID, d, lw=1.5, label="/".join(map(str, keys)) + f" (n = {lx.size})")
            n_lines += 1
    if has_threshold:
        ax2.axvline(thr_log, color="C3", ls="--", lw=1.5)
    ax2.set_title("per " + "/".join(cols) if cols else "per group", fontsize=10)
    _format_axis(ax2, result)
    if n_lines:
        ax2.legend(fontsize=8, loc="upper left")
    else:
        ax2.text(0.5, 0.5, "no group with >= 20 candidates", ha="center", va="center",
                 transform=ax2.transAxes, fontsize=9, color="0.4")

    fig.suptitle(
        "Size criterion of the bud heuristic - threshold source: " + result.source,
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=180)
    plt.close(fig)


def run_bud_size_threshold(
    cells: pd.DataFrame,
    out_dir: Path,
    params: Optional[LineageParams] = None,
    fallback: Optional[float] = None,
    plausible: Sequence[float] = (0.15, 0.9),
    min_candidates: int = 50,
) -> float:
    """Schwelle ableiten, Tabellen und Abbildung schreiben, Schwelle zurueckgeben
    (inf = kein Groessenfilter)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    kwargs = dict(fallback=fallback, plausible=plausible, min_candidates=min_candidates)

    candidates = bud_size_candidates(cells, params)
    fractions = candidates["bud_area_fraction"] if not candidates.empty else []
    result = resolve_bud_size_threshold(fractions, **kwargs)
    table = bud_size_threshold_table(candidates, result, **kwargs)

    candidates.to_csv(out_dir / "20_bud_size_at_appearance.csv", index=False)
    table.to_csv(out_dir / "20_bud_size_threshold.csv", index=False)
    plot_bud_size_distribution(candidates, result, out_dir / "20_bud_size_at_appearance.pdf")

    logger.info(
        "Groessenkriterium: Schwelle %.3f (%s) aus %d Kandidaten; Gipfel %.2f / %.2f, "
        "Tal-Tiefe %.2f. Tabellen: 20_bud_size_threshold.csv, 20_bud_size_at_appearance.csv",
        result.threshold, result.source, result.n_candidates,
        result.peak_low, result.peak_high, result.valley_depth,
    )
    if not np.isfinite(result.threshold):
        logger.warning(
            "Groessenkriterium: KEIN Groessenfilter aktiv - %s. Die Verteilung von "
            "bud_area / mother_area traegt keine Schwelle; die Events bleiben ungefiltert. "
            "Diagnose-Spalten in 20_bud_size_at_appearance.csv.", result.source,
        )
    elif result.source != "antimode":
        logger.warning(
            "Groessenkriterium: fester Rueckfallwert %.2f in Kraft (%s) - das ist ein "
            "willkuerlicher Schnitt, vor der Interpretation 20_bud_size_at_appearance.pdf "
            "ansehen.", result.threshold, result.source,
        )
    return float(result.threshold)


def load_bud_size_threshold(output_dir: Path) -> Optional[float]:
    """Die von run_analysis.py angewandte Schwelle aus 20_bud_size_threshold.csv,
    None wenn die Datei fehlt."""
    path = Path(output_dir) / "20_bud_size_threshold.csv"
    if not path.exists():
        logger.info("Keine %s gefunden - keine abgeleitete Groessenschwelle verfuegbar.", path)
        return None
    table = pd.read_csv(path)
    if "threshold" not in table.columns:
        return None
    rows = table[table["scope"].astype(str).str.startswith("global")] if "scope" in table.columns else table
    if rows.empty:
        return None
    threshold = float(rows["threshold"].iloc[0])
    source = rows["source"].iloc[0] if "source" in rows.columns else "?"
    logger.info("Groessenkriterium aus %s: Schwelle %.3f (%s)", path.name, threshold, source)
    return threshold
