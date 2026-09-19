"""
pko_comparison.py
=================
Pullulan-Produzenten gegen PKO: verhalten sich die Kontrollen ohne Pullulan?

HINTERGRUND
-----------
Die Kontrollkammern (PosCtrl = durchgehend Feast, NegCtrl = durchgehend
Starvation) verhalten sich im Wildtyp nicht so, wie sie sollten. Die
Arbeitshypothese: Pullulan - das Exopolysaccharid, das der WT und die
Biosensor-Staemme ausscheiden - setzt die Chip-Strukturen zu und erzeugt
unerwartete Stroemungsprofile. Der PKO-Stamm produziert kein Pullulan.

DIE VERSUCHSSTRUKTUR, AUF DER DAS HIER STEHT
--------------------------------------------
Pro (Stamm, osc_type, Periode) gibt es EINEN Chip (experiment_units.py). Auf
ihm liegen 3 PosCtrl- und 3 NegCtrl-Kammern auf verschiedenen Arrays - nominell
identisch, gleiche Kultur, gleicher Tag, gleiches Medium. Ihre Uneinigkeit ist
deshalb der reinste verfuegbare Ausdruck dessen, was der CHIP mit einer Kammer
macht: der hydraulische Anteil. Genau den sollte eine Verstopfung aufblaehen.

PKO ist EIN Chip (Periode 3). Die Vergleichsgruppe:

  * PRIMAER (Entscheidung des Autors): die Chips der Produzenten-Staemme bei
    DERSELBEN Periode - WT/Glc/3 und die vier Biosensor-Staemme bei 3.
  * HINTERGRUND: alle Produzenten-Chips aller Perioden. Kontrollkammern sind
    periodenunabhaengig (PosCtrl ist durchgehend Feast, egal in welchem Batch
    sie mitlief), also ist jeder Chip ein gueltiger Vergleichspunkt. Die
    Verteilung ueber ~50 Chips zeigt, ob der PKO-Wert innerhalb oder
    ausserhalb dessen liegt, was Pullulan-Produzenten ueberhaupt tun.

DREI GROESSEN PRO CHIP, KEINE DAVON LINEAGE-ABHAENGIG
----------------------------------------------------
  1. Kammer-Uebereinstimmung: CV der Kammer-Medianflaeche ueber die 3 PosCtrl-
     (bzw. 3 NegCtrl-)Kammern EINES Chips. Das ist die Clogging-Signatur.
  2. Zeitliche Stabilitaet: Streuung der Residuen des Frame-Medians um einen
     linearen Zeittrend, relativ zum Niveau - pro Kontrollkammer, dann pro
     Chip gemittelt. Der glatte Trend (PosCtrl waechst) ist Biologie und wird
     abgezogen; erratische Spruenge und Drift bleiben. Konstantes Medium,
     also kein Aliasing.
  3. Kontroll-Bracket: trennt PosCtrl von NegCtrl ueberhaupt? Cliff's Delta
     auf mu_area (alle Zellen, alle Fits - siehe compute_control_bracket()).

Zell-Ausbeute ueber die Zeit gibt es NICHT mehr: Blastokonidien werden
zwischen Kammern gespuelt, die Zellzahl einer Kammer sagt deshalb nichts ueber
Verstopfung.

STATISTIK: WAS DIESES DESIGN TRAEGT
-----------------------------------
PKO ist EIN Chip. Es gibt keinen Test PKO-gegen-Produzenten; die Aussage ist,
wo der PKO-Chip in der Verteilung der Produzenten-Chips liegt. Innerhalb eines
Chips sind die 3 Kontrollkammern technische Wiederholungen - deshalb ist der
CV ueber sie die Kennzahl, nicht ein Test.

CONFOUND, DER IN DIE DISKUSSION GEHOERT
---------------------------------------
PKO ist eine Mutante mit veraendertem Stoffwechsel, nicht "WT ohne
Verstopfung". Unterscheidbar sind hydraulische und metabolische Erklaerung nur
ueber das MUSTER: hydraulisch heisst hohe Kammer-Streuung und zeitliche Drift
bei Produzenten, metabolisch heisst Niveau-Verschiebung bei erhaltener
Uebereinstimmung. Deshalb stehen hier Streuungen, nicht Mittelwerte.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from growth_rate import classify_condition_type

logger = logging.getLogger(__name__)

CONTROL_ORDER = ["NegCtrl", "PosCtrl"]
CONTROL_COLORS = {"NegCtrl": "#4C78A8", "PosCtrl": "#E45756"}
GROUP_PRODUCER = "producer"
GROUP_PKO = "PKO"
GROUP_COLORS = {GROUP_PRODUCER: "#6E6E6E", GROUP_PKO: "#54A24B"}


def _minutes_to_hours(minutes: float) -> float:
    """*_min-Angaben sind MINUTEN, 'time_h' ist in STUNDEN."""
    return minutes / 60.0


# ==============================================================================
# Vorbereitung
# ==============================================================================

def add_pullulan_group(df: pd.DataFrame, pko_biosensor: str) -> pd.DataFrame:
    """'pullulan' = 'PKO' fuer den PKO-Stamm, sonst 'producer' (WT und Biosensoren).

    Nicht 'WT' als Gruppenname: WT ist ein eigener Stamm in den Daten, und die
    Biosensor-Staemme sind ebenfalls Pullulan-Produzenten.
    """
    out = df.copy()
    out["pullulan"] = np.where(out["biosensor"].astype(str) == str(pko_biosensor), GROUP_PKO, GROUP_PRODUCER)
    return out


def add_condition_type(df: pd.DataFrame) -> pd.DataFrame:
    if "condition_type" in df.columns:
        return df
    out = df.copy()
    out["condition_type"] = out["condition"].map(classify_condition_type)
    return out


def matched_periods(cells_pko: pd.DataFrame) -> list[str]:
    periods = sorted(cells_pko["osc_freq"].dropna().astype(str).unique().tolist())
    logger.info("PKO deckt die Periode(n) %s ab - die Produzenten-Chips dieser Periode(n) sind der "
                "primaere Vergleich, alle uebrigen Produzenten-Chips der Hintergrund.", periods)
    return periods


def control_cells(cells: pd.DataFrame, pko_biosensor: str) -> pd.DataFrame:
    """Nur PosCtrl-/NegCtrl-Zeilen, mit chip, condition_type, pullulan."""
    if "chip" not in cells.columns:
        raise ValueError("control_cells(): 'chip' fehlt - experiment_units.add_experiment_units() zuerst.")
    df = add_condition_type(add_pullulan_group(cells, pko_biosensor))
    df = df[df["condition_type"].isin(CONTROL_ORDER)].copy()
    if df.empty:
        logger.warning(
            "control_cells(): keine PosCtrl-/NegCtrl-Zeilen. classify_condition_type() vergleicht "
            "exakt auf die Endungen '_PosCtrl'/'_NegCtrl'."
        )
    return df


# ==============================================================================
# 1 + 2. Kammer-Uebereinstimmung und zeitliche Stabilitaet
# ==============================================================================

def per_chamber_control_values(
    controls: pd.DataFrame,
    value_col: str = "area",
    analysis_start_min: float = 120.0,
) -> pd.DataFrame:
    """Pro Kontrollkammer: Median-Niveau und zeitlicher CV (beide nach der Vorkonditionierung).

    Stufe 1: Median ueber die Zellen eines Frames (robust gegen Segmentierungs-
    Ausreisser). Stufe 2 fuer das Niveau: Median ueber die Frames. Fuer die
    zeitliche Stabilitaet: SD/Mittelwert der Frame-Mediane ueber die Zeit.
    """
    if controls.empty:
        return pd.DataFrame()
    missing = [c for c in ("exp_id", "time_h", value_col) if c not in controls.columns]
    if missing:
        logger.warning("per_chamber_control_values(): Spalten %s fehlen.", missing)
        return pd.DataFrame()
    meta = [c for c in ["pullulan", "biosensor", "osc_type", "osc_freq", "chip", "condition",
                        "condition_type", "replicate", "chamber", "exp_id"] if c in controls.columns]
    per_time = (controls.dropna(subset=[value_col])
                .groupby(meta + ["time_h"], dropna=False)[value_col]
                .median().reset_index(name="frame_median"))
    late = per_time[per_time["time_h"] >= _minutes_to_hours(analysis_start_min)]
    if late.empty:
        logger.warning("per_chamber_control_values(): keine Frames ab %.0f min.", analysis_start_min)
        return pd.DataFrame()

    def _detrended_cv(grp: pd.DataFrame) -> float:
        # Zeitliche Stabilitaet = Streuung der RESIDUEN um einen linearen Trend,
        # relativ zum Niveau. Der glatte Trend (PosCtrl waechst durchgehend)
        # ist Biologie und gehoert nicht hinein; erratische Spruenge und Drift
        # - das, was eine Stroemungsstoerung macht - bleiben uebrig.
        y = grp["frame_median"].to_numpy(float); t = grp["time_h"].to_numpy(float)
        if len(y) < 4 or not np.isfinite(y).all() or abs(y.mean()) < 1e-9:
            return np.nan
        slope, intercept = np.polyfit(t, y, 1)
        resid = y - (slope * t + intercept)
        return float(resid.std(ddof=2) / abs(y.mean()))

    out = (late.groupby(meta, dropna=False)
           .apply(lambda g: pd.Series({"level": g["frame_median"].median(),
                                       "temporal_cv": _detrended_cv(g),
                                       "n_timepoints": len(g)}), include_groups=False)
           .reset_index())
    return out


def within_chip_agreement(per_chamber: pd.DataFrame) -> pd.DataFrame:
    """Pro Chip und Kontrollart: CV des Niveaus ueber die Kontrollkammern, mittlerer zeitlicher CV.

    Braucht mindestens zwei Kammern derselben Kontrollart auf dem Chip - sonst
    gibt es keine Kammer-zu-Kammer-Streuung, und die Zeile wird mit
    chamber_cv = NaN gefuehrt (nicht 0: 'nicht messbar' ist keine gute
    Uebereinstimmung).
    """
    if per_chamber.empty:
        return pd.DataFrame()
    key = [c for c in ["pullulan", "biosensor", "osc_type", "osc_freq", "chip", "condition_type"]
           if c in per_chamber.columns]
    out = (per_chamber.groupby(key, dropna=False)
           .agg(n_chambers=("level", "count"),
                level_mean=("level", "mean"),
                chamber_cv=("level", lambda s: s.std(ddof=1) / abs(s.mean()) if len(s) > 1 and s.mean() else np.nan),
                temporal_cv=("temporal_cv", "mean"))
           .reset_index())
    single = out["n_chambers"] < 2
    if single.any():
        logger.warning("within_chip_agreement(): %d Chip/Kontrollart-Paare haben nur eine Kammer - "
                       "chamber_cv dort NaN.", int(single.sum()))
    return out


# ==============================================================================
# 3. Kontroll-Bracket
# ==============================================================================

def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """Cliff's Delta P(a > b) - P(a < b) in [-1, 1], aus Mann-Whitney-U."""
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    if a.size == 0 or b.size == 0:
        return float("nan")
    u = stats.mannwhitneyu(a, b, alternative="two-sided").statistic
    return float(2.0 * u / (a.size * b.size) - 1.0)


def compute_control_bracket(area_table: pd.DataFrame, value_col: str = "mu_area") -> pd.DataFrame:
    """Pro Chip: trennt das PosCtrl/NegCtrl-Bracket? Cliff's Delta auf mu_area.

    Alle Zellen, ALLE Fits: fit_is_reliable heisst R^2 >= 0.5, und eine flache
    Bedingung (NegCtrl soll flach sein) erklaert per Konstruktion kaum Varianz -
    der Filter behielte dort nur die Tracks, in denen Rauschen wie ein Trend
    aussieht, und kippte das Vorzeichen (auf synthetischen Daten passiert:
    11 von 119 NegCtrl-Fits ueberlebten, mit hoeherem Median als PosCtrl).
    Der Anteil zuverlaessiger Fits wird stattdessen als Diagnose mitgefuehrt.
    cell_type wird nicht gefiltert - die Mutter/Bud-Heuristik ist unvalidiert.
    """
    if area_table is None or area_table.empty or "chip" not in area_table.columns:
        return pd.DataFrame()
    df = add_condition_type(area_table).dropna(subset=[value_col])
    df = df[df["condition_type"].isin(CONTROL_ORDER)]
    if df.empty:
        return pd.DataFrame()
    key = [c for c in ["pullulan", "biosensor", "osc_type", "osc_freq", "chip"] if c in df.columns]
    rows = []
    for keys, grp in df.groupby(key, dropna=False):
        pos = grp[grp["condition_type"] == "PosCtrl"]; neg = grp[grp["condition_type"] == "NegCtrl"]
        rec = dict(zip(key, keys if isinstance(keys, tuple) else (keys,)))
        rec["n_posctrl_cells"] = len(pos); rec["n_negctrl_cells"] = len(neg)
        rec["median_posctrl"] = float(pos[value_col].median()) if len(pos) else np.nan
        rec["median_negctrl"] = float(neg[value_col].median()) if len(neg) else np.nan
        rec["bracket_delta"] = cliffs_delta(pos[value_col].to_numpy(), neg[value_col].to_numpy())
        for lab, part in (("posctrl", pos), ("negctrl", neg)):
            rec[f"frac_reliable_{lab}"] = float(part["fit_is_reliable"].mean()) if "fit_is_reliable" in part.columns and len(part) else np.nan
        rows.append(rec)
    return pd.DataFrame(rows)


# ==============================================================================
# Abbildung
# ==============================================================================

def _strain_order(df: pd.DataFrame) -> list[str]:
    """Produzenten alphabetisch, WT zuletzt davon, PKO ganz rechts."""
    prod = sorted(df.loc[df["pullulan"] == GROUP_PRODUCER, "biosensor"].astype(str).unique())
    prod = [s for s in prod if s != "WT"] + (["WT"] if "WT" in prod else [])
    pko = sorted(df.loc[df["pullulan"] == GROUP_PKO, "biosensor"].astype(str).unique())
    return prod + pko


def _panel(ax, df, value_col, strain_order, matched_chips, ylabel, title, by_control=True, seed=0):
    """Ein Panel: pro Stamm die Chips als Punkte - passende Periode gefuellt, Rest blass."""
    rng = np.random.default_rng(seed)
    controls = CONTROL_ORDER if by_control and "condition_type" in df.columns else [None]
    for x, strain in enumerate(strain_order):
        sub_s = df[df["biosensor"].astype(str) == strain]
        for offset, ct in zip(np.linspace(-0.18, 0.18, len(controls)), controls):
            sub = sub_s if ct is None else sub_s[sub_s["condition_type"] == ct]
            sub = sub.dropna(subset=[value_col])
            if sub.empty:
                continue
            is_pko = (sub["pullulan"] == GROUP_PKO)
            matched = sub["chip"].isin(matched_chips) | is_pko
            color = CONTROL_COLORS[ct] if ct is not None else GROUP_COLORS[GROUP_PRODUCER]
            for sel, alpha, size, z in ((~matched, 0.28, 18, 2), (matched, 0.95, 46, 4)):
                part = sub[sel]
                if part.empty:
                    continue
                jitter = rng.uniform(-0.05, 0.05, size=len(part))
                col = GROUP_COLORS[GROUP_PKO] if (ct is None and part["pullulan"].eq(GROUP_PKO).all()) else color
                ax.scatter(np.full(len(part), x + offset) + jitter, part[value_col], s=size, alpha=alpha,
                           color=col, edgecolor="white", linewidth=0.5, zorder=z)
            bg = sub[~matched][value_col]
            if len(bg) >= 2:  # Spanne der Hintergrund-Chips als duenner Balken
                ax.plot([x + offset, x + offset], [bg.quantile(0.1), bg.quantile(0.9)],
                        color=color, alpha=0.35, linewidth=3, solid_capstyle="butt", zorder=1)
    ax.set_xticks(range(len(strain_order))); ax.set_xticklabels(strain_order, rotation=25, ha="right")
    ax.set_ylabel(ylabel, fontsize=9); ax.set_title(title, fontsize=10, loc="left")
    ax.grid(axis="y", alpha=0.25, linewidth=0.6)


def plot_pko_control_agreement(
    agreement: pd.DataFrame,
    bracket: pd.DataFrame,
    matched_chips: set,
    out_path: Path,
    periods: Optional[list[str]] = None,
) -> None:
    panels = []
    if agreement is not None and not agreement.empty:
        panels.append((agreement, "chamber_cv", "CV of chamber median area\n(3 control chambers, one chip)",
                       "a) Do nominally identical control chambers on one chip agree?", True))
        panels.append((agreement, "temporal_cv", "Detrended temporal CV\n(residuals around linear trend, chip mean)",
                       "b) Are the control chambers stable over the run (beyond their growth trend)?", True))
    if bracket is not None and not bracket.empty:
        panels.append((bracket, "bracket_delta", "PosCtrl vs NegCtrl\nCliff's δ on µ_area",
                       "c) Does the control bracket separate at all?", False))
    if not panels:
        logger.warning("plot_pko_control_agreement(): keine Daten - uebersprungen.")
        return
    strain_order = _strain_order(pd.concat([p[0] for p in panels], ignore_index=True))
    fig, axes = plt.subplots(len(panels), 1, figsize=(1.35 * len(strain_order) + 3.6, 3.3 * len(panels)),
                             squeeze=False, sharex=True)
    for ax, (df, col, ylab, title, by_ct) in zip(axes[:, 0], panels):
        _panel(ax, df, col, strain_order, matched_chips, ylab, title, by_control=by_ct)
        if col == "bracket_delta":
            ax.axhline(0.0, color="black", linewidth=0.9, linestyle="--", alpha=0.7)
    axes[-1, 0].set_xlabel("Strain (pullulan producers, WT last; PKO right)")
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker="o", linestyle="", color=CONTROL_COLORS["NegCtrl"], label="NegCtrl"),
               Line2D([], [], marker="o", linestyle="", color=CONTROL_COLORS["PosCtrl"], label="PosCtrl"),
               Line2D([], [], marker="o", linestyle="", color="#555555", markersize=8,
                      label=f"chip at the matched period ({', '.join(periods or [])} min) — primary comparison"),
               Line2D([], [], marker="o", linestyle="", color="#555555", alpha=0.3, markersize=5,
                      label="other chips of that strain (all periods) — backdrop; bar = 10–90 % range")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=8, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Pullulan producers vs PKO: control-chamber behaviour, one point per chip", y=0.995)
    fig.text(0.5, -0.09,
             "Each point is one chip (= one period batch, one culture). PKO is a single chip, so there is no "
             "test — read where it falls relative to the producers.\nControl chambers are constant-medium, so "
             "every producer chip is a valid comparator regardless of its period.",
             ha="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=180)
    plt.close(fig)
    logger.info("Plot gespeichert: %s", out_path.name)


# ==============================================================================
# Orchestrierung
# ==============================================================================

def run_pko_comparison(
    cells_producers: pd.DataFrame,
    cells_pko: pd.DataFrame,
    area_producers: pd.DataFrame,
    area_pko: pd.DataFrame,
    output_dir: Path,
    pko_biosensor: str = "PKO",
    value_col: str = "area",
    analysis_start_min: float = 120.0,
) -> None:
    """Der komplette Produzenten-gegen-PKO-Vergleich. Schreibt nur nach output_dir."""
    output_dir.mkdir(parents=True, exist_ok=True)
    # Ausgaben frueherer Versionen, die INHALTLICH ungueltig sind, entfernen:
    # die Zell-Ausbeute (Blastokonidien werden zwischen Kammern gespuelt) und
    # die 'per_replicate'-Tabellen (Rep ist kein Replikat). Eine Abbildung,
    # die als Evidenz verworfen wurde, darf nicht liegen bleiben.
    for stale in ("60_pko_cell_yield_slopes.csv", "60_pko_cell_yield_timeseries.csv", "62_pko_cell_yield.pdf",
                  "60_pko_control_bracket_per_replicate.csv", "60_pko_control_bracket_per_strain.csv",
                  "60_pko_chamber_agreement_per_replicate.csv", "60_pko_chamber_agreement_per_strain.csv"):
        path = output_dir / stale
        if path.exists():
            path.unlink()
            logger.info("Veraltete PKO-Ausgabe entfernt: %s", stale)
    if cells_pko.empty or cells_producers.empty:
        logger.warning("run_pko_comparison(): PKO- oder Produzenten-Daten fehlen - uebersprungen.")
        return
    periods = matched_periods(cells_pko)
    joint = pd.concat([control_cells(cells_producers, pko_biosensor),
                       control_cells(cells_pko, pko_biosensor)], ignore_index=True)
    if joint.empty:
        return
    matched_chips = set(joint.loc[joint["osc_freq"].astype(str).isin(periods), "chip"].unique())

    per_chamber = per_chamber_control_values(joint, value_col=value_col, analysis_start_min=analysis_start_min)
    if not per_chamber.empty:
        per_chamber.to_csv(output_dir / "60_pko_control_chambers.csv", index=False)
    agreement = within_chip_agreement(per_chamber)
    if not agreement.empty:
        agreement["matched_period"] = agreement["chip"].isin(matched_chips)
        agreement.to_csv(output_dir / "60_pko_within_chip_agreement.csv", index=False)
        show = agreement[agreement["matched_period"] | (agreement["pullulan"] == GROUP_PKO)]
        logger.info("Kammer-Uebereinstimmung auf den Chips der passenden Periode (CV ueber Kontrollkammern):\n%s",
                    show[["pullulan", "biosensor", "condition_type", "n_chambers", "chamber_cv", "temporal_cv"]]
                    .round(3).to_string(index=False))

    area_frames = [add_pullulan_group(t, pko_biosensor) for t in (area_producers, area_pko)
                   if t is not None and not t.empty and "chip" in t.columns]
    bracket = compute_control_bracket(pd.concat(area_frames, ignore_index=True)) if area_frames else pd.DataFrame()
    if not bracket.empty:
        bracket["matched_period"] = bracket["chip"].isin(matched_chips) | (bracket["pullulan"] == GROUP_PKO)
        bracket.to_csv(output_dir / "60_pko_control_bracket.csv", index=False)
        logger.info("Kontroll-Bracket (Cliff's delta, Chips der passenden Periode):\n%s",
                    bracket[bracket["matched_period"]][["pullulan", "biosensor", "bracket_delta",
                                                        "frac_reliable_posctrl", "frac_reliable_negctrl"]]
                    .round(3).to_string(index=False))

    plot_pko_control_agreement(agreement, bracket, matched_chips,
                               output_dir / "61_pko_control_agreement.pdf", periods=periods)
    logger.warning(
        "PKO-VERGLEICH: PKO ist EIN Chip; es gibt keinen Test. Die Aussage ist, wo der PKO-Chip in der "
        "Verteilung der Produzenten-Chips liegt (primaer: gleiche Periode; Hintergrund: alle). Unterscheidbar "
        "sind hydraulische und metabolische Erklaerung nur ueber das Muster (Streuung, Drift), nicht ueber "
        "einen Niveau-Unterschied. Siehe pko_comparison.py."
    )
