"""
endpoint_trends.py
==================
Kumulativer Endzustand gegen die Zyklusperiode - plus der Trendtest, der dazu
gehoert.

WARUM DIESE DATEI EXISTIERT
---------------------------
Bei MIN_PER_FRAME = 10 liegt die kuerzeste aufloesbare Periode bei 20 min;
fuenf der sechs Bedingungen liegen darunter (siehe config.py und den
README-Abschnitt zur Abtastung). Ein einzelner Zyklus ist damit NICHT
beobachtbar, und eine scheinbare Periodizitaet in einer Zeitreihe waere ein
Alias-Artefakt. Interpretierbar ist nur die KUMULATIVE Wirkung ueber Stunden.

Genau diese Auswertung fehlte: die Pipeline erzeugte Zeitreihen (10_/30_/31_)
und Varianzmasse (40_), aber keine einzige Abbildung, die den ENDZUSTAND einer
Bedingung gegen die Periode stellt. Fuer die Sensor-Ratios gab es sie
ueberhaupt nicht - 50_summary_statistics.csv bekommt nur intensity_cols
uebergeben, die ratio_*-Spalten erscheinen dort also nie.

Ebenso fehlte der Test, der zur eigentlichen Vorhersage passt. Die Erwartung
ist MONOTON in der Periode (lange Famine-Halbzyklen belasten mehr, siehe
README), und Monotonie prueft man mit einer Rangkorrelation. Kruskal-Wallis -
der einzige Test, der bisher verdrahtet war, und das nur fuer die Kontrollen -
prueft "irgendeine Gruppe unterscheidet sich" und ist dafuer das falsche
Werkzeug.

WIE DER ENDZUSTAND GEBILDET WIRD
--------------------------------
Dreistufig, damit nichts pseudorepliziert wird:

    1. pro Kammer: Mittelwert ueber die Zellen in den letzten
       ENDPOINT_LAST_FRACTION der Frames DIESER Kammer
    2. pro Replikat: Mittelwert ueber die Kammern des Replikats
       (Kammern sind technische Messungen desselben Chips)
    3. ueber Replikate: Mittelwert + SEM, n = echte Replikatzahl

Das Fenster ist relativ zur jeweiligen Kammer und nicht absolut, weil Kammern
unterschiedlich lang aufgenommen sein koennen - ein fester Frame-Bereich
wuerde sonst bei kurzen Kammern ins Leere greifen.

SPEARMAN, NICHT PEARSON
-----------------------
Geprueft wird Monotonie, nicht Linearitaet: die Perioden sind geometrisch
gestuft (0.75 ... 24, jeweils Faktor 2), und es gibt keinen Grund, einen
linearen Zusammenhang in der Periode zu erwarten. Getestet wird auf
REPLIKAT-Ebene (ein Wert pro Replikat und Periode), nicht ueber Zellen -
sonst haengt das p fast nur an der Zellzahl.

Kontrollen (PosCtrl/NegCtrl) gehen NICHT in die Korrelation ein: sie haben
keine Periode, der Ordnername ihres Batches ist keine Behandlung. Sie
erscheinen in der Abbildung als Referenzbaender, was der eigentliche Zweck
der Kontrollen ist.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from scipy import stats

from growth_rate import classify_condition_type
from experiment_units import summarise_hierarchical

logger = logging.getLogger(__name__)

from plot_style import (CONTROL_LINESTYLES, INK, INK_MUTED, INK_SOFT, SURFACE, control_handles, edge_color,
                        errorbar_kwargs, finish, legend_below, marker_kwargs, ordered_strains, panel_title,
                        strain_color, strain_handles)
from matplotlib.lines import Line2D

CONTROL_TYPES = ["NegCtrl", "PosCtrl"]

# Spearman braucht mindestens so viele verschiedene Perioden, damit eine
# Rangkorrelation ueberhaupt etwas aussagt. Bei zwei Punkten ist rho immer
# +-1 - das ist keine Evidenz, sondern Arithmetik.
MIN_PERIODS_FOR_TREND = 3


def add_condition_type(df: pd.DataFrame) -> pd.DataFrame:
    """Ergaenzt 'condition_type', falls noch nicht vorhanden."""
    if "condition_type" in df.columns:
        return df
    if "condition" not in df.columns:
        raise ValueError("add_condition_type() braucht eine 'condition'-Spalte.")
    out = df.copy()
    out["condition_type"] = out["condition"].map(classify_condition_type)
    return out


def period_minutes(values: pd.Series) -> pd.Series:
    """'osc_freq' -> Periode in Minuten als float, sonst NaN.

    Die Spalte heisst 'osc_freq', enthaelt aber die PERIODE in Minuten
    (config.OSC_FREQ_IS_PERIOD_IN_MINUTES). Statische Bedingungen
    ('static_ypd') und alles andere Nicht-Numerische werden zu NaN und fallen
    damit aus jeder Trendrechnung, statt eine Scheinzahl zu erzeugen.
    """
    return pd.to_numeric(values, errors="coerce")


def summarise_per_replicate(
    df: pd.DataFrame,
    value_col: str,
    group_cols: Optional[Sequence[str]] = None,
    replicate_col: str = "chip",
    chamber_col: str = "exp_id",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Zelle -> Kammer -> Chip -> Bedingung; siehe experiment_units.summarise_hierarchical().

    Der Name bleibt aus Kompatibilitaet, die Einheit ist jetzt der CHIP: bei
    den Oszillationsdaten gibt es davon einen pro Bedingung, der Fehlerbalken
    ist dann der Kammer-Fehler dieses Chips und in 'error_unit' auch so
    beschriftet. Rueckgabe: (pro Chip, pro Bedingung) - die Chip-Ebene ist die
    Grundlage des Spearman-Tests (ein Wert je Chip und Periode), die
    Bedingungs-Ebene die des Plots.
    """
    _, per_chip, per_condition = summarise_hierarchical(
        df, value_col, condition_cols=group_cols, chamber_col=chamber_col, chip_col=replicate_col,
    )
    return per_chip, per_condition


def detect_saturation_frame(
    cells: pd.DataFrame,
    level: float = 0.90,
    smooth_frames: int = 5,
    min_growth: float = 1.5,
    chamber_col: str = "exp_id",
    cell_id_col: str = "cell_uid",
) -> pd.DataFrame:
    """Pro Kammer: ab welchem Frame ist die Zellzahl gesaettigt?

    Saettigung = erster Frame, ab dem die (rollend geglaettete) Zahl
    verfolgbarer Zellen >= level x ihres Maximums liegt - aber NUR fuer
    Kammern, die ueberhaupt gewachsen sind (Maximum / Startwert >= min_growth).
    Eine flache Kammer liegt von Anfang an bei 90 % ihres Maximums und wuerde
    sonst als 'sofort gesaettigt' das Fenster auf die ersten Stunden ziehen;
    sie bekommt stattdessen ihren letzten Frame und begrenzt nichts. Eine
    Kammer, die bis zum Ende waechst, saettigt erst nahe ihrem letzten Frame -
    gewollt: begrenzend ist die Kammer, die FRUEH voll ist (ypd).

    Rueckgabe: exp_id, n_frames, n_cells_max, saturation_frame, saturated_early
    (True, wenn die Saettigung vor 75 % der Aufnahme liegt).
    """
    if cells is None or cells.empty or chamber_col not in cells.columns:
        return pd.DataFrame()
    counts = cells.groupby([chamber_col, "frame"])[cell_id_col].nunique().reset_index(name="n")
    records = []
    for exp_id, grp in counts.groupby(chamber_col):
        grp = grp.sort_values("frame")
        smooth = grp["n"].rolling(smooth_frames, center=True, min_periods=1).median()
        peak = float(smooth.max())
        start = float(smooth.iloc[0]) if smooth.iloc[0] > 0 else float("nan")
        growth = peak / start if start and start == start else float("inf")
        if growth >= min_growth:
            hit = grp.loc[smooth >= level * peak, "frame"]
            sat = int(hit.iloc[0]) if not hit.empty else int(grp["frame"].max())
        else:
            sat = int(grp["frame"].max())  # nicht gewachsen -> begrenzt das Fenster nicht
        n_frames = int(grp["frame"].nunique())
        records.append({chamber_col: exp_id, "n_frames": n_frames, "n_cells_max": int(grp["n"].max()),
                        "growth_factor": round(growth, 2), "saturation_frame": sat,
                        "saturated_early": growth >= min_growth and sat < 0.75 * (grp["frame"].max() + 1)})
    out = pd.DataFrame(records)
    meta_cols = [c for c in ["biosensor", "osc_type", "osc_freq", "condition", "medium", "chip_family",
                             "chip", "replicate"] if c in cells.columns]
    if meta_cols:
        out = out.merge(cells[[chamber_col] + meta_cols].drop_duplicates(chamber_col), on=chamber_col, how="left")
    n_early = int(out["saturated_early"].sum())
    logger.info(
        "Saettigung bestimmt fuer %d Kammern; %d davon saettigen frueh (< 75 %% der Aufnahme). "
        "Frueheste Saettigung: Frame %d.", len(out), n_early, int(out["saturation_frame"].min()),
    )
    return out


def static_endpoint_window(
    saturation: pd.DataFrame,
    last_fraction: float = 0.25,
    min_frame: int = 0,
) -> tuple[int, int]:
    """Absolutes Endfenster fuer den statischen Zweig.

    Endet an der FRUEHESTEN Saettigung ueber alle Kammern (damit keine
    ueberwachsene Kammer im Fenster liegt) und ist last_fraction dieser
    Spanne breit; beginnt nie vor min_frame (Ende der Vorkonditionierung).
    """
    hi = int(saturation["saturation_frame"].min())
    lo = max(min_frame, int(round(hi - last_fraction * hi)))
    if hi - lo < 3:
        logger.warning(
            "static_endpoint_window(): das Fenster %d-%d ist sehr schmal - die frueheste "
            "Saettigung liegt nahe der Vorkonditionierung. STATIC_SATURATION_LEVEL pruefen.", lo, hi,
        )
    return lo, hi


def compute_endpoint_per_replicate(
    cells: pd.DataFrame,
    value_col: str,
    last_fraction: float = 0.25,
    group_cols: Optional[Sequence[str]] = None,
    frame_window: Optional[tuple[int, int]] = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Endzustand einer Bedingung.

    Ohne frame_window: die letzten `last_fraction` der Frames PRO KAMMER
    (relativ zu deren eigenem Frame-Bereich - Kammern koennen unterschiedlich
    lang aufgenommen sein). Mit frame_window=(lo, hi): ein ABSOLUTES Fenster,
    dasselbe fuer alle Kammern - fuer den statischen Zweig, wo ypd ueberwaechst
    und das relative Fenster Unvergleichbares vergleichen wuerde.

    Anschliessend laeuft dieselbe hierarchische Aggregation wie in
    summarise_per_replicate().
    """
    if cells is None or cells.empty or value_col not in cells.columns:
        return pd.DataFrame(), pd.DataFrame()
    if not 0.0 < last_fraction <= 1.0:
        raise ValueError(f"last_fraction muss in (0, 1] liegen, ist {last_fraction}.")
    if "frame" not in cells.columns:
        logger.warning("compute_endpoint_per_replicate(): Spalte 'frame' fehlt - uebersprungen.")
        return pd.DataFrame(), pd.DataFrame()

    work = cells.dropna(subset=[value_col]).copy()
    if work.empty:
        logger.info(
            "compute_endpoint_per_replicate(): keine gueltigen '%s'-Werte - uebersprungen.",
            value_col,
        )
        return pd.DataFrame(), pd.DataFrame()

    chamber_col = "exp_id" if "exp_id" in work.columns else None
    if chamber_col is None:
        logger.warning("compute_endpoint_per_replicate(): 'exp_id' fehlt - uebersprungen.")
        return pd.DataFrame(), pd.DataFrame()

    if frame_window is not None:
        lo, hi = frame_window
        endpoint_rows = work[(work["frame"] >= lo) & (work["frame"] <= hi)]
    else:
        fmin = work.groupby(chamber_col)["frame"].transform("min")
        fmax = work.groupby(chamber_col)["frame"].transform("max")
        cutoff = fmax - (fmax - fmin + 1) * last_fraction
        endpoint_rows = work[work["frame"] > cutoff]
    if endpoint_rows.empty:
        logger.warning(
            "compute_endpoint_per_replicate(): das Endfenster (letzte %.0f%% der Frames) ist "
            "fuer '%s' leer - uebersprungen.", 100 * last_fraction, value_col,
        )
        return pd.DataFrame(), pd.DataFrame()

    return summarise_per_replicate(endpoint_rows, value_col, group_cols=group_cols)


def spearman_against_period(
    per_replicate: pd.DataFrame,
    group_cols: Optional[Sequence[str]] = None,
    freq_col: str = "osc_freq",
    value_col: str = "value",
    min_periods: int = MIN_PERIODS_FOR_TREND,
) -> pd.DataFrame:
    """Spearman-Rangkorrelation gegen die Periode, pro Stamm/Oszillationstyp.

    Eingabe ist die CHIP-Ebene (ein Wert je Chip). Bei den Oszillationsdaten
    ist das EIN Wert pro Periode - n ist die Zahl der Perioden (6 bei Glc, 4
    bei pH). Bei n = 6 braucht p < 0.05 ein |rho| >= 0.83: rho ist hier eine
    Effektstaerke, die man berichtet, kein Test, auf den man sich stuetzt.
    Ueber Kammern gerechnet saehe n groesser aus, waere aber Pseudoreplikation -
    die 5 Kammern einer Periode sind eine Kultur auf einem Chip.

    Kontrollen fallen heraus: PosCtrl/NegCtrl haben keine Periode.
    Nicht-numerische Bedingungen (statisch) ebenfalls, ueber period_minutes().

    p_value = NaN, wenn weniger als `min_periods` verschiedene Perioden
    vorliegen - explizit markiert, nicht stillschweigend uebersprungen.
    """
    if per_replicate is None or per_replicate.empty:
        return pd.DataFrame()
    if value_col not in per_replicate.columns or freq_col not in per_replicate.columns:
        return pd.DataFrame()

    df = add_condition_type(per_replicate) if "condition" in per_replicate.columns else per_replicate.copy()
    if "condition_type" in df.columns:
        df = df[~df["condition_type"].isin(CONTROL_TYPES)]

    df = df.assign(_period=period_minutes(df[freq_col])).dropna(subset=["_period", value_col])
    if df.empty:
        logger.info(
            "spearman_against_period(): keine numerischen Perioden (Spalte '%s') - "
            "Trendtest entfaellt (z.B. bei statischen Daten).", freq_col,
        )
        return pd.DataFrame()

    if group_cols is None:
        group_cols = [c for c in ["value_col", "biosensor", "osc_type"] if c in df.columns]
    group_cols = list(group_cols)

    records = []
    for keys, grp in df.groupby(group_cols, dropna=False) if group_cols else [((), df)]:
        keys = keys if isinstance(keys, tuple) else (keys,)
        record = dict(zip(group_cols, keys))
        n_periods = int(grp["_period"].nunique())
        record["n_periods"] = n_periods
        record["n_chips"] = int(len(grp))
        if n_periods < min_periods:
            record["spearman_rho"] = float("nan")
            record["p_value"] = float("nan")
            logger.warning(
                "Spearman fuer %s: nur %d verschiedene Perioden (mindestens %d noetig) - "
                "Trendtest nicht durchfuehrbar.", record, n_periods, min_periods,
            )
        else:
            rho, p = stats.spearmanr(grp["_period"], grp[value_col])
            record["spearman_rho"] = float(rho)
            record["p_value"] = float(p)
        records.append(record)

    out = pd.DataFrame(records)
    if not out.empty and out["p_value"].notna().any():
        n_sig = int((out["p_value"] < 0.05).sum())
        logger.info(
            "Spearman gegen die Periode: %d Gruppen getestet, %d mit p < 0.05 "
            "(monotoner Zusammenhang mit der Periode).", int(out["p_value"].notna().sum()), n_sig,
        )
    return out


def control_trend_check(
    per_chip: pd.DataFrame,
    value_col: str = "value",
    freq_col: str = "osc_freq",
    group_cols: Optional[Sequence[str]] = None,
    min_periods: int = MIN_PERIODS_FOR_TREND,
    strong: float = 0.6,
) -> pd.DataFrame:
    """Laufen die KONTROLLEN einer Struktur mit deren Periode?

    Jede Periode ist eine Struktur (im Code 'chip') mit eigenen PosCtrl- und
    NegCtrl-Kammern in konstantem Medium. Die Kontrollen koennen auf die
    Periode nicht reagieren. Aendern sie sich ueber die Strukturen einer Serie
    trotzdem monoton mit der Periode, traegt die Struktur selbst den Trend
    (Position auf dem physischen Chip, Beladungsreihenfolge, Kulturalter am
    Tag, Stroemung) - und der Trend der Oszillationskammern ist erst dann ein
    Periodeneffekt, wenn er ueber die Kontrollen HINAUSgeht
    (rho_osc_minus_ctrl). Das ist die Differenz-Variante des Bracket-Scores,
    die auch dann etwas sagt, wenn PosCtrl und NegCtrl nicht auseinanderliegen.

    Eine Zeile je (value_col, biosensor, osc_type): rho gegen die Periode fuer
    die Oszillationskammern, PosCtrl, NegCtrl, das Kontrollmittel und
    Osc - Kontrollmittel, mit n_periods und verdict. Ein 'period effect'
    verlangt drei Dinge zugleich: die Oszillationskammern trenden
    (|rho| >= strong), KEINE Kontrollart trendet gleichsinnig, und die
    Differenz Osc - Kontrollen trendet ebenfalls (sonst sind die
    Kontrollschwankungen zwischen den Strukturen so gross wie der Trend).
    """
    if per_chip is None or per_chip.empty or value_col not in per_chip.columns:
        return pd.DataFrame()
    df = add_condition_type(per_chip) if "condition_type" not in per_chip.columns else per_chip.copy()
    df = df.assign(_period=period_minutes(df[freq_col])).dropna(subset=["_period", value_col])
    if df.empty:
        return pd.DataFrame()
    if group_cols is None:
        group_cols = [c for c in ["value_col", "biosensor", "osc_type"] if c in df.columns]
    group_cols = list(group_cols)

    def _rho(x: pd.Series, y: pd.Series) -> float:
        ok = x.notna() & y.notna()
        if ok.sum() < min_periods or y[ok].nunique() < 2:
            return np.nan
        return float(stats.spearmanr(x[ok], y[ok])[0])

    rows = []
    for keys, grp in df.groupby(group_cols, dropna=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        wide = grp.pivot_table(index="_period", columns="condition_type", values=value_col, aggfunc="mean")
        for col in ("Oscillation", "PosCtrl", "NegCtrl"):
            if col not in wide.columns:
                wide[col] = np.nan
        wide["ctrl_mean"] = wide[["PosCtrl", "NegCtrl"]].mean(axis=1)
        wide["osc_minus_ctrl"] = wide["Oscillation"] - wide["ctrl_mean"]
        period = pd.Series(wide.index.to_numpy(dtype=float), index=wide.index)
        rec = dict(zip(group_cols, keys))
        rec.update({
            "n_periods": int(wide["Oscillation"].notna().sum()),
            "rho_osc": _rho(period, wide["Oscillation"]),
            "rho_posctrl": _rho(period, wide["PosCtrl"]),
            "rho_negctrl": _rho(period, wide["NegCtrl"]),
            "rho_ctrl_mean": _rho(period, wide["ctrl_mean"]),
            "rho_osc_minus_ctrl": _rho(period, wide["osc_minus_ctrl"]),
        })
        r_osc, r_diff = rec["rho_osc"], rec["rho_osc_minus_ctrl"]
        # Der staerkste Kontrolltrend zaehlt: schon EINE Kontrollart, die mit
        # der Periode laeuft, belegt einen Struktureffekt - Kontrollen in
        # konstantem Medium koennen auf die Periode nicht reagieren.
        ctrl_rhos = [rec["rho_posctrl"], rec["rho_negctrl"], rec["rho_ctrl_mean"]]
        ctrl_rhos = [r for r in ctrl_rhos if not np.isnan(r)]
        r_ctrl = max(ctrl_rhos, key=abs) if ctrl_rhos else np.nan
        rec["rho_ctrl_strongest"] = r_ctrl
        same_sign_diff = (not np.isnan(r_diff)) and abs(r_diff) >= strong and np.sign(r_diff) == np.sign(r_osc)
        if np.isnan(r_osc) or np.isnan(r_ctrl):
            verdict = "fewer than 3 periods with oscillation and control values"
        elif abs(r_osc) < strong:
            verdict = "no monotone trend of the oscillation chambers"
        elif abs(r_ctrl) >= strong and np.sign(r_ctrl) == np.sign(r_osc):
            verdict = ("structure effect: a constant-medium control trends with the period like the "
                       "oscillation chambers" + ("; Osc - controls still trends, residual period effect on top"
                                                  if same_sign_diff else "; no period effect beyond the controls"))
        elif same_sign_diff:
            verdict = "period effect: oscillation chambers trend, their controls do not, and the difference trends"
        else:
            verdict = ("not robust: oscillation chambers trend, but not after subtracting their controls "
                       "(control swings between structures as large as the trend)")
        rec["verdict"] = verdict
        rows.append(rec)
    out = pd.DataFrame(rows)
    if not out.empty:
        n_struct = int(out["verdict"].str.startswith("structure effect").sum())
        n_period = int(out["verdict"].str.startswith("period effect").sum())
        n_weak = int(out["verdict"].str.startswith("not robust").sum())
        logger.info(
            "Kontroll-Trend-Check (%s): %d Serien mit Struktureffekt (eine Kontrolle laeuft mit der Periode), "
            "%d mit robustem Periodeneffekt (Osc und Osc - Kontrollen trenden, Kontrollen nicht), "
            "%d nicht robust (Osc trendet, die Differenz nicht), %d ohne Osc-Trend.",
            value_col if "value_col" not in out.columns else ", ".join(sorted(out["value_col"].astype(str).unique())),
            n_struct, n_period, n_weak, int(out["verdict"].str.startswith("no monotone").sum()),
        )
    return out


_READOUT_LABELS = {
    "area": "endpoint area", "eccentricity": "endpoint eccentricity",
    "mu_area": "µ_area", "budding_rate_per_h": "budding rate (sparse window)",
    "mu_bud": "µ_bud (births per cell-hour)", "immigration_per_cell_h": "immigration per cell-hour",
}


def plot_control_trend_summary(ctrl_trend: pd.DataFrame, out_path: Path, strong: float = 0.6) -> None:
    """EINE Abbildung fuer den Befund: pro Readout und Serie der Spearman der Oszillationskammern
    gegen die Periode (x) und der der staerksten Kontrolle derselben Strukturen (y). Punkte nahe der
    Diagonale: die Struktur traegt beide. Farbe = Stamm, Marker = Readout, hohl = 'not robust'.
    Die Zonen tragen die Verdict-Regel: dunkleres Grau in den Ecken = |rho_osc| >= strong UND die
    staerkste Kontrolle gleichsinnig >= strong (Struktureffekt); helles Band = |rho_osc| >= strong
    ohne gleichsinnige Kontrolle (Periodeneffekt oder 'not robust', je nach Differenz).
    Eine Facette je Oszillationstyp."""
    from matplotlib.patches import Rectangle

    if ctrl_trend is None or ctrl_trend.empty:
        return
    df = ctrl_trend.dropna(subset=["rho_osc", "rho_ctrl_strongest"]).copy()
    if df.empty:
        logger.warning("plot_control_trend_summary(): keine Serie mit >= 3 Perioden - uebersprungen.")
        return
    if "osc_type" not in df.columns:
        df["osc_type"] = "all"
    osc_types = sorted(df["osc_type"].dropna().astype(str).unique())
    readouts = list(dict.fromkeys(df["value_col"].astype(str)))
    marker_cycle = ["o", "s", "^", "D", "v", "P", "X", "*", "<", ">"]
    markers = {r: marker_cycle[i % len(marker_cycle)] for i, r in enumerate(readouts)}
    zone_structure, zone_period = "#e3e3e3", "#f3f3f3"

    fig, axes = plt.subplots(1, len(osc_types), figsize=(3.6 * len(osc_types) + 0.4, 3.9), squeeze=False, sharey=True)
    for ax, ot in zip(axes[0], osc_types):
        sub = df[df["osc_type"].astype(str) == ot]
        ax.grid(False)
        for sx in (1, -1):
            y0 = -1.05 if sx > 0 else -strong
            ax.add_patch(Rectangle((min(sx * strong, sx * 1.05), y0), 1.05 - strong, 1.05 + strong,
                                   color=zone_period, lw=0, zorder=0))
        for sx, sy in ((1, 1), (-1, -1)):
            ax.add_patch(Rectangle((min(sx * strong, sx * 1.05), min(sy * strong, sy * 1.05)),
                                   1.05 - strong, 1.05 - strong, color=zone_structure, lw=0, zorder=0))
        ax.text(0.83, 1.0, "structure\neffect", ha="center", va="top", fontsize=6.5, color=INK_MUTED, zorder=1)
        ax.text(0.83, -0.02, "period\neffect", ha="center", va="top", fontsize=6.5, color=INK_MUTED, zorder=1)
        ax.plot([-1.05, 1.05], [-1.05, 1.05], color=INK_MUTED, lw=0.7, ls="--", zorder=1)
        ax.axhline(0, color="#d0d0d0", lw=0.6, zorder=1)
        ax.axvline(0, color="#d0d0d0", lw=0.6, zorder=1)
        for _, r in sub.iterrows():
            color = strain_color(r.get("biosensor", ""))
            hollow = str(r.get("verdict", "")).startswith("not robust")
            ax.scatter(r["rho_osc"], r["rho_ctrl_strongest"], marker=markers[str(r["value_col"])], s=46,
                       facecolor=SURFACE if hollow else color, edgecolor=edge_color(color), linewidth=0.9, zorder=3)
        ax.set_xlim(-1.05, 1.05)
        ax.set_ylim(-1.05, 1.05)
        ax.set_aspect("equal")
        ax.set_xticks([-1, -0.5, 0, 0.5, 1]); ax.set_yticks([-1, -0.5, 0, 0.5, 1])
        panel_title(ax, f"{ot}  (n = {len(sub)} readout × strain series)")
        ax.set_xlabel("Spearman ρ vs period: oscillation chambers")
    axes[0][0].set_ylabel("Spearman ρ vs period:\nstrongest control of the same structures")

    handles = [Line2D([], [], marker=markers[r], color=INK_SOFT, ls="", markersize=6, label=_READOUT_LABELS.get(r, r))
               for r in readouts]
    handles += strain_handles(df["biosensor"].unique()) if "biosensor" in df.columns else []
    handles.append(Line2D([], [], marker="o", mfc=SURFACE, mec=INK_SOFT, ls="", markersize=6,
                          label="not robust: trend vanishes after subtracting the controls"))
    legend_below(fig, handles, ncol=3, y=0.0)
    finish(fig, out_path, logger)


def within_culture_trend(
    per_chip: pd.DataFrame,
    value_col: str = "value",
    freq_col: str = "osc_freq",
    culture_col: str = "culture",
    group_cols: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    """Aenderung von der kuerzesten zur laengsten Periode INNERHALB einer Kultur.

    Kultur = physischer Chip = Vorkultur = Datum; er traegt 2-3 Strukturen mit
    je einer Periode. Zwischen Kulturen ist ein Periodenvergleich mit der
    Kultur konfundiert, innerhalb einer Kultur nicht. Eine Zeile je
    (value_col, biosensor, osc_type, culture): die Perioden, die Werte der
    Oszillationskammern und ihrer Kontrollen an der kuerzesten und laengsten
    Periode und die Differenzen (delta_* < 0 = faellt mit der Periode).
    """
    if per_chip is None or per_chip.empty or value_col not in per_chip.columns or culture_col not in per_chip.columns:
        return pd.DataFrame()
    df = add_condition_type(per_chip) if "condition_type" not in per_chip.columns else per_chip.copy()
    df = df.assign(_period=period_minutes(df[freq_col])).dropna(subset=["_period", value_col, culture_col])
    if df.empty:
        return pd.DataFrame()
    if group_cols is None:
        group_cols = [c for c in ["value_col", "biosensor", "osc_type"] if c in df.columns]
    group_cols = list(group_cols) + [culture_col]
    rows = []
    for keys, grp in df.groupby(group_cols, dropna=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        wide = grp.pivot_table(index="_period", columns="condition_type", values=value_col, aggfunc="mean").sort_index()
        if "Oscillation" not in wide.columns or wide["Oscillation"].notna().sum() < 2:
            continue
        for col in ("PosCtrl", "NegCtrl"):
            if col not in wide.columns:
                wide[col] = np.nan
        wide["ctrl_mean"] = wide[["PosCtrl", "NegCtrl"]].mean(axis=1)
        osc = wide["Oscillation"].dropna()
        first, last = osc.index.min(), osc.index.max()
        rec = dict(zip(group_cols, keys))
        rec.update({
            "periods": " ".join(f"{p:g}" for p in osc.index),
            "n_structures": int(len(osc)),
            "osc_shortest": float(osc.loc[first]), "osc_longest": float(osc.loc[last]),
            "delta_osc": float(osc.loc[last] - osc.loc[first]),
            "delta_ctrl": float(wide.loc[last, "ctrl_mean"] - wide.loc[first, "ctrl_mean"]),
        })
        rec["delta_osc_minus_ctrl"] = rec["delta_osc"] - rec["delta_ctrl"]
        rows.append(rec)
    out = pd.DataFrame(rows)
    if not out.empty:
        n = int(out["delta_osc"].notna().sum())
        logger.info(
            "Innerhalb der Kulturen (%d Kulturen mit >= 2 Perioden): Osc faellt zur laengeren Periode in %d, "
            "Kontrollen in %d, Osc - Kontrollen in %d.",
            n, int((out["delta_osc"] < 0).sum()), int((out["delta_ctrl"] < 0).sum()),
            int((out["delta_osc_minus_ctrl"] < 0).sum()),
        )
    return out


def bracket_normalise(
    per_chip: pd.DataFrame,
    degenerate_k: float = 2.0,
) -> pd.DataFrame:
    """Endzustand jeder Oszillationsbedingung relativ zu den Kontrollen IHRES Chips.

    score = (osc - NegCtrl) / (PosCtrl - NegCtrl)
        0 = wie durchgehend Starvation, 1 = wie durchgehend Feast.

    Warum das die richtige Groesse ist: jede Periode ist ein eigener Chip aus
    einer eigenen Vorkultur. Der Rohwert einer Periode enthaelt damit den
    Chip-/Tages-/Kultur-Effekt. Die Kontrollen liegen auf DEMSELBEN Chip und
    tragen denselben Effekt - die Normierung entfernt ihn, ohne dass man den
    Chip kennen muesste. Das ist die gepaarte Auswertung, die dieses Design
    umsonst mitliefert.

    Entartung: wenn |PosCtrl - NegCtrl| kleiner ist als degenerate_k mal die
    Kammer-Streuung der Kontrollen auf diesem Chip, trennt das Bracket nichts
    und der Score ist Rauschen geteilt durch Rauschen. Solche Chips werden mit
    bracket_degenerate=True markiert, im Plot hohl gezeichnet und fallen aus
    dem Trendtest - und genau diese Chips sind der Befund von Abschnitt 4.
    """
    if per_chip is None or per_chip.empty:
        return pd.DataFrame()
    df = add_condition_type(per_chip) if "condition_type" not in per_chip.columns else per_chip.copy()
    key = [c for c in ["value_col", "biosensor", "osc_type", "osc_freq", "chip"] if c in df.columns]
    rows = []
    for keys, grp in df.groupby(key, dropna=False):
        by = grp.set_index("condition_type")
        rec = dict(zip(key, keys if isinstance(keys, tuple) else (keys,)))
        for ct, name in (("Oscillation", "osc"), ("PosCtrl", "pos"), ("NegCtrl", "neg")):
            rec[name] = float(by.loc[ct, "value"]) if ct in by.index else float("nan")
            rec[f"{name}_sd_chamber"] = float(by.loc[ct, "sd_chamber"]) if ct in by.index and "sd_chamber" in by.columns else float("nan")
        span = rec["pos"] - rec["neg"]
        noise = np.nanmean([rec["pos_sd_chamber"], rec["neg_sd_chamber"]])
        rec["bracket_span"] = span
        rec["bracket_degenerate"] = bool(np.isnan(span) or (noise == noise and abs(span) < degenerate_k * noise))
        rec["value"] = (rec["osc"] - rec["neg"]) / span if span and span == span and not rec["bracket_degenerate"] else float("nan")
        rec["condition"] = "bracket_score"
        rec["condition_type"] = "Oscillation"
        rows.append(rec)
    out = pd.DataFrame(rows)
    if not out.empty:
        n_deg = int(out["bracket_degenerate"].sum())
        if n_deg:
            logger.warning(
                "bracket_normalise(): %d von %d Chips haben ein entartetes Bracket "
                "(|PosCtrl - NegCtrl| < %.1f x Kammer-Streuung) - Score dort NaN, aus dem "
                "Trendtest ausgeschlossen. Das ist der Befund aus Abschnitt 4, kein Fehler.",
                n_deg, len(out), degenerate_k,
            )
    return out


def plot_endpoint_vs_period(
    per_chip: pd.DataFrame,
    out_path: Path,
    value_col: str,
    trend: Optional[pd.DataFrame] = None,
    score: Optional[pd.DataFrame] = None,
    score_trend: Optional[pd.DataFrame] = None,
    ylabel: Optional[str] = None,
    strain_order: Optional[Sequence[str]] = None,
    title: Optional[str] = None,
) -> None:
    """Endzustand gegen die Periode: EINE Struktur pro Periode, mit IHREN Kontrollen.

    Obere Reihe: pro Periode der Mittelwert der Oszillationskammern (Fehlerbalken = Kammern dieser
    Struktur, technisch) als gefuellter Kreis mit Linie in der Stammfarbe; an derselben x-Position
    die Feast- (Dreieck hoch, gefuellt) und Famine-Kontrolle (Dreieck runter, hohl) DERSELBEN
    Struktur. Duenne graue Referenzlinien: das ueber alle Strukturen des Stamms gepoolte Mittel je
    Kontrollart (gestrichelt Feast, gepunktet Famine) - wie weit die Kontrollen selbst wandern,
    zeigen die Dreiecke.

    Untere Reihe: der bracket-normierte Score (bracket_normalise()), 0 = Famine-Kontrolle,
    1 = Feast-Kontrolle; Strukturen mit entartetem Bracket hohl bei 0.5.

    Eine Facette pro Stamm in der Stammfarbe (config.STRAIN_COLORS); Spearman-Wert in der Facette,
    n = Zahl der Strukturen = Zahl der Perioden. Schritt 20/24 nutzen dieselbe Abbildung fuer die
    Knospungsrate und µ_bud (title/ylabel).
    """
    if per_chip is None or per_chip.empty:
        logger.warning("plot_endpoint_vs_period(): keine Daten fuer '%s' - uebersprungen.", value_col)
        return
    df = add_condition_type(per_chip) if "condition_type" not in per_chip.columns else per_chip.copy()
    df = df.assign(_period=period_minutes(df["osc_freq"]))
    osc = df[(df["condition_type"] == "Oscillation")].dropna(subset=["_period", "value"])
    if osc.empty:
        logger.warning("plot_endpoint_vs_period(): keine Oszillations-Chips mit numerischer Periode "
                       "fuer '%s' - uebersprungen.", value_col)
        return
    strains = [b for b in (strain_order or ordered_strains(osc["biosensor"].unique())) if b in set(osc["biosensor"])]
    has_score = score is not None and not score.empty
    n_rows = 2 if has_score else 1
    periods_all = sorted(osc["_period"].unique())
    fig, axes = plt.subplots(n_rows, len(strains), figsize=(2.7 * len(strains) + 0.5, 2.7 * n_rows + 0.4),
                             squeeze=False, sharex=True)
    any_degenerate = False

    for j, strain in enumerate(strains):
        color = strain_color(strain)
        ax = axes[0][j]
        sub = df[df["biosensor"] == strain]
        # Gepooltes Kontrollmittel je Art als duenne Referenzlinie (Feast gestrichelt, Famine gepunktet).
        for ct in CONTROL_TYPES:
            vals = sub.loc[sub["condition_type"] == ct, "value"].dropna()
            if len(vals) >= 2:
                ax.axhline(vals.mean(), color=INK_MUTED, linewidth=0.8, linestyle=CONTROL_LINESTYLES[ct], zorder=1)
        # Kontrollen DIESER Struktur an der x-Position ihrer Periode.
        for ct in ("NegCtrl", "PosCtrl"):
            c = sub[sub["condition_type"] == ct].dropna(subset=["_period", "value"]).sort_values("_period")
            if not c.empty:
                ax.scatter(c["_period"], c["value"], **marker_kwargs(ct, color, size=30))
        o = sub[sub["condition_type"] == "Oscillation"].dropna(subset=["_period", "value"]).sort_values("_period")
        yerr = o["sd_chamber"].fillna(0.0) if "sd_chamber" in o.columns else None
        ax.errorbar(o["_period"], o["value"], yerr=yerr, **errorbar_kwargs("Oscillation", color))
        _log_period_axis(ax, periods_all)
        panel_title(ax, f"{strain}  (n = {o['_period'].nunique()} structures)")
        if j == 0:
            ax.set_ylabel(ylabel or f"{value_col}\n(endpoint)")
        _annotate_trend(ax, trend, strain, value_col)

        if has_score:
            ax2 = axes[1][j]
            sc = score[(score["biosensor"] == strain)].assign(_period=lambda d: period_minutes(d["osc_freq"]))
            sc = sc.dropna(subset=["_period"]).sort_values("_period")
            good = sc[~sc["bracket_degenerate"]].dropna(subset=["value"])
            bad = sc[sc["bracket_degenerate"]]
            ax2.axhline(0, color=INK_MUTED, linewidth=0.8, linestyle=CONTROL_LINESTYLES["NegCtrl"], zorder=1)
            ax2.axhline(1, color=INK_MUTED, linewidth=0.8, linestyle=CONTROL_LINESTYLES["PosCtrl"], zorder=1)
            if not good.empty:
                ax2.plot(good["_period"], good["value"], marker="o", markersize=5.5, linewidth=1.4, color=color,
                         markerfacecolor=color, markeredgecolor=edge_color(color), markeredgewidth=0.8, zorder=3)
            if not bad.empty:
                any_degenerate = True
                ax2.scatter(bad["_period"], np.full(len(bad), 0.5), marker="o", s=34, facecolor=SURFACE,
                            edgecolor=edge_color(color), linewidth=0.9, zorder=3)
            _log_period_axis(ax2, periods_all)
            ax2.set_ylim(-0.3, 1.3)
            ax2.set_xlabel("feast/famine cycle period [min]")
            if j == 0:
                ax2.set_ylabel("bracket score\n0 = famine control, 1 = feast control")
            _annotate_trend(ax2, score_trend, strain, value_col)
        else:
            ax.set_xlabel("feast/famine cycle period [min]")

    handles = control_handles(INK_SOFT)
    handles[0].set_label("oscillation chambers (mean ± SD of the structure's chambers)")
    handles += [Line2D([], [], color=INK_MUTED, linewidth=0.9, linestyle=CONTROL_LINESTYLES["PosCtrl"],
                       label="feast controls, mean over structures"),
                Line2D([], [], color=INK_MUTED, linewidth=0.9, linestyle=CONTROL_LINESTYLES["NegCtrl"],
                       label="famine controls, mean over structures")]
    if any_degenerate:
        handles.append(Line2D([], [], marker="o", linestyle="", markersize=6, markerfacecolor=SURFACE,
                              markeredgecolor=INK_SOFT, label="bracket undefined (controls do not separate)"))
    legend_below(fig, handles, ncol=3, y=0.0)
    if title:
        fig.suptitle(title, fontsize=9.5, y=1.0)
    finish(fig, out_path, logger)


def _log_period_axis(ax, periods) -> None:
    """Log-x mit expliziten Ticks auf den Perioden; Minor-Ticks aus (sonst '2 x 10^0')."""
    ax.set_xscale("log")
    ax.set_xticks(periods)
    ax.set_xticklabels([f"{p:g}" for p in periods])
    ax.xaxis.set_minor_locator(mticker.NullLocator())
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax.tick_params(axis="x", labelsize=8)


def _annotate_trend(ax, trend, strain, value_col) -> None:
    if trend is None or trend.empty:
        return
    mask = pd.Series(True, index=trend.index)
    for col_name, want in (("biosensor", strain), ("value_col", value_col)):
        if col_name in trend.columns:
            mask &= trend[col_name] == want
    hit = trend[mask]
    if len(hit) == 1 and pd.notna(hit.iloc[0]["spearman_rho"]):
        rho, p, n = hit.iloc[0]["spearman_rho"], hit.iloc[0]["p_value"], int(hit.iloc[0].get("n_chips", 0))
        ax.text(0.03, 0.97, f"ρ = {rho:+.2f}, p = {p:.2f}, n = {n}", transform=ax.transAxes,
                va="top", ha="left", fontsize=7.5, color=INK_SOFT, zorder=6,
                bbox=dict(boxstyle="round,pad=0.2", facecolor=SURFACE, edgecolor="none", alpha=0.85))


