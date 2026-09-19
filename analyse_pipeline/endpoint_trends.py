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

CONTROL_TYPES = ["NegCtrl", "PosCtrl"]
CONTROL_COLORS = {"NegCtrl": "#4C78A8", "PosCtrl": "#E45756"}

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

    Eingabe ist die REPLIKAT-Ebene (ein Wert je Replikat und Periode) - nicht
    die Zellebene und nicht die bereits gemittelte Summary. Ueber Zellen
    gerechnet haengt das p fast nur an der Zellzahl; auf der gemittelten
    Summary bliebe pro Periode ein einziger Punkt und die Streuung zwischen
    den Replikaten waere unsichtbar.

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
        record["n_replicate_values"] = int(len(grp))
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


def plot_endpoint_vs_period(
    summary: pd.DataFrame,
    out_path: Path,
    value_col: str,
    trend: Optional[pd.DataFrame] = None,
    ylabel: Optional[str] = None,
) -> None:
    """Endzustand gegen die Periode, mit PosCtrl/NegCtrl als Referenzbaender.

    Die x-Achse ist LOGARITHMISCH: die Perioden sind geometrisch gestuft
    (jeweils Faktor 2 von 0.75 bis 24 min, also ein 32-facher Dosisbereich).
    Linear dargestellt draengen sich fuenf der sechs Bedingungen links
    zusammen, und die Dosis-Wirkung waere nicht ablesbar.

    Die Kontrollen erscheinen als waagerechte Baender (Mittelwert +- SEM ueber
    ihre Replikate) statt als Punkte auf der Periodenachse: sie HABEN keine
    Periode, und der Ordnername ihres Batches ist keine Behandlung. Genau so
    sind sie auch gemeint - als Bezugsrahmen, gegen den die
    Oszillationsbedingungen gelesen werden.
    """
    if summary is None or summary.empty:
        logger.warning("plot_endpoint_vs_period(): keine Daten fuer '%s' - uebersprungen.", value_col)
        return

    df = add_condition_type(summary) if "condition" in summary.columns else summary.copy()
    df = df.assign(_period=period_minutes(df["osc_freq"]))

    is_control = df["condition_type"].isin(CONTROL_TYPES) if "condition_type" in df.columns \
        else pd.Series(False, index=df.index)
    osc = df[~is_control].dropna(subset=["_period", "mean"])
    controls = df[is_control]

    if osc.empty:
        logger.warning(
            "plot_endpoint_vs_period(): keine Oszillationsbedingungen mit numerischer Periode "
            "fuer '%s' - uebersprungen (bei statischen Daten erwartet).", value_col,
        )
        return

    osc_types = sorted(osc["osc_type"].dropna().unique())
    biosensors = sorted(osc["biosensor"].dropna().unique())
    fig, axes = plt.subplots(
        len(osc_types), len(biosensors),
        figsize=(4.2 * len(biosensors), 3.6 * len(osc_types)),
        squeeze=False, sharex=True,
    )

    for i, osc_type in enumerate(osc_types):
        for j, biosensor in enumerate(biosensors):
            ax = axes[i][j]
            sub = osc[(osc["osc_type"] == osc_type) & (osc["biosensor"] == biosensor)].sort_values("_period")

            # Kontrollbaender zuerst, damit die Datenpunkte darueber liegen.
            csub = controls[(controls["osc_type"] == osc_type) & (controls["biosensor"] == biosensor)]
            for _, row in csub.iterrows():
                ctype = row["condition_type"]
                col = CONTROL_COLORS.get(ctype, "#999999")
                spread = row.get("sem", 0.0)
                spread = 0.0 if pd.isna(spread) else spread
                ax.axhspan(row["mean"] - spread, row["mean"] + spread,
                           color=col, alpha=0.16, zorder=1)
                ax.axhline(row["mean"], color=col, linewidth=1.2, linestyle="--",
                           alpha=0.9, zorder=2,
                           label=f"{ctype} (n={int(row.get('n_units', 0))} {row.get('error_unit', '')})")

            if not sub.empty:
                ax.errorbar(
                    sub["_period"], sub["mean"],
                    yerr=sub["sem"].fillna(0.0) if "sem" in sub.columns else None,
                    marker="o", markersize=5.5, linewidth=1.5, capsize=3,
                    color="#333333", zorder=3, label="Oscillation",
                )
            # Log-Achse mit EXPLIZITEN Ticks auf den echten Perioden. Die
            # Minor-Ticks muessen dabei abgeschaltet werden: matplotlib
            # beschriftet sie auf einer Log-Achse per Default in
            # Exponentialschreibweise ("2 x 10^0"), was sich mit den
            # Major-Labels zu unlesbarem Text ueberlagert.
            periods_present = sorted(osc["_period"].unique())
            ax.set_xscale("log")
            ax.set_xticks(periods_present)
            ax.set_xticklabels([f"{p:g}" for p in periods_present])
            ax.xaxis.set_minor_locator(mticker.NullLocator())
            ax.xaxis.set_minor_formatter(mticker.NullFormatter())
            ax.tick_params(axis="x", labelsize=8)
            ax.grid(alpha=0.22, linewidth=0.6)
            ax.set_title(f"{osc_type} | {biosensor}", fontsize=10)

            # Trendtest in die Facette schreiben, zu der er gehoert - so kann
            # die Abbildung nicht von ihrer Statistik getrennt werden.
            if trend is not None and not trend.empty:
                mask = pd.Series(True, index=trend.index)
                for col_name, want in (("biosensor", biosensor), ("osc_type", osc_type),
                                       ("value_col", value_col)):
                    if col_name in trend.columns:
                        mask &= trend[col_name] == want
                hit = trend[mask]
                if len(hit) == 1 and pd.notna(hit.iloc[0]["spearman_rho"]):
                    rho = hit.iloc[0]["spearman_rho"]
                    p = hit.iloc[0]["p_value"]
                    ax.annotate(
                        f"Spearman ρ = {rho:+.2f}\np = {p:.3g}",
                        xy=(0.03, 0.97), xycoords="axes fraction", va="top", ha="left",
                        fontsize=8, bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                                              edgecolor="#BBBBBB", alpha=0.9),
                    )

            if i == len(osc_types) - 1:
                ax.set_xlabel("Feast/famine cycle period [min]")
            if j == 0:
                ax.set_ylabel(ylabel or f"{value_col}\n(endpoint, mean over replicates)")

    handles, labels = axes[0][0].get_legend_handles_labels()
    if handles:
        seen, uh, ul = set(), [], []
        for h, l in zip(handles, labels):
            if l not in seen:
                seen.add(l); uh.append(h); ul.append(l)
        fig.legend(uh, ul, loc="lower center", ncol=min(len(ul), 4), bbox_to_anchor=(0.5, -0.06),
                   frameon=False, fontsize=9)

    fig.suptitle(
        f"{value_col}: cumulative endpoint vs cycle period\n"
        "(point = mean over biological replicates ± SEM; bands = constant-medium controls)",
        y=1.02,
    )
    fig.text(
        0.5, -0.13,
        "Individual cycles are below the sampling limit and are NOT resolved — only the "
        "cumulative endpoint after the full run is interpretable.\nControls have no period and "
        "are drawn as reference bands. Spearman tests monotonicity across periods on "
        "replicate-level values.",
        ha="center", fontsize=8,
    )
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=180)
    plt.close(fig)
    logger.info("Plot gespeichert: %s", out_path.name)
