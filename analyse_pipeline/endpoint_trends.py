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
    replicate_col: str = "replicate",
    chamber_col: str = "exp_id",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Dreistufige Aggregation: Kammer -> Replikat -> Bedingung.

    Das ist die Aggregation, die die vorhandenen Summary-Funktionen NICHT
    machen: summarise_growth_rate() und summarise_area_growth() poolen
    einzelne Zellen, und aggregate_robustness_over_replicates() kommt per
    Default nur bis auf Kammer-Ebene (exp_id enthaelt 'chamber', siehe
    data_loading.py). Hier ist die Einheit das biologische Replikat.

    Rueckgabe: (pro Replikat, pro Bedingung). Die Replikat-Ebene ist die
    Grundlage fuer den Spearman-Test, die Bedingungs-Ebene die fuer den Plot.
    """
    if df is None or df.empty or value_col not in df.columns:
        return pd.DataFrame(), pd.DataFrame()

    if group_cols is None:
        group_cols = [c for c in ["biosensor", "osc_type", "osc_freq", "condition"]
                      if c in df.columns]
    group_cols = list(group_cols)

    work = df.dropna(subset=[value_col])
    if work.empty:
        return pd.DataFrame(), pd.DataFrame()

    # Stufe 1: pro Kammer
    if chamber_col in work.columns:
        per_chamber = (
            work.groupby(group_cols + [c for c in (replicate_col, chamber_col)
                                       if c in work.columns], dropna=False)[value_col]
            .mean().reset_index()
        )
    else:
        per_chamber = work

    # Stufe 2: pro Replikat
    if replicate_col not in per_chamber.columns:
        logger.warning(
            "summarise_per_replicate(): Spalte '%s' fehlt - es kann nicht auf biologische "
            "Replikate aggregiert werden; die Fehlerbalken beschreiben dann Kammern.",
            replicate_col,
        )
        per_replicate = per_chamber
        rep_key = [c for c in (chamber_col,) if c in per_chamber.columns]
    else:
        per_replicate = (
            per_chamber.groupby(group_cols + [replicate_col], dropna=False)[value_col]
            .mean().reset_index()
        )
        rep_key = [replicate_col]

    # Stufe 3: ueber Replikate
    summary = (
        per_replicate.groupby(group_cols, dropna=False)[value_col]
        .agg(mean="mean", sd="std",
             sem=lambda s: s.std(ddof=1) / np.sqrt(len(s)) if len(s) > 1 else 0.0,
             n_replicates="count")
        .reset_index()
    )
    if "condition" in summary.columns:
        summary = add_condition_type(summary)
    if "condition" in per_replicate.columns:
        per_replicate = add_condition_type(per_replicate)
    per_replicate = per_replicate.rename(columns={value_col: "value"})
    per_replicate["value_col"] = value_col
    summary["value_col"] = value_col
    logger.debug("summarise_per_replicate(%s): Replikat-Einheit = %s", value_col, rep_key)
    return per_replicate, summary


def compute_endpoint_per_replicate(
    cells: pd.DataFrame,
    value_col: str,
    last_fraction: float = 0.25,
    group_cols: Optional[Sequence[str]] = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Endzustand einer Bedingung: die letzten `last_fraction` der Frames.

    Das Fenster wird PRO KAMMER bestimmt (relativ zu deren eigenem
    Frame-Bereich), nicht global - Kammern koennen unterschiedlich lang
    aufgenommen sein.

    Anschliessend laeuft dieselbe dreistufige Aggregation wie in
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
                           alpha=0.9, zorder=2, label=f"{ctype} (n={int(row.get('n_replicates', 0))})")

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
