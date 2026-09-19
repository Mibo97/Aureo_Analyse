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
) -> None:
    """Endzustand gegen die Periode: EIN Chip pro Periode, mit SEINEN Kontrollen.

    Obere Reihe (Rohwert): pro Periode der Chip-Mittelwert ueber die
    Oszillationskammern (Fehlerbalken = Kammern dieses Chips, technisch), und
    an derselben x-Position die PosCtrl-/NegCtrl-Mittelwerte DESSELBEN Chips
    als kleine Marker. Dahinter, blass, das ueber alle Chips gepoolte
    Kontrollband (Mittelwert +- SD ueber Chips) als Bezugsrahmen. So sieht man
    beides: wo die Periode relativ zu ihren eigenen Kontrollen liegt, und wie
    stark die Kontrollen selbst von Chip zu Chip wandern.

    Untere Reihe: der bracket-normierte Score (bracket_normalise()), 0 = wie
    Starvation, 1 = wie Feast. Entartete Chips hohl.

    Eine Facette pro Stamm - die Staemme bleiben getrennt (n = 1 Serie je
    Stamm). Der Spearman-Wert steht in der Facette, zu der er gehoert, mit
    n = Zahl der Chips = Zahl der Perioden.
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
    strains = [b for b in (strain_order or sorted(osc["biosensor"].unique()))
               if b in set(osc["biosensor"])]
    has_score = score is not None and not score.empty
    n_rows = 2 if has_score else 1
    degenerate_labelled = False  # Legendeneintrag an der ERSTEN Facette, die einen hohlen Punkt hat
    fig, axes = plt.subplots(n_rows, len(strains), figsize=(3.9 * len(strains), 3.4 * n_rows),
                             squeeze=False, sharex=True)

    for j, strain in enumerate(strains):
        ax = axes[0][j]
        sub = df[df["biosensor"] == strain]
        # Gepooltes Kontrollband ueber alle Chips dieses Stamms (Bezugsrahmen).
        for ct in CONTROL_TYPES:
            vals = sub.loc[sub["condition_type"] == ct, "value"].dropna()
            if len(vals) >= 2:
                ax.axhspan(vals.mean() - vals.std(), vals.mean() + vals.std(),
                           color=CONTROL_COLORS[ct], alpha=0.10, zorder=0)
                ax.axhline(vals.mean(), color=CONTROL_COLORS[ct], linewidth=0.9, linestyle=":",
                           alpha=0.7, zorder=1)
        # Kontrollen DIESES Chips an der x-Position seiner Periode.
        for ct, marker in (("NegCtrl", "v"), ("PosCtrl", "^")):
            c = sub[sub["condition_type"] == ct].dropna(subset=["_period", "value"]).sort_values("_period")
            if not c.empty:
                ax.scatter(c["_period"], c["value"], marker=marker, s=34, color=CONTROL_COLORS[ct],
                           edgecolor="white", linewidth=0.5, zorder=3,
                           label=f"{ct} of the same chip" if j == 0 else None)
        o = sub[sub["condition_type"] == "Oscillation"].dropna(subset=["_period", "value"]).sort_values("_period")
        ax.errorbar(o["_period"], o["value"],
                    yerr=o["sd_chamber"].fillna(0.0) if "sd_chamber" in o.columns else None,
                    marker="o", markersize=5.5, linewidth=1.4, capsize=3, color="#333333", zorder=4,
                    label="Oscillation (mean of chambers on the chip;\nbar = chamber SD, technical)" if j == 0 else None)
        _log_period_axis(ax, sorted(osc["_period"].unique()))
        ax.set_title(f"{strain}  (n = {o['_period'].nunique()} chips)", fontsize=10)
        ax.grid(alpha=0.22, linewidth=0.6)
        if j == 0:
            ax.set_ylabel(ylabel or f"{value_col}\n(endpoint)")
        _annotate_trend(ax, trend, strain, value_col)

        if has_score:
            ax2 = axes[1][j]
            sc = score[(score["biosensor"] == strain)].assign(_period=lambda d: period_minutes(d["osc_freq"]))
            sc = sc.dropna(subset=["_period"]).sort_values("_period")
            good = sc[~sc["bracket_degenerate"]].dropna(subset=["value"])
            bad = sc[sc["bracket_degenerate"]]
            ax2.axhspan(0, 1, color="#999999", alpha=0.06, zorder=0)
            ax2.axhline(0, color=CONTROL_COLORS["NegCtrl"], linewidth=0.9, linestyle=":", alpha=0.8)
            ax2.axhline(1, color=CONTROL_COLORS["PosCtrl"], linewidth=0.9, linestyle=":", alpha=0.8)
            if not good.empty:
                ax2.plot(good["_period"], good["value"], marker="o", markersize=5.5, linewidth=1.4,
                         color="#333333", zorder=3)
            if not bad.empty:
                ax2.scatter(bad["_period"], np.full(len(bad), 0.5), marker="o", s=40, facecolor="white",
                            edgecolor="#333333", linewidth=1.2, zorder=3,
                            label=None if degenerate_labelled else "bracket degenerate on this chip (score undefined)")
                degenerate_labelled = True
            _log_period_axis(ax2, sorted(osc["_period"].unique()))
            ax2.set_ylim(-0.3, 1.3)
            ax2.grid(alpha=0.22, linewidth=0.6)
            ax2.set_xlabel("Feast/famine cycle period [min]")
            if j == 0:
                ax2.set_ylabel("Bracket score\n0 = NegCtrl, 1 = PosCtrl (same chip)")
            _annotate_trend(ax2, score_trend, strain, value_col)
        else:
            ax.set_xlabel("Feast/famine cycle period [min]")

    handles, labels = [], []
    for row in axes:
        for ax_ in row:
            h, l = ax_.get_legend_handles_labels()
            for hh, ll in zip(h, l):
                if ll not in labels:
                    handles.append(hh); labels.append(ll)
    if handles:
        fig.legend(handles, labels, loc="lower center", ncol=min(len(labels), 3),
                   bbox_to_anchor=(0.5, -0.10 if has_score else -0.16), frameon=False, fontsize=8)
    fig.suptitle(
        f"{value_col}: cumulative endpoint vs cycle period — one chip per period\n"
        "(shaded = PosCtrl/NegCtrl pooled over chips, mean ± SD; markers = controls of that chip)",
        y=1.02,
    )
    fig.text(
        0.5, -0.20 if has_score else -0.26,
        "Each period is one chip from one preculture: n = 1 biological replicate per point. Error bars "
        "are chambers on that chip (technical).\nIndividual cycles are below the sampling limit; only "
        "the cumulative endpoint is interpretable. Spearman on chip means: an effect size, not a test "
        "to lean on at n ≤ 6.",
        ha="center", fontsize=8,
    )
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=180)
    plt.close(fig)
    logger.info("Plot gespeichert: %s", out_path.name)


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
        ax.annotate(f"ρ = {rho:+.2f}, p = {p:.2g}\n(n = {n} chips)",
                    xy=(0.03, 0.97), xycoords="axes fraction", va="top", ha="left", fontsize=8,
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#BBBBBB", alpha=0.92))
