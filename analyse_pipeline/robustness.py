"""
robustness.py
==============
Robustness-Quantifizierung nach Blöbaum et al. 2024 (Microbial Cell
Factories), basierend auf der Fano-Factor-Formel aus Trivellin et al. 2022:

    R = -sigma^2 / x_bar * 1/m                                      (Eq. 1)

wobei sigma/x_bar Standardabweichung/Mittelwert EINER bestimmten
Verteilung sind (je nach R(t) oder R(p) unterschiedlich definiert, siehe
unten), und m der Mittelwert der Funktion ÜBER ALLE Zeitpunkte/Bedingungen
ist (ein globaler Normalisierungsfaktor).

WICHTIG: R ist eine RELATIVE Größe - der Wert ändert sich, wenn weitere
Replikate oder Bedingungen hinzukommen (weil sich 'm' ändert). R ist NICHT
dazu gedacht, in Isolation interpretiert zu werden, sondern im Vergleich
zwischen Bedingungen INNERHALB EIN UND DESSELBEN Analyse-Laufs. Wird der
Datensatz erweitert (neue Bedingung, neues Replikat), ändern sich ALLE
zuvor berechneten R-Werte erneut - sie sind nicht stabil/reproduzierbar
über verschiedene Analyse-Läufe mit unterschiedlichem Datenumfang hinweg.

ZWEI VARIANTEN, ZWEI BEDEUTUNGEN (siehe Paper Methods):

  R(t) - Robustheit ÜBER DIE ZEIT ("wie stabil ist diese Funktion für eine
         Bedingung/Zelle über den Beobachtungszeitraum?")
         Zwei Auflösungsstufen:
         - Populationsebene: pro Kammer x Zeitpunkt wird zunächst über alle
           Zellen gemittelt, DANN wird sigma/x_bar über die Zeitpunkte
           (innerhalb dieser Kammer) berechnet.
         - Einzelzell-Ebene: pro Zelle wird sigma/x_bar über ALLE Zeitpunkte
           berechnet, an denen genau diese Zelle beobachtet wurde.

  R(p) - Robustheit ÜBER DIE POPULATION ("wie homogen ist die Population
         zu einem gegebenen Zeitpunkt?")
         sigma/x_bar wird über ALLE Zellen einer Kammer zu EINEM Zeitpunkt
         berechnet - eine Zahl pro Kammer x Zeitpunkt.

In allen Fällen ist 'm' der Mittelwert über ALLE Zeitpunkte UND alle
Bedingungen, die in der jeweiligen Berechnung berücksichtigt werden sollen
(siehe `m_scope` Parameter unten) - NICHT pro Kammer einzeln.
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
import pandas as pd

from growth_rate import classify_condition_type

logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def _fano_robustness(sigma: pd.Series, x_bar: pd.Series, m: float) -> pd.Series:
    """Eq. 1: R = -sigma^2/x_bar * 1/m. x_bar nahe 0 -> NaN (keine sinnvolle Robustheit definierbar)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        r = -(sigma ** 2) / x_bar * (1.0 / m)
    r = r.where(x_bar.abs() > 1e-12, np.nan)
    return r


def compute_rt_population(
    cells: pd.DataFrame,
    value_col: str,
    group_cols: Optional[list[str]] = None,
) -> pd.DataFrame:
    """
    R(t) auf POPULATIONS-(Kammer-)Ebene.

    Schritt 1: pro Kammer x Zeitpunkt über alle Zellen mitteln (Populations-
               Mittelwert von value_col zu diesem Zeitpunkt).
    Schritt 2: pro Kammer über die Zeitpunkte hinweg sigma/x_bar berechnen.
    'm' ist der Mittelwert der Populations-Mittelwerte über ALLE Kammern und
    Zeitpunkte im übergebenen Datensatz (globaler Normalisierungsfaktor -
    ändert sich, wenn ihr den Datensatz erweitert, siehe Modul-Docstring).

    Parameters
    ----------
    cells : Datensatz mit exp_id, frame, value_col (eine Zeile pro Zelle x Frame)
    group_cols : zusätzliche Hierarchie-Spalten, die pro exp_id konstant
                 sind und in der Ausgabe erhalten bleiben sollen
                 (Standard: biosensor, osc_type, osc_freq, condition,
                 replicate, chamber - je nach Verfügbarkeit)

    Returns
    -------
    DataFrame mit einer Zeile pro exp_id: sigma, x_bar, m, R_t_population
    """
    if group_cols is None:
        group_cols = [c for c in ["biosensor", "osc_type", "osc_freq", "condition", "replicate", "chamber"]
                      if c in cells.columns]

    pop_mean_per_t = (
        cells.dropna(subset=[value_col])
        .groupby(["exp_id", "frame"])[value_col]
        .mean()
        .rename("pop_mean")
        .reset_index()
    )

    m = pop_mean_per_t["pop_mean"].mean()

    stats = (
        pop_mean_per_t.groupby("exp_id")["pop_mean"]
        .agg(sigma="std", x_bar="mean", n_timepoints="count")
        .reset_index()
    )
    stats["m"] = m
    stats["R_t_population"] = _fano_robustness(stats["sigma"], stats["x_bar"], m)

    if group_cols:
        meta = cells[["exp_id"] + group_cols].drop_duplicates("exp_id")
        stats = stats.merge(meta, on="exp_id", how="left")

    logger.info(
        "R(t) Populationsebene für '%s' berechnet: %d Kammern, m=%.4g.",
        value_col, len(stats), m,
    )
    return stats


def compute_rt_single_cell(
    cells: pd.DataFrame,
    value_col: str,
    group_cols: Optional[list[str]] = None,
    cell_id_col: str = "cell_uid",
) -> pd.DataFrame:
    """
    R(t) auf EINZELZELL-Ebene.

    Pro Zelle (cell_uid) wird sigma/x_bar von value_col über alle Zeitpunkte
    berechnet, an denen diese Zelle beobachtet wurde. 'm' ist der Mittelwert
    von value_col über ALLE Zellen und Zeitpunkte im übergebenen Datensatz.

    Returns
    -------
    DataFrame mit einer Zeile pro cell_uid: sigma, x_bar, m, R_t_single_cell, n_timepoints
    """
    if group_cols is None:
        group_cols = [c for c in ["biosensor", "osc_type", "osc_freq", "condition", "replicate", "chamber", "exp_id"]
                      if c in cells.columns]

    valid = cells.dropna(subset=[value_col])
    m = valid[value_col].mean()

    stats = (
        valid.groupby(cell_id_col)[value_col]
        .agg(sigma="std", x_bar="mean", n_timepoints="count")
        .reset_index()
    )
    stats["m"] = m
    stats["R_t_single_cell"] = _fano_robustness(stats["sigma"], stats["x_bar"], m)

    if group_cols:
        meta = cells[[cell_id_col] + group_cols].drop_duplicates(cell_id_col)
        stats = stats.merge(meta, on=cell_id_col, how="left")

    logger.info(
        "R(t) Einzelzell-Ebene für '%s' berechnet: %d Zellen, m=%.4g.",
        value_col, len(stats), m,
    )
    return stats


def compute_rp(
    cells: pd.DataFrame,
    value_col: str,
    group_cols: Optional[list[str]] = None,
) -> pd.DataFrame:
    """
    R(p) - Robustheit ÜBER DIE POPULATION zu jedem Zeitpunkt.

    Pro Kammer x Zeitpunkt wird sigma/x_bar von value_col ÜBER ALLE ZELLEN
    DIESER KAMMER ZU DIESEM ZEITPUNKT berechnet (Maß für Homogenität der
    Population). 'm' ist der Mittelwert von value_col über ALLE Zellen,
    Zeitpunkte und Kammern im übergebenen Datensatz.

    Returns
    -------
    DataFrame mit einer Zeile pro exp_id x frame: sigma, x_bar, m, R_p, n_cells
    """
    if group_cols is None:
        group_cols = [c for c in ["biosensor", "osc_type", "osc_freq", "condition", "replicate", "chamber"]
                      if c in cells.columns]

    valid = cells.dropna(subset=[value_col])
    m = valid[value_col].mean()

    stats = (
        valid.groupby(["exp_id", "frame"])[value_col]
        .agg(sigma="std", x_bar="mean", n_cells="count")
        .reset_index()
    )
    stats["m"] = m
    stats["R_p"] = _fano_robustness(stats["sigma"], stats["x_bar"], m)

    n_single_cell = (stats["n_cells"] <= 1).sum()
    if n_single_cell > 0:
        logger.warning(
            "%d von %d Kammer-Zeitpunkten haben <= 1 Zelle - sigma ist dort "
            "undefiniert (NaN), R_p kann für diese Zeitpunkte nicht berechnet werden.",
            n_single_cell, len(stats),
        )

    if group_cols:
        meta = cells[["exp_id"] + group_cols].drop_duplicates("exp_id")
        stats = stats.merge(meta, on="exp_id", how="left")

    logger.info(
        "R(p) für '%s' berechnet: %d Kammer-Zeitpunkte, m=%.4g.",
        value_col, len(stats), m,
    )
    return stats


def aggregate_robustness_over_replicates(
    robustness_df: pd.DataFrame,
    value_col: str,
    group_cols: Optional[list[str]] = None,
    replicate_id_col: str = "replicate",
    chamber_id_col: str = "exp_id",
) -> pd.DataFrame:
    """
    Mittelt einen Robustness-Wert (R_t_population, R_t_single_cell oder R_p)
    über biologische Replikate, analog zu den Punkt+Errorbar-Plots im Paper
    (Fig. 6, 7b: Mittelwert ± SD über Triplikate).

    DREISTUFIGE Aggregation:
      1) pro Kammer (chamber_id_col, Standard 'exp_id') mitteln.
      2) pro Replikat (replicate_id_col, Standard 'replicate') mitteln.
      3) über Replikate: Mittelwert + SD, n_replicates = echte Replikat-Anzahl.

    Schritt 1 ist für R_t_population/R_t_single_cell ein No-Op (compute_rt_*()
    liefert dort bereits genau eine Zeile pro Kammer). Für R_p dagegen ist er
    NOTWENDIG: compute_rp() liefert eine Zeile pro Kammer x Frame, d.h. ohne
    diesen Zwischenschritt würde jeder FRAME als eigenes "Replikat" gezählt -
    n_replicates wäre "Kammern x Frames", und die ausgewiesene SD würde Frame-
    und Kammer-Varianz vermischen.

    WARUM SCHRITT 2 DAZUGEKOMMEN IST
    --------------------------------
    Vorher war replicate_id_col='exp_id' der Default, und damit endete die
    Aggregation auf KAMMER-Ebene: 'exp_id' wird in data_loading.py aus
    biosensor + osc_type + osc_freq + condition + replicate + chamber gebaut,
    ist also eine Kammer-ID und keine Replikat-ID. Jede Kammer zählte danach
    als eigenes biologisches Replikat - bei 3 Replikaten mit je 2 Kammern
    stand in n_replicates eine 6, und die ausgewiesene SD mischte technische
    (Kammer-zu-Kammer) mit biologischer (Replikat-zu-Replikat) Varianz.
    Kammern eines Replikats sind technische Messungen desselben Chips, nicht
    unabhängige biologische Wiederholungen.

    Das ändert die Zahlen in allen 40_*_aggregated.csv: 'sd' und
    'n_replicates' fallen, und 'mean' verschiebt sich leicht, weil Kammern
    jetzt gleich gewichtet werden statt nach ihrer Zeilenzahl. Die Änderung
    ist gewollt - die Fehlerbalken beschreiben ab jetzt biologische Replikate.
    Wer die alte Kammer-Ebene braucht, ruft mit replicate_id_col='exp_id' auf.

    Parameters
    ----------
    robustness_df : Ergebnis von compute_rt_population/compute_rt_single_cell/compute_rp
    value_col : welche Spalte gemittelt wird (z.B. 'R_t_population')
    group_cols : Spalten, die EINE Bedingung definieren (ohne 'replicate'/'chamber').
                 Standard: biosensor, osc_type, osc_freq, condition - falls vorhanden.
    replicate_id_col : Spalte, die EIN biologisches Replikat identifiziert
                 (Standard 'replicate'). Fehlt sie, wird mit Warnung auf die
                 Kammer-Ebene zurückgefallen.
    chamber_id_col : Spalte, die EINE Kammer identifiziert (Standard 'exp_id').
    """
    if group_cols is None:
        group_cols = [c for c in ["biosensor", "osc_type", "osc_freq", "condition"]
                      if c in robustness_df.columns]

    df = robustness_df.dropna(subset=[value_col])

    # Stufe 1: pro Kammer - entfernt die Frame-Pseudoreplikation bei R_p.
    if chamber_id_col in df.columns:
        keys = group_cols + [c for c in (replicate_id_col, chamber_id_col) if c in df.columns]
        per_chamber = df.groupby(keys, dropna=False)[value_col].mean().reset_index()
    else:
        per_chamber = df

    # Stufe 2: pro Replikat - entfernt die Kammer-Pseudoreplikation.
    if replicate_id_col in per_chamber.columns:
        per_replicate = (
            per_chamber.groupby(group_cols + [replicate_id_col], dropna=False)[value_col]
            .mean()
            .reset_index()
        )
    else:
        logger.warning(
            "aggregate_robustness_over_replicates(): Spalte '%s' fehlt - es wird nur bis auf "
            "die Kammer-Ebene ('%s') aggregiert. 'n_replicates' zaehlt dann Kammern statt "
            "biologische Replikate, und 'sd' mischt technische mit biologischer Varianz.",
            replicate_id_col, chamber_id_col,
        )
        per_replicate = per_chamber

    agg = (
        per_replicate.groupby(group_cols)[value_col]
        .agg(mean="mean", sd="std", n_replicates="count")
        .reset_index()
    )

    # Siehe growth_rate.classify_condition_type(): dieselbe Kollision (mehrere
    # Zeilen pro osc_freq, weil PosCtrl/NegCtrl unter demselben osc_freq-Ordner
    # abgelegt sind) tritt hier ebenso auf wie bei summarise_growth_rate().
    if "condition" in agg.columns:
        agg["condition_type"] = agg["condition"].apply(classify_condition_type)

    return agg
