"""
control_consistency.py
=======================
Validierung der Kontrollbedingungen (PosCtrl/NegCtrl) ÜBER die verschiedenen
Oszillationsfrequenz-Batches hinweg.

HINTERGRUND
-----------
PosCtrl/NegCtrl-Kammern laufen in jedem osc_freq-Batch mit (siehe
growth_rate.classify_condition_type()) und sollten - wenn die
Mikrofluidik-Strömungsverhältnisse den Organismus nicht relevant
beeinflussen - über alle Frequenz-Batches hinweg VERGLEICHBARE Werte
liefern. Dieses Modul macht diesen Vergleich explizit, getrennt nach
PosCtrl/NegCtrl (und, falls mehrere Biosensoren vorhanden sind, auch
getrennt nach Biosensor):

  1. VISUELL: Punkt+Errorbar-Plot der Kontrollwerte über osc_freq
     (Wiederverwendung von summary_plots.plot_point_errorbar), Farbe =
     Kontroll-Kategorie (+ Biosensor, falls mehrere vorhanden).
  2. STATISTISCH: Kruskal-Wallis-Test pro (biosensor, osc_type,
     condition_type): unterscheiden sich die Werte signifikant zwischen den
     osc_freq-Batches? Eine SIGNIFIKANTE Differenz bei einer Kontrolle ist
     ein Hinweis auf einen Strömungs-/Batch-Effekt, der unabhängig von der
     eigentlichen Oszillationsbedingung ist.

WICHTIG: Kruskal-Wallis testet nur "unterscheidet sich MINDESTENS eine
Gruppe von den anderen", nicht WELCHE - für Post-hoc-Paarvergleiche wäre
z.B. Dunn's Test nötig. Das ist hier bewusst nicht eingebaut, um die
explorative Fragestellung ("gibt es überhaupt einen Batch-Effekt?") nicht
unnötig zu verkomplizieren.

WICHTIG zu den Eingabedaten: die statistischen Tests brauchen die Werte AUF
REPLIKAT-/EINZELZELL-EBENE (z.B. mu_table aus compute_specific_growth_rate(),
oder rt_pop/rp aus robustness.py) - NICHT die bereits über Replikate
gemittelten Summary-Tabellen, sonst hätte jeder osc_freq-Batch nur einen
einzigen Wert und Kruskal-Wallis wäre nicht durchführbar. Für den visuellen
Plot ist dagegen die aggregierte Summary-Tabelle (mit mean/sd-Spalten)
richtig - genau wie bei plot_point_errorbar() sonst auch.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import pandas as pd
from scipy import stats

from growth_rate import classify_condition_type
from summary_plots import plot_point_errorbar

logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def add_condition_type(df: pd.DataFrame) -> pd.DataFrame:
    """Ergänzt 'condition_type' (Oscillation/PosCtrl/NegCtrl/Unknown), falls noch nicht vorhanden."""
    if "condition_type" in df.columns:
        return df
    if "condition" not in df.columns:
        raise ValueError("add_condition_type() braucht eine 'condition' Spalte.")
    df = df.copy()
    df["condition_type"] = df["condition"].apply(classify_condition_type)
    return df


def filter_controls(df: pd.DataFrame) -> pd.DataFrame:
    """Reduziert auf PosCtrl/NegCtrl-Zeilen (schließt 'Oscillation' und 'Unknown' aus)."""
    df = add_condition_type(df)
    return df[df["condition_type"].isin(["PosCtrl", "NegCtrl"])].copy()


def test_control_consistency_across_freq(
    df: pd.DataFrame,
    value_col: str,
    group_cols: Optional[list[str]] = None,
    freq_col: str = "osc_freq",
    min_n_per_group: int = 2,
) -> pd.DataFrame:
    """
    Kruskal-Wallis-Test pro (group_cols)-Kombination: unterscheiden sich die
    value_col-Werte EINER Kontrolle (PosCtrl oder NegCtrl) signifikant
    ZWISCHEN den osc_freq-Batches, in denen sie mitgelaufen ist?

    Parameters
    ----------
    df : Rohdaten AUF REPLIKAT-/EINZELZELL-EBENE (siehe Modul-Docstring).
         Muss 'condition' (oder bereits 'condition_type') und freq_col
         enthalten. Wird intern auf PosCtrl/NegCtrl gefiltert.
    value_col : zu testende numerische Spalte (z.B. 'mu', 'R_t_population', 'R_p')
    group_cols : definieren EINE "Kontroll-Bedingung", innerhalb derer über
                 die osc_freq-Batches getestet wird.
                 Standard: biosensor, osc_type, condition_type - je nach Verfügbarkeit.
    freq_col : Spalte mit dem Frequenz-Batch (Standard 'osc_freq')
    min_n_per_group : osc_freq-Batches mit weniger als so vielen Werten
                 werden aus dem Test ausgeschlossen (Kruskal-Wallis braucht
                 mind. 1 Wert/Gruppe, ist aber mit sehr kleinen Gruppen
                 statistisch kaum aussagekräftig - Standard 2).

    Returns
    -------
    DataFrame mit einer Zeile pro group_cols-Kombination:
        group_cols..., n_freq_batches, n_total, h_statistic, p_value
    p_value = NaN, falls weniger als 2 osc_freq-Batches mit ausreichend
    Werten übrig bleiben (Test nicht durchführbar - explizit markiert, nicht
    stillschweigend übersprungen).
    """
    df = filter_controls(df)
    df = df.dropna(subset=[value_col, freq_col])

    if df.empty:
        logger.warning(
            "test_control_consistency_across_freq(): keine PosCtrl/NegCtrl-Zeilen mit "
            "gültigem '%s' gefunden - Test übersprungen.", value_col,
        )
        return pd.DataFrame()

    if group_cols is None:
        group_cols = [c for c in ["biosensor", "osc_type", "condition_type"] if c in df.columns]

    results = []
    for keys, grp in df.groupby(group_cols):
        keys = keys if isinstance(keys, tuple) else (keys,)
        freq_groups = [
            g[value_col].to_numpy()
            for _, g in grp.groupby(freq_col)
            if len(g) >= min_n_per_group
        ]
        n_freq_batches = len(freq_groups)

        record = dict(zip(group_cols, keys))
        record["n_freq_batches"] = n_freq_batches
        record["n_total"] = sum(len(g) for g in freq_groups)

        if n_freq_batches < 2:
            record["h_statistic"] = float("nan")
            record["p_value"] = float("nan")
            logger.warning(
                "Kruskal-Wallis für %s: nur %d osc_freq-Batch(es) mit >= %d Werten - "
                "Test nicht durchführbar (braucht mindestens 2).",
                dict(zip(group_cols, keys)), n_freq_batches, min_n_per_group,
            )
        else:
            h, p = stats.kruskal(*freq_groups)
            record["h_statistic"] = h
            record["p_value"] = p

        results.append(record)

    out = pd.DataFrame(results)
    if not out.empty:
        n_significant = int((out["p_value"] < 0.05).sum())
        logger.info(
            "Kontroll-Konsistenz-Test ('%s' über osc_freq-Batches): %d Gruppen getestet, "
            "%d mit p < 0.05 (signifikanter Unterschied zwischen Frequenz-Batches - "
            "möglicher Strömungs-/Batch-Effekt).",
            value_col, len(out), n_significant,
        )
    return out


def plot_control_consistency(
    summary: pd.DataFrame,
    value_col: str,
    out_path: Path,
    facet_col: str = "osc_type",
    x_col: str = "osc_freq",
    x_order: Optional[list] = None,
    ylabel: Optional[str] = None,
    title: Optional[str] = None,
) -> None:
    """
    Punkt+Errorbar-Plot NUR der Kontrollbedingungen (PosCtrl/NegCtrl) über
    osc_freq - visuelles Gegenstück zu test_control_consistency_across_freq().

    Erwartet die bereits über Replikate AGGREGIERTE Summary-Tabelle (z.B.
    mu_summary aus growth_rate.summarise_growth_rate(), oder rt_pop_agg/
    rp_agg aus robustness.aggregate_robustness_over_replicates()) - siehe
    Modul-Docstring zum Unterschied zur statistischen Testfunktion.

    Farbe = Biosensor (konsistent mit den übrigen Plots der Pipeline, z.B.
    09_specific_growth_rate.pdf), Symbol-Form = PosCtrl/NegCtrl - so bleibt
    "welcher Biosensor" und "welche Kontrollart" auf einen Blick trennbar,
    statt beides in eine kombinierte Farbe zu packen.
    """
    controls = filter_controls(summary)
    if controls.empty:
        logger.warning("plot_control_consistency(): keine Kontroll-Zeilen (PosCtrl/NegCtrl) gefunden - Plot übersprungen.")
        return

    plot_point_errorbar(
        controls, value_col=value_col, out_path=out_path,
        x_col=x_col, facet_col=facet_col,
        color_col="biosensor" if "biosensor" in controls.columns else None,
        style_col="condition_type",
        x_order=x_order, ylabel=ylabel,
        title=title or "Control consistency across oscillation-frequency batches\n(flat = control robust, trend = possible flow/batch effect)",
    )
