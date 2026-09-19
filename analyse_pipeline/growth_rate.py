"""
growth_rate.py
==============
Spezifische Wachstumsrate µ pro Zelle, nach dem Prinzip von Eq. 2 aus
Blöbaum et al. 2024 (Microbial Cell Factories), VERALLGEMEINERT auf
Reproduktionsereignisse mit MEHREREN gleichzeitigen Tochterzellen
(z.B. mehrere Blastokonidien pro Knospungs-Ereignis):

    µ = ln((n + k) / n) / t

wobei:
    n = Zellzahl unmittelbar VOR dem Ereignis (bei Eq. 2 implizit n=1,
        da das Paper von klassischem Einzel-Budding ausgeht)
    k = Anzahl Tochterzellen, die GLEICHZEITIG in diesem Ereignis entstehen
        (k=1 im klassischen Fall -> reduziert sich exakt auf Eq. 2: ln(2)/t)
    t = Zeit (h) zwischen diesem und dem VORHERIGEN Reproduktionsereignis
        derselben Mutterzelle

HERLEITUNG: Eq. 2 selbst stammt aus dem exponentiellen Wachstumsgesetz
N(t) = N0 * 2^(t/tau) für eine Verdopplung (Mutter -> Mutter + 1 Bud, d.h.
Faktor 2 pro Ereignis). Entstehen pro Ereignis k Töchter statt 1, wächst die
lokale Zellzahl um den Faktor (n+k)/n statt um den Faktor 2 - die Formel
wird entsprechend verallgemeinert, bleibt aber für k=1 identisch zu Eq. 2.

GRUPPIERUNG GLEICHZEITIGER BUDS: Alle Budding-Events derselben Mutterzelle
mit IDENTISCHEM budding_frame werden zu EINEM Ereignis mit k=Anzahl dieser
Events zusammengefasst, BEVOR die Zeitintervalle zwischen Ereignissen
gebildet werden. Das ist der biologisch korrekte Umgang mit Organismen, die
routinemäßig mehrere Tochterzellen pro Zyklus bilden (z.B. multipolares
Budding) - im Unterschied zu klassischem Hefe-Budding, wo Gleichzeitigkeit
eher ein Tracking-Artefakt wäre.

Produziert eine Mutter z.B. 3 Ereignisse (unabhängig davon wie viele
Töchter jedes Ereignis hatte), ergeben sich 2 Zeitintervalle und damit
2 µ-Werte für diese Zelle (NICHT 3 - das erste Intervall würde von
Kammerstart bis zum ersten Ereignis reichen und enthält die unbekannte
Vorgeschichte der Zelle vor Beobachtungsbeginn; das ist nicht mit den
späteren, vollständig beobachteten Intervallen vergleichbar und wird daher
bewusst ausgeschlossen).

ARTEFAKT-FILTER (Paper, Methods): µ-Werte über einer Schwelle (Standard:
0.6 h^-1) werden markiert (nicht gelöscht). Im Paper wird das Beispiel
genannt, dass eine Zelle, die 30 min nach Inokulation der Kammer knospt,
einen unplausiblen Wert von µ=1.3 h^-1 ergäbe - das ist aber genau der Fall
des ERSTEN Intervalls (Kammerstart bis 1. Ereignis), den wir hier bereits
strukturell ausschließen. Der Schwellenwert-Filter bleibt trotzdem als
zusätzliche Sicherung bestehen, für den Fall von Tracking-Artefakten.

ANNAHME ZUR ZELLZAHL "n" VOR DEM EREIGNIS: dieses Modul setzt n=1 für JEDES
Ereignis (d.h. es betrachtet die Mutterzelle isoliert, nicht die gesamte
Kammer-Population). Das entspricht der Paper-Logik (Eq. 2 ist explizit
EINZELZELL-bezogen, im Unterschied zur Budding Ratio in Eq. 3, die auf
Kammer-Ebene rechnet). Falls eure Mutterzelle zum Zeitpunkt des Ereignisses
selbst aus einem vorherigen Multi-Bud-Ereignis stammt, wird das NICHT
rückwirkend in n eingerechnet - n bleibt 1 (die eine Mutter), nur k variiert
je nach Anzahl ihrer eigenen, gleichzeitig entstandenen Töchter.

UNABHÄNGIG VON DER BEDINGUNG: µ wird für jede Zelle berechnet, die
mindestens 2 Reproduktionsereignisse hat - unabhängig davon, ob die Zelle
aus einer Oszillationsbedingung, der PosCtrl (durchgehend Feast) oder der
NegCtrl (durchgehend Starvation, siehe Hinweis in budding_ratio_timeseries.py)
stammt. PosCtrl/NegCtrl sind hier einfach zwei weitere Werte der
'condition'-Spalte, für die dieselbe Formel gilt.
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def classify_condition_type(condition: str) -> str:
    """
    Ordnet einen rohen 'condition'-String (z.B. 'Osc0.75_NegCtrl') seiner
    Kontroll-Kategorie zu. Folgt dem festen Namensmuster aus
    extract_metadata_from_filename() (cellpose_pipeline_v9.py): Kontrollkammern
    enden auf '_NegCtrl'/'_PosCtrl', alles andere ist die eigentliche
    Oszillationsbedingung.

    Wird u.a. von plot_point_errorbar() (summary_plots.py) als style_col
    genutzt, um Oscillation/PosCtrl/NegCtrl-Zeilen sichtbar zu trennen statt
    sie beim Reindex auf einen gemeinsamen osc_freq-Wert kollidieren zu lassen.
    """
    if pd.isna(condition):
        return "Unknown"
    if condition.endswith("NegCtrl"):
        return "NegCtrl"
    if condition.endswith("PosCtrl"):
        return "PosCtrl"
    if condition.startswith("St."):
        # Statische (nicht-oszillierende) Bedingungen, z.B. 'St.omlp'/'St.ypd'
        # (siehe Data/<Biosensor>/static/static_<medium>/...). Diese sollen
        # als eigene Kategorie sichtbar bleiben statt fälschlich unter
        # 'Oscillation' zu laufen - für sie gibt es ja keine Oszillation.
        return condition
    return "Oscillation"


def compute_specific_growth_rate(
    lineage_events: pd.DataFrame,
    cells: pd.DataFrame,
    min_per_frame: float,
    mu_max_threshold: float = 0.6,
) -> pd.DataFrame:
    """
    Berechnet µ für jede Mutterzelle mit >= 2 Reproduktionsereignissen.
    Gleichzeitige Buds (identischer budding_frame) werden vorher zu einem
    Ereignis mit Töchterzahl k zusammengefasst (siehe Modul-Docstring).

    Parameters
    ----------
    lineage_events : Ergebnis von lineage.classify_mother_bud() - braucht
                     Spalten mother_cell_uid, budding_frame
    cells : voller Zelldatensatz (für Hierarchie-Metadaten)
    min_per_frame : Minuten pro Frame (Zeitkalibrierung, z.B. MIN_PER_FRAME
                    aus run_analysis.py) - zur Umrechnung von Frame-Differenz in Stunden
    mu_max_threshold : µ-Werte darüber werden als Artefakt markiert
                    (Standard 0.6 h^-1, siehe Paper Methods: maximales
                    realistisches µ für Hefe liegt bei ca. 0.5 h^-1; passt
                    diesen Wert ggf. an euren Organismus an)

    Returns
    -------
    DataFrame mit einer Zeile pro µ-Wert (eine Mutterzelle kann mehrere
    Zeilen haben, wenn sie mehrere Ereignis-Intervalle hat):
        exp_id, mother_cell_uid, interval_index, frame_start, frame_end,
        delta_t_h, k_daughters, mu, mu_is_artefact (+ Hierarchie-Metadaten)

    k_daughters bezieht sich auf das Ereignis am ENDE des Intervalls
    (frame_end) - also wie viele Töchter zu DIESEM Zeitpunkt gleichzeitig
    entstanden sind.

    Zeilen mit mu_is_artefact=True sind in der Tabelle ENTHALTEN (nicht
    stillschweigend gelöscht), aber klar markiert - so bleibt nachvollziehbar,
    wie viele Werte gefiltert wurden und warum. Für die übliche Auswertung
    filtert man einfach auf mu_is_artefact == False.
    """
    required = {"mother_cell_uid", "budding_frame", "exp_id"}
    missing = required - set(lineage_events.columns)
    if missing:
        raise ValueError(
            f"compute_specific_growth_rate() fehlen Spalten in 'lineage_events': {missing}. "
            f"Stellt sicher, dass classify_mother_bud() mit einem Datensatz aufgerufen wurde, "
            f"der eine 'cell_uid' Spalte enthält (siehe lineage.py)."
        )

    if lineage_events.empty:
        logger.warning("lineage_events ist leer - compute_specific_growth_rate() gibt eine leere Tabelle zurück.")
        return pd.DataFrame()

    results = []
    n_single_event_mothers = 0
    n_multi_bud_events_total = 0

    for mother_uid, grp in lineage_events.groupby("mother_cell_uid"):
        exp_id = grp["exp_id"].iloc[0]

        # Schritt 1: gleichzeitige Buds (identischer budding_frame) zu EINEM
        # Ereignis mit Töchterzahl k zusammenfassen.
        events_by_frame = (
            grp.groupby("budding_frame")
            .size()
            .rename("k_daughters")
            .reset_index()
            .sort_values("budding_frame")
        )
        n_multi_bud_events_total += int((events_by_frame["k_daughters"] > 1).sum())

        frames = events_by_frame["budding_frame"].to_numpy()
        k_values = events_by_frame["k_daughters"].to_numpy()

        if len(frames) < 2:
            n_single_event_mothers += 1
            continue

        # Schritt 2: µ aus den Intervallen ZWISCHEN aufeinanderfolgenden
        # Ereignissen berechnen, verallgemeinerte Formel ln((n+k)/n)/t mit n=1.
        for i in range(1, len(frames)):
            frame_start, frame_end = frames[i - 1], frames[i]
            k = int(k_values[i])  # Töchterzahl des Ereignisses AM ENDE des Intervalls
            delta_t_h = (frame_end - frame_start) * min_per_frame / 60.0

            if delta_t_h <= 0:
                # Sollte nach der Gruppierung in Schritt 1 nicht mehr vorkommen
                # (gleichzeitige Buds sind bereits zusammengefasst), außer bei
                # widersprüchlichen/doppelten Frame-Einträgen in den Rohdaten.
                logger.warning(
                    "Mutter '%s': unerwartetes delta_t<=0 zwischen Frame %d und %d "
                    "nach Gruppierung gleichzeitiger Buds - wird übersprungen. "
                    "Mögliche Ursache: doppelte/widersprüchliche Einträge in lineage_events.",
                    mother_uid, frame_start, frame_end,
                )
                continue

            n = 1  # Mutterzelle isoliert betrachtet, siehe Modul-Docstring
            mu = np.log((n + k) / n) / delta_t_h
            is_artefact = mu > mu_max_threshold

            results.append({
                "exp_id": exp_id,
                "mother_cell_uid": mother_uid,
                "interval_index": i,
                "frame_start": frame_start,
                "frame_end": frame_end,
                "delta_t_h": delta_t_h,
                "k_daughters": k,
                "mu": mu,
                "mu_is_artefact": is_artefact,
            })

    out = pd.DataFrame(results)

    if n_single_event_mothers > 0:
        logger.info(
            "%d Mutterzellen hatten nur 1 Reproduktionsereignis (kein Intervall berechenbar) - "
            "übersprungen, nicht als Artefakt gezählt.",
            n_single_event_mothers,
        )
    if n_multi_bud_events_total > 0:
        logger.info(
            "%d Ereignisse mit mehreren gleichzeitigen Buds (k > 1) erkannt und als "
            "EIN Ereignis mit entsprechender Töchterzahl in µ eingerechnet.",
            n_multi_bud_events_total,
        )

    if out.empty:
        logger.warning("Keine µ-Werte berechnet (keine Mutter mit >= 2 Reproduktionsereignissen).")
        return out

    n_artefacts = int(out["mu_is_artefact"].sum())
    if n_artefacts > 0:
        logger.warning(
            "%d von %d µ-Werten liegen über mu_max_threshold=%.2f h^-1 und sind als "
            "Artefakt markiert (mu_is_artefact=True). Für die Auswertung auf "
            "mu_is_artefact == False filtern.",
            n_artefacts, len(out), mu_max_threshold,
        )

    meta_cols = [c for c in ["biosensor", "osc_type", "osc_freq", "condition", "replicate", "chamber",
                             "chip", "chip_family", "medium", "date"]
                 if c in cells.columns]
    if meta_cols:
        meta = cells[["exp_id"] + meta_cols].drop_duplicates("exp_id")
        out = out.merge(meta, on="exp_id", how="left")

    logger.info(
        "Spezifische Wachstumsrate berechnet: %d µ-Werte aus %d Mutterzellen "
        "(%d Artefakte, %d gültig).",
        len(out), out["mother_cell_uid"].nunique(), n_artefacts, len(out) - n_artefacts,
    )

    return out


def summarise_growth_rate(
    mu_table: pd.DataFrame,
    group_cols: Optional[list[str]] = None,
    exclude_artefacts: bool = True,
) -> pd.DataFrame:
    """
    Aggregiert µ-Werte pro Bedingung (z.B. für Punkt+Errorbar-Plots wie
    Paper Fig. 3b, "Specific Growth Rate" Panel).

    Parameters
    ----------
    mu_table : Ergebnis von compute_specific_growth_rate()
    group_cols : Spalten, die eine Bedingung definieren. Standard: biosensor,
                 osc_type, osc_freq, condition (ohne replicate/chamber, damit
                 über Replikate gemittelt wird - siehe exclude_artefacts).
    exclude_artefacts : ob mu_is_artefact==True Zeilen vor der Aggregation
                 ausgeschlossen werden (Standard True - empfohlen für jede
                 normale Auswertung).
    """
    if mu_table.empty:
        return pd.DataFrame()

    df = mu_table
    if exclude_artefacts:
        df = df[~df["mu_is_artefact"]]

    if group_cols is None:
        group_cols = [c for c in ["biosensor", "osc_type", "osc_freq", "condition",
                                   "medium", "chip_family"] if c in df.columns]

    summary = (
        df.groupby(group_cols)
        .agg(
            mean_mu=("mu", "mean"),
            sd_mu=("mu", "std"),
            n_values=("mu", "count"),
            mean_k_daughters=("k_daughters", "mean"),
        )
        .reset_index()
    )

    # Grobe Kontroll-Kategorie (Oscillation/PosCtrl/NegCtrl) als eigene Spalte,
    # damit Plots (siehe summary_plots.plot_point_errorbar, style_col) die
    # Kontrollen sichtbar von der eigentlichen Oszillationsbedingung trennen
    # können, statt beim Reindex auf denselben osc_freq-Wert zu kollidieren.
    if "condition" in summary.columns:
        summary["condition_type"] = summary["condition"].apply(classify_condition_type)

    return summary
