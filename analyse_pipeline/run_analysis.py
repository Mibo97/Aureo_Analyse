"""
run_analysis.py
===============
Hauptskript der Analyse-Pipeline. Ausführen mit:

    python run_analysis.py

Erwartete Ordnerstruktur unter DATA_ROOT:

    Data/<Biosensor>/<Oszillationstyp>/<Oszillationsfrequenz>/03_results/Combined_Results.{csv,parquet}

z.B.:
    Data/iGlucoSnFR/Glucose/5min/03_results/Combined_Results.csv
    Data/pHRed/pH/10min/03_results/Combined_Results.csv

Ablauf (Kernreihe 00->50, danach Anhang 90/95):
    1. Alle Combined_Results einlesen (gecacht als .parquet für schnelles Neuladen)
    2. QC-Exclusion-Tabelle anwenden (nicht-destruktiv, siehe qc_exclusions.py)
    3. Mutter/Bud-Klassifikation vorab berechnen (Datenabhängigkeit: wird von
       Sektion 10 UND Sektion 20 benötigt, siehe Kommentar dort)
    00. Übersicht / Sanity-Check
    10. Zellmorphologie & Wachstum (Fläche, µ_event, µ_area)
    20. Lineage/Budding (Klassifikation-Output, Budding Ratio, Panel A, Lineage-Tiefe)
    30. Sensor-/Ratio-Zeitverläufe (unterstützend, Grundlage für Robustness)
    40. Robustness R(t)/R(p) - Synthese, inkl. Sensor-Ratio-Spalten,
        µ_event (R(t) Einzelzelle) und µ_area (R(p)); siehe Kommentare in
        dieser Sektion, warum jeweils nur EINE der beiden Varianten
        sinnvoll berechenbar ist
    50. Zusammenfassungstabelle
    90. Anhang: Morphologie-Scatter, Einzelzell-Trajektorien, stabile Mutter-Trajektorien
    95. Anhang: Sensor-QC/Controls (PosCtrl-vs-NegCtrl-Validierung pro Biosensor)
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from data_loading import load_all_results
from qc_exclusions import (
    init_qc_exclusions,
    read_qc_exclusions,
    summarise_qc_exclusions,
    summarise_track_merges,
    apply_track_merges,
    apply_qc_exclusions,
)
from sensors import compute_ratios, SENSOR_CONFIG
from lineage import classify_mother_bud, compute_budding_ratio, identify_mothers
from mother_trajectories import plot_stable_mother_per_group, build_lineage_tree, summarise_lineage_depth
from budding_ratio_timeseries import compute_budding_ratio_timeseries, aggregate_budding_ratio_over_replicates
from growth_rate import compute_specific_growth_rate, summarise_growth_rate, classify_condition_type
from area_growth import compute_area_growth_rate, summarise_area_growth, plot_mu_event_vs_mu_area
from robustness import (
    compute_rt_population,
    compute_rt_single_cell,
    compute_rp,
    aggregate_robustness_over_replicates,
)
from control_consistency import test_control_consistency_across_freq, plot_control_consistency
from queen_controls import (
    prepare_sensor_controls,
    summarise_sensor_controls,
    plot_sensor_control_timeseries,
    plot_sensor_control_comparison,
    plot_sensor_raw_channel_timeseries,
    summarise_sensor_controls_by_chamber,
    plot_control_chamber_comparison,
)
from violin_plots import plot_panel_a
from summary_plots import (
    plot_budding_ratio_timeseries,
    plot_point_errorbar,
    plot_rt_vs_rp_quadrant,
    plot_rt_single_cell_distribution,
)
from analysis import (
    pretty_label,
    add_time_column,
    find_intensity_columns,
    data_overview,
    plot_n_tracks_overview,
    plot_metric_over_time_by_frequency,
    plot_morphology_scatter,
    plot_single_cell_trajectories,
    summary_statistics,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def exclude_controls(df: pd.DataFrame) -> pd.DataFrame:
    """
    Entfernt PosCtrl-/NegCtrl-Zeilen aus einer Tabelle für die 'normalen'
    Plots - die Kontrollen selbst werden bereits eigenständig über
    test_control_consistency_across_freq()/plot_control_consistency()
    ausgewertet (Suffix '_control_consistency_*'), daher sollen sie in den
    übrigen Plots nicht mehr doppelt auftauchen.

    Nutzt bevorzugt eine vorhandene 'condition_type'-Spalte (z.B. aus
    summarise_growth_rate() oder den Robustness-Aggregaten); ist die nicht
    vorhanden, wird sie aus 'condition' on-the-fly berechnet (siehe
    growth_rate.classify_condition_type()). Fehlen beide Spalten, wird die
    Tabelle unverändert zurückgegeben (mit Warnung), statt fälschlich alles
    zu verwerfen. Wirkt NUR auf das, was geplottet wird - CSV-Exporte und
    die Kontroll-Konsistenz-Berechnungen selbst bleiben unverändert und
    behalten weiterhin PosCtrl/NegCtrl.
    """
    if df is None or df.empty:
        return df
    if "condition_type" in df.columns:
        ctype = df["condition_type"]
    elif "condition" in df.columns:
        ctype = df["condition"].apply(classify_condition_type)
    else:
        logger.warning(
            "exclude_controls(): weder 'condition' noch 'condition_type' in der Tabelle "
            "vorhanden - wird ungefiltert (inkl. evtl. Kontrollen) zurückgegeben."
        )
        return df
    return df[~ctype.isin(["PosCtrl", "NegCtrl"])].copy()

# ==============================================================================
# 0. KONFIGURATION
# ==============================================================================
# Alle Pfade und Parameter stehen in config.py - dort (oder über die
# Umgebungsvariablen AUREO_DATA_ROOT / AUREO_OUTPUT_DIR) anpassen, NICHT hier.
# inspect_lineage.py liest dieselbe Datei, damit Kalibrierung und Auswertung
# garantiert mit denselben Schwellen laufen.
from config import (
    DATA_ROOT,
    OUTPUT_DIR,
    OUTPUT_DIR_STATIC,
    CACHE_PATH,
    QC_EXCLUSIONS_PATH,
    FORCE_RELOAD,
    MIN_PER_FRAME,
    FREQ_ORDER,
    STATIC_ORDER,
    OSCILLATION_START_MIN,
    PANEL_A_GROUP_COL,
    PANEL_A_FACET_COL,
    CONTROL_CONCENTRATION_LABELS,
    LINEAGE_PARAMS,
    STABLE_MOTHER_MIN_COVERAGE,
    STABLE_MOTHER_GROUP_COLS,
    STABLE_MOTHER_BASE_VALUE_COLS,
    FLUX_CONFIG,
    MU_MAX_THRESHOLD,
    ROBUSTNESS_VALUE_COLS,
    log_active_configuration,
)


def resolve_x_order(cells: pd.DataFrame, col: str, preferred_order: list[str] | None) -> list[str]:
    """
    Bestimmt die x-Achsen-Reihenfolge für einen bestimmten Datensatz robust
    aus dem, was TATSÄCHLICH in 'col' vorkommt - statt blind eine hart
    codierte Liste (z.B. FREQ_ORDER/STATIC_ORDER) an x_order durchzureichen.

    Hintergrund: plot_point_errorbar() nutzt x_order typischerweise, um die
    Spalte in eine kategoriale Reihenfolge zu bringen - Werte, die NICHT in
    x_order stehen, fallen dabei meist als NaN raus und verschwinden
    kommentarlos aus dem Plot (z.B. weil der reale Ordnername anders
    geschrieben ist als in FREQ_ORDER/STATIC_ORDER angenommen). Damit das
    nie mehr unbemerkt passiert, werden hier nur die Werte aus
    preferred_order übernommen, die auch WIRKLICH in den Daten vorkommen,
    und alles Weitere, was in den Daten ist aber nicht in preferred_order
    steht, hinten angehängt (statt verworfen).
    """
    present = sorted(cells[col].dropna().unique().tolist(), key=str)
    if not preferred_order:
        return present
    ordered = [v for v in preferred_order if v in present]
    extra = [v for v in present if v not in preferred_order]
    if extra:
        logger.warning(
            "resolve_x_order(): folgende Werte in '%s' stehen nicht in der konfigurierten "
            "Reihenfolge %s und werden hinten angehängt statt verworfen: %s. "
            "Falls das nicht erwartet ist (z.B. Tippfehler/Groß-Kleinschreibung), "
            "bitte FREQ_ORDER/STATIC_ORDER anpassen.",
            col, preferred_order, extra,
        )
    return ordered + extra


def main() -> None:
    # Aktive Konfiguration inkl. aller Abweichungen vom dokumentierten Standard
    # ganz am Anfang ins Log - siehe config.log_active_configuration().
    log_active_configuration()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR_STATIC.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Daten laden
    # ------------------------------------------------------------------
    cells = load_all_results(DATA_ROOT, cache_path=CACHE_PATH, force_reload=FORCE_RELOAD)

    # ------------------------------------------------------------------
    # 1b. Ratiometrische Biosensoren: Verhältnis-Spalten berechnen
    # ------------------------------------------------------------------
    cells = compute_ratios(cells)  # nutzt SENSOR_CONFIG aus sensors.py
    if not SENSOR_CONFIG:
        logger.info(
            "Hinweis: SENSOR_CONFIG in sensors.py ist noch leer. Falls ihr "
            "ratiometrische Biosensoren habt, dort die Kanal-Paare eintragen."
        )

    # Spaltennamen stehen jetzt fest (QC entfernt nur Zeilen, keine Spalten) -
    # daher hier schon einmal bestimmen, statt wie zuvor erst mitten in der
    # Sensor-Plot-Sektion. Wird von Sektion 30 (Plots), 40 (Robustness) und
    # 50/90 (Summary/Supplement) weiterverwendet.
    intensity_cols = find_intensity_columns(cells)
    ratio_cols = [c for c in cells.columns if c.startswith("ratio_")]

    # ------------------------------------------------------------------
    # 2. QC: Track-Merges DANN Exclusions anwenden (nicht-destruktiv)
    # ------------------------------------------------------------------
    # WICHTIG: Track-Merges müssen VOR den Exclusions UND vor jeglicher
    # Lineage-Klassifikation laufen - sonst sieht classify_mother_bud() weiterhin
    # zwei getrennte Tracks und ein Tracking-Bruch (z.B. durch starke
    # Zellbewegung) wird fälschlich als neuer Bud interpretiert.
    init_qc_exclusions(QC_EXCLUSIONS_PATH)  # legt leere Datei an, falls noch keine existiert
    exclusions = read_qc_exclusions(QC_EXCLUSIONS_PATH)

    merge_summary = summarise_track_merges(exclusions)
    if not merge_summary.empty:
        logger.info("Konfigurierte Track-Merges:\n%s", merge_summary.to_string(index=False))
    cells = apply_track_merges(cells, exclusions)

    qc_summary = summarise_qc_exclusions(exclusions)
    if not qc_summary.empty:
        logger.info("QC-Exclusions pro Experiment:\n%s", qc_summary.to_string(index=False))

    cells = apply_qc_exclusions(cells, exclusions, mode="remove")
    cells = add_time_column(cells, MIN_PER_FRAME)

    logger.info("=== Daten nach QC (gesamt, vor Trennung Oszillation/statisch) ===")
    logger.info("  Zeilen total:                    %d", len(cells))
    logger.info("  Eindeutige Tracks (cell_uid):    %d", cells["cell_uid"].nunique())
    logger.info("  Eindeutige Experimente (exp_id): %d", cells["exp_id"].nunique())
    logger.info("  Biosensoren:            %s", sorted(cells["biosensor"].dropna().unique()))
    logger.info("  Oszillationstypen:      %s", sorted(cells["osc_type"].dropna().unique()))
    logger.info("  Oszillationsfrequenzen: %s", sorted(cells["osc_freq"].dropna().unique()))

    # ------------------------------------------------------------------
    # 2b. Statische Daten (Data/<Biosensor>/static/static_<medium>/...,
    #     condition St.omlp/St.ypd) von den Oszillationsdaten trennen - ab
    #     hier laufen beide komplett getrennt durch dieselbe Auswertung
    #     (run_pipeline()), landen aber in eigenen Output-Ordnern/Dateien.
    # ------------------------------------------------------------------
    # .str.strip().str.lower() statt striktem "==" Vergleich: robust gegen
    # Groß-/Kleinschreibung ("Static"/"STATIC") oder versehentliche
    # Leerzeichen im Ordnernamen - sonst landen statische Zeilen fälschlich
    # (und unbemerkt) in cells_osc statt in cells_static, siehe die Warnung
    # 'static_omlp'/'static_ypd' stehen nicht in FREQ_ORDER.
    is_static = cells["osc_type"].astype(str).str.strip().str.lower() == "static"
    cells_osc = cells[~is_static].copy()
    cells_static = cells[is_static].copy()

    logger.info(
        "Aufgeteilt: %d Zeilen Oszillationsdaten, %d Zeilen statische Daten (osc_type == 'static').",
        len(cells_osc), len(cells_static),
    )

    run_pipeline(
        cells_osc, output_dir=OUTPUT_DIR, freq_order=resolve_x_order(cells_osc, "osc_freq", FREQ_ORDER),
        intensity_cols=intensity_cols, ratio_cols=ratio_cols,
        run_sensor_controls=True,
    )

    if cells_static.empty:
        logger.warning("Keine statischen Daten (osc_type == 'static') gefunden - Schritt für statische Daten übersprungen.")
    else:
        # Keine Sensor-/Ratio-Auswertung für statische Daten: hier wird nur
        # der Wildtyp kultiviert, es gibt also keine Fluoreszenzkanäle
        # (siehe Rücksprache) - daher intensity_cols/ratio_cols leer und
        # run_sensor_controls=False, statt auf durchgehend NaN-Spalten zu
        # rechnen.
        run_pipeline(
            cells_static, output_dir=OUTPUT_DIR_STATIC,
            freq_order=resolve_x_order(cells_static, "osc_freq", STATIC_ORDER),
            intensity_cols=[], ratio_cols=[],
            run_sensor_controls=False,
        )

    logger.info("=== Fertig! Oszillationsdaten: %s | statische Daten: %s ===", OUTPUT_DIR, OUTPUT_DIR_STATIC)


def run_pipeline(
    cells: pd.DataFrame,
    output_dir: Path,
    freq_order: list[str] | None,
    intensity_cols: list[str],
    ratio_cols: list[str],
    run_sensor_controls: bool,
) -> None:
    """
    Kompletter Analyse-Kern (Schritte 00-95), unabhängig davon ob er auf den
    Oszillationsdaten oder auf den statischen Kontrolldaten (St.omlp/St.ypd)
    läuft - siehe main(), wo diese Funktion für beide Teilmengen separat
    aufgerufen wird, mit jeweils eigenem output_dir/freq_order.

    'freq_order' steuert die x-Achsen-Reihenfolge der Plots über 'osc_freq'
    (früher die globale Konstante FREQ_ORDER) - bei den statischen Daten
    steht in 'osc_freq' stattdessen 'static_omlp'/'static_ypd' (siehe
    STATIC_ORDER), was hier direkt der gewünschten Vergleichsgruppe
    St.omlp vs. St.ypd entspricht.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Für die 'normalen' Zeitverlaufsplots (10, 30) werden PosCtrl/NegCtrl
    # ausgeblendet - die laufen jetzt eigenständig über die
    # control_consistency-Plots (siehe exclude_controls()). Lineage-
    # Klassifikation, µ/Robustness-Tabellen, Kontroll-Konsistenz-Tests und
    # alle CSV-Exporte laufen weiterhin auf dem VOLLEN 'cells' (mit
    # Kontrollen) - nur die Sichtbarkeit in den Plots ändert sich.
    cells_plot = exclude_controls(cells)

    # ------------------------------------------------------------------
    # 3. Mutter/Bud-Klassifikation VORAB berechnen (nur Berechnung, keine
    #    Ausgabe hier). Notwendig, weil Sektion 10 (µ_event/µ_area) diese
    #    Daten braucht - die zugehörigen Output-Dateien (lineage_events,
    #    Budding Ratio, Panel A, Lineage-Tree) werden aber erst unter den
    #    20er-Präfixen in Sektion 20 geschrieben, damit die Datei-Reihenfolge
    #    der gewünschten inhaltlichen Reihenfolge entspricht.
    #
    #    HEURISTIK ohne echtes Lineage-Tracking - siehe lineage.py Docstring.
    #    Bitte mit inspect_classification() an echten QC-Overlays kalibrieren,
    #    bevor den Ergebnissen für eine Publikation vertraut wird.
    # ------------------------------------------------------------------
    lineage_events = classify_mother_bud(cells, LINEAGE_PARAMS, flux_config=FLUX_CONFIG)
    # ALLE Mutterzellen (auch ohne Budding-Event) - notwendig, damit ruhende
    # Mütter mit budding_ratio=0 in Panel A erscheinen, siehe
    # lineage.compute_budding_ratio()-Docstring.
    mothers = identify_mothers(cells, LINEAGE_PARAMS)

    # ==================================================================
    # 00. Übersicht / Sanity-Check
    # ==================================================================
    overview = data_overview(cells)
    overview.to_csv(output_dir / "00_data_overview.csv", index=False)
    logger.info("Übersicht gespeichert: 00_data_overview.csv")
    plot_n_tracks_overview(overview, output_dir / "00_n_tracks_overview.pdf", freq_order=freq_order)

    # ==================================================================
    # 10. Zellmorphologie & Wachstum (Fläche, µ_event, µ_area)
    # ==================================================================
    if "area" in cells.columns:
        plot_metric_over_time_by_frequency(
            cells_plot, "area", output_dir / "10_cell_area_over_time.pdf",
            freq_order=freq_order, ylabel="Cell area [px²]",
        )

    # -- µ_event: spezifische Wachstumsrate nach Eq. 2 (Blöbaum et al. 2024):
    #    µ = ln(2)/t, aus der Zeit zwischen aufeinanderfolgenden Budding-
    #    Events EINER Mutterzelle. Läuft identisch für alle Bedingungen
    #    inkl. PosCtrl (durchgehend Feast) und NegCtrl (durchgehend
    #    Starvation) - µ ist pro Zelle definiert, nicht pro Bedingung.
    #
    #    style_col="condition_type" wird hier NICHT mehr gesetzt: seit
    #    exclude_controls() PosCtrl/NegCtrl vorher rausfiltert, bleibt in
    #    diesem Plot ohnehin nur noch "Oscillation" übrig (bzw. bei den
    #    statischen Daten St.omlp/St.ypd, was schon über die x-Achse
    #    sichtbar ist) - eine zusätzliche Marker-Form dafür wäre redundant.
    mu_table = compute_specific_growth_rate(
        lineage_events, cells, min_per_frame=MIN_PER_FRAME, mu_max_threshold=MU_MAX_THRESHOLD,
    )
    if not mu_table.empty:
        mu_table.to_csv(output_dir / "11_specific_growth_rate.csv", index=False)
        logger.info("Tabelle gespeichert: 11_specific_growth_rate.csv")

        mu_summary = summarise_growth_rate(mu_table)
        mu_summary.to_csv(output_dir / "11_specific_growth_rate_summary.csv", index=False)
        logger.info("Tabelle gespeichert: 11_specific_growth_rate_summary.csv")

        plot_point_errorbar(
            exclude_controls(mu_summary), value_col="mean_mu", out_path=output_dir / "11_specific_growth_rate.pdf",
            x_col="osc_freq", facet_col="osc_type", color_col=PANEL_A_GROUP_COL,
            x_order=freq_order,
            ylabel="Specific growth rate µ [h⁻¹]",
        )

        # Kontroll-Konsistenz: driften PosCtrl/NegCtrl über die Frequenz-Batches
        # hinweg, obwohl sie eigentlich unabhängig von der Oszillation sein
        # sollten? Test auf Replikat-Ebene (mu_table), Plot auf aggregierter
        # Ebene (mu_summary) - siehe control_consistency.py Docstring.
        mu_control_test = test_control_consistency_across_freq(mu_table, value_col="mu")
        if not mu_control_test.empty:
            mu_control_test.to_csv(output_dir / "11_control_consistency_mu_kruskal.csv", index=False)
            logger.info("Tabelle gespeichert: 11_control_consistency_mu_kruskal.csv")

        plot_control_consistency(
            mu_summary, value_col="mean_mu", out_path=output_dir / "11_control_consistency_mu.pdf",
            x_order=freq_order, ylabel="Specific growth rate µ [h⁻¹] (controls)",
        )
    else:
        logger.warning("Keine µ-Werte berechnet - benötigt Mutterzellen mit >= 2 Budding-Events.")
        mu_summary = pd.DataFrame()  # leer, damit spätere Prüfung (µ_area-Scatter) nicht mit NameError abbricht

    # -- µ_area: flächenbasierte Wachstumsrate
    area_table = compute_area_growth_rate(
        cells, min_per_frame=MIN_PER_FRAME,
        lineage_events=lineage_events,
        mothers=mothers,
    )
    if not area_table.empty:
        area_table.to_csv(output_dir / "12_area_growth_rate.csv", index=False)

        # Separate Summaries für Mütter und Knospen
        area_summary_all = summarise_area_growth(area_table)
        area_summary_mother = summarise_area_growth(area_table, cell_type="mother")
        area_summary_bud = summarise_area_growth(area_table, cell_type="bud")

        area_summary_all.to_csv(output_dir / "12_area_growth_rate_summary_all.csv", index=False)
        area_summary_mother.to_csv(output_dir / "12_area_growth_rate_summary_mother.csv", index=False)
        area_summary_bud.to_csv(output_dir / "12_area_growth_rate_summary_bud.csv", index=False)

        # Plot: µ_area für Mütter (primär interessant) - ohne PosCtrl/NegCtrl,
        # siehe exclude_controls().
        plot_point_errorbar(
            exclude_controls(area_summary_mother), value_col="mean_mu_area",
            out_path=output_dir / "12_area_growth_rate_mother.pdf",
            x_col="osc_freq", facet_col="osc_type", color_col=PANEL_A_GROUP_COL,
            x_order=freq_order,
            ylabel="µ_area, mother cells [h⁻¹]",
        )

        # Scatter: µ_event vs. µ_area (ebenfalls ohne Kontrollen)
        if not mu_summary.empty:
            plot_mu_event_vs_mu_area(
                exclude_controls(mu_summary), exclude_controls(area_summary_mother),
                output_dir / "12_mu_event_vs_mu_area.pdf",
                label_col="osc_freq", facet_col=PANEL_A_GROUP_COL,
                mu_event_table=exclude_controls(mu_table), mu_area_table=exclude_controls(area_table),
            )

    # ==================================================================
    # 20. Lineage/Budding (Ausgabe der in Schritt 3 vorab berechneten
    #     Mutter/Bud-Klassifikation, plus Budding Ratio, Panel A und
    #     Lineage-Tiefe)
    # ==================================================================
    if not lineage_events.empty:
        lineage_events.to_csv(output_dir / "20_budding_events.csv", index=False)
        logger.info("Tabelle gespeichert: 20_budding_events.csv (%d Events)", len(lineage_events))
    else:
        logger.warning("Keine Budding-Events erkannt - Budding Ratio besteht dann ausschließlich aus ruhenden Müttern.")

    # -- Budding Ratio "pro Mutter über die gesamte Beobachtungsdauer"
    #    (NICHT die Paper-Eq.-3-Zeitreihe, siehe weiter unten dafür)
    per_mother, per_experiment = compute_budding_ratio(lineage_events, cells, mothers)
    if not per_mother.empty:
        per_mother.to_csv(output_dir / "21_budding_ratio_per_mother.csv", index=False)
        per_experiment.to_csv(output_dir / "21_budding_ratio_per_experiment.csv", index=False)
        logger.info("Tabellen gespeichert: 21_budding_ratio_per_mother.csv / _per_experiment.csv")

        plot_panel_a(
            cells_plot, exclude_controls(per_mother), output_dir / "21_panel_a_violin.pdf",
            group_col=PANEL_A_GROUP_COL, facet_col=PANEL_A_FACET_COL,
        )
    else:
        logger.warning("compute_budding_ratio() lieferte keine per_mother-Tabelle - Panel A wird übersprungen.")

    # -- Budding Ratio als ZEITREIHE nach Eq. 3 (Blöbaum et al. 2024):
    #    n_buds(t) / n_cells(t-1), analog zu Fig. 3a im Paper.
    #    HINWEIS: bei euch ist die Negativkontrolle durchgehend Starvation
    #    (im Paper: durchgehend Feast) - die Berechnung selbst ist medien-
    #    unabhängig, aber die INTERPRETATION im Vergleich zum Paper ist es
    #    nicht. Siehe budding_ratio_timeseries.py Docstring.
    budding_ts = compute_budding_ratio_timeseries(cells, lineage_events)
    budding_ts.to_csv(output_dir / "22_budding_ratio_timeseries.csv", index=False)
    logger.info("Tabelle gespeichert: 22_budding_ratio_timeseries.csv")

    budding_ts_agg = aggregate_budding_ratio_over_replicates(budding_ts)
    budding_ts_agg.to_csv(output_dir / "22_budding_ratio_timeseries_aggregated.csv", index=False)
    logger.info("Tabelle gespeichert: 22_budding_ratio_timeseries_aggregated.csv")

    plot_budding_ratio_timeseries(
        exclude_controls(budding_ts), output_dir / "22_budding_ratio_timeseries.pdf",
        group_col=PANEL_A_GROUP_COL, facet_col=PANEL_A_FACET_COL, freq_order=freq_order,
    )

    # -- Lineage-Baum & Generationstiefe
    if not lineage_events.empty:
        lineage_tree = build_lineage_tree(lineage_events)
        if not lineage_tree.empty:
            lineage_tree.to_csv(output_dir / "23_lineage_tree.csv", index=False)
            logger.info("Tabelle gespeichert: 23_lineage_tree.csv")

            lineage_depth = summarise_lineage_depth(lineage_tree)
            lineage_depth.to_csv(output_dir / "23_lineage_depth_summary.csv", index=False)
            logger.info(
                "Tabelle gespeichert: 23_lineage_depth_summary.csv (mittlere max. "
                "Generationstiefe=%.1f über %d Kammern)",
                lineage_depth["max_generation"].mean(), len(lineage_depth),
            )

    # ==================================================================
    # 30. Sensor-/Ratio-Zeitverläufe (unterstützend - Grundlage für die
    #     Robustness-Auswertung in Sektion 40)
    # ==================================================================
    if not intensity_cols:
        logger.warning("Keine 'mean_<kanal>' Spalten gefunden - Sensor-Intensitäts-Plots werden übersprungen.")

    for col in intensity_cols:
        plot_metric_over_time_by_frequency(
            cells_plot, col, output_dir / f"30_{col}_over_time.pdf", freq_order=freq_order,
        )

    if not ratio_cols:
        logger.warning("Keine 'ratio_*' Spalten gefunden - Ratio-Plots werden übersprungen.")

    for col in ratio_cols:
        plot_metric_over_time_by_frequency(
            cells_plot, col, output_dir / f"31_{col}_over_time.pdf",
            freq_order=freq_order, ylabel=pretty_label(col),
        )

    # ==================================================================
    # 40. Robustness-Quantifizierung R(t) und R(p) (Eq. 1, Trivellin et al.
    #     2022 / Blöbaum et al. 2024) für die in ROBUSTNESS_VALUE_COLS
    #     konfigurierten Spalten PLUS die zur Laufzeit erkannten
    #     ratio_*-Spalten (Sektion 30) - dadurch fließen die Sensor-Daten
    #     hier in die zusammenfassende Auswertung mit ein.
    #
    #     Zusätzlich (siehe unten, nach der Haupt-Schleife): µ_event und
    #     µ_area aus Sektion 10, JEWEILS NUR IN DER SINNVOLLEN VARIANTE
    #     (µ_event -> R(t) Einzelzelle, µ_area -> R(p)) - die jeweils
    #     andere Variante ist mit der Datenstruktur dieser Tabellen nicht
    #     sauber definierbar, siehe Kommentare dort.
    # ==================================================================
    # dict.fromkeys() statt set(): erhält die Reihenfolge und entfernt
    # Duplikate, falls eine Spalte versehentlich in beiden Listen auftaucht.
    robustness_value_cols = list(dict.fromkeys(ROBUSTNESS_VALUE_COLS + ratio_cols))
    logger.info("Robustness R(t)/R(p) wird berechnet für: %s", robustness_value_cols)

    for value_col in robustness_value_cols:
        if value_col not in cells.columns:
            logger.warning("Robustness: Spalte '%s' nicht in den Daten - übersprungen.", value_col)
            continue

        rt_pop = compute_rt_population(cells, value_col)
        rt_pop.to_csv(output_dir / f"40_Rt_population_{value_col}.csv", index=False)

        rt_cell = compute_rt_single_cell(cells, value_col)
        rt_cell.to_csv(output_dir / f"40_Rt_single_cell_{value_col}.csv", index=False)
        plot_rt_single_cell_distribution(
            exclude_controls(rt_cell), output_dir / f"40_Rt_single_cell_{value_col}.pdf",
            value_col="R_t_single_cell", facet_col="osc_freq", freq_order=freq_order,
        )

        rp = compute_rp(cells, value_col)
        rp.to_csv(output_dir / f"40_Rp_{value_col}.csv", index=False)

        rt_pop_agg = aggregate_robustness_over_replicates(rt_pop, "R_t_population")
        rt_pop_agg.to_csv(output_dir / f"40_Rt_population_{value_col}_aggregated.csv", index=False)

        rp_agg = aggregate_robustness_over_replicates(rp, "R_p")
        rp_agg.to_csv(output_dir / f"40_Rp_{value_col}_aggregated.csv", index=False)

        # Kontrollen raus aus den 'normalen' R(t)/R(p)-Plots - die laufen
        # jetzt separat über die control_consistency-Plots weiter unten.
        rt_pop_agg_plot = exclude_controls(rt_pop_agg)
        rp_agg_plot = exclude_controls(rp_agg)

        plot_point_errorbar(
            rt_pop_agg_plot, value_col="mean", out_path=output_dir / f"40_Rt_population_{value_col}.pdf",
            x_col="osc_freq", facet_col="osc_type", color_col=PANEL_A_GROUP_COL,
            x_order=freq_order,
            ylabel=f"R(t) — {value_col}", title=f"R(t) population level — {value_col}",
        )
        plot_point_errorbar(
            rp_agg_plot, value_col="mean", out_path=output_dir / f"40_Rp_{value_col}.pdf",
            x_col="osc_freq", facet_col="osc_type", color_col=PANEL_A_GROUP_COL,
            x_order=freq_order,
            ylabel=f"R(p) — {value_col}", title=f"R(p) — {value_col}",
        )
        plot_rt_vs_rp_quadrant(
            rt_pop_agg_plot, rp_agg_plot, output_dir / f"40_Rt_vs_Rp_{value_col}.pdf",
            label_col="osc_freq", facet_col=PANEL_A_GROUP_COL,
        )

        # Kontroll-Konsistenz für R(t)/R(p), analog zur Wachstumsrate oben.
        rt_control_test = test_control_consistency_across_freq(rt_pop, value_col="R_t_population")
        if not rt_control_test.empty:
            rt_control_test.to_csv(output_dir / f"40_control_consistency_Rt_population_{value_col}_kruskal.csv", index=False)
        plot_control_consistency(
            rt_pop_agg, value_col="mean", out_path=output_dir / f"40_control_consistency_Rt_population_{value_col}.pdf",
            x_order=freq_order, ylabel=f"R(t) — {value_col} (controls)",
        )

        rp_control_test = test_control_consistency_across_freq(rp, value_col="R_p")
        if not rp_control_test.empty:
            rp_control_test.to_csv(output_dir / f"40_control_consistency_Rp_{value_col}_kruskal.csv", index=False)
        plot_control_consistency(
            rp_agg, value_col="mean", out_path=output_dir / f"40_control_consistency_Rp_{value_col}.pdf",
            x_order=freq_order, ylabel=f"R(p) — {value_col} (controls)",
        )

        logger.info("Robustness R(t)/R(p) für '%s' berechnet und gespeichert.", value_col)

    # -- µ_event -> R(t) Einzelzelle: wie stabil ist die Reproduktionsrate
    #    EINER Mutter über ihre eigenen Budding-Intervalle?
    #    mu_table hat KEINE 'frame'-Spalte (nur frame_start/frame_end pro
    #    Intervall) und ist damit NICHT mit compute_rt_population()/
    #    compute_rp() kompatibel - die brauchen eine gemeinsame Zeitachse
    #    über mehrere Zellen hinweg, die es bei Event-Daten mit ihren pro
    #    Mutter unterschiedlichen Intervall-Zeitpunkten nicht sinnvoll gibt
    #    (Intervall-Enden verschiedener Mütter fallen kaum je zusammen).
    #    compute_rt_single_cell() braucht dagegen KEINE 'frame'-Spalte,
    #    nur mehrere Werte PRO ZELLE - das ist hier gegeben (mehrere
    #    Intervalle pro Mutter), daher funktioniert NUR diese Variante.
    if not mu_table.empty:
        mu_valid = mu_table[~mu_table["mu_is_artefact"]]
        rt_cell_mu = compute_rt_single_cell(mu_valid, value_col="mu", cell_id_col="mother_cell_uid")
        rt_cell_mu.to_csv(output_dir / "40_Rt_single_cell_mu_event.csv", index=False)
        plot_rt_single_cell_distribution(
            exclude_controls(rt_cell_mu), output_dir / "40_Rt_single_cell_mu_event.pdf",
            value_col="R_t_single_cell", facet_col="osc_freq", freq_order=freq_order,
        )
        rt_cell_mu_agg = aggregate_robustness_over_replicates(rt_cell_mu, "R_t_single_cell")
        rt_cell_mu_agg.to_csv(output_dir / "40_Rt_single_cell_mu_event_aggregated.csv", index=False)
        logger.info("R(t) Einzelzelle für µ_event berechnet und gespeichert (%d Mütter).", len(rt_cell_mu))
    else:
        logger.warning("mu_table ist leer - R(t) Einzelzelle für µ_event wird übersprungen.")

    # -- µ_area -> R(p): wie homogen ist die flächenbasierte Wachstumsrate
    #    ÜBER DIE ZELLEN einer Kammer?
    #    area_table hat pro Zelle nur EINEN Wert (der ganze Track wurde
    #    bereits zu einer Steigung verdichtet) - eine Zeitachse gibt es
    #    hier nicht mehr, R(t) ist für µ_area daher NICHT definierbar
    #    (weder Populations- noch Einzelzell-Variante: beide bräuchten
    #    mehrere Zeitpunkte pro Zelle bzw. pro Kammer).
    #    compute_rp() gruppiert technisch nach exp_id x frame - da wir
    #    keine Zeit haben, setzen wir 'frame' konstant auf 0. Dadurch
    #    entspricht jede Kammer genau einer Gruppe, und R(p) misst wie
    #    beabsichtigt die Homogenität ÜBER DIE ZELLEN einer Kammer (statt
    #    über Zeitpunkte). Nur zuverlässige Fits (fit_is_reliable) gehen
    #    ein - unzuverlässige Fits würden R(p) sonst mit Rauschen aus
    #    schlecht bestimmten Steigungen aufblähen. Mutter- und Knospen-
    #    Tracks werden NICHT getrennt (analog zu 'area'/'eccentricity'
    #    oben, die ebenfalls nicht nach cell_type filtern).
    if not area_table.empty:
        area_reliable = area_table[area_table["fit_is_reliable"]].copy()
        area_reliable["frame"] = 0
        rp_mu_area = compute_rp(area_reliable, value_col="mu_area")
        rp_mu_area.to_csv(output_dir / "40_Rp_mu_area.csv", index=False)

        rp_mu_area_agg = aggregate_robustness_over_replicates(rp_mu_area, "R_p")
        rp_mu_area_agg.to_csv(output_dir / "40_Rp_mu_area_aggregated.csv", index=False)

        plot_point_errorbar(
            exclude_controls(rp_mu_area_agg), value_col="mean", out_path=output_dir / "40_Rp_mu_area.pdf",
            x_col="osc_freq", facet_col="osc_type", color_col=PANEL_A_GROUP_COL,
            x_order=freq_order,
            ylabel="R(p) — µ_area", title="R(p) — µ_area (homogeneity across cells)",
        )

        rp_mu_area_control_test = test_control_consistency_across_freq(rp_mu_area, value_col="R_p")
        if not rp_mu_area_control_test.empty:
            rp_mu_area_control_test.to_csv(output_dir / "40_control_consistency_Rp_mu_area_kruskal.csv", index=False)
        plot_control_consistency(
            rp_mu_area_agg, value_col="mean", out_path=output_dir / "40_control_consistency_Rp_mu_area.pdf",
            x_order=freq_order, ylabel="R(p) — µ_area (controls)",
        )
        logger.info("R(p) für µ_area berechnet und gespeichert (%d zuverlässige Tracks).", len(area_reliable))
    else:
        logger.warning("area_table ist leer - R(p) für µ_area wird übersprungen.")

    # ==================================================================
    # 50. Zusammenfassungstabelle
    # ==================================================================
    summary = summary_statistics(cells, intensity_cols)
    summary.to_csv(output_dir / "50_summary_statistics.csv", index=False)
    logger.info("Tabelle gespeichert: 50_summary_statistics.csv")

    # ==================================================================
    # 90. ANHANG: Morphologie-Scatter, Einzelzell-Trajektorien, stabile
    #     Mutter-Trajektorien - deskriptive Plausibilitäts-Checks ohne
    #     eigene quantitative Kennzahl, siehe Diskussion zur Restrukturierung.
    # ==================================================================
    plot_morphology_scatter(cells, output_dir / "90_morphology_scatter.pdf")

    # Dieser Plot zeigt 'area', hing aber an intensity_cols - dadurch fehlte er
    # bei den statischen Daten komplett (Wildtyp, keine Fluoreszenzkanäle).
    if "area" in cells.columns:
        plot_single_cell_trajectories(cells, "area", output_dir / "91_single_cell_trajectories.pdf")
    else:
        logger.warning("Spalte 'area' fehlt - 91_single_cell_trajectories.pdf übersprungen.")

    # Stabil getrackte Mütter im Detail: pro Bedingungs-Kombination
    # (STABLE_MOTHER_GROUP_COLS) EINE stabile Mutter, mit Budding-Markern
    # und nur ihrem eigenen Ratio-Kanal (siehe mother_trajectories.py) -
    # ergänzt die aggregierten Auswertungen oben (Panel A, Fig. 3a) um
    # eine Einzelzell-Ebene zur visuellen Plausibilitätsprüfung.
    stable_base_cols = [c for c in STABLE_MOTHER_BASE_VALUE_COLS if c in cells.columns]
    if not mothers.empty and stable_base_cols:
        plot_stable_mother_per_group(
            cells, lineage_events, mothers,
            out_path=output_dir / "92_stable_mother_trajectories.pdf",
            group_cols=STABLE_MOTHER_GROUP_COLS, base_value_cols=stable_base_cols,
            min_coverage=STABLE_MOTHER_MIN_COVERAGE, min_per_frame=MIN_PER_FRAME,
        )
        logger.info("Plot gespeichert: 92_stable_mother_trajectories.pdf")
    else:
        logger.warning("Keine stabilen Mütter oder keine base_value_cols - 92_stable_mother_trajectories.pdf übersprungen.")

    # ==================================================================
    # 95. ANHANG: Sensor-Controls (PosCtrl-vs-NegCtrl-Validierung pro
    #     Biosensor). Reine Kalibrierungs-/QC-Auswertung des Sensors selbst,
    #     nicht Teil der biologischen Kernaussage - siehe Diskussion zur
    #     Restrukturierung. Zellwerte werden vorab pro Replikat summiert,
    #     um Pseudoreplikation zu vermeiden.
    #
    #     Wird für die statischen Daten übersprungen (run_sensor_controls=
    #     False): dort wird nur der Wildtyp kultiviert, es gibt also keine
    #     Fluoreszenzkanäle/Sensor-Ratios auszuwerten.
    # ==================================================================
    if not run_sensor_controls:
        logger.info("Sensor-Control-Validierung (Schritt 95) übersprungen (run_sensor_controls=False).")
        SENSOR_CONFIG_ITEMS = {}
    else:
        SENSOR_CONFIG_ITEMS = SENSOR_CONFIG

    for biosensor, sensor_cfg in SENSOR_CONFIG_ITEMS.items():
        ratio_col = sensor_cfg.ratio_name
        if ratio_col not in cells.columns or not cells[ratio_col].notna().any():
            logger.info("Sensor-control validation: no usable '%s' values for %s; skipped.", ratio_col, biosensor)
            continue

        sensor_data = cells[cells["biosensor"] == biosensor]
        for osc_type in sorted(sensor_data["osc_type"].dropna().unique()):
            analysis_start_min = OSCILLATION_START_MIN
            controls = prepare_sensor_controls(
                cells, value_col=ratio_col, biosensor=biosensor, osc_type=osc_type,
            )
            control_time, control_summary = summarise_sensor_controls(
                controls, value_col=ratio_col, analysis_start_min=analysis_start_min,
            )
            if control_time.empty:
                continue

            output_stem = f"95_{biosensor}_{osc_type}_{ratio_col}"
            control_time.to_csv(output_dir / f"{output_stem}_timeseries_per_replicate.csv", index=False)
            control_summary.to_csv(output_dir / f"{output_stem}_summary_per_replicate.csv", index=False)
            sensor_label = ratio_col.replace("ratio_", "")
            plot_sensor_control_timeseries(
                control_time, output_dir / f"{output_stem}_timeseries.pdf",
                sensor_label=sensor_label,
                preconditioning_end_min=OSCILLATION_START_MIN,
            )
            plot_sensor_control_comparison(
                control_summary, output_dir / f"{output_stem}_comparison.pdf",
                sensor_label=sensor_label, analysis_start_min=analysis_start_min,
                control_labels=CONTROL_CONCENTRATION_LABELS.get(osc_type),
            )

            # A flat ratio can arise because both raw channels shift together
            # or because neither channel contains a detectable control effect.
            # The two raw-channel panels make this distinction inspectable.
            channel_a_controls = prepare_sensor_controls(
                cells, value_col=sensor_cfg.channel_a, biosensor=biosensor, osc_type=osc_type,
            )
            channel_b_controls = prepare_sensor_controls(
                cells, value_col=sensor_cfg.channel_b, biosensor=biosensor, osc_type=osc_type,
            )
            channel_a_time, _ = summarise_sensor_controls(
                channel_a_controls, value_col=sensor_cfg.channel_a, analysis_start_min=analysis_start_min,
            )
            channel_b_time, _ = summarise_sensor_controls(
                channel_b_controls, value_col=sensor_cfg.channel_b, analysis_start_min=analysis_start_min,
            )
            plot_sensor_raw_channel_timeseries(
                channel_a_time, channel_b_time, output_dir / f"{output_stem}_raw_channels.pdf",
                sensor_label=sensor_label, channel_a_label=sensor_cfg.channel_a,
                channel_b_label=sensor_cfg.channel_b,
                preconditioning_end_min=OSCILLATION_START_MIN,
            )

            chamber_summary = summarise_sensor_controls_by_chamber(
                controls, value_col=ratio_col, analysis_start_min=analysis_start_min,
            )
            if not chamber_summary.empty:
                chamber_summary.to_csv(output_dir / f"{output_stem}_control_chambers.csv", index=False)
                plot_control_chamber_comparison(
                    chamber_summary, output_dir / f"{output_stem}_control_chambers.pdf",
                    sensor_label=sensor_label,
                )

    logger.info("=== Fertig! Alle Plots & Tabellen in: %s ===", output_dir)


if __name__ == "__main__":
    main()
