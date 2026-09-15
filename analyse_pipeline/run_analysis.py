"""
run_analysis.py
===============
Hauptskript der Analyse-Pipeline.

    python run_analysis.py                    # alles
    python run_analysis.py --list-steps       # welche Schritte gibt es?
    python run_analysis.py --steps 13 40      # nur diese neu rechnen

Erwartete Ordnerstruktur unter DATA_ROOT:

    Data/<Biosensor>/<Oszillationstyp>/<Oszillationsperiode>/03_results/Combined_Results.{csv,parquet}

z.B.:
    Data/BSG/Glc/1.5/03_results/Combined_Results.csv
    Data/BSA/pH/24/03_results/Combined_Results.csv

Dieses Modul macht nur noch das Drumherum: laden, QC anwenden, Oszillations-
von statischen Daten trennen, und dann die Schritte aus pipeline_steps.py
ausfuehren. Die Auswertung selbst steht dort, ein Schritt pro Funktion.

Ablauf:
    1. Alle Combined_Results einlesen (gecacht als .parquet)
    2. Track-Merges DANN QC-Exclusions anwenden (nicht-destruktiv)
    3. Oszillations- und statische Daten trennen, Morphospace-Referenz bilden
    4. pipeline_steps.STEPS auf beide Teilmengen anwenden

Konfiguration: config.py (oder AUREO_DATA_ROOT / AUREO_OUTPUT_DIR).
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

# ==============================================================================
# 0. KONFIGURATION
# ==============================================================================
# Alle Pfade und Parameter stehen in config.py - dort (oder über die
# Umgebungsvariablen AUREO_DATA_ROOT / AUREO_OUTPUT_DIR) anpassen, NICHT hier.
# inspect_lineage.py und validate_lineage.py lesen dieselbe Datei, damit
# Kalibrierung, Validierung und Auswertung mit denselben Schwellen laufen.
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
    MORPHOLOGY_REFERENCE_CONDITIONS,
    MORPHOLOGY_PERCENTILE,
    log_active_configuration,
)
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
from morphology import derive_thresholds_from_reference
from analysis import add_time_column, find_intensity_columns
from pipeline_steps import STEPS, PipelineContext

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


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


def select_steps(patterns: list[str] | None):
    """Waehlt Schritte anhand von Teilstrings ihres Schluessels aus.

    '13' trifft '13_morphology', 'morph' ebenso. Ein Muster ohne Treffer ist
    ein Fehler und kein stilles Ueberspringen - sonst laeuft man versehentlich
    eine leere Pipeline und haelt das Ergebnis fuer aktuell.
    """
    if not patterns:
        return list(STEPS)

    selected, unknown = [], []
    for pattern in patterns:
        hits = [s for s in STEPS if pattern.lower() in s.key.lower()]
        if not hits:
            unknown.append(pattern)
        selected.extend(h for h in hits if h not in selected)

    if unknown:
        raise SystemExit(
            f"Unbekannte Schritt(e): {unknown}\n"
            f"Verfuegbar: {', '.join(s.key for s in STEPS)}\n"
            f"(--list-steps zeigt die Beschreibungen)"
        )
    # Immer in Pipeline-Reihenfolge ausfuehren, egal wie die Muster stehen.
    return [s for s in STEPS if s in selected]


def run_steps(ctx: PipelineContext, steps) -> None:
    """Fuehrt die gewaehlten Schritte auf einem Kontext aus.

    Ein Fehler in einem Schritt bricht den Rest NICHT ab: bei einem Lauf ueber
    Nacht ist es deutlich nuetzlicher, die uebrigen Auswertungen zu haben und
    am Ende eine Liste der gescheiterten Schritte zu lesen, als alles zu
    verlieren, weil eine Abbildung eine leere Facette hatte. Der Traceback
    steht vollstaendig im Log, und der Exit-Code am Ende ist ungleich 0.
    """
    ctx.output_dir.mkdir(parents=True, exist_ok=True)
    failed = []
    for step in steps:
        logger.info("--- Schritt %s: %s ---", step.key, step.title)
        try:
            step.run(ctx)
        except Exception:
            logger.exception("Schritt '%s' fehlgeschlagen - wird uebersprungen.", step.key)
            failed.append(step.key)
    if failed:
        logger.error("Fehlgeschlagene Schritte in %s: %s", ctx.output_dir, ", ".join(failed))
    ctx.failed_steps = failed


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Analyse-Pipeline fuer Einzelzell-Mikroskopiedaten. "
                    "Pfade und Parameter kommen aus config.py bzw. aus "
                    "AUREO_DATA_ROOT / AUREO_OUTPUT_DIR.",
    )
    parser.add_argument(
        "--steps", nargs="+", metavar="MUSTER", default=None,
        help="Nur diese Schritte ausfuehren, als Teilstring des Schluessels "
             "(z.B. '13 40' oder 'morphology'). Standard: alle.",
    )
    parser.add_argument(
        "--list-steps", action="store_true",
        help="Verfuegbare Schritte auflisten und beenden.",
    )
    parser.add_argument(
        "--skip-static", action="store_true",
        help="Die statischen Daten (osc_type == 'static') nicht auswerten.",
    )
    args = parser.parse_args(argv)

    if args.list_steps:
        print("Verfuegbare Schritte (Reihenfolge = Ausfuehrungsreihenfolge):\n")
        for s in STEPS:
            print(f"  {s.key:<20} {s.title}")
        print("\nBeispiel:  python run_analysis.py --steps 13 40")
        return 0

    steps = select_steps(args.steps)

    # Aktive Konfiguration inkl. aller Abweichungen vom dokumentierten Standard
    # ganz am Anfang ins Log - siehe config.log_active_configuration().
    log_active_configuration()
    if args.steps:
        logger.info("Teillauf: nur %s", ", ".join(s.key for s in steps))

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
    # daher hier schon einmal bestimmen. Wird von Sektion 30 (Plots), 40
    # (Robustness) und 50/90 (Summary/Supplement) weiterverwendet.
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
    logger.info("  Oszillationsperioden:   %s", sorted(cells["osc_freq"].dropna().unique()))

    # ------------------------------------------------------------------
    # 2b. Statische Daten (Data/<Biosensor>/static/static_<medium>/...,
    #     condition St.omlp/St.ypd) von den Oszillationsdaten trennen - ab
    #     hier laufen beide komplett getrennt durch dieselben Schritte,
    #     landen aber in eigenen Output-Ordnern/Dateien.
    # ------------------------------------------------------------------
    # .str.strip().str.lower() statt striktem "==" Vergleich: robust gegen
    # Groß-/Kleinschreibung ("Static"/"STATIC") oder versehentliche
    # Leerzeichen im Ordnernamen - sonst landen statische Zeilen fälschlich
    # (und unbemerkt) in cells_osc statt in cells_static.
    is_static = cells["osc_type"].astype(str).str.strip().str.lower() == "static"
    cells_osc = cells[~is_static].copy()
    cells_static = cells[is_static].copy()

    logger.info(
        "Aufgeteilt: %d Zeilen Oszillationsdaten, %d Zeilen statische Daten (osc_type == 'static').",
        len(cells_osc), len(cells_static),
    )

    # ------------------------------------------------------------------
    # 2c. Zell-Positionen für die QC-Kalibrierung exportieren
    # ------------------------------------------------------------------
    # plot_qc_lineage_overlay.py und validate_lineage.py brauchen den
    # Zelldatensatz NACH Track-Merges und Exclusions - genau den Stand, auf dem
    # classify_mother_bud() gelaufen ist. Der Parquet-Cache ist dagegen der
    # Stand VOR QC. Metadaten mit exportieren: validate_lineage.py gruppiert
    # die Erkennungsrate danach (osc_freq/osc_type/...), und ohne diese Spalten
    # faellt genau der Konfundierungs-Test aus, um den es dort geht.
    qc_position_cols = [c for c in
                        ["exp_id", "cell_uid", "track_id", "frame", "centroid_x", "centroid_y",
                         "area", "filename",
                         "biosensor", "osc_type", "osc_freq", "condition", "replicate", "chamber"]
                        if c in cells.columns]
    qc_positions_path = OUTPUT_DIR / "00_cell_positions.parquet"
    cells[qc_position_cols].to_parquet(qc_positions_path, index=False)
    logger.info(
        "Tabelle gespeichert: %s (%d Zeilen) - Eingabe für plot_qc_lineage_overlay.py "
        "und validate_lineage.py.",
        qc_positions_path.name, len(cells),
    )

    # ------------------------------------------------------------------
    # 2d. Morphospace-Normalfenster EINMAL aus dem Gesamtdatensatz ableiten
    # ------------------------------------------------------------------
    # Bewusst hier und nicht pro Teil-Pipeline: die Referenz (PosCtrl) gibt es
    # nur bei den Oszillationsdaten. Wuerde jede Teilmenge ihre eigenen
    # Schwellen bilden, waeren die aberranten Anteile von Oszillations- und
    # statischen Daten nicht mehr miteinander vergleichbar.
    morph_thresholds = derive_thresholds_from_reference(
        cells,
        reference_condition_types=MORPHOLOGY_REFERENCE_CONDITIONS,
        percentile=MORPHOLOGY_PERCENTILE,
    )

    # ------------------------------------------------------------------
    # 3. Schritte ausfuehren - erst Oszillation, dann statisch
    # ------------------------------------------------------------------
    failed = []

    osc_ctx = PipelineContext(
        cells=cells_osc, output_dir=OUTPUT_DIR,
        freq_order=resolve_x_order(cells_osc, "osc_freq", FREQ_ORDER),
        intensity_cols=intensity_cols, ratio_cols=ratio_cols,
        run_sensor_controls=True, morph_thresholds=morph_thresholds,
    )
    run_steps(osc_ctx, steps)
    failed += [f"Oszillation/{k}" for k in osc_ctx.failed_steps]

    if args.skip_static:
        logger.info("Statische Daten uebersprungen (--skip-static).")
    elif cells_static.empty:
        logger.warning(
            "Keine statischen Daten (osc_type == 'static') gefunden - Schritt für "
            "statische Daten übersprungen."
        )
    else:
        # Keine Sensor-/Ratio-Auswertung für statische Daten: hier wird nur
        # der Wildtyp kultiviert, es gibt also keine Fluoreszenzkanäle -
        # daher intensity_cols/ratio_cols leer und run_sensor_controls=False,
        # statt auf durchgehend NaN-Spalten zu rechnen.
        static_ctx = PipelineContext(
            cells=cells_static, output_dir=OUTPUT_DIR_STATIC,
            freq_order=resolve_x_order(cells_static, "osc_freq", STATIC_ORDER),
            intensity_cols=[], ratio_cols=[],
            run_sensor_controls=False, morph_thresholds=morph_thresholds,
        )
        run_steps(static_ctx, steps)
        failed += [f"statisch/{k}" for k in static_ctx.failed_steps]

    if failed:
        logger.error("=== Mit Fehlern beendet. Fehlgeschlagen: %s ===", ", ".join(failed))
        return 1

    logger.info("=== Fertig! Oszillationsdaten: %s | statische Daten: %s ===",
                OUTPUT_DIR, OUTPUT_DIR_STATIC)
    return 0


if __name__ == "__main__":
    sys.exit(main())
