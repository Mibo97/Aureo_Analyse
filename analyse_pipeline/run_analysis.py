"""
run_analysis.py
===============
Hauptskript der Analyse-Pipeline.

    python run_analysis.py                    # alles
    python run_analysis.py --list-steps       # welche Schritte gibt es?
    python run_analysis.py --steps 20 40      # nur diese neu rechnen

Erwartete Ordnerstruktur unter DATA_ROOT:

    Data/<Biosensor>/<Oszillationstyp>/<Oszillationsperiode>/03_results/Combined_Results.{csv,parquet}

z.B.:
    Data/BSG/Glc/1.5/03_results/Combined_Results.csv
    Data/BSA/pH/24/03_results/Combined_Results.csv

Dieses Modul macht nur noch das Drumherum: laden, QC anwenden, die drei
Zweige trennen, und dann die Schritte aus pipeline_steps.py ausfuehren. Die
Auswertung selbst steht dort, ein Schritt pro Funktion.

DREI UNABHAENGIGE ZWEIGE, drei Output-Ordner:

    Oszillation   analysis_output/          osc_type != 'static', biosensor != 'PKO'
    statisch      analysis_output/static/   osc_type == 'static'
    PKO           analysis_output/pko/      biosensor == 'PKO'  (Data/PKO/...)

Sie laufen durch DIESELBEN Schritte, teilen aber keine Zahlen: jeder Zweig
bekommt seinen eigenen PipelineContext. Der PKO-Zweig kommt zusaetzlich mit
einem WT-gegen-PKO-Vergleich (pko_comparison.py) fuer die Clogging-Hypothese -
der braucht beide Teilmengen und laeuft deshalb ausserhalb der Schritt-Schleife.

Ablauf:
    1. Alle Combined_Results einlesen (gecacht als .parquet)
    2. Track-Merges DANN QC-Exclusions anwenden (nicht-destruktiv)
    3. PKO abtrennen, dann Oszillations- und statische Daten trennen
    4. pipeline_steps.STEPS auf alle drei Teilmengen anwenden
    5. WT gegen PKO vergleichen, falls PKO-Daten vorhanden sind

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
    OUTPUT_DIR_PKO,
    PKO_BIOSENSOR_NAME,
    CACHE_PATH,
    QC_EXCLUSIONS_PATH,
    FORCE_RELOAD,
    MIN_PER_FRAME,
    FREQ_ORDER,
    STATIC_CHIP_LABELS,
    STATIC_MEDIUM_PREFIX,
    STATIC_MEDIUM_ORDER,
    PANEL_A_GROUP_COL_STATIC,
    OSCILLATION_START_MIN,
    LINEAGE_PARAMS,
    BUD_MAX_AREA_FRACTION_FALLBACK,
    BUD_SIZE_PLAUSIBLE_RANGE,
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
    qc_batches,
    find_qc_conflicts,
)
from qc_comparison import compare_qc_runs, qc_exclusion_inventory
from sensors import compute_ratios, SENSOR_CONFIG
from analysis import add_time_column, find_intensity_columns
from pipeline_steps import STEPS, PipelineContext
from pko_comparison import run_pko_comparison
from experiment_units import add_experiment_units, chip_overview
from bud_size import run_bud_size_threshold

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

    '40' trifft '40_robustness', 'robust' ebenso. Ein Muster ohne Treffer ist
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
             "(z.B. '20 40' oder 'robustness'). Standard: alle.",
    )
    parser.add_argument(
        "--list-steps", action="store_true",
        help="Verfuegbare Schritte auflisten und beenden.",
    )
    parser.add_argument(
        "--skip-static", action="store_true",
        help="Die statischen Daten (osc_type == 'static') nicht auswerten.",
    )
    parser.add_argument(
        "--skip-qc-comparison", action="store_true",
        help="Den Lauf auf den Rohdaten (ohne Merges/Exclusions) und den Vergleich "
             "mit/ohne QC fuer die QC-beruehrten Batches weglassen.",
    )
    parser.add_argument(
        "--skip-pko", action="store_true",
        help=f"Die PKO-Daten (biosensor == '{PKO_BIOSENSOR_NAME}') nicht auswerten, "
             "inklusive des WT-gegen-PKO-Vergleichs.",
    )
    args = parser.parse_args(argv)

    if args.list_steps:
        print("Verfuegbare Schritte (Reihenfolge = Ausfuehrungsreihenfolge):\n")
        for s in STEPS:
            print(f"  {s.key:<20} {s.title}")
        print("\nBeispiel:  python run_analysis.py --steps 20 40")
        return 0

    steps = select_steps(args.steps)

    # Aktive Konfiguration inkl. aller Abweichungen vom dokumentierten Standard
    # ganz am Anfang ins Log - siehe config.log_active_configuration().
    log_active_configuration()
    if args.steps:
        logger.info("Teillauf: nur %s", ", ".join(s.key for s in steps))

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR_STATIC.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR_PKO.mkdir(parents=True, exist_ok=True)

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
    # 2. QC: Exclusions, Track-Merges, Exclusions (nicht-destruktiv)
    # ------------------------------------------------------------------
    # WICHTIG: Track-Merges müssen VOR jeglicher
    # Lineage-Klassifikation laufen - sonst sieht classify_mother_bud() weiterhin
    # zwei getrennte Tracks und ein Tracking-Bruch (z.B. durch starke
    # Zellbewegung) wird fälschlich als neuer Bud interpretiert.
    init_qc_exclusions(QC_EXCLUSIONS_PATH)  # legt leere Datei an, falls noch keine existiert
    exclusions = read_qc_exclusions(QC_EXCLUSIONS_PATH)

    # Rohdaten VOR jedem QC festhalten - fuer den Vergleich mit/ohne QC. Die
    # Ratio-Spalten sind schon da (nicht-destruktiv), die Zeitspalte und die
    # Versuchseinheiten kommen unten noch dazu.
    cells_raw = cells.copy()

    conflicts = find_qc_conflicts(exclusions)
    if not conflicts.empty:
        (OUTPUT_DIR / "qc_comparison").mkdir(parents=True, exist_ok=True)
        conflicts.to_csv(OUTPUT_DIR / "qc_comparison" / "70_qc_exclusions_conflicts.csv", index=False)

    qc_summary = summarise_qc_exclusions(exclusions)
    if not qc_summary.empty:
        logger.info("QC-Exclusions pro Experiment:\n%s", qc_summary.to_string(index=False))

    # Ausschluesse ZWEIMAL: vor den Merges auf die urspruengliche cell_uid (sonst
    # verliert ein Track, der gemergt UND ausgeschlossen ist, seinen Ausschluss,
    # weil der Merge ihn umbenennt), und nach den Merges noch einmal, damit ein
    # Ausschluss auf dem ZIEL-Track auch die hineingemergten Frames trifft.
    # Ausschluesse sind idempotent - der zweite Durchlauf kostet nichts Falsches.
    cells = apply_qc_exclusions(cells, exclusions, mode="remove")

    merge_summary = summarise_track_merges(exclusions)
    if not merge_summary.empty:
        logger.info("Konfigurierte Track-Merges:\n%s", merge_summary.to_string(index=False))
    cells = apply_track_merges(cells, exclusions)

    cells = apply_qc_exclusions(cells, exclusions, mode="remove")
    cells = add_time_column(cells, MIN_PER_FRAME)

    # Versuchseinheiten: chip / chip_family / medium / date. Erst hier, nach
    # dem QC, damit die Kammerzahlen pro Chip den ausgewerteten Stand zeigen.
    # WICHTIG fuer alles Weitere: 'replicate' ist ein Array-Index, kein
    # Replikat - siehe experiment_units.py.
    cells = add_experiment_units(cells, STATIC_CHIP_LABELS, static_medium_prefix=STATIC_MEDIUM_PREFIX)
    overview_chips = chip_overview(cells)
    overview_chips.to_csv(OUTPUT_DIR / "00_chip_overview.csv", index=False)
    logger.info("Tabelle gespeichert: 00_chip_overview.csv (%d Chips)", len(overview_chips))

    failed: list[str] = []

    # ------------------------------------------------------------------
    # 2a. Groessenkriterium der Mutter/Bud-Heuristik: EINE Schwelle aus
    #     allen Daten nach QC, fuer alle Zweige dieselbe (bud_size.py).
    # ------------------------------------------------------------------
    # Angespuelte Blastokonidien tauchen "neu" neben sitzenden Zellen auf und
    # bestehen die raeumliche Zuordnung wie eine Knospe - sind aber beim ersten
    # Auftreten etwa so gross wie die vermeintliche Mutter. Die Schwelle ist
    # der Antimodus der zweigipfligen Verteilung von bud_area / mother_area;
    # Tabelle, Abbildung und die gewaehlte Quelle ('source') landen in
    # 20_bud_size_threshold.csv / 20_bud_size_at_appearance.*.
    try:
        bud_size_threshold = run_bud_size_threshold(
            cells, OUTPUT_DIR, params=LINEAGE_PARAMS,
            fallback=BUD_MAX_AREA_FRACTION_FALLBACK, plausible=BUD_SIZE_PLAUSIBLE_RANGE,
        )
    except Exception:
        logger.exception(
            "Groessenkriterium konnte nicht aus den Daten abgeleitet werden - "
            "Rueckfallwert %.2f wird verwendet.", BUD_MAX_AREA_FRACTION_FALLBACK,
        )
        bud_size_threshold = BUD_MAX_AREA_FRACTION_FALLBACK
        failed.append("bud_size_threshold")

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
    # PKO ZUERST abtrennen, VOR der statisch/Oszillations-Trennung: PKO liegt
    # auf der Biosensor-Ebene (Data/PKO/...) und landet damit in 'biosensor',
    # nicht in 'osc_type'. Ohne diesen Schnitt liefe PKO durch is_static
    # hindurch in cells_osc und erschiene - weil PANEL_A_GROUP_COL == 'biosensor'
    # auch color_col ist - als zusaetzliche Farbe in JEDEM bestehenden
    # Oszillations-Plot und als zusaetzliches Violin in Panel A. Die
    # vorhandenen Ergebnisse wuerden sich also aendern, was der Sinn eines
    # dritten, unabhaengigen Zweigs gerade nicht ist. Siehe config.PKO_BIOSENSOR_NAME.
    is_pko = cells["biosensor"].astype(str).str.strip().str.upper() == PKO_BIOSENSOR_NAME.upper()
    cells_pko = cells[is_pko].copy()
    cells_rest = cells[~is_pko]

    is_static = cells_rest["osc_type"].astype(str).str.strip().str.lower() == "static"
    cells_osc = cells_rest[~is_static].copy()
    cells_static = cells_rest[is_static].copy()

    logger.info(
        "Aufgeteilt: %d Zeilen Oszillationsdaten, %d Zeilen statische Daten "
        "(osc_type == 'static'), %d Zeilen PKO-Daten (biosensor == '%s').",
        len(cells_osc), len(cells_static), len(cells_pko), PKO_BIOSENSOR_NAME,
    )
    if cells_pko.empty:
        logger.info(
            "Keine PKO-Daten gefunden - der PKO-Zweig und der WT-gegen-PKO-Vergleich "
            "entfallen. Das ist der normale Zustand, solange unter %s/%s/ keine Daten "
            "liegen; die uebrigen Auswertungen sind davon unberuehrt.",
            DATA_ROOT.name, PKO_BIOSENSOR_NAME,
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
                         "biosensor", "osc_type", "osc_freq", "condition", "replicate", "chamber",
                         "chip", "chip_family", "medium", "date"]
                        if c in cells.columns]
    qc_positions_path = OUTPUT_DIR / "00_cell_positions.parquet"
    cells[qc_position_cols].to_parquet(qc_positions_path, index=False)
    logger.info(
        "Tabelle gespeichert: %s (%d Zeilen) - Eingabe für plot_qc_lineage_overlay.py "
        "und validate_lineage.py.",
        qc_positions_path.name, len(cells),
    )

    # ------------------------------------------------------------------
    # 3. Schritte ausfuehren - erst Oszillation, dann statisch, dann PKO
    # ------------------------------------------------------------------
    osc_freq_order = resolve_x_order(cells_osc, "osc_freq", FREQ_ORDER)
    osc_ctx = PipelineContext(
        cells=cells_osc, output_dir=OUTPUT_DIR,
        freq_order=osc_freq_order,
        intensity_cols=intensity_cols, ratio_cols=ratio_cols,
        run_sensor_controls=True,
        run_control_consistency=len(osc_freq_order) >= 2,
        bud_size_threshold=bud_size_threshold,
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
        static_freq_order = resolve_x_order(cells_static, "medium", STATIC_MEDIUM_ORDER)
        static_ctx = PipelineContext(
            cells=cells_static, output_dir=OUTPUT_DIR_STATIC,
            freq_order=static_freq_order,
            # Medium auf der x-Achse, Chip-Familie (W109/W65) als Facette.
            x_col="medium", facet_col="chip_family", panel_a_facet_col="chip_family",
            intensity_cols=[], ratio_cols=[],
            run_sensor_controls=False,
            run_control_consistency=len(static_freq_order) >= 2,
            # Panel A muss hier nach 'osc_freq' gruppieren: bei den statischen
            # Daten steht die Vergleichsgruppe (static_omlp/static_ypd) dort,
            # und plot_panel_a() sieht nur group_col/facet_col. Mit dem
            # Default 'biosensor' landeten BEIDE Medien in EINEM Violin, und
            # die Abbildung konnte die Frage des statischen Experiments
            # ("komplexes vs. minimales Medium") gar nicht beantworten.
            panel_a_group_col=PANEL_A_GROUP_COL_STATIC,
            bud_size_threshold=bud_size_threshold,
        )
        run_steps(static_ctx, steps)
        failed += [f"statisch/{k}" for k in static_ctx.failed_steps]

    # ------------------------------------------------------------------
    # 3b. PKO: dritter, unabhaengiger Zweig (Clogging-Hypothese)
    # ------------------------------------------------------------------
    # Laeuft durch DIESELBEN Schritte wie die statischen Daten und landet in
    # einem eigenen Output-Ordner. Wie bei static: keine Fluoreszenzkanaele,
    # also intensity_cols/ratio_cols leer und run_sensor_controls=False -
    # die Schritte 30/31/95 fallen damit von selbst aus.
    #
    # Zusaetzlich run_control_consistency=False, sobald PKO nur EINE Periode
    # abdeckt: test_control_consistency_across_freq() braucht mindestens zwei
    # osc_freq-Batches und liefert sonst p = NaN, und plot_control_consistency()
    # zeichnet dann einen einzelnen Punkt pro Kontrollart - eine trivial flache
    # Linie, die wie "PKO-Kontrollen sind konsistent" aussieht, obwohl sie
    # keine Information enthaelt. Der eigentliche Vergleich laeuft stattdessen
    # ueber pko_comparison.run_pko_comparison() (siehe dort).
    if args.skip_pko:
        logger.info("PKO-Daten uebersprungen (--skip-pko).")
    elif cells_pko.empty:
        pass  # oben schon gemeldet
    else:
        pko_freq_order = resolve_x_order(cells_pko, "osc_freq", FREQ_ORDER)
        pko_ctx = PipelineContext(
            cells=cells_pko, output_dir=OUTPUT_DIR_PKO,
            freq_order=pko_freq_order,
            intensity_cols=[], ratio_cols=[],
            run_sensor_controls=False,
            run_control_consistency=len(pko_freq_order) >= 2,
            bud_size_threshold=bud_size_threshold,
        )
        run_steps(pko_ctx, steps)
        failed += [f"PKO/{k}" for k in pko_ctx.failed_steps]

        # WT gegen PKO: braucht BEIDE Teilmengen und laeuft deshalb ausserhalb
        # der Schritt-Schleife. Die area_table-Objekte kommen aus den bereits
        # gerechneten Kontexten (lazy properties), damit die ln(Flaeche)-Fits
        # nicht zweimal laufen.
        try:
            run_pko_comparison(
                cells_producers=cells_osc, cells_pko=cells_pko,
                area_producers=osc_ctx.area_table, area_pko=pko_ctx.area_table,
                output_dir=OUTPUT_DIR_PKO,
                pko_biosensor=PKO_BIOSENSOR_NAME,
                analysis_start_min=OSCILLATION_START_MIN,
            )
        except Exception:
            logger.exception("WT-gegen-PKO-Vergleich fehlgeschlagen - wird uebersprungen.")
            failed.append("PKO/wt_vs_pko_comparison")

    # ------------------------------------------------------------------
    # 4. Mit QC vs. ohne QC - nur fuer die Batches, die QC gesehen haben
    # ------------------------------------------------------------------
    # Ein dritter Lauf auf den ROHDATEN (keine Merges, keine Exclusions), auf
    # die QC-beruehrten Batches beschraenkt, in analysis_output/no_qc/. Danach
    # stellt qc_comparison.py Kammer fuer Kammer gegenueber, was das QC an den
    # Kennzahlen geaendert hat. Beschraenkt, weil "mit QC" ueberall sonst mit
    # "ohne QC" identisch ist - der Vergleich saehe dort wie "kein Effekt" aus.
    batches = qc_batches(exclusions)
    if args.skip_qc_comparison:
        logger.info("QC-Vergleich uebersprungen (--skip-qc-comparison).")
    elif batches.empty:
        logger.info("Keine QC-Exclusions vorhanden - kein Vergleich mit/ohne QC.")
    else:
        logger.info("QC-beruehrte Batches (Vergleich mit/ohne QC):\n%s", batches.to_string(index=False))
        keys = set(map(tuple, batches[["biosensor", "osc_type", "osc_freq"]].astype(str).values))
        raw = add_time_column(cells_raw, MIN_PER_FRAME)
        raw = add_experiment_units(raw, STATIC_CHIP_LABELS, static_medium_prefix=STATIC_MEDIUM_PREFIX)
        in_batches = pd.Series(
            list(map(tuple, raw[["biosensor", "osc_type", "osc_freq"]].astype(str).values)),
            index=raw.index,
        ).isin(keys)
        raw = raw[in_batches].copy()
        # Nur der Oszillationszweig ist QC-beruehrt; statische/PKO-Batches im
        # QC-File wuerden hier ebenfalls mitlaufen, mit denselben Kontext-Regeln.
        is_static_raw = raw["osc_type"].astype(str).str.strip().str.lower() == "static"
        raw_osc = raw[~is_static_raw]
        if raw_osc.empty:
            logger.warning("QC-Vergleich: die QC-beruehrten Batches enthalten keine Oszillationsdaten.")
        else:
            out_noqc = OUTPUT_DIR / "no_qc"
            noqc_freq_order = resolve_x_order(raw_osc, "osc_freq", FREQ_ORDER)
            noqc_ctx = PipelineContext(
                cells=raw_osc, output_dir=out_noqc,
                freq_order=noqc_freq_order,
                intensity_cols=intensity_cols, ratio_cols=ratio_cols,
                run_sensor_controls=False,
                run_control_consistency=len(noqc_freq_order) >= 2,
                bud_size_threshold=bud_size_threshold,
            )
            run_steps(noqc_ctx, steps)
            failed += [f"no_qc/{k}" for k in noqc_ctx.failed_steps]
            try:
                qc_exclusion_inventory(exclusions, OUTPUT_DIR / "qc_comparison")
                compare_qc_runs(OUTPUT_DIR, out_noqc, batches, OUTPUT_DIR / "qc_comparison")
            except Exception:
                logger.exception("QC-Vergleich fehlgeschlagen - wird uebersprungen.")
                failed.append("qc_comparison")

    if failed:
        logger.error("=== Mit Fehlern beendet. Fehlgeschlagen: %s ===", ", ".join(failed))
        return 1

    logger.info("=== Fertig! Oszillationsdaten: %s | statische Daten: %s | PKO: %s ===",
                OUTPUT_DIR, OUTPUT_DIR_STATIC, OUTPUT_DIR_PKO)
    return 0


if __name__ == "__main__":
    sys.exit(main())
