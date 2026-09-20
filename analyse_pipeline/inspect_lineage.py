"""
inspect_lineage.py
===================
Interaktives Hilfsskript zur Kalibrierung der Lineage-Parameter
(mother_min_frames, bud_max_frames, tolerance_px in lineage.LineageParams).

Macht dieselben Lade-/QC-Schritte wie run_analysis.py (nutzt denselben
Parquet-Cache, ist also nach dem ersten Mal schnell), klassifiziert dann
Mutter/Bud-Events und lässt dich mit lineage.inspect_classification() pro
Kammer nachschauen, welche Budding-Events erkannt wurden - zum Abgleich
gegen das QC-Overlay-TIF derselben Kammer (siehe Hinweise unten).

AUSFÜHRUNG
----------
Am bequemsten INTERAKTIV, damit `cells`/`lineage_events` im Speicher bleiben
und du exp_id / Parameter mehrfach durchprobieren kannst, ohne neu zu laden:

    ipython -i inspect_lineage.py
    # oder: python -i inspect_lineage.py
    # oder Zeile für Zeile in eine Jupyter-Zelle kopieren

Als reines Skript (python inspect_lineage.py) geht's auch, druckt dann nur
die Beispiel-Kammer am Ende und beendet sich.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from data_loading import load_all_results
from qc_exclusions import read_qc_exclusions, apply_track_merges, apply_qc_exclusions
from sensors import compute_ratios
from lineage import classify_mother_bud, inspect_classification, LineageParams
from analysis import add_time_column
from bud_size import load_bud_size_threshold
from relink import gap_close_tracks, detect_sparse_window, flag_lineage_window

# ==============================================================================
# Konfiguration: exakt dieselbe wie in run_analysis.py, weil beide aus config.py
# lesen. Vorher standen die Werte hier ein zweites Mal und waren auseinander
# gelaufen (anderer DATA_ROOT, anderer CACHE_PATH, mother_min_frames 20 statt
# 10) - damit wurde hier faktisch ein anderer Datenstand kalibriert als der,
# den die Pipeline auswertet.
# ==============================================================================
from config import (
    DATA_ROOT,
    OUTPUT_DIR,
    CACHE_PATH,
    QC_EXCLUSIONS_PATH,
    MIN_PER_FRAME,
    LINEAGE_PARAMS,
    BUD_MAX_AREA_FRACTION_FALLBACK,
    LINEAGE_SPARSE_MAX_OBJECTS,
    LINEAGE_SPARSE_SMOOTH_FRAMES,
    LINEAGE_SPARSE_MIN_FRAMES,
    RELINK_MAX_GAP_FRAMES,
    RELINK_MAX_DISTANCE_PX,
    RELINK_MAX_AREA_RATIO,
    log_active_configuration,
)

# Zum Ausprobieren anderer Schwellen einfach hier überschreiben, z.B.:
#     PARAMS = LineageParams(mother_min_frames=20, bud_max_frames=7, tolerance_px=30.0)
# Der Default ist bewusst der PRODUKTIVE Wert aus config.py.
PARAMS: LineageParams = LINEAGE_PARAMS


def load_prepared_cells() -> pd.DataFrame:
    """Lädt & bereitet die Zelldaten exakt wie run_analysis.py auf (bis vor der Lineage-Klassifikation)."""
    cells = load_all_results(DATA_ROOT, cache_path=CACHE_PATH, force_reload=False)
    cells = compute_ratios(cells)

    exclusions = read_qc_exclusions(QC_EXCLUSIONS_PATH)
    # Dieselbe Reihenfolge wie run_analysis.py: Ausschluesse vor UND nach den
    # Merges (ein Track, der gemergt und ausgeschlossen ist, verloere sonst
    # seinen Ausschluss, weil der Merge ihn umbenennt).
    cells = apply_qc_exclusions(cells, exclusions, mode="remove")
    cells = apply_track_merges(cells, exclusions)
    cells = apply_qc_exclusions(cells, exclusions, mode="remove")
    # Wie run_analysis.py: Gap Closing, dann nur das Sparse-Phase-Fenster
    # (relink.py) - die Heuristik sieht nichts anderes.
    cells, _, _ = gap_close_tracks(cells, max_gap=RELINK_MAX_GAP_FRAMES,
                                   max_distance_px=RELINK_MAX_DISTANCE_PX, max_area_ratio=RELINK_MAX_AREA_RATIO)
    cells = add_time_column(cells, MIN_PER_FRAME)
    window = detect_sparse_window(cells, max_objects=LINEAGE_SPARSE_MAX_OBJECTS,
                                  smooth_frames=LINEAGE_SPARSE_SMOOTH_FRAMES, min_frames=LINEAGE_SPARSE_MIN_FRAMES)
    cells = flag_lineage_window(cells, window)
    return cells[cells["in_lineage_window"]].copy()


def qc_overlay_path_for(cells: pd.DataFrame, exp_id: str) -> Path:
    """
    Findet den Pfad zum passenden QC-Overlay-TIF für eine exp_id, damit du
    inspect_classification() direkt gegen das Bild abgleichen kannst.

    exp_id (aus data_loading.py) und der QC-Dateiname (aus
    cellpose_pipeline_v9.py, basiert auf dem NB2-Dateinamen) sind NICHT
    derselbe String - die Brücke ist die Spalte 'filename' in cells (der
    ursprüngliche Datei-Stem). QC-Overlays liegen im Bildverarbeitungs-
    output_dir/QC/, NICHT im hiesigen OUTPUT_DIR - Pfad ggf. anpassen.
    """
    rows = cells.loc[cells["exp_id"] == exp_id, "filename"]
    if rows.empty:
        raise ValueError(f"exp_id '{exp_id}' kommt in cells nicht vor.")
    stem = rows.iloc[0]
    return Path(f"{stem}_QC_overlay.tif")  # ggf. Ordnerpfad (output_dir/QC/) voranstellen


if __name__ == "__main__":
    log_active_configuration()
    print("Lade Zelldaten (nutzt Parquet-Cache falls vorhanden)...")
    cells = load_prepared_cells()

    # Groessenkriterium: dieselbe Schwelle wie im letzten Pipeline-Lauf
    # (20_bud_size_threshold.csv); ohne die Datei der Rueckfallwert, damit
    # hier nicht eine andere Heuristik kalibriert wird als die produktive.
    BUD_SIZE_THRESHOLD = load_bud_size_threshold(OUTPUT_DIR)
    if BUD_SIZE_THRESHOLD is None:
        BUD_SIZE_THRESHOLD = (float("inf") if BUD_MAX_AREA_FRACTION_FALLBACK is None
                              else float(BUD_MAX_AREA_FRACTION_FALLBACK))
        print(f"Keine abgeleitete Groessenschwelle gefunden - Rueckfall {BUD_SIZE_THRESHOLD} (inf = kein Filter).")
    print(f"Klassifiziere Mutter/Bud-Events mit: {PARAMS}, Groessenschwelle {BUD_SIZE_THRESHOLD} (inf = kein Filter)")
    lineage_events = classify_mother_bud(cells, PARAMS, bud_size_threshold=BUD_SIZE_THRESHOLD)

    exp_ids_with_events = sorted(lineage_events["exp_id"].unique()) if not lineage_events.empty else []
    print(f"\n{len(exp_ids_with_events)} Kammern mit mindestens einem erkannten Budding-Event.")
    print("Erste 10 davon (zum Reinkopieren unten):")
    for e in exp_ids_with_events[:10]:
        print(" ", e)

    if exp_ids_with_events:
        example_exp_id = exp_ids_with_events[0]
        print(f"\n--- Beispiel: inspect_classification() für '{example_exp_id}' ---")
        print(inspect_classification(lineage_events, cells, example_exp_id).to_string(index=False))
        print(f"\nZugehöriges QC-Overlay (Dateiname, Pfad ggf. anpassen): "
              f"{qc_overlay_path_for(cells, example_exp_id).name}")

    print(
        "\nInteraktiv weitermachen (falls mit -i gestartet):\n"
        "  inspect_classification(lineage_events, cells, 'DEINE_EXP_ID')\n"
        "  qc_overlay_path_for(cells, 'DEINE_EXP_ID')\n"
        "Parameter ändern & neu klassifizieren, ohne neu zu laden:\n"
        "  PARAMS = LineageParams(mother_min_frames=20, bud_max_frames=10, tolerance_px=30.0)\n"
        "  lineage_events = classify_mother_bud(cells, PARAMS, bud_size_threshold=BUD_SIZE_THRESHOLD)"
    )
