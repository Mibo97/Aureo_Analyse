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

# ==============================================================================
# Dieselbe Konfiguration wie in run_analysis.py - bei Änderungen dort auch hier
# nachziehen, sonst vergleichst du gegen einen anderen Datenstand.
# ==============================================================================
DATA_ROOT = Path(r"/prj/microfluidic/ma_mimorde/Data")
OUTPUT_DIR = DATA_ROOT.parent / "analysis_output"
CACHE_PATH = OUTPUT_DIR / "combined_results_cache.parquet"
QC_EXCLUSIONS_PATH = OUTPUT_DIR / "qc_exclusions.csv"
MIN_PER_FRAME = 10.0

# Aktuell in run_analysis.py verwendete Lineage-Parameter - hier zum Testen anpassen.
PARAMS = LineageParams(
    mother_min_frames=20,
    bud_max_frames=7,
    tolerance_px=30.0,
)


def load_prepared_cells() -> pd.DataFrame:
    """Lädt & bereitet die Zelldaten exakt wie run_analysis.py auf (bis vor der Lineage-Klassifikation)."""
    cells = load_all_results(DATA_ROOT, cache_path=CACHE_PATH, force_reload=False)
    cells = compute_ratios(cells)

    exclusions = read_qc_exclusions(QC_EXCLUSIONS_PATH)
    cells = apply_track_merges(cells, exclusions)
    cells = apply_qc_exclusions(cells, exclusions, mode="remove")
    cells = add_time_column(cells, MIN_PER_FRAME)
    return cells


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
    print("Lade Zelldaten (nutzt Parquet-Cache falls vorhanden)...")
    cells = load_prepared_cells()

    print(f"Klassifiziere Mutter/Bud-Events mit: {PARAMS}")
    lineage_events = classify_mother_bud(cells, PARAMS)

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
        "  lineage_events = classify_mother_bud(cells, PARAMS)"
    )
