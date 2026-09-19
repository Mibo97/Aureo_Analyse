"""
inspect_runs.py
===============
Welche Kammer stammt aus welchem Chip-Lauf? Vorlage fuer run_map.csv erzeugen.

WARUM
-----
Die Spalte 'replicate' (Rep1, Rep2, ...) ist ein Index ueber die ausgewerteten
Kammern, KEIN biologisches Replikat: mehrere Reps koennen vom selben Chip am
selben Tag stammen. Die biologische Einheit ist der Chip-Lauf, und die einzige
Stelle, an der er in den Daten steckt, ist das Datum am Anfang von 'filename'
(z.B. '260616_Osc1.5_NegCtrl_Rep1_ChamA13' -> 260616).

Zwei Laeufe koennen aber am selben Tag stattfinden. Das Datum allein reicht
also nicht - es braucht eine 'run'-Spalte, die im Zweifel aus dem Laborbuch
kommt. Dieses Skript schreibt dafuer eine VORLAGE mit vorbelegtem Datum und
markiert die Faelle, in denen ein Datum sicher mehr als einen Lauf enthaelt:
dieselbe Kammerposition (z.B. ChamA13) zweimal am selben Tag in derselben
Bedingung kann nicht ein Lauf sein.

    python inspect_runs.py                # Pfade aus config.py
    -> analysis_output/run_map_template.csv

Danach: die Spalte 'run' pruefen bzw. ergaenzen (z.B. '260616a', '260616b'
bei zwei Laeufen an einem Tag), die Datei als analysis_output/run_map.csv
speichern. run_analysis.py liest sie ein und aggregiert dann auf Chip-Laeufe
statt auf den Kammer-Index.
"""

from __future__ import annotations

import logging
import re
import sys

import pandas as pd

from config import OUTPUT_DIR, CACHE_PATH

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

KEY = ["biosensor", "osc_type", "osc_freq", "condition", "replicate", "chamber"]
DATE_RE = re.compile(r"^(\d{6,8})[_-]")


def main() -> int:
    positions = OUTPUT_DIR / "00_cell_positions.parquet"
    source = positions if positions.exists() else CACHE_PATH
    if not source.exists():
        logger.error("Weder %s noch %s gefunden - erst run_analysis.py laufen lassen.", positions, CACHE_PATH)
        return 1
    logger.info("Lese %s", source)
    cells = pd.read_parquet(source)

    if "filename" not in cells.columns:
        logger.error(
            "Die Spalte 'filename' fehlt in den Daten. Dann steckt der Chip-Lauf NIRGENDS in den "
            "Tabellen, und run_map.csv muss vollstaendig von Hand aus dem Laborbuch kommen. "
            "Die Vorlage wird trotzdem geschrieben - ohne Datum."
        )
        cells["filename"] = ""

    missing = [c for c in KEY if c not in cells.columns]
    if missing:
        logger.error("Spalten fehlen: %s", missing)
        return 1

    per_chamber = (
        cells.groupby(KEY, dropna=False)["filename"]
        .agg(lambda s: " | ".join(sorted(str(v) for v in s.dropna().unique())))
        .reset_index()
    )
    per_chamber["date"] = per_chamber["filename"].map(
        lambda f: (DATE_RE.match(str(f)).group(1) if DATE_RE.match(str(f)) else "")
    )
    per_chamber["run"] = per_chamber["date"]  # Vorbelegung: ein Lauf pro Tag - bitte korrigieren

    # Markieren, wo ein Datum SICHER mehr als einen Lauf enthaelt.
    per_chamber["condition_type"] = per_chamber["condition"].astype(str).str.extract(
        r"_(PosCtrl|NegCtrl)$"
    )[0].fillna("Oscillation")
    batch = ["biosensor", "osc_type", "osc_freq", "date", "condition_type", "chamber"]
    dup = per_chamber.groupby(batch)["replicate"].transform("count")
    per_chamber["flag"] = ""
    per_chamber.loc[(dup > 1) & (per_chamber["date"] != ""), "flag"] = (
        "same chamber position twice on this date -> at least 2 runs, split 'run' by hand"
    )
    per_chamber.loc[per_chamber["date"] == "", "flag"] = "no date in filename -> fill 'run' by hand"

    out = OUTPUT_DIR / "run_map_template.csv"
    cols = KEY + ["filename", "date", "run", "flag"]
    per_chamber[cols].sort_values(KEY, key=lambda s: s.astype(str)).to_csv(out, index=False)

    n_dates = per_chamber.loc[per_chamber["date"] != "", "date"].nunique()
    n_flag = int((per_chamber["flag"] != "").sum())
    logger.info("Vorlage geschrieben: %s (%d Kammern, %d verschiedene Daten, %d markiert)",
                out, len(per_chamber), n_dates, n_flag)
    logger.info("Beispiele fuer 'filename':\n  %s",
                "\n  ".join(per_chamber["filename"].head(5).astype(str)))
    logger.info(
        "Kammern pro Datum und Bedingungs-Batch (Median / Max): %.0f / %d",
        per_chamber.groupby(["biosensor", "osc_type", "osc_freq", "date"]).size().median(),
        per_chamber.groupby(["biosensor", "osc_type", "osc_freq", "date"]).size().max(),
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
