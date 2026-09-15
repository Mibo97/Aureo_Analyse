#!/usr/bin/env python3
"""
plot_qc_lineage_overlay.py
============================
Eigenständiges CLI-Skript (unabhängig von run_analysis.py aufrufbar): nimmt
das bereits von der Cellpose-Pipeline erzeugte QC-Overlay-TIFF für eine
Kammer (Segmentierungskonturen bereits eingezeichnet) und ergänzt es um die
Mutter/Bud-Verknüpfungen aus classify_mother_bud() (lineage.py) - Mutter- und
Bud-Zentroid werden markiert und mit einer Linie verbunden, dazu ein
Text-Label mit den cell_uids. Kein Neuzeichnen der Segmentierung nötig, das
bestehende QC-Overlay wird einfach wiederverwendet (siehe Frage im Chat).

VERWENDUNG
----------
--exp-id ist NICHT mehr nötig - der Pfad zum QC-Overlay enthält bereits alle
sechs Bestandteile von exp_id (biosensor/osc_type/osc_freq im Ordnerpfad,
condition/replicate/chamber im Dateinamen, siehe data_loading.py), das
Skript leitet exp_id automatisch daraus ab (derive_exp_id()) und prüft das
Ergebnis gegen cells['exp_id'], bevor annotiert wird.

Variante A - Dateipfad direkt angeben (immer zuverlässig, empfohlen):
    python plot_qc_lineage_overlay.py \
        --cells cells.parquet --lineage-events lineage_events.csv \
        --qc-tif "/Data/StrainX/Glc/1.5/03_results/QC/260616_Osc1.5_NegCtrl_Rep1_ChamA13_QC_overlay.tif" \
        --output ChamA13_lineage_overlay.tif

Variante B - automatisch unter --data-root suchen, über den entscheidenden
Dateinamen-Teilstring (Datum kann variieren):
    python plot_qc_lineage_overlay.py \
        --cells cells.parquet --lineage-events lineage_events.csv \
        --data-root /Data --match "_Osc1.5_NegCtrl_Rep1_ChamA13_"

--exp-id bleibt als Override verfügbar, falls euer Namensschema doch einmal
abweicht oder die automatische Ableitung fehlschlägt.

Standardmäßig werden nur Frames mit einem erkannten Budding-Event
ausgegeben (--mode events - das ist der eigentliche Kalibrierungs-Zweck).
Mit --mode all bzw. --mode range lässt sich das erweitern.

ANNAHMEN (bei Abweichung bitte anpassen)
-----------------------------------------
- 'cells' und 'lineage_events' sind die Tabellen aus eurer Pipeline
  (lineage.classify_mother_bud() mit cell_uid-Spalte, siehe lineage.py) -
  als .parquet oder .csv einlesbar.
- Das QC-Overlay-TIFF ist ein Multi-Page-RGB-Stack, eine Seite pro Frame,
  Seitenindex 0-basiert == 'frame'-Spalte in 'cells' (Standard-Konvention
  eurer restlichen Pipeline, siehe budding_ratio_timeseries.py etc.).
- lineage_events enthält 'mother_cell_uid'/'bud_cell_uid' (nicht nur
  track_id) - diese Spalten existieren nur, wenn classify_mother_bud() mit
  einem Datensatz aufgerufen wurde, der eine 'cell_uid'-Spalte hatte.
"""

from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

MOTHER_COLOR = (0, 200, 140)     # teal
BUD_COLOR = (255, 90, 40)        # orange-rot
LINE_COLOR = (255, 220, 0)       # gelb
MARKER_RADIUS = 10

# Dateinamen-Muster eurer Cellpose-QC-Overlays, siehe Beispiele im Chat:
#   260616_Osc1.5_NegCtrl_Rep1_ChamA13_QC_overlay.tif  (Kontroll-Kammer)
#   260616_Osc1.5_Rep3_ChamA6_QC_overlay.tif           (Oszillations-Kammer, kein _NegCtrl/_PosCtrl-Segment)
# 'condition' ist alles zwischen Datum und '_Rep' - dadurch werden beide
# Fälle abgedeckt, ohne _NegCtrl/_PosCtrl explizit im Muster zu benötigen.
_QC_FILENAME_RE = re.compile(
    r"^\d+_(?P<condition>.+?)_Rep(?P<replicate>\d+)_Cham(?P<chamber>[A-Za-z0-9]+)_QC_overlay\.tif$"
)


# ==============================================================================
# 0. exp_id automatisch aus dem QC-Pfad ableiten (statt redundant per --exp-id)
# ==============================================================================

def parse_qc_filename(name: str) -> tuple[str, str, str]:
    """Extrahiert (condition, replicate, chamber) aus dem QC-Overlay-Dateinamen."""
    m = _QC_FILENAME_RE.match(name)
    if not m:
        raise ValueError(
            f"Dateiname '{name}' passt nicht auf das erwartete Muster "
            f"'<Datum>_<condition>_Rep<N>_Cham<X>_QC_overlay.tif'. "
            f"Bitte --exp-id explizit angeben, wenn euer Namensschema abweicht."
        )
    return m.group("condition"), m.group("replicate"), m.group("chamber")


def find_hierarchy_before_results_dir(
    qc_path: Path, results_dir_name: str = "03_results",
) -> tuple[str, str, str]:
    """
    Liest (biosensor, osc_type, osc_freq) aus den DREI Ordnern direkt VOR
    'results_dir_name' im Pfad - exakt dieselbe Hierarchie-Bedeutung wie in
    data_loading._parse_hierarchy(), aber ankerbasiert (sucht '03_results'
    irgendwo im Pfad) statt an einer festen Tiefe relativ zu data_root. Das
    funktioniert unabhängig davon, wie viele Unterordner (z.B. 'QC') nach
    '03_results' noch folgen, und braucht kein --data-root.
    """
    parts = qc_path.resolve().parts
    if results_dir_name not in parts:
        raise ValueError(
            f"'{results_dir_name}' kommt im Pfad '{qc_path}' nicht vor - kann Biosensor/"
            f"Oszillationstyp/-frequenz nicht ableiten. Bitte --exp-id explizit angeben "
            f"oder --results-dir-name anpassen, falls euer Ordner anders heißt."
        )
    idx = parts.index(results_dir_name)
    if idx < 3:
        raise ValueError(
            f"Vor '{results_dir_name}' stehen im Pfad '{qc_path}' weniger als 3 Ordner "
            f"(biosensor/osc_type/osc_freq erwartet) - bitte --exp-id explizit angeben."
        )
    biosensor, osc_type, osc_freq = parts[idx - 3:idx]
    return biosensor, osc_type, osc_freq


def derive_exp_id(
    qc_path: Path,
    cells: pd.DataFrame,
    results_dir_name: str = "03_results",
) -> str:
    """
    Baut exp_id EXAKT nach der Konvention aus data_loading.load_all_results()
    ("__".join([biosensor, osc_type, osc_freq, condition, replicate, chamber]))
    aus dem QC-Overlay-Pfad zusammen - der Pfad enthält bereits alle sechs
    Bestandteile, --exp-id separat anzugeben ist daher redundant.

    Wird gegen cells['exp_id'] verifiziert, damit eine falsche Ableitung
    NICHT still die falsche Kammer annotiert:
      1. exakter Treffer -> wird verwendet.
      2. kein exakter Treffer, aber GENAU EIN exp_id stimmt in allen Teilen
         außer 'condition' überein (z.B. weil euer condition-String anders
         aufgebaut ist als hier angenommen) -> wird mit Warnung verwendet.
      3. sonst -> klarer Fehler mit Vorschlag, --exp-id explizit zu setzen.
    """
    biosensor, osc_type, osc_freq = find_hierarchy_before_results_dir(qc_path, results_dir_name)
    condition, replicate, chamber = parse_qc_filename(qc_path.name)

    candidate = "__".join([biosensor, osc_type, osc_freq, condition, replicate, chamber])
    known_ids = set(cells["exp_id"].astype(str))

    if candidate in known_ids:
        logger.info("exp_id aus Pfad abgeleitet: '%s'", candidate)
        return candidate

    key_no_condition = (biosensor, osc_type, osc_freq, replicate, chamber)
    fallback_matches = [
        eid for eid in known_ids
        if len((p := eid.split("__"))) == 6 and (p[0], p[1], p[2], p[4], p[5]) == key_no_condition
    ]
    if len(fallback_matches) == 1:
        logger.warning(
            "exp_id '%s' (aus Pfad abgeleitet) nicht exakt in cells['exp_id'] gefunden, aber "
            "genau ein eindeutiger Treffer über Biosensor/Osz.typ/-frequenz/Replikat/Kammer "
            "(condition weicht ab) - verwende '%s'. Falls das falsch ist: --exp-id explizit setzen.",
            candidate, fallback_matches[0],
        )
        return fallback_matches[0]

    example_ids = sorted(known_ids)[:8]
    raise ValueError(
        f"Konnte exp_id nicht sicher aus dem Pfad ableiten (erwartet '{candidate}', "
        f"{'mehrere' if len(fallback_matches) > 1 else 'keine'} Kandidaten über "
        f"Biosensor/Osz.typ/-frequenz/Replikat/Kammer gefunden: {fallback_matches}). "
        f"Bitte --exp-id explizit angeben. Beispiele vorhandener exp_ids: {example_ids}"
    )


# ==============================================================================
# 1. QC-Overlay-Datei finden bzw. laden
# ==============================================================================

def find_qc_overlay(data_root: Path, match: str, chamber: Optional[str] = None) -> Path:
    """
    Sucht rekursiv unter data_root nach *_QC_overlay.tif, deren Dateiname
    'match' enthält (der von euch als entscheidend genannte Teil, z.B.
    '_Osc1.5_NegCtrl_Rep1_' - Datum/Kammer-Kürzel am Anfang/Ende der Datei
    dürfen dabei variieren). Falls 'chamber' gesetzt ist, wird zusätzlich
    'Cham{chamber}' im Namen verlangt, um bei mehreren Kammern derselben
    Bedingung eindeutig zu bleiben.

    Bricht mit einer klaren Fehlermeldung ab, wenn 0 oder >1 Dateien
    gefunden werden - Ratelogik bei Bildpfaden ist gefährlich, deshalb lieber
    laut scheitern und die Kandidaten auflisten, als das falsche Bild zu nehmen.
    """
    candidates = [
        p for p in data_root.rglob("*QC_overlay.tif")
        if match in p.name and (chamber is None or f"Cham{chamber}" in p.name)
    ]
    if len(candidates) == 0:
        raise FileNotFoundError(
            f"Kein QC-Overlay unter '{data_root}' gefunden, das '*{match}*QC_overlay.tif' "
            f"entspricht{f' (mit Cham{chamber})' if chamber else ''}. "
            f"Prüft den --match-Teilstring, oder übergebt den Pfad direkt mit --qc-tif."
        )
    if len(candidates) > 1:
        listing = "\n".join(f"  - {p}" for p in sorted(candidates))
        raise FileNotFoundError(
            f"{len(candidates)} QC-Overlays passen auf '*{match}*QC_overlay.tif' - nicht eindeutig:\n"
            f"{listing}\n"
            f"Bitte --match enger fassen (z.B. Kammer ergänzen) oder --qc-tif direkt angeben."
        )
    logger.info("QC-Overlay gefunden: %s", candidates[0])
    return candidates[0]


def load_stack_as_rgb_frames(path: Path) -> np.ndarray:
    """
    Lädt ein Multi-Page-TIFF und normalisiert es auf die Form
    ``(n_frames, height, width, 3)`` mit dtype ``uint8``.

    Pillow wird bewusst statt ``tifffile`` verwendet: Pillow ist bereits
    eine Abhängigkeit dieses Skripts und kann mehrseitige RGB- und
    Graustufen-TIFFs lesen. Damit bleibt das QC-Werkzeug ohne zusätzliche
    TIFF-spezifische Python-Abhängigkeit ausführbar.
    """
    frames: list[np.ndarray] = []
    try:
        with Image.open(path) as stack:
            for frame_index in range(getattr(stack, "n_frames", 1)):
                stack.seek(frame_index)
                # convert("RGB") behandelt Graustufen-, Palette-, RGBA- und
                # 16-bit-Darstellungen einheitlich für die Overlay-Ausgabe.
                frames.append(np.asarray(stack.convert("RGB"), dtype=np.uint8).copy())
    except (OSError, ValueError) as exc:
        raise ValueError(f"QC-Overlay-TIFF konnte nicht gelesen werden: '{path}' ({exc})") from exc

    if not frames:
        raise ValueError(f"QC-Overlay-TIFF enthält keine Frames: '{path}'.")
    if len({frame.shape for frame in frames}) != 1:
        raise ValueError(
            f"QC-Overlay-TIFF enthält Frames mit unterschiedlichen Bildgrößen: '{path}'."
        )
    return np.stack(frames, axis=0)


# ==============================================================================
# 2. Frame-Auswahl je nach --mode
# ==============================================================================

def select_frames(
    n_frames_in_stack: int,
    events: pd.DataFrame,
    mode: str,
    frame_range: Optional[tuple[int, int]],
    context: int,
) -> list[int]:
    if mode == "all":
        return list(range(n_frames_in_stack))

    if mode == "range":
        if frame_range is None:
            raise ValueError("--mode range erfordert --frame-range START-ENDE (z.B. 0-50).")
        lo, hi = frame_range
        return [f for f in range(lo, hi + 1) if 0 <= f < n_frames_in_stack]

    if mode == "events":
        if events.empty:
            logger.warning("Keine Budding-Events für dieses exp_id - Ausgabe wäre leer, gebe alle Frames zurück.")
            return list(range(n_frames_in_stack))
        base_frames = sorted(events["budding_frame"].unique().tolist())
        frames = set()
        for f in base_frames:
            frames.update(range(max(0, f - context), min(n_frames_in_stack - 1, f + context) + 1))
        return sorted(frames)

    raise ValueError(f"Unbekannter --mode '{mode}' (erlaubt: events, all, range).")


# ==============================================================================
# 3. Marker + Verbindungen einzeichnen
# ==============================================================================

def _centroid(cells: pd.DataFrame, exp_id: str, frame: int, cell_uid: str) -> Optional[tuple[float, float]]:
    row = cells[(cells["exp_id"] == exp_id) & (cells["frame"] == frame) & (cells["cell_uid"] == cell_uid)]
    if row.empty:
        return None
    return float(row["centroid_x"].iloc[0]), float(row["centroid_y"].iloc[0])


def draw_lineage_annotations(
    frame_img: np.ndarray,
    exp_id: str,
    frame: int,
    cells: pd.DataFrame,
    events_at_frame: pd.DataFrame,
    font: ImageFont.ImageFont,
) -> np.ndarray:
    """
    Zeichnet für ALLE Budding-Events, deren budding_frame == frame ist, die
    Mutter- und Bud-Zentroide plus Verbindungslinie und Label auf ein
    einzelnes Frame-Bild. Fehlt der Zentroid einer Zelle in diesem Frame
    (sollte am budding_frame selbst nicht vorkommen, außer bei QC-Lücken),
    wird das Event übersprungen und eine Warnung geloggt statt abzustürzen.
    """
    img = Image.fromarray(frame_img)
    draw = ImageDraw.Draw(img)
    draw.text((6, 4), f"frame={frame}", fill=(255, 255, 255), font=font)

    for _, ev in events_at_frame.iterrows():
        mom_uid, bud_uid = ev["mother_cell_uid"], ev["bud_cell_uid"]
        mom_xy = _centroid(cells, exp_id, frame, mom_uid)
        bud_xy = _centroid(cells, exp_id, frame, bud_uid)
        if mom_xy is None or bud_xy is None:
            logger.warning(
                "Frame %d: Zentroid für Event Mutter='%s'/Bud='%s' fehlt in 'cells' - Event übersprungen.",
                frame, mom_uid, bud_uid,
            )
            continue

        mx, my = mom_xy
        bx, by = bud_xy
        draw.line([(mx, my), (bx, by)], fill=LINE_COLOR, width=2)
        draw.ellipse(
            [mx - MARKER_RADIUS, my - MARKER_RADIUS, mx + MARKER_RADIUS, my + MARKER_RADIUS],
            outline=MOTHER_COLOR, width=2,
        )
        draw.ellipse(
            [bx - MARKER_RADIUS, by - MARKER_RADIUS, bx + MARKER_RADIUS, by + MARKER_RADIUS],
            outline=BUD_COLOR, width=2,
        )
        draw.text((mx + MARKER_RADIUS + 2, my - MARKER_RADIUS - 4), f"M:{mom_uid}", fill=MOTHER_COLOR, font=font)
        draw.text((bx + MARKER_RADIUS + 2, by + MARKER_RADIUS - 8), f"B:{bud_uid}", fill=BUD_COLOR, font=font)

    return np.array(img)


# ==============================================================================
# 4. Hauptablauf
# ==============================================================================

def build_lineage_overlay(
    cells: pd.DataFrame,
    lineage_events: pd.DataFrame,
    exp_id: str,
    qc_tif_path: Path,
    output_path: Path,
    mode: str = "events",
    frame_range: Optional[tuple[int, int]] = None,
    context: int = 0,
) -> Path:
    required_cells = {"exp_id", "frame", "cell_uid", "centroid_x", "centroid_y"}
    missing_cells = required_cells - set(cells.columns)
    if missing_cells:
        raise ValueError(f"'cells' fehlen Spalten: {missing_cells}")

    required_events = {"exp_id", "mother_cell_uid", "bud_cell_uid", "budding_frame"}
    missing_events = required_events - set(lineage_events.columns)
    if missing_events:
        raise ValueError(
            f"'lineage_events' fehlen Spalten: {missing_events}. Stellt sicher, dass "
            f"classify_mother_bud() (lineage.py) mit einem 'cell_uid'-Datensatz lief."
        )

    cells_sub = cells[cells["exp_id"] == exp_id]
    events_sub = lineage_events[lineage_events["exp_id"] == exp_id]
    if cells_sub.empty:
        raise ValueError(f"exp_id '{exp_id}' nicht in 'cells' gefunden.")
    logger.info(
        "exp_id='%s': %d Zellen-Frame-Zeilen, %d Budding-Events in lineage_events.",
        exp_id, len(cells_sub), len(events_sub),
    )

    stack = load_stack_as_rgb_frames(qc_tif_path)
    n_frames_in_stack = stack.shape[0]
    logger.info("QC-Overlay geladen: %d Frames, Bildgröße %dx%d.", n_frames_in_stack, stack.shape[2], stack.shape[1])

    frames_to_render = select_frames(n_frames_in_stack, events_sub, mode, frame_range, context)
    if not frames_to_render:
        raise ValueError("Keine Frames zum Rendern ausgewählt (leere Auswahl) - --mode/--frame-range prüfen.")

    try:
        font = ImageFont.truetype("DejaVuSans.ttf", 12)
    except OSError:
        font = ImageFont.load_default()

    out_frames = []
    for f in frames_to_render:
        if f >= n_frames_in_stack:
            logger.warning("Frame %d liegt außerhalb des QC-Overlay-Stacks (%d Frames) - übersprungen.", f, n_frames_in_stack)
            continue
        events_at_frame = events_sub[events_sub["budding_frame"] == f]
        annotated = draw_lineage_annotations(stack[f], exp_id, f, cells_sub, events_at_frame, font)
        out_frames.append(annotated)

    out_stack = np.stack(out_frames, axis=0)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    # Pillow schreibt einen TIFF-Stack, ohne dass tifffile installiert sein
    # muss. Die einzelnen Seiten sind bereits RGB-uint8-Arrays.
    output_images = [Image.fromarray(frame, mode="RGB") for frame in out_stack]
    output_images[0].save(
        output_path,
        format="TIFF",
        save_all=True,
        append_images=output_images[1:],
        compression="tiff_deflate",
    )
    logger.info(
        "Fertig: %d Frames mit Mutter/Bud-Markern (%d Events insgesamt markiert) -> %s",
        len(out_frames), len(events_sub), output_path,
    )
    return output_path


# ==============================================================================
# 5. CLI
# ==============================================================================

def _read_table(path: Path) -> pd.DataFrame:
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path)


def _parse_frame_range(s: Optional[str]) -> Optional[tuple[int, int]]:
    if s is None:
        return None
    lo, hi = s.split("-")
    return int(lo), int(hi)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Zeichnet Mutter/Bud-Verknüpfungen (lineage.classify_mother_bud) auf ein "
                    "bestehendes Cellpose-QC-Overlay-TIFF für eine Kammer.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--cells", required=True, type=Path, help="Pfad zu cells (.parquet oder .csv)")
    parser.add_argument("--lineage-events", required=True, type=Path, help="Pfad zu lineage_events (.parquet oder .csv)")
    parser.add_argument(
        "--exp-id", default=None,
        help="exp_id der Kammer. Optional - wird per Default automatisch aus dem QC-Pfad "
             "abgeleitet (biosensor/osc_type/osc_freq aus dem Ordnerpfad, condition/replicate/"
             "chamber aus dem Dateinamen, siehe derive_exp_id()). Nur nötig, wenn euer "
             "Namensschema von der Konvention aus data_loading.py abweicht.",
    )
    parser.add_argument(
        "--results-dir-name", default="03_results",
        help="Name des Ergebnis-Ordners, der Biosensor/Osz.typ/-frequenz von condition/"
             "replicate/chamber trennt (Default '03_results', siehe data_loading.py)",
    )

    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--qc-tif", type=Path, help="Direkter Pfad zum bestehenden QC-Overlay-TIFF (immer zuverlässig)")
    src.add_argument("--data-root", type=Path, help="Wurzelverzeichnis, unter dem rekursiv gesucht wird (mit --match)")

    parser.add_argument(
        "--match", default=None,
        help="Entscheidender Dateinamen-Teilstring für die Suche unter --data-root, "
             "z.B. '_Osc1.5_NegCtrl_Rep1_' (Datum/Kammer-Kürzel dürfen variieren). "
             "Erforderlich, wenn --data-root statt --qc-tif verwendet wird.",
    )
    parser.add_argument("--chamber", default=None, help="Optional: engt die Suche zusätzlich auf 'Cham{chamber}' ein")

    parser.add_argument("--mode", choices=["events", "all", "range"], default="events",
                         help="events (Default) = nur Frames mit Budding-Event, all = alle Frames, "
                              "range = --frame-range")
    parser.add_argument("--frame-range", default=None, help="Nur bei --mode range: z.B. '0-50' (inklusive)")
    parser.add_argument("--context", type=int, default=0,
                         help="Nur bei --mode events: zusätzlich N Frames vor/nach jedem Event mit ausgeben (Default 0)")

    parser.add_argument("--output", type=Path, default=None,
                         help="Ausgabe-TIFF-Pfad (Default: ./{exp_id}_lineage_overlay.tif)")

    args = parser.parse_args()

    if args.data_root is not None and args.match is None:
        parser.error("--match ist erforderlich, wenn --data-root statt --qc-tif verwendet wird.")

    cells = _read_table(args.cells)
    lineage_events = _read_table(args.lineage_events)

    qc_tif_path = args.qc_tif or find_qc_overlay(args.data_root, args.match, args.chamber)

    exp_id = args.exp_id or derive_exp_id(qc_tif_path, cells, args.results_dir_name)

    output_path = args.output or Path(f"{exp_id}_lineage_overlay.tif")
    frame_range = _parse_frame_range(args.frame_range)

    build_lineage_overlay(
        cells, lineage_events, exp_id, qc_tif_path, output_path,
        mode=args.mode, frame_range=frame_range, context=args.context,
    )


if __name__ == "__main__":
    main()
