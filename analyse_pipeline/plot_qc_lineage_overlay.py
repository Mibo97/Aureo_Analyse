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

Beide Eingabetabellen erzeugt run_analysis.py in OUTPUT_DIR:
    --cells           analysis_output/00_cell_positions.parquet
    --lineage-events  analysis_output/20_budding_events.csv

Variante A - Dateipfad direkt angeben (immer zuverlässig, empfohlen):
    python plot_qc_lineage_overlay.py \
        --cells   ../analysis_output/00_cell_positions.parquet \
        --lineage-events ../analysis_output/20_budding_events.csv \
        --qc-tif "/Data/BSG/Glc/1.5/03_results/QC/260616_Osc1.5_NegCtrl_Rep1_ChamA13_QC_overlay.tif" \
        --output ChamA13_lineage_overlay.tif

Variante B - automatisch unter --data-root suchen, über den entscheidenden
Dateinamen-Teilstring (Datum kann variieren):
    python plot_qc_lineage_overlay.py \
        --cells   ../analysis_output/00_cell_positions.parquet \
        --lineage-events ../analysis_output/20_budding_events.csv \
        --data-root /Data --match "_Osc1.5_NegCtrl_Rep1_ChamA13_"

--exp-id bleibt als Override verfügbar, falls euer Namensschema doch einmal
abweicht oder die automatische Ableitung fehlschlägt.

Standardmäßig werden die Frames ausgegeben, in denen ein Budding-Event
erkannt wurde ODER ein Bud-Kandidat KEINER Mutter zugeordnet werden konnte
(--mode events - das ist der eigentliche Kalibrierungs-Zweck). Mit --mode all
bzw. --mode range lässt sich das erweitern.

LEGENDE DER MARKER
------------------
    tuerkiser Kreis (klein)  Mutter-Zentroid
    tuerkiser Kreis (gross)  adaptiver Suchradius dieser Mutter
    oranger Kreis            zugeordneter Bud, Label 'd/r' = Distanz / Suchradius
                             (nahe 1.0 = gerade eben noch akzeptiert -> pruefen)
    gelbe Linie              Mutter-Bud-Zuordnung
    hellblauer Kreis '?'     Bud-Kandidat OHNE Mutter - faellt aus Budding Ratio,
                             µ_event und Stammbaum heraus

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
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

MOTHER_COLOR = (0, 200, 140)     # teal
BUD_COLOR = (255, 90, 40)        # orange-rot
LINE_COLOR = (255, 220, 0)       # gelb
MISS_COLOR = (170, 170, 255)     # hellblau: Bud-Kandidat OHNE zugeordnete Mutter
MARKER_RADIUS = 10


def short_uid(cell_uid: str) -> str:
    """'BSG__Glc__1.5__Osc1.5_NegCtrl__Rep1__ChamA13__track7' -> 'track7'.

    Innerhalb EINES Bildes ist die exp_id konstant. Der volle cell_uid ist
    schnell >50 Zeichen lang, lief über den Bildrand hinaus und überlagerte
    bei mehreren Events im selben Frame alle anderen Labels unlesbar.
    """
    return str(cell_uid).rsplit("__", 1)[-1]

# Dateinamen-Muster eurer Cellpose-QC-Overlays, siehe Beispiele im Chat:
#   260616_Osc1.5_NegCtrl_Rep1_ChamA13_QC_overlay.tif  (Kontroll-Kammer)
#   260616_Osc1.5_Rep3_ChamA6_QC_overlay.tif           (Oszillations-Kammer, kein _NegCtrl/_PosCtrl-Segment)
# 'condition' ist alles zwischen Datum und '_Rep' - dadurch werden beide
# Fälle abgedeckt, ohne _NegCtrl/_PosCtrl explizit im Muster zu benötigen.
_QC_FILENAME_RE = re.compile(
    r"^\d+_(?P<condition>.+?)_Rep(?P<replicate>\d+)_Cham(?P<chamber>[A-Za-z0-9]+)_QC_overlay\.tiff?$",
    re.IGNORECASE,
)

# Alles vor '_QC_overlay.tif[f]' ist der Datei-Stem, den die Bildverarbeitung
# auch in die Spalte 'filename' von Combined_Results schreibt - das ist die
# zuverlaessigste Bruecke zwischen Bild und Tabelle (siehe derive_exp_id()).
_QC_STEM_RE = re.compile(r"^(?P<stem>.+?)_QC_overlay\.tiff?$", re.IGNORECASE)


def qc_stem(name: str) -> str:
    """'..._ChamA13_QC_overlay.tif' -> '..._ChamA13' (der Stem aus 'filename')."""
    m = _QC_STEM_RE.match(name)
    if not m:
        raise ValueError(
            f"Dateiname '{name}' endet nicht auf '_QC_overlay.tif' / '_QC_overlay.tiff' - "
            f"bitte --exp-id explizit angeben."
        )
    return m.group("stem")


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
    stem = qc_stem(qc_path.name)
    known_ids = sorted({str(e) for e in cells["exp_id"].dropna().unique()})
    if not known_ids:
        raise ValueError("'cells' enthält keine exp_id-Werte.")

    # --- Weg 1: über die Spalte 'filename' (exakt, konventionsunabhängig).
    # Die Bildverarbeitung schreibt denselben Datei-Stem in 'filename' und in
    # den QC-Overlay-Namen - damit ist die Zuordnung eindeutig, egal wie
    # condition/replicate/chamber im Detail geschrieben sind.
    if "filename" in cells.columns:
        hits = sorted({
            str(e) for e in cells.loc[cells["filename"].astype(str) == stem, "exp_id"].dropna().unique()
        })
        if len(hits) == 1:
            logger.info("exp_id über cells['filename'] == '%s' bestimmt: '%s'", stem, hits[0])
            return hits[0]
        if len(hits) > 1:
            raise ValueError(
                f"Datei-Stem '{stem}' gehört zu mehreren exp_ids {hits} - bitte --exp-id explizit angeben."
            )
        logger.info(
            "Kein cells['filename'] == '%s' - versuche Zuordnung über Pfad-Hierarchie und "
            "Dateinamen-Bestandteile.", stem,
        )

    # --- Weg 2: Hierarchie aus dem Pfad + Teilstring-Abgleich der restlichen
    # exp_id-Bestandteile gegen den Dateinamen. Bewusst KEIN Zusammenbauen
    # eines Kandidaten-Strings mehr: die frühere Version zerlegte 'Rep1' zu
    # '1' und 'ChamA13' zu 'A13' und konnte deshalb NIE auf eine echte exp_id
    # treffen (die 'Rep1'/'ChamA13' enthält).
    biosensor, osc_type, osc_freq = find_hierarchy_before_results_dir(qc_path, results_dir_name)
    same_hierarchy = [
        eid for eid in known_ids if eid.split("__")[:3] == [biosensor, osc_type, osc_freq]
    ]
    if not same_hierarchy:
        raise ValueError(
            f"Keine exp_id mit Hierarchie '{biosensor}__{osc_type}__{osc_freq}' in 'cells'. "
            f"Passt --results-dir-name ('{results_dir_name}') zum Pfad '{qc_path}'? "
            f"Beispiele vorhandener exp_ids: {known_ids[:8]}"
        )

    def _match_score(eid: str) -> int:
        """Wie viele der restlichen exp_id-Bestandteile kommen im Dateinamen vor?"""
        return sum(1 for part in eid.split("__")[3:] if part and part in stem)

    scored = [(_match_score(eid), eid) for eid in same_hierarchy]
    best = max(s for s, _ in scored)
    winners = [eid for s, eid in scored if s == best]

    if best > 0 and len(winners) == 1:
        logger.info(
            "exp_id über Pfad-Hierarchie + Dateiname bestimmt: '%s' (%d passende Bestandteile).",
            winners[0], best,
        )
        return winners[0]

    raise ValueError(
        f"Konnte exp_id nicht eindeutig aus '{qc_path.name}' ableiten "
        f"({len(winners)} gleich gute Kandidaten bei {best} passenden Bestandteilen: {winners[:5]}). "
        f"Bitte --exp-id explizit angeben. Kandidaten dieser Hierarchie: {same_hierarchy[:8]}"
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
        p for p in data_root.rglob("*QC_overlay.tif*")
        if p.suffix.lower() in (".tif", ".tiff")
        and match in p.name and (chamber is None or f"Cham{chamber}" in p.name)
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
    miss_frames: Sequence[int] = (),
) -> list[int]:
    if mode == "all":
        return list(range(n_frames_in_stack))

    if mode == "range":
        if frame_range is None:
            raise ValueError("--mode range erfordert --frame-range START-ENDE (z.B. 0-50).")
        lo, hi = frame_range
        return [f for f in range(lo, hi + 1) if 0 <= f < n_frames_in_stack]

    if mode == "events":
        # Nicht zugeordnete Bud-Kandidaten sind genauso Teil der Kalibrierung
        # wie die erkannten Events - ihre Frames kommen daher mit in die Auswahl.
        base_frames = sorted(set(events["budding_frame"].tolist() if not events.empty else []) | set(miss_frames))
        if not base_frames:
            logger.warning(
                "Weder Budding-Events noch nicht zugeordnete Kandidaten für dieses exp_id - "
                "gebe alle Frames zurück."
            )
            return list(range(n_frames_in_stack))
        frames = set()
        for f in base_frames:
            frames.update(range(max(0, f - context), min(n_frames_in_stack - 1, f + context) + 1))
        return sorted(frames)

    raise ValueError(f"Unbekannter --mode '{mode}' (erlaubt: events, all, range).")


# ==============================================================================
# 3. Marker + Verbindungen einzeichnen
# ==============================================================================

def find_unassigned_bud_candidates(
    cells_sub: pd.DataFrame, events_sub: pd.DataFrame,
) -> dict[str, int]:
    """cell_uid -> erster Frame, für alle Bud-KANDIDATEN OHNE zugeordnete Mutter.

    Spiegelt die Kandidaten-Definition aus lineage.classify_mother_bud():
    jede Zelle, die NACH dem ersten Frame der Kammer neu auftaucht, ist ein
    Kandidat. Wer davon in lineage_events nicht als 'bud_cell_uid' vorkommt,
    wurde keiner Mutter zugeordnet und fehlt damit in Budding Ratio, µ_event
    und im Stammbaum.

    classify_mother_bud() zählt diese Fälle nur als eine Gesamtzahl über ALLE
    Kammern ins Log. Für die Kalibrierung ist aber genau interessant, WO im
    Bild sie auftreten: ein Kandidat mitten im Feld neben einer klar
    erkennbaren Mutter bedeutet 'tolerance_px zu klein', einer am Bildrand
    oder auf Debris bedeutet 'ist gar kein Bud'.
    """
    if cells_sub.empty:
        return {}
    first_frame = cells_sub.groupby("cell_uid")["frame"].min()
    chamber_start = cells_sub["frame"].min()
    candidates = first_frame[first_frame > chamber_start]
    assigned = (
        set(events_sub["bud_cell_uid"].astype(str))
        if not events_sub.empty and "bud_cell_uid" in events_sub.columns else set()
    )
    return {str(uid): int(f) for uid, f in candidates.items() if str(uid) not in assigned}


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
    misses_at_frame: Optional[Sequence[str]] = None,
) -> tuple[np.ndarray, int]:
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

    n_drawn = 0
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

        # Suchradius der Mutter mitzeichnen: so ist im Bild direkt sichtbar,
        # ob der Bud knapp oder komfortabel innerhalb der Toleranz lag - und
        # ob der Radius fremde Nachbarzellen mit einschließt.
        radius = ev.get("adaptive_radius_px")
        if radius is not None and pd.notna(radius):
            draw.ellipse([mx - radius, my - radius, mx + radius, my + radius],
                         outline=MOTHER_COLOR)

        draw.line([(mx, my), (bx, by)], fill=LINE_COLOR, width=2)
        draw.ellipse(
            [mx - MARKER_RADIUS, my - MARKER_RADIUS, mx + MARKER_RADIUS, my + MARKER_RADIUS],
            outline=MOTHER_COLOR, width=2,
        )
        draw.ellipse(
            [bx - MARKER_RADIUS, by - MARKER_RADIUS, bx + MARKER_RADIUS, by + MARKER_RADIUS],
            outline=BUD_COLOR, width=2,
        )
        draw.text((mx + MARKER_RADIUS + 2, my - MARKER_RADIUS - 4),
                  f"M:{short_uid(mom_uid)}", fill=MOTHER_COLOR, font=font)

        # d/r = Distanz relativ zum Suchradius: 0 = Bud sitzt auf dem
        # Mutter-Zentroid, ~1 = gerade eben noch akzeptiert. Werte dicht an 1
        # sind die Zuordnungen, die man beim Kalibrieren prüfen will.
        bud_label = f"B:{short_uid(bud_uid)}"
        dist = ev.get("distance_px")
        if dist is not None and pd.notna(dist) and radius is not None and pd.notna(radius) and radius > 0:
            bud_label += f" d/r={dist / radius:.2f}"
        draw.text((bx + MARKER_RADIUS + 2, by + MARKER_RADIUS - 8), bud_label, fill=BUD_COLOR, font=font)
        n_drawn += 1

    for miss_uid in misses_at_frame or []:
        xy = _centroid(cells, exp_id, frame, miss_uid)
        if xy is None:
            continue
        cx, cy = xy
        draw.ellipse(
            [cx - MARKER_RADIUS, cy - MARKER_RADIUS, cx + MARKER_RADIUS, cy + MARKER_RADIUS],
            outline=MISS_COLOR, width=2,
        )
        draw.text((cx + MARKER_RADIUS + 2, cy - MARKER_RADIUS - 4),
                  f"?:{short_uid(miss_uid)}", fill=MISS_COLOR, font=font)

    return np.array(img), n_drawn


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

    required_events = ["exp_id", "mother_cell_uid", "bud_cell_uid", "budding_frame"]
    if lineage_events.empty:
        # Kein einziges Event im gesamten Datensatz: weitermachen statt
        # abbrechen - das Overlay zeigt dann nur die nicht zugeordneten
        # Bud-Kandidaten, und genau das ist hier die gesuchte Information.
        logger.warning(
            "'lineage_events' enthaelt keine Zeilen - es werden nur nicht zugeordnete "
            "Bud-Kandidaten markiert. Das ist der erwartete Ablauf, wenn die Heuristik "
            "gar nichts erkannt hat (z.B. tolerance_px zu klein)."
        )
        lineage_events = pd.DataFrame(columns=required_events)
    missing_events = set(required_events) - set(lineage_events.columns)
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

    # Annahme dieses Skripts: Seitenindex im TIFF == Spalte 'frame'. Stimmt das
    # nicht (z.B. weil 'frame' bei 1 beginnt), landen ALLE Marker eine Seite
    # daneben - das faellt beim Draufschauen kaum auf, macht die Kalibrierung
    # aber wertlos. Deshalb hier laut pruefen statt still annehmen.
    frame_min, frame_max = int(cells_sub["frame"].min()), int(cells_sub["frame"].max())
    if frame_min != 0 or frame_max != n_frames_in_stack - 1:
        logger.warning(
            "'frame' laeuft in dieser Kammer von %d bis %d, der TIFF-Stack hat %d Seiten "
            "(erwartet: 0 bis %d). Falls die Zaehlung versetzt ist, sitzen alle Marker auf "
            "der falschen Seite - bitte an einem Event mit sichtbarem Knospungsereignis pruefen.",
            frame_min, frame_max, n_frames_in_stack, n_frames_in_stack - 1,
        )

    misses = find_unassigned_bud_candidates(cells_sub, events_sub)
    if misses:
        logger.info(
            "%d von %d Bud-Kandidaten dieser Kammer wurden KEINER Mutter zugeordnet - "
            "sie werden hellblau mit '?' markiert.",
            len(misses), len(misses) + len(events_sub),
        )
    misses_by_frame: dict[int, list[str]] = {}
    for uid, f in misses.items():
        misses_by_frame.setdefault(f, []).append(uid)

    frames_to_render = select_frames(
        n_frames_in_stack, events_sub, mode, frame_range, context,
        miss_frames=sorted(misses_by_frame),
    )
    if not frames_to_render:
        raise ValueError("Keine Frames zum Rendern ausgewählt (leere Auswahl) - --mode/--frame-range prüfen.")

    try:
        font = ImageFont.truetype("DejaVuSans.ttf", 12)
    except OSError:
        font = ImageFont.load_default()

    out_frames = []
    n_events_drawn = 0
    for f in frames_to_render:
        if f >= n_frames_in_stack:
            logger.warning("Frame %d liegt außerhalb des QC-Overlay-Stacks (%d Frames) - übersprungen.", f, n_frames_in_stack)
            continue
        events_at_frame = events_sub[events_sub["budding_frame"] == f]
        annotated, n_drawn = draw_lineage_annotations(
            stack[f], exp_id, f, cells_sub, events_at_frame, font,
            misses_at_frame=misses_by_frame.get(f, []),
        )
        out_frames.append(annotated)
        n_events_drawn += n_drawn

    if not out_frames:
        raise ValueError(
            f"Kein einziger ausgewählter Frame liegt im QC-Overlay-Stack ({n_frames_in_stack} Seiten) - "
            f"passen 'frame'-Spalte und TIFF-Seiten zusammen?"
        )

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
        "Fertig: %d Frames gerendert, %d von %d Events gezeichnet, %d nicht zugeordnete "
        "Kandidaten markiert -> %s",
        len(out_frames), n_events_drawn, len(events_sub), len(misses), output_path,
    )
    return output_path


# ==============================================================================
# 5. CLI
# ==============================================================================

def _read_table(path: Path) -> pd.DataFrame:
    """Liest .parquet oder .csv; eine LEERE CSV ergibt einen leeren DataFrame.

    classify_mother_bud() liefert eine leere Tabelle, wenn in einem Datensatz
    kein einziges Event erkannt wurde - to_csv() schreibt dann eine Datei ganz
    ohne Kopfzeile, an der pd.read_csv() mit 'No columns to parse from file'
    abbricht. Genau dieser Fall ist aber der wichtigste fuer das QC-Overlay:
    man will ja sehen, WARUM nichts erkannt wurde.
    """
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        logger.warning("'%s' ist leer - wird als Tabelle ohne Zeilen behandelt.", path)
        return pd.DataFrame()


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
