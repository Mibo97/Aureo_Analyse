"""cell_filter.py - was zaehlt als Zelle?

Die Bild-Pipeline liefert jedes segmentierte Objekt ab min_size_px, und der Tracker verbindet, was
er kann. Uebrig bleiben Spuren, die keine Zellen sind: Schmutz und Halo-Stuecke, die im Fluss
vorbeitreiben, Segmentierungsflackern von einem Frame Dauer. Auf dem re-getrackten QC-Batch
(WT/pH/6, 0.0733 um/px) dauern 26 %% der Spuren einen Frame, ihre groesste Flaeche liegt im
Median bei 580 px2 (3 um2), waehrend Spuren ab 5 Frames zu 95 %% ueber 2,100 px2 kommen; eine
Blastokonidie hat >= 20 um2 (3,700 px2). Zwei Regeln auf SPUR-Ebene trennen das sauber:

    is_cell = (Frames der Spur >= min_frames) & (groesste Flaeche der Spur >= min_max_area_px)

Die groesste Flaeche statt der ersten, damit eine Knospe, die klein beginnt und waechst, Zelle
bleibt. Mit min_frames = 2 und min_max_area_px = 1,500 (8 um2) fallen im QC-Batch 33 %% der
Spuren, aber nur 4 %% der Objekt-Frames. Die Regel gilt fuer alle Tabellen (v11 und re-getrackt)
und laeuft NACH dem manuellen QC; 00_cell_filter.csv zeigt je Kammer, was sie entfernt hat.
"""
from __future__ import annotations

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def flag_cells(cells: pd.DataFrame, min_frames: int = 2, min_max_area_px: float = 1500.0) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Spalte is_cell an jeder Zeile plus Bericht je exp_id (Spuren/Objekt-Frames entfernt)."""
    if cells.empty or not {"cell_uid", "frame", "area"}.issubset(cells.columns):
        return cells.assign(is_cell=True), pd.DataFrame()
    per_track = cells.groupby("cell_uid").agg(n_frames=("frame", "nunique"), max_area=("area", "max"))
    ok = (per_track["n_frames"] >= min_frames) & (per_track["max_area"] >= min_max_area_px)
    out = cells.copy()
    out["is_cell"] = out["cell_uid"].map(ok).fillna(False).astype(bool)
    key = "exp_id" if "exp_id" in out.columns else None
    if key:
        rep = out.groupby(key).agg(
            n_tracks=("cell_uid", "nunique"),
            n_object_frames=("cell_uid", "size"),
            n_tracks_removed=("cell_uid", lambda s: s[~out.loc[s.index, "is_cell"]].nunique()),
            n_object_frames_removed=("is_cell", lambda s: int((~s).sum())),
        ).reset_index()
        rep["share_tracks_removed"] = rep["n_tracks_removed"] / rep["n_tracks"].clip(lower=1)
        rep["share_object_frames_removed"] = rep["n_object_frames_removed"] / rep["n_object_frames"].clip(lower=1)
        rep["min_frames"] = min_frames
        rep["min_max_area_px"] = min_max_area_px
    else:
        rep = pd.DataFrame()
    n_tr = int((~ok).sum()); n_of = int((~out["is_cell"]).sum())
    logger.info(
        "Zellfilter (>= %d Frames, groesste Flaeche >= %.0f px2): %d von %d Spuren (%.0f %%) und %d von %d "
        "Objekt-Frames (%.1f %%) sind keine Zellen und werden entfernt (00_cell_filter.csv).",
        min_frames, min_max_area_px, n_tr, len(per_track), 100 * n_tr / max(len(per_track), 1),
        n_of, len(out), 100 * n_of / max(len(out), 1),
    )
    return out, rep
