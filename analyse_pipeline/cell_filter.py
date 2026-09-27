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

KONTRASTREGEL (nur Tabellen mit phase_mean/phase_std, Pipeline v12): tote Zellen und Zelltruemmer
verlieren im Phasenkontrast ihren Kontrast (kein Dichteunterschied zum Medium mehr). Je Spur wird
der Median von phase_std / phase_mean gebildet. Auf dem truemmerreichen Chip BSG/pH/6 ist diese
Groesse zweigipflig (Moden 0.05 und 0.28, Tal 0.11-0.16), auf WT/pH/6 eingipflig bei 0.2-0.3;
Zellen ab 3,000 px2 liegen in beiden Experimenten nie unter 0.18 (5 %%-Quantil). Eine Spur unter
min_phase_cv ist keine Zelle. Mit 0.12 entfernt die Regel auf BSG/pH/6 ZUSAETZLICH zur
Groessenregel 44 Spuren mit 20 %% der Objekt-Frames (runde, schrumpfende Objekte von 1,300-1,700 px2
ueber ~30 Frames, 70 %% davon mit beruehrenden "Knospen"), auf WT/pH/6 4 Spuren (0.1 %%). Ohne die
Spalten oder mit min_phase_cv = None bleibt die Regel aus; der Bericht traegt dann NaN.
"""
from __future__ import annotations

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def track_phase_cv(cells: pd.DataFrame) -> pd.Series | None:
    """Median von phase_std / phase_mean je cell_uid (None ohne die Spalten)."""
    if not {"phase_mean", "phase_std"}.issubset(cells.columns):
        return None
    mean = pd.to_numeric(cells["phase_mean"], errors="coerce")
    std = pd.to_numeric(cells["phase_std"], errors="coerce")
    cv = std / mean.where(mean > 0)
    return cv.groupby(cells["cell_uid"]).median()


def flag_cells(cells: pd.DataFrame, min_frames: int = 2, min_max_area_px: float = 1500.0,
               min_phase_cv: float | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Spalte is_cell an jeder Zeile plus Bericht je exp_id (Spuren/Objekt-Frames entfernt).

    min_phase_cv: Kontrastregel (siehe Modul-Docstring); None oder fehlende phase_*-Spalten = aus.
    """
    if cells.empty or not {"cell_uid", "frame", "area"}.issubset(cells.columns):
        return cells.assign(is_cell=True), pd.DataFrame()
    per_track = cells.groupby("cell_uid").agg(n_frames=("frame", "nunique"), max_area=("area", "max"))
    ok_size = (per_track["n_frames"] >= min_frames) & (per_track["max_area"] >= min_max_area_px)
    phase_cv = track_phase_cv(cells) if min_phase_cv is not None else None
    contrast_active = phase_cv is not None
    if contrast_active:
        # NaN (phase_mean fehlt oder 0) verwirft keine Spur: nur ein gemessener niedriger Kontrast zaehlt.
        low_contrast = (phase_cv.reindex(per_track.index) < min_phase_cv).fillna(False)
        ok = ok_size & ~low_contrast
    else:
        low_contrast = pd.Series(False, index=per_track.index)
        ok = ok_size
    out = cells.copy()
    out["is_cell"] = out["cell_uid"].map(ok).fillna(False).astype(bool)
    by_contrast = ok_size & low_contrast          # nur die Kontrastregel entfernt sie
    out["_by_contrast"] = out["cell_uid"].map(by_contrast).fillna(False).astype(bool)
    key = "exp_id" if "exp_id" in out.columns else None
    if key:
        rep = out.groupby(key).agg(
            n_tracks=("cell_uid", "nunique"),
            n_object_frames=("cell_uid", "size"),
            n_tracks_removed=("cell_uid", lambda s: s[~out.loc[s.index, "is_cell"]].nunique()),
            n_object_frames_removed=("is_cell", lambda s: int((~s).sum())),
            n_tracks_removed_by_contrast=("cell_uid", lambda s: s[out.loc[s.index, "_by_contrast"]].nunique()),
            n_object_frames_removed_by_contrast=("_by_contrast", lambda s: int(s.sum())),
        ).reset_index()
        rep["share_tracks_removed"] = rep["n_tracks_removed"] / rep["n_tracks"].clip(lower=1)
        rep["share_object_frames_removed"] = rep["n_object_frames_removed"] / rep["n_object_frames"].clip(lower=1)
        rep["share_object_frames_removed_by_contrast"] = (
            rep["n_object_frames_removed_by_contrast"] / rep["n_object_frames"].clip(lower=1))
        rep["min_frames"] = min_frames
        rep["min_max_area_px"] = min_max_area_px
        rep["min_phase_cv"] = min_phase_cv if contrast_active else np.nan
        if not contrast_active:
            for c in ("n_tracks_removed_by_contrast", "n_object_frames_removed_by_contrast",
                      "share_object_frames_removed_by_contrast"):
                rep[c] = np.nan
    else:
        rep = pd.DataFrame()
    out = out.drop(columns="_by_contrast")
    n_tr = int((~ok).sum()); n_of = int((~out["is_cell"]).sum())
    logger.info(
        "Zellfilter (>= %d Frames, groesste Flaeche >= %.0f px2): %d von %d Spuren (%.0f %%) und %d von %d "
        "Objekt-Frames (%.1f %%) sind keine Zellen und werden entfernt (00_cell_filter.csv).",
        min_frames, min_max_area_px, n_tr, len(per_track), 100 * n_tr / max(len(per_track), 1),
        n_of, len(out), 100 * n_of / max(len(out), 1),
    )
    if contrast_active:
        n_ct = int(by_contrast.sum()); n_cf = int((out["cell_uid"].map(by_contrast).fillna(False)).sum())
        logger.info(
            "  davon durch die Kontrastregel (Median phase_std/phase_mean < %.2f, zusaetzlich zur Groessenregel): "
            "%d Spuren, %d Objekt-Frames (%.1f %%) - tote Zellen / Zelltruemmer ohne Phasenkontrast.",
            min_phase_cv, n_ct, n_cf, 100 * n_cf / max(len(out), 1),
        )
    elif min_phase_cv is not None:
        logger.info("  Kontrastregel aus: keine phase_mean/phase_std-Spalten (Tabellen vor Pipeline v12).")
    return out, rep
