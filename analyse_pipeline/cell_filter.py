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
Groesse zweigipflig (Moden 0.05 und 0.28, Tal 0.11-0.16), auf WT/pH/6 eingipflig bei 0.2-0.3.
Die Schwelle ist RELATIV zur Struktur (Aufnahmesitzung: biosensor/osc_type/osc_freq): eine Spur ist
keine Zelle, wenn ihr Kontrast unter min_phase_cv_rel x Median der Spuren liegt, die die Groessen-
regel derselben Struktur bestehen. Warum relativ: die W109-Filme (statisch, Minimalmedium) sind
vierfach dunkler aufgenommen (Zellmittel ~480 statt 1,500-2,300 Zaehlwerte, Hintergrund ~455), dort
liegen ALLE Zellen bei 0.07-0.12 - eine feste Schwelle von 0.12 entfernte 59 %% ihrer Objekt-Frames,
lebende Zellen. Mit 0.45 x Median: BSG/pH/6 dieselben 44 Spuren (20 %% der groessenbestandenen
Objekt-Frames: runde, schrumpfende Objekte von 1,300-1,700 px2 ueber ~30 Frames, 70 %% davon mit
beruehrenden "Knospen"), WT/pH/6 2 Spuren, W65 statisch 13 Spuren (0.3 %%), W109 keine. Ohne die
Spalten oder mit min_phase_cv_rel = None bleibt die Regel aus; der Bericht traegt dann NaN. Die
Referenz je Struktur steht im Bericht (phase_cv_reference, min_phase_cv).
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
               min_phase_cv_rel: float | None = None,
               structure_cols: tuple[str, ...] = ("biosensor", "osc_type", "osc_freq")) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Spalte is_cell an jeder Zeile plus Bericht je exp_id (Spuren/Objekt-Frames entfernt).

    min_phase_cv_rel: Kontrastregel relativ zum Median der groessenbestandenen Spuren derselben Struktur
    (structure_cols; fehlen die Spalten, gilt die ganze Tabelle als eine Struktur). None oder fehlende
    phase_*-Spalten = aus.
    """
    if cells.empty or not {"cell_uid", "frame", "area"}.issubset(cells.columns):
        return cells.assign(is_cell=True), pd.DataFrame()
    per_track = cells.groupby("cell_uid").agg(n_frames=("frame", "nunique"), max_area=("area", "max"))
    ok_size = (per_track["n_frames"] >= min_frames) & (per_track["max_area"] >= min_max_area_px)
    phase_cv = track_phase_cv(cells) if min_phase_cv_rel is not None else None
    contrast_active = phase_cv is not None
    threshold_per_track = pd.Series(np.nan, index=per_track.index)
    reference_per_track = pd.Series(np.nan, index=per_track.index)
    if contrast_active:
        cols = [c for c in structure_cols if c in cells.columns]
        first = cells.drop_duplicates("cell_uid").set_index("cell_uid")
        struct = (first[cols].astype(str).agg("__".join, axis=1) if cols
                  else pd.Series("all", index=first.index)).reindex(per_track.index)
        pcv = phase_cv.reindex(per_track.index)
        # Referenz je Struktur: Median des Kontrasts der Spuren, die die Groessenregel bestehen.
        ref = pcv[ok_size].groupby(struct[ok_size]).median()
        reference_per_track = struct.map(ref)
        threshold_per_track = float(min_phase_cv_rel) * reference_per_track
        # NaN (phase_mean fehlt oder 0, oder keine Referenz) verwirft keine Spur: nur ein gemessener
        # niedriger Kontrast zaehlt.
        low_contrast = (pcv < threshold_per_track).fillna(False)
        ok = ok_size & ~low_contrast
    else:
        low_contrast = pd.Series(False, index=per_track.index)
        ok = ok_size
    out = cells.copy()
    out["is_cell"] = out["cell_uid"].map(ok).fillna(False).astype(bool)
    by_contrast = ok_size & low_contrast          # nur die Kontrastregel entfernt sie
    out["_by_contrast"] = out["cell_uid"].map(by_contrast).fillna(False).astype(bool)
    out["_thr"] = out["cell_uid"].map(threshold_per_track)
    out["_ref"] = out["cell_uid"].map(reference_per_track)
    key = "exp_id" if "exp_id" in out.columns else None
    if key:
        rep = out.groupby(key).agg(
            n_tracks=("cell_uid", "nunique"),
            n_object_frames=("cell_uid", "size"),
            n_tracks_removed=("cell_uid", lambda s: s[~out.loc[s.index, "is_cell"]].nunique()),
            n_object_frames_removed=("is_cell", lambda s: int((~s).sum())),
            n_tracks_removed_by_contrast=("cell_uid", lambda s: s[out.loc[s.index, "_by_contrast"]].nunique()),
            n_object_frames_removed_by_contrast=("_by_contrast", lambda s: int(s.sum())),
            phase_cv_reference=("_ref", "first"),
            min_phase_cv=("_thr", "first"),
        ).reset_index()
        rep["share_tracks_removed"] = rep["n_tracks_removed"] / rep["n_tracks"].clip(lower=1)
        rep["share_object_frames_removed"] = rep["n_object_frames_removed"] / rep["n_object_frames"].clip(lower=1)
        rep["share_object_frames_removed_by_contrast"] = (
            rep["n_object_frames_removed_by_contrast"] / rep["n_object_frames"].clip(lower=1))
        rep["min_frames"] = min_frames
        rep["min_max_area_px"] = min_max_area_px
        rep["min_phase_cv_rel"] = min_phase_cv_rel if contrast_active else np.nan
        if not contrast_active:
            for c in ("n_tracks_removed_by_contrast", "n_object_frames_removed_by_contrast",
                      "share_object_frames_removed_by_contrast", "phase_cv_reference", "min_phase_cv"):
                rep[c] = np.nan
    else:
        rep = pd.DataFrame()
    out = out.drop(columns=["_by_contrast", "_thr", "_ref"])
    n_tr = int((~ok).sum()); n_of = int((~out["is_cell"]).sum())
    logger.info(
        "Zellfilter (>= %d Frames, groesste Flaeche >= %.0f px2): %d von %d Spuren (%.0f %%) und %d von %d "
        "Objekt-Frames (%.1f %%) sind keine Zellen und werden entfernt (00_cell_filter.csv).",
        min_frames, min_max_area_px, n_tr, len(per_track), 100 * n_tr / max(len(per_track), 1),
        n_of, len(out), 100 * n_of / max(len(out), 1),
    )
    if contrast_active:
        n_ct = int(by_contrast.sum()); n_cf = int((out["cell_uid"].map(by_contrast).fillna(False)).sum())
        ref_vals = reference_per_track.dropna()
        logger.info(
            "  davon durch die Kontrastregel (Median phase_std/phase_mean < %.2f x Strukturmedian; Referenzen %.3f-%.3f, "
            "Schwellen %.3f-%.3f): %d Spuren, %d Objekt-Frames (%.1f %%) - tote Zellen / Zelltruemmer ohne Phasenkontrast.",
            min_phase_cv_rel, float(ref_vals.min()) if len(ref_vals) else float("nan"),
            float(ref_vals.max()) if len(ref_vals) else float("nan"),
            float(threshold_per_track.min()) if threshold_per_track.notna().any() else float("nan"),
            float(threshold_per_track.max()) if threshold_per_track.notna().any() else float("nan"),
            n_ct, n_cf, 100 * n_cf / max(len(out), 1),
        )
        if not rep.empty:
            heavy = rep[rep["share_object_frames_removed_by_contrast"] > 0.25]
            if not heavy.empty:
                logger.warning(
                    "  Kontrastregel entfernt in %d Kammer(n) mehr als 25 %% der Objekt-Frames - viele tote Zellen, oder "
                    "die Struktur-Referenz ist selbst von Truemmern dominiert: %s",
                    len(heavy),
                    ", ".join(f"{r.exp_id} ({100 * r.share_object_frames_removed_by_contrast:.0f} %)"
                              for r in heavy.sort_values("share_object_frames_removed_by_contrast", ascending=False)
                              .head(8).itertuples()),
                )
    elif min_phase_cv_rel is not None:
        logger.info("  Kontrastregel aus: keine phase_mean/phase_std-Spalten (Tabellen vor Pipeline v12).")
    return out, rep


