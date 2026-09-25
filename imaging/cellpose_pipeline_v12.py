"""cellpose_pipeline_v12.py - Segmentierung + Tracking je Film, Tracking VOR dem Filtern.

Nachfolger von ../cellpose_pipeline_v11.py fuer einen Film je Aufruf (Array-Jobs, interaktive Schleifen).
Was sich gegenueber v11 aendert (Begruendung: ../docs/tracking_diagnosis.md, ../docs/tracking_plan.md):

  * Kein Objekt wird vor dem Tracking verworfen. Die roi_filter-Regeln (Rand, Flaeche, Solidity,
    Exzentrizitaet) werden zu Spalten: at_border, below_min_area, above_max_area, low_solidity,
    high_eccentricity. Die Analyse entscheidet, was sie ausschliesst.
  * Rohe Cellpose-Labels je Frame nach labels_<stem>.zarr; Tracking (track_labels.py: Gedaechtnis, Kosten
    statt IoU-Tor, Merge/Split, parent_track_id) laeuft danach auf diesem Stack, ohne Cellpose erneut.
  * Rotation (falls konfiguriert) auf den GANZEN Frame vor der Kammererkennung und auf JEDEN Kanal.
  * Pixelgroesse (um/px) aus den nd2-Metadaten, Cellpose-Version und Einstellungen in jeder Zeile.
  * Formmerkmale fuer die Zelltypen nach Rensink et al. 2026: Umfang, Zirkularitaet, Achsen, und der
    Phasenkontrast in der Maske (phase_mean, phase_std); Skelettlaenge optional (--skeleton).
  * Ausgabe je Film: <results>/Single-Cell-Results_<stem>.csv mit allen v11-Spalten (frame, track_id, roi_id,
    area, centroid_y/x, solidity, eccentricity, <stat>_<kanal>, date, condition, replicate, chamber, filename)
    plus den neuen; merge_results.py baut daraus Combined_Results.csv je Experiment.

Aufruf:
  python cellpose_pipeline_v12.py --config <EXPERIMENT>.yaml --settings pipeline_template_v12.yaml --file FILM.nd2
      [--exp-dir DIR] [--processed-subdir 02_processed_v12] [--results-subdir 03_results_v12]
      [--frames a-b] [--no-track] [--no-overlay] [--skeleton]

  --config   die Experiment-YAML von v11 (Kanaele, Kammererkennung, Messkanaele). Ihre io-Pfade werden
             ignoriert, wenn --exp-dir gegeben ist oder aus dem Film abgeleitet wird (<exp>/01_raw_data/FILM).
  --settings ueberschreibt die Abschnitte segmentation, roi_filter, tracking, measurements (eine Datei fuer
             alle Experimente, damit alle mit derselben Einstellung laufen).
"""
from __future__ import annotations

import argparse
import json
import logging
import re
import sys
import time
from dataclasses import asdict
from pathlib import Path
from typing import Callable, Optional

import numpy as np
import pandas as pd
import yaml
from scipy import ndimage as ndi
from skimage import exposure, filters, measure, restoration, transform
from skimage.segmentation import relabel_sequential

sys.path.insert(0, str(Path(__file__).resolve().parent))
from track_labels import LabelTracker, TrackParams, frame_features, relabel_stack, save_stack, tracking_summary  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger("cellpose_v12")

PIPELINE_VERSION = "12.0"
FLAG_COLS = ["at_border", "below_min_area", "above_max_area", "low_solidity", "high_eccentricity"]
_NEW_MODELS = {"cpsam", "cpsam_v2", "cpdino", "cpdino-vitb"}


# ----------------------------------------------------------------------------------------------
# Konfiguration, Metadaten, Kammer (aus v11 uebernommen)
# ----------------------------------------------------------------------------------------------
def load_config(config_path: Path, settings_path: Optional[Path] = None) -> dict:
    cfg = yaml.safe_load(Path(config_path).read_text())
    for key in ("preprocessing", "segmentation", "channels"):
        if key not in cfg:
            raise ValueError(f"Abschnitt '{key}' fehlt in {config_path}")
    cfg.setdefault("roi_filter", {}); cfg.setdefault("tracking", {}); cfg.setdefault("measurements", {})
    if settings_path is not None:
        over = yaml.safe_load(Path(settings_path).read_text()) or {}
        for key in ("segmentation", "roi_filter", "tracking", "measurements"):
            if key in over:
                cfg[key] = {**cfg.get(key, {}), **over[key]}
        cfg["settings_file"] = str(settings_path)
    return cfg


def extract_metadata_from_filename(stem: str) -> dict:
    m = re.match(r"^(\d{6})_(.+?)_(Rep\d+)_(Cham[A-Za-z0-9]+)", stem)
    if m:
        return {"date": m.group(1), "condition": m.group(2), "replicate": m.group(3), "chamber": m.group(4)}
    logger.warning("Konnte Metadaten nicht aus '%s' extrahieren.", stem)
    return {"date": "Unknown", "condition": "Unknown", "replicate": "Unknown", "chamber": "Unknown"}


def detect_chamber(frame: np.ndarray, cfg: dict) -> Optional[tuple[int, int, int, int]]:
    """Iterative Kammererkennung wie in v11 (Otsu, dann Schwelle anheben, groesste Region, nach innen schrumpfen)."""
    img = exposure.rescale_intensity(frame.astype(float), out_range=(0.0, 1.0))
    img = 1.0 - filters.gaussian(img, sigma=1)
    base = filters.threshold_otsu(img); step = 25 / 255.0
    min_area = cfg["min_area_px"]; shrink = cfg.get("shrink_px", 10)
    for it in range(cfg.get("max_iterations", 100)):
        thresh = base + it * step
        if thresh >= 1.0:
            break
        labeled = measure.label(img > thresh)
        props = [p for p in measure.regionprops(labeled) if p.area >= min_area]
        if props:
            r0, c0, r1, c1 = max(props, key=lambda p: p.area).bbox
            logger.info("Kammer gefunden (Iter %d): y=%d:%d, x=%d:%d", it, r0 + shrink, r1 - shrink, c0 + shrink, c1 - shrink)
            return (r0 + shrink, c0 + shrink, r1 - shrink, c1 - shrink)
    return None


def rolling_ball_subtract(frame: np.ndarray, radius: float, downscale_size: int = 256) -> np.ndarray:
    f = frame.astype(np.float32); h, w = f.shape; longest = max(h, w)
    if downscale_size <= 0 or longest <= downscale_size:
        bg = restoration.rolling_ball(f, radius=radius)
    else:
        scale = downscale_size / longest
        small = transform.resize(f, (max(1, round(h * scale)), max(1, round(w * scale))), preserve_range=True, anti_aliasing=True).astype(np.float32)
        bg = transform.resize(restoration.rolling_ball(small, radius=max(radius * scale, 1.0)), f.shape, preserve_range=True).astype(np.float32)
    return np.clip(f - bg, 0, None)


# ----------------------------------------------------------------------------------------------
# Filme: nd2 oder ein Array (Tests)
# ----------------------------------------------------------------------------------------------
class Nd2Movie:
    def __init__(self, path: Path):
        import nd2
        self.path = Path(path); self._f = nd2.ND2File(str(path))
        self.sizes = dict(self._f.sizes); self.axis_order = tuple(self.sizes.keys())
        self._dask = self._f.to_dask()
        assert tuple(self._dask.shape) == tuple(self.sizes.values()), f"Shape mismatch: {self._dask.shape} vs {self.sizes}"
        self.n_frames = self.sizes.get("T", 1)
        try:
            v = self._f.voxel_size(); self.um_per_px = float(v.x)
        except Exception:
            self.um_per_px = float("nan")

    def get_frame(self, t: int, c: int, z: int = 0) -> np.ndarray:
        idx = []
        for ax in self.axis_order:
            idx.append(t if ax == "T" else c if ax == "C" else z if ax == "Z" else slice(None))
        return np.asarray(self._dask[tuple(idx)].compute())

    def close(self):
        self._f.close()


class ArrayMovie:
    """(T, C, Y, X)-Array als Film - fuer Tests ohne nd2."""
    def __init__(self, arr: np.ndarray, um_per_px: float = 0.0733):
        self.arr = arr; self.n_frames = arr.shape[0]; self.um_per_px = um_per_px
    def get_frame(self, t, c, z=0):
        return self.arr[t, c]
    def close(self):
        pass


# ----------------------------------------------------------------------------------------------
# Cellpose
# ----------------------------------------------------------------------------------------------
def make_segment_fn(seg_cfg: dict) -> tuple[Callable[[np.ndarray], np.ndarray], str]:
    from cellpose import models
    model_type = seg_cfg["model_type"]; use_gpu = seg_cfg.get("use_gpu", True)
    is_new = model_type in _NEW_MODELS
    model = (models.CellposeModel(pretrained_model=model_type, gpu=use_gpu) if is_new
             else models.Cellpose(model_type=model_type, gpu=use_gpu))
    try:
        import cellpose; version = str(getattr(cellpose, "version", getattr(cellpose, "__version__", "?")))
    except Exception:
        version = "?"
    niter = seg_cfg.get("niter") or None
    def segment(frame: np.ndarray) -> np.ndarray:
        norm = exposure.rescale_intensity(frame.astype(float), out_range=(0.0, 1.0))
        if is_new:
            res = model.eval(norm, cellprob_threshold=seg_cfg["cellprob_threshold"], flow_threshold=seg_cfg.get("flow_threshold", 0.4),
                             min_size=seg_cfg["min_size_px"], niter=niter)
        else:
            res = model.eval(norm, diameter=seg_cfg.get("diameter_px"), channels=[0, 0], cellprob_threshold=seg_cfg["cellprob_threshold"],
                             flow_threshold=seg_cfg.get("flow_threshold", 0.4), min_size=seg_cfg["min_size_px"], niter=niter)
        return np.asarray(res[0]).astype(np.int32)
    logger.info("Cellpose %s, Modell %s, flow %.2f, cellprob %.2f, min_size %d, niter %s", version, model_type,
                seg_cfg.get("flow_threshold", 0.4), seg_cfg["cellprob_threshold"], seg_cfg["min_size_px"], niter)
    return segment, version


# ----------------------------------------------------------------------------------------------
# Objekte je Frame: Merkmale, Flags, Intensitaeten
# ----------------------------------------------------------------------------------------------
_STAT_FUNCS = {"mean": ndi.mean, "std": ndi.standard_deviation, "min": ndi.minimum, "max": ndi.maximum, "median": ndi.median}


def frame_objects(labels: np.ndarray, phase: np.ndarray, channels: dict[str, np.ndarray], filt_cfg: dict,
                  intensity_stats: list[str], skeleton: bool = False) -> pd.DataFrame:
    objs = frame_features(labels, skeleton=skeleton)
    if objs.empty:
        return objs
    H, W = labels.shape; bnd = int(filt_cfg.get("exclude_boundary_px", 0))
    objs["at_border"] = ((objs.bbox_r0 <= bnd) | (objs.bbox_c0 <= bnd) | ((H - objs.bbox_r1) <= bnd) | ((W - objs.bbox_c1) <= bnd)).astype(bool)
    objs["below_min_area"] = (objs.area < float(filt_cfg.get("min_area_px", 0))).astype(bool)
    objs["above_max_area"] = (objs.area > float(filt_cfg.get("max_area_px", np.inf))).astype(bool)
    objs["low_solidity"] = ((objs.solidity < float(filt_cfg["min_solidity"])) if filt_cfg.get("min_solidity") else False)
    objs["high_eccentricity"] = ((objs.eccentricity > float(filt_cfg["max_eccentricity"])) if filt_cfg.get("max_eccentricity") else False)
    idx = objs.label.values.astype(int)
    objs["phase_mean"] = ndi.mean(phase.astype(np.float64), labels, idx)
    objs["phase_std"] = ndi.standard_deviation(phase.astype(np.float64), labels, idx)
    for name, img in channels.items():
        img64 = img.astype(np.float64)
        for stat in intensity_stats:
            fn = _STAT_FUNCS.get(stat)
            if fn is None:
                continue
            objs[f"{stat}_{name}"] = fn(img64, labels, idx)
    return objs


# ----------------------------------------------------------------------------------------------
# Ein Film
# ----------------------------------------------------------------------------------------------
def track_params_from_cfg(tr_cfg: dict, skeleton: bool) -> TrackParams:
    keys = {f.name for f in TrackParams.__dataclass_fields__.values()}
    kw = {k: v for k, v in (tr_cfg or {}).items() if k in keys}
    kw["skeleton"] = skeleton
    return TrackParams(**kw)


def process_movie(nd2_path: Path, cfg: dict, processed_dir: Path, results_dir: Path, frames: Optional[list[int]] = None,
                  segment_fn: Optional[Callable] = None, movie=None, track: bool = True, qc_overlay: bool = True,
                  skeleton: bool = False, cellpose_version: str = "?") -> Optional[pd.DataFrame]:
    t0 = time.time(); stem = Path(nd2_path).stem; meta = extract_metadata_from_filename(stem)
    processed_dir.mkdir(parents=True, exist_ok=True); results_dir.mkdir(parents=True, exist_ok=True)
    logger.info("=== %s | %s ===", Path(nd2_path).name, meta["condition"])
    movie = movie or Nd2Movie(nd2_path)
    chan = {ch["name"]: {"index": ch["index"], "measure": ch.get("measure_intensity", False)} for ch in cfg["channels"]}
    phase_idx = chan["phase_contrast"]["index"]
    measure_chans = [(n, c["index"]) for n, c in chan.items() if n != "phase_contrast" and c["measure"]]
    rot = cfg["preprocessing"].get("rotation_angle_deg")
    def full_frame(t, c):
        fr = movie.get_frame(t, c, 0)
        if rot is not None:   # Rotation auf den GANZEN Frame, fuer jeden Kanal - vor dem Crop
            fr = transform.rotate(fr, rot, resize=False, preserve_range=True)
        return fr
    ref = full_frame(0, phase_idx)
    roi = detect_chamber(ref, cfg["preprocessing"]["chamber_detection"])
    if roi is None:
        if cfg["preprocessing"]["chamber_detection"].get("fallback_to_full_image", True):
            logger.warning("%s: Kammer nicht gefunden - ganzes Bild.", stem); roi = (0, 0, ref.shape[0], ref.shape[1])
        else:
            logger.error("%s: Kammer nicht gefunden - uebersprungen.", stem); return None
    r0, c0, r1, c1 = roi; H, W = r1 - r0, c1 - c0
    seg_cfg = cfg["segmentation"]; filt_cfg = cfg.get("roi_filter", {}); meas_cfg = cfg.get("measurements", {})
    rb_radius = meas_cfg.get("rolling_ball_radius_px"); rb_down = meas_cfg.get("rolling_ball_downscale_size", 256)
    stats = list(meas_cfg.get("intensity_stats", ["mean", "std"]))
    if segment_fn is None:
        segment_fn, cellpose_version = make_segment_fn(seg_cfg)
    frames = list(frames) if frames is not None else list(range(movie.n_frames))
    T = len(frames)
    labels_path = processed_dir / f"labels_{stem}.zarr"
    import zarr
    kw = dict(mode="w", shape=(T, H, W), chunks=(1, H, W), dtype=np.uint16)
    try:
        labels_z = zarr.open(str(labels_path), zarr_format=2, **kw)
    except TypeError:
        labels_z = zarr.open(str(labels_path), **kw)
    tracker = LabelTracker(track_params_from_cfg(cfg.get("tracking"), skeleton)) if track else None
    per_frame = []; phase_keep = [] if qc_overlay else None
    t_batch = time.perf_counter()
    for i, t in enumerate(frames):
        t_fr = time.perf_counter()
        phase = full_frame(t, phase_idx)[r0:r1, c0:c1]
        chans = {}
        for name, cidx in measure_chans:
            img = full_frame(t, cidx)[r0:r1, c0:c1]
            chans[name] = rolling_ball_subtract(img, rb_radius, rb_down) if rb_radius is not None else img
        labels = segment_fn(phase).astype(np.int32)
        if labels.max() > 0:
            labels = relabel_sequential(labels)[0].astype(np.int32)   # 1..N, wie Cellpose; sicher gegen Luecken
        if labels.max() > 65535:
            raise ValueError(f"{stem} Frame {t}: mehr als 65535 Objekte")
        labels_z[i] = labels.astype(np.uint16)
        objs = frame_objects(labels, phase, chans, filt_cfg, stats, skeleton=skeleton)
        if tracker is not None:
            objs = tracker.update(i, labels.astype(np.int64), objs if not objs.empty else None)
            objs = objs.drop(columns=["frame"])
        objs.insert(0, "frame", t)
        per_frame.append(objs)
        if phase_keep is not None:
            phase_keep.append(phase)
        el = time.perf_counter() - t_batch
        logger.info("Frame %3d/%d | Objekte %3d | %.1f s | verbleibend ~%.1f min", i + 1, T, len(objs), time.perf_counter() - t_fr, el / (i + 1) * (T - i - 1) / 60)
    movie.close()
    table = pd.concat([p for p in per_frame if not p.empty], ignore_index=True) if any(not p.empty for p in per_frame) else pd.DataFrame()
    if table.empty:
        logger.warning("%s: keine Objekte.", stem); return None
    if tracker is not None:
        events = pd.DataFrame(tracker.events)
        summary = tracking_summary(table.rename(columns={}), events); summary["stem"] = stem
        (events if not events.empty else pd.DataFrame(columns=["frame", "type", "track_id", "other_track_id", "gap_frames"])
         ).to_csv(processed_dir / f"{stem}_events.csv", index=False)
        json.dump(summary, open(processed_dir / f"{stem}_tracking_summary.json", "w"), indent=2)
        # tracks_<stem>.zarr: Stack-Index i, Spalte frame = echter Frame; fuer relabel_stack braucht es den Stack-Index
        tab_idx = table.copy(); tab_idx["frame"] = tab_idx["frame"].map({t: i for i, t in enumerate(frames)})
        save_stack(processed_dir / f"tracks_{stem}.zarr", relabel_stack(labels_z, tab_idx))
        logger.info("%s: %d Spuren, %.3f neue IDs je Objekt-Frame, Median %s Frames", stem, summary["n_tracks"],
                    summary["new_ids_per_object_frame"], summary["median_track_frames"])
    else:
        table["track_id"] = -1
    # Tabelle im v11-Layout plus neue Spalten
    table = table.rename(columns={"label": "roi_id"})
    for k, v in meta.items():
        table[k] = v
    table["filename"] = stem
    table["um_per_px"] = float(getattr(movie, "um_per_px", np.nan))
    table["pipeline_version"] = PIPELINE_VERSION; table["cellpose_version"] = cellpose_version
    table["model_type"] = seg_cfg["model_type"]; table["flow_threshold"] = seg_cfg.get("flow_threshold", 0.4)
    table["cellprob_threshold"] = seg_cfg["cellprob_threshold"]
    lead = ["frame", "track_id", "roi_id", "area", "centroid_y", "centroid_x", "solidity", "eccentricity"]
    table = table[lead + [c for c in table.columns if c not in lead and c != "r"]]
    out_path = results_dir / f"Single-Cell-Results_{stem}.csv"
    table.to_csv(out_path, index=False)
    if qc_overlay and phase_keep is not None:
        try:
            write_overlay(phase_keep, labels_z, table, frames, results_dir / "QC" / f"{stem}_QC_overlay.tif")
        except Exception as e:  # noqa: BLE001
            logger.warning("QC-Overlay fehlgeschlagen: %s", e)
    logger.info("%s: fertig in %.0f s | %d Zeilen, %d Spuren -> %s", stem, time.time() - t0, len(table), table["track_id"].nunique(), out_path.name)
    return table


def write_overlay(phase_frames, labels_z, table: pd.DataFrame, frames: list[int], out_path: Path) -> None:
    """Konturen + Track-IDs auf Phasenkontrast als ImageJ-TIF (wie v11)."""
    import tifffile
    from PIL import Image, ImageDraw
    from skimage.segmentation import find_boundaries
    out_path.parent.mkdir(parents=True, exist_ok=True)
    by_frame = {f: sub for f, sub in table.groupby("frame")}
    stack = []
    for i, t in enumerate(frames):
        img = exposure.rescale_intensity(phase_frames[i].astype(float), out_range=(0, 255)).astype(np.uint8)
        rgb = np.stack([img] * 3, axis=-1)
        lab = np.asarray(labels_z[i]); rgb[find_boundaries(lab, mode="outer")] = (255, 60, 60)
        pil = Image.fromarray(rgb); draw = ImageDraw.Draw(pil)
        for r in by_frame.get(t, pd.DataFrame()).itertuples():
            draw.text((int(r.centroid_x) + 1, int(r.centroid_y) + 1), str(int(r.track_id)), fill=(0, 0, 0))
            draw.text((int(r.centroid_x), int(r.centroid_y)), str(int(r.track_id)), fill=(255, 255, 255))
        stack.append(np.array(pil))
    arr = np.stack(stack)
    try:
        tifffile.imwrite(str(out_path), arr, compression="lzw", imagej=True, metadata={"axes": "TYXS"}, photometric="rgb")
    except Exception:  # LZW braucht imagecodecs; ohne Kompression geht es immer
        tifffile.imwrite(str(out_path), arr, imagej=True, metadata={"axes": "TYXS"}, photometric="rgb")


def parse_frames(s: Optional[str]) -> Optional[list[int]]:
    if not s:
        return None
    a, b = s.split("-"); return list(range(int(a), int(b) + 1))


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, required=True, help="Experiment-YAML (Kanaele, Kammer, Messung)")
    ap.add_argument("--settings", type=Path, default=None, help="globale Einstellungen (segmentation/roi_filter/tracking/measurements)")
    ap.add_argument("--file", type=Path, required=True, help="ein nd2-Film")
    ap.add_argument("--exp-dir", type=Path, default=None, help="Experimentordner; Standard: Elternordner von 01_raw_data")
    ap.add_argument("--processed-subdir", default="02_processed_v12"); ap.add_argument("--results-subdir", default="03_results_v12")
    ap.add_argument("--frames", default=None, help="a-b (inklusive), nur fuer Tests")
    ap.add_argument("--no-track", action="store_true"); ap.add_argument("--no-overlay", action="store_true"); ap.add_argument("--skeleton", action="store_true")
    a = ap.parse_args(argv)
    cfg = load_config(a.config, a.settings)
    exp_dir = a.exp_dir or a.file.resolve().parent.parent
    process_movie(a.file, cfg, exp_dir / a.processed_subdir, exp_dir / a.results_subdir, frames=parse_frames(a.frames),
                  track=not a.no_track, qc_overlay=not a.no_overlay, skeleton=a.skeleton)


if __name__ == "__main__":
    main()
