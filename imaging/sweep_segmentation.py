"""sweep_segmentation.py - Cellpose-Einstellungen auf wenigen aufeinanderfolgenden Frames vergleichen.

Laeuft auf dem Cluster (GPU). Fuer jede Einstellung des Gitters werden die gewaehlten Frames eines nd2-Films
segmentiert (Kammer-Crop wie in der Pipeline) und ohne Hand-Labels bewertet:

  consistency      neue IDs je Objekt-Frame nach track_labels (weniger = die Objekte werden von Frame zu
                   Frame wiedergefunden); dropout = Anteil Verknuepfungen ueber eine Luecke
  merge / split    aus der Maskenueberlappung aufeinanderfolgender Frames
  large_split      Aufspaltungen von Objekten ueber dem Doppelten der Medianflaeche (septierte
                   Schwellzellen sollen EIN Objekt bleiben)
  n_obj_cv, area_cv  Streuung der Objektzahl und der Gesamtmaskenflaeche ueber die Frames
  s_per_frame      Rechenzeit

Dazu je Einstellung Overlays (Konturen auf Phasenkontrast) fuer die Sichtpruefung.

    python sweep_segmentation.py --config CONFIG.yaml --nd2 FILM.nd2 --frames 40-52 --out DIR \
        [--models cpsam,cpsam_v2] [--flow 0.4,0.6,0.8] [--cellprob -1,0,1] [--niter 0,500] [--min-size 100]

Mehrere --nd2/--frames-Paare sind erlaubt (gleiche Reihenfolge). --frames a-b ist inklusive.
"""
from __future__ import annotations
import argparse, itertools, json, sys, time
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent))
from track_labels import LabelTracker, TrackParams, frame_features, overlap_counts, tracking_summary  # noqa: E402


def load_frames(nd2_path: Path, config_path: Path, frames: list[int]):
    """Phasenkontrast-Frames, gecroppt auf die Kammer wie in cellpose_pipeline_v11."""
    import nd2
    from skimage import transform
    from cellpose_pipeline_v11 import MicroscopyPipeline
    pipe = MicroscopyPipeline(str(config_path))
    out = []
    with nd2.ND2File(str(nd2_path)) as f:
        dims = f.sizes; axis_order = tuple(dims.keys()); dask = f.to_dask()
        phase_idx = pipe.channel_map["phase_contrast"]["index"]
        def get(t, c):
            idx = []
            for ax in axis_order:
                idx.append(t if ax == "T" else c if ax == "C" else 0 if ax == "Z" else slice(None))
            return np.asarray(dask[tuple(idx)].compute())
        ref = get(0, phase_idx)
        if pipe.rotation_angle is not None:
            ref = transform.rotate(ref, pipe.rotation_angle, resize=False, preserve_range=True)
        roi = pipe.detect_chamber(ref) or (0, 0, ref.shape[0], ref.shape[1])
        r0, c0, r1, c1 = roi
        for t in frames:
            fr = get(t, phase_idx)
            if pipe.rotation_angle is not None:
                fr = transform.rotate(fr, pipe.rotation_angle, resize=False, preserve_range=True)
            out.append(fr[r0:r1, c0:c1])
        try:
            vox = f.voxel_size(); px = float(vox.x)
        except Exception:
            px = float("nan")
    return out, pipe, px


def segment(model, frame, flow, cellprob, niter, min_size):
    from skimage import exposure
    norm = exposure.rescale_intensity(frame.astype(float), out_range=(0.0, 1.0))
    res = model.eval(norm, flow_threshold=flow, cellprob_threshold=cellprob, min_size=min_size, niter=(niter or None))
    return np.asarray(res[0]).astype(np.int32)


def consistency_metrics(labels: list[np.ndarray]) -> dict:
    tr = LabelTracker(TrackParams(sparse_n=0))
    feats = []
    for t, lab in enumerate(labels):
        f = frame_features(lab); feats.append(f); tr.update(t, lab.astype(np.int64), f)
    tab = tr.table(); ev = pd.DataFrame(tr.events); s = tracking_summary(tab, ev)
    n_merge = n_split = n_large_split = 0
    med = np.median(np.concatenate([f.area.values for f in feats if len(f)])) if any(len(f) for f in feats) else np.nan
    for t in range(1, len(labels)):
        oc = overlap_counts(labels[t - 1].astype(np.int64), labels[t].astype(np.int64))
        prev_a = dict(zip(feats[t - 1].label.astype(int), feats[t - 1].area)); cur_a = dict(zip(feats[t].label.astype(int), feats[t].area))
        by_prev, by_cur = {}, {}
        for (a, b), c in oc.items():
            by_prev.setdefault(a, []).append((b, c)); by_cur.setdefault(b, []).append((a, c))
        for a, lst in by_prev.items():   # ein Vorgaenger -> zwei Nachfolger, beide substanziell
            big = [b for b, c in lst if c >= 0.3 * cur_a[b]]
            if len(big) >= 2:
                n_split += 1; n_large_split += int(prev_a[a] > 2 * med)
        for b, lst in by_cur.items():    # ein Nachfolger deckt zwei Vorgaenger ab
            if sum(c >= 0.5 * prev_a[a] for a, c in lst) >= 2:
                n_merge += 1
    n_obj = np.array([len(f) for f in feats]); area = np.array([f.area.sum() if len(f) else 0 for f in feats])
    n_links = int(len(tab) - (tab.frame == 0).sum())
    return dict(new_ids_per_object_frame=s["new_ids_per_object_frame"], gap_links_share=s["n_gap_links"] / max(n_links, 1),
                merges_per_frame=n_merge / max(len(labels) - 1, 1), splits_per_frame=n_split / max(len(labels) - 1, 1),
                large_splits_per_frame=n_large_split / max(len(labels) - 1, 1),
                n_obj_mean=float(n_obj.mean()), n_obj_cv=float(n_obj.std() / max(n_obj.mean(), 1e-9)),
                area_cv=float(area.std() / max(area.mean(), 1e-9)), median_area=float(med))


def save_overlays(frames, labels, out_dir: Path, tag: str, every: int = 4):
    from PIL import Image
    from skimage import exposure
    from skimage.segmentation import find_boundaries
    out_dir.mkdir(parents=True, exist_ok=True)
    for t in range(0, len(frames), every):
        img = exposure.rescale_intensity(frames[t].astype(float), out_range=(0, 255)).astype(np.uint8)
        rgb = np.stack([img] * 3, axis=-1)
        b = find_boundaries(labels[t], mode="outer"); rgb[b] = (255, 60, 60)
        Image.fromarray(rgb).save(out_dir / f"{tag}_frame{t:03d}.png")


def parse_frames(s: str) -> list[int]:
    a, b = s.split("-"); return list(range(int(a), int(b) + 1))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, required=True); ap.add_argument("--nd2", type=Path, action="append", required=True)
    ap.add_argument("--frames", action="append", required=True, help="a-b je --nd2")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--models", default="cpsam,cpsam_v2"); ap.add_argument("--flow", default="0.4,0.6,0.8")
    ap.add_argument("--cellprob", default="-1,0,1"); ap.add_argument("--niter", default="0,500"); ap.add_argument("--min-size", type=int, default=100)
    ap.add_argument("--gpu", action="store_true", default=True); ap.add_argument("--overlays-every", type=int, default=4)
    a = ap.parse_args()
    from cellpose import models
    a.out.mkdir(parents=True, exist_ok=True)
    grid = list(itertools.product(a.models.split(","), [float(x) for x in a.flow.split(",")],
                                  [float(x) for x in a.cellprob.split(",")], [int(x) for x in a.niter.split(",")]))
    print(f"{len(grid)} Einstellungen x {len(a.nd2)} Filme")
    movies = []
    for nd2_path, fr in zip(a.nd2, a.frames):
        frames, pipe, px = load_frames(nd2_path, a.config, parse_frames(fr))
        movies.append((nd2_path.stem, frames)); print(f"{nd2_path.name}: {len(frames)} Frames, {frames[0].shape}, {px} um/px")
    rows = []; loaded = {}
    for model_name, flow, cellprob, niter in grid:
        if model_name not in loaded:
            loaded[model_name] = models.CellposeModel(pretrained_model=model_name, gpu=a.gpu)
        model = loaded[model_name]; tag = f"{model_name}_f{flow:g}_p{cellprob:g}_n{niter}"
        for stem, frames in movies:
            t0 = time.perf_counter(); labels = [segment(model, f, flow, cellprob, niter, a.min_size) for f in frames]
            dt = (time.perf_counter() - t0) / len(frames)
            m = consistency_metrics(labels); m.update(model=model_name, flow=flow, cellprob=cellprob, niter=niter, movie=stem, s_per_frame=dt)
            rows.append(m); print(json.dumps(m))
            np.savez_compressed(a.out / f"labels_{stem}_{tag}.npz", labels=np.stack(labels))
            save_overlays(frames, labels, a.out / "overlays" / stem, tag, every=a.overlays_every)
        pd.DataFrame(rows).to_csv(a.out / "sweep_results.csv", index=False)
    df = pd.DataFrame(rows)
    summary = df.groupby(["model", "flow", "cellprob", "niter"])[["new_ids_per_object_frame", "gap_links_share", "merges_per_frame",
                                                                 "splits_per_frame", "large_splits_per_frame", "n_obj_mean", "n_obj_cv",
                                                                 "area_cv", "median_area", "s_per_frame"]].mean()
    summary.to_csv(a.out / "sweep_summary.csv"); print(summary.round(3).to_string())


if __name__ == "__main__":
    main()
