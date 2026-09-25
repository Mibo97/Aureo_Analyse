"""probe_env.py - Was kann der Rechenknoten? Einmal interaktiv oder per sbatch laufen lassen.

    python probe_env.py [--nd2 FILM.nd2] [--config CONFIG.yaml] [--frames 3]

Gibt aus: Host, GPU (nvidia-smi, torch), Paketversionen (cellpose, zarr, nd2, skimage), ob der
Projektpfad gemountet ist, und - mit --nd2 - die nd2-Metadaten (Achsen, Kanaele, Pixelgroesse in um)
sowie die Cellpose-Zeit je Frame mit der Config-Einstellung. Damit sind die offenen Cluster-Fragen aus
docs/tracking_plan.md beantwortet.
"""
from __future__ import annotations
import argparse, importlib, os, platform, shutil, subprocess, sys, time
from pathlib import Path


def ver(name):
    try:
        m = importlib.import_module(name); return getattr(m, "__version__", getattr(m, "version", "?"))
    except Exception as e:
        return f"MISSING ({type(e).__name__})"


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--nd2", type=Path); ap.add_argument("--config", type=Path)
    ap.add_argument("--frames", type=int, default=3); ap.add_argument("--project", default="/prj/microfluidic/ma_mimorde")
    a = ap.parse_args()
    print("host:", platform.node(), "| python:", sys.version.split()[0], "| cpus:", os.cpu_count())
    for k in ("SLURM_JOB_ID", "SLURM_JOB_PARTITION", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE", "CUDA_VISIBLE_DEVICES", "TMPDIR"):
        print(f"  {k} = {os.environ.get(k)}")
    if shutil.which("nvidia-smi"):
        print(subprocess.run(["nvidia-smi", "--query-gpu=name,memory.total,driver_version", "--format=csv"], capture_output=True, text=True).stdout)
    else:
        print("nvidia-smi: nicht gefunden")
    for m in ("torch", "cellpose", "zarr", "nd2", "skimage", "scipy", "numpy", "pandas", "tifffile", "yaml"):
        print(f"  {m:<10} {ver(m)}")
    try:
        import torch
        print("torch.cuda.is_available():", torch.cuda.is_available(),
              "| device:", torch.cuda.get_device_name(0) if torch.cuda.is_available() else "-",
              "| VRAM GB:", round(torch.cuda.get_device_properties(0).total_memory / 1e9, 1) if torch.cuda.is_available() else "-")
    except Exception as e:
        print("torch:", e)
    print("project path mounted:", Path(a.project).exists(), a.project)
    if a.nd2:
        import nd2, numpy as np
        with nd2.ND2File(str(a.nd2)) as f:
            print("nd2 sizes:", dict(f.sizes))
            try:
                v = f.voxel_size(); print(f"voxel size: x={v.x} y={v.y} z={v.z} um/px")
            except Exception as e:
                print("voxel_size:", e)
            try:
                chans = [c.channel.name for c in f.metadata.channels]; print("channels:", chans)
            except Exception as e:
                print("channels:", e)
            try:
                exp = f.experiment; print("experiment loops:", [(type(l).__name__, getattr(l, 'count', None)) for l in exp])
            except Exception as e:
                print("experiment:", e)
        if a.config:
            sys.path.insert(0, str(Path(__file__).resolve().parent.parent)); sys.path.insert(0, str(Path(__file__).resolve().parent))
            from sweep_segmentation import load_frames, segment
            from cellpose import models
            frames, pipe, px = load_frames(a.nd2, a.config, list(range(a.frames)))
            seg = pipe.config["segmentation"]
            model = models.CellposeModel(pretrained_model=seg["model_type"], gpu=seg.get("use_gpu", True))
            segment(model, frames[0], seg.get("flow_threshold", 0.4), seg["cellprob_threshold"], seg.get("niter"), seg["min_size_px"])  # warm-up
            t0 = time.perf_counter(); n = 0
            for fr in frames:
                lab = segment(model, fr, seg.get("flow_threshold", 0.4), seg["cellprob_threshold"], seg.get("niter"), seg["min_size_px"]); n += int(lab.max())
            dt = (time.perf_counter() - t0) / len(frames)
            print(f"cellpose {seg['model_type']}: {dt:.2f} s/frame on {frames[0].shape}, {n / len(frames):.0f} objects/frame")


if __name__ == "__main__":
    main()
