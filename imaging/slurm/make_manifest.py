"""make_manifest.py - Liste aller Experimente und Filme unter dem Datenbaum.

    python imaging/slurm/make_manifest.py /prj/microfluidic/ma_mimorde/Data --out manifests

Schreibt manifests/experiments.csv (eine Zeile je <strain>/<osc_type>/<period>: Pfade zu 01_raw_data,
02_processed, 03_results, Zahl der nd2-Filme und der masks_*.zarr) und manifests/movies.csv (eine Zeile je
nd2). Die Zeilennummer (0-basiert, ohne Kopf) ist der Index fuer die Array-Jobs.
"""
import argparse, csv
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("data_root", type=Path); ap.add_argument("--out", type=Path, default=Path("manifests"))
    a = ap.parse_args(); a.out.mkdir(parents=True, exist_ok=True)
    exps, movies = [], []
    for raw in sorted(a.data_root.glob("*/*/*/01_raw_data")):
        exp = raw.parent; rel = exp.relative_to(a.data_root).parts
        nd2s = sorted(p for p in raw.glob("*.nd2") if "overview" not in p.name.lower())
        zarrs = sorted((exp / "02_processed").glob("masks_*.zarr"))
        combined = exp / "03_results" / "Combined_Results.csv"
        exps.append(dict(strain=rel[0], osc_type=rel[1], period=rel[2], exp_dir=str(exp), raw_dir=str(raw),
                         processed_dir=str(exp / "02_processed"), results_dir=str(exp / "03_results"),
                         n_nd2=len(nd2s), n_zarr=len(zarrs), combined_exists=combined.exists()))
        for p in nd2s:
            movies.append(dict(strain=rel[0], osc_type=rel[1], period=rel[2], nd2=str(p), stem=p.stem, exp_dir=str(exp),
                               zarr_exists=(exp / "02_processed" / f"masks_{p.stem}.zarr").exists()))
    for name, rows in (("experiments.csv", exps), ("movies.csv", movies)):
        with open(a.out / name, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["empty"]); w.writeheader(); w.writerows(rows)
    print(f"{len(exps)} Experimente, {len(movies)} Filme -> {a.out}/experiments.csv, {a.out}/movies.csv")
    print(f"  mit masks_*.zarr: {sum(m['zarr_exists'] for m in movies)} Filme; Combined_Results vorhanden: {sum(e['combined_exists'] for e in exps)} Experimente")


if __name__ == "__main__":
    main()
