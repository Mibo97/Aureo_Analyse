"""segment_all.py - alle Filme eines Datenbaums durch cellpose_pipeline_v12.py schicken, wieder aufrufbar.

    python imaging/segment_all.py /prj/.../Data --config-dir /pfad/zu/den/yamls --settings imaging/pipeline_template_v12.yaml \
        [--workers 3] [--gpus 0,1,2] [--config-pattern "{strain}_{osc_upper}_{period}config.yaml"] [--dry-run]

Je Film ein Unterprozess; fertige Filme (Single-Cell-Results_<stem>.csv in 03_results_v12) werden uebersprungen.
Die Experiment-YAML wird je Experiment ueber --config-pattern gesucht (Platzhalter strain, osc, osc_upper,
period), ersatzweise ueber ein Muster *<strain>*<period>*.yaml im Ordner; ohne Treffer wird das Experiment
mit Hinweis ausgelassen. Am Ende merge_results.py je Experiment.

Array-Jobs: --manifest manifests/movies.csv --index $SLURM_ARRAY_TASK_ID verarbeitet genau einen Film.
"""
from __future__ import annotations
import argparse, csv, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent


def find_config(config_dir: Path, pattern: str, strain: str, osc: str, period: str) -> Path | None:
    cand = config_dir / pattern.format(strain=strain, osc=osc, osc_upper=osc.upper(), osc_lower=osc.lower(), period=period)
    if cand.exists():
        return cand
    hits = sorted(p for p in config_dir.glob("*.y*ml") if strain.lower() in p.name.lower() and f"{period}" in p.name and osc.lower() in p.name.lower())
    return hits[0] if hits else None


def movies_of(data_root: Path):
    for raw in sorted(data_root.glob("*/*/*/01_raw_data")):
        exp = raw.parent; strain, osc, period = exp.relative_to(data_root).parts
        for nd2 in sorted(p for p in raw.glob("*.nd2") if "overview" not in p.name.lower()):
            yield dict(strain=strain, osc=osc, period=period, exp_dir=exp, nd2=nd2)


def run_one(m: dict, args) -> tuple[str, int]:
    cfg = find_config(args.config_dir, args.config_pattern, m["strain"], m["osc"], m["period"])
    if cfg is None:
        return (f"{m['nd2'].name}: keine Experiment-YAML fuer {m['strain']}/{m['osc']}/{m['period']} in {args.config_dir}", 2)
    out = m["exp_dir"] / args.results_subdir / f"Single-Cell-Results_{m['nd2'].stem}.csv"
    if out.exists() and not args.force:
        return (f"{m['nd2'].name}: schon fertig", 0)
    cmd = [sys.executable, str(HERE / "cellpose_pipeline_v12.py"), "--config", str(cfg), "--settings", str(args.settings),
           "--file", str(m["nd2"]), "--exp-dir", str(m["exp_dir"]), "--processed-subdir", args.processed_subdir,
           "--results-subdir", args.results_subdir] + (["--no-overlay"] if args.no_overlay else []) + (["--skeleton"] if args.skeleton else [])
    if args.dry_run:
        return (" ".join(cmd), 0)
    env = dict(os.environ)
    if args.gpus:
        gpus = args.gpus.split(","); env["CUDA_VISIBLE_DEVICES"] = gpus[m["slot"] % len(gpus)]
    Path("logs").mkdir(exist_ok=True)
    with open(Path("logs") / f"segment_{m['nd2'].stem}.log", "w") as f:
        rc = subprocess.call(cmd, stdout=f, stderr=subprocess.STDOUT, env=env)
    return (f"{m['nd2'].name}: {'ok' if rc == 0 else 'FEHLER (rc ' + str(rc) + ')'}  Log: logs/segment_{m['nd2'].stem}.log", rc)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data_root", type=Path); ap.add_argument("--config-dir", type=Path, required=True)
    ap.add_argument("--settings", type=Path, default=HERE / "pipeline_template_v12.yaml")
    ap.add_argument("--config-pattern", default="{strain}_{osc_upper}_{period}config.yaml")
    ap.add_argument("--processed-subdir", default="02_processed_v12"); ap.add_argument("--results-subdir", default="03_results_v12")
    ap.add_argument("--workers", type=int, default=1); ap.add_argument("--gpus", default="", help="z.B. 0,1,2 - Worker werden reihum auf die GPUs verteilt")
    ap.add_argument("--no-overlay", action="store_true"); ap.add_argument("--skeleton", action="store_true")
    ap.add_argument("--force", action="store_true"); ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--manifest", type=Path, default=None); ap.add_argument("--index", type=int, default=None)
    ap.add_argument("--no-merge", action="store_true")
    args = ap.parse_args()
    movies = list(movies_of(args.data_root))
    if args.manifest is not None and args.index is not None:   # Array-Job: genau ein Film
        rows = list(csv.DictReader(open(args.manifest))); r = rows[args.index]
        movies = [m for m in movies if str(m["nd2"]) == r["nd2"]]
        if not movies:
            raise SystemExit(f"Film aus Manifest nicht gefunden: {r['nd2']}")
    for i, m in enumerate(movies):
        m["slot"] = i
    print(f"{len(movies)} Filme, {args.workers} Worker, GPUs {args.gpus or 'wie Umgebung'}", flush=True)
    n_err = 0
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        for msg, rc in ex.map(lambda m: run_one(m, args), movies):
            n_err += int(rc != 0); print(msg, flush=True)
    if not args.no_merge and not args.dry_run and args.index is None:
        subprocess.call([sys.executable, str(HERE / "merge_results.py"), "--all", str(args.data_root), "--results-subdir", args.results_subdir])
    print(f"fertig, {n_err} Fehler/uebersprungene Experimente", flush=True)


if __name__ == "__main__":
    main()
