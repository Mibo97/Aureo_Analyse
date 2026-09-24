#!/bin/bash
# Re-Tracking ALLER Experimente, nacheinander oder parallel, ohne sbatch. Wieder aufrufbar: fertige
# Experimente (Combined_Results_retracked.csv vorhanden) werden uebersprungen.
#   bash imaging/run_track_all.sh /prj/microfluidic/ma_mimorde/Data 4      # 4 Experimente gleichzeitig
# Laeuft lange (Stunden). Im slterm-Fenster mit nohup starten, damit ein geschlossenes Fenster nichts abbricht:
#   nohup bash imaging/run_track_all.sh /prj/.../Data 4 > logs/track_all.log 2>&1 &
#   tail -f logs/track_all.log
set -euo pipefail
ROOT="${1:?Datenwurzel, z.B. /prj/microfluidic/ma_mimorde/Data}"; PAR="${2:-2}"
mkdir -p logs manifests
python imaging/slurm/make_manifest.py "$ROOT" --out manifests
python - "$PAR" <<'PY'
import csv, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
par = int(sys.argv[1])
rows = list(csv.DictReader(open("manifests/experiments.csv")))
todo = [r for r in rows if int(r["n_zarr"]) > 0 and r["combined_exists"] == "True"
        and not os.path.exists(os.path.join(r["results_dir"], "Combined_Results_retracked.csv"))]
print(f"{len(rows)} Experimente, {len(todo)} noch zu tracken, {par} parallel", flush=True)
def run(r):
    exp = r["exp_dir"]; log = f"logs/track_{r['strain']}_{r['osc_type']}_{r['period']}.log"
    with open(log, "w") as f:
        rc = subprocess.call(["python", "imaging/track_labels.py", "--batch", os.path.join(exp, "02_processed"),
                              "--combined", os.path.join(exp, "03_results", "Combined_Results.csv"),
                              "--out", os.path.join(exp, "04_tracking")] + os.environ.get("TRACK_ARGS", "").split(), stdout=f, stderr=subprocess.STDOUT)
    print(("ok   " if rc == 0 else "FEHLER ") + exp + f"  (Log: {log})", flush=True)
with ThreadPoolExecutor(max_workers=par) as ex:
    list(ex.map(run, todo))
print("fertig", flush=True)
PY
