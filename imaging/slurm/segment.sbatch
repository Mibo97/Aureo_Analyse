#!/bin/bash
# Pipeline v12 als Array-Job, ein Film je Task (vom Login-Knoten abschicken):
#   python imaging/slurm/make_manifest.py /prj/.../Data --out manifests
#   N=$(($(wc -l < manifests/movies.csv) - 1)); CONFIG_DIR=/pfad/zu/den/yamls sbatch --array=0-$((N-1))%3 imaging/slurm/segment.sbatch
#   danach: python imaging/merge_results.py --all /prj/.../Data
#SBATCH -p cuda
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH -o logs/segment_%A_%a.out
set -euo pipefail
mkdir -p logs
# conda: erst "conda" auf dem PATH, sonst bekannte Installationen (CONDA_BASE setzt den Pfad fest)
if command -v conda >/dev/null 2>&1; then eval "$(conda shell.bash hook)"; else
  for d in "${CONDA_BASE:-}" /prj/microfluidic/ma_mimorde/Miniconda3/Miniconda3 "$HOME/miniforge3" "$HOME/miniconda3"; do
    [ -n "$d" ] && [ -f "$d/etc/profile.d/conda.sh" ] && source "$d/etc/profile.d/conda.sh" && break; done; fi
conda activate "${CONDA_ENV:-cellpose}"
cd "${REPO_DIR:-$SLURM_SUBMIT_DIR}"
python imaging/segment_all.py "${DATA_ROOT:-/prj/microfluidic/ma_mimorde/Data}" --config-dir "${CONFIG_DIR:?CONFIG_DIR setzen}" \
    --manifest manifests/movies.csv --index "$SLURM_ARRAY_TASK_ID" --workers 1
