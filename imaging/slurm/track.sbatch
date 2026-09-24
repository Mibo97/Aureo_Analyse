#!/bin/bash
# Re-Tracking aller Experimente aus den vorhandenen masks_*.zarr (CPU, kein Cellpose).
#   python imaging/slurm/make_manifest.py /prj/.../Data --out manifests
#   N=$(($(wc -l < manifests/experiments.csv) - 1)); sbatch --array=0-$((N-1))%16 imaging/slurm/track.sbatch manifests/experiments.csv
# Schreibt je Experiment 03_results/Combined_Results_retracked.csv und 04_tracking/ (Tabellen, Ereignisse, tracks_*.zarr).
#SBATCH -p cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH -o logs/track_%A_%a.out
set -euo pipefail
mkdir -p logs
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV:-cellpose}"
cd "${REPO_DIR:-$SLURM_SUBMIT_DIR}"
MANIFEST="${1:?manifests/experiments.csv}"
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 2))p" "$MANIFEST")
EXP_DIR=$(echo "$LINE" | python -c "import csv,sys; r=next(csv.reader(sys.stdin)); print(r[3])")
echo "task $SLURM_ARRAY_TASK_ID: $EXP_DIR"
python imaging/track_labels.py --batch "$EXP_DIR/02_processed" --combined "$EXP_DIR/03_results/Combined_Results.csv" \
    --out "$EXP_DIR/04_tracking" ${TRACK_ARGS:-}
