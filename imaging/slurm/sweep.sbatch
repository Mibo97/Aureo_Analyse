#!/bin/bash
# Segmentierungs-Sweep auf wenigen Frames (Phase A). Beispiel:
#   sbatch imaging/slurm/sweep.sbatch WT_GLC_0.75config.yaml sweep_out \
#       /prj/.../260805_Osc6_Rep1_ChamA4.nd2 40-52  /prj/.../260805_Osc6_PosCtrl_Rep1_ChamA1.nd2 100-112
#SBATCH -p cuda
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH -o logs/sweep_%j.out
set -euo pipefail
mkdir -p logs
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV:-cellpose}"
cd "${REPO_DIR:-$SLURM_SUBMIT_DIR}"
CONFIG="$1"; OUT="$2"; shift 2
ARGS=()
while [ $# -ge 2 ]; do ARGS+=(--nd2 "$1" --frames "$2"); shift 2; done
python imaging/sweep_segmentation.py --config "$CONFIG" --out "$OUT" "${ARGS[@]}" \
    --models "${MODELS:-cpsam,cpsam_v2}" --flow "${FLOW:-0.4,0.6,0.8}" --cellprob "${CELLPROB:--1,0,1}" --niter "${NITER:-0,500}"
