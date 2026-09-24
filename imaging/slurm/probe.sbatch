#!/bin/bash
# Einmal laufen lassen: beantwortet GPU-Typ, VRAM, Versionen, Mount, Pixelgroesse, Sekunden je Frame.
#   sbatch imaging/slurm/probe.sbatch /prj/.../01_raw_data/<film>.nd2 WT_GLC_0.75config.yaml
#SBATCH -p cuda
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH -o logs/probe_%j.out
set -euo pipefail
mkdir -p logs
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV:-cellpose}"     # Name der Umgebung mit cellpose 4.1.1, nd2, zarr
cd "${REPO_DIR:-$SLURM_SUBMIT_DIR}"
python imaging/probe_env.py --nd2 "${1:?nd2 file}" --config "${2:?config yaml}" --frames 3
