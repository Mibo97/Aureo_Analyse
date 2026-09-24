#!/bin/bash
# Re-Tracking EINES Experiments aus 02_processed/masks_*.zarr (kein Cellpose, keine GPU noetig):
#   bash imaging/run_track_one.sh /prj/microfluidic/ma_mimorde/Data/WT/pH/6
# Ergebnis: <exp>/03_results/Combined_Results_retracked.csv und <exp>/04_tracking/
set -euo pipefail
EXP="${1:?Experiment-Ordner, z.B. .../Data/WT/pH/6}"
mkdir -p logs
LOG="logs/track_$(basename "$(dirname "$(dirname "$EXP")")")_$(basename "$(dirname "$EXP")")_$(basename "$EXP").log"
python imaging/track_labels.py --batch "$EXP/02_processed" --combined "$EXP/03_results/Combined_Results.csv" \
    --out "$EXP/04_tracking" ${TRACK_ARGS:-} 2>&1 | tee "$LOG"
echo "-> Log: $LOG"
