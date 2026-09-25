#!/bin/bash
# Interaktiv (im slterm -p cuda Fenster, Umgebung aktiviert, im Repo-Ordner):
#   bash imaging/run_probe.sh /prj/.../01_raw_data/<film>.nd2 WT_GLC_0.75config.yaml
set -euo pipefail
mkdir -p logs
LOG="logs/probe_$(date +%Y%m%d_%H%M).log"
python imaging/probe_env.py --nd2 "${1:?nd2 file}" --config "${2:?config yaml}" --frames 3 2>&1 | tee "$LOG"
echo "-> gespeichert in $LOG"
