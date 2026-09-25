#!/bin/bash
# Segmentierungs-Sweep interaktiv (im slterm -p cuda Fenster):
#   bash imaging/run_sweep.sh WT_GLC_0.75config.yaml sweep_out \
#       /prj/.../260805_Osc6_Rep1_ChamA4.nd2 40-52 /prj/.../260805_Osc6_PosCtrl_Rep1_ChamA1.nd2 100-112
# Kleiner anfangen (schneller Test, 2 Einstellungen):  MODELS=cpsam_v2 FLOW=0.8 CELLPROB=0 NITER=0,500 bash imaging/run_sweep.sh ...
set -euo pipefail
CONFIG="${1:?config yaml}"; OUT="${2:?output dir}"; shift 2
ARGS=(); while [ $# -ge 2 ]; do ARGS+=(--nd2 "$1" --frames "$2"); shift 2; done
mkdir -p logs
python imaging/sweep_segmentation.py --config "$CONFIG" --out "$OUT" "${ARGS[@]}" \
    --models="${MODELS:-cpsam,cpsam_v2}" --flow="${FLOW:-0.4,0.6,0.8}" --cellprob="${CELLPROB:--1,0,1}" --niter="${NITER:-0,500}" \
    2>&1 | tee "logs/sweep_$(date +%Y%m%d_%H%M).log"
