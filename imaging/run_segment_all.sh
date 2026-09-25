#!/bin/bash
# Alle Filme mit Pipeline v12 (wieder aufrufbar, fertige werden uebersprungen). Drei GPUs auf dem Knoten:
#   nohup bash imaging/run_segment_all.sh /prj/microfluidic/ma_mimorde/Data /pfad/zu/den/yamls 3 0,1,2 > logs/segment_all.log 2>&1 &
#   tail -f logs/segment_all.log
# Danach liegen je Experiment 03_results_v12/Combined_Results.csv; Analyse mit AUREO_RESULTS_SUBDIR=03_results_v12.
set -euo pipefail
ROOT="${1:?Datenwurzel}"; CFGDIR="${2:?Ordner mit den Experiment-YAMLs}"; WORKERS="${3:-1}"; GPUS="${4:-}"
mkdir -p logs
python imaging/segment_all.py "$ROOT" --config-dir "$CFGDIR" --settings imaging/pipeline_template_v12.yaml \
    --workers "$WORKERS" ${GPUS:+--gpus "$GPUS"} ${SEGMENT_ARGS:-}
