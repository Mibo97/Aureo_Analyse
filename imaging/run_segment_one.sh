#!/bin/bash
# Einen Film mit Pipeline v12 segmentieren + tracken (im slterm -p cuda Fenster, cellpose-Umgebung aktiv):
#   bash imaging/run_segment_one.sh /prj/.../Data/WT/pH/6/01_raw_data/<film>.nd2 /pfad/WT_PH_6config.yaml [--frames 0-12]
# Ausgabe: <exp>/02_processed_v12/ (labels_*.zarr, tracks_*.zarr, Ereignisse) und <exp>/03_results_v12/Single-Cell-Results_<stem>.csv
set -euo pipefail
FILE="${1:?nd2}"; CFG="${2:?Experiment-YAML}"; shift 2
mkdir -p logs
python imaging/cellpose_pipeline_v12.py --config "$CFG" --settings imaging/pipeline_template_v12.yaml --file "$FILE" "$@" \
    2>&1 | tee "logs/segment_$(basename "${FILE%.nd2}").log"
