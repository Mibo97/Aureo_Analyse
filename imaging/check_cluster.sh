#!/bin/bash
# Was bietet der Cluster? Einmal im slterm-Fenster ausfuehren und die Ausgabe schicken.
#   bash imaging/check_cluster.sh
echo "== host / user / dir";  hostname; whoami; pwd; date
echo "== slurm commands (leer = nicht vorhanden)"
for c in sbatch srun salloc squeue sinfo scancel slterm; do printf "  %-8s %s\n" "$c" "$(command -v $c 2>/dev/null || echo '-')"; done
echo "== what is slterm?"; type slterm 2>/dev/null | head -3; f=$(command -v slterm 2>/dev/null); [ -n "$f" ] && file "$f" 2>/dev/null && head -c 1500 "$f" 2>/dev/null | grep -v '^$' | head -25
echo "== current allocation (inside slterm these are set)"
for k in SLURM_JOB_ID SLURM_JOB_PARTITION SLURM_CPUS_ON_NODE SLURM_MEM_PER_NODE SLURM_JOB_NODELIST CUDA_VISIBLE_DEVICES TMPDIR; do echo "  $k=${!k}"; done
command -v squeue >/dev/null && squeue -u "$USER" 2>/dev/null | head -5
command -v sinfo  >/dev/null && sinfo -o "%P %a %l %D %G" 2>/dev/null | head -12
echo "== cpu / memory / gpu"; nproc; free -g | head -2; command -v nvidia-smi >/dev/null && nvidia-smi -L || echo "  nvidia-smi: -"
echo "== python env"; command -v python; python -c "import sys; print(sys.version.split()[0])"
python - <<'PY'
import importlib
for m in ("cellpose","torch","nd2","zarr","skimage","scipy","pandas","numpy","yaml"):
    try: print(f"  {m:<9} {getattr(importlib.import_module(m),'__version__','?')}")
    except Exception as e: print(f"  {m:<9} MISSING")
try:
    import torch; print("  cuda:", torch.cuda.is_available(), torch.cuda.get_device_name(0) if torch.cuda.is_available() else "")
except Exception as e: print("  torch:", e)
PY
echo "== project path"; ls -d /prj/microfluidic/ma_mimorde/Data 2>/dev/null || echo "  /prj/microfluidic/ma_mimorde/Data nicht sichtbar"
echo "== tools"; for c in git tmux screen nohup; do printf "  %-6s %s\n" "$c" "$(command -v $c 2>/dev/null || echo '-')"; done
