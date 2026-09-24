# imaging/ — Werkzeuge um die Cellpose-Pipeline herum

Die Bild-Pipeline selbst ist `../cellpose_pipeline_v11.py` (Segmentierung + Tracking im selben Loop,
Konfiguration z.B. `../WT_GLC_0.75config.yaml`). Hier liegen die Werkzeuge, die auf ihren Ausgaben
arbeiten oder ihre Einstellungen pruefen. Befund und Plan: `../docs/tracking_diagnosis.md`,
`../docs/tracking_plan.md`.

| Datei | Zweck | Laeuft wo |
| --- | --- | --- |
| `track_labels.py` | Neues Tracking auf gespeicherten Label-Stacks (`masks_<stem>.zarr`), ohne Cellpose: Gedaechtnis, Kosten statt IoU-Tor, Merge/Split aus der Ueberlappung, `parent_track_id` fuer beruehrende neue Objekte. `--batch` verarbeitet ein Experiment und schreibt `Combined_Results_retracked.csv` (alte Spalten + neue IDs). | CPU, Minuten je Film |
| `synthetic_tracking_test.py` | Ground-Truth-Test des Trackers gegen einen v11-artigen Tracker auf einem synthetischen Film | ueberall |
| `sweep_segmentation.py` | Gitter ueber Cellpose-Einstellungen auf wenigen aufeinanderfolgenden Frames; Bewertung ohne Hand-Labels (Konsistenz, Merge/Split, Aufspaltung grosser Objekte) plus Overlays | GPU |
| `probe_env.py` | Was kann der Knoten: GPU, Versionen, Mount, nd2-Metadaten (Pixelgroesse), Sekunden je Frame | GPU |
| `slurm/make_manifest.py` | Liste aller Experimente/Filme fuer Array-Jobs | ueberall |
| `slurm/probe.sbatch`, `slurm/sweep.sbatch`, `slurm/track.sbatch` | Jobskripte (Partition `cuda` bzw. `cpu`; Umgebungsname ueber `CONDA_ENV`, Standard `cellpose`) | Cluster |

## Reihenfolge auf dem Cluster

```bash
# 0. einmalig: Repo auf den Cluster, Umgebung (conda) mit cellpose 4.1.1, nd2, zarr, scikit-image
cd /pfad/zum/Aureo_Analyse && mkdir -p logs

# 1. Knoten und Daten pruefen (GPU, Versionen, Pixelgroesse, s/Frame)
sbatch imaging/slurm/probe.sbatch /prj/.../01_raw_data/<film>.nd2 WT_GLC_0.75config.yaml
cat logs/probe_*.out

# 2. Re-Tracking aller Experimente aus den vorhandenen Masken (kein Cellpose)
python imaging/slurm/make_manifest.py /prj/microfluidic/ma_mimorde/Data --out manifests
N=$(($(wc -l < manifests/experiments.csv) - 1))
sbatch --array=0-$((N-1))%16 imaging/slurm/track.sbatch manifests/experiments.csv
#    -> je Experiment 03_results/Combined_Results_retracked.csv und 04_tracking/

# 3. Analyse auf den neuen Tabellen
export AUREO_RESULTS_PATTERN="Combined_Results_retracked.*"
python analyse_pipeline/run_analysis.py

# 4. Segmentierungs-Sweep (Phase A), vier Filme, je 12 aufeinanderfolgende Frames
sbatch imaging/slurm/sweep.sbatch WT_GLC_0.75config.yaml sweep_out \
    /prj/.../260805_Osc6_Rep1_ChamA4.nd2 40-52 /prj/.../260805_Osc6_PosCtrl_Rep1_ChamA1.nd2 100-112 \
    /prj/.../260805_Osc6_NegCtrl_Rep2_ChamA14.nd2 20-32 /prj/.../WT/Glc/6/01_raw_data/<film>.nd2 60-72
#    -> sweep_out/sweep_summary.csv und sweep_out/overlays/
```

Ein einzelner Film zum Ausprobieren, ohne Slurm:

```bash
python imaging/track_labels.py /prj/.../02_processed/masks_260805_Osc6_Rep1_ChamA4.zarr --out /tmp/tracktest
```

Partition `cpu` in `track.sbatch` ggf. an den Cluster anpassen (`sinfo` zeigt die Namen). Fehlgeschlagene
Array-Indizes werden mit `--array=3,17` einzeln wiederholt.

## Was `track_labels.py` in die Tabelle schreibt

`track_id` (neu), `track_id_v11` (alt = Label im Stack), `parent_track_id` (−1 = keins),
`parent_area_ratio`, `link_type` (`continued`, `gap`, `long_range`, `split`, `new_touching`, `new`),
`gap_frames`, `overlap_prev_frac`, `n_candidates`, `axis_major`, `axis_minor`, `perimeter`, `circularity`,
optional `skeleton_px` (`--skeleton`). Alle alten Spalten (Fluoreszenz, Metadaten) bleiben erhalten.
Ereignisse (merge, split, gap, long_range) stehen je Film in `04_tracking/<stem>_events.csv`.

## Synthetischer Test (`synthetic_tracking_test.py`, zwei Seeds, 80 Frames)

| | v11-artig (IoU-Tor, kein Gedaechtnis) | `track_labels.py` |
| --- | --- | --- |
| Verknuepfungen ueber 1 Frame wiedergefunden | 0.97 | 0.97 |
| Verknuepfungen ueber 2–3 Frames (Objekt zwischendurch nicht segmentiert) | 0.06–0.08 | 0.90–0.92 |
| Spuren je Objekt (Objekte mit >= 5 Auftritten) | 6.4–6.7 | 2.7–2.9 |
| Identitaetswechsel (Spur auf zwei Objekten) | 0–1 | 1–2 |
| Knospen mit richtigem Elternteil | 0 | 0.95–1.00 |
