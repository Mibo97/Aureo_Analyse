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

## Reihenfolge auf dem Cluster, interaktiv (ohne sbatch)

Alles hier laeuft in dem Fenster, das `slterm -p cuda` oeffnet (fuer das Re-Tracking reicht ein
Fenster ohne GPU, falls es eine CPU-Partition gibt). Vorher wie gewohnt die Umgebung aktivieren
(miniforge, `conda activate <env>`), dann in den Repo-Ordner wechseln.

```bash
# 0. einmalig: Repo auf den Cluster holen (Branch mit den Werkzeugen) und Pakete pruefen
git clone -b claude/adoring-allen-70cicr https://github.com/Mibo97/Aureo_Analyse.git
cd Aureo_Analyse
pip install -r imaging/requirements.txt          # in der aktivierten Umgebung; cellpose ist schon da
bash imaging/check_cluster.sh                    # Ausgabe schicken: zeigt, was der Cluster kann

# 1. Knoten und Daten pruefen (GPU, Versionen, Pixelgroesse, Sekunden je Frame) - ca. 2 Minuten
bash imaging/run_probe.sh /prj/microfluidic/ma_mimorde/Data/WT/pH/6/01_raw_data/<film>.nd2 WT_GLC_0.75config.yaml

# 2. Re-Tracking des QC-Batches (WT/pH/6) aus den vorhandenen Masken - Minuten, keine GPU
bash imaging/run_track_one.sh /prj/microfluidic/ma_mimorde/Data/WT/pH/6
#    -> .../WT/pH/6/03_results/Combined_Results_retracked.csv und .../WT/pH/6/04_tracking/
python analyse_pipeline/diagnose_tracking.py /prj/.../WT/pH/6/03_results/Combined_Results_retracked.csv --out /prj/.../WT/pH/6/04_tracking
#    schicken: 04_tracking/tracking_summary.csv, tracking_new_id_categories.csv, tracking_relink_prototype.csv

# 3. Re-Tracking aller Experimente (Stunden; wieder aufrufbar, fertige werden uebersprungen)
nohup bash imaging/run_track_all.sh /prj/microfluidic/ma_mimorde/Data 4 > logs/track_all.log 2>&1 &
tail -f logs/track_all.log                        # Strg+C beendet nur die Anzeige, nicht den Lauf

# 4. Analyse auf den neuen Tabellen
export AUREO_RESULTS_PATTERN="Combined_Results_retracked.*"
python analyse_pipeline/run_analysis.py

# 5. Segmentierungs-Sweep (Phase A, GPU), erst klein testen, dann das ganze Gitter
MODELS=cpsam_v2 FLOW=0.8 CELLPROB=0 NITER=0,500 bash imaging/run_sweep.sh WT_GLC_0.75config.yaml sweep_test \
    /prj/.../WT/pH/6/01_raw_data/260805_Osc6_Rep1_ChamA4.nd2 40-52
bash imaging/run_sweep.sh WT_GLC_0.75config.yaml sweep_out \
    /prj/.../260805_Osc6_Rep1_ChamA4.nd2 40-52 /prj/.../260805_Osc6_PosCtrl_Rep1_ChamA1.nd2 100-112 \
    /prj/.../260805_Osc6_NegCtrl_Rep2_ChamA14.nd2 20-32 /prj/.../WT/Glc/6/01_raw_data/<film>.nd2 60-72
#    -> sweep_out/sweep_summary.csv und sweep_out/overlays/
```

Wichtig bei interaktiven Sitzungen: die Sitzung hat ein Zeitlimit (steht in `squeue -u $USER`, falls
vorhanden). Schritt 3 deshalb mit `nohup` starten und bei Abbruch einfach erneut aufrufen. Wenn `git` auf
dem Cluster fehlt: das Repo als ZIP von GitHub laden und den Ordner `imaging/` sowie
`analyse_pipeline/diagnose_tracking.py` hochladen.

## Mit sbatch (optional, wenn `sbatch` vorhanden ist)

`imaging/slurm/*.sbatch` sind dieselben Schritte als Batch-Jobs (Probe, Sweep, Re-Tracking als
Array-Job ueber `manifests/experiments.csv`). Auf diesem Cluster (Slurm 19, Partitionen `cebitec01`,
`kurs`, `cuda`, `interactive`, alle ohne Zeitlimit) gibt es keine GRES-Definition fuer GPUs: GPU-Jobs
laufen einfach mit `-p cuda`, `--gres=gpu:1` wuerde abgelehnt. CPU-Jobs gehen an `cebitec01`.
Umgebungsname ueber `CONDA_ENV` (Standard `cellpose`), conda-Installation ueber `CONDA_BASE`.
Abschicken vom Login-Knoten aus (nicht aus dem slterm-Fenster): `sbatch imaging/slurm/probe.sbatch ...`;
Stand mit `squeue -u $USER`, Ausgabe in `logs/`.

## Was `track_labels.py` in die Tabelle schreibt

`track_id` (neu), `track_id_v11` (alt = Label im Stack), `parent_track_id` (−1 = keins),
`parent_area_ratio`, `link_type` (`continued`, `gap`, `long_range`, `unmerge`, `split`, `new_touching`, `new`),
`gap_frames`, `overlap_prev_frac`, `n_candidates`, `axis_major`, `axis_minor`, `perimeter`, `circularity`,
optional `skeleton_px` (`--skeleton`). Alle alten Spalten (Fluoreszenz, Metadaten) bleiben erhalten.
Ereignisse (merge, split, gap, long_range) stehen je Film in `04_tracking/<stem>_events.csv`.

## Synthetischer Test (`synthetic_tracking_test.py`, zwei Seeds, 80 Frames)

| | v11-artig (IoU-Tor, kein Gedaechtnis) | `track_labels.py` |
| --- | --- | --- |
| Verknuepfungen ueber 1 Frame wiedergefunden | 0.97 | 0.97 |
| Verknuepfungen ueber 2–3 Frames (Objekt zwischendurch nicht segmentiert) | 0.06–0.08 | 0.95–0.96 |
| Spuren je Objekt (Objekte mit >= 5 Auftritten) | 6.4–6.7 | 2.3 |
| Identitaetswechsel (Spur auf zwei Objekten) | 0–1 | 1–3 |
| Knospen mit richtigem Elternteil | 0 | 0.95–1.00 |
| dito mit 12 % verschmolzenen Masken je Objekt-Frame (`--p-merge 0.12`): Spuren je Objekt | 8.2–8.6 | 2.6–2.7 |
