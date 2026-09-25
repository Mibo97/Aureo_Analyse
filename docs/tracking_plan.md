# Plan: segmentation tuning, tracking before filtering, cell types, batch jobs

Status: proposal for discussion. Nothing in the pipeline is changed by this document. Facts come from
`WT_GLC_0.75config.yaml`, `cellpose_pipeline_v11.py`, `docs/tracking_diagnosis.md` and Rensink et al. 2026
(Fungal Biology 130:101754).

## 0. What the config and the QC batch tell us

- `rotation_angle_deg: null`, so the rotation/crop bug is dormant. It is still worth fixing, it is three lines.
- The ROI filters are loose (`min_solidity 0.40`, `max_eccentricity 0.99`). Only the border band (10 px) and
  `min_area_px 200` actually drop objects. In the QC batch, new IDs of the kinds "isolated" and "missing for
  one frame" sit near the crop edge 3× more often than objects in general (22 % and 23 % vs 8 %), and 7 % of
  new IDs are just above 200 px versus 1 % of all objects. Tracking before filtering removes exactly these.
- The tracking values differ from the code defaults (`small_object_radius_px 35`, `max_area_growth 7`). Under
  these values 66 % of the same-object losses would still fail the IoU gate and 28 % the 50 px distance gate,
  so the tracker rewrite stands.
- The largest group of new IDs, an object appearing next to a neighbour whose mask did not change (51 %), is
  neither border nor mask splitting. Cellpose does not find the same objects in every frame, or cells arrive.
  A segmentation setting is therefore judged by its frame-to-frame consistency, not only by single frames.
- `flow_threshold 0.8` and `niter 500` were chosen so that constricted, figure-eight objects are not split.
  If those objects are mother-and-bud pairs, this setting delays bud detection, which is what the lineage
  needs to see early. If they are septated swollen cells, the setting is right. This has to be decided by eye.
- The pixel size is not in the table. It is needed for areas in µm² and for cell-type size classes; nd2 files
  carry it (`ND2File.voxel_size()`), the pipeline can write it into every row.

## 1. Phase A: segmentation sweep on a small validation set (GPU, hours)

Validation set: three movies of the QC batch (one oscillation chamber, one PosCtrl at high density, one
NegCtrl at low density) plus one Glc movie, 12 consecutive frames each from a dense and a sparse phase.

Grid: `model_type` {cpsam, cpsam_v2}, `flow_threshold` {0.4, 0.6, 0.8}, `cellprob_threshold` {−1, 0, +1},
`niter` {null, 500}; `min_size_px` 100 for all (objects below 200 px are only flagged later).

Metrics without hand labels, on the consecutive frames:
1. consistency: new-ID categories per object-frame from `diagnose_tracking.py` (fewer "next to an unchanged
   neighbour" and "missing for one frame" is better);
2. merge/split rate from mask overlap between consecutive frames (one mask at t over two at t−1 and vice versa);
3. object count and total mask area per frame (stability).

Then, for the three best settings, overlays of five frames each for a by-eye count of merged, split, missed and
false objects, about 30 minutes of annotation. If no setting is consistent enough, the fallback is fine-tuning
the Cellpose model on 20–40 corrected frames, which is a separate step.

Checkpoint 1: table of the sweep, overlays, one chosen setting.

## 2. Phase B: pipeline changes

1. `segment_frame` keeps every Cellpose object at or above `min_size_px` and adds flags instead of dropping:
   `at_border`, `below_min_area`, `above_max_area`, `low_solidity`, `high_eccentricity`. A config switch
   `roi_filter.mode: flag | drop` keeps the old behaviour available.
2. Features per object: area, centroid, bbox, perimeter, circularity, major and minor axis, eccentricity,
   solidity, mean and std of the phase-contrast signal inside the mask (contrast), optionally skeleton length.
3. Raw label masks are saved as `labels_<stem>.zarr` (Cellpose labels), track masks as `tracks_<stem>.zarr`.
4. Tracking becomes a separate pass `track_labels.py` on CPU that reads the label stack and the feature table:
   Hungarian assignment per frame with a memory of three frames; cost from centroid distance in radii, area
   ratio and mask overlap; gates on distance (1.5 r·√gap + 10 px) and area ratio (2.5), no IoU gate; a second
   pass in sparse frames that links leftover pairs over a large radius only when they are mutually unique;
   merge and split detected by overlap with two predecessors or two successors; a new object whose mask touches
   an existing track gets `parent_track_id`; an empty frame does not reset. Output columns: `track_id`,
   `parent_track_id`, `link_type` (continued, gap, long_range, new), `gap_frames`.
5. Rotation applied to the full frame before chamber detection and to every channel; `um_per_px` from the nd2
   metadata in every row; Cellpose version and config in the log.
6. Command line: `--file` for one movie, `--config-template` with io paths derived from the file's location,
   `--merge` to build Combined_Results per experiment. Needed for array jobs.

Checkpoint 2, in two steps that need no full GPU run:
- the new tracker on the existing `02_processed/masks_*.zarr` stacks of the QC batch (unique values per frame
  are objects), scored as in the diagnosis: fragmentation, manual merges reproduced, lineage events in the
  sparse window;
- the QC batch re-segmented with the chosen setting, same scoring.
Only then the full data set.

## 3. Phase C: analysis side

- `data_loading` accepts the new columns; which flags exclude an object is set in `config.py`, not in the
  imaging pipeline.
- `lineage.py` uses `parent_track_id` when present and falls back to the heuristic otherwise; a bud counts only
  if it persists two frames.
- `relink.py` keeps the prototype linker for tables without the new columns.
- Validation as before: `validate_lineage.py` against the manual QC, `lv_03`, window events.

Checkpoint 3: QC batch through the whole chain, old versus new.

## 4. Cell types after Rensink et al. 2026 (Table 2)

What the paper does: objects from an oCelloScope (brightfield, low magnification, static wells, PDB pH 5,
26 °C) are typed by four rules: area (≤ 275 vs > 275 px), circularity (> 0.85 vs ≤ 0.85), thinned length
(≤ 30 vs > 30 px, hyphae) and contrast (> 0.4, melanised chlamydospores). Blastoconidia become swollen cells
after about 6 h, swollen cells become hyphae after 10–13 h, and after 10 h the new blastoconidia are produced
by swollen cells (3.8–5.8 each) and hyphae (1–4.9 each), no longer by blastoconidia. Li et al. 2009, cited
there: blastoconidia turn into swollen cells at pH 4.5 but stay blastoconidia at pH 6.

The thresholds are in their pixels and their optics and do not transfer. The features do. Plan:

- Phase B already adds circularity, axes, contrast and optionally skeleton length. Calibrate the four rules on
  our data: area bimodality for blastoconidia versus swollen cells in µm², circularity and skeleton length for
  hyphae, dark contrast for chlamydospores. Validate on about 100 hand-labelled objects (confusion matrix).
- For tracking: gates per type. Hyphae grow along their axis, so allow a larger area ratio and a displacement
  along the major axis; swollen cells hardly move; blastoconidia are the ones that get flushed, so the
  long-range rule applies to them only. A bud is a blastoconidium by definition, so `parent_track_id` is only
  assigned to small round objects touching a swollen cell, hypha or blastoconidium.
- For interpretation: the composition per frame (fractions of blastoconidia, swollen cells, hyphae over 22 h)
  is a readout that needs no tracks at all, and it turns the eccentricity endpoint into a hyphal fraction. The
  budding rate can be split by mother type, which is how Rensink report division. For the pH series the
  Li 2009 result predicts a composition shift between pH 4.5 and 6; whether our feast and famine pH values
  bracket that is a question below.
- Caveat: their timings are for static PDB at pH 5; in a chamber with flow and a different medium they are
  orientation, not expectation.

## 5. Batch jobs on the Slurm cluster (partition cuda)

Files to add under `slurm/`:

- `make_manifest.py`: walks the data tree and writes `manifest.csv`, one nd2 per line with strain, osc_type,
  period, replicate, chamber and the derived output directories.
- `segment.sbatch`: array job, one task per movie.

```
#!/bin/bash
#SBATCH -p cuda
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --array=0-564%8
#SBATCH -o logs/segment_%A_%a.out
module load cuda            # or: source activate cellpose-env
FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 2))p" manifest.csv | cut -d, -f1)
python cellpose_pipeline.py --config-template template.yaml --file "$FILE" --scratch "$TMPDIR"
```

- `track.sbatch`: CPU array job with `--dependency=afterok:<segment job id>`, runs `track_labels.py` per movie.
- `merge.sbatch`: one task per experiment, builds Combined_Results.

Notes: `%8` limits the number of GPUs used at once; failed tasks are re-run with `--array=12,57`; the zarr
stacks are written to the node's `$TMPDIR` and copied to `02_processed` at the end; a first test is done
interactively with `srun -p cuda --gres=gpu:1 --pty bash` on one movie to measure time per frame, which sets
`--time`. At roughly two seconds per frame a 132-frame movie takes about ten minutes including I/O, so
565 movies on eight GPUs take about twelve hours.

## 6. Answers so far and what they change

1. **Pixel size:** probably in the nd2 metadata. `imaging/probe_env.py` prints `voxel_size()`; the
   sweep script prints it per movie. Until confirmed, all areas stay in pixels.
2. **Figure-eight objects are septated swollen cells.** The segmentation setting (`flow_threshold 0.8`,
   `niter 500`) is right and stays the base of the sweep; the sweep scores splits of large objects as
   errors (`large_splits_per_frame`).
3. **`masks_*.zarr` exist for all experiments.** The tracker runs on them on CPU for every experiment
   before any GPU time; only border-band objects and objects below 200 px are missing from those stacks.
4. **Cluster:** miniconda3, Cellpose 4.1.1. GPU type, memory, wall time, mounts are answered by
   `sbatch imaging/slurm/probe.sbatch <film.nd2> <config.yaml>`.
5. **pH series: feast pH 3, famine pH 8.** Li et al. 2009 (cited by Rensink): chlamydospores from
   swollen cells below pH 3, blastoconidia to swollen cells at pH 4.5, stable at pH 6. Feast sits at the
   chlamydospore edge and famine outside the reported range, so a composition readout per frame is worth
   having for this series in particular.

## 6a. Built at this checkpoint

- `imaging/track_labels.py`: the tracker of B4 on label stacks (memory, cost instead of IoU gate,
  sparse-frame unique links, merge/split by overlap, `parent_track_id` by touching), with `--batch` writing
  `Combined_Results_retracked.csv` next to the old table (old columns kept, `track_id_v11` = old ID).
- `imaging/synthetic_tracking_test.py`: ground truth with jumps, buds, dropouts, false splits, merged masks,
  transient cells and an empty frame. Two seeds: links across a 2–3 frame gap recovered 0.90–0.92 versus
  0.06–0.08 for the v11-like tracker; tracks per object 2.7–2.9 versus 6.4–6.7; bud parent correct in
  0.95–1.00 versus none; identity switches 1–2 versus 0–1.
- `imaging/sweep_segmentation.py` and `imaging/slurm/sweep.sbatch`: Phase A, scored by consistency,
  merge/split, large-object splits, count stability, time per frame, with overlays.
- `imaging/probe_env.py`, `imaging/slurm/probe.sbatch`, `imaging/slurm/make_manifest.py`,
  `imaging/slurm/track.sbatch`: cluster probe and array jobs for the re-tracking of all experiments.
- `analyse_pipeline`: `AUREO_RESULTS_PATTERN="Combined_Results_retracked.*"` switches the analysis to the
  re-tracked tables (own cache file).

## 6b. First cluster results (probe, re-tracking of WT/Glc/0.75, mini sweep)

- Node wohlrose: NVIDIA L40S 48 GB (plus L4 and a second L40S), 168 CPUs, no time limits; Cellpose 4.1.1,
  torch 2.12 cu130, zarr 3.2.1, nd2 0.11.3. Cellpose `cpsam_v2` takes 0.7–1.4 s per frame at 1139 × 1041
  px, so one 133-frame movie is 2–3 minutes and the whole data set 20–30 GPU-hours on one card.
- **Pixel size 0.0733 µm/px** (nd2 metadata). The chamber crop is 83 × 76 µm. The median object of the
  QC batch, 3,730 px², is 20 µm² (equivalent diameter 5.1 µm), the largest 29,000 px² are 160 µm². In
  Rensink's terms a blastoconidium of 9–11 × 3–6.5 µm is 20–56 µm² (3,700–10,500 px²) and a swollen cell of
  12 × 9 to 15 × 11 µm is 85–130 µm² (16,000–24,000 px²); the size classes are calibrated on our own
  distribution, since Cellpose masks exclude the phase halo.
- Re-tracking of WT/Glc/0.75 (11 chambers, 3–61 objects per frame) from the existing masks: losses of an
  object that is still in the table fell from 31 % of new IDs (v11 tracker on the QC batch) to 2.5 %.
  The remaining new IDs are small and short-lived: median 714 px² (3.8 µm²), 67 % below 1,000 px²
  (5.4 µm²) while a blastoconidium is at least 20 µm²; 47 % appear isolated, in the sparsest frames 78 %.
  Debris and halo fragments passing through the chamber under flow, and buds at their first frames.
  Consequence for the analysis: an object counts as a cell above a size floor and after two frames; the
  table keeps everything.
- A flaw found on the way and fixed: a track absorbed into a neighbour's mask was ended, so every
  merge-and-split flicker of a mother-and-bud mask produced a new ID. Absorbed tracks now stay linkable
  for `memory_merged` frames while the host lives. Synthetic test with 12 % merged masks: tracks per
  object 2.6–2.7 versus 8.2–8.6 for the v11-like tracker.
- Mini sweep (two settings on a movie with 2–3 cells): identical results, `niter` makes no difference
  there. The sweep needs a chamber with tens of cells to say anything.

## 6c. Checkpoint 2: the QC batch re-tracked, the full sweep on one dense chamber

Re-tracking of WT/pH/6 from the existing masks (11 chambers, means):

| | v11 tracker | `track_labels.py` |
| --- | --- | --- |
| new IDs per object-frame | 0.137 | 0.079 |
| median track length (frames) | 3 | 7 |
| tracks of at least 10 frames | 21 % | 43 % |
| tracks per chamber | 824 | 461 |
| manual merges reproduced, all 244 | | 22 % |
| manual merges reproduced, gap 1 to 3 frames (116) | | 34 % |

Gap links 2,846, merges 1,340, splits 1,417, new objects touching a tracked mask 1,612. The manual
merges that stay unreproduced are the large jumps (median 2.2 radii across one frame) and the long
absences (92 of 244 links span more than five frames); a linker without the images cannot make those
safely, and the centroid-only prototype reached the same 23 %. Whether that matters for the result is
decided at checkpoint 3 by `validate_lineage.py` on the re-tracked table, not by the link count.

Sweep on PosCtrl_Rep1_ChamA1, frames 100 to 112 (dense), 36 settings:

- `cpsam_v2` and `cpsam` give identical masks: Cellpose 4.1.1 resolves both names to the same default
  model (the probe log says so). `niter` 0 versus 500 changes nothing measurable.
- `flow_threshold` 0.4 beats 0.8 on every metric at the same `cellprob`: new IDs −27 %, merges −13 %,
  splits −10 %, large-object splits −12 %. The setting chosen to protect the septated swollen cells does
  not protect them better.
- `cellprob_threshold` is a trade-off: −1 gives the fewest new IDs (0.044) but the most merges (7.6 per
  frame); +1 the fewest merges and splits (5.3 and 24.6 per frame, large splits −22 %) but more new IDs
  (0.058). 0 sits between (0.051, 7.3, 30.2).
- Candidates for the by-eye check: (0.4, −1), (0.4, 0), (0.4, +1) against the current (0.8, 0), which
  (0.4, 0) dominates on all metrics. The sweep is one chamber and 13 frames; two more movies with the
  reduced grid (flow × cellprob, `niter` 500, `cpsam`) come before a decision.

## 6d. Checkpoint 3: the analysis on the re-tracked QC batch

Built (analysis side, `analyse_pipeline/`):

- `cell_filter.py`: a track is a cell with at least `CELL_MIN_FRAMES` = 2 frames and a largest area of at
  least `CELL_MIN_MAX_AREA_PX` = 1,500 px² (8 µm²). On the QC batch it removes 31–33 % of the tracks and
  3.5–6 % of the object-frames, for v11 and re-tracked tables alike; report `00_cell_filter.csv`.
- `lineage.classify_mother_bud_measured()`: on re-tracked tables the mother is the track whose mask the
  bud's first mask touches (`parent_track_id`, `link_type` new_touching or split), with the same rules as
  the heuristic (mother established, bud persists `bud_min_frames` = 2, size criterion). Hybrid: candidates
  whose mask touches nothing go through the radius heuristic; `method` says which rule made the event.
- `qc_exclusions.translate_exclusions_to_retracked()`: the manual QC file is mapped through
  `track_id_v11` onto the new IDs, frame-exact; 503 rows became 594, 66 merges were already made by the
  tracker, 13 old tracks were no longer in the data. Written to `qc_exclusions_retracked.csv`.
- No gap closing in the analysis on re-tracked tables; `AUREO_RESULTS_PATTERN` selects the tables.

Results on the QC batch, sparse window, after manual QC:

| | v11 tables, radius heuristic | re-tracked tables, hybrid |
| --- | --- | --- |
| median track length after QC (frames) | 6 | 12 |
| new tracks per object-frame | 0.093 | 0.055 |
| candidates (new tracks in the window) | 250 | 252 |
| events | 102 | 116 (70 measured, 46 radius) |
| the same bud in both sets | 77 | 77 |
| bud growth over 3 frames, median | ×4.8 | ×3.2 |
| contact ratio, median | 1.08 | 1.08 |
| measured mother = nearest established cell | | 75 % |

The two methods agree on three quarters of the events and on the bud signatures. They do not agree on
the chambers: per chamber the counts are 2–16 events, mean 10.5, and their spread (CV 0.29–0.40) is the
Poisson noise of such counts (CV 0.31). Spearman between the methods' chamber rates is −0.08. A budding
rate per chamber in the window is therefore a noise-limited number; it becomes a readout only pooled
over the chambers of a structure, and the story's chip-level treatment already does that.

Earlier note, now done: the QC batch re-tracked on the cluster (`track.sbatch` on one experiment), scored with
`analyse_pipeline/diagnose_tracking.py` and `validate_lineage.py` against the manual QC, plus the sweep
table and overlays. Pipeline v12 (flags instead of drops, raw label stacks, `--file`) follows once the
setting is chosen.

## 7. Order and checkpoints

Phase A (sweep) and the re-tracking on the existing zarr stacks run in parallel on the cluster. B1–B3,
B5–B6 (pipeline v12) after checkpoint 1. Phase C after checkpoint 2. The full re-run last.
