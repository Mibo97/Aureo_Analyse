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

## 6e. Checkpoint 4: full run with the new analysis, the sweep decided

Correction: the per-chamber track counts in `00_cell_filter.csv` (9,063 tracks for WT/pH/6) show that this
full run used the **v11 tables** with the new analysis (cell filter, persistence), not the re-tracked tables.
The re-tracking of all experiments was still running. The comparison below therefore isolates the effect of
the analysis changes; the run on the re-tracked tables follows.

Full run, v11 tables, new analysis (565 chambers):

- Cell filter: median 38 % of the tracks removed per chamber (IQR 31–46 %) but only 7 % of the
  object-frames (`00_cell_filter.csv`). A handful of nearly empty chambers lose most of their
  object-frames because they held little besides debris.
- Chip-level budding rates: 146 condition means, Spearman 0.95 against the v11-based run, mean
  0.36 → 0.30 per mother-hour (persistence and the cell filter remove events). The chip level is robust to
  these analysis changes; whether it is also robust to the tracking is checked on the re-tracked run.
- Control-trend summary: 21 no trend, 19 structure effect, 5 period effect, 3 not robust; 44 of the 48
  verdicts unchanged, ρ_osc between the runs 0.94. The five period-effect rows are area BSG/pH, µ_area
  BSG/Glc and BSG/pH (opposite signs), budding rate BSO/Glc and BSPH/pH, each with the strongest control
  trending the other way or flat. The story of `docs/data_story.md` stands under the new analysis rules.

Sweep on two more movies (about 20 objects per frame, reduced grid), together with the dense chamber:

- Flow threshold 0.4 beats 0.8 on every metric on all three movies.
- Cell probability +1 gives the fewest merges and splits (large-object splits −20 % against 0) and −1
  the fewest new IDs in dense frames (−12 % against 0); 0 sits between and is what the config has.
- Overlays at flow 0.4: clusters segmented cleanly, septated swollen cells one mask each, attached
  daughters separate.

Decision proposed: `flow_threshold 0.4`, `cellprob_threshold 0.0`, `model_type cpsam` (cpsam_v2 is the
same model), `niter 500` or null (no effect). Re-segmenting the existing data with it buys 10–20 % on
the consistency metrics on top of what the re-tracking already gave; the thesis can proceed on the
re-tracked tables while pipeline v12 (tracking before filtering, raw label stacks, single-file command
line, rotation fix, µm per pixel, cell-type features) is built for the re-run and for future
experiments.

## 6f. Pipeline v12 built

`imaging/cellpose_pipeline_v12.py` (segmentation env): one movie per call, ROI filters become the columns
`at_border`, `below_min_area`, `above_max_area`, `low_solidity`, `high_eccentricity`; raw Cellpose labels
saved as `labels_<stem>.zarr`; tracking with `track_labels.py` right after segmentation (`tracks_<stem>.zarr`,
events, summary); rotation on the full frame before chamber detection and on every channel; `um_per_px`,
Cellpose version and settings in every row; perimeter, circularity, axes, phase-contrast mean and std inside
the mask. Output goes to `02_processed_v12/` and `03_results_v12/` next to the v11 folders, which stay
untouched. `pipeline_template_v12.yaml` carries the decided setting; the experiment YAMLs of v11 keep
providing channels and chamber detection. `segment_all.py` runs all movies, resumable, workers spread over
the GPUs; `merge_results.py` builds `Combined_Results.csv`. The analysis selects the folder with
`AUREO_RESULTS_SUBDIR=03_results_v12`, drops flagged rows (`FLAG_EXCLUDE_ROWS`) and uses the measured
parent; the manual QC file, which names v11 IDs, is not applied to v12 tables. Tested end to end here with a
stubbed Cellpose on a synthetic movie.

## 6g. Checkpoint 5: the full run on the v12 tables

547 of 565 chambers came through v12 (18 movies still to segment or merge). Compared with the v11 tables:

| | v11 original | v11, new analysis | v12 |
| --- | --- | --- | --- |
| new tracks per object-frame, median chamber | | 0.093 (QC batch) | 0.070 (Glc 0.080, pH 0.057) |
| tracks per chamber, median | | 178 | 125 |
| object-frames per chamber, median | | 1,235 | 1,168 |
| cell filter: tracks removed / object-frames removed | | 38 % / 7.2 % | 33 % / 5.3 % |
| control-trend verdicts: no trend / structure / period / not robust | 19 / 18 / 7 / 4 | 21 / 19 / 5 / 3 | 23 / 20 / 3 / 2 |
| chip-level budding rate, Spearman against v11 original (146 conditions) | | 0.95 | 0.77 (Glc 0.66, pH 0.79) |
| mean budding rate per mother-hour | 0.36 | 0.30 | 0.31 |

- The three remaining period-effect rows are area BSG/pH, µ_area BSO/Glc and budding rate WT/Glc, each
  a different readout and strain, none recurring. The count fell with every improvement of the tracking,
  7 → 5 → 3.
- The per-series Spearman of the budding rate against the period does not survive the change of tracking:
  across the ten series it correlates at 0.05 between v11 and v12, while the chip-level rates themselves
  correlate at 0.77–0.82. Which series looked like a period trend was decided by the tracker, not by the
  cells. That is the methods argument for treating the budding rate per structure as the readout and its
  trend against the period as noise.
- BSG pH 6 min keeps losing 60–80 % of its object-frames to the cell filter under both segmentations; the
  small objects there are dead-cell debris (confirmed by eye), so the filter does what it should.

Earlier note, now done: the QC batch re-tracked on the cluster (`track.sbatch` on one experiment), scored with
`analyse_pipeline/diagnose_tracking.py` and `validate_lineage.py` against the manual QC, plus the sweep
table and overlays. Pipeline v12 (flags instead of drops, raw label stacks, `--file`) follows once the
setting is chosen.

## 6h. Checkpoint 6: the complete v12 run (565 chambers) and the dead-cell rule

All 565 chambers, the 18 static ones included, are through pipeline v12 and the analysis
(`analysis_output_v12`, `RESULTS_VERSION = "v12"`, no environment variables). The numbers of checkpoint 5
barely moved with the last 18 movies: new tracks per object-frame 0.070 in the median chamber (q10 0.048,
q90 0.111; Glc 0.080, pH 0.057), median track 5 frames, 78 % of the object-frames in tracks of 10 frames or
more; the cell filter removes 34 % of the tracks and 6.2 % of the object-frames (median chamber 5.3 %,
BSG/pH/6 52 %); 11,351 budding events in 530 chambers, 74 % of them with the mother from the mask contact;
control-trend verdicts 23 no trend / 20 structure effect / 3 period effect / 2 not robust. `docs/data_story.md`
is rewritten on these tables; the earlier numbers stay only where they are marked as v11.

Two things changed on the analysis side with the measured parent:

- **The size criterion is active.** In the v11 tables the ratio of bud area to mother area at first detection
  had one broad mode and the criterion was switched off. With the parent from the mask contact the
  distribution is bimodal (modes 0.07 and 0.45, antimode 0.32, valley depth 0.25) and `bud_size.py` applies
  the global threshold: 2,160 of 13,938 candidates (15.5 %) are rejected as mother-sized objects that appear
  next to a mother, washed-in cells and masks split in two. Per strain and oscillation type the valley is too
  shallow in 11 of 12 groups, so the one global threshold stays (`20_bud_size_threshold.csv`).
- **Dead cells and debris lose their phase contrast.** On the two v12 tables uploaded for this purpose,
  BSG/pH/6 (debris-rich by eye) and WT/pH/6, the track median of `phase_std / phase_mean` is bimodal in
  BSG/pH/6 (modes 0.05 and 0.28, valley 0.11–0.16, 71 % of the tracks in the low mode) and unimodal at
  0.2–0.3 in WT/pH/6; objects with a median area of 3,000 px² or more never fall below 0.18 (5 % quantile) in
  either table, although the illumination differs by a factor of 1.8 between the two experiments (the ratio
  cancels it). Most low-contrast objects are already removed by the size rule; what the contrast rule adds
  is the persistent debris that passes it: on BSG/pH/6 44 tracks with 1,432 object-frames, 20 % of the
  object-frames that survive the size rule, round (circularity 0.91), shrinking (area ratio last/first 0.74),
  1,300–1,700 px², tracked for a median of 31 frames, 70 % of them with a touching "bud" that the lineage
  would have counted; concentrated in three chambers (NegCtrl Rep1/Rep2, Osc Rep3 with 27 % of its frames).
  On WT/pH/6 it removes 4 tracks and 30 object-frames (0.1 %). Implemented as a third rule of
  `cell_filter.py`: `CELL_MIN_PHASE_CV = 0.12` (made relative to the structure in 6k), applied only when the phase columns exist, reported per
  chamber in `00_cell_filter.csv` (`n_tracks_removed_by_contrast`, `n_object_frames_removed_by_contrast`).
  A chamber whose objects all sit below the threshold (all dead, or focus lost) drops out entirely, which the
  report makes visible. The v12 outputs quoted in the data story predate the rule; the next analysis run
  applies it. Physically the rule reads the loss of the refractive-index difference of a lysed cell; a dying
  cell that keeps its contrast is not caught, and nothing in these tables distinguishes dead from dormant
  cells that still have contrast.

## 6i. Checkpoint 7: second v12 run (contrast rule, validation), W65/W109, growth rate from budding, restyle

**Second run** (564 chambers: W65 minimal Rep5 dropped, see below; contrast rule active). The contrast
rule removed 1,951 tracks and 13,321 object-frames (1.3 %) over the whole run, but unevenly: 49 chambers
lose more than 5 % of their object-frames and 17 more than 20 %, and the heaviest losses are the W109 static
chambers (minimal medium Rep1 to Rep4: 46, 39, 49 and 68 %; complex medium Rep2 and Rep4: 65 and 53 %), far
beyond anything the two calibration chips showed. The W109 numbers moved accordingly (minimal-medium area
5,489 → 8,199 px², budding rate 0.33 → 0.13 per mother-hour; complex medium 0.05 → 0.009). Open: are those
chambers full of dead cells, or do the W109 movies have a different phase contrast so that live cells fall
below 0.12? That needs the W109 `Combined_Results.csv` tables; until then the static section is provisional.
The control-trend verdicts moved from 23/20/3/2 to 26/18/2/2: area BSG/pH, area WT/pH and µ_area WT/pH went
from ρ −0.8 to −0.4, one rank swap in a four-period series each.

**Validation on v12** (`lineage_validation/`): the detection rate of the mother assignment differs between
the structures of a series in 4 of 10 series (Kruskal p < 0.05; v11: 7 of 10) and trends with the period in
no consistent direction (BSA/Glc −0.77, BSPH/Glc −0.77, BSG/pH −0.80, BSPH/pH +0.80, BSO/Glc +0.60);
median assignment rate 0.67 per chamber, median distance over search radius 0.93, 4 % ambiguous candidates.

**W65 minimal medium**: Rep1 and Rep5 are the same stage position recorded twice (110 and 109 tracks, 37 and
35 objects, the same 51,000 px² cell); Rep5 is dropped through `config.EXCLUDED_CHAMBERS`, four chambers
remain. The overlay of frames 85 to 88 shows the lower of the two swollen cells releasing a ring of ten
blastoconidia within three frames (tracks 3 to 14): synchronous multipolar budding, the thesis illustration
of a swollen cell turning into blastoconidia.

**W109** is one chip with three structures (lab book), so it joins `STATIC_SINGLE_CHIP_FAMILIES`: the four
movies per medium are four distinct chambers of that chip (186, 225, 369 and 410 tracks in minimal medium,
no duplicates), the unit is the chamber and n is one culture for both static families.

**Growth rate from budding** (`growth_from_budding.py`, outputs `24_*`): µ_bud = births per cell-hour in the
sparse window, births = accepted budding events, cell-hours = cells present per frame summed over the window
times 10 min. In balanced growth each birth adds one cell, so the ratio is the specific growth rate of the
population; bursts count every bud and no interval is needed. The interbud table `11_specific_growth_rate`
stays but is an interbud rate, not µ (medians of 1.6 to 2.9 h⁻¹ per condition; its summary table had been
empty because of a pandas groupby default, fixed). Alongside: immigration = new tracks without a parent mask
per cell-hour. Check on the WT/pH/6 v12 table: µ_bud 0.09 to 0.22 h⁻¹ in the oscillation chambers (doubling
3 to 7 h), immigration 0.15 to 0.28 per cell-hour, i.e. washed-in cells arrive faster than cells are born.
Same tables and figure as the budding rate (per chamber, per chip, summary, Spearman, bracket, control trend,
within culture, vs period), a row in `50_control_trend_summary`, and `24_mu_bud_vs_mu_area.pdf`.

**Restyle** (`plot_style.py`, `config.STRAIN_COLORS`): colour = strain everywhere (WT #2a78d6, BSA #eb6834,
BSO #1baf7a, BSG #eda100, BSPH #e87ba4, PKO grey); controls by marker and fill in the strain colour
(oscillation filled circle, feast control filled up-triangle, famine control hollow down-triangle), periods
as a light-to-dark ramp of the strain hue, static media by fill; white background, faint horizontal grid,
two spines, 8 to 9 pt text, dark marker edges (yellow, green and pink are below 3:1 contrast on white),
editable PDF text. The five colours pass the colour-vision check (weakest pair BSG/BSO, ΔE 9.1); the old
control red clashed with BSA orange (ΔE 5.9), hence markers instead of colours for the controls. The
explanatory footers inside the figures are gone; they belong in the captions. Rendered on the synthetic data.

## 6j. Checkpoint 8: third run (W109 one chip, growth rate from budding)

Same tables as the second run (cell filter and events identical), W109 as one chip (two structure groups of
four chambers, unit = chamber) and the `24_` family. Full run: 11,211 births in 49,657 cell-hours over 535
chambers; µ_bud median 0.21 per cell-hour (q10 0.10, q90 0.37), doubling time 3.3 h; immigration median 0.19
per cell-hour and above the birth rate in 37 % of the chambers. Chamber type without effect (oscillation 0.22,
famine control 0.20, feast control 0.19; bracket degenerate on 42 of 49 structures), Glc series 0.22 to 0.30
against pH series 0.14 to 0.18. Control-trend rows for µ_bud: 4 structure effects, 6 no trend, no period
effect; the summary now has 58 rows, 32 / 22 / 2 / 2. µ_bud against µ_area per chip and chamber type:
Spearman −0.4 (n = 146). Static: W109 minimal 0.11, complex 0.01 (the complex-medium chambers keep cells in
only 83, 28, 77 and 19 frames after the contrast rule, two without any event: the W109 question of 6i
stands); W65 0.11 and 0.14. PKO famine controls 0.32 to 0.44, feast controls 0.05 to 0.17.
`docs/data_story.md` carries these numbers (2.8, 4, 5, 7).

The restyled figures were checked on the real data (172 PDFs of this run). Fixed after that check: a bracket
score beyond the axis is drawn as a triangle at the panel edge with its value printed instead of a line leaving
the panel; the Spearman text sits in the panel title, not on the data; the violin figure compares each strain
against WT only (four brackets instead of ten); the point figures and the validation figure show only the
periods present in each facet; the size-criterion figure keeps its group legend below the panels; the sensor
time course keeps its "2 h" mark inside the axes; the sparse-window figure carries a strain legend.

## 6k. The contrast rule made relative (W109 resolved)

The W109 minimal-medium table shows a different recording, not dead cells: cell mean about 480 counts
(background about 455) against 1,500 to 2,300 on every other chip, and with that little signal above the
camera offset the ratio phase_std / phase_mean of every live cell lies at 0.07 to 0.12, unimodal at 0.095 in
all four chambers and at every size (tracks of 3,000 px² or more: median 0.095, 86 % below 0.12). The fixed
threshold of 0.12 removed 59 % of the size-passing object-frames there. W65 static (the other uploaded table)
looks like the calibration chips (reference 0.31, 0.2 % below 0.12).

Rule now: a track is low-contrast when its ratio is below `CELL_MIN_PHASE_CV_REL` = 0.45 times the median of
the size-passing tracks of the same structure (biosensor/osc_type/osc_freq, i.e. one imaging session). On the
four tables: BSG/pH/6 threshold 0.122, the same 44 tracks (20 % of the size-passing frames) as before;
WT/pH/6 0.109, 2 tracks; W65 static 0.139, 13 tracks (0.3 %); W109 0.043, none. The reference and the
effective threshold stand per chamber in `00_cell_filter.csv`; the log still warns when a chamber loses more
than a quarter of its frames, which is where a debris-dominated reference would show. The oscillation series
are unaffected by the change (same tracks removed), so the third-run numbers of the data story stand; only the
W109 static values return to the first-run numbers with the next run.

## 6l. Checkpoint 9: final run with the relative contrast rule

564 chambers. The contrast rule now removes 432 tracks and 0.6 % of all object-frames (fixed threshold:
1,951 and 1.3 %); 6 chambers lose more than 20 % (17 before), W109 nothing (structure references 0.095 and
0.127 against 0.19 to 0.35 elsewhere). W109 static is back at the first-run values: minimal medium 5,489 px²
against complex 19,105 (Welch p 0.004 over four chambers each), budding 0.33 against 0.05 per mother-hour,
µ_bud 0.18 against 0.03 per cell-hour. Because the relative threshold differs from 0.12 wherever a structure's
reference is far from 0.27, a few oscillation rows moved as well: verdicts 31 / 23 / 2 / 2 (was 32 / 22 / 2 / 2),
the two period effects unchanged (µ_area BSO/Glc, budding rate WT/Glc); µ_bud median 0.21 per cell-hour,
immigration 0.19, 11,361 births in 49,764 cell-hours; events 11,361. `docs/data_story.md` carries the final
numbers throughout; the figure fixes of 6j are in the PDFs of this run.

## 7. Order and checkpoints

Phase A (sweep) and the re-tracking on the existing zarr stacks run in parallel on the cluster. B1–B3,
B5–B6 (pipeline v12) after checkpoint 1. Phase C after checkpoint 2. The full re-run last.
