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

## 6. Open questions

1. Pixel size and objective, or confirmation that the nd2 metadata carries it.
2. The figure-eight objects that `flow_threshold 0.8` protects: septated swollen cells or mother-and-bud pairs?
3. Do the `02_processed/masks_*.zarr` stacks still exist for all experiments?
4. Cluster: GPU model and memory, wall-time limit, whether `/prj/microfluidic` is mounted on the compute nodes,
   how the Python environment is provided (module or conda), and the installed Cellpose version.
5. The pH values of feast and famine in the pH series.

## 7. Order and checkpoints

Phase A and the tracker on the existing zarr stacks (B4 on old masks) can run in parallel. B1–B3, B5–B6 after
checkpoint 1. Phase C after checkpoint 2. The full re-run last.
