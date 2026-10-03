# The pipeline, step by step

An overview for explaining the analysis to a supervisor. Written 2026-10-01 for the state of the final v12 run.
Details and the reasoning behind every rule are in `docs/data_story.md` (the argument), `docs/tracking_plan.md`
(the history of the tracking rebuild) and `README.md` (every module and output prefix). Parameter values are
those of `analyse_pipeline/config.py`, `imaging/pipeline_template_v12.yaml` and `imaging/track_labels.py`.

The pipeline has two halves. The imaging half runs on the cluster, once per movie, and turns a time-lapse into a
table of segmented, tracked objects. The analysis half runs on one machine, reads all those tables at once and
turns them into the readouts, tests and figures of the thesis. Nothing is decided by hand in between.

---

## Part A: from a movie to an object table (cluster, one job per movie)

**A1. Input.** One .nd2 file per chamber: phase contrast plus the fluorescence channels of the strain, one frame
every 10 min, 133 frames in a standard 22 h run. The chamber is detected and cropped automatically (iterative
thresholding, morphological filtering); the pixel size (0.0733 µm) is read from the metadata.

**A2. Segmentation.** Cellpose-SAM (`cpsam`, flow threshold 0.4, cell probability 0) on every phase-contrast
frame. Every segmented object is kept. Whether an object touches the chamber border or lies outside the area,
solidity or eccentricity limits is written as a flag column, not applied; the analysis decides what to exclude.

**A3. Measurements.** Per object: area, centroid, solidity, eccentricity, perimeter, circularity, axis lengths,
the mean and standard deviation of the phase-contrast signal inside the mask (used later to recognise dead
cells), and the fluorescence intensities per channel after a rolling-ball background subtraction.

**A4. Tracking** (`imaging/track_labels.py`, on the stored label stacks, so it can be rerun without
re-segmenting). Objects of consecutive frames are matched by linear assignment on a cost made of centroid
distance in cell radii, log area ratio and skipped frames, credited by mask overlap; gates: 1.5 radii per square
root of the gap plus 10 px, area ratio at most 2.5. A track that is missing in a frame stays linkable for 3
frames (15 frames if it was absorbed into a neighbour's mask). In sparse frames (at most 20 objects) a second
pass links unambiguous pairs over a larger radius. Merges and splits are read from the overlap with the previous
frame. A new object whose mask, dilated by 2 px, touches a tracked mask gets that track as its parent. Every
row records how it was linked (`link_type`: continued, gap, long_range, unmerge, split, new_touching, new) and
its parent track.

**A5. One table per experiment.** `merge_results.py` joins the per-movie tables into one
`Combined_Results.csv` per strain, oscillation type and period (`03_results_v12/`), with QC overlays in
`03_results_v12/QC/`.

---

## Part B: the analysis run (`analyse_pipeline/run_analysis.py`)

### B1. Loading and the experimental units

- All `Combined_Results.csv` of the data tree are read into one table (a parquet cache beside the output
  folder makes the second run fast). Strain, oscillation type and period come from the folder path; chamber
  and array index from the file name. Every row gets a chamber id and a global cell id.
- Rows flagged at the border or outside the area limits are dropped from every measurement. One static recording
  that duplicated another stage position (W65 minimal medium, Rep5) is dropped by name.
- Sensor ratios are computed per cell (sensor channel over reference channel).
- The experimental units are derived: a **structure** is one period with its 5 oscillation chambers and its 3
  feast (PosCtrl) and 3 to 4 famine (NegCtrl) chambers; a **culture** is the physical chip, inoculated from one
  preculture on one day, carrying 2 or 3 structures. The pipeline calls a structure `chip` in its tables. For
  the static data W109 is one chip with three structures and W65 one chip with one. Tables
  `00_chip_overview.csv` and `00_chip_run_order.csv` (were the periods run in date order? they were, in day
  blocks) record this.

### B2. The cell filter (replaces the manual quality control)

A track counts as a cell if it has at least 2 frames and a largest area of at least 1,500 px² (8 µm²). That
removes dirt, halo pieces and single-frame flicker (35 % of the tracks, 7 % of the object-frames). A second
rule removes dead cells and debris, which lose their phase contrast: a track whose median ratio of
phase-contrast standard deviation to mean is below 45 % of the median of the size-passing tracks of the same
structure is dropped. The threshold is relative because the recordings differ in exposure. `00_cell_filter.csv`
reports per chamber what went.

### B3. The sparse-phase window

Objects per frame are counted per chamber. The window is the set of frames before the five-frame rolling
median exceeds 20 objects; chambers with fewer than 20 such frames have no window. Reason: in a dense chamber
every new object touches an established cell and most of them are fragments of touching masks, not buds
(`00_new_objects_vs_density.pdf`, `00_track_fragmentation.csv`, `20_lineage_window.pdf`). Everything that
needs a mother-bud assignment is evaluated only inside the window. Endpoint morphology, µ_area and the
robustness metrics use the whole run.

### B4. Mother and bud

Inside the window every newly appearing object is a bud candidate. If the tracker gave it a parent from mask
contact, that parent is the mother, provided the parent had been tracked for at least 3 frames and the
candidate persists for at least 2 frames (74 % of the events). Candidates without a touching mask go through a
distance heuristic (mother within 30 px, tracked for at least 10 frames, candidate persists 2 to 7 frames as a
separate object). A size criterion removes washed-in cells and split masks: the ratio of candidate area to
mother area at first detection is bimodal over all candidates, and candidates above the antimode (0.32) are
rejected. The threshold is determined once from all data and applied everywhere. Output: `20_budding_events.csv`.

### B5. The analysis steps (each writes tables and figures with its prefix)

| step | what it computes | main outputs |
| --- | --- | --- |
| 00 overview | chambers, frames, objects per frame, tracks; density figure; fragmentation table | `00_*` |
| 10 growth | cell area over time; µ_event per mother from interbud intervals (reported, not used as a growth rate); µ_area per track as the slope of ln(area), tracks of at least 10 frames, fitted before the first breakpoint; per-chip summaries, trend and control-trend tests for µ_area | `10_*`, `11_*`, `12_*` |
| 13 endpoint | mean over the last 25 % of a chamber's frames for area, eccentricity and sensor ratios (static: an absolute window ending at the earliest saturation); chamber, structure and condition levels; Spearman against the period; bracket score; control-trend classification; within-culture change; figures against the period with the controls | `13_*` |
| 20 lineage | budding events; budding ratio per mother and per chamber; buds per mother-hour in the window with the same trend and control tests; Panel A violins; budding-ratio time series; lineage depth; µ_bud (births per cell-hour) and immigration (new tracks without parent per cell-hour); µ_bud against µ_area | `20_*` to `24_*` |
| 30 sensors | channel intensities and ratios over time, after the 2 h preconditioning, as drift over hours (cycles lie below the sampling limit) | `30_*`, `31_*` |
| 40 robustness | R(t) at population and single-cell level and R(p) for area, eccentricity and the sensor ratios; R(p) of µ_area; R(t) of µ_event per mother; R(p) of the budding rate per mother; each metric per chamber, then the same chip logic, tests and figures as the readouts; control consistency across structures (Kruskal-Wallis) | `40_*` |
| 50 summary | summary statistics; control-trend summary over all readouts (one point per readout and series, oscillation ρ against control ρ); the same for the robustness metrics; the paired comparison of oscillation chambers with the controls of their own structure for all readouts and metrics | `50_*`, `51_*` |
| 90 morphology and appendix | morphology per cell (mean area against mean eccentricity, one point per cell, guide lines for large round cells, share in the title), single-cell trajectories, one stable mother per group with its budding marks | `90_*` to `92_*` |
| 95 sensor controls | feast against famine control chambers per sensor, per structure and over time (does a sensor respond to feast versus famine at all?) | `95_*` |

### B6. The three branches and the validation

- **Oscillation branch** (`analysis_output_v12/`): all steps, x axis = period, one structure per period.
- **Static branch** (`static/`): the same steps with x axis = medium (complex versus minimal) and one panel per
  chip family; Welch's t-test over the chambers of a chip; endpoint in a common window before the complex-medium
  chambers saturate.
- **PKO branch** (`pko/`): the pullulan-knockout structure against the producer structures, in addition the
  within-structure agreement of the control chambers (`60_*`, `61_*`), the test of the clogging hypothesis.
- **no_qc branch** and `validate_lineage.py`: the first compares with and without the manual curation of the
  old pipeline (now historical), the second measures how many bud candidates get a mother per chamber and
  whether that assignment rate differs between structures (`lineage_validation/`).

---

## Part C: the statistical logic

1. **What is a replicate.** The five oscillation chambers of a period are chambers of one structure on one chip
   from one preculture: technical replicates. Every error bar of the oscillation figures is a chamber error bar.
   A comparison between periods is a comparison between structures, so the Spearman correlation against the
   period has n equal to the number of periods (4 to 6) and is reported as an effect size, not as a test.
2. **The controls are the handle on the structure.** The feast and famine chambers of a structure receive
   constant medium and cannot respond to the period. If they trend with the period, the structure carries the
   trend. For every readout and series the pipeline computes the Spearman of the oscillation chambers, of each
   control type and of the difference, and classifies: no trend; structure effect (a control trends the same
   way); not robust (the trend vanishes after subtracting the controls); period effect (oscillation chambers
   trend, no control does, the difference trends). The thresholds are |ρ| ≥ 0.6.
3. **Bracket score.** Every oscillation value is expressed relative to the two controls of its structure
   (0 = famine, 1 = feast); where the two controls do not separate beyond the chamber scatter the bracket is
   degenerate and marked hollow.
4. **Within cultures.** The change from the shortest to the longest period inside one physical chip, where the
   preculture cannot differ.
5. **Oscillation against constant medium.** The oscillation mean of a structure against the mean of its two
   controls, paired over the 49 structures (Wilcoxon signed-rank), as ratio or difference, per oscillation type
   and strain, and whether the effect depends on the period.
6. **Tests used.** Spearman (trends), Wilcoxon signed-rank (paired over structures), Welch's t-test (static
   chambers, BioLector wells), Kruskal-Wallis (controls across structures, assignment rate). No multiple-testing
   correction; the number of tested combinations is reported.

---

## Part D: what came out, in five lines

- Every apparent dose response against the period (budding, shape, three sensor ratios) is a structure effect:
  the constant-medium controls on the same structures trend the same way. Of 58 readout-by-series combinations,
  2 remain as period effects, which is what chance produces; 23 are structure effects.
- Under glucose oscillation cells bud about a quarter more often, are 14 % larger and grow 12 % slower in area
  than under constant feast or famine on the same structure, in every strain and independent of the period;
  under pH oscillation nothing differs from the controls.
- The robustness metrics say the same: cell size is slightly less stable and more heterogeneous under
  oscillation, the growth rate more homogeneous, the sensor signals unchanged; against the period they behave
  like the readouts.
- In static cultivation on chip W109, cells in complex medium end 3.5 times larger and bud a sixth as often as
  in minimal medium; chip W65 does not repeat the minimal-medium side, largely because of one chamber in which a
  swollen cell released ten blastoconidia at once.
- Rebuilding the tracker halved the fragmentation and removed the budding trends the old tracker had produced:
  the per-series trends of a lineage readout were tracker properties.

---

## Part E: running it

```
# cluster, per experiment (imaging/README.md): segment_all.py -> track_labels.py --batch -> merge_results.py
cd analyse_pipeline
python run_analysis.py          # RESULTS_VERSION = "v12" in config.py -> analysis_output_v12/ with static/ and pko/
python validate_lineage.py      # lineage_validation/
```

`config.py` holds every threshold named above and logs them, with the method caveats, at the start of every
run. The thresholds that matter most: cell filter 2 frames and 1,500 px², contrast rule 0.45 of the structure
median, sparse limit 20 objects and 20 frames, mother established 3 frames, bud persistence 2 frames, size
criterion 0.32, endpoint last 25 % of the frames, µ_area 10 frames, trend threshold |ρ| 0.6.
