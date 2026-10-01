# Replacement text for the methods chapter: Section 2.3 (growth parameters) and Section 2.5

Draft of 2026-09-29. The current Sections 2.5.1 to 2.5.6 describe the earlier pipeline (manual quality control,
the track-merging routine, the distance heuristic with 10 and 7 frames and 30 px as the only mother-bud rule,
µ_event as growth rate, Kruskal-Wallis as the main test). The results chapter rests on the final pipeline, which
the text below describes in the style of the existing chapter. The paragraph for Section 2.3 is an addition at the
end of the BioLector paragraph (after "... absolute OD₆₀₀ values."); Section 2.5 is replaced as a whole. Sentences
in *[italic square brackets]* are notes to the author, not thesis text. Parameter values are those of
`analyse_pipeline/config.py`, `imaging/pipeline_template_v12.yaml` and `imaging/track_labels.py` of the final run.

---

## Addition to 2.3 Biosensor Details and Validation

Growth parameters were estimated from the scattered-light signal of every well with a custom Python script. The
constant background of medium and plate was estimated per well as the offset for which the early part of the
curve, up to half of its total increase, was most linear on a logarithmic scale, and was subtracted. The maximum
specific growth rate µ_max was determined as the steepest slope of the natural logarithm of the corrected signal
in a sliding window of 2 h (at least five readings, coefficient of determination R² ≥ 0.95, signal above 5 % of
the range of the well), the doubling time as ln 2 / µ_max, the lag time by the tangent method as the intersection
of the µ_max tangent with the logarithm of the initial signal (mean of the first three readings), and y_max as the
maximum of the corrected signal. The parameters of the five wells of a condition were averaged (mean ± SD), and
each biosensor strain was compared with the wild type in the same medium by Welch's t-test. The five wells of a
condition were wells of one plate inoculated from one preculture and represent technical replicates.

---

## 2.5 Data Processing, Image Analysis, and Statistics

### 2.5.1 Image Preprocessing, Segmentation, and Tracking

Automated image processing and single-cell tracking were performed using a custom Python pipeline. All
processing steps were controlled via hierarchical YAML configuration files to ensure full reproducibility and
parameter transparency. Raw time-lapse data were loaded lazily to limit memory usage, and the microfluidic
cultivation chambers were detected and cropped automatically using iterative thresholding and morphological
filtering. The pixel size (0.0733 µm) was read from the metadata of the .nd2 files.

Single-cell segmentation was performed on the phase-contrast images with the Cellpose framework using the
Segment Anything Model architecture (Cellpose-SAM, model `cpsam`, flow threshold 0.4, cell probability threshold
0). Every segmented object was retained, and the morphological limits (contact with the chamber border, area
below or above the admissible range, low solidity, high eccentricity) were recorded as flag columns instead of
being applied before tracking, so that the decision to exclude an object was taken in the analysis. For every
object the projected area, centroid, solidity, eccentricity, perimeter, circularity and axis lengths were
measured, together with the mean and standard deviation of the phase-contrast signal inside the mask. Fluorescence
intensities were extracted from the segmented masks after a rolling-ball background subtraction of the
measurement channels to account for spatially varying background and uneven illumination.

Tracking was performed on the stored label stacks as a separate step (Blöbaum et al., 2024 for the cultivation
system; tracking after the design described in `docs/tracking_diagnosis.md` *[replace by a citation or drop]*).
Objects of consecutive frames were assigned by linear assignment minimising a cost composed of the centroid
distance in units of the cell radius, the logarithm of the area ratio and the number of skipped frames, credited by
the overlap of the two masks; assignments were admitted only within a distance of 1.5 radii per square root of the
gap plus 10 px and an area ratio of at most 2.5. A track that was not detected in a frame remained linkable for
three frames (15 frames for a track that had been absorbed into the mask of a neighbouring cell, as long as that
neighbour was tracked). In sparsely populated frames (at most 20 objects) a second pass linked remaining pairs over
a larger radius if the pair was unambiguous. Merging and splitting of masks were detected from the overlap with the
previous frame: when one mask covered two predecessors the larger track was continued and the smaller one was
marked as absorbed; when a mask fell apart again the separated piece was returned to the absorbed track or, if
none existed, received a new identity with the host as parent. A newly appearing object whose mask, dilated by
two pixels, touched a tracked mask received that track as its parent. For every object the type of link
(continued, gap, long-range, unmerge, split, new touching a mask, new without contact) and the parent track were
stored. The per-movie tables were merged into one table per experiment.

### 2.5.2 Object Filtering and Quality Control

Objects flagged as touching the chamber border or as lying outside the admissible area range were tracked but
excluded from every measurement. On the level of tracks, an object was counted as a cell only if its track
comprised at least two frames and its largest area was at least 1,500 px² (8 µm²); these limits were calibrated
on the manually curated experiment WT/pH/6 min, in which single-frame objects had a median area of 580 px²
whereas 95 % of the tracks of at least five frames exceeded 2,100 px². In addition, dead cells and cell debris,
which lose their phase contrast, were removed by the ratio of the standard deviation to the mean of the
phase-contrast signal inside the mask: a track whose median ratio was below 45 % of the median of the
size-passing tracks of the same structure was excluded. The threshold was defined relative to the structure
because the exposure differed between recordings. The manual curation of the earlier pipeline version, which had
been performed for one experiment (503 corrections), was not applied to the final tables; it served to calibrate
the tracker and the filters. One static recording that duplicated another stage position was excluded.

### 2.5.3 Experimental Units and Aggregation

The curated single-cell data were organised by biosensor strain, oscillation type and cycle period. One
microfluidic structure carried one cycle period with five oscillation chambers, three constant-feast (PosCtrl) and
three to four constant-famine (NegCtrl) chambers distributed over several arrays; the array index and the chamber
position were taken from the file names. On every array the oscillation chambers occupied the positions A3 to
A12, the feast controls A1 and A2 and the famine controls A13 and A14. A physical chip carried two or three structures, was inoculated from
one preculture and was imaged on one day; it is referred to as a culture. For the static cultivations, chip W109
carried three structures with four chambers per medium and chip W65 one structure with five chambers in complex
and four in minimal medium. All readouts were aggregated hierarchically from cells to chambers (mean per chamber),
from chambers to structures (mean over chambers, with the standard deviation and standard error over chambers)
and from structures to conditions. Because the chambers of a structure were technical replicates of one culture,
every comparison between cycle periods was a comparison between structures with n equal to the number of periods,
and every error bar of the oscillation data is a chamber error bar of one structure. Comparisons between the two
media of a static chip were made over the chambers of that chip.

### 2.5.4 Sparse-Phase Window and Lineage Classification

With rising cell density the number of newly appearing objects in a chamber rose in proportion to the number of
objects, and most of them touched an established cell at their first detection, because the masks of neighbouring
cells touched, split and re-merged from frame to frame (Section 3.2). Lineage readouts were therefore evaluated
only in the sparse phase of every chamber, defined as the frames before the rolling median (five frames) of the
number of objects per frame exceeded 20. Chambers with fewer than 20 such frames were excluded from the lineage
readouts. All other readouts used the entire cultivation.

Within the window, every newly appearing object was a bud candidate. If the object had a parent track from the
mask contact at its first detection, the parent was accepted as the mother when it had been tracked for at least
three frames and the candidate persisted for at least two frames. Candidates without a touching mask were assigned
by a distance heuristic: the mother was the nearest cell tracked for at least ten frames within 30 px, and the
candidate had to persist for at least two and at most seven frames as a separate object. The assignment method was
recorded for every event. Washed-in cells and masks that had split in two were removed by a size criterion: the
ratio of the candidate's area to the mother's area at first detection was bimodal over all candidates, and
candidates above the antimode of this distribution (0.32) were rejected. The threshold was determined once from
all data and applied to all strains and conditions. From the accepted events the budding ratio of every mother
(fraction of its observed frames spent in the budding phase) and the budding rate per mother-hour of every chamber
(accepted events per hour of tracked mother time in the window) were calculated. The lineage assignment was
validated by the share of candidates that could be assigned to a mother and the share of ambiguous candidates per
chamber, compared across the structures of a series.

### 2.5.5 Growth, Morphology, and Sensor Readouts

The endpoint of a chamber was the mean over its cells in the last 25 % of its frames. For the static cultivations,
in which the complex-medium chambers filled up early, an absolute window common to all chambers was used instead,
ending at the earliest frame at which the smoothed cell number of any chamber that had grown at least 1.5-fold
reached 90 % of its maximum, and spanning the last 25 % of the frames before it. Endpoint readouts were the
projected cell area, the eccentricity and, for the biosensor strains, the ratiometric sensor signal (intensity of
the sensor channel divided by the reference channel per cell). Sensor readouts excluded the first 120 min of every
cultivation, the period before the environmental perturbations started.

Three specific growth rates were computed. The area-based growth rate µ_area of a single cell was the slope of
the natural logarithm of its projected area against time for tracks of at least ten frames, fitted only before the
first breakpoint of the track (a sudden drop or jump of the area or a loss of solidity, which indicate a division,
a cluster or a segmentation error); mother cells and buds were fitted separately, and the share of reliable fits
was reported per condition. The population birth rate µ_bud of a chamber was the number of accepted budding events
in the sparse-phase window divided by the cell-hours of the window (the number of cells in every frame summed over
the frames and multiplied by the frame interval); in balanced growth every birth adds one cell, so births per
cell-hour are the specific growth rate, and no interval between events is required. Next to it, the immigration
rate was the number of tracks appearing without a parent mask after the first frame of the window per cell-hour.
The event-based growth rate µ_event from the interval between two buddings of the same mother (Blöbaum et al.,
2024) was computed for comparison but not used as a growth rate, because the generation time of the daughters
does not enter it.

The robustness metrics R(t) and R(p) were adapted from Trivellin et al. (2022) and Blöbaum et al. (2024).
Temporal robustness R(t) quantifies the stability of a readout over time, at the population level as the
variation of the chamber mean over the frames and at the single-cell level as the variation of a cell's value
over its frames (cells tracked for at least ten frames); population robustness R(p) quantifies the cell-to-cell heterogeneity of a readout within a
chamber at each time point. Both were applied to the cell area, the eccentricity and the sensor ratios, R(p) also
to µ_area. The robustness of budding was quantified at the single-cell level as R(t) of the event-based rate
µ_event over the budding intervals of a mother (mothers with at least three intervals) and as R(p) of the
budding rate per mother across the mothers of a chamber. Every robustness metric was averaged per chamber and
evaluated against the cycle period and against the controls of the structure in the same way as the readouts
(Section 2.5.6). For the pullulan knockout strain, the agreement of the control chambers within a structure was
quantified as the coefficient of variation of the chamber-median area over the three feast (or famine) chambers,
as the temporal coefficient of variation of the frame median after removal of the linear time trend, and as the
ratio of the feast to the famine level, and was compared with the distribution of the same quantities over the
structures of the pullulan-producing strains.

### 2.5.6 Statistical Analysis

All aggregations, tests and figures were produced within the Python pipeline (pandas, NumPy, SciPy,
Matplotlib). Trends against the cycle period were quantified by the Spearman rank correlation ρ over the
structures of a series, with n equal to the number of periods (four to six); ρ is reported as an effect size, since
with six periods a |ρ| of 0.83 is required for p < 0.05. To separate an effect of the period from an effect of
the structure, the constant-medium control chambers of the same structures, which cannot respond to the period,
were evaluated in the same way. For every readout and series the correlation with the period was computed for the
oscillation chambers, for the feast and famine controls and their mean, and for the difference between the
oscillation chambers and the control mean, and the combination was classified: a period effect required that the
oscillation chambers trended (|ρ| ≥ 0.6), that no control type trended in the same direction, and that the
difference trended as well; a structure effect was recorded when a control type trended in the same direction as
the oscillation chambers; a trend that vanished after subtraction of the controls was recorded as not robust. In
addition, every oscillation value was expressed relative to the two controls of its structure as a bracket score
(0 = famine control, 1 = feast control); the bracket was regarded as degenerate when the two controls differed by
less than twice the chamber standard deviation of the controls, and degenerate structures were excluded from the
trend test of the score. Within cultures that carried two or three periods, the change from the shortest to the
longest period was recorded for the oscillation chambers, their controls and the difference.

The oscillation chambers of a structure were compared with the mean of its two controls, paired over the
structures of all series, by the Wilcoxon signed-rank test, and the ratio of the two was reported per oscillation
type and strain; the robustness metrics were classified against the period in the same way as the readouts.
Feast and famine controls were compared over structures by the Wilcoxon signed-rank test. The two media of a
static chip and the growth parameters of the BioLector cultivations were compared by Welch's t-test over chambers
and wells, respectively. The consistency of the control chambers across the structures of a series was tested by
the Kruskal-Wallis H-test on chamber-level values, and the assignment rate of the lineage classification across
structures likewise. No correction for multiple testing was applied; the number of tested combinations is
reported with the results. Static control conditions were processed through the identical pipeline but evaluated
independently.

*[The current Sections 2.5.2 ("track-merging routine", manual curation of out-of-focus and overlapping cells)
and 2.5.3 (mother ≥ 10 frames, bud ≤ 7 frames, 30 px as the sole rule) are superseded by 2.5.2 and 2.5.4 above;
the description of the QUEEN-2m normalisation in the current 2.5.5 is kept in 2.5.5 above.]*
