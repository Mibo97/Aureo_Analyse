# 3. Results

Draft of 2026-09-29 for the thesis, updated 2026-10-05. Numbers come from the final v12 run of the analysis pipeline (564 chambers,
`analysis_output_v12/`, `docs/data_story.md`) and from the BioLector cultivations (`growth_parameters_summary.csv`).
Figure and table numbers follow the chapter PDF of 2026-10-05 (Figures 2 to 17, Tables 1 to 6); the remaining
appendix figures keep the A-numbers of the notes table until the appendix is numbered. The pipeline file behind
every figure and table stands in square brackets in its caption and is to be removed before submission. Citation
keys follow the methods chapter (`2024_Blöbaum`). Strain names follow Section 2.3 (BSpH is labelled BSPH in the
pipeline figures). Sentences in *[italic square brackets]* are notes to the author, not thesis text. The robustness results of 3.4
and 3.5 (oscillation against control chambers, R(t) and R(p) against the controls and the period) were computed
from the per-chamber tables of the final run with the pipeline's own functions (`docs/scratch/robust3.py`); from
the next run on they are pipeline outputs of steps 40 and 50 (`40_*_per_chip.csv`, `40_*_control_trend.csv`,
`50_robustness_control_trend_summary.*`, `51_osc_vs_controls*`).

---

This chapter first describes the growth of the wild type and the four biosensor strains in the four media that
formed the poles of the feast/famine oscillations, measured in a conventional microtiter cultivation (3.1). The
microfluidic data set and the limit that the cell density set to single-cell tracking are described next, because
they define which single-cell readouts could be evaluated and over which part of a cultivation (3.2). The
comparison of complex and minimal medium under static conditions follows (3.3), then the oscillation
experiments: growth, morphology and their robustness against the half-cycle period (3.4), the biosensor readouts
(3.5), and the
behaviour of the constant-medium control chambers on the same structures (3.6). The chapter closes with the
pullulan knockout strain (3.7).

## 3.1 Strain performance in the four reservoir media

The two glucose variants of the oMLP medium (50 g/L and 0 g/L glucose, pH 5.5) and the two pH variants (20 g/L
glucose, pH 3.0 and pH 8.0), i.e. the four media that were later switched against each other in the microfluidic
chip (Section 2.1), were used to compare the growth of the wild type and the four biosensor strains in a
BioLector cultivation (Section 2.3). Each strain was cultivated in five wells per medium on one plate from one
preculture; the five wells are therefore technical replicates of one culture. Growth was followed as scattered
light for 24 to 26 h (Figure 2), and the maximum specific growth rate µ_max, the doubling time t_d, the lag time
and the maximum background-corrected signal y_max were estimated for every well (Table 1).

*[Figure 2: the existing `growth_curves.pdf`, one panel per medium, mean of five wells with a ±1 SD band, strain
colours as in all other figures.]*

**Figure 2: Growth of the wild type and the four biosensor strains in the four reservoir media.** Scattered
light (a.u.) over time in a BioLector cultivation at 30 °C, 1000 rpm; oMLP with 50 g/L or 0 g/L glucose at pH 5.5
(top) and oMLP with 20 g/L glucose at pH 3.0 or pH 8.0 (bottom). Lines are the mean of five wells of one plate,
bands ±1 SD. [`growth_curves.pdf`]

**Table 1: Growth parameters in the four reservoir media.** Mean ± SD of five wells per strain and medium.
µ_max: steepest slope of ln(signal) in a sliding 2 h window after background subtraction; t_d = ln 2 / µ_max; lag
time by the tangent method; y_max: maximum background-corrected scattered light. [`growth_parameters_summary.csv`]

| medium | strain | µ_max (h⁻¹) | t_d (h) | lag time (h) | y_max (a.u.) |
| --- | --- | --- | --- | --- | --- |
| oMLP 50 g/L glucose | WT | 0.50 ± 0.01 | 1.4 | 5.8 ± 0.5 | 114 ± 12 |
| | BSA | 0.36 ± 0.02 | 2.0 | 7.5 ± 2.2 | 52 ± 2 |
| | BSO | 0.43 ± 0.05 | 1.6 | 6.6 ± 1.8 | 80 ± 4 |
| | BSG | 0.33 ± 0.03 | 2.1 | 2.8 ± 1.2 | 123 ± 6 |
| | BSpH | 0.45 ± 0.03 | 1.5 | 8.6 ± 0.9 | 133 ± 5 |
| oMLP 0 g/L glucose | WT | 0.09 ± 0.03 | 8.6 | 12.0 ± 1.2 | 59 ± 13 |
| | BSA | 0.02 ± 0.00 | 42 | 9.6 ± 5.6 | 38 ± 6 |
| | BSO | 0.02 ± 0.01 | 60 | 10.3 ± 6.5 | 37 ± 3 |
| | BSG | 0.02 ± 0.00 | 35 | 3.3 ± 0.2 | 39 ± 2 |
| | BSpH | 0.02 ± 0.01 | 36 | 2.9 ± 1.7 | 38 ± 2 |
| oMLP pH 3 | WT | 0.36 ± 0.09 | 2.0 | 11.1 ± 5.2 | 143 ± 16 |
| | BSA | 0.19 ± 0.06 | 4.1 | 19.6 ± 0.6 | 73 ± 17 |
| | BSO | 0.34 ± 0.01 | 2.0 | 6.2 ± 0.4 | 53 ± 2 |
| | BSG | 0.34 ± 0.01 | 2.1 | 3.8 ± 2.8 | 87 ± 18 |
| | BSpH | 0.30 ± 0.02 | 2.3 | 3.9 ± 1.1 | 69 ± 10 |
| oMLP pH 8 | WT | 0.34 ± 0.03 | 2.0 | 3.5 ± 0.7 | 130 ± 4 |
| | BSA | 0.32 ± 0.02 | 2.2 | 7.5 ± 3.7 | 57 ± 1 |
| | BSO | 0.31 ± 0.01 | 2.3 | 4.9 ± 0.2 | 44 ± 2 |
| | BSG | 0.29 ± 0.01 | 2.4 | 2.4 ± 2.1 | 125 ± 7 |
| | BSpH | 0.23 ± 0.02 | 3.1 | 1.8 ± 0.6 | 98 ± 12 |

In oMLP with 50 g/L glucose all five strains grew exponentially after a lag of 3 to 9 h. The wild type reached
the highest maximum specific growth rate (0.50 h⁻¹, doubling time 1.4 h). The four biosensor strains grew more
slowly, BSpH at 90 %, BSO at 85 %, BSA at 71 % and BSG at 66 % of the wild-type rate (Welch's t-test against the
wild type, p ≤ 0.024 for each). The maximum signal did not follow the same order: BSpH (133 a.u.) and BSG
(123 a.u.) reached the wild-type level (114 a.u.) or exceeded it, whereas BSO reached 71 % and BSA 45 % of it
(p ≤ 0.002). The curves were not single exponentials. The wild type rose steeply between 12 and 20 h and reached a
plateau at 22 h, BSG and BSO showed a second rise after 18 h, and BSpH reached its maximum at 20.5 h, after which
the scattered light declined. BSG had the shortest lag (2.8 h) and BSpH the longest (8.6 h).

Without glucose only the wild type increased its signal, from about 32 to 65 a.u. within 24 h, corresponding to a
maximum specific growth rate of 0.09 h⁻¹ after a lag of 12 h. The four biosensor strains stayed at their initial
level (µ_max ≤ 0.02 h⁻¹, y_max 37 to 39 a.u.).

In oMLP with 20 g/L glucose at pH 3.0 the growth rates of all strains were lower than at 50 g/L glucose and
pH 5.5. The wild type grew at 0.36 ± 0.09 h⁻¹ with a large well-to-well variation; BSO, BSG and BSpH (0.30 to
0.34 h⁻¹) were not distinguishable from it (p 0.23 to 0.62), BSA was slower (0.19 h⁻¹, p 0.010). The wild type
reached the highest maximum signal of all cultivations (143 a.u., plateau from 22 h), and every biosensor strain
stayed below it (37 to 61 % of the wild-type value, p ≤ 0.001). The BSA wells started at an elevated signal of
about 75 a.u. that declined over the first 10 h before growth began; the cause is unknown, and the lag time of
19.6 h estimated for BSA in this medium reflects this initial decline rather than a growth delay. A smaller initial
decline was present in the wild-type wells.

At pH 8.0 the wild type grew at 0.34 h⁻¹; BSA (0.32 h⁻¹) and BSO (0.31 h⁻¹) were not distinguishable from it
(p 0.15 and 0.052), BSG (0.29 h⁻¹, p 0.014) and BSpH (0.23 h⁻¹, p < 0.001) were slower. The lag times were the
shortest of the four media (1.8 to 7.5 h), and none of the cultures had reached a plateau after 24 h. BSG matched
the wild-type maximum signal (125 against 130 a.u., p 0.25), BSpH reached 76 % of it, BSA 44 % and BSO 34 %.

Across the four media the ranking of the biosensor strains against the wild type therefore depended on the
medium. BSA had the lowest maximum specific growth rate in three of the four media and the lowest or second
lowest maximum signal in all four; BSO reached wild-type growth rates in the pH media but the lowest maximum
signal at pH 8; BSG and BSpH combined a reduced growth rate with a wild-type-like maximum signal at 50 g/L
glucose, and BSG also at pH 8.

*[The methods say that one well per condition was measured for OD₆₀₀ at the end of the cultivation to calibrate
the scattered light. The figure and the table are in a.u. because the calibration values are not in the uploaded
files; µ_max and the lag time are unaffected by a linear calibration, y_max scales with it. If you send the OD₆₀₀
readings with the corresponding scattered-light values, I add the calibration to the plot script and the y axis
and y_max column change to OD₆₀₀.]*

## 3.2 The microfluidic data set, the cells in it, and the limit set by cell density

### 3.2.1 Size of the data set

The oscillation experiments comprised 50 microfluidic structures, each carrying one switching interval of the
medium with five oscillation chambers and its own constant-feast (PosCtrl) and constant-famine (NegCtrl) chambers, on 20 physical
chips (Table 2). Every structure was one interval of one strain and one oscillation type; a physical chip carried
two or three structures and was inoculated from one preculture on one day. The interval is the half-cycle
period: the medium was switched every 0.75 to 24 min, so a full feast/famine cycle lasted 1.5 to 48 min. It is
called the period of a structure below. The five oscillation chambers of a
period therefore were technical replicates of one culture, and every comparison between periods was a comparison
between structures. The static cultivations comprised 17 chambers on two chips. Images were acquired every
10 min for 22 h (133 frames in the standard run; one chamber broke off after 13 frames), with a pixel size of
0.0733 µm, so that a blastoconidium of 20 µm² covered about 3,700 px².

**Table 2: The microfluidic data set.** One structure carries one switching interval (half-cycle period); a culture is one physical chip
inoculated from one preculture on one day. [`00_chip_overview.csv`]

| strain | oscillation type | half-cycle periods (min) | structures | chambers | cultures |
| --- | --- | --- | --- | --- | --- |
| WT | Glc | 0.75, 1.5, 3, 6, 12 | 5 | 55 | 2 |
| WT | pH | 0.75, 1.5, 6, 24 | 4 | 44 | 2 |
| BSA | Glc | 0.75 to 24 | 6 | 67 | 2 |
| BSA | pH | 0.75, 1.5, 6, 24 | 4 | 44 | 2 |
| BSO | Glc | 0.75 to 24 | 6 | 66 | 2 |
| BSO | pH | 0.75, 1.5, 6, 24 | 4 | 44 | 2 |
| BSG | Glc | 0.75 to 24 | 6 | 66 | 2 |
| BSG | pH | 0.75, 1.5, 6, 24 | 4 | 43 | 2 |
| BSpH | Glc | 0.75 to 24 | 6 | 66 | 1 |
| BSpH | pH | 0.75, 1.5, 6, 24 | 4 | 41 | 2 |
| pullulan KO | Glc | 3 | 1 | 11 | 1 |
| WT, static | oMLP and YPD | none | W109: 3, W65: 1 | 17 | 2 |
| total | | | 54 | 564 | 22 |

After exclusion of the objects at the chamber border, segmentation and tracking yielded 96,720 object tracks
with 1,019,226 object-frames in the 564 chambers. The
cell filter (at least two frames and a largest area of at least 1,500 px², i.e. 8 µm², and a phase contrast of at
least 45 % of the structure median) removed 35 % of the tracks but only 6.8 % of the object-frames (median
chamber 5.7 %), since the removed tracks were mostly single-frame objects of a few µm²; the phase-contrast rule
accounted for 432 tracks and 0.6 % of the object-frames. The debris-rich structure BSG/pH/6 min lost 52 % of its
object-frames. 63,215 cell tracks with 950,203 object-frames remained. The number of cell tracks per chamber
ranged from 30 to 236 (10th to 90th percentile over the 536 chambers of the oscillation experiments, median 83);
feast-control chambers held more tracks (median 121) than oscillation (71) and famine-control chambers (69), the
WT/pH structures the most (up to 333 per chamber on average) and BSG/pH/6 min the fewest (24) (Figure 3).

**Figure 3: Number of cell tracks per chamber.** Cell tracks after the cell filter, mean ± SD over the chambers of
each structure, against the half-cycle period; bars: oscillation chambers, feast control, famine control; top row
Glc, bottom row pH series. [`00_n_tracks_overview_summary.pdf`; numbers from `00_data_overview.csv`]


### 3.2.2 Observed morphology

The cells in the chambers took some of the forms described for *A. pullulans* by (2024_Rensink): ovoid
yeast-like cells, large round swollen cells, pseudohyphae and the blastoconidia that bud from them (Figure 5);
chlamydospores and germ tubes were not observed. Figure 4 shows every tracked cell as one point, its mean
projected area over its track against its mean eccentricity, for the five strains of the oscillation experiments
and for the chips and media of the static cultivations.

The typical cell was an ovoid yeast-like cell. Over the 22,090 cells tracked for at least ten frames in the
oscillation experiments the median mean area was 22 µm² (10th to 90th percentile 11 to 56 µm²) and the median
eccentricity 0.75, which corresponds to an axis ratio of 1.5; 90 % of the cells had an axis ratio below 1.9.
Cells in the Glc series were slightly larger and more elongated than in the pH series (median 24 against
22 µm², eccentricity 0.80 against 0.72). Among the strains the wild type had the smallest cells (median 18 µm²)
and BSpH the largest (31 µm²), with BSA, BSO and BSG in between (24 to 26 µm²; Figure 6).

A minority of the cells were large and round, the swollen cells of (2026_Rensink). They form the second cluster
in Figure 4, separate from the band of yeast-like cells: 5 % of the tracked cells combined a mean area of at
least 30 µm² with an eccentricity below 0.6 (WT 3 %, BSA 5 %, BSO 6 %, BSG 8 %, BSpH 9 %), 34 % reached
30 µm² and 1.2 % exceeded 100 µm². In the static cultivations the share depended on medium and chip. In complex
medium on chip W109, 82 % of the cells were at least 30 µm² and the median cell measured 82 µm²; in minimal
medium on the same chip the cells were small (median 21 µm²) and the most elongated of all groups (eccentricity
0.85, axis ratio 1.9). On chip W65 the minimal-medium chambers held the largest share of large round cells
(25 % of 51 tracked cells), among them the two cells of about 270 µm² that released blastoconidia (3.3), whereas
its complex-medium chambers resembled the oscillation chambers (median 28 µm², 1 % large and round).

The buds were blastoconidia released from large cells. Accepted buds measured 4.7 µm² at first detection (10th to
90th percentile 2.3 to 9.9 µm²), and their mothers 71 µm² (32 to 120 µm²), more than twice the median cell. This
holds for all bud candidates, not only the accepted ones (median mother 66 µm², 4 % of the mothers below
15 µm²), so it is not a consequence of the size criterion of 3.2.5. Pseudohyphal chains of elongated cells
(Figure 5 A) were seen in complex medium on chip W65 *[your text says "only in the W109 chip in complex medium";
Figure 5 A is from W65: on which chips and media did you see them?]*; the segmentation split them into their
individual cells (among the objects of the four structures with complete shape data, WT/pH/6 min, BSG/pH/6 min,
W65 and W109, about 125,000 objects of at least 8 µm², the median solidity was 0.98, fewer than 0.3 % had an axis
ratio above 3 and 0.1 % were flagged as highly eccentric by the segmentation), so they enter Figure 4 as single
elongated cells, not as filaments.

**Figure 4: Morphology of the tracked cells.** One point per cell tracked for at least ten frames: mean projected
area (log scale) against mean eccentricity, one panel per strain over all oscillation experiments (top row) and
per chip and medium of the static cultivations (bottom row). The guide lines at 30 µm² and eccentricity 0.6 mark
the quadrant of large round cells, the swollen cells of (2026_Rensink); their share stands in the panel title.
[`90_morphology_scatter.pdf`, `static/90_morphology_scatter.pdf`, from the next pipeline run]

**Figure 5: Phase-contrast examples of the observed morphologies.** (A) A pseudohyphal chain of elongated cells,
chip W65, complex medium; (B) two swollen cells with the blastoconidia they released, chip W65, minimal medium,
the chamber described in 3.3 (frame *[n]*). Red outlines: segmentation masks; scale bar 10 µm *[check]*. [from the
QC overlays in `03_results_v12/QC/`]

**Figure 6: Cell area and eccentricity per strain.** Violins over one value per cell (mean over the track, tracks
of at least ten frames), left Glc and right pH series, dot = mean; stars: Mann-Whitney U over structures (median
per structure, each strain against the wild type; ns p ≥ 0.05). [`21_panel_a_violin.pdf`, area and eccentricity
rows, from the next pipeline run; the budding-ratio row of the file is not shown]

### 3.2.3 How well the cells were tracked

Tracking links the segmented objects of consecutive frames into tracks. A new track starts whenever an object
cannot be linked to one of the previous frame: when a bud appears, when a cell is washed in, or when the tracker
loses a cell and restarts it under a new identity. The last case, fragmentation, decides how long the tracks are
on which growth and budding can be measured. Over all chambers, new tracks began at a median rate of 0.070 per
object and frame (10th to 90th percentile over chambers 0.047 to 0.110), about one per 14 object-frames, buds and
washed-in cells included. The median track lasted 5 frames (50 min), but 78 % of all object-frames belonged to
tracks of at least 10 frames, so most of the data lay in long tracks. Structures of the Glc series fragmented more
than those of the pH series (median 0.080 against 0.057 new tracks per object-frame). The previous tracker
(pipeline v11: overlap threshold, no memory) had fragmented twice as much, with 0.14 new tracks per object-frame,
a median track length of 3 frames and 31 % single-frame tracks (`docs/tracking_diagnosis.md`).

### 3.2.4 Why budding events can be counted only in the sparse phase

The chambers were inoculated with one to three cells and filled to as many as 180 objects per frame within 22 h
(Figure 7 A). Whether a new object is a bud is decided from its contact with a tracked cell: the tracker records
for every new object whether its mask touched an existing cell when it was first detected. In a sparse chamber a
new object that touches a cell is a bud emerging from its mother. In a crowded chamber, however, the masks of
neighbouring cells touch, merge and split from one frame to the next, and every fragment of such a mask is also a
new object touching a tracked cell. Figure 7 B quantifies how the two cases scale with density. Each point is one
chamber during one block of 22 frames (3.7 h); its position gives the number of objects in the chamber during
that block on the x axis and the number of new objects appearing per frame on the y axis. Dark points are new
objects that touched an existing cell, light points new objects that touched nothing. The dark series rises along
a straight line of slope one on the logarithmic axes: the number of touching new objects is proportional to the
number of objects present. The right panel shows the same data divided by the number of objects. From about
8 objects per frame upward every object in the chamber gave rise to a new touching object at a constant rate of
0.04 to 0.06 per frame, i.e. once every 17 to 25 frames, whether the chamber held 10 or 200 objects. The light
series, the washed-in cells, saturated at about two new objects per frame, because the flow, not the density,
delivers them. In the manually curated structure WT/pH/6 min this amounted to 0.5 touching new objects per
block at 2 objects per frame and 130 per block, about six per frame, at 132 objects per frame, while the free new
objects stayed at 40 to 50 per block. Above 20 objects per frame, 55 to 63 % of all new objects touched a mask on
the three structures that reached that density. *[Numbers of Figure 7 B from the four structures examined in
detail, 35 chambers of WT/pH/6 min, BSG/pH/6 min, W65 and W109; to be replaced by those of the full-run figure.]*

At the high densities the touching new objects are therefore dominated by fragments of touching masks, which the
contact rule cannot tell from buds. Three further observations support this reading: the masks visibly touch and
re-split in the dense phase; the manual curation of the previous tracker found a third of its new identities to be
losses of tracked cells; and the budding trends against the period changed with the tracker (3.2.6). Lineage
readouts were therefore restricted to the sparse phase of every chamber, defined as the frames before the rolling
median of objects per frame exceeded 20 (Section 2.5.4). 563 of the 564 chambers had such a window of at least
20 frames; its median length was 117 of 133 frames (range 26 to 136), and 245 chambers never exceeded 20 objects
and were used whole. Endpoint morphology, the area growth rate and the robustness metrics used the whole
cultivation.

**Figure 7: The limit set by cell density.** (A) Objects per frame over time for 80 of the 564 chambers, with the
sparse limit of 20 objects (dashed), and the distribution of the sparse-phase window length over all chambers.
(B) Left: new objects per frame against objects per frame, one point per chamber and block of 22 frames; dark:
new objects that touched an existing cell at first detection (bud candidates or fragments of touching masks),
light: new objects that touched nothing (washed-in cells); lines and bands: median and quartiles per density
class (powers of two); dashed line: sparse limit. Right: the dark series divided by the number of objects, i.e.
touching new objects per existing object and frame; a flat line means proportionality to density.
[`20_lineage_window.pdf`; `00_new_objects_vs_density.pdf` from the next pipeline run, currently drawn for the
four structures WT/pH/6, BSG/pH/6, W65 and W109]

### 3.2.5 Budding events, the size criterion and validation

Within the sparse-phase windows 11,361 budding events were accepted in 531 chambers (median 18 per chamber, 10th
to 90th percentile 8 to 39). For 74 % of them the mother was the cell whose mask the bud touched at first
detection, for 26 % it was assigned by the distance heuristic because the bud touched no mask; 72 % of the mothers
were tracked for 30 frames or longer. The ratio of bud area to mother area at first detection was bimodal, with
modes at 0.07 and 0.45 and an antimode at 0.32 (valley depth 0.25; Figure 8). Objects above 0.32 were
mother-sized objects appearing next to a mother, i.e. washed-in cells and masks that had split in two, and were
rejected: 2,162 of 13,940 candidates (15.5 %). Per strain and oscillation type the valley was too shallow to
support separate thresholds in 10 of 12 groups, so the global threshold was applied throughout.

**Figure 8: Size criterion of the mother-bud assignment.** Area of a new object at first detection relative to the
area of its assigned mother, all 13,940 candidates. Left: all branches pooled, histogram and density on a log
scale, modes at 0.07 (buds) and 0.45 (washed-in cells and split masks), antimode 0.32 (dashed) with the search
range (grey). Right: the density per strain and oscillation type (n = candidates). [`20_bud_size_at_appearance.pdf`]

The share of bud candidates that could be assigned to a mother (assignment rate) was 0.67 per chamber (median),
with 4 % of the candidates ambiguous between two mothers. The assignment rate differed between the structures of
a series in 4 of the 10 strain-by-oscillation-type series (Kruskal-Wallis, p < 0.05) and correlated with the
period in no consistent direction (Spearman ρ −0.77 for BSA/Glc, −0.66 for BSpH/Glc, −0.80 for BSG/pH,
+0.60 for BSO/Glc; Figure 9).

**Figure 9: Assignment rate of the bud candidates.** Fraction of the bud candidates of a chamber that received a
mother, one point per chamber, colour = strain, bar = median per half-cycle period; Glc series (left) and pH
series (right). [`lineage_validation/lv_02_detection_rate.pdf`, without the static panel from the next
validation run]

### 3.2.6 Dependence of the budding trends on the tracking

Rebuilding the segmentation and the tracking changed which series showed a budding-rate trend against the period.
With the previous tracker, eight of the nine series with a non-zero correlation had a negative one (WT/Glc
ρ −0.90, BSpH/Glc −0.89, BSG/Glc and BSO/Glc −0.66). With the rebuilt tracker the same four series gave −0.60,
−0.26, −0.20 and +0.77. Across the ten series the per-series correlations of the two trackings correlated with
each other at ρ 0.05, whereas the budding rates of the individual structures correlated at ρ 0.77. The
morphology and sensor trends, which do not depend on tracking, were the same in both runs (3.4, 3.5).

## 3.3 Static cultivation: complex against minimal medium

The wild type was cultivated without medium switching in complex medium (YPD) and in minimal medium (oMLP, 50 g/L
glucose) on two chips (Section 2.4). Chip W109 carried three structures with four chambers per medium; chip W65
carried five chambers in complex and, after one stage position that had been recorded twice was counted once,
four chambers in minimal medium. Both chips were inoculated from one preculture each, so the comparison within
each chip is a comparison between chambers of one culture. The endpoint was evaluated in an absolute frame window
ending at the earliest saturation of any chamber (frames 58 to 77, i.e. 9.7 to 12.8 h after the start), because
the complex-medium chambers of chip W109 filled up early (Appendix Figure A4).

**Table 3: Complex against minimal medium on two static chips.** Mean over chambers; endpoint area and
eccentricity in the common endpoint window, budding rate per mother-hour and specific birth rate µ_bud in the
sparse-phase window. p: Welch's t-test over the chambers of one chip (four against four on W109, four against
five on W65). [`static/13_endpoint_summary.csv`, `static/21_budding_rate_summary.csv`,
`static/24_growth_from_budding_summary.csv`]

| readout | W109 minimal | W109 complex | p | W65 minimal | W65 complex | p |
| --- | --- | --- | --- | --- | --- | --- |
| endpoint cell area (µm²) | 29.5 | 102.7 | 0.004 | 112.5 | 107.3 | 0.89 |
| endpoint eccentricity | 0.80 | 0.75 | 0.23 | 0.43 | 0.72 | 0.006 |
| buds per mother-hour | 0.33 | 0.05 | 0.12 | 0.14 | 0.18 | 0.66 |
| µ_bud (births per cell-hour) | 0.18 | 0.03 | 0.099 | 0.11 | 0.14 | 0.68 |

On chip W109 cells in complex medium ended 3.5 times larger than in minimal medium (102.7 against 29.5 µm²,
Welch's t-test over four chambers each, p 0.004), were marginally rounder (eccentricity 0.75 against 0.80,
p 0.23), and budded a sixth as often per mother-hour (0.05 against 0.33; the four minimal-medium chambers gave
0.45, 0.58, 0.27 and 0.00 buds per mother-hour, the complex-medium chambers 0.03 to 0.08). The population birth
rate showed the same difference (0.03 against 0.18 births per cell-hour). In complex medium the growth of the
population therefore went into cell size rather than into blastoconidia (Figure 10).

Chip W65 did not repeat the medium effect. Its cells reached the same endpoint area in complex medium as those of
W109 (107.3 against 102.7 µm²), but in minimal medium they were 3.8 times larger than on W109 (112.5 µm²), and the
budding rates of the two media were alike (0.14 and 0.18 buds per mother-hour). The eccentricity differed in the
opposite direction to W109 (0.43 in minimal against 0.72 in complex medium, p 0.006). The four minimal-medium
chambers of W65 were not one population: one chamber ended at 215 µm², the other three at 70 to 86 µm², a
chamber-to-chamber spread of 61 % against 4 to 24 % in the other groups. That chamber contained two swollen
cells of about 270 µm² each, and the lower one released a ring of ten blastoconidia within three frames (frames
85 to 88; Figure 5 B), the transition from swollen cell to blastoconidia described by (2026_Rensink), which set the chamber's
endpoint area and budding rate apart from the other three. With one culture per chip, the medium effect is a
result of chip W109, with W65 as a second single chip that showed the same cell size in complex medium and a
different picture in minimal medium.

**Figure 10: Static cultivation in complex and minimal medium.** (A) Endpoint cell area, (B) buds per mother-hour
in the sparse-phase window and (C) births per cell-hour, per chip family; mean ± SEM over chambers, every chamber
as a small point beside the mean, complex medium filled and minimal medium hollow.
[`static/13_endpoint_vs_medium_area.pdf`, `static/21_budding_rate_vs_medium.pdf`,
`static/24_growth_from_budding_vs_medium.pdf`; the per-chamber points appear with the next pipeline run]

## 3.4 Oscillations: growth, morphology and their robustness against the half-cycle period

The half-cycle periods of 0.75 to 6 min are shorter than the 10-min frame interval, so single cycles are not
resolved there; the 12-min half-cycle (full period 24 min) is at the sampling limit and the 24-min half-cycle
(48 min) is resolved with 4.8 frames per cycle. All readouts below are therefore cumulative over hours and do
not describe the course of a cycle. Robustness under oscillation was read in three ways: from the stability of
the readouts across the half-cycle periods
against the trends of the constant-medium controls, from the comparison of the oscillation chambers with the
controls of their own structure, and from the temporal and population robustness metrics R(t) and R(p)
(Section 2.5.5).

For every oscillation series (strain by oscillation type) one value per structure was compared across the cycle
periods, with the constant-medium control chambers of the same structure plotted as separate markers
(Figure 11, Appendix Figures A5 to A7). Because a structure carried one period, the
Spearman rank correlation of a readout with the period was computed over structures, with n equal to the number
of periods (four to six); with n = 6 a |ρ| of 0.83 corresponds to p 0.05, and ρ is reported as an effect size, not
as a test. In addition, every oscillation value was expressed relative to the two controls of its structure as a
bracket score (0 = as under constant famine, 1 = as under constant feast; `*_bracket_score.csv`), which is
undefined where the two controls did not separate beyond the chamber scatter.

**Figure 11: Readouts of the Glc oscillation series against the half-cycle period.** (A) Buds per mother-hour in the
sparse-phase window, (B) births per cell-hour (µ_bud) and (C) endpoint cell area, one panel per strain; filled
circles: oscillation chambers (mean ± SD over the chambers of the structure), triangles: feast control (filled, up) and famine control
(hollow, down) of the same structure. ρ, p and n in the panel titles: Spearman correlation over structures. The pH series are shown in Appendix Figure A5.
*[Your Figure 11 shows B and C; the budding rate per mother-hour (A) then needs a place in the text or the
appendix, because the paired comparison of 3.4 rests on it.]*
[`21_budding_rate_vs_period_Glc.pdf`, `24_growth_from_budding_vs_period_Glc.pdf`,
`13_endpoint_vs_period_area_Glc.pdf`]

**Levels.** The endpoint cell area of the oscillation chambers lay between 21 and 59 µm² across all series (wild
type 21 to 40 µm², BSG/pH 40 to 59 µm²), the eccentricity between 0.57 and 0.82. The budding rate in the
sparse-phase window ranged from 0.22 to 0.97 buds per mother-hour in the Glc series and from 0.08 to 0.53 in the
pH series. The specific birth rate of the population, µ_bud, had a median of 0.21 births per cell-hour over all
535 chambers with a window (10th to 90th percentile 0.10 to 0.37, doubling time 3.3 h); the Glc series ran at
0.22 to 0.30 and the pH series at 0.16 to 0.18 births per cell-hour (medians per series). At the same time cells
entered the chambers with the flow at a median rate of 0.19 per cell-hour (new tracks without a parent mask), and
immigration exceeded births in 37 % of the chambers, so that the object count of a chamber was not a growth
curve (Appendix Figure A6). Per structure and chamber type, µ_bud and the area growth rate of single cells
µ_area were not coupled in one direction (Figure 12): pooled over both oscillation types they
correlated negatively (Spearman ρ −0.41, n = 146), but only because the Glc structures combined a higher µ_bud
with a lower µ_area than the pH structures (medians 0.25 against 0.16 births per cell-hour and 0.11 against
0.18 h⁻¹); within the Glc series the correlation was −0.24 (n = 87) and within the pH series +0.39 (n = 59).

**Figure 12: Population growth from budding against single-cell area growth.** µ_area (mean over the tracks of
the chambers) against µ_bud, one point per structure and chamber type (circles: oscillation chambers; triangles:
feast control filled, famine control hollow), colour = strain; left Glc, right pH series; dashed line: equal
rates; ρ: Spearman over the points of the panel. [`24_mu_bud_vs_mu_area.pdf`]

**Oscillation against constant medium.** Robustness to the oscillation itself was read from the comparison of
the oscillation chambers with the constant-medium controls of the same structure, paired over structures
(Wilcoxon signed-rank test of the oscillation mean against the mean of the two controls, n = 49). Under glucose
oscillation the cells budded more often than under constant feast or famine: the budding rate of the oscillation
chambers exceeded the control mean on 27 of the 29 Glc structures (median ratio 1.25, p < 0.001), and µ_bud on
23 of 29 (ratio 1.15, p < 0.001). The cells were also larger (endpoint area above the control mean on 25 of 29
structures, ratio 1.14, p < 0.001), marginally rounder (eccentricity ratio 0.98, p < 0.001) and grew more slowly
in area (µ_area below the control mean on 22 of 29 structures, ratio 0.88, p 0.010). Under pH oscillation none
of the readouts differed from the controls (budding rate ratio 1.09, p 0.09; µ_bud 1.05, p 0.12; area 1.00,
p 0.65; µ_area 0.97, p 0.33; eccentricity 0.98, p 0.45). Over both oscillation types the budding rate of the
oscillation chambers lay above both controls of their structure on 26 of the 49 structures and below both on 4,
and the difference was present in every strain (median ratio 1.07 for BSG to 1.37 for BSA). The ratio of the
oscillation chambers to their controls did not depend on the period: it correlated with the period at |ρ| ≥ 0.6
in two to five of the ten series per readout, with both signs. The immigration rate, which reflects the flow
through a chamber, did not differ between oscillation and control chambers in either oscillation type (ratio
1.02 and 1.10, p 0.90 and 0.45). The oscillation chambers occupied the positions A3 to A12 of an array and the
controls the positions A1, A2, A13 and A14 (Section 2.5.3), so the comparison is also one between chamber
positions (Appendix Figure A13).

**Trends against the period.** Monotone trends of the oscillation chambers with |ρ| ≥ 0.6 appeared in every
readout, but in different series (Table 5, columns "oscillation chambers"). The budding rate fell with the period
in WT/pH (ρ −1.00), WT/Glc, BSA/pH and BSG/pH (−0.60) and rose in BSO/Glc (+0.77); µ_bud rose in BSO/Glc (+1.00)
and fell in BSA/pH, BSG/Glc and BSG/pH (−0.60). The endpoint area fell with the period in BSpH/Glc (−0.89) and
WT/pH (−0.80) and rose in BSA/Glc (+0.71); µ_area fell in BSA/Glc (−0.94), BSO/Glc (−0.83) and WT/pH (−0.80). The
eccentricity trended in seven of the ten series (|ρ| 0.60 to 0.83), positively in six of them. No trend recurred
across strains within an oscillation type, and none recurred in both oscillation types of one strain.

**Feast against famine controls.** The constant-feast and constant-famine chambers of a structure barely
separated. For the budding rate the bracket was degenerate on 44 of 49 structures, and the feast control exceeded
the famine control on 26 of 48 structures (Wilcoxon signed-rank test, p 0.61); for µ_bud on 43 of 49 and 23 of 48
(p 0.91). Cells under constant famine budded in the sparse phase as often as cells under constant feast (0.05 to
0.63 against 0.00 to 0.94 buds per mother-hour). The endpoint area and eccentricity had a direction, feast cells
being larger on 31 and more eccentric on 33 of 48 structures (p 0.01 each), but the bracket was degenerate on 43
and 39 of 49 structures.

**Within cultures.** Inside the 19 cultures that carried two or three periods, where the preculture could not
differ, the endpoint area fell from the shortest to the longest period in 8 cultures for the oscillation chambers
and in 13 for their controls, and the difference between the two fell in 6; for µ_bud the counts were 8, 8 and
11 of 19. Neither readout changed with the period in a consistent direction within cultures.

**Temporal and population robustness.** The robustness metrics R(t) and R(p) (Section 2.5.5) were computed for
the cell area, the eccentricity, the area growth rate µ_area (R(p) across the cells of a chamber) and the budding
of the mothers of a chamber (R(p) of the budding rate per mother across the mothers, and single-cell R(t) of the
event-based rate µ_event over the intervals of mothers with at least three buddings), and the oscillation
chambers were compared with the controls of their structure in the same way as above (Table 4; Appendix
Figure A11). Under oscillation the mean cell area of a chamber was less stable over time than under
constant medium (R(t) at the population level −0.127 against −0.096 for the feast and −0.123 for the famine
controls; below the control mean on 31 of 49 structures, p 0.013), and so was the area of the individual cell
(single-cell R(t) −0.106 against −0.100 and −0.077, below the control mean on 35 of 49, p 0.007). The
population was also more heterogeneous in cell area (R(p) −0.671 against −0.571 and −0.613, on 34 of 49
structures, p 0.003) and in eccentricity (−0.048 against −0.041 and −0.042, p 0.003), but more homogeneous in
its area growth rate (R(p) of µ_area −0.333 against −0.412 and −0.453, above the control mean on 33 of 49
structures, p 0.008). The heterogeneity of cell area and that of the growth rate were unrelated across chambers
(Spearman ρ 0.06, n = 497). The reproductive output of the mothers was as heterogeneous under oscillation as
under constant medium (R(p) of the budding rate per mother −2.02 against −1.98 and −1.92, p 0.95), whereas the
budding rhythm of the individual mother was less stable under glucose oscillation (single-cell R(t) of µ_event
−1.11 against −1.13 and −1.00; below the control mean on 29 of 47 structures, p 0.022; Glc p 0.043, pH
p 0.28). Against the period the robustness metrics behaved like the readouts themselves: of the 90 combinations
of the nine metrics and ten series, 51 showed no monotone trend of the oscillation chambers, 23 a trend shared
by a control of the same structures, 8 a trend that vanished after subtraction of the controls, and 8 met the
conditions of a period effect (population R(t) of the area in WT/Glc, single-cell R(t) of the area in BSG/pH,
R(p) of the area in BSA/Glc and BSpH/pH, and R(p) of the eccentricity, R(p) of µ_area, R(p) of the budding rate
per mother and single-cell R(t) of µ_event in BSA/Glc). Five of the eight fell in the series BSA/Glc, in which
the population became more homogeneous in area, in µ_area and in budding and the budding rhythm of the mothers
more stable, but the population more heterogeneous in eccentricity, with longer periods.

**Table 4: Temporal and population robustness of growth and morphology under oscillation and under constant
medium.** R(t) and R(p) after (2022_Trivellin; 2024_Blöbaum), R ≤ 0 with 0 for a perfectly stable readout; medians
over structures of the chamber means, single-cell R(t) over cells tracked for at least ten frames and, for
µ_event, over mothers with at least three budding intervals; R(p) of the budding rate per mother across the
mothers of a chamber. "Less robust": structures on which R of the oscillation chambers lay below the mean of the
two controls; p: Wilcoxon signed-rank test over the structures; period effects: classification of 3.6 over the
ten series. [`40_<metric>_<readout>_per_chip.csv`, `51_osc_vs_controls.csv`,
`40_<metric>_<readout>_control_trend.csv`]

| metric | readout | oscillation | feast control | famine control | less robust (of 49) | p | period effects (of 10) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| R(t), population | cell area | −0.127 | −0.096 | −0.123 | 31 | 0.013 | 1 |
| R(t), population | eccentricity | −0.015 | −0.012 | −0.014 | 33 | 0.019 | 0 |
| R(t), single cell | cell area | −0.106 | −0.100 | −0.077 | 35 | 0.007 | 1 |
| R(t), single cell | eccentricity | −0.022 | −0.020 | −0.018 | 30 | 0.010 | 0 |
| R(p) | cell area | −0.671 | −0.571 | −0.613 | 34 | 0.003 | 2 |
| R(p) | eccentricity | −0.048 | −0.041 | −0.042 | 33 | 0.003 | 1 |
| R(p) | µ_area | −0.333 | −0.412 | −0.453 | 16 | 0.008 | 1 |
| R(p) | budding rate per mother | −2.018 | −1.984 | −1.921 | 25 | 0.95 | 1 |
| R(t), single cell | µ_event (interbud rate) | −1.115 | −1.132 | −0.999 | 29 of 47 | 0.022 | 1 |

## 3.5 Biosensor readouts against the half-cycle period

The four biosensor strains reported their ratiometric signal (sensor over reference channel) at the endpoint of
every chamber (Figure 14, Appendix Figure A8). The ratios lay at 0.10 to 0.16 for the ATP sensor (BSA,
QUEEN-2m), 0.01 to 0.04 for the glucose-flux sensor (BSG, Gly-RNA; one feast-control structure at 1.00), 0.00 to
0.02 for the oxidative-stress sensor (BSO, OxPro) and 0.63 to 0.96 for the pH sensor (BSpH, sfpHluorin).
Figure 13 shows one cell of every biosensor strain in the merged image, the phase contrast and the two fluorescence
channels that form the ratio.

**Figure 13: Example cells of the four biosensor strains.** From left to right: merged image, phase contrast,
first and second fluorescence channel of the ratio; BSA (QUEEN-2m) and BSpH (sfpHluorin): uvGFP and GFP channel
*[excitation 405 and 488 nm? check]*, BSO (OxPro): YFP and mCherry, BSG (Gly-RNA): mTurquoise and mCherry.
*[chamber and time point; scale bar]* [your microscopy images]

None of the four sensors showed a consistent difference between its constant-feast and constant-famine chambers:
the bracket was degenerate on 6 to 9 of the 10 structures per sensor, and the pairwise comparison of the two
control types within every structure (Appendix Figure A8) gave no direction for any sensor. The dynamic range
of the sensors under the cultivation conditions in the chip is therefore not established by these data.

Against the period, the OxPro ratio of BSO rose in the oscillation chambers of both oscillation types (ρ +0.77
in Glc, +1.00 in pH), and so did the pHluorin ratio of BSpH (+0.71 and +0.60) and the Gly-RNA ratio of BSG in
the Glc series (+0.89); the QUEEN-2m ratio of BSA showed no trend. In every one of these cases the famine-control
chambers of the same structures rose with the period as well (OxPro +0.94 and +1.00, pHluorin +0.89 and +1.00,
Table 5), with the exception of Gly-RNA, whose trend did not survive the subtraction of its controls (3.6).

**Figure 14: Biosensor ratios of the Glc oscillation series against the half-cycle period.** (A) OxPro ratio of BSO,
(B) sfpHluorin ratio of BSpH; layout as in Figure 11. [`13_endpoint_vs_period_ratio_OxPro_Glc.pdf`,
`13_endpoint_vs_period_ratio_pHluorin_Glc.pdf`]

*[Optional, not in your chapter:]* **Robustness of the sensor signals.** The temporal and population robustness of the four ratios differed by
two orders of magnitude between the sensors (Table S1): the sfpHluorin ratio was the most stable over time and
the most homogeneous across cells (R(t) and R(p) of −0.01), the OxPro ratio, whose values lie close to zero, the
least (R(t) −0.21 and R(p) −0.73 in the oscillation chambers). For none of the four sensors did the oscillation
chambers differ from the controls of their structure in R(t) or R(p) (Wilcoxon signed-rank test over ten
structures, p 0.06 to 0.92). Against the period, 12 of the 24 combinations of sensor, metric and series were
structure effects, among them both series of OxPro in all three metrics, in which the robustness of the ratio
fell with the period in the oscillation and in the famine chambers alike; one combination met the conditions of
a period effect (single-cell R(t) of the Gly-RNA ratio in BSG/pH), two were not robust and nine showed no trend.

**Table S1: Temporal and population robustness of the four sensor ratios.** *[Optional, not in your chapter.]*
Layout as in Table 4; ten structures per sensor. [`40_<metric>_ratio_<sensor>_per_chip.csv`, `51_osc_vs_controls.csv`]

| sensor | metric | oscillation | feast control | famine control | less robust (of 10) | p |
| --- | --- | --- | --- | --- | --- | --- |
| QUEEN-2m (BSA) | R(t), population | −0.140 | −0.113 | −0.137 | 4 | 0.77 |
| | R(t), single cell | −0.037 | −0.034 | −0.040 | 5 | 0.92 |
| | R(p) | −0.049 | −0.074 | −0.058 | 2 | 0.28 |
| Gly-RNA (BSG) | R(t), population | −0.015 | −0.009 | −0.022 | 3 | 0.28 |
| | R(t), single cell | −0.042 | −0.044 | −0.023 | 5 | 0.92 |
| | R(p) | −0.065 | −0.128 | −0.117 | 5 | 0.77 |
| OxPro (BSO) | R(t), population | −0.212 | −0.173 | −0.096 | 5 | 0.92 |
| | R(t), single cell | −0.249 | −0.292 | −0.127 | 4 | 0.56 |
| | R(p) | −0.728 | −1.032 | −0.380 | 4 | 0.62 |
| sfpHluorin (BSpH) | R(t), population | −0.010 | −0.010 | −0.010 | 7 | 0.32 |
| | R(t), single cell | −0.012 | −0.010 | −0.010 | 8 | 0.06 |
| | R(p) | −0.010 | −0.005 | −0.015 | 5 | 0.43 |

## 3.6 The constant-medium controls on the same structures

The control chambers of a structure received constant medium and could not respond to the period. Their values
were nevertheless compared across the structures of each series in the same way as the oscillation chambers
(Section 2.5.6). For every readout and series, the Spearman correlation of the oscillation chambers, of the feast
and famine controls, and of the difference between the oscillation chambers and the control mean was computed
against the period, and each of the 58 readout-by-series combinations was classified: a period effect required
that the oscillation chambers trended (|ρ| ≥ 0.6), that no control type trended in the same direction, and that
the difference trended as well.

**Table 5: Trends of the oscillation chambers and of their controls against the half-cycle period.** Spearman ρ
over structures for all readout-by-series combinations in which the oscillation chambers trended with |ρ| ≥ 0.6;
"strongest control" is the control type (feast, famine or their mean) with the largest |ρ|; "difference" is
oscillation chambers minus control mean. Residual: the difference trends in the same direction with |ρ| ≥ 0.6
although a control trends too. [`50_control_trend_summary.csv`, `24_growth_from_budding_control_trend.csv`]

| readout | series | oscillation chambers | strongest control | difference | classification |
| --- | --- | --- | --- | --- | --- |
| buds per mother-hour | WT/pH | −1.00 | −1.00 (feast) | +0.40 | structure effect |
| | BSO/Glc | +0.77 | +0.66 (famine) | +0.77 | structure effect, residual |
| | BSA/pH | −0.60 | −1.00 (feast) | +0.20 | structure effect |
| | BSG/pH | −0.60 | −0.80 (feast) | −0.20 | structure effect |
| | WT/Glc | −0.60 | 0.00 | −0.80 | period effect |
| µ_bud | BSO/Glc | +1.00 | +0.60 (famine) | +0.09 | structure effect |
| | BSA/pH | −0.60 | −1.00 (feast) | +0.20 | structure effect |
| | BSG/Glc | −0.60 | −0.77 (famine) | +0.14 | structure effect |
| | BSG/pH | −0.60 | −0.80 (feast) | −0.20 | structure effect |
| OxPro ratio | BSO/Glc | +0.77 | +0.94 (famine) | +0.26 | structure effect |
| | BSO/pH | +1.00 | +1.00 (famine) | +0.40 | structure effect |
| sfpHluorin ratio | BSpH/Glc | +0.71 | +0.89 (famine) | −0.09 | structure effect |
| | BSpH/pH | +0.60 | +1.00 (famine) | −0.40 | structure effect |
| Gly-RNA ratio | BSG/Glc | +0.89 | +0.43 | −0.26 | not robust |
| endpoint area | BSA/Glc | +0.71 | +0.94 (mean) | +0.60 | structure effect, residual |
| | BSpH/Glc | −0.89 | −0.89 (famine) | −0.26 | structure effect |
| | WT/pH | −0.80 | −0.80 (famine) | +0.40 | structure effect |
| µ_area | BSA/Glc | −0.94 | −0.60 (famine) | −0.60 | structure effect, residual |
| | BSO/Glc | −0.83 | −0.49 | −0.77 | period effect |
| | WT/pH | −0.80 | −1.00 (famine) | +0.60 | structure effect |
| eccentricity | BSO/Glc | +0.83 | +0.77 | +0.14 | structure effect |
| | WT/pH | +0.80 | +0.80 | −0.60 | structure effect |
| | BSpH/Glc | +0.71 | +0.77 | +0.49 | structure effect |
| | BSG/Glc | +0.66 | +0.77 | +0.60 | structure effect, residual |
| | BSA/pH | +0.60 | +1.00 | 0.00 | structure effect |
| | BSG/pH | −0.60 | −0.60 | −0.40 | structure effect |
| | BSpH/pH | +0.60 | +0.50 | −0.20 | not robust |

**Table 6: Classification of the 58 readout-by-series combinations** (endpoint area, eccentricity, four sensor
ratios, µ_area, budding rate and µ_bud over ten series). [`50_control_trend_summary.csv`]

| classification | combinations |
| --- | --- |
| no monotone trend of the oscillation chambers | 31 |
| structure effect: a control trends with the period in the same direction | 23 |
| not robust: the trend vanishes after subtracting the controls | 2 |
| period effect: oscillation chambers trend, controls do not, difference trends | 2 |

In 23 of the 27 combinations in which the oscillation chambers trended with the period, at least one control type
of the same structures trended in the same direction, in 15 of them with a |ρ| of 0.8 or more (Table 5). Two
trends vanished after the controls were subtracted (eccentricity BSpH/pH, Gly-RNA ratio BSG/Glc). Two combinations
met all three conditions of a period effect: µ_area in BSO/Glc (ρ −0.83, difference −0.77) and the budding rate in
WT/Glc (ρ −0.60, difference −0.80). Neither recurred in the other oscillation type of the same strain (µ_area
BSO/pH ρ −0.20; budding rate WT/pH a structure effect). Figure 15 shows all 58 combinations at once: the correlation of the oscillation chambers on the x axis and that of the strongest control on the y axis, where a
point on the diagonal is a structure effect; 23 points lie in the two diagonal corners, 31 in the central band,
and the two period effects and the two non-robust trends lie off the diagonal.

**Figure 15: Trends of the oscillation chambers against trends of their controls.** One point per readout and
series (58 combinations): Spearman ρ of the oscillation chambers against the period on the x axis, ρ of the
strongest control type on the y axis, marker = readout, colour = strain, hollow = not robust; grey zones mark
|ρ| < 0.6. [`50_control_trend_summary.pdf`; `_growth.pdf` and `_sensors.pdf` are the same figure split for the
appendix]

The control chambers also differed between the structures of a series in their cell-to-cell heterogeneity: the
population robustness R(p) of the cell area differed between structures for every control type and series
(Kruskal-Wallis over chambers, p < 0.05 in 20 of 20 tests), R(p) of the eccentricity in 19 of 20 and R(p) of the
sensor ratios in 16 of 16, whereas R(p) of µ_area differed in 2 of 20 and the temporal robustness R(t) of area
and eccentricity in 5 and 4 of 20 (Appendix Figure A9).

The periods of a series had been run in day blocks (Section 2.4.3): the Glc periods as three plus three on two
days, in ascending order for WT and BSA (Spearman of period against date +0.87 and +0.88), in descending order for
BSG and BSO (−0.88), and all six on one day for BSpH; the pH periods as two plus two, in ascending order for
every strain (+0.89). The shortest period of a series usually occupied the same structure position on the chip
*[author statement; no run sheet of the positions exists]*. Period, structure position and cultivation day were
therefore aligned within every series.

## 3.7 The pullulan knockout strain

The pullulan knockout strain was cultivated on one structure of the Glc series at a period of 3 min, to test
whether the control chambers of a structure agree better in the absence of the exopolysaccharide. The
chamber-to-chamber variation of the cell-area level among the three famine-control chambers of the knockout
structure was 0.23 (coefficient of variation) and among the three feast-control chambers 0.37, which placed the
knockout at the 67th and 96th percentile of the 49 producer structures (Figure 16). The temporal variation of the
area after removing the linear trend lay at the 92nd and 61st percentile. Famine-control cells of the knockout
were 1.5 times larger than its feast-control cells, where the producer structures had a median ratio of 0.98
(knockout at the 92nd percentile); the µ_area bracket of feast over famine was 0.10, the producers' median. The
famine controls of the knockout also budded faster than its feast controls (0.32 to 0.43 against 0.18 to 0.22
births per cell-hour). The single knockout structure therefore did not show a better agreement of its control
chambers than the pullulan-producing strains.

**Figure 16: Agreement of the control chambers on the pullulan-knockout structure against the producer
structures.** Chamber-to-chamber coefficient of variation of the cell-area level and detrended temporal
coefficient of variation for the feast and famine controls of every structure; the knockout structure is marked.
[`pko/61_pko_control_agreement.pdf`]

---

## Notes to the author (not thesis text)

**Figure and table sources.** Every pipeline file is relative to `analysis_output_v12/` of the final run.

| thesis figure | pipeline file(s) | status |
| --- | --- | --- |
| Figure 2 | `growth_curves.pdf` (BioLector script) | exists; a.u., OD₆₀₀ calibration pending |
| Figure 3 | `00_n_tracks_overview_summary.pdf` | exists; x label "half-cycle period" with the next run |
| Figure 4 | `90_morphology_scatter.pdf`, `static/90_morphology_scatter.pdf` | with the next pipeline run (per-cell version) |
| Figure 5 | QC overlays from `03_results_v12/QC/` (your TIFFs) | exists |
| Figure 6 | `21_panel_a_violin.pdf` (area and eccentricity rows) | with the next run (one value per cell, test over structures, µm²) |
| Figure 7 A | `20_lineage_window.pdf` | exists |
| Figure 7 B | `00_new_objects_vs_density.pdf` | with the next pipeline run |
| Figure 8 | `20_bud_size_at_appearance.pdf` | exists |
| Figure 9 | `lineage_validation/lv_02_detection_rate.pdf` | with the next validation run (no static panel) |
| Figure 10 | `static/13_endpoint_vs_medium_area.pdf`, `static/21_budding_rate_vs_medium.pdf`, `static/24_growth_from_budding_vs_medium.pdf` | with the next run (per-chamber points, µm²) |
| Figure 11 | `21_budding_rate_vs_period_Glc.pdf`, `24_growth_from_budding_vs_period_Glc.pdf`, `13_endpoint_vs_period_area_Glc.pdf` | with the next run (µm², half-cycle label) |
| Figure 12 | `24_mu_bud_vs_mu_area.pdf` | exists |
| Figure 13 | your microscopy images | exists |
| Figure 14 | `13_endpoint_vs_period_ratio_OxPro_Glc.pdf`, `13_endpoint_vs_period_ratio_pHluorin_Glc.pdf` | with the next run (label) |
| Figure 15 | `50_control_trend_summary.pdf` | with the next run (combined figure restored 2026-10-05) |
| Figure 16 | `pko/61_pko_control_agreement.pdf` | exists |
| Figure 17 (appendix) | `95_<strain>_<osc_type>_ratio_<sensor>_comparison.pdf` | exists |
| Appendix A4 | `static/10_cell_area_over_time.pdf`, `static/13_endpoint_vs_medium_eccentricity.pdf` | with the next run (µm²) |
| Appendix A5 | `13_endpoint_vs_period_area_pH.pdf`, `24_growth_from_budding_vs_period_pH.pdf`, `21_budding_rate_vs_period_Glc.pdf` and `_pH.pdf` | with the next run |
| Appendix A6 | `24_immigration_vs_period_Glc.pdf`, `24_immigration_vs_period_pH.pdf` | exist |
| Appendix A7 | `13_endpoint_vs_period_eccentricity_Glc.pdf`, `_pH.pdf`, `12_area_growth_rate_all.pdf` | exist |
| Appendix A8 | `13_endpoint_vs_period_ratio_<sensor>_<Glc,pH>.pdf` (the remaining sensor endpoint figures) | exist |
| Appendix A9 | `40_control_consistency_*.pdf` | exist |
| Appendix A10 | `pko/13_endpoint_vs_period_area_Glc.pdf` | exists |
| Appendix A11 | `40_<metric>_<readout>_vs_period_<osc_type>.pdf` (R(t) population, R(t) single cell, R(p) of area, eccentricity, µ_area, budding rate per mother, µ_event) | with the next pipeline run; layout of Figure 11 with the controls |
| Appendix A12 | `50_robustness_control_trend_summary.pdf` | with the next pipeline run |
| Appendix A13 | `51_osc_vs_controls.pdf`, `51_osc_vs_controls_robustness.pdf`, `51_osc_vs_controls_sensors.pdf` | with the next pipeline run |

**Open items.**

- OD₆₀₀ calibration of the BioLector data (3.1): send the end-point OD₆₀₀ readings and I add the calibration to
  `biolector_plot.py`.
- The W65 overlay of frames 85 to 88 (3.3) could become an appendix figure from your TIFF.
- The methods chapter must describe the final pipeline before this chapter can cite it; the replacement text for
  Sections 2.3 (growth parameters) and 2.5 is in `docs/methods_2_5_rewrite.md`.
- Conversion used throughout: 1 px² = 0.00537 µm² (0.0733 µm per pixel).
- Pipeline outputs for the robustness results (built 2026-10-01, in the next run): step 40 writes every robustness
  metric per chamber through the chip logic of the readouts (`40_<metric>_<readout>_per_chip.csv`, `_spearman`,
  `_bracket_score`, `_control_trend`, `_vs_period_<osc_type>.pdf` with the controls), step 50 collects them into
  `50_robustness_control_trend_summary.*` and writes the paired comparison of the oscillation chambers with the
  controls of their structure for all readouts and metrics (`51_osc_vs_controls_per_structure.csv`,
  `51_osc_vs_controls.csv`, `51_osc_vs_controls_vs_period.csv`, three figures). Until that run, the numbers of
  3.4, 3.5 and Tables 8 and 9 come from `docs/scratch/robust3.py`, the same functions on the final-run tables.
