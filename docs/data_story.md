# The story of the data

Argumentation of the thesis results, with the figures and tables of the pipeline that carry each step.
Draft of 2026-09-22, rewritten 2026-09-27 on the complete v12 run: segmentation and tracking rebuilt
(`docs/tracking_plan.md`), output `analysis_output_v12/` (paths below are relative to it). Numbers are from the
final run of 2026-09-28 (564 chambers: structure-relative contrast rule, W65 minimal Rep5 dropped as a
duplicate, W109 one chip, growth rate from budding).
Numbers marked *WT/pH/6* come from the v12 `Combined_Results.csv` of that batch, the batch that carries the
manual QC of the v11 run; *BSG/pH/6* is the debris-rich batch used for the dead-cell rule. Numbers marked *v11*
come from the earlier tables and are kept only where the v12 run has no counterpart yet. Statements marked
**(open)** need the author's input, see section 9.

---

## 0. In one paragraph

An oscillation experiment across five strains, two media and six feast/famine periods produced apparent
dose responses in the budding rate, in cell shape and in three biosensor ratios. Every one of them turned
out to be a property of the microfluidic structure, not of the period: the constant-medium controls on the
same structures drift identically, and after subtracting them nothing is left that recurs across strains.
Rebuilding the segmentation and the tracking halved the fragmentation of the cell tracks, gave every bud a
measured mother and removed the budding trends the old tracker had produced; what survived the rebuild are
the structure effects. Two properties of the design made this inevitable, and one decision made it visible.
Each period sat on its own structure, usually in the same position on the chip, so period and position were
aligned; the periods were shorter than the sampling interval, so only cumulative effects could ever be read;
and every structure carried its own feast and famine controls, which is the only reason the drift could be
told from a treatment effect. In a crowded chamber every new object touches a cell, so budding is still
counted only in the sparse first hours of a chamber. The static comparison of complex versus minimal medium
on the W109 chips shows a coherent medium effect, cells three times larger that bud a sixth as often in
complex medium; whether its four replicates are four cultures or one is still open, and W65, one chip and
one culture, does not repeat it.

---

## 1. What was measured, and what the units are

### 1.1 Design

| element | value |
| --- | --- |
| organism, imaging | *Aureobasidium pullulans*, microfluidic chambers, one frame per 10 min, 133 frames in the standard run (about 22 h; 13 to 136 in the tables, one chamber broke off after 13 frames); 0.0733 µm per px, a blastoconidium of 20 µm² is 3,700 px² |
| segmentation, tracking | Cellpose `cpsam`, flow 0.4, cell probability 0 (`imaging/pipeline_template_v12.yaml`); tracking on the raw labels before any filtering, the mother of a new object from the mask it touches (`imaging/cellpose_pipeline_v12.py`, `imaging/track_labels.py`) |
| strains | WT and four biosensor strains: BSA (`ratio_Queen-2m`), BSG (`ratio_GlyRNA`), BSO (`ratio_OxPro`), BSPH (`ratio_pHluorin`); PKO = no pullulan |
| oscillation types | Glc and pH switching; the "period" of a structure is the switching interval of the medium, i.e. the half-cycle period (`osc_freq`): 0.75, 1.5, 3, 6, 12, 24 min (Glc; WT lacks 24) and 0.75, 1.5, 6, 24 min (pH); a full feast/famine cycle lasts twice as long (1.5 to 48 min) |
| one structure | one period; 5 oscillation chambers, 3 constant-feast (PosCtrl) and 3 (sometimes 4) constant-famine (NegCtrl) chambers on several arrays; `replicate` is the array index, `chamber` the position (A1/A2 feast, A13/A14 famine, A3 to A12 switching) |
| one physical chip | 2 or 3 structures = 2 or 3 periods, one pre-culture, one day (`culture`) |
| oscillation and PKO total | 50 structures, 547 chambers (8 to 12 per structure), 20 cultures |
| static | 17 chambers on two chips. W109: one chip with three structures, four chambers per medium; W65: one chip, one pre-culture, four chambers in minimal and five in complex medium; `chamber` always A0 |

The pipeline calls a structure `chip` in every table and figure. `culture` is the physical chip.

- `00_chip_overview.csv`: one row per structure with its chambers per control type, its culture and `n_structures_in_culture`.
- `00_chip_run_order.csv`: for each strain and oscillation type, the periods in run-date order and Spearman(period, date). The periods were run in day blocks: Glc as 3 + 3 periods on two days (BSPH all six on one day), pH as 2 + 2.
- `00_data_overview.csv`, `00_n_tracks_overview*.pdf`: chambers, frames and objects per frame.

### 1.2 What the design can and cannot separate

- **Within a period there is no biological replication.** Five oscillation chambers of a period are five
  chambers of one structure. Every error bar in the oscillation figures is a chamber error bar; the tables
  say so in `error_unit` and count `n_units` accordingly.
- **Period and structure are aligned.** One structure per period, and the shortest period usually on the same
  structure position **(author statement)**. Whatever varies with position (flow path, loading order, focus)
  varies with the period.
- **Period and culture are aligned in blocks.** Short periods on one day, long periods on another. For Glc the
  block order differs between strains (BSA and WT short first, BSG and BSO long first, BSPH all six on one
  day); for pH every series ran the short periods first.
- **The controls sit on the structure.** PosCtrl and NegCtrl chambers of every structure receive constant
  medium and cannot respond to the period. They are the only handle on the structure and culture effects.

### 1.3 The sampling limit

Ten-minute frames resolve full periods of 20 min and longer. A full feast/famine cycle is twice the switching
interval, so the intervals 0.75 to 6 min (full periods 1.5 to 12 min), four of the six Glc intervals and three of
the four pH intervals, lie below that; the 12-min interval is at the limit (full period 24 min, 2.4 frames per
cycle) and the 24-min interval is resolved (48 min, 4.8 frames per cycle). For the short intervals single cycles
are not observable, and apparent periodicity in their time series is an alias. The pipeline therefore reads
cumulative quantities for all intervals alike: the endpoint of a chamber (`13_`), the budding rate over the sparse
window (`21_`), µ_bud (`24_`) and the robustness metrics (`40_`). The time series (`10_`, `30_`, `31_`) are shown
as drift over hours, never as cycles. (`config.py` logs the unresolved and marginal intervals at every run.)

---

## 2. The tracking limit and the sparse-phase window

### 2.1 The tracker fragments, half as much as before

The v11 tracker matched masks frame to frame by an overlap threshold, had no memory and reset on empty
frames; every object received a new track ID with a probability of about 12 % per frame
(`docs/tracking_diagnosis.md`). The rebuilt tracker (assignment by distance, area and overlap, a memory of
three frames, merges and splits by overlap) halves that.

| `00_track_fragmentation.csv`, per chamber, after the cell filter | q10 | median | q90 |
| --- | --- | --- | --- |
| new tracks per object-frame | 0.047 | 0.070 | 0.110 |
| median track length, frames | 4 | 5 | 12 |
| share of object-frames in tracks of 10 frames or more | 0.60 | 0.78 | 0.89 |

*v11*, same table before gap closing: 0.09 / 0.14 / 0.22 new tracks per object-frame, median track 1 / 3 / 4
frames, 31 % of the tracks one frame long. *WT/pH/6*: 4,894 tracks in 11 chambers in the raw v12 table, 1,173
of them one frame long (v11: 9,063 and 2,712). Glc structures fragment more than pH structures (median 0.080
against 0.057). The fragmentation rises with the period in two of the ten series (BSA/Glc ρ 0.89, BSO/Glc
0.71) and falls in one (BSG/pH −0.60).

### 2.2 New objects follow the density

Chambers start with 1 to 3 objects and fill to as many as 180. In the dense phase most new objects touch an
established cell and are bud candidates. *WT/pH/6*, v12 table, means over the 11 chambers by 22-frame block
(tracks present in the first frame are not counted as new):

| block (frames) | 0–21 | 22–43 | 44–65 | 66–87 | 88–109 | 110–131 |
| --- | --- | --- | --- | --- | --- | --- |
| objects per frame, median | 2.0 | 5.2 | 12.8 | 32.9 | 83.6 | 131.7 |
| new tracks touching or split from a tracked mask | 0.5 | 5.3 | 8.2 | 33.2 | 83.3 | 129.8 |
| new tracks touching nothing | 2.5 | 5.9 | 13.5 | 33.7 | 50.1 | 42.5 |

237 touching new tracks per chamber over the run, 130 of them in the last block, where their number equals
the object count: at 130 objects per frame the masks touch, split and re-merge, and a "bud" is any fragment.
The measured parent does not remove this; only the sparse window does (2.4). Evidence in the pipeline:
`20_lineage_window.pdf` (objects per frame over time for every chamber) and `00_track_fragmentation.csv`
(`objects_first_frame`, `objects_last_frame`).

`00_new_objects_vs_density.pdf` draws this table for every chamber (`00_new_objects_vs_density.csv`,
`relink.new_objects_vs_density()`, block length `DENSITY_BLOCK_FRAMES` = 22): one point per chamber and
block, new tracks per frame against objects per frame, separately for new tracks that touch a tracked mask or
split from one and for new tracks that touch nothing, with median and quartiles per density class; the right
panel gives the touching new tracks per object and frame. On the four v12 tables at hand (35 chambers of
WT/pH/6, BSG/pH/6, W65 and W109) the touching new tracks rise in proportion to the density, 0.04 to 0.06 per
object and frame between 8 and 256 objects per frame, while the new tracks touching nothing saturate at about
two per frame; above the sparse limit 55 to 63 % of the new tracks touch a mask on the three growing chips and
35 % on W109, which never gets dense. The figure over all 564 chambers comes with the next run.

### 2.3 What replaced the manual quality control

The manual QC of WT/pH/6 (503 rows: 248 merges after tracking breaks, 214 background objects, 16 mother
cells on the edge, 13 dead cells, 11 debris) names v11 track IDs and is not applied to the v12 tables. It
served twice: it showed that the v11 tracker, not the cells, produced the breaks (31 % of the new IDs were
losses of a tracked object; `docs/tracking_diagnosis.md`), and it calibrated the new tracker (34 % of the
manual merges fall inside its memory; `docs/tracking_plan.md` 6c). Its four kinds of correction are now rules:

- tracking breaks: the tracker's memory and its merge/split handling (2.1);
- background objects and debris: the cell filter (`cell_filter.py`, `00_cell_filter.csv`): a track counts as
  a cell with at least 2 frames and a largest area of at least 1,500 px² (8 µm²). Full run: 35 % of the tracks
  and 6.8 % of the object-frames removed (median chamber 5.7 %), of which the contrast rule below accounts for
  0.6 % of all object-frames; BSG/pH/6 loses 52 % of its object-frames, and the objects there are dead-cell
  debris by eye;
- dead cells: the contrast rule of the same filter (2.7);
- cells on the edge: the `at_border` flag of pipeline v12, tracked but excluded from every measurement
  (13 % of the rows of WT/pH/6).

The v11 comparison with and without the manual file (`qc_comparison/70_qc_effect.pdf`, *v11*) had shown that
it changed 2.6 % of the raw events over the run and 19 % inside the sparse window; repeating it for 50 batches
would not have changed the lineage readouts.

### 2.4 Sparse-phase window and events

- **Window** (`20_lineage_window.csv`, `20_lineage_window.pdf`): per chamber, the frames before the rolling
  median of objects per frame exceeds 20; chambers with fewer than 20 such frames drop out. Full run: 563 of
  564 chambers qualify (the one that broke off after 13 frames does not), median window 117 of 133 frames,
  range 26 to 136; 245 chambers never exceed 20 objects and are used whole. Every lineage output (`20_` to
  `23_`, `21_budding_rate_*`, µ_event, mother/bud labels) uses only this window; endpoint, µ_area and
  robustness still see the whole run.
- **Gap closing** no longer runs in the analysis: the tracker's memory (3 frames, 15 for a track absorbed by
  a neighbour's mask) does the same work on the masks, and every row records how its object was linked
  (`link_type`: continued, gap, long_range, unmerge, split, new_touching, new).
- **Mother from the mask** (`20_budding_events.csv`, column `method`): a new object that touches a tracked
  mask at first detection has that mask's track as `parent_track_id`; the lineage takes the parent as the
  mother when it is an established cell and the bud persists for 2 frames or more. Candidates without a
  touching mask go through the old radius heuristic. Full run: 11,361 events in 531 chambers, median 18 per
  chamber (q10 8, q90 39); 74 % measured, 26 % heuristic; 72 % of the mothers were tracked for 30 frames or
  more. *WT/pH/6*: 93 events in 10 chambers, 4 to 14 per chamber, 57 measured and 36 heuristic (v11 with gap
  closing: 136).
- **Validation** (`lineage_validation/`, `validate_lineage.py` on the v12 tables): the share of bud candidates
  in the window that are assigned to a mother differs between the structures of a series in 4 of 10 series
  (Kruskal p < 0.05; v11: 7 of 10) and trends with the period in no consistent direction (BSA/Glc −0.77,
  BSPH/Glc −0.66, BSG/pH −0.80, BSO/Glc +0.60; `lv_03_detection_rate_kruskal.csv`). Median
  assignment rate 0.67 per chamber, median distance over search radius 0.93, 4 % of the candidates ambiguous
  (`lv_02_detection_rate_per_chamber.csv`). Inflowing cells are new tracks without a mother; a washed-in cell
  that lands on a mother still counts, which is what the size criterion (2.5) removes.

### 2.5 The size criterion, now active

A bud should be small at first detection and a washed-in cell mother-sized. In the v11 tables the ratio of
bud area to mother area had one broad mode around 0.4 to 0.5 and the criterion was inactive. With the mother
from the mask contact the distribution is bimodal: modes at 0.07 and 0.45, antimode 0.32, valley depth 0.25
(`20_bud_size_threshold.csv`, `20_bud_size_at_appearance.pdf`). The pipeline applies the global threshold of
0.32 and rejects 2,162 of 13,940 candidates (15.5 %) as mother-sized objects that appear next to a mother:
washed-in cells and masks that split in two (*WT/pH/6*: 30 % of the touching new objects and 29 % of the
splits lie above 0.32). Per strain and oscillation type the valley is too shallow in 10 of 12 groups, so the
one global threshold stays. The mother's mask does not change when a bud is first segmented **(author
statement)**, so the area balance of the mother is not a criterion; the diagnostic columns stay in
`20_bud_size_at_appearance.csv`.

### 2.6 What remains trustworthy

Readouts that need no lineage: endpoint area and eccentricity per chamber (`13_`), population R(t) and R(p)
(`40_Rt_population_*`, `40_Rp_*`), sensor ratios (`30_`, `31_`, `95_`). µ_area (`12_`) fits ln(area) per track
and requires 10 frames. The lineage readouts describe the sparse phase only, and their per-series trends
against the period changed with the tracking (section 4).

### 2.7 Dead cells lose their phase contrast

Pipeline v12 measures the mean and the standard deviation of the phase-contrast signal inside every mask
(`phase_mean`, `phase_std`). Their ratio per track separates dead from live objects: on *BSG/pH/6* it is
bimodal (modes 0.05 and 0.28, valley 0.11 to 0.16, 71 % of the tracks in the low mode), on *WT/pH/6*
unimodal at 0.2 to 0.3; objects with a median area of 3,000 px² or more never fall below 0.18 in either
table, although the illumination differs by a factor of 1.8 between the two experiments (the ratio cancels
it). A lysed cell has no refractive-index difference to the medium left; the debris the author sees in
BSG/pH/6 is the low mode. Most of it is already removed by the size rule. What the contrast rule adds is the
persistent debris that passes it: on BSG/pH/6, 44 tracks with 20 % of the object-frames that survive the
size rule, round (circularity 0.91), shrinking (area ratio last/first 0.74), 1,300 to 1,700 px², tracked for
a median of 31 frames, 70 % of them with a touching "bud" that the lineage would have counted; on WT/pH/6,
4 tracks and 0.1 % of the frames. It cannot tell a dead cell that kept its contrast from a live one.

A fixed threshold failed on W109. The second and third runs applied 0.12 and removed 39 to 68 % of the
object-frames of the W109 chambers, and its complex-medium chambers kept cells in only 83, 28, 77 and 19
frames. The W109 tables show why: those movies were recorded four times darker (cell mean about 480 counts
against 1,500 to 2,300 elsewhere, background about 455), and with that little signal above the camera offset
every live cell sits at a ratio of 0.07 to 0.12, unimodal at 0.095 across all four chambers and all sizes. The
rule is therefore relative: a track is low-contrast when its ratio is below 0.45 times the median of the
size-passing tracks of the same structure (`CELL_MIN_PHASE_CV_REL` in `config.py`; the reference and the
effective threshold per chamber stand in `00_cell_filter.csv`). On the calibration chips this removes the
same tracks as the fixed threshold (the 44 of BSG/pH/6, 2 of WT/pH/6), on W65 static 13 tracks (0.3 % of the
frames), on W109 none. Its limit is a structure whose cells are mostly dead: the reference itself would be
debris, and the log warns when a chamber loses more than a quarter of its frames. In the final run the rule
removes 432 tracks and 0.6 % of all object-frames; 6 chambers lose more than 20 %, the largest share in
BSG/pH/6 (Osc Rep3, 27 %), and W109 loses nothing.

### 2.8 Specific growth rate from budding

With the mother measured, the budding events give a growth rate of the population (`growth_from_budding.py`,
`24_growth_from_budding_*`): µ_bud = births per cell-hour in the sparse window, births being the accepted
events and cell-hours the cells present in every frame of the window times 10 min. In balanced growth each
birth adds one cell, so births per cell-hour are the specific growth rate; a burst of ten buds counts ten, and
no interval between events is needed. The interbud table (`11_specific_growth_rate.csv`, µ = ln(1 + k)/Δt per
mother after Blöbaum 2024) is not a growth rate here: the median interval is 1.3 h, one in ten is a single
frame, and the daughters' generation time never enters, so it comes out at 1.6 to 2.9 h⁻¹. Next to µ_bud the
same table carries the immigration rate, new tracks without a parent mask per cell-hour. On the *WT/pH/6*
table the five oscillation chambers give µ_bud 0.09 to 0.22 h⁻¹ (doubling times 3 to 7 h) and an immigration
rate of 0.15 to 0.28 per cell-hour: washed-in cells arrive faster than cells are born, which is why the count
of objects in a chamber is not a growth curve. µ_bud is the specific birth rate, an upper bound of the net
growth rate; deaths are not countable (a dead cell stays as debris until the cell filter removes it). The
figure against the period (`24_growth_from_budding_vs_period_<osc_type>.pdf`) has the layout of the budding
rate figure, and `24_mu_bud_vs_mu_area.pdf` sets the population rate against the single-cell area growth.

Full run, 535 chambers with a window: 11,361 births in 49,764 cell-hours.

| `24_growth_from_budding_per_chamber.csv` | q10 | median | q90 |
| --- | --- | --- | --- |
| µ_bud, births per cell-hour | 0.10 | 0.21 | 0.37 |
| doubling time, h | | 3.3 | |
| immigration, new tracks without a parent per cell-hour | 0.08 | 0.19 | 0.34 |

- The chamber type does not matter: oscillation chambers 0.22, famine controls 0.21, feast controls 0.19 per
  cell-hour (medians); the feast control exceeds the famine control on 23 of 48 structures (Wilcoxon p 0.91),
  and the bracket is degenerate on 43 of 49. Constant famine does not slow the birth rate in the sparse window.
- The oscillation type does: the Glc series run at 0.22 to 0.30 per cell-hour, the pH series at 0.16 to 0.18
  (medians per series), which is a difference between run days and media as much as between treatments.
- Immigration exceeds births in 37 % of the chambers.
- Against the period (`24_growth_from_budding_control_trend.csv`): no period effect in any series; BSO/Glc
  rises with the period (ρ +1.00) together with its famine control (+0.60), BSA/pH, BSG/Glc and BSG/pH fall
  (−0.60) together with a control (−0.77 to −1.00), the other six series show no trend.
- Population against single cell: µ_bud and µ_area per chip and chamber type correlate at Spearman −0.41
  pooled (n = 146, `24_mu_bud_vs_mu_area.pdf`), but the pooled value is a between-type contrast: the Glc
  structures have the higher µ_bud (median 0.25 against 0.16) and the lower µ_area (0.11 against 0.18 h⁻¹).
  Within Glc the correlation is −0.24 (n = 87; oscillation chambers alone −0.43), within pH +0.39 (n = 59;
  oscillation chambers +0.48). The trade-off "growth into size or into blastoconidia" holds between the two
  oscillation types and in the static comparison of section 3, not within the pH series.
- Static (`static/24_growth_from_budding_summary.csv`): W109 minimal medium 0.18, complex medium 0.03 per
  cell-hour; W65 0.11 and 0.14, with the burst chamber at 0.28. PKO: its famine controls bud at 0.32 to 0.43
  per cell-hour, its feast controls at 0.18 to 0.22 (`pko/24_growth_from_budding_per_chamber.csv`).

---

## 3. Section 1: static medium

Both static families are one chip each (`STATIC_SINGLE_CHIP_FAMILIES`). W109 is one chip with three
structures; its four movies per medium are four distinct chambers of that chip (no duplicates: 186 to 410
tracks in minimal medium), so its p-values below describe chamber scatter within one culture, not cultures.
W65 was one chip with one pre-culture; of its five minimal-medium recordings, Rep1 and Rep5 are the same
stage position recorded twice (`EXCLUDED_CHAMBERS` drops Rep5), so four chambers remain. Neither family has
biological replication; the medium effect rests on chambers of one chip in each family.

| W109, complex (ypd) vs minimal (omlp) | omlp | ypd | Welch p over the four chambers |
| --- | --- | --- | --- |
| endpoint cell area, px² (`static/13_endpoint_summary.csv`) | 5,489 | 19,105 | 0.004 |
| endpoint eccentricity | 0.80 | 0.75 | 0.23 |
| buds per mother-hour, sparse window (`static/21_budding_rate_summary.csv`) | 0.33 | 0.05 | 0.12 |

Cells in complex medium end up 3.5 times larger, marginally rounder, and bud a sixth as often per
mother-hour (0.45, 0.58, 0.27 and 0.00 against 0.03 to 0.08; the fourth minimal-medium chamber had no bud in
its window): growth goes into cell size rather than into blastoconidia, and the population birth rate says the
same (0.18 against 0.03 per cell-hour, section 2.8). The v11 tables gave the same picture (5,049 against
14,451 px², 0.47 against 0.09 per mother-hour). Two intermediate runs with a fixed contrast threshold had
halved the effect because that threshold removed live cells from the darker W109 recording (section 2.7); the
structure-relative rule leaves W109 untouched.

The W65 chip does not repeat it. With the duplicate dropped, the area ratio ypd/omlp is 0.95, the
eccentricity 0.72 against 0.43 (Welch p 0.006 over chambers, the opposite direction to W109), the budding rate
0.18 against 0.14 per mother-hour. The four W65 minimal-medium chambers are not one population: Rep1 ends at
40,000 px² (215 µm²), Rep2 to Rep4 at 13,000 to 16,000 px², a spread of 61 % against 4 to 24 % in the other
groups (`static/13_endpoint_per_chamber.csv`). Rep1 holds two swollen cells of about 50,000 px²
each, and the overlay of frames 85 to 88 shows the lower one releasing a ring of ten blastoconidia within
three frames (tracks 3 to 14): the swollen-cell-to-blastoconidia transition of Rensink et al. 2026, caught
in one chamber, and the reason the chamber's endpoint area and budding rate stand apart. In complex medium
the two chip families agree (W65 19,963 px², W109 19,105); in minimal medium W65 is 3.8 times larger, and
that difference is largely the one chamber. With one culture per family, the medium effect
is a W109 result with W65 as a second, single chip that shows the same cell size in complex medium and a
different picture in minimal medium. The static figures (`static/13_endpoint_vs_medium_*.pdf`,
`static/21_budding_rate_vs_medium.pdf`, `static/24_growth_from_budding_vs_medium.pdf`,
`static/12_area_growth_rate_all.pdf`) draw every chamber as a small point next to its mean, so the W65 spread
and the burst chamber are visible in the figure itself; the per-chamber tables behind them are
`static/13_endpoint_per_chamber.csv`, `static/21_budding_ratio_per_experiment.csv` and
`static/24_growth_from_budding_per_chamber.csv`.

- Main text: `static/13_endpoint_vs_medium_area.pdf`, `static/13_endpoint_vs_medium_eccentricity.pdf`,
  `static/21_budding_rate_vs_medium.pdf` (mean ± SEM, one panel per chip family; the error unit per family
  stands in `error_unit`: chips for W109, chambers for W65).
- Appendix: `static/10_cell_area_over_time.pdf` (the ypd chambers saturate, which fixes the endpoint window,
  `static/13_endpoint_saturation_per_chamber.csv`), `static/21_panel_a_violin.pdf`,
  `static/12_area_growth_rate_all.pdf`, `static/40_Rp_*.pdf`.

---

## 4. Section 2: oscillations, the apparent dose responses

The pipeline reads one value per structure and plots it against the period with the controls of the same
structure as markers (`13_endpoint_vs_period_<readout>_<osc_type>.pdf`, `21_budding_rate_vs_period_<osc_type>.pdf`;
since 2026-10-05 without connecting lines, pooled-control lines and the bracket-score row, `PERIOD_FIGURE_*` in
`config.py`; the bracket score, 0 = famine control, 1 = feast control, undefined when the two controls do not
separate, stays in `*_bracket_score.csv`). Spearman over structures, n = number of periods (4 to 6), is
an effect size, not a test (`13_endpoint_spearman.csv`, `21_budding_rate_spearman.csv`,
`12_area_growth_rate_spearman.csv`).

What appears in the v12 run:

| readout | series | Spearman vs period |
| --- | --- | --- |
| buds per mother-hour | WT/pH | −1.00 |
| | WT/Glc, BSA/pH, BSG/pH | −0.60 |
| | BSO/Glc | +0.77 |
| µ_bud, births per cell-hour | BSO/Glc | +1.00 |
| | BSA/pH, BSG/Glc, BSG/pH | −0.60 |
| `ratio_OxPro` | BSO/pH, BSO/Glc | +1.00, +0.77 |
| `ratio_pHluorin` | BSPH/Glc, BSPH/pH | +0.71, +0.60 |
| `ratio_GlyRNA` | BSG/Glc | +0.89 |
| area | BSPH/Glc, WT/pH | −0.89, −0.80 |
| | BSA/Glc | +0.71 |
| µ_area | BSA/Glc, BSO/Glc, WT/pH | −0.94, −0.83, −0.80 |
| eccentricity | 7 of 10 series positive | 0.60 to 0.83 in six of them |

The budding trend of the v11 run is gone. There, eight of the nine budding-rate series with a non-zero ρ
were negative (WT/Glc −0.90, BSPH/Glc −0.89, BSG/Glc and BSO/Glc −0.66); the v12 run gives the same series
−0.60, −0.26, −0.20 and +0.77, and across the ten series the per-series ρ of the two runs correlate at 0.05,
while the budding rates of the structures themselves correlate at 0.77 (`docs/tracking_plan.md` 6g). Which
series looked like a dose response was decided by the tracker, whose fragmentation depends on the density
and shape of the cells in each structure, not by the cells. The sensor and shape trends, which need no
tracking, are the same in both runs; what they are is the subject of section 5.

---

## 5. The pivot: the controls carry the trends

Constant-medium chambers cannot respond to the period. On the same structures they trend with it:

| readout | series | Osc vs period | strongest control vs period | Osc minus controls |
| --- | --- | --- | --- | --- |
| buds per mother-hour | WT/pH | −1.00 | −1.00 (PosCtrl) | +0.40 |
| | BSO/Glc | +0.77 | +0.66 (NegCtrl) | +0.77 |
| | BSA/pH | −0.60 | −1.00 (PosCtrl) | +0.20 |
| | BSG/pH | −0.60 | −0.80 (PosCtrl) | −0.20 |
| | WT/Glc | −0.60 | 0.00 | −0.80 |
| µ_bud, births per cell-hour | BSO/Glc | +1.00 | +0.60 (NegCtrl) | +0.09 |
| | BSA/pH | −0.60 | −1.00 (PosCtrl) | +0.20 |
| | BSG/Glc | −0.60 | −0.77 (NegCtrl) | +0.14 |
| | BSG/pH | −0.60 | −0.80 (PosCtrl) | −0.20 |
| `ratio_OxPro` | BSO/Glc | +0.77 | +0.94 (NegCtrl) | +0.26 |
| `ratio_OxPro` | BSO/pH | +1.00 | +1.00 (NegCtrl) | +0.40 |
| `ratio_pHluorin` | BSPH/Glc | +0.71 | +0.89 (NegCtrl) | −0.09 |
| `ratio_pHluorin` | BSPH/pH | +0.60 | +1.00 (NegCtrl) | −0.40 |
| `ratio_GlyRNA` | BSG/Glc | +0.89 | +0.43 | −0.26 |
| area | BSA/Glc | +0.71 | +0.94 (mean of both) | +0.60 |
| area | BSPH/Glc | −0.89 | −0.89 (NegCtrl) | −0.26 |
| area | WT/pH | −0.80 | −0.80 (NegCtrl) | +0.40 |
| µ_area | BSA/Glc | −0.94 | −0.60 (NegCtrl) | −0.60 |
| µ_area | BSO/Glc | −0.83 | −0.49 | −0.77 |
| µ_area | WT/pH | −0.80 | −1.00 (NegCtrl) | +0.60 |

- `13_endpoint_control_trend.csv`, `12_area_growth_rate_control_trend.csv`, `21_budding_rate_control_trend.csv`
  (collected in `50_control_trend_summary.csv`): per readout and series, the Spearman of the oscillation
  chambers, of PosCtrl, NegCtrl and their mean, of the difference, and a verdict. A `period effect` needs all
  three: the oscillation chambers trend (|ρ| ≥ 0.6), no control trends the same way, and the difference
  trends. Over the 58 readout × series rows of the full run (area, eccentricity, four sensor ratios, µ_area,
  budding rate, µ_bud):

  | verdict | rows | which |
  | --- | --- | --- |
  | no monotone trend of the oscillation chambers | 31 | |
  | structure effect: a control trends the same way | 23 | 4 of them with a residual difference: area and µ_area BSA/Glc, eccentricity BSG/Glc, budding rate BSO/Glc |
  | not robust: the trend vanishes after subtracting the controls | 2 | eccentricity BSPH/pH; GlyRNA BSG/Glc |
  | period effect: all three conditions | 2 | µ_area BSO/Glc; budding rate WT/Glc |

  The two are scattered: two readouts, two strains, and neither recurs in the other oscillation type of the
  same strain (µ_area BSO/pH −0.20, budding rate WT/pH a structure effect). With four periods a |ρ| of 0.8
  arises by chance in one series of three, and the count of period effects fell with every change of the
  pipeline: 7 in the v11 run, 5 with the new analysis rules on the v11 tables, 3 in the first v12 run, 2 once
  the contrast rule removed the dead cells (area BSG/pH went from ρ −0.8 to −0.4, one rank swap among four
  periods); three further rows changed verdict between the fixed and the relative contrast rule, none of them
  a period effect. What is left is what 58 tries produce by chance. The one recurrence of the v11 run, µ_area falling
  with the period in all five pH series, did not survive the rebuild (ρ −0.4 to +0.4).
- `13_endpoint_within_culture.csv`, `24_growth_from_budding_within_culture.csv`: the change from the shortest
  to the longest period inside one culture, where the pre-culture cannot differ. Area: the oscillation
  chambers fall in 8 of 19 cultures, the controls in 13 of 19, the difference in 6 of 19; µ_bud: 8, 8 and 11
  of 19. No direction.
- `13_endpoint_bracket_score.csv`, `21_budding_rate_bracket_score.csv`, `24_growth_from_budding_bracket_score.csv`:
  the feast and famine controls barely separate. For the budding rate the bracket is degenerate on 44 of 49
  structures, and PosCtrl exceeds NegCtrl on 26 of 48 (Wilcoxon p 0.61); for µ_bud on 43 of 49 and 23 of 48
  (p 0.91). Cells under constant famine bud in the sparse window as often as under constant feast, 0.05 to
  0.63 against 0.00 to 0.94 per mother-hour. Area and eccentricity have a direction, feast cells larger and
  more eccentric on 31 and 33 of 48 structures (Wilcoxon p 0.01 each), but not a separation beyond the
  chamber scatter: the bracket is degenerate on 43 and 39 of 49 structures.
- `40_control_consistency_*.pdf`, `40_control_consistency_*_kruskal.csv`, `11_control_consistency_mu_kruskal.csv`:
  the same controls compared across the structures of a series; they differ.

The mechanism: the structure sets the readout for all its chambers, and the shortest period usually sat on the
same structure position, so position and period were confounded within every culture. No run sheet records the
position, so it stays inferred from the period order. The day blocks add a culture effect on top, balanced
across strains for Glc, aligned with the period for pH (`00_chip_run_order.csv`). Without the controls every
drift in section 4 would have been reported as a dose response.

`50_control_trend_summary_growth.pdf` and `_sensors.pdf` (with `50_control_trend_summary.csv`) show the whole
finding in two figures, the growth and morphology readouts and the sensor ratios:
one point per readout and strain series, the oscillation Spearman on x, the strongest control Spearman on y,
coloured by verdict; points along the diagonal are the structure effects. It collects the control-trend tables
of the endpoint (`13_`), of µ_area (`12_area_growth_rate_control_trend.csv`) and of the budding rate (`21_`).
On the full run, 23 of the 58 points sit in the structure-effect corners on the diagonal and 31 in the
middle band; 2 period-effect points and 2 not-robust points remain off the diagonal.

---

## 5b. Robustness: the oscillation against constant medium, and R(t) / R(p)

The thesis question is robustness under oscillation (introduction 1.4). Two readings of it. Since 2026-10-01 the
pipeline writes both: step 40 runs every robustness metric per chamber through the chip logic of the readouts
(`40_<metric>_<readout>_per_chip.csv`, `_spearman`, `_bracket_score`, `_control_trend`, `_vs_period_<osc_type>.pdf`
with the controls of the structure), step 50 collects the verdicts (`50_robustness_control_trend_summary.*`) and
writes the paired comparison for all readouts and metrics (`51_osc_vs_controls_per_structure.csv`,
`51_osc_vs_controls.csv`, `51_osc_vs_controls_vs_period.csv`, `51_osc_vs_controls*.pdf`; `osc_vs_controls.py`).
The numbers below come from `docs/scratch/robust3.py`, the same functions on the final-run tables, until the next
run:

**Oscillation chambers against the controls of their structure**, paired over the 49 structures (Wilcoxon of the
oscillation mean against the control mean; ratio Osc / control mean):

| readout | Glc (29 structures): Osc > ctrl | ratio | p | pH (20): Osc > ctrl | ratio | p |
| --- | --- | --- | --- | --- | --- | --- |
| buds per mother-hour | 27 | 1.25 | < 0.001 | 13 | 1.09 | 0.09 |
| µ_bud | 23 | 1.15 | < 0.001 | 14 | 1.05 | 0.12 |
| endpoint area | 25 | 1.14 | < 0.001 | 10 | 1.00 | 0.65 |
| eccentricity | 7 | 0.98 | < 0.001 | 7 | 0.98 | 0.45 |
| µ_area | 7 | 0.88 | 0.010 | 8 | 0.97 | 0.33 |
| immigration | 16 | 1.02 | 0.90 | 13 | 1.10 | 0.45 |

Under Glc oscillation the cells bud a quarter more often, are 14 % larger, marginally rounder and grow 12 % slower
in area than under constant feast or famine on the same structure; under pH oscillation nothing differs. The
budding rate lies above both controls on 26 of 49 structures and below both on 4; every strain shows it (ratio
1.07 BSG to 1.37 BSA). The ratio does not depend on the period (|ρ| ≥ 0.6 in 2 to 5 of 10 series per readout,
both signs). Immigration is the same in oscillation and control chambers, so the flow is comparable; but the
oscillation chambers are the positions A3 to A12 and the controls A1/A2 and A13/A14, so position is confounded
with treatment here too. Note also that the famine controls (0 g/L) bud as often as the feast controls in the
sparse window (section 5), while in the BioLector the biosensor strains do not grow at 0 g/L within 24 h.

**R(t) and R(p)** (Trivellin 2022 / Blöbaum 2024, R = −σ²/x̄ · 1/m ≤ 0), medians over structures, oscillation
against feast / famine controls, structures with the oscillation chambers less robust than the control mean:

| metric | readout | Osc | feast | famine | less robust of 49 | p | period effects of 10 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| R(t) population | area | −0.127 | −0.096 | −0.123 | 31 | 0.013 | 1 (WT/Glc) |
| R(t) population | eccentricity | −0.015 | −0.012 | −0.014 | 33 | 0.019 | 0 |
| R(t) single cell (≥ 10 frames) | area | −0.106 | −0.100 | −0.077 | 35 | 0.007 | 1 (BSG/pH) |
| R(t) single cell | eccentricity | −0.022 | −0.020 | −0.018 | 30 | 0.010 | 0 |
| R(p) | area | −0.671 | −0.571 | −0.613 | 34 | 0.003 | 2 (BSA/Glc, BSPH/pH) |
| R(p) | eccentricity | −0.048 | −0.041 | −0.042 | 33 | 0.003 | 1 (BSA/Glc) |
| R(p) | µ_area | −0.333 | −0.412 | −0.453 | 16 | 0.008 | 1 (BSA/Glc) |
| R(p) | budding rate per mother (across the mothers of a chamber) | −2.018 | −1.984 | −1.921 | 25 | 0.95 | 1 (BSA/Glc) |
| R(t) single cell (≥ 3 intervals) | µ_event per mother | −1.115 | −1.132 | −0.999 | 29 of 47 | 0.022 | 1 (BSA/Glc) |

Under oscillation the cell area is less stable over time and more heterogeneous across cells, the growth rate more
homogeneous; the two heterogeneities are unrelated across chambers (ρ 0.06, n 497). µ_bud itself has one value
per chamber and no R; its robustness is carried by the two budding rows: the mothers of a chamber are as
heterogeneous in their budding rate under oscillation as under constant medium (p 0.95), the rhythm of the single
mother is less stable under Glc oscillation (p 0.043; pH 0.28). Sensors: no difference between oscillation and
control chambers for any of the four ratios (p 0.06 to 0.92); pHluorin the most robust (R −0.01), OxPro the
least (R(t) −0.21, R(p) −0.73, values near zero). The control-trend classification over all 114 robustness
combinations (21 metrics x series): 60 no trend, 35 structure effects, 10 not robust, 9 period effects, five of
them in BSA/Glc (R(p) of area, eccentricity, µ_area and budding rate per mother, R(t) of µ_event: with longer
periods that series gets more homogeneous in size, growth and budding, more stable in the budding rhythm, and more
heterogeneous in shape), the same picture as the 58 readout rows of section 5. The R(t) caveat: the cycle itself
contributes to the temporal variance of the oscillation chambers, and differently per interval (a random phase of
the cycle at intervals up to 6 min, near the Nyquist limit at 12 min, a resolved cycle at 24 min), so an R(t)
trend of the oscillation chambers alone is not interpretable; against the controls, which see no cycle, it stays
readable.

Static context: R(t) of the area −0.18 in minimal and −0.39 in complex medium (the growing cells), R(p) −0.47 and
−0.51.

---

## 6. Section 3: biosensors

On the real data none of the sensors shows a clear feast versus famine difference between its PosCtrl and
NegCtrl chambers (`95_*_comparison.pdf`, author's reading; bracket degenerate on 6 to 9 of 10 structures per
sensor, `13_endpoint_bracket_score.csv`), so the sensors' dynamic range in this setup is not established.
Beyond that, the sensors report differences between structures, and their famine controls report the same
differences. `ratio_OxPro` rises with the period in the NegCtrl chambers exactly as in the oscillation
chambers (Glc +0.94 against +0.77, pH +1.00 against +1.00), `ratio_pHluorin` likewise (Glc +0.89 against
+0.71, pH +1.00 against +0.60); `ratio_GlyRNA` rises in the oscillation chambers of BSG/Glc (+0.89) but not
after its controls are subtracted; `ratio_Queen-2m` shows no trend. Nothing the sensors show is attributable
to the period.

- Main text: `13_endpoint_vs_period_ratio_<sensor>_Glc.pdf` (and `_pH.pdf`), the ratio rows of
  `13_endpoint_control_trend.csv`.
- Whether a sensor responds to feast versus famine at all, independent of the period:
  `95_<strain>_<osc_type>_ratio_<sensor>_comparison.pdf` (PosCtrl against NegCtrl per structure),
  `95_..._control_chambers.pdf` (every control chamber within every structure), `95_..._timeseries.pdf`,
  `95_..._raw_channels.pdf`, with `95_..._summary_per_replicate.csv`: no clear difference for any sensor.
- Appendix: `30_<channel>_over_time_<osc_type>.pdf`, `31_ratio_<sensor>_over_time_<osc_type>.pdf` (drift over
  hours; single cycles are not resolved for intervals up to 6 min), `40_Rp_ratio_*.pdf`, `40_Rt_population_ratio_*.pdf`.

---

## 7. PKO

One structure at period 3, built for the clogging hypothesis: without pullulan the control chambers of a
structure should agree better than the producers' do. They do not (`pko/60_pko_within_chip_agreement.csv`,
`pko/61_pko_control_agreement.pdf`): the PKO chamber-to-chamber CV of the area level is 0.23 for famine and
0.37 for feast, the 67th and 96th percentile of all producer structures; the detrended temporal CV sits at the
92nd and 61st percentile. One observation, pointing the wrong way for the hypothesis. PKO famine cells are 1.5
times larger than PKO feast cells, where the producers' median ratio is 0.98 (PKO at the 92nd percentile,
`level_mean` in the same table); the µ_area bracket of feast over famine (`pko/60_pko_control_bracket.csv`)
is 0.10, the producers' median. Its famine controls also bud faster than its feast controls (0.32 to 0.43
against 0.18 to 0.22 births per cell-hour, `pko/24_growth_from_budding_per_chamber.csv`). Unexplained at
n = 1.

---

## 8. Figure plan

| thesis section | main-text figures | appendix figures | evidence tables |
| --- | --- | --- | --- |
| Methods: units and sampling | none (a schematic of chip, structures, arrays and chambers is the author's) | `00_n_tracks_overview_summary.pdf` | `00_chip_overview.csv`, `00_chip_run_order.csv` |
| 3.2.2 Observed morphology | `90_morphology_scatter.pdf` (one cell = one point, mean area against mean eccentricity, guides at 30 µm² and eccentricity 0.6, share of large round cells in the title; static per chip family and medium) | | `40_Rt_single_cell_area.csv` / `_eccentricity.csv` (per-cell means), `20_bud_size_at_appearance.csv` |
| Methods: tracking and sparse window | `20_lineage_window.pdf`, `00_new_objects_vs_density.pdf` | `20_bud_size_at_appearance.pdf`, `qc_comparison/70_qc_effect.pdf` and `lineage_validation/lv_02_detection_rate.pdf` (v11 run) | `00_track_fragmentation.csv`, `00_cell_filter.csv`, `20_bud_size_threshold.csv`, `00_new_objects_vs_density.csv` (the block table of 2.2 for all chambers) |
| 1 Static medium | `static/13_endpoint_vs_medium_area.pdf`, `static/21_budding_rate_vs_medium.pdf` | `static/13_endpoint_vs_medium_eccentricity.pdf`, `static/10_cell_area_over_time.pdf`, `static/21_panel_a_violin.pdf` | `static/13_endpoint_summary.csv`, `static/13_endpoint_per_chamber.csv`, `static/21_budding_rate_summary.csv` |
| 2 Oscillations | `21_budding_rate_vs_period_Glc.pdf`, `24_growth_from_budding_vs_period_Glc.pdf`, `13_endpoint_vs_period_area_Glc.pdf` | the `_pH.pdf` counterparts, `24_immigration_vs_period_*.pdf`, `24_mu_bud_vs_mu_area.pdf`, `13_endpoint_vs_period_eccentricity_*.pdf`, `10_cell_area_over_time_*.pdf`, `12_area_growth_rate_all.pdf` (with the control chambers), `40_Rp_*.pdf` | `13_endpoint_spearman.csv`, `21_budding_rate_spearman.csv`, `24_growth_from_budding_spearman.csv`, `12_area_growth_rate_spearman.csv` |
| 3 Biosensors | `13_endpoint_vs_period_ratio_OxPro_Glc.pdf`, `13_endpoint_vs_period_ratio_pHluorin_Glc.pdf` | `95_*_comparison.pdf`, `31_ratio_*_over_time_*.pdf` | ratio rows of `13_endpoint_control_trend.csv`, `95_*_summary_per_replicate.csv` |
| 4 Controls (the pivot) | `50_control_trend_summary_growth.pdf`, `50_control_trend_summary_sensors.pdf`; `13_endpoint_vs_period_ratio_OxPro_Glc.pdf` as the worked example | `40_control_consistency_*.pdf` | `*_bracket_score.csv` of `13_`/`21_`/`24_`, `50_control_trend_summary.csv`, `13_endpoint_control_trend.csv`, `12_area_growth_rate_control_trend.csv`, `21_budding_rate_control_trend.csv`, `24_growth_from_budding_control_trend.csv`, `13_endpoint_within_culture.csv`, `*_bracket_score.csv`, `00_chip_run_order.csv` |
| PKO | `pko/61_pko_control_agreement.pdf` | `pko/13_endpoint_vs_period_area_Glc.pdf` | `pko/60_pko_within_chip_agreement.csv`, `pko/60_pko_control_bracket.csv` |

Suggested order of the results chapter: static first (the clean result), then the tracking limit as a short
methods-in-results block, then oscillations and biosensors together with the control pivot, then PKO. Section 4
is not a separate finding after the others; it is what the others turn into. The change of the budding trends
between the two trackings (section 4) belongs in the methods block: it is the direct evidence that the
per-series trends of a lineage readout are tracker properties.

All figures share one style (`plot_style.py`, colours in `config.py`): colour means strain (WT blue, BSA
orange, BSO green, BSG yellow, BSPH pink, PKO grey), the control chambers are markers in the strain colour
(feast control filled up-triangle, famine control hollow down-triangle, oscillation filled circle), periods in
the time-series figures run from light (short) to dark (long) within the strain hue, and the static media are
filled (complex) versus lighter (minimal), with every chamber as a small point beside its mean. The figures
carry no explanatory footers; what a band or a marker means goes into the caption.

---

## 9. Open points and proposed additions

Answered so far: the thesis is in English and the four sections stay in the results chapter; W65 was one chip
with one pre-culture, and its minimal-medium Rep5 is a second recording of Rep1; W109 is one chip with three
structures; the sensors show no clear feast/famine difference on the real data; there is no run sheet of
structure positions, so the outlook asks for one; the small objects of BSG/pH/6 are dead-cell debris; the W65
outlier chamber holds two swollen cells, one of which releases ten blastoconidia at once (section 3).

Resolved since: the contrast rule on W109 (the W109 movies are a darker recording, the rule is now relative to
the structure, section 2.7, and the final run carries it). Still open: nothing that needs the author's input.

Built since the first draft: the control-trend summary figure (`50_control_trend_summary.pdf`), the control-trend
check for µ_area, the rebuilt segmentation and tracking with the measured mother, the cell filter, the
size criterion on a bimodal distribution, the new-objects-against-density table (2.2, from the v12 table of
WT/pH/6), the dead-cell rule (2.7).

Not built, still possible:

- **Nothing pending on the pipeline side.** The robustness outputs of section 5b (steps 40 and 50, built
  2026-10-01) appear with the next run, together with the density figure, the per-chamber points and the per-cell
  morphology figure (`90_morphology_scatter.pdf`, rebuilt 2026-10-03 for results section 3.2.2: 22,090 cells of
  at least 10 frames, median 22 µm², eccentricity 0.75; 5 % large and round; buds 4.7 µm² at first detection on
  mothers of 71 µm²; no filamentous objects: solidity 0.98, axis ratio > 3 in < 0.3 %). Every number in this document comes from the final run. The two
  figure additions of 2026-09-29, the per-chamber points on the static figures (section 3) and the block table
  of 2.2 as a figure over all chambers (`00_new_objects_vs_density.pdf`), are in the code and appear with the
  next run; the density numbers quoted in 2.2 are from the four tables at hand, not from the full run.

---

## 10. Reproduction

```
cd analyse_pipeline
python run_analysis.py            # RESULTS_VERSION = "v12" in config.py: analysis_output_v12/ with static/ and pko/
python validate_lineage.py        # lineage_validation/ (needs the full run first; pending on v12)
```

Segmentation, tracking and the cluster route: `imaging/README.md` (`segment_all.py`, `merge_results.py`).
`config.py` holds every threshold named above (`CELL_MIN_FRAMES`, `CELL_MIN_MAX_AREA_PX`, `CELL_MIN_PHASE_CV_REL`,
`EXCLUDED_CHAMBERS`, `STATIC_SINGLE_CHIP_FAMILIES`, `LINEAGE_PARAMS`, `LINEAGE_SPARSE_*`, `DENSITY_BLOCK_FRAMES`,
`AREA_GROWTH_MIN_FRAMES`,
`BUD_MAX_AREA_FRACTION_FALLBACK`, `FLAG_EXCLUDE_ROWS`) and the figure colours (`STRAIN_COLORS`), and logs them,
together with the method caveats, at the start of every run. `README.md` describes each module and output prefix.
