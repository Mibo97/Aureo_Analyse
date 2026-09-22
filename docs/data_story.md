# The story of the data

Argumentation of the thesis results, with the figures and tables of the pipeline that carry each step.
Draft of 2026-09-22 (English; can be translated). Paths are relative to `analysis_output/`.
Numbers marked *QC batch* come from the manually curated batch WT/pH/6 and its raw `Combined_Results.csv`;
everything else comes from the full run. Statements marked **(open)** need the author's input, see section 9.

---

## 0. In one paragraph

An oscillation experiment across five strains, two media and six feast/famine periods produced
dose responses in the budding rate, in cell shape and in three biosensor ratios. Every one of them
turned out to be a property of the microfluidic structure, not of the period: the constant-medium
controls on the same structures drift identically, and after subtracting them nothing is left that
recurs across strains. Two properties of the design made this inevitable, and one decision made it
visible. Each period sat on its own structure, usually in the same position on the chip, so period
and position were aligned; the periods were shorter than the sampling interval, so only cumulative
effects could ever be read; and every structure carried its own feast and famine controls, which is
the only reason the drift could be told from a treatment effect. Underneath, the tracker fragments
every cell track into pieces of a few frames, so budding can only be counted in the sparse first
hours of a chamber. The one section with independent cultures, complex versus minimal medium on
the W109 chips, shows a coherent medium effect: larger, rounder, less budding cells in complex medium.

---

## 1. What was measured, and what the units are

### 1.1 Design

| element | value |
| --- | --- |
| organism, imaging | *Aureobasidium pullulans*, microfluidic chambers, one frame per 10 min, 130 to 136 frames (about 22 h), Cellpose segmentation with tracking |
| strains | WT and four biosensor strains: BSA (`ratio_Queen-2m`), BSG (`ratio_GlyRNA`), BSO (`ratio_OxPro`), BSPH (`ratio_pHluorin`); PKO = no pullulan |
| oscillation types | Glc and pH switching; periods 0.75, 1.5, 3, 6, 12, 24 min (Glc; WT lacks 24) and 0.75, 1.5, 6, 24 min (pH) |
| one structure | one period; 5 oscillation chambers, 3 constant-feast (PosCtrl) and 3 constant-famine (NegCtrl) chambers on several arrays; `replicate` is the array index, `chamber` the position (A1/A2 feast, A13/A14 famine, A3 to A12 switching) |
| one physical chip | 2 or 3 structures = 2 or 3 periods, one pre-culture, one day (`culture`) |
| oscillation and PKO total | 50 structures, 547 chambers, 24 cultures |
| static | W109: 4 chips per medium; W65: 5 chips, each with both media; `chamber` always A0; every chip its own culture |

The pipeline calls a structure `chip` in every table and figure. `culture` is the physical chip.

- `00_chip_overview.csv`: one row per structure with its chambers per control type, its culture and `n_structures_in_culture`.
- `00_chip_run_order.csv`: for each strain and oscillation type, the periods in run-date order and Spearman(period, date). The periods were run in day blocks: Glc as 3 + 3 periods on two days, pH as 2 + 2.
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

Ten-minute frames resolve periods of 20 min and longer. Five of the six Glc periods and three of the four pH
periods lie below that. Single cycles are not observable; apparent periodicity in the time series is an
alias. The pipeline therefore reads only cumulative quantities: the endpoint of a chamber (`13_`), the budding
rate over the sparse window (`21_`), and the robustness R(p) (`40_Rp_*`). The time series (`10_`, `30_`,
`31_`) are shown as drift over hours, never as cycles. (`config.py` logs the unresolved periods at every run.)

---

## 2. The tracking limit and the sparse-phase window

### 2.1 The tracker fragments

Every object receives a new track ID with a probability of about 12 % per frame, on every chip.

| `00_track_fragmentation.csv`, per chamber, before gap closing | q10 | median | q90 |
| --- | --- | --- | --- |
| new tracks per object and frame | 0.09 | 0.14 | 0.22 |
| median track length, frames | 1 | 3 | 4 |
| share of single-frame tracks | 0.16 | 0.31 | 0.51 |
| share of object-frames in tracks of 10 frames or more | 0.48 | 0.64 | 0.77 |

*QC batch*: 9 063 tracks in 11 chambers, 2 712 of them one frame long.

### 2.2 The heuristic budding events follow the density

Chambers start with 1 to 3 objects and fill to as many as 180. Because every new track next to an
established cell is a candidate bud, the number of events follows the object count.

*QC batch*, raw events per chamber by 22-frame block: 0.7, 4.1, 13.7, 51.8, 148.4, 243.4, in total 462 per
chamber. Median objects per frame in the same blocks: 6, 8, 18, 49, 106, 149. These are fragment statistics,
not biology. Evidence in the pipeline: `20_lineage_window.pdf` (objects per frame over time for every chamber)
and `00_track_fragmentation.csv` (`objects_first_frame`, `objects_last_frame`). A figure of events against
density over the full run is proposed in section 9.

### 2.3 Manual quality control does not fix it

*QC batch*: 503 rows of manual QC, 248 merges after tracking breaks ("cell movement"), 214 background objects,
16 mother cells on the edge, 13 dead cells, 11 debris. Applied, they change 2.6 % of the raw events over the
run and 19 % inside the sparse window. The 248 merges reassemble 96 cells from a median of 3 fragments each,
out of thousands of fragments.

- `qc_comparison/70_qc_effect.pdf`, `70_qc_effect_summary.csv`, `70_qc_effect_per_chamber.csv`: with QC against
  without QC, chamber by chamber, for endpoint area, µ_area, R(p) and the budding ratio.
- `qc_comparison/70_qc_exclusion_inventory.csv`, `70_qc_exclusions_conflicts.csv`: what the QC file contains.

Consequence for the thesis: repeating this QC for 50 batches would not change the lineage readouts. The
automatic steps below replace the bulk of it; the manual file keeps precedence where it exists.

### 2.4 Sparse-phase window and gap closing

- **Window** (`20_lineage_window.csv`, `20_lineage_window.pdf`): per chamber, the frames before the rolling
  median of objects per frame exceeds 20; chambers with fewer than 20 such frames drop out. Full run: all
  565 chambers qualify, median window 110 of 133 frames, 211 chambers never exceed 20 objects, range 27 to
  133 frames. Every lineage output (`20_` to `23_`, `21_budding_rate_*`, µ_event, mother/bud labels) uses
  only this window; endpoint, µ_area and robustness still see the whole run.
- **Gap closing** (`00_track_relinks.csv`, `00_track_relinks_per_chamber.csv`): a new track is joined to a track
  that ended at most 2 frames earlier when the distance is within 150 px, the area ratio within 0.5 to 2, and
  the match is unique in both directions. Manual merges run first and win. *QC batch*: 40 % of the 70 manual
  links inside the window are recovered, the rest are refused as ambiguous (same-cell jumps of median 97 px
  with several candidates). Wider radii or longer gaps recover fewer, because ambiguity grows faster than hits.
- **Events in the window**, *QC batch*: 158 raw, 136 with gap closing, 128 with manual QC, 116 with both;
  118 of the 136 have a mother tracked for 30 frames or more. Ten to fifteen per chamber is the number the
  growth in the window (2 to 20 objects in about 11 h) allows.
- **Validation** (`lineage_validation/lv_02_detection_rate.pdf`, `lv_02_detection_rate_per_chamber.csv`,
  `lv_03_detection_rate_kruskal.csv`): the share of new tracks in the window that are assigned to a mother
  differs between the structures of a series (Kruskal p < 0.05 in 7 of 10 series) and falls with the period in
  BSG/Glc (rho −0.94) and BSPH/Glc (−0.77). Inflowing cells are new tracks without a mother; a washed-in cell
  that lands next to a mother still counts. `lv_01_d_over_r_distribution.pdf`, `lv_03_assignment_ambiguity.pdf`
  and `lv_04_tolerance_sweep.pdf` document the spatial rule.

### 2.5 The size criterion that was tried and dropped

A bud should be small at first detection and a washed-in cell mother-sized. On the real data the ratio of bud
area to mother area has one broad mode around 0.4 to 0.5 and no second peak (`20_bud_size_at_appearance.pdf`,
`20_bud_size_threshold.csv`, column `source`). A fixed cut at 0.5 removed half of the events and made the
detection rate more, not less, structure dependent. The criterion is inactive; the mother's mask does not
change when a bud is first segmented **(author statement)**, so the area balance of the mother is not a
criterion either. The diagnostic columns stay in `20_bud_size_at_appearance.csv`.

### 2.6 What remains trustworthy

Readouts that need no tracking: endpoint area and eccentricity per chamber (`13_`), population R(t) and R(p)
(`40_Rt_population_*`, `40_Rp_*`), sensor ratios (`30_`, `31_`, `95_`). µ_area (`12_`) fits ln(area) per track
and now requires 10 frames instead of 2. The lineage readouts describe the sparse phase only.

---

## 3. Section 1: static medium, the section with replication

Four (W109) and five (W65) independent chips per medium; the unit is the chip, the error bar a chip error bar.

| W109, complex (ypd) vs minimal (omlp) | omlp | ypd | Welch p over chip means |
| --- | --- | --- | --- |
| endpoint cell area, px² (`static/13_endpoint_summary.csv`) | 5 049 | 14 451 | 0.004 |
| endpoint eccentricity | 0.79 | 0.74 | 0.06 |
| buds per mother-hour, sparse window (`static/21_budding_rate_summary.csv`) | 0.47 | 0.09 | 0.09 |

Cells in complex medium end up almost three times larger, slightly rounder, and bud a fifth as often per
mother-hour: growth goes into cell size rather than into blastoconidia.

The W65 chips do not confirm it. Every direction flips (area ratio ypd/omlp 0.73, p 0.45; eccentricity 0.45
vs 0.69, p 0.04; budding 0.18 vs 0.24, p 0.20), and the W65 minimal-medium chips are not one population: area
varies by 71 % between chips and eccentricity by 42 %, against 14 to 19 % in every other group. W65 cells are
also three to five times larger than W109 cells in both media, so the two families cannot be pooled. **(open)**:
the five W65 omlp chips need a per-chip look (`static/13_endpoint_per_chip.csv`) before this section is written.

- Main text: `static/13_endpoint_vs_medium_area.pdf`, `static/13_endpoint_vs_medium_eccentricity.pdf`,
  `static/21_budding_rate_vs_medium.pdf` (mean ± SEM over chips, one panel per chip family).
- Appendix: `static/10_cell_area_over_time.pdf` (the ypd chambers saturate, which fixes the endpoint window,
  `static/13_endpoint_saturation_per_chamber.csv`), `static/21_panel_a_violin.pdf`,
  `static/12_area_growth_rate_all.pdf`, `static/40_Rp_*.pdf`.

---

## 4. Section 2: oscillations, the apparent dose responses

The pipeline reads one value per structure and plots it against the period with the controls of the same
structure as markers and the pooled control band behind (`13_endpoint_vs_period_<readout>_<osc_type>.pdf`,
`21_budding_rate_vs_period_<osc_type>.pdf`; lower row: bracket score, 0 = famine control, 1 = feast control,
hollow when the two controls do not separate). Spearman over structures, n = number of periods, is an effect
size, not a test (`13_endpoint_spearman.csv`, `21_budding_rate_spearman.csv`).

What appeared:

| readout | series | Spearman vs period |
| --- | --- | --- |
| buds per mother-hour | WT/Glc | −0.90 (p 0.04) |
| | BSPH/Glc | −0.89 (p 0.02) |
| | BSG/Glc, BSO/Glc | −0.66 |
| | pH series | −0.4 to −0.8 |
| `ratio_OxPro` | BSO/Glc | +0.94 (p 0.005) |
| `ratio_pHluorin` | BSPH/Glc | +0.90 (p 0.04) |
| `ratio_GlyRNA` | BSG/Glc | +0.89 (p 0.02) |
| eccentricity | 7 of 10 series positive |

Eight of the nine budding-rate series with a non-zero rho are negative. Tracking cannot have produced this:
fragmentation rises with the period in only 3 of 10 series (`00_track_fragmentation.csv` per structure), and
where it rises it adds false buds and shortens mothers, which pushes the rate up, not down.

---

## 5. The pivot: the controls carry the trends

Constant-medium chambers cannot respond to the period. On the same structures they trend with it:

| readout | series | Osc vs period | strongest control vs period | Osc minus controls |
| --- | --- | --- | --- | --- |
| buds per mother-hour | BSPH/Glc | −0.89 | −0.94 | −0.31 |
| | BSG/Glc | −0.66 | −0.83 | +0.77 |
| | BSO/Glc | −0.66 | +0.83 | −0.60 |
| | WT/Glc | −0.90 | −0.30 | −0.20 |
| `ratio_OxPro` | BSO/Glc | +0.94 | +0.94 (NegCtrl) | +0.03 |
| `ratio_OxPro` | BSO/pH | +1.00 | +1.00 (NegCtrl) | +0.80 |
| `ratio_pHluorin` | BSPH/Glc | +0.90 | +0.90 (NegCtrl) | +0.30 |
| `ratio_GlyRNA` | BSG/Glc | +0.89 | +0.43 | −0.26 |

- `13_endpoint_control_trend.csv`, `21_budding_rate_control_trend.csv`: per readout and series, the Spearman
  of the oscillation chambers, of PosCtrl, NegCtrl and their mean, of the difference, and a verdict. A
  `period effect` needs all three: the oscillation chambers trend, no control trends the same way, and the
  difference trends. Of 28 endpoint combinations, 12 show no oscillation trend, 12 are shared by a control,
  2 vanish once the controls are subtracted, and 2 survive without recurring in another strain
  (eccentricity BSG/Glc, area BSG/pH with four points).
- `13_endpoint_within_culture.csv`, `21_budding_rate_within_culture.csv`: the change from the shortest to the
  longest period inside one culture, where the pre-culture cannot differ. Budding rate: the oscillation
  chambers fall in 9 of 19 cultures, the controls in 6 of 18, the difference in 9 of 19.
- `13_endpoint_bracket_score.csv`, `21_budding_rate_bracket_score.csv`: the feast and famine controls barely
  separate. For the budding rate the bracket is degenerate on 42 of 48 structures, and PosCtrl exceeds NegCtrl
  on 26 of 48 (Wilcoxon p 0.53). Cells under constant famine bud in the sparse window as often as under
  constant feast, 0.18 to 0.77 per mother-hour.
- `40_control_consistency_*.pdf`, `40_control_consistency_*_kruskal.csv`, `11_control_consistency_mu_kruskal.csv`:
  the same controls compared across the structures of a series; they differ.

The mechanism: the structure sets the readout for all its chambers, and the shortest period usually sat on the
same structure position, so position and period were confounded within every culture. The day blocks add a
culture effect on top, balanced across strains for Glc, aligned with the period for pH (`00_chip_run_order.csv`).
Without the controls every drift in section 4 would have been reported as a dose response.

---

## 6. Section 3: biosensors

The sensors report differences between structures, and their famine controls report the same differences.
`ratio_OxPro` and `ratio_pHluorin` rise with the period in the NegCtrl chambers exactly as in the oscillation
chambers; `ratio_GlyRNA` rises in the oscillation chambers but not after its controls are subtracted;
`ratio_Queen-2m` shows no trend. Nothing the sensors show is attributable to the period.

- Main text: `13_endpoint_vs_period_ratio_<sensor>_Glc.pdf` (and `_pH.pdf`), the ratio rows of
  `13_endpoint_control_trend.csv`.
- Whether a sensor responds to feast versus famine at all, independent of the period:
  `95_<strain>_<osc_type>_ratio_<sensor>_comparison.pdf` (PosCtrl against NegCtrl per structure),
  `95_..._control_chambers.pdf` (every control chamber within every structure), `95_..._timeseries.pdf`,
  `95_..._raw_channels.pdf`, with `95_..._summary_per_replicate.csv`. **(open)**: not yet reviewed on the real data.
- Appendix: `30_<channel>_over_time_<osc_type>.pdf`, `31_ratio_<sensor>_over_time_<osc_type>.pdf` (drift over
  hours; cycles are below the sampling limit), `40_Rp_ratio_*.pdf`, `40_Rt_population_ratio_*.pdf`.

---

## 7. PKO

One structure at period 3, built for the clogging hypothesis: without pullulan the control chambers of a
structure should agree better than the producers' do. They do not (`pko/60_pko_within_chip_agreement.csv`,
`pko/61_pko_control_agreement.pdf`): the PKO chamber-to-chamber CV of the area level is 0.24 for famine and
0.36 for feast, the 58th and 88th percentile of all producer structures; the detrended temporal CV sits at the
81st and 43rd percentile. One observation, pointing the wrong way for the hypothesis. PKO feast cells are 2.7
times smaller than PKO famine cells, where producers show a ratio of 1.0 (`pko/60_pko_control_bracket.csv`);
unexplained at n = 1.

---

## 8. Figure plan

| thesis section | main-text figures | appendix figures | evidence tables |
| --- | --- | --- | --- |
| Methods: units and sampling | none (a schematic of chip, structures, arrays and chambers is the author's) | `00_n_tracks_overview_summary.pdf` | `00_chip_overview.csv`, `00_chip_run_order.csv` |
| Methods: tracking and sparse window | `20_lineage_window.pdf` | `qc_comparison/70_qc_effect.pdf`, `lineage_validation/lv_02_detection_rate.pdf`, `20_bud_size_at_appearance.pdf` | `00_track_fragmentation.csv`, `00_track_relinks_per_chamber.csv`, `70_qc_effect_summary.csv`, `lv_03_detection_rate_kruskal.csv` |
| 1 Static medium | `static/13_endpoint_vs_medium_area.pdf`, `static/21_budding_rate_vs_medium.pdf` | `static/13_endpoint_vs_medium_eccentricity.pdf`, `static/10_cell_area_over_time.pdf`, `static/21_panel_a_violin.pdf` | `static/13_endpoint_summary.csv`, `static/13_endpoint_per_chip.csv`, `static/21_budding_rate_summary.csv` |
| 2 Oscillations | `21_budding_rate_vs_period_Glc.pdf`, `13_endpoint_vs_period_area_Glc.pdf` | the `_pH.pdf` counterparts, `13_endpoint_vs_period_eccentricity_*.pdf`, `10_cell_area_over_time_*.pdf`, `12_area_growth_rate_all.pdf`, `40_Rp_*.pdf` | `13_endpoint_spearman.csv`, `21_budding_rate_spearman.csv` |
| 3 Biosensors | `13_endpoint_vs_period_ratio_OxPro_Glc.pdf`, `13_endpoint_vs_period_ratio_pHluorin_Glc.pdf` | `95_*_comparison.pdf`, `31_ratio_*_over_time_*.pdf` | ratio rows of `13_endpoint_control_trend.csv`, `95_*_summary_per_replicate.csv` |
| 4 Controls (the pivot) | proposed summary figure (section 9); `13_endpoint_vs_period_ratio_OxPro_Glc.pdf` as the worked example | `40_control_consistency_*.pdf`, bracket rows of the `13_`/`21_` figures | `13_endpoint_control_trend.csv`, `21_budding_rate_control_trend.csv`, `*_within_culture.csv`, `*_bracket_score.csv`, `00_chip_run_order.csv` |
| PKO | `pko/61_pko_control_agreement.pdf` | `pko/13_endpoint_vs_period_area_Glc.pdf` | `pko/60_pko_within_chip_agreement.csv`, `pko/60_pko_control_bracket.csv` |

Suggested order of the results chapter: static first (the clean result), then the tracking limit as a short
methods-in-results block, then oscillations and biosensors together with the control pivot, then PKO. Section 4
is not a separate finding after the others; it is what the others turn into.

---

## 9. Open points and proposed additions

Questions for the author:

1. **Language and chapter structure.** Is the thesis written in German, and are the four sections the results
   chapters as planned, or does the control pivot become part of the oscillation chapter?
2. **W65 minimal medium.** Which of the five omlp chips carry the 71 % spread? Aggregates, a failed run, or real?
3. **Sensor dynamic range.** Do the `95_*` figures show a feast/famine difference for each sensor on the real
   data? Section 3 needs that sentence before the structure argument.
4. **Structure position.** "Most of the time" the shortest period sat on the same structure: is there a run sheet
   that says which structure each period had? A column `structure_position` per structure would let the
   pipeline test position directly instead of inferring it.

Proposed additions to the pipeline, none of them started:

- A **control-trend summary figure**: one point per readout and series, oscillation Spearman on x, strongest
  control Spearman on y; points on the diagonal are structure effects. This is the single figure for section 4.
- An **events-against-density figure** over the full run for the QC batch, to show that raw budding events follow
  the object count; it costs one full-run lineage pass and would replace the numbers quoted in 2.2.
- **Per-chip points** on the static figures, so the W65 heterogeneity is visible in the figure itself.
- The **control-trend check for µ_area** (`12_area_growth_rate_per_chip.csv`), and the controls back into
  `12_area_growth_rate_all.pdf`, which currently shows oscillation conditions only.

---

## 10. Reproduction

```
cd analyse_pipeline
python run_analysis.py            # all branches: analysis_output/, static/, pko/, no_qc/, qc_comparison/
python validate_lineage.py        # lineage_validation/ (needs the full run first)
```

`config.py` holds every threshold named above (`LINEAGE_SPARSE_*`, `RELINK_*`, `AREA_GROWTH_MIN_FRAMES`,
`BUD_MAX_AREA_FRACTION_FALLBACK`) and logs them, together with the method caveats, at the start of every run.
`README.md` describes each module and output prefix.
