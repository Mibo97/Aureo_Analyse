# Review of the results chapter (PDF "3_results", 20 pages, received 2026-10-05)

Checklist for editing the Word file. The reference text for every passage that has to be updated is
`docs/results_chapter.md` (draft state of 2026-10-05, half-cycle wording, nine robustness metrics). Figure and
table numbers below are those of the PDF; "draft" means the numbering of `docs/results_chapter.md`.

## 0. Answers of 2026-10-05 and what followed

| question | answer | done |
| --- | --- | --- |
| 1 Figure 5 | A: chip W65, complex medium; B: chip W65, minimal medium, the chamber of 3.3, frame 90; scale bars 10 µm | caption in the draft; 3.2.2 keeps pseudohyphae (W65, Figure 5 A) and filamentous growth (true hyphae, only W109 in complex medium) apart |
| 2 Figure 6 stars | pipeline figure `21_panel_a_violin.pdf`: the stars were a Mann-Whitney U over object-frames (every frame of every cell), so every pair was "****" by n alone | pipeline changed: one value per cell, test over structures (median per structure, ≥ 3 per group), stated on the figure; caption in the draft |
| 3 Figure 13 | merged, phase contrast, signal 1 (uvGFP for BSA/BSpH, YFP for OxPro, mTurquoise for Gly-RNA), signal 2 (GFP for BSA/BSpH, mCherry for BSO/BSG); time point not needed, the figure shows that both channels are present in every strain | caption in the draft; a scale bar is still missing |
| 4 Figure 15 | keep the combined figure | pipeline writes `50_control_trend_summary.pdf` again (growth readouts first, then sensors); the split pair stays for the appendix |
| 5 Figure 9 | drop the static panel | `validate_lineage.py` plots the oscillation series only |
| 6 Table 2 | drop the caption sentences | draft |
| 7 Rensink | both papers exist | keys kept as they are; check each citation against the right paper |
| 8 area axes | µm² | every area figure in µm², tables in px² plus `*_um2` columns in `13_endpoint_*.csv` |
| 9 sensor robustness | drop the sentence in 3.6 | draft; the 3.5 paragraph and Table 9 marked optional |

The draft `docs/results_chapter.md` now carries the thesis numbering (Figures 2 to 17, Tables 1 to 6) and
captions for Figures 3, 5, 6, 8, 9, 12 and 13; the mapping table in section 2 below is therefore only needed for
the Word file.

## 1. Coherence in short

The argument of the chapter is intact and in the right order. What is not coherent yet is the apparatus around
it: two Word cross-references are broken, the subsection numbers of 3.2 repeat, almost every figure and table
reference in the running text still carries the draft numbering, four paragraphs describe figures or tables of
an earlier draft (bands, bracket rows, lines, seven instead of nine robustness metrics), one sentence is
contradicted by its own figure (µ_bud against µ_area), and the appendix holds one of the twelve figures the
text refers to. Seven of the new figures are never cited in the text.

## 2. Broken references and numbering

- "Fehler! Verweisquelle konnte nicht gefunden werden." appears twice: p. 10 (should be Table 1) and p. 12
  (should be Table 2). Re-insert the two cross-references and update all fields (Ctrl+A, F9).
- Subsections: the second "3.2.3 Budding events, the size criterion and validation" and the second "3.2.4
  Dependence of the budding trends on the tracking" must be 3.2.5 and 3.2.6. The text already refers to
  "the size criterion of 3.2.5" (3.2.2) and "changed with the tracker (3.2.6)" (3.2.4).
- Figure references in the text still use the draft numbering. Mapping draft → PDF:

  | where in the text | reads | should read |
  | --- | --- | --- |
  | 3.2.2, twice ("Figure 3 shows every tracked cell", "second cluster in Figure 3") | Figure 3 | Figure 4 |
  | 3.2.4 | Figure 4 A, Figure 4 B | Figure 7 A, Figure 7 B |
  | 3.2.5 "(valley depth 0.25; Appendix Figure A2)" | Appendix Figure A2 | Figure 8 |
  | 3.3, end of the W109 paragraph | Figure 4 | Figure 10 |
  | 3.4, second paragraph | Figure 5, Appendix Figures A5 to A7 | Figure 11, Appendix … |
  | 3.5, first sentence | Figure 6, Appendix Figure A8 | Figure 14, Appendix Figure 17 (or A-number) |
  | 3.6 "Figure 7 shows all 58 combinations" | Figure 7 | Figure 15 |
  | 3.7 "(Figure 8)" | Figure 8 | Figure 16 |

- Table references: "Table 6, columns 'oscillation chambers'" (3.4) → Table 5; "(Table 8; Appendix Figure
  A11)" (3.4) → Table 4; "in 15 of them with a |ρ| of 0.8 or more (Table 6)" (3.6) → Table 5; "(Tables 8 and
  9)" (3.6) → Table 4 and the sensor-robustness table, which does not exist yet (section 3 below).
- Figures 3, 5, 6, 8, 9, 12 and 13 are never cited. Every figure needs one sentence that sends the reader to
  it (suggestions in section 4).
- Appendix: the text refers to A2, A4, A5, A7, A8 and A11; the appendix contains one figure, numbered
  "Figure 17". Give the appendix its own caption label (Figure A1, A2, …) so that "Appendix Figure A8" can
  be a cross-reference, and add the figures listed in section 5.

## 3. Passages that predate the current draft

Take the wording from `docs/results_chapter.md`; the numbers there are reproduced by the pipeline.

- **Half-cycle period.** "cycle period" still stands in the chapter intro, 3.2.1, the Table 2 caption and
  header ("periods (min)"), the headings of 3.4 and 3.5, the Table 5 caption and the captions of Figures 9,
  11 and 14. The defining sentence of draft 3.2.1 is missing ("The interval is the half-cycle period: the
  medium was switched every 0.75 to 24 min, so a full feast/famine cycle lasted 1.5 to 48 min"), and so is
  the sampling sentence that opens draft 3.4 (intervals of 0.75 to 6 min unresolved at 10 min per frame,
  12 min at the limit, 24 min resolved with 4.8 frames per cycle; all readouts cumulative over hours).
- **3.4, second paragraph** describes the old figure: "the pooled control range as a band" and "In the bottom
  row of every panel the oscillation chambers are expressed relative to the two controls (bracket score …)".
  Figure 11 has neither. Replace by the draft paragraph (bracket score as a table quantity, undefined where
  the controls do not separate).
- **Figure 11 caption**: "dashed and dotted lines: pooled control means" (no lines in the figure) and "mean ±
  SEM over chambers" (the figure legend says mean ± SD of the structure's chambers, which is what the
  pipeline draws; the draft caption said SEM too and is corrected). Same SD wording for Figure 14.
- **Robustness paragraph in 3.4 and Table 4** are the seven-metric version (70 combinations: 39/18/7/6, three
  of six period effects in BSA/Glc). The current draft has nine metrics, with the two budding metrics you
  asked for: R(p) of the budding rate per mother (−2.018 oscillation, −1.984 feast, −1.921 famine; below the
  control mean on 25 of 49, p 0.95; 1 period effect) and single-cell R(t) of µ_event (−1.115, −1.132, −0.999;
  29 of 47, p 0.022, Glc p 0.043, pH p 0.28; 1 period effect). The sentence on the period then reads 90
  combinations: 51 no trend, 23 structure effects, 8 not robust, 8 period effects, five of them in BSA/Glc.
  Table 4 needs the two extra rows and the caption addition on µ_event (mothers with at least three
  intervals) and on R(p) of the budding rate.
- **3.5 lacks the "Robustness of the sensor signals" paragraph and Table 9** (sensor robustness, 12 rows),
  although 3.6 refers to "Tables 8 and 9" and reports "94 combinations … 48/30/9/7". With the nine metrics
  the draft says 114 combinations: 60 no trend, 35 structure effects, 10 not robust, 9 period effects.
  Either add the paragraph and the table (recommended; robustness is the question of the thesis) or drop
  the sentence and the reference to Table 9.
- **Table 2 caption** promises "Cell tracks are counted after the cell filter" and "Budding events are the
  accepted mother-bud assignments", but the table has no such columns. Either add the two columns (I can
  produce them per row from `00_cell_filter.csv` and `20_budding_events.csv`) or delete the two sentences.
- **Table 3 (static)**: the p column is incomplete. Welch's t-test over the chambers gives W109 µ_bud
  p 0.099, W65 area p 0.89, W65 buds per mother-hour p 0.66, W65 µ_bud p 0.68; replace "n.s." by the values
  (updated in the draft).
- **µ_bud against µ_area (3.4, before Figure 12).** "µ_bud and µ_area were negatively correlated (ρ −0.4,
  n = 146): where cells grew in area, they budded less" is contradicted by Figure 12 itself, which
  annotates ρ −0.24 for Glc and +0.39 for pH. The pooled value exists only because the Glc structures
  combine a higher µ_bud with a lower µ_area than the pH structures (medians 0.25 against 0.16 births per
  cell-hour, 0.11 against 0.18 h⁻¹). Corrected sentence (now in the draft): "Per structure and chamber type,
  µ_bud and µ_area were not coupled in one direction (Figure 12): pooled over both oscillation types they
  correlated negatively (Spearman ρ −0.41, n = 146), but only because the Glc structures combined a higher
  µ_bud with a lower µ_area than the pH structures (…); within the Glc series the correlation was −0.24
  (n = 87) and within the pH series +0.39 (n = 59)." The same correction is in `docs/data_story.md`.
- **Citation keys**: "(2024_Rensink)", "(2026_Rensink)" (three times in 3.2.2 and in the Figure 4 caption)
  and "(2022_Trivellin; 2024_Blöbaum)" (Table 4 caption) are placeholders from the draft. The Rensink paper
  is cited as 2024 in 3.2.1 and 3.3 and as 2026 elsewhere; settle on one year.

## 4. Figures: content and captions

The pipeline figures in the PDF come from different pipeline states (Figure 4 from the per-cell morphology
version, Figure 11 without lines, Figure 15 still the combined summary, x labels "cycle period"). After the
final run replace all pipeline figures in one go; their x labels then read "half-cycle period [min]".

- **Figure 2** (BioLector): fine. Open: OD₆₀₀ calibration (needs your end-point OD readings).
- **Figure 3** (tracks per strain and period, `00_n_tracks_overview_summary.pdf`): not cited. Caption:
  "Number of cell tracks per chamber after the cell filter, mean ± SD over the chambers of each structure,
  against the half-cycle period; bars: oscillation chambers, feast control, famine control; top row Glc,
  bottom row pH series." Sentence for 3.2.1, after "63,215 cell tracks … remained": "The number of cell
  tracks per chamber ranged from 30 to 236 (10th to 90th percentile over the 536 chambers of the oscillation
  experiments, median 83); feast-control chambers held more tracks (median 121) than oscillation (71) and
  famine-control chambers (69), the WT/pH structures the most (up to 333 per chamber on average) and
  BSG/pH/6 min the fewest (24) (Figure 3)." (Numbers from `00_data_overview.csv` of the v12 run.)
- **Figure 4** (morphology scatter): caption fine. The text says "for the WT strain of the oscillation
  experiments"; the figure shows all five strains. The static panels (`static/90_morphology_scatter.pdf`,
  bottom row in the draft) are missing although the text quotes their shares (82 % on W109 complex, 25 % of
  51 cells on W65); add them or write "not shown".
- **Figure 5** (phase-contrast examples): good addition. The caption needs chip, medium, strain and frame
  for A and B, "pseudohyphae" (not "pseudohyphea"), "red outlines: segmentation masks" (not "segmentate
  ROIs") and the scale-bar length. If B is the W65 minimal-medium chamber of 3.3 (two cells of about 270 µm²
  and the ring of ten blastoconidia, frames 85 to 88), cite "Figure 5 B" there. If A is from W109 in complex
  medium, cite it in the last sentence of 3.2.2.
- **Morphology text**: the chapter now says pseudohyphae were observed (first sentence of 3.2.2, last
  sentence "Filamentous growth was only observed in the W109 chip in complex medium"), the draft said "no
  filamentous growth" from the mask statistics. Both hold: a pseudohyphal chain is segmented into its single
  elongated cells (Figure 5 A shows four masks). Suggested wording for the last sentence: "Pseudohyphal
  chains of elongated cells (Figure 5 A) were seen only on chip W109 in complex medium; the segmentation
  split them into their individual cells (median solidity 0.98, axis ratio above 3 in fewer than 0.3 % of
  the objects of the four structures with complete shape data), so they enter Figure 4 as single elongated
  cells, not as filaments."
- **Figure 6** (violins with stars): three problems. It is not cited and sits inside 3.2.4 although it
  belongs to 3.2.2 next to Figure 4. "Stars mark significancy" names no test and no unit; if the test ran
  over cells (thousands per strain) every comparison is significant and the stars say nothing about the
  strains, which are one culture per structure. The area axis is in px² while the text is in µm². The
  per-strain medians are already in the text and the distributions in Figure 4. Recommendation: drop the
  figure, or keep it without stars (or with a test over structures, n = 4 to 6 per strain, named in the
  caption) and move it to 3.2.2.
- **Figure 7** (density): caption fine; "chambers(B)" → "chambers. (B)".
- **Figure 8** (bud size, `20_bud_size_at_appearance.pdf`): caption: "Size criterion of the mother-bud
  assignment. Area of a new object at first detection relative to the area of its assigned mother, all
  13,940 candidates. Left: all branches pooled, histogram and density on a log scale, modes at 0.07 (buds)
  and 0.45 (washed-in cells and split masks), antimode 0.32 (dashed) with the search range (grey). Right:
  the density per strain and oscillation type (n = candidates)." Cite it in 3.2.5 instead of A2.
- **Figure 9** (assignment rate, `lineage_validation/lv_02_detection_rate.pdf`): caption: "Assignment rate
  of the bud candidates. Fraction of the bud candidates of a chamber that received a mother, one point per
  chamber, colour = strain, bar = median per period; Glc series (left), static chips (middle), pH series
  (right)." The static panel labels (static_W65 / static_omlp / static_ypd) mix chip and medium; the text
  does not discuss the static assignment rate, so drop the panel, or I relabel the groups (W65 both media,
  W109 minimal, W109 complex) in `validate_lineage.py`.
- **Figure 10** (static): caption matches (A area, B µ_bud). The draft had buds per mother-hour as a third
  panel (`static/21_budding_rate_vs_medium.pdf`); the text discusses that readout at length, so consider
  adding it.
- **Figure 11**: caption as in section 3 (no lines, SD, half-cycle). The budding rate per mother-hour
  (`21_budding_rate_vs_period_Glc.pdf`), on which the paired comparison rests (27 of 29 structures), is in
  neither the text nor the appendix; add it to one of them.
- **Figure 12** (`24_mu_bud_vs_mu_area.pdf`): caption: "Population growth from budding against single-cell
  area growth. µ_area (mean over the tracks of the chambers) against µ_bud, one point per structure and
  chamber type (circles: oscillation chambers; triangles: feast control filled, famine control hollow),
  colour = strain; left Glc, right pH series; dashed line: equal rates; ρ: Spearman over the points of the
  panel." Text: "In Figure 12 the both growth rates are compared" → "Figure 12 compares the two growth
  rates."
- **Figure 13** (microscopy of the biosensor strains): caption needs the panel order and the channel of each
  image (merged, phase contrast, reference channel, sensor channel, with the fluorophore or filter per
  strain), the condition (chamber, time point) and a scale bar. "Microscopy-footage" → "Microscopy images".
  Cite it where the ratio is defined (first sentence of 3.5): already there, keep.
- **Figure 14**: caption: "(A) OxPro ratio of BSO and (B) sfpHluorin ratio of BSpH at the endpoint against
  the half-cycle period, Glc series; markers as in Figure 11." ("BSPH" → BSpH.)
- **Figure 15**: the combined summary of the earlier run. The final pipeline writes two figures
  (`50_control_trend_summary_growth.pdf`, `_sensors.pdf`, draft Figure 8 A/B). Keeping the combined one is
  fine (same data), but it cannot be regenerated by the final run. Caption fine.
- **Figure 16** (PKO): the figure has panels a, b and c; the caption covers a and b. Add: "(c) Cliff's δ of
  the feast against the famine controls for µ_area (does the control bracket separate at all); the knockout
  structure at the producers' median."
- **Figure 17** (appendix, `95_*_comparison.pdf`): caption: "Sensor ratios of the control chambers. Ratio
  per chamber (median over the time points after the 2 h preconditioning), famine control against feast
  control, one point per chamber, bar = median and interquartile range over chambers; top row Glc series
  (0 against 50 g/L glucose), bottom row pH series (pH 8 against pH 3); one panel per sensor."

## 5. What is missing in the appendix

Referenced in the text or needed for a claim; pipeline file names in brackets, "next run" = produced by the
final pipeline run.

| appendix | content | files | state |
| --- | --- | --- | --- |
| A4 | static: objects per frame over time (why the window ends at frame 77); endpoint eccentricity per chip | `static/10_cell_area_over_time.pdf`, `static/13_endpoint_vs_medium_eccentricity.pdf` | exist |
| A5 | pH counterparts of Figure 11, and the budding rate per mother-hour for both series | `13_endpoint_vs_period_area_pH.pdf`, `24_growth_from_budding_vs_period_pH.pdf`, `21_budding_rate_vs_period_Glc.pdf`, `_pH.pdf` | exist; relabelled x axis with the next run |
| A6 | immigration against the period | `24_immigration_vs_period_Glc.pdf`, `_pH.pdf` | exist |
| A7 | eccentricity against the period, µ_area overview | `13_endpoint_vs_period_eccentricity_Glc.pdf`, `_pH.pdf`, `12_area_growth_rate_all.pdf` | exist |
| A8 | sensor controls (= Figure 17) plus the pH counterparts of Figure 14 and the QUEEN-2m and Gly-RNA endpoint figures | `13_endpoint_vs_period_ratio_<sensor>_<Glc,pH>.pdf` | exist |
| A9 | control consistency across structures (Kruskal-Wallis paragraph of 3.6) | `40_control_consistency_*.pdf` | exist |
| A10 | PKO endpoint area against the producers (optional) | `pko/13_endpoint_vs_period_area_Glc.pdf` | exists |
| A11 | robustness metrics against the period with the controls | `40_<metric>_<readout>_vs_period_<osc_type>.pdf` | next run |
| A12 | robustness control-trend summary | `50_robustness_control_trend_summary.pdf` | next run |
| A13 | paired comparison oscillation against controls | `51_osc_vs_controls.pdf`, `_robustness.pdf`, `_sensors.pdf` | next run |
| Table | sensor robustness (draft Table 9) if it does not go into 3.5 | `40_*_ratio_*_per_chip.csv`, `51_osc_vs_controls.csv` | next run |

## 6. Small things

- Typos and spacing: "Table 1:Growth"; "mediansover"; "0.6although"; "thesparse-phase"; "combinations**"
  (stray asterisks, Table 6 caption); "significancy"; "pseudohyphea" (twice); "segmentate ROIs";
  "Diskussion"; "diƯerence" (twice on p. 26, a broken ff ligature; check the font or the PDF export);
  "chambers(B)" (Figure 7 caption); "Microscopy-footage".
- Table 2 header "periods (min)" → "half-cycle periods (min)".
- Justified paragraph on p. 22 ("marginally rounder (eccentricity ratio 0.98, p < 0.001) and grew more
  slowly") is spaced out; a non-breaking space in "p < 0.001" fixes it.
- Area axes of Figures 6, 10 and 11 are in px² while the text is in µm². The pipeline can plot in µm²
  (UM_PER_PX = 0.0733 is in the config); say if you want that.
- Strain name: "BSPH" in pipeline titles and the Figure 14 caption, "BSpH" in the text.

## 7. Questions

1. Figure 5: which chip, medium, strain and frame are A and B? Is B the W65 chamber of 3.3?
2. Figure 6: which test and which unit produced the stars? Keep the figure at all?
3. Figure 13: the channel of each panel per strain, and which chamber and time point?
4. Figure 15: keep the combined figure or switch to the two figures of the final run?
5. Figure 9: drop the static panel, or should I relabel its groups?
6. Table 2: add the "cell tracks" and "budding events" columns or delete the two caption sentences?
7. Rensink: 2024 or 2026?
8. Area axes in µm² in the pipeline figures?
9. Sensor robustness: add the paragraph and Table 9 to 3.5, or drop the 114-combinations sentence in 3.6?
