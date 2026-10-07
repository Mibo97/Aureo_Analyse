# Literature passages for the discussion (4.1 and 4.4 on 2026-10-06; 4.5 to 4.7 added 2026-10-07)

Extracted 2026-10-06 from the twelve PDFs you sent (`literatur-4.1.zip`, `literatur-4.4.zip`). Page numbers are
PDF pages (page 1 = first page of the file), not the printed journal pages. Quotes are verbatim from the text
layer, shortened with "…". Each passage names the claim of `docs/discussion_outline.md` it supports and how to
use it. Where a paper contradicts or qualifies a claim, that is said.

## What each paper is good for

| key | paper | use in the discussion |
| --- | --- | --- |
| 2022_Trivellin | Quantification of microbial robustness in yeast (ACS Synth Biol) | definition of R, normalisation by m, robustness is function-specific, R is relative to the strains tested (4.1, 4.3) |
| 2022_Olsson | Robustness: linking strain design to viable bioprocesses (Trends Biotechnol) | robustness against tolerance, phenotype stability under perturbations, heterogeneity from gradients (4.1, 4.3, 4.7) |
| 2024_Blöbaum | Quantifying microbial robustness in dynamic environments using dMSCC (Microb Cell Fact) | the yeast chip design, replicate = chamber, R(t) and R(p) definitions and limits, the yeast responses to the same intervals (4.1, 4.2, 4.3) |
| 2023_Blöbaum | Protocol dMSCC of C. glutamicum (STAR Protocols) | control zones and the switching zone, chamber exchange to 95 % (4.1, 4.4, 4.8) |
| 2020_Täuber | dMSCC platform (bioRxiv) | exchange time with fluorescein (95 % after 5 s), C. glutamicum growth and cell length against the interval, time scales of reactors (4.1, 4.2, 4.4, 4.7) |
| 2010_Lazic | Pseudoreplication (BMC Neurosci) | technical against biological replicates, experimental unit, inflated false positives, remedies (4.1) |
| 2025_Pianale | Physiology and robustness of yeasts in dynamic pH and glucose environments (Biotechnol Bioeng) | the closest comparison: substrate against pH oscillations in three yeasts at 0.75 to 48 min; sensor response times (4.2, 4.3, 4.4) |
| 2022_Pianale | Biosensor toolbox (Front Microbiol) | the four sensors, their channels and normalisation, calibration by permeabilised cells, no growth burden in S. cerevisiae (4.4, 4.7) |
| 2014_Yaginuma | QUEEN ATP indicator (Sci Rep) | principle, Kd, dynamic range, pH and temperature dependence of the QUEEN signal (4.4) |
| 1998_Miesenböck | pHluorins (Nature) | ratiometric pHluorin principle and calibration in buffers (4.4) |
| 2024_Haala | DoE medium optimisation for liamocins in A. pullulans (Front Bioeng Biotechnol) | oMLP origin, same strain NRRL 62042, storage lipids under high C/N, growth rates (4.4, 4.5) |
| 2024_Zhang | UV mutagenesis and fluorescence screening for pullulan (Fermentation) | not a burden paper (see end); usable for morphotypes and melanin in 4.5 |

---

## 4.1 The period did not matter; the structure did

### Claim: robustness means phenotype stability across a perturbation space, distinct from tolerance

- 2022_Trivellin p.1: "Robustness is defined for a specific function (or phenotype) and set of perturbations.
  Robustness is therefore different from tolerance, which specifically describes stable growth or survival to
  various perturbations". Use: the opening definition; our perturbation space is the set of switching
  intervals.
- 2022_Trivellin p.1 (abstract): "The methodology validated that robustness is function-specific and
  characterized by positive and negative function-specific trade-offs." Use: why size, growth rate and
  budding can answer differently (4.3).
- 2022_Olsson p.1: "microbial robustness … describes the stability of a phenotype (e.g., titer, production rate,
  and yield) when challenged by different disturbances or perturbations". p.3: "although tolerance only
  considers viability and/or specific growth rate, microbial robustness encompasses the stability of the
  production and growth of a microbe in different conditions evaluated using different measures." Use: the
  BioLector growth (3.1) is tolerance, the chip data are robustness.
- 2022_Olsson p.1: "Single-cell analysis can elucidate the roles of subpopulations and population dynamics in
  microbial robustness. Quantification of robustness is suggested as a tool for guiding strain and
  bioprocess development." Use: why single-cell data and R(p) belong in the answer.

### Claim: what the design can and cannot separate; one culture per structure, technical replicates

- 2010_Lazic p.1: "Pseudoreplication occurs when observations are not statistically independent, but treated
  as if they are. This can occur when there are multiple observations on the same subjects, when samples are
  nested or hierarchically organised, or when measurements are correlated in time or space."
- 2010_Lazic p.1, quoting Hurlbert: "use of inferential statistics to test for treatment effects with data
  from experiments where either treatments are not replicated (though samples maybe) or replicates are not
  statistically independent".
- 2010_Lazic p.2: "An experimental unit is defined as the smallest entity that can be randomly assigned to a
  different treatment condition". Use: the structure is the experimental unit of the period; the five
  chambers are not.
- 2010_Lazic p.2: "The difference here is between biological replicates which are independent … and technical
  replicates which are not independent (e.g. dividing a blood sample from a single rat into two sub-samples
  …)". Use: chambers of one structure are the sub-samples.
- 2010_Lazic p.2: "measuring the height of one person ten times provides only one independent piece of
  information about the population." Use: five chambers, one piece of information about the period.
- 2010_Lazic p.4: "When there is a positive within group correlation (the more common situation), the Type I
  error probability (α) will be greater than 0.05, and the greater the correlation the greater the number of
  false positives." Use: why chamber-level tests against the period would have been wrong.
- 2010_Lazic p.4, Table 1: hierarchical or nested data, solutions "Include random effects" or "Average over
  observations"; p.11: "In the first case, average all the x and y observations on each subject, and then use
  these for analyses. In the second case, analysis of the slopes or mixed models can be used." Use: our
  structure means and the Spearman over structures are the "average over observations" remedy; mixed models
  are the outlook.

### Claim: the dMSCC designs in the literature use the same unit; replicate = chamber

- 2024_Blöbaum p.3: "The microfluidic structure … consists of six connected cultivation structures, that
  include one dedicated positive control and five oscillation-structures. Within each structure, six arrays
  of 23 monolayer-growth chambers were exposed to the dynamic condition, with one additional array acting as
  negative control". p.6: "Further development of the dMSCC chip allowed the simultaneous assessment of five
  independent media oscillation frequency … rather than the previous three". Use: one frequency per structure
  is the published design; our structures carry one interval for the same reason.
- 2024_Blöbaum p.6: "Error bars generally represented the standard deviation for the average robustness or
  performance of each replicate chamber." p.17 (supplement captions): "Dispersion of the data corresponds
  to the standard deviation across triplicates (three chambers)." Use: in the published yeast work, as here,
  the replicate is the chamber; the constant-medium chambers on every structure are what let us see what this
  unit hides. Say it as a shared limitation, not as a fault.
- 2023_Blöbaum p.10: "The chip can therefore be divided into control zones, which are always exposed to the
  same liquid and a switching zone, across which the laminar boundary layer moves when changing the flow
  profile." p.9: "The laminar boundary layer can be shifted between each array pair … This allows to adjust,
  for example, the size of the control zones." Use: where the feast and famine controls sit and why they are
  at the ends of the array; the position confound of 4.2 follows from the laminar layout.
- 2020_Täuber p.15: control experiments for the switching itself: "a switch between BHI and BHI was performed in
  10 second oscillation intervals, which resulted in a colony growth rate of µcolony = (0.92 ± 0.02) h⁻¹.
  These findings show that pressure variation which might result from the medium oscillations have no
  significant effect on the colony growth rate." Use: the pressure switching as such is not the stressor; a
  control we did not run and can cite instead.

### Claim: the medium reaches the cells within seconds, so a 45-s half-cycle is transmitted (also 4.4, 4.8)

- 2020_Täuber p.13: "95% is reached in the cultivation chamber after 5 seconds." and "medium switches up to
  5 seconds can be performed, guaranteeing an almost complete medium exchange within the cultivation chamber.
  At oscillation intervals of 2 seconds (0.5 Hz) the fluorescence signal is 82% so that a full exchange of
  nutrients cannot be guaranteed anymore". Also p.13: "Depending on the position of the cultivation chamber
  and thus the flow regime during the switching zone, the exchange rates of fluorescein can slightly differ."
  Use: with 45 s as the shortest half-cycle, the chambers see the full amplitude; the position remark is the
  literature hook for the position confound. Caveat to write: measured in the C. glutamicum chip (chamber
  height below 1 µm), not in the yeast-type chip used here, and under their flow; hence the fluorescein test
  in 4.8.
- 2023_Blöbaum p.9: "With this procedure, you can validate, that the medium inside the chamber is exchanged
  to 95% even during oscillation." and the note that the measurement "is only necessary when your experiments
  include very fast oscillations (<10 s) between media conditions or if you test new chamber heights." Use:
  the authors themselves treat intervals above 10 s as safely transmitted; our chamber height differs, which
  is the one reason to still recommend the test.

### Claim: R is a relative measure; its level depends on the data set

- 2022_Trivellin p.2: "For each function i, a strain S, and a perturbation space P, R was calculated as
  … −(Fano factor)·(1/m) … To allow comparison of R between functions, we normalized the different Fano
  factors to the mean of the functions they describe (m) across all strains … We set the upper limit for R to
  0 (highest robustness)". p.3: "its suitability could be limited by the normalization with m (the mean of the
  function performance among the strains) … it further adjusts R to the strains tested in the study." Use:
  why only comparisons within our data set (oscillation against controls of the same structure) are
  interpretable, not the R levels themselves; also why the unit of the area does not change R.
- 2024_Blöbaum p.8: "According to Eq. 1, R(t) changes upon addition of more replicates or conditions, making
  it a relative and not an absolute term."

---

## 4.2 (passages found while reading for 4.1 and 4.4)

- 2024_Blöbaum p.1 (abstract): "We observed a decrease in specific growth rate but an increase in
  intracellular ATP levels with longer oscillation intervals. Cells subjected to 48 min oscillations
  exhibited the highest average ATP content, but the lowest stability over time and the highest heterogeneity
  within the population." Use: in S. cerevisiae the interval mattered; in A. pullulans it did not. Note that
  their frequencies also sat one per structure with chamber replicates (above), so the comparison is between
  designs of the same kind.
- 2024_Blöbaum p.6: "Even though oscillation frequencies differed, the total time spent under starvation or
  feast conditions was the same for all cells." Use: the time-averaged medium is the same for all intervals,
  so an interval effect cannot be a dose effect; supports reading 2 of 4.2 against the average.
- 2024_Blöbaum p.8 and p.14: "Cell circularity varied substantially with 1.5 and 6 min oscillations, owing to a
  shift towards pseudohyphal growth … Pseudohyphal growth is a known response to nitrogen starvation and
  stress"; "pseudohyphal growth was triggered by glucose shifts every 1.5 and 6 min". Use: oscillation can
  change morphology in yeast; our pseudohyphae appeared in complex medium, not with the interval (4.5).
- 2025_Pianale p.1 (abstract): "All strains showed higher sensitivity to substrate than pH oscillations."
  p.4: "S. cerevisiae is notoriously resistant to low pH …, explaining the absence of growth impairment by any
  pH oscillation … in the pH-dynamic setup, substrate (i.e. glucose) was not a limiting factor as it was
  sufficient to maintain internal homeostasis and stable growth." Use: the same Glc-against-pH asymmetry as
  ours, with the same explanation offered: glucose present throughout the pH cycle.
- 2025_Pianale p.4: "the budding ratio declined the most with oscillations of 1.5 and 6 min in CEN.PK113-7D
  and of 6 and 24 min in industrial strains". p.5: "Overall, substrate oscillations resulted in bigger and
  rounder cells than pH oscillations". Use: interval effects in yeast were strain-specific and
  non-monotone, which a Spearman over four to six periods would not capture either; and the direction of the
  size effect (bigger under substrate oscillation) matches our Glc result.
- 2025_Pianale p.5: "Relative ATP levels dropped significantly (to < 1) when cells were exposed to substrate
  oscillations while remaining similar to the control (~1) in pH oscillations". Use: what a responsive ATP
  sensor shows under glucose oscillation; our QUEEN-2m ratio did not separate even feast from famine (4.4).
- 2020_Täuber p.15 to p.17: C. glutamicum under BHI/PBS oscillation: µ 0.32 to 0.37 h⁻¹ at 1 h to 30 min,
  0.29 to 0.30 h⁻¹ at 5 to 15 min, 0.49 to 0.58 h⁻¹ at 1 min to 10 s, constant BHI 0.90 h⁻¹; p.18: average
  cell length changed with the interval (3.4 µm at 1 h, 2.0 to 2.6 µm at shorter intervals). Use: the bacterial
  precedent for interval-dependent growth and size, against which the flat A. pullulans response stands out.
- 2020_Täuber p.18: "fitness costs resulting from gene regulatory expression may be reduced after the first
  stimulus is applied … so that the cells do not have to constantly adapt to the new conditions". Use: one
  mechanism offered for insensitivity to fast switching.

---

## 4.3 (passages for the robustness metrics)

- 2024_Blöbaum p.5: definitions. R(p): "σ and x refer to the standard deviation and mean, respectively, of a
  function … across all cells at each time point. Instead, m refers to the mean of a function across all time
  points and conditions. Therefore, R(p) describes how homogeneous a function is across a cell population."
  R(t) at population and single-cell level as in our methods 2.5.5.
- 2024_Blöbaum p.16: "R(t) is not able to differentiate between oscillating and steadily changing functions."
  Use: our caveat that R(t) of the oscillation chambers alone is not interpretable at unresolved intervals.
- 2024_Blöbaum p.12: "slower feast-starvation oscillations led to greater heterogeneity in intracellular ATP
  levels … as manifested by a decrease in R(p) with longer oscillations." Use: contrast with our sensors,
  whose R(p) showed no interval dependence beyond the controls.
- 2025_Pianale p.8: "All strains exhibited less population heterogeneity when exposed to substrate than pH
  oscillations." Use: heterogeneity under oscillation is organism- and perturbation-specific; ours was larger
  under glucose oscillation than under constant medium for area, not for budding output.
- 2022_Olsson p.5: "Phenotypic heterogeneity may occur owing to intrinsic factors (e.g., individual
  physiological state, stochastic gene expression/noise, and cell-cycle phases) and/or heterogeneous
  extrinsic factors (e.g., chemical mixing gradients and variation in cell density)". Use: the two sources
  behind R(p); cell density varies with the sparse-phase window.

---

## 4.4 Famine did not starve, and the sensors did not separate

### Claim: the famine medium in the chamber; what the published chips transmit

- 2020_Täuber p.13 and 2023_Blöbaum p.9 as under 4.1 (95 % exchange after 5 s). Use: the published value we
  rely on; the famine state should have reached the cells, so the budding readout's insensitivity to famine is
  biological or a window effect, not a transport effect. Keep the fluorescein test as future work because of
  the different chamber height and the passive 10 mbar supply of the control arrays (methods 2.4.3).

### Claim: reserves and the preculture; blastoconidia formation under limitation

- 2024_Haala p.6: "A. pullulans form intracellular storage lipids under unbalanced cultivation conditions,
  like high carbon (C) and low nitrogen (N) concentrations (high C/N ratio)." p.2: "storage molecules, which
  are present as intracellular lipids, primarily as C16:n or C18:n, and can reach 65% of the cell dry weight
  (Xue et al., 2018)". Use: the preculture in oMLP with 50 g/L glucose and 1.9 g/L NH4NO3 (C/N 125 in the
  production medium, p.10) is a lipid-storing condition; cells enter the chip with reserves that can carry
  22 h of famine. Same strain: p.13 "A. pullulans NRRL 62042".
- 2024_Haala p.13: "The maximal growth rate was increased from 0.11 h⁻¹ to 0.18 h⁻¹ … only a few samples were
  taken during the growth phase, and the growth rate could, therefore, be inaccurate." Use: the bioreactor
  growth rates of the same strain in oMLP are far below our BioLector 0.50 h⁻¹; note the different glucose
  (208 g/L) and sampling before comparing.
- 2025_Pianale p.5: "More than 500 genes have been found to be differentially expressed between starving and
  calorie-restricted (i.e., enough substrate for metabolic activity, but not reproduction) cells (Gulli et al.
  2019)." and the authors' own change from 0 g/L to 10 mg/L glucose as the famine condition. Use: "famine"
  is not one state; our 0 g/L oMLP still contains the other components; a sentence on what the famine medium
  is and is not.
- For conidiation under nutrient limitation in A. pullulans none of the twelve papers has a passage; use
  2024_Rensink / 2026_Rensink (cell-type dynamics) and look for a dedicated source.

### Claim: the sensors are expressed, but their dynamic range in the chip is not established

- 2022_Pianale p.4, Table 1: QUEEN-2m, ATP, λex 410 and 480, λem 520; sfpHluorin, intracellular pH, λex 390
  and 470, λem 512; GlyRNA, glycolytic flux via fructose-bisphosphate, mTurquoise2 with mCherry normaliser;
  OxPro, oxidative stress response via YAP1 activation, ymYPET with mCherry normaliser. p.7: "a constitutively
  expressed fluorescent reporter … was added to normalize the biosensor output … As both QUEEN-2m and
  sfpHluorin are ratiometric probes, they did not require this addition." Use: the exact meaning of the four
  ratios in 3.5 and of the channels in Figure 13.
- 2022_Pianale p.6: intracellular pH calibration: digitonin-permeabilised cells "added to citric acid/Na2HPO4
  buffer, whose pH ranged from 4.5 to 8 … Fluorescence was plotted against pH and calibration curves were
  generated." Use: the calibration that exists for the yeast toolbox and that our methods 2.3 reproduced in
  the BioLector; what is missing is the same curve measured in the chip and a clamp for the other three
  sensors.
- 2025_Pianale p.5: "While the ATP biosensor (QUEEN-2m) had a response time of a few seconds …, the oxidative
  stress and glycolytic flux biosensor (GlyOx) had a response time in the order of tens of minutes (being
  based on transcription, translation, and maturation of fluorescent proteins …). However, these biosensors
  were still useful for grasping general trends … biosensors with response times slower than environmental
  oscillations are still able to spot differences as the biosensor response shifts to a new steady state".
  Use: why OxPro and Gly-RNA cannot follow a 45-s half-cycle and why an endpoint ratio is the right readout;
  also why a flat endpoint ratio can still mean "no shift of the steady state".
- 2025_Pianale p.5: "For all biosensors, relative changes to the positive control were computed … and
  represent the fold increase/decrease for each parameter". Use: the published way to report sensor data is
  relative to a constant-medium control on the same chip, which is what our bracket does; our bracket was
  degenerate.
- 2014_Yaginuma p.2 to p.3: "QUEEN-2m was designed as a low affinity ATP indicator so it is suitable for
  measuring physiological ATP levels"; "QUEEN-2m showed a Kd of 4.5 mM for ATP at 25 °C"; "the in vivo dynamic
  range of QUEEN-2m (~3.0) … was larger than that of ATeam1.03YEMK (~1.8)"; p.1: "QUEEN was apparently
  insensitive to bacteria growth rate changes". Use: the sensor can resolve physiological ATP; a ratio range
  of 0.10 to 0.16 across all our BSA chambers against a threefold in vivo range says the cells sat in one ATP
  state or the sensor did not respond.
- 2014_Yaginuma p.2 (Fig. 1 legend) and p.3: "Response of QUEEN-2m to ATP concentration at different pH values"
  (pH 6.2 to 8.6) and "to clarify the effect of diversity in intracellular pH on the QUEEN-2m signal,
  ratiometric pHluorin was expressed". Use: the QUEEN ratio is pH-dependent, so under pH oscillation an ATP
  reading needs the pH sensor beside it; a flat QUEEN ratio under pH oscillation is consistent with a
  constant intracellular pH, which is what BSpH reported.
- 2014_Yaginuma p.2: Kd of QUEEN-7m differs between 25 and 37 °C. Use: calibration must be done at the
  cultivation temperature (30 °C), one more reason the BioLector curve does not transfer to the chip.
- 1998_Miesenböck p.2: ratiometric pHluorin "displays a reversible excitation ratio change between pH 7.5
  and 5.5 … with a response time of < 20 ms"; calibration "Cells expressing GPI-anchored ratiometric pHluorin
  at their surface … were imaged in buffers adjusted to pH values between 5.28 and 7.8." Use: the principle
  and the in situ calibration by imaging cells in buffers of known pH, the counterpart of what is missing in
  the chip; the working range covers the cytosolic pH of fungi.
- 2022_Pianale p.14: "we showed a tight connection between ATP concentration and intracellular pH". Use: the
  two ratiometric sensors are coupled; read BSA and BSpH together.

### Claim: the sensor strains grew slower than the wild type in the BioLector (3.1)

- 2022_Pianale p.1 and p.8 to p.9: "Even when combined together, the biosensors did not significantly affect
  key physiological parameters, such as specific growth rate and product yields."; "Neither the growth curves
  … nor the maximum specific growth rates … of strains bearing single biosensors showed any significant
  difference with respect to the parental strain."; "Overall, the selected biosensors did not affect yeast
  growth". Use: the contrast. In S. cerevisiae with genome integration at a defined site the sensors were
  neutral; in A. pullulans the four sensor strains grew at 66 to 90 % of the wild type in oMLP 50 g/L. The
  discussion can offer integration site, copy number, expression level and clonal variation as candidates
  and say that an isogenic empty-vector control is missing.
- 2024_Zhang does not support this claim (see below). A source for fluorescent-protein or heterologous
  expression burden in fungi or yeasts is still needed; alternatively keep the sentence descriptive and cite
  only the contrast with 2022_Pianale.

---

## Papers that fit something else than planned

- **2024_Zhang** is a UV-mutagenesis and screening paper: "we designed a high-throughput screening system
  utilizing a unique fluorescent protein specific to pullulan … flow cytometry allowed for single-cell
  screening" (p.1); the fluorescence is a bimolecular-complementation reporter of pullulan, not an expressed
  sensor whose cost is measured. It says nothing about growth burden. It is usable in 4.5 for the morphotype
  list, "It manifests five distinct cell shapes, including yeast-like cells, chlamydospores, blastospores,
  swollen blastospores, and mycelium. A. pullulans is also called 'black yeast' because it also produces
  melanin in the later stages of fermentation" (p.2), and for "a 10:1 carbon/nitrogen ratio is the most
  favorable condition for pullulan production" (p.2, citing Kumar et al.).
- **2022_Olsson** p.9 gives Kitano's integral definition of robustness; only needed if you want the formal
  background before Trivellin's R.

## Added 2026-10-06: 2004_Campbell and 2026_Malat

- **2004_Campbell** (FEMS Microbiol Lett 232:225): p.1 (abstract): "only swollen cells and chlamydospores, and
  neither hyphae nor unicellular blastospores, often held responsible for pullulan formation, appeared to
  produce pullulan-like material"; p.4: "Only the multi-celled chlamydospores and swollen cells were coated with
  silver grains … None of the blastospores …, germ tubes arising from swollen cells … or hyphae … showed this
  silver deposition"; p.1: "Many factors affect this complex life cycle, with the nitrogen source (both organic
  and inorganic) and medium pH particularly influential"; p.2: "An inverse relationship between intracellular
  glycogen and extracellular pullulan formation in A. pullulans [7] suggested a link between their production."
  Use: 4.5 and 4.7 (morphotype decides the product; a shift toward blastoconidia under glucose heterogeneity
  is a shift away from the producing forms); 4.4 (glycogen as a reserve).
- **2026_Malat** (Sci Rep 16:16672): p.1: fluorescence decline "may be associated with reduced dye penetration
  into the compact, melanized matrix rather than a reduction in biomass"; p.1: morphological forms "ranging from
  yeast-like blastoconidia to filamentous hyphae and melanized arthroconidia"; p.8: quantitative fluorescence
  "may not resolve the mechanical or structural properties of adhering biomass, particularly when EPS
  accumulation limits dye penetration". Use: 4.4 (melanisation as a possible attenuator of fluorescence
  readouts, hypothesis), 4.5 (surface growth and EPS as the context of the chamber biofilm).

## 4.5 to 4.7 (added 2026-10-07: 2025_Pachitariu, 2021_Nadal-Rey, 2022_Losoi, plus passages of the earlier papers)

### 4.5 Growth mode and morphology
- 2026_Rensink p.1 to p.2 (compiled dimensions): yeast cells "with a size of 9-11 x 3-6.5 μm, while the hyphae
  have a width of 2-16 μm (de Hoog et al., 2000; Samson et al., 2019)"; "The globular and ellipsoid swollen cells
  can either or not be septated with an average size of 15 × 11 μm and 12 × 9 μm, respectively, and have a thick
  cell wall (Campbell et al., 2004; Li et al., 2009; Pechak and Crang, 1977)"; chlamydospores "have an average
  size of 13 × 12 μm". Use: map the size classes of 3.2.2 onto the literature (yeast cells about 20 to 55 µm²,
  swollen cells 85 to 130 µm², chlamydospores about 120 µm² as projected areas).
- 2026_Rensink p.1 (abstract): the cell-type cycle and "blastoconidia were no longer the most dominant cell
  type after 10 h of culturing, although they represented more than 75 % of the cells at the moment of
  inoculation"; model prediction "57 %, 32 %, and 11 %" blastoconidia, swollen cells and hyphae after 72 h.
- 2004_Campbell p.1: nitrogen source and medium pH "particularly influential" on the life cycle; abstract and p.4
  on swollen cells and chlamydospores as the pullulan producers (see above).
- 2024_Blöbaum p.8: "Cell circularity varied substantially with 1.5 and 6 min oscillations, owing to a shift
  towards pseudohyphal growth … Pseudohyphal growth is a known response to nitrogen starvation and stress".
  Use: a different trigger of pseudohyphae than ours (complex medium).
- 2024_Haala p.6 (storage lipids under high C/N) and p.10 (oMLP composition, C/N 125); 2024_Zhang p.2 ("a 10:1
  carbon/nitrogen ratio is the most favorable condition for pullulan production").

### 4.6 Measuring a blastoconidia-forming fungus in a chip
- 2025_Pachitariu p.1: "developers may specifically focus on methods that can be proved to generalize well
  out-of-distribution"; "We increase generalization performance further by making the model robust to channel
  shuffling, cell size, shot noise, downsampling, isotropic and anisotropic blur. The new model can be readily
  adopted into the Cellpose ecosystem which includes finetuning, human-in-the-loop training, image restoration
  and 3D segmentation approaches."
- 2025_Pachitariu p.5: "Cellpose-SAM can run out-of-the-box on images that have been acquired with varying
  levels of image degradation, at different pixel sizes or in arbitrary channel order, substantially simplifying
  the logistics typically associated with setting up an image segmentation pipeline."
- 2025_Pachitariu p.6: "the model is especially good at generalizing with zero-shot or limited data."
- 2025_Pachitariu p.11 to p.12: training data include phase-contrast yeast (YeaZ, "16 2D images from the phase
  contrast dataset") and phase-contrast bacteria (Omnipose). Use: yeast in phase contrast is inside its
  training distribution, a black yeast with 270-µm² cells and pseudohyphal chains at its edge.
- 2024_Blöbaum p.3: "In perfusion-based microfluidic systems, a maximum of 150–1000 microbial cells can be
  cultivated and trapped in one monolayer-growth chambers". Use: the density the yeast chambers reach against
  the 180 objects at most in ours.
- 2020_Täuber p.13 (position-dependent exchange) and p.15 (identical-media switch) as under 4.1 and 4.2.

### 4.7 Implications for bioprocess design
- 2021_Nadal-Rey p.2: mixing time "defined as the time needed for the concentration of an inert tracer to reach
  a pre-defined degree of homogeneity … with values being of the order 10–1000 s depending on the reactor design
  and operating conditions (Sweere et al., 1987)"; "Use of the circulation time will give a lower value of τ
  than the mixing time (Doran, 1995), typically by a factor of three to five."; "gradients are more likely to
  occur at high biomass concentrations, and/or for microorganisms which have high substrate and/or oxygen uptake
  rates."
- 2021_Nadal-Rey p.3: "The changes that can occur as a response to such fluctuations vary over a wide range of
  timescales (from seconds to hours). Changes in the metabolome … can occur at the second time scale."; "It has
  been estimated that exposure to fluctuating conditions increases the ATP demand for maintenance by 40–50% in
  E. coli (Löffler et al., 2016). The practical consequence of this is that less carbon is available for product
  formation, leading to a reduction in the product yield."; "Using regime analysis, substrate gradients are
  expected in the majority of cases".
- 2021_Nadal-Rey p.7 and p.9 (Table 3): S. cerevisiae, 20 to 22 m³, four impellers, mixing time 147 to 183 s;
  "A substrate concentration peak from approximately 40 to 80 mg L−1 was observed, being 80 mg L−1 the maximum
  local concentration of glucose reported." Use: the amplitude of real fed-batch gradients against the 50 to
  0 g/L alternation of the chip.
- 2021_Nadal-Rey p.10: scale-down simulators apply pulses "at such frequency and magnitude that the simulator
  is able to mimic the oscillations that the cells experience at a particular location in a large-scale
  reactor"; "Many of these issues are amplified when working with viscous fermentation broths (e.g. filamentous
  fungi broth)"; "Tip speed affects the shear rate, causing changes in the rheology and morphology of
  filamentous fungi". Use: viscosity of pullulan broths as the scale-down problem the chip avoids.
- 2022_Losoi p.1: "insufficient macromixing of feeds, leading to heterogeneities in pH, substrate, and oxygen";
  multipoint feeds "reduced the mixing time substantially by more than a minute and mitigated gradients of pH,
  substrate, and oxygen"; "the phenotypical heterogeneity of the biomass population was diminished".
- 2022_Losoi p.9: "the pH gradients 10 s after the carbonate pulse … were similar to the substrate gradients."
- 2022_Losoi p.10: "Both R4 and B13 achieved mixing times of 10 s (Table 3), which are common in
  laboratory-scale reactors"; p.14: "Simulated mixing times were reduced from the scale of minutes to the scale
  of 10 s, which mitigated substrate gradients and restored ideal homogeneous reactor performance."; p.17:
  "locally measured mixing times of 165 s in R4 … and 124 s in R1".
- 2025_Pianale p.3: oscillations of 24 and 48 min "included to assess the physiological effects of longer
  oscillations in the order of large-reactor mixing times (Lara et al. 2006)".
- 2022_Olsson p.1: "Quantification of robustness is suggested as a tool for guiding strain and bioprocess
  development."

## Still missing after these seventeen

- Conidiation or blastoconidia release under nutrient limitation in A. pullulans (4.2, 4.4): Rensink 2024
  and 2026 are the fallback; a physiological source would be better.
- Fluorescent-protein or heterologous-expression burden in a fungus (4.4, 4.7).
- Mixing and circulation times of large reactors with numbers (4.7): your intro refs 2022_Losoi,
  2021_Nadal-Ray, 2025_Arulrajah; 2020_Täuber p.3 gives only "Typical time scales are in second to minutes
  range based on the size and geometry of the bioreactor".
