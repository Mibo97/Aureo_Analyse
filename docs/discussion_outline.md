# 4. Discussion: outline and evidence

Working plan of 2026-10-06 for the discussion chapter; all sections drafted in `docs/discussion_chapter.md` by 2026-10-07. Every section lists the claims, the evidence in the
results chapter (numbers from the final v12 run, `docs/data_story.md`), the interpretations to weigh, and what
the paragraph needs from the literature. "Intro refs" are keys of the introduction's reference list; "to find" are
papers not yet in the list. The chapter has no length limit; the Conclusion is a separate chapter.

Decisions so far: the chamber medium-exchange time is taken from Blöbaum 2024 (no fluorescein test of our own),
and the exchange time and the position within the array are discussed as future work, not analysed further.
The biosensors were validated in the BioLector only; in the chip they are uncalibrated.

---

## 0. The question and the one-paragraph answer

Introduction 1.4: "How robust is *A. pullulans* under controlled glucose concentration and pH oscillations, and
what are the implications for bioprocess design in a circular bioeconomy context?" Plus the two method problems
named there: mother-bud assignment with moving blastoconidia, and pullulan in a flow device.

Answer to carry through the chapter: at the hour scale, growth, morphology and sensor signals of all five strains
did not depend on the half-cycle period between 0.75 and 24 min; every apparent dose response was carried by the
structure. Glucose oscillation itself shifted the population toward blastoconidia formation and larger cells
compared with constant feast or famine on the same structure; pH oscillation changed nothing. The robustness
metrics confirm this. The static comparison shows that the readouts respond when the medium differs strongly
(complex against minimal). The sensors are expressed and detected but their dynamic range in the chip is not
established. Pullulan does not explain the disagreement of control chambers.

---

## 4.1 The period did not matter; the structure did

**Claims.**
1. No readout changed with the half-cycle period beyond what the constant-medium controls of the same structures
   did (3.4, 3.5, 3.6: 58 readout-by-series combinations, 31 no trend, 23 structure effects, 2 not robust,
   2 period effects; Table 5, Table 6, Figure 15).
2. The two period effects (µ_area in BSO/Glc, budding rate in WT/Glc) did not recur in the other oscillation type
   of the same strain, and their number lies below what the classification produces by chance (3.5 on average
   when the periods are shuffled among the structures of a series; section "Numbers"), whereas the 23 structure
   effects exceed the chance level (13.5 on average, p 0.002): the controls trend with the real period order more
   often than with a random one, which is the alignment of period, position and day made visible.
3. Structure, cultivation day and chip position were aligned with the period (3.6 last paragraph: day blocks,
   Spearman of period against date ±0.87 to ±0.89, shortest period on the same structure position).
4. Within cultures that carried two or three periods, neither area nor µ_bud changed in a consistent direction
   (8, 13, 6 of 19 and 8, 8, 11 of 19).
5. The controls of a structure differed from the controls of the next structure more than the period moved
   anything: R(p) of the area differed between structures in 20 of 20 Kruskal-Wallis tests (3.6).

**Interpretation.**
- What a structure effect is physically: one preculture, one loading, one day, one position on the chip, one
  flow path. Oscillation and control chambers share it; the period does not.
- This is robustness in the sense of Trivellin et al.: the phenotype is stable across the perturbation space
  (switching intervals). Say clearly what it is not: the design cannot detect period effects smaller than the
  between-structure variation, and the feast and famine controls did not span a range (bracket degenerate on 44
  of 49 structures for budding), so effect sizes cannot be expressed relative to the feast-famine difference.
- What the finding teaches about dMSCC experiments: without the control chambers on the same structure, 27
  trends would have been reported as dose responses. Randomise periods over structures and days, run every
  period on at least two cultures, keep constant-medium chambers on every structure.

**Literature need.**
- Definition of microbial robustness and of R: intro refs 2022_Trivellin, 2024_Blöbaum; to find: Olsson et al.
  review on robustness in bioprocesses (same group as Trivellin).
- How dMSCC studies replicated and controlled: 2024_Blöbaum, 2023_Blöbaum; to find: Täuber et al. 2020
  (dMSCC, *C. glutamicum*, oscillation periods).
- Pseudoreplication and technical replicates: to find: Hurlbert 1984, Lazic 2010 (or any statistics text you
  have on nested designs).

---

## 4.2 Glucose oscillation shifted growth toward blastoconidia; pH oscillation changed nothing

**Claims.**
1. Paired over structures, oscillation chambers budded more than the control mean on 27 of 29 Glc structures
   (ratio 1.25), had higher µ_bud (1.15), larger cells (1.14), marginally rounder cells (0.98) and slower area
   growth (0.88); in every strain; independent of the period; immigration equal (3.4).
2. The oscillation chambers lay above both controls on 26 of 49 structures and below both on 4: not an average
   of feast and famine.
3. Under pH oscillation no readout differed from the controls (ratios 0.97 to 1.09, p 0.09 to 0.65).
4. In the BioLector the wild type grew at pH 3.0 and pH 8.0 at nearly the same rate (0.36 and 0.34 h⁻¹) (3.1).

**Interpretation, three readings to weigh.**
- Transitions as a trigger: every feast-to-famine switch is a nutrient downshift; in fungi nutrient limitation
  induces conidiation, and in *A. pullulans* blastoconidia are released from swollen cells (Rensink). Cells that
  see 20 to 700 downshifts in 20 h may release blastoconidia more often than cells in constant medium of either
  kind. Fits "above both controls".
- Time-averaged medium: the oscillation chambers see 25 g/L glucose on average; but an average effect would place
  them between the controls, which they are not (claim 2). Weak.
- Position: the oscillation chambers occupy A3 to A12, the controls A1, A2, A13 and A14 of the same array
  (methods 2.5.3). Flow, seeding density and focus may differ between the centre and the ends of an array. The
  design cannot separate this from the oscillation. Future work: constant medium in central chambers, or
  oscillation at the ends.
- pH: the external pH alternated between 3 and 8 and nothing moved, not even the pH sensor ratio (4.4). Either
  intracellular pH homeostasis (HOG pathway, cation transporters, intro 1.2) or an unresponsive sensor; the
  growth readouts argue for homeostasis, since the BioLector growth rates at the two pH values were alike.

**Literature need.**
- Responses of other organisms to glucose oscillations in the same chip type: 2024_Blöbaum (*S. cerevisiae*),
  to find: Täuber 2020 (*C. glutamicum*); scale-down reactors: to find: Limberg 2017, Löffler 2016, Lara 2006
  "Living with heterogeneities", Neubauer & Junne 2010.
- Conidiation under nutrient limitation in *A. pullulans* or related fungi; blastoconidia release from swollen
  cells: 2024_Rensink, 2026_Rensink; to find: Ramos & García Acha 1975 (cited by Rensink).
- pH tolerance and homeostasis: 2014_Gostinčar, 2019_Gostinčar, 2009_Schoch; to find: effect of pH on
  *A. pullulans* morphology and pullulan production (optimum around pH 6.5; morphology shifts at low pH).

---

## 4.3 What the robustness metrics add

**Claims.**
1. Under oscillation the mean cell area of a chamber was less stable over time and the population more
   heterogeneous in area and eccentricity, but more homogeneous in area growth rate; the budding output of the
   mothers was as heterogeneous as under constant medium, the budding rhythm of the individual mother less stable
   under glucose oscillation (Table 4; 3.4).
2. Against the period the metrics behaved like the readouts: 90 combinations, 51 no trend, 23 structure effects,
   8 not robust, 8 period effects, five of them in BSA/Glc.
3. Heterogeneity of cell area and of growth rate were unrelated across chambers (ρ 0.06, n = 497).
4. The sensor ratios did not differ in R(t) or R(p) between oscillation and controls (Table 5).

**Interpretation.**
- More budding events mean more small new objects and more mothers in transition, so a wider size distribution
  and a less stable chamber mean are the expected echo of claim 1 of 4.2, not a separate effect.
- A more homogeneous area growth rate under oscillation can be read as synchronisation by the medium switches;
  say that this is a hypothesis.
- What R measures and what it cannot: R is the Fano factor over the global mean, dimensionless and comparable
  across readouts, but R(t) depends on the sampling (10-min frames alias the cycles differently per interval), so
  only the comparison against the controls of the same structure is interpretable, not the level.
- Five of eight period effects in one series (BSA/Glc) is a series property, not a metric property.
- Checked 2026-10-06, confirmed by the cluster run of 2026-10-08 (`51_osc_vs_controls_cv.csv`): as coefficient of
  variation, the growth-rate homogeneity (p 0.002), the eccentricity heterogeneity (p 0.001) and the temporal
  instability of the mean area (p 0.042) and of the mean eccentricity (p 0.005) hold; the area heterogeneity
  does not (p 0.13), it is the larger mean size expressed as a Fano factor. All differences are Glc-only.

**Literature need.** 2022_Trivellin (robustness is phenotype-specific; how R was compared across conditions),
2024_Blöbaum (use of R(t) and R(p) in the chip); to find: phenotypic heterogeneity reviews (Ackermann 2015) if
you want the heterogeneity angle.

---

## 4.4 Famine did not starve, and the sensors did not separate

This is the central caveat of the chapter; give it its own section.

**Claims.**
1. Feast and famine controls barely separated: bracket degenerate on 44 of 49 structures for the budding rate,
   feast above famine on 26 of 48 (p 0.61); cells under constant famine budded in the sparse phase as often as
   under constant feast (0.05 to 0.63 against 0.00 to 0.94 buds per mother-hour); only area and eccentricity had a
   direction (feast larger on 31 of 48, more eccentric on 33 of 48, p 0.01) (3.4).
2. In the BioLector, oMLP without glucose supported almost no growth (biosensor strains µ_max ≤ 0.02 h⁻¹, wild
   type 0.09 h⁻¹ after a 12 h lag) (3.1).
3. No sensor separated its feast from its famine chambers (degenerate on 6 to 9 of 10 structures per sensor,
   no direction in the pairwise comparison; Figure 17) (3.5).
4. The sensors are expressed and both channels are detected in every strain (Figure 13); the Gly-RNA ratio of
   one feast structure at 1.00 is an outlier to explain or exclude.

**Interpretation.**
- Reserves and the preculture: the sparse-phase window covers most of the 22 h, and cells came from a glucose
  preculture; blastoconidia formation may draw on reserves (lipids, glycogen) and is itself a starvation
  response, so the budding readout cannot separate feast from famine within 22 h. Cell size can, and did.
- Did the famine medium reach the cells? The control arrays were supplied at a passive 10 mbar against 90 mbar
  on the oscillation line (methods 2.4.3); the exchange time of the chambers is taken from Blöbaum 2024. Name
  this as an untested assumption and as future work (fluorescein switching test, a glucose indicator in the
  chamber).
- Sensors: without an in-chip calibration (ionophore clamp for pHluorin, 2-deoxyglucose or azide for ATP
  depletion, H₂O₂ for OxPro) the absence of a feast-famine difference cannot be read as physiology. The BioLector
  validation shows the strains grow and fluoresce; it does not show that the ratio moves in the chip. Therefore
  the sensor results are negative results about the measurement, not about the cells.
- Consequence for the thesis question: robustness is demonstrated for growth and morphology against the period
  and against oscillation itself; for intracellular physiology it is not demonstrated, only not contradicted.

**Literature need.** 2022_Pianale, 2025_Pianale (sensor dynamic ranges and calibration in yeast), 2024_Blöbaum
(sensor use in the chip, exchange time); to find: Yaginuma 2014 (QUEEN), Miesenböck 1998 (pHluorin), a paper on
storage compounds or starvation survival in *A. pullulans* (lipid accumulation, liamocins: 2024_Haala), and on
fluorescent-protein burden on growth (for 4.7).

---

## 4.5 Growth mode and morphology

**Claims.**
1. Typical cell 22 µm², axis ratio 1.5; 5 % large and round; Glc series slightly larger and more elongated than
   pH; wild type smallest, BSpH largest (3.2.2).
2. Complex medium made swollen cells: on W109 82 % of cells ≥ 30 µm², median 82 µm²; cells 3.5 times larger and
   budding a sixth as often as in minimal medium (3.3). W65 did not repeat the minimal-medium side; one chamber
   with two 270 µm² cells releasing ten blastoconidia in three frames (Figure 5 B).
3. Buds are blastoconidia of 4.7 µm² released from mothers of 71 µm², more than twice the median cell; the
   mother size holds for all candidates, so it is not made by the size criterion.
4. µ_bud and µ_area are not coupled in one direction; the Glc structures combine more budding with slower area
   growth than the pH structures (Figure 12).
5. Pseudohyphal chains in complex medium on W65 (Figure 5 A); true hyphae only on W109 in complex medium.

**Interpretation.**
- Two growth modes, size and blastoconidia: complex medium drives swelling and hyphae, minimal medium
  blastoconidia; glucose oscillation moves the population toward blastoconidia and larger mothers at once.
- Map the size classes onto the literature dimensions: yeast cells 9 to 11 × 3 to 6.5 µm correspond to about
  20 to 55 µm² (our median 22 µm²), Rensink's swollen cells 15 × 11 µm to about 130 µm², chlamydospores 13 ×
  12 µm to about 120 µm²; our median mother at bud release (71 µm²) is a 9.5 µm sphere, the W65 cells of 270 µm²
  are 18.5 µm, larger than the published swollen cells.
- Rensink's cell-type dynamics (blastoconidia divide, differentiate into swollen cells, swollen cells release
  blastoconidia and form hyphae) give the frame; our data add the medium dependence and the oscillation shift.
- The W109-W65 disagreement is the culture problem of 4.1 in miniature: one culture per chip.
- Bridge to the product: which morphotype produces pullulan (Campbell 2004) decides whether a morphology shift
  under glucose heterogeneity matters for the process.

**Literature need.** 2024_Rensink, 2026_Rensink (cell types, dimensions, dynamics; both are in the list, cite
each for what it contains); to find: Campbell et al. 2004 FEMS Microbiol Lett 232:225 (morphotype and pullulan),
de Hoog 2000 (dimensions, cited by Rensink), medium and C/N effects on *A. pullulans* morphology, 2024_Haala
(oMLP and liamocin production context).

---

## 4.6 Measuring a blastoconidia-forming fungus in a chip

**Claims.**
1. The tracker fragments at 0.070 new tracks per object-frame (half the previous version); median track 5 frames
   but 78 % of object-frames in tracks of ≥ 10 frames (3.2.3).
2. Above about 8 to 20 objects per frame, new objects touching a cell scale with density and are mostly mask
   fragments; lineage readouts are therefore restricted to the sparse phase (median window 117 of 133 frames;
   245 chambers never left it) (3.2.4, Figure 7).
3. The size criterion from a bimodal ratio (modes 0.07 and 0.45, antimode 0.32) rejected 15.5 % of candidates;
   assignment rate 0.67 per chamber; 74 % of mothers by mask contact (3.2.5).
4. The budding trends against the period depended on the tracker (per-series correlations of the two trackings
   ρ 0.05; per-structure rates ρ 0.77) (3.2.6).
5. Immigration (0.19 per cell-hour) rivalled births (0.21); 37 % of chambers had more immigrants than births, so
   object counts are not growth curves (3.4).
6. The pullulan knockout did not agree better among its control chambers than the producers (3.7, Figure 16).

**Interpretation.**
- What the platform can measure for this organism: hour-scale cumulative readouts (endpoint morphology,
  sparse-phase budding, µ_bud, robustness metrics), not lineages over generations and not intra-cycle dynamics
  at 10-min frames. Both intro problems are answered: mother-bud assignment by mask contact plus a size
  criterion, valid in the sparse phase; pullulan is not what makes control chambers disagree.
- Blastoconidia wash-out makes µ_bud a lower bound and lineage depth shallow; the open chamber design that lets
  blastoconidia leave is also what keeps the density low for 22 h.
- Rule-based filtering replaced the manual curation; the rules are documented and the thresholds are logged with
  every run (methods 2.5).
- The tracker dependence of per-series trends is a warning for any lineage readout from automated tracking:
  report structure-level agreement, not only series-level trends.
- Alternatives for the control disagreement after the knockout result: seeding density, position, flow asymmetry
  between arrays; future work: fluorescein exchange test, position test, seeding-density control.

**Literature need.** Cellpose and Cellpose-SAM (to find: Stringer 2021, Pachitariu 2025), yeast trackers for
comparison (to find: YeaZ, DeLTA, TrackMate), dMSCC chamber design and exchange time (2024_Blöbaum, to find:
Täuber 2020), 2019_Chen (AmAgs2 knockout), pullulan rheology if you want to argue viscosity.

---

## 4.7 Implications for bioprocess design in a circular bioeconomy

**Claims to connect.**
1. Mixing and circulation times of large reactors (intro refs 2022_Losoi, 2021_Nadal-Ray, 2025_Arulrajah) lie at
   the short end of the half-cycles tested; the 0.75 to 24 min intervals bracket feeding and circulation cycles.
2. At the hour scale the organism's growth and morphology did not depend on the switching frequency, in all five
   strains; glucose heterogeneity shifted the morphology distribution; pH heterogeneity between 3 and 8 was
   tolerated.
3. The biosensor strains paid for the sensor: 66 to 90 % of the wild-type growth rate in oMLP 50 g/L, and the
   ranking depended on the medium (3.1).

**Interpretation.**
- For a circular-feedstock process, insensitivity to glucose and pH fluctuation frequency is the asset the
  introduction asked for; the morphology shift under glucose heterogeneity is the cost to watch if the
  morphotype decides the product (pullulan, liamocins, enzymes).
- Scale-down studies with *A. pullulans* should vary the amplitude and the medium, not only the frequency, and
  replicate per culture.
- Monitoring strains with intracellular sensors are possible in this organism but carry a growth burden and need
  in-chip calibration before they report physiology.

**Literature need.** 2022_Losoi, 2021_Nadal-Ray, 2025_Arulrajah, 2023_Zhang (gradients and *A. pullulans*),
2021_He, 2022_Wan, 2011_Manitchotpisit (products), 2024_von_Braun, 2021_Lange (bioeconomy frame); to find:
Haringa 2016 (lifelines), scale-down reviews.

---

## 4.8 Limitations and future work (consolidated)

- Biological replication: every period on at least two or three cultures; randomise periods over structures,
  positions and days.
- Chamber medium exchange: fluorescein switching test for every interval, especially 0.75 and 1.5 min; a glucose
  indicator to verify the famine state in the chamber.
- Position within the array: constant medium in central chambers or oscillation at the ends, or a dedicated
  position test with identical medium in all chambers.
- Sampling: 1-min frames for a subset of chambers to see the cycles of the 12- and 24-min intervals and the
  damping of the short ones.
- Sensors: in-chip calibration of all four sensors; exclude or explain the Gly-RNA outlier structure.
- Static branch: more chips per medium; the W65 minimal-medium disagreement.
- BioLector: OD₆₀₀ calibration of the scattered light; biological replicates of the growth parameters.
- Pipeline: bud retention (trap geometry) to deepen lineages; mixed models for structure-level inference.

---

## Conclusion (separate chapter, one paragraph plan)

Question, design in one sentence, the three findings (period did not matter beyond the structure; glucose
oscillation shifted growth toward blastoconidia and larger cells, pH did not; sensors and static comparison),
the two method answers (mother-bud assignment, pullulan), the central caveat (famine and sensor range), and the
one recommendation for scale-down work with this organism.

---

## Literature workflow (how to find the fitting parts)

Done for 4.1 and 4.4 on 2026-10-06: `docs/literature_4.1_4.4.md` lists the passages of the twelve papers you
sent, per claim, with PDF page and quote, plus the papers that fit something else and what is still missing.

For every claim above the "literature need" names what the sentence must be able to cite. Two ways to get there:

1. Upload the PDFs of the candidates (or the whole library as a zip). I extract the passages that match each
   claim with page numbers and the exact wording, so that you decide what to cite from a short list instead of
   re-reading papers.
2. Alternatively export the reference list (BibTeX or RIS) with abstracts; I mark which references fit which
   claim and which claims have no source yet.

Search keys per section, for your own library search: 4.1 "robustness" "Fano" "pseudoreplication" "technical
replicate"; 4.2 "oscillation" "feast famine" "scale-down" "conidiation" "nutrient limitation" "intracellular pH";
4.3 "phenotypic heterogeneity" "robustness metric"; 4.4 "biosensor calibration" "dynamic range" "pHluorin"
"QUEEN" "starvation" "storage lipid"; 4.5 "swollen cell" "blastoconidia" "chlamydospore" "pullulan morphology"
"yeast-like mycelial"; 4.6 "Cellpose" "cell tracking" "microfluidic single-cell" "medium exchange"; 4.7 "mixing
time" "circulation time" "gradient" "lifeline" "scale-down".

---

## Numbers

- Permutation null of the control-trend classification (`docs/scratch/trend_null.py` on the final-run tables,
  1,000 shuffles of the periods among the structures of each series; a structure keeps its oscillation chambers
  and controls together, so chip effects survive and only the period order is destroyed):

  | verdict | observed | null mean | null 5th to 95th percentile | P(null ≥ observed) |
  | --- | --- | --- | --- | --- |
  | no monotone trend | 31 | 38.3 | 32 to 44 | 0.98 |
  | structure effect | 23 | 13.5 | 9 to 19 | 0.002 |
  | not robust | 2 | 2.7 | 0 to 5 | 0.78 |
  | period effect | 2 | 3.5 | 1 to 7 | 0.88 |

  Reading: two period effects are fewer than chance produces, so nothing in the 58 combinations needs a period
  to explain it. The excess of structure effects says that the structure drifts (day, position, loading order)
  were monotone in the real period order more often than in a random order; this is the alignment of 3.6
  quantified, and the strongest argument for randomising periods over structures and days.
- Position within the array and the chamber medium-exchange time: not computed, future work (4.8).
- OD₆₀₀ calibration of the BioLector data: pending your end-point readings.
