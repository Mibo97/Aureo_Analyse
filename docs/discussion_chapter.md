# 4. Discussion

Draft of 2026-10-06, sections 4.1 and 4.4; the other sections follow the plan in `docs/discussion_outline.md`.
Cross-references use the thesis numbering of the results chapter (Tables 1 to 6, Figures 2 to 17, Sections 3.x)
and of the methods chapter (2.x). Citation keys are those of your reference list; the page and the full passage
behind every quotation are in `docs/literature_4.1_4.4.md`. Direct quotations can be paraphrased where the
thesis style prefers it. Sentences in *[italic square brackets]* are notes to the author, not thesis text.

This chapter answers the question of the introduction, how robust *A. pullulans* is under glucose and pH
oscillations and what follows for bioprocess design, in the order in which the results allow it: first what the
period did and did not do and what carried the trends (4.1), then the one effect of the oscillation itself
(4.2) and the robustness metrics (4.3), then the two caveats that limit the answer, the controls that did not
separate and the sensors whose range in the chip is unknown (4.4). Morphology and growth mode (4.5), the
method (4.6) and the implications for bioprocess design (4.7) follow, and the limitations are collected in 4.8.

## 4.1 The period did not matter; the structure did

Across half-cycle periods of 0.75 to 24 min, that is full feast/famine cycles of 1.5 to 48 min, no growth,
morphology or sensor readout of the five strains changed with the period once the constant-medium controls of
the same structures were taken into account. Of the 58 readout-by-series combinations, 31 showed no monotone
trend of the oscillation chambers, 23 a trend that a control of the same structures shared, 2 a trend that
vanished after the controls were subtracted and 2 a trend that met all conditions of a period effect (Table 6,
Figure 15). Shuffling the periods among the structures of each series produced 3.5 period effects on average,
so the two observed ones lie below chance level, and neither recurred in the other oscillation type of the same
strain. The 23 structure effects, by contrast, exceeded the 13 that shuffled periods produced: the control
chambers trended with the real period order more often than with a random one. The result is as negative as it
is clear. Whatever the period did to the cells over 20 h was smaller than what the structure did.

In the terms of Trivellin et al., who define robustness "for a specific function (or phenotype) and set of
perturbations" and separate it from tolerance, which "specifically describes stable growth or survival to
various perturbations" (2022_Trivellin), the perturbation space of this thesis is the set of switching
intervals, and the functions are cell size and shape, budding, area growth and the sensor ratios. Every one of
these functions was stable across the perturbation space. The BioLector cultivations of 3.1, in which the
strains grew at pH 3.0 and pH 8.0 and with and without glucose, describe tolerance: "tolerance only considers
viability and/or specific growth rate", whereas robustness "encompasses the stability of the production and
growth of a microbe in different conditions evaluated using different measures" (2022_Olsson). The chip data
describe robustness.

What a structure effect is, physically, follows from how the experiments were made. A structure carried one
period, one loading of cells from one preculture on one day, one position on the physical chip and one flow
path. Its oscillation chambers and its constant-medium chambers share all of this; the period they do not
share. When the feast and famine controls of a series rose or fell with the period as the oscillation chambers
did, the common factor was the structure, and the period order was aligned with it: the periods were run in day
blocks (Spearman of period against date ±0.87 to ±0.89), and the shortest period usually occupied the same
structure position (3.6). The permutation result makes this alignment measurable. The structure drifts were
monotone in the real period order more often than in a random one.

The unit of the experiment follows from the same fact. Lazic defines the experimental unit as "the smallest
entity that can be randomly assigned to a different treatment condition" and separates "biological replicates
which are independent" from "technical replicates which are not independent" (2010_Lazic). The entity that
could be assigned a period was the structure. Its five oscillation chambers are sub-samples of one culture, and
"measuring the height of one person ten times provides only one independent piece of information about the
population" (2010_Lazic). A test over chambers against the period would have treated correlated observations
as independent, and "the greater the correlation the greater the number of false positives" (2010_Lazic). This
is why the analysis compared one value per structure across the periods, why the Spearman correlation with n of
four to six is reported as an effect size and not as a test, why the comparison of oscillation and constant
medium was paired within structures (3.4), and why the change within cultures that carried two or three periods
was examined separately (3.4). Averaging the observations of a unit and random-effects models are the two
remedies Lazic lists (2010_Lazic); the first was used here, the second needs more cultures per period than
this data set holds.

The published dMSCC studies share the unit. The yeast chip of Blöbaum et al. "consists of six connected
cultivation structures, that include one dedicated positive control and five oscillation-structures", each
frequency on its own structure, and the dispersion of their data "corresponds to the standard deviation across
triplicates (three chambers)" (2024_Blöbaum); Pianale et al. report each frequency over five replicate chambers
(2025_Pianale); and Täuber et al. validated the switching itself with an oscillation between two identical
media, which left the growth rate unchanged (2020_Täuber). None of this is a fault of those studies or of this
one. One period per structure is what the laminar layout of the chip allows, in which "control zones, which are
always exposed to the same liquid" flank "a switching zone, across which the laminar boundary layer moves"
(2023_Blöbaum). The difference here is that every structure carried constant-medium chambers of both kinds, so
that the drift of the structure could be read off and subtracted. Without them, the 27 monotone trends of the
oscillation chambers in Table 5 would have been reported as dose responses, among them the rise of the OxPro
and sfpHluorin ratios with the period in both oscillation types (3.5). The control-trend classification of 3.6
is therefore the methodological result of this thesis that carries beyond *A. pullulans*: in a dMSCC
experiment the constant-medium chambers on the same structure are the only handle on the structure, and the
period can be credited only with what exceeds them.

Two things the result does not say. First, it does not exclude period effects smaller than the variation
between structures. The feast and famine chambers of one structure differed from those of the next structure
more than any period moved the oscillation chambers; the population robustness of the cell area differed
between the structures of a series in 20 of 20 Kruskal-Wallis tests (3.6). This variation is the detection
limit of the design. Second, the controls did not span a range. For the budding rate the bracket was degenerate
on 44 of 49 structures, so no effect could be expressed as a fraction of the feast-famine difference (4.4). The
robustness metrics add the same caveat in their own terms. R is "a relative and not an absolute term" that
"changes upon addition of more replicates or conditions" (2024_Blöbaum), normalised by the mean of the function
in the data set, which "adjusts R to the strains tested in the study" (2022_Trivellin); only comparisons within
this data set, oscillation chambers against the controls of their structure, are interpretable, not the levels
of R. And R(t) "is not able to differentiate between oscillating and steadily changing functions"
(2024_Blöbaum), which at 10 min per frame and half-cycles of 0.75 to 24 min means that an R(t) difference
between periods would not have been interpretable in any case; against the controls it is.

Whether the cells saw the short half-cycles at all is a question these data cannot answer, but the chip
literature can. In the *C. glutamicum* chip of the same family, "95% is reached in the cultivation chamber
after 5 seconds" of a medium switch, and "medium switches up to 5 seconds can be performed, guaranteeing an
almost complete medium exchange within the cultivation chamber" (2020_Täuber); the protocol regards a
measurement of the exchange as necessary only for oscillations faster than 10 s or for new chamber heights
(2023_Blöbaum). The shortest half-cycle used here, 45 s, lies an order of magnitude above that limit, so the
cells received the full amplitude at every interval. The chamber height of the yeast-type chip and the passive
supply of the control arrays (2.4.3) differ from the measured configuration, which is why the fluorescein test
stays in the outlook (4.8).

Against this background the flat response of *A. pullulans* stands out. In *S. cerevisiae* exposed to the
same intervals, "a decrease in specific growth rate but an increase in intracellular ATP levels" was observed
"with longer oscillation intervals", and the longest cycles gave the least stable and most heterogeneous ATP
levels (2024_Blöbaum). In *C. glutamicum* the colony growth rate rose from about 0.3 h⁻¹ at intervals of 5 min
to 1 h to 0.5 to 0.6 h⁻¹ at intervals of 1 min and shorter, and the mean cell length changed with the interval
(2020_Täuber). Both organisms responded to how often the medium switched. *A. pullulans*, over 20 h and
against its own controls, did not. Section 4.2 shows that it did respond to whether the medium switched at all.

The design lessons are concrete. Periods should be assigned to structures and days at random and not in
blocks; every period should run on at least two cultures; constant-medium chambers of both kinds belong on
every structure; and the statistics belong at the level of the structure. A control with constant medium in
central chambers, or a switch between two identical media as in Täuber et al. (2020_Täuber), would separate
the chamber position and the pressure switching from the oscillation (4.2). With these, the variation between
structures stops being the detection limit and becomes a measured quantity.

## 4.2 Glucose oscillation shifted growth toward blastoconidia; pH oscillation changed nothing

*[To be drafted; plan and evidence in `docs/discussion_outline.md`, passages in `docs/literature_4.1_4.4.md`
under 4.2.]*

## 4.3 What the robustness metrics add

*[To be drafted; plan in the outline, passages under 4.3 of the literature file.]*

## 4.4 Famine did not starve, and the sensors did not separate

The weakest point of the design is not the period but the controls. The constant-feast and constant-famine
chambers of a structure, which received 50 or 0 g/L glucose, or pH 3.0 or pH 8.0, for the 22 h of a run, barely
separated in the readouts that carry the argument. For the budding rate the bracket was degenerate on 44 of 49
structures, and the feast control exceeded the famine control on 26 of 48; cells under constant famine budded in
the sparse phase as often as cells under constant feast (3.4). Only the morphology had a direction, feast cells
being larger on 31 and more eccentric on 33 of 48 structures. None of the four sensors separated its feast from
its famine chambers (3.5, Figure 17). The same medium without glucose supported almost no growth in the
BioLector: the biosensor strains stayed at their initial signal, and the wild type grew at 0.09 h⁻¹ after a lag
of 12 h (3.1). In the chip, 22 h of the same medium did not stop blastoconidia formation. Three readings have to
be weighed, and they lead to different consequences for the oscillation results.

The first reading is reserves. The cells entered the chip from a preculture in oMLP *[confirm: with 50 g/L
glucose]*, a medium derived from the liamocin production medium of Haala et al.; with 50 g/L glucose and
0.46 g/L ammonium nitrate it has the carbon-to-nitrogen ratio of 125 of that production medium (2024_Haala).
Under such conditions, "high carbon (C) and low nitrogen (N) concentrations (high C/N ratio)", *A. pullulans*
forms "intracellular storage lipids" (2024_Haala), which in oleaginous *Aureobasidium* strains "can reach 65%
of the cell dry weight" (2024_Haala, citing Xue et al. 2018); the strain of that study is the wild type used
here, NRRL 62042. An inverse relation between intracellular glycogen and extracellular pullulan has been
described as well (2004_Campbell). Cells that arrive with such reserves can bud for many hours without external
carbon, and blastoconidia formation is in any case the response that *A. pullulans* showed whenever growth in
size stopped: in the static chambers, complex medium produced large cells and few blastoconidia and minimal
medium small cells and many (3.3), and under glucose oscillation the cells budded more and grew in area more
slowly than under constant feast (3.4). Under this reading the famine medium did reach the cells, but the
readout that was expected to fall, the budding rate over the sparse phase, is the one readout that a starving
*A. pullulans* population keeps up, whereas cell size, which needs new material, responded in the expected
direction. The oscillation results of 3.4 then stand, and the missing separation of the controls is a property
of the organism and of the 22-h window, not of the measurement.

The second reading is that the famine medium was not famine. The medium without glucose still contained the
salts, the nitrogen source, the vitamins and the trace elements (2.1), and the preculture medium entered the
chip with the cells. Pianale et al. changed their own starvation condition from 0 g/L to 10 mg/L glucose
between two studies and note that "more than 500 genes have been found to be differentially expressed between
starving and calorie-restricted (i.e., enough substrate for metabolic activity, but not reproduction) cells"
(2025_Pianale): what famine means for a cell depends on residual carbon in a range that the chip experiments did
not measure. A glucose indicator in the chamber or a sampled outlet would settle this.

The third reading is transport. The control arrays were supplied at a passive 10 mbar against 90 mbar on the
oscillation line (2.4.3), and the exchange time of the chambers was not measured here. The published values
argue against a transport limitation: in the *C. glutamicum* chip the chamber reaches 95 % of a new medium
within 5 s (2020_Täuber), and the protocol considers the measurement necessary only below 10 s or for new chamber
heights (2023_Blöbaum). The chamber height of the yeast-type chip is such a new height, so the fluorescein test
is listed in 4.8. But a chamber that is perfused for 22 h with glucose-free medium, however slowly, does not keep
50 g/L glucose; the transport reading can delay the famine, not remove it.

The sensors add an independent reason for caution. All four were expressed, and both channels of every ratio
were detected in every strain (Figure 13). The ratios were computed as in the toolbox they come from: the
excitation ratios of QUEEN-2m and sfpHluorin, and the sensor channel over a constitutively expressed mCherry for
GlyRNA and OxPro, because "as both QUEEN-2m and sfpHluorin are ratiometric probes, they did not require this
addition" (2022_Pianale). But the ratios lay in narrow bands, 0.10 to 0.16 for QUEEN-2m and 0.63 to 0.96 for
sfpHluorin over all chambers of all structures (3.5), and no sensor separated constant feast from constant
famine. In *S. cerevisiae* the same QUEEN-2m reported relative ATP levels that "dropped significantly (to < 1)
when cells were exposed to substrate oscillations while remaining similar to the control (~1) in pH
oscillations" (2025_Pianale), and in *E. coli* its in vivo dynamic range is about threefold (2014_Yaginuma). A
sensor that works in a population whose ATP level does not change between 50 and 0 g/L glucose, and a sensor
that does not respond in this host, give the same picture in these data. The reserves of the first reading
would produce exactly a constant ATP level and a constant intracellular pH. The decision needs a calibration
in the chip: a pH clamp of permeabilised cells in buffers of known pH, as the toolbox did for sfpHluorin in the
plate reader (2022_Pianale) and as Miesenböck et al. did by imaging cells in buffers between pH 5.3 and 7.8
(1998_Miesenböck), and a defined depletion of ATP for QUEEN-2m, as Yaginuma et al. achieved with cyanide
(2014_Yaginuma). The validation of the methods chapter was made in the BioLector and not in the chip, at a
different cell state and without the optics of the microscope, and the affinity of QUEEN depends on temperature
(2014_Yaginuma); a calibration at 30 °C under the microscope is the one that would count.

Two properties of the sensors shape what could have been seen even with a calibration. QUEEN-2m responds
within seconds, whereas the promoter-based OxPro and GlyRNA sensors have "a response time in the order of tens
of minutes (being based on transcription, translation, and maturation of fluorescent proteins)"
(2025_Pianale). At half-cycles of 0.75 to 6 min the promoter sensors integrate over many cycles and can only
move to a new steady state, which "biosensors with response times slower than environmental oscillations are
still able to spot" (2025_Pianale). The endpoint ratio of 3.5 is therefore the appropriate readout, and a flat
endpoint ratio means that no new steady state was reached, not that nothing happened within a cycle. And the
QUEEN-2m ratio depends on pH, which is why Yaginuma et al. expressed pHluorin beside it "to clarify the effect
of diversity in intracellular pH on the QUEEN-2m signal" (2014_Yaginuma), and why Pianale et al. found "a tight
connection between ATP concentration and intracellular pH" (2022_Pianale). Under pH oscillation BSA and BSpH
have to be read together. Both ratios were flat, which is consistent with a constant intracellular pH between
external pH 3.0 and 8.0 and with the pH tolerance that the BioLector growth rates showed (3.1), and inconsistent
only with a large change in ATP.

A property of the organism may add to this. *A. pullulans* is a black yeast that "produces melanin in the
later stages of fermentation" (2024_Zhang), and in surface cultures the fluorescence of a cell-wall dye fell as
the layer matured, which the authors attribute to "reduced dye penetration into the compact, melanized matrix
rather than a reduction in biomass" (2026_Malat). Melanin absorbs more strongly at shorter wavelengths, so an
excitation ratio measured at 390 to 410 nm against 470 to 480 nm would drift with melanisation independently of
pH or ATP. Whether the chambers melanised within 22 h was not assessed. This is a hypothesis to test with the
autofluorescence of the wild type, which carries no sensor. *[Hypothesis; drop if you do not want it in the
thesis.]*

The growth of the sensor strains sets a last caveat. In oMLP with 50 g/L glucose the four biosensor strains grew
at 66 to 90 % of the wild-type rate, and their ranking against the wild type depended on the medium (3.1). In
*S. cerevisiae* the same sensors, integrated at a defined site, "did not significantly affect key physiological
parameters, such as specific growth rate and product yields" (2022_Pianale). In *A. pullulans* they did, or the
transformants differ from the wild type in more than the sensor: integration site, copy number and expression
level were not characterised, and no empty-vector transformant was available as an isogenic control. The
consequence for the chip data is limited, because every comparison of 3.4 and 3.5 is made within a strain,
oscillation chambers against the controls of the same structure. A comparison between strains in the chip,
however, is a comparison between cultures that differ in growth rate, and the strain differences in cell size
(3.2.2) should be read with this in mind.

Taken together, the robustness of growth and morphology against the period and against the oscillation itself is
established by the controls on the same structures. The robustness of the intracellular physiology is not
established, and it is not contradicted. The sensor data are a negative result about the measurement, not
about the cells, until a calibration in the chip shows the range in which the ratios can move; and the famine
control is a negative result about the readout, not about the medium, until the carbon state of the chamber has
been measured.

## 4.5 Growth mode and morphology

*[To be drafted. New passage for the product link: Campbell et al. found with a gold-conjugated pullulanase
probe that "only swollen cells and chlamydospores, and neither hyphae nor unicellular blastospores, often held
responsible for pullulan formation, appeared to produce pullulan-like material" (2004_Campbell, p.1; p.4:
"Only the multi-celled chlamydospores and swollen cells were coated with silver grains"), and that "the nitrogen
source (both organic and inorganic) and medium pH" are "particularly influential" on the life cycle
(2004_Campbell, p.1).]*

## 4.6 Measuring a blastoconidia-forming fungus in a chip

*[To be drafted; plan in the outline.]*

## 4.7 Implications for bioprocess design in a circular bioeconomy

*[To be drafted; plan in the outline.]*

## 4.8 Limitations and future work

*[To be drafted from the consolidated list in the outline: biological replication and randomisation,
fluorescein exchange test and famine verification, position within the array, 1-min frames for a subset,
in-chip sensor calibration, static replication, OD₆₀₀ calibration, bud retention and mixed models.]*
