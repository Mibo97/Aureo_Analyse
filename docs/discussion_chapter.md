# 4. Discussion

Draft of 2026-10-06/07, revised 2026-10-08: plausibility and redundancy pass, most quotations paraphrased,
section titles renamed, the conclusion folded in as 4.9; plan in `docs/discussion_outline.md`. Cross-references
use the thesis numbering of the results chapter (Tables 1 to 7, Figures 2 to 17, Sections 3.x) and of the methods
chapter (2.x). Citation keys are those of your reference list; the page and the full passage behind every
quotation and paraphrase are in `docs/literature_4.1_4.4.md`, and the quotations that remain verbatim are listed
in the notes at the end of this file. Sentences in *[italic square brackets]* are notes to the author, not
thesis text.

This chapter answers the question of the introduction, how robust *A. pullulans* is under glucose and pH
oscillations and what follows for bioprocess design, in the order in which the results allow it: first what the
period did and did not do and what carried the trends (4.1), then the one effect of the oscillation itself
(4.2) and the robustness metrics (4.3), then the two caveats that limit the answer, the controls that did not
separate and the sensors whose range in the chip is unknown (4.4). Morphology and growth mode (4.5), the
method (4.6) and the implications for bioprocess design (4.7) follow; the limitations are collected in 4.8 and
the conclusion is drawn in 4.9.

## 4.1 Robustness against the oscillation period and the structure as the experimental unit

Across half-cycle periods of 0.75 to 24 min, that is full feast/famine cycles of 1.5 to 48 min, no growth,
morphology or sensor readout of the five strains changed with the period once the constant-medium controls of
the same structures were taken into account. Of the 58 readout-by-series combinations, 31 showed no monotone
trend of the oscillation chambers, 23 a trend that a control of the same structures shared, 2 a trend that
vanished after the controls were subtracted and 2 a trend that met all conditions of a period effect (Table 7,
Figure 15). Shuffling the periods among the structures of each series produced 3.5 period effects on average,
so the two observed ones lie below chance level, and neither recurred in the other oscillation type of the same
strain. The 23 structure effects, by contrast, exceeded the 13.5 that shuffled periods produced: the control
chambers trended with the real period order more often than with a random one. The result is as negative as it
is clear. Whatever the period did to the cells over the 20 h of oscillation was smaller than what the structure
did.

Trivellin et al. define robustness "for a specific function (or phenotype) and set of perturbations" and
separate it from tolerance, which "specifically describes stable growth or survival to various perturbations"
(2022_Trivellin). In these terms the perturbation space of this thesis is the set of switching intervals, and
the functions are cell size and shape, budding, area growth and the sensor ratios; every one of them was stable
across the perturbation space. The BioLector cultivations of 3.1, in which the strains grew at pH 3.0 and
pH 8.0 and with and without glucose, describe tolerance in the sense of Olsson et al., who restrict that term to
viability and growth rate and reserve robustness for the stability of growth and production across conditions,
assessed with more than one measure (2022_Olsson). The chip data describe robustness.

What a structure effect is, physically, follows from how the experiments were made. A structure carried one
period, one loading of cells from one preculture on one day, one position on the physical chip and one flow
path. Its oscillation chambers and its constant-medium chambers share all of this; the period they do not
share. When the feast and famine controls of a series rose or fell with the period as the oscillation chambers
did, the common factor was the structure, and the period order was aligned with it: the periods were run in day
blocks (Spearman of period against date ±0.87 to ±0.89), and the shortest period usually occupied the same
structure position (3.6). The permutation result makes this alignment measurable: the structure drifts were
monotone in the real period order more often than in a random one.

The unit of the experiment follows from the same fact. Lazic defines the experimental unit as "the smallest
entity that can be randomly assigned to a different treatment condition" and separates biological replicates,
which are independent, from technical replicates, which are sub-samples of one unit and are not (2010_Lazic).
The entity that could be assigned a period was the structure. Its five oscillation chambers are sub-samples of
one culture and contribute one independent observation about the period, in the way that repeated
measurements of one person contribute one observation about a population (2010_Lazic). A test over chambers
against the period would have treated correlated observations as independent, which inflates the number of
false positives in proportion to the correlation (2010_Lazic). This is why the analysis compared one value per
structure across the periods, why the Spearman correlation with n of four to six is reported as an effect size
and not as a test, why the comparison of oscillation and constant medium was paired within structures (3.4),
and why the change within cultures that carried two or three periods was examined separately (3.4). Of the two
remedies Lazic lists for nested data, averaging the observations of a unit and random-effects models, the first
was used here; the second needs more cultures per period than this data set holds (4.8).

The published dMSCC studies share the unit. The yeast chip of Blöbaum et al. carries one positive-control
structure and five oscillation structures, one frequency per structure, and the dispersion of their data is the
standard deviation over three replicate chambers (2024_Blöbaum); Pianale et al. report each frequency over five
replicate chambers (2025_Pianale); and Täuber et al. validated the switching itself with an oscillation between
two identical media, which left the growth rate unchanged (2020_Täuber). None of this is a fault of those
studies or of this one. One period per structure is what the laminar layout of the chip allows: the arrays at
the ends of a structure lie in control zones that always receive the same medium, and only the arrays in the
switching zone, across which the laminar boundary moves, see the alternation (2023_Blöbaum). The difference
here is that every structure carried constant-medium chambers of both kinds, so that the drift of the structure
could be read off and subtracted. Without them, the 27 monotone trends of the oscillation chambers in Table 6
would have been reported as dose responses, among them the rise of the OxPro and sfpHluorin ratios with the
period in both oscillation types (3.5). The control-trend classification of 3.6 is therefore the methodological
result of this thesis that carries beyond *A. pullulans*: in a dMSCC experiment the constant-medium chambers on
the same structure are the only handle on the structure, and the period can be credited only with what exceeds
them.

Two things the result does not say. First, it does not exclude period effects smaller than the variation
between structures. The feast and famine chambers of one structure differed from those of the next structure
more than any period moved the oscillation chambers; the population robustness of the cell area differed
between the structures of a series in 20 of 20 Kruskal-Wallis tests (3.6). This variation is the detection
limit of the design. Second, the controls did not span a range. For the budding rate the bracket was degenerate
on 44 of 49 structures, so no effect could be expressed as a fraction of the feast-famine difference (4.4). The
robustness metrics carry the same two limits in their own terms, together with a sampling caveat of their own
(4.3).

Whether the cells saw the short half-cycles at all is a question these data cannot answer, but the chip
literature can. In the *C. glutamicum* chip of the same family, a fluorescein switch reached 95 % of the new
medium in the chamber within 5 s, so that switches of 5 s and longer exchange the medium almost completely
(2020_Täuber), and the protocol regards a measurement of the exchange as necessary only for oscillations faster
than 10 s or for new chamber heights (2023_Blöbaum). The shortest half-cycle used here, 45 s, lies an order of
magnitude above that limit, so the cells received the full amplitude at every interval. The chamber height of
the yeast-type chip and the passive supply of the control arrays (2.4.3) differ from the measured
configuration, which is why the fluorescein test stays in the outlook (4.8).

Against this background the flat response of *A. pullulans* stands out. In *S. cerevisiae* exposed to the
same intervals, the specific growth rate fell and the intracellular ATP level rose with longer intervals, and
the cells under 48-min cycles had the highest ATP content but the least stable and most heterogeneous one
(2024_Blöbaum). In *C. glutamicum* the colony growth rate rose from about 0.3 h⁻¹ at intervals of 5 min to 1 h
to 0.5 to 0.6 h⁻¹ at intervals of 1 min and shorter, and the mean cell length changed with the interval
(2020_Täuber). Both organisms responded to how often the medium switched. *A. pullulans*, over 20 h and
against its own controls, did not. Section 4.2 shows that it did respond to whether the medium switched at all.
The design consequences, random assignment of periods to structures and days, two cultures per period,
constant-medium chambers on every structure and a control for the chamber position, are collected in 4.8.

## 4.2 The response to medium alternation: glucose and pH oscillation

The period did not matter, but the oscillation did. Paired over the 29 structures of the Glc series, cells in
the oscillation chambers budded more often than cells in the constant-medium chambers of the same structure: the
budding rate lay above the control mean on 27 of 29 structures (median ratio 1.25) and µ_bud on 23 of 29
(1.15); the cells were larger (endpoint area above the control mean on 25 of 29 structures, ratio 1.14),
marginally rounder (eccentricity 0.98) and grew more slowly in area (µ_area below the control mean on 22 of 29
structures, ratio 0.88) (3.4). The shift appeared in every strain, with median budding ratios from 1.07 for BSG
to 1.37 for BSA, and it did not depend on the period: the ratio of oscillation to control correlated with the
period at |ρ| ≥ 0.6 in two to five of the ten series per readout, with both signs. Under pH oscillation none of
the readouts differed from the controls (ratios 0.97 to 1.09, p 0.09 to 0.65). The glucose result is the one
positive finding of the oscillation experiments. It is a finding about the alternation of the medium, not about
its frequency, and three readings of it have to be weighed.

The first is that the oscillation chambers saw a medium of intermediate strength. Blöbaum et al. note for the
same kind of design that the total time spent under feast and under starvation is the same at every frequency
(2024_Blöbaum); here every oscillation chamber spent half of its 20 h in 50 g/L and half in 0 g/L glucose,
whatever the interval, which is a time-average of 25 g/L. A medium of 25 g/L would place the oscillation
chambers between the controls. They did not lie between them: on 26 of the 49 structures the oscillation
chambers budded more than both controls and on 4 less than both (3.4), and the two controls themselves barely
differed in budding (4.4), so there was no gradient between feast and famine on which an average could lie. The
time-average reading explains the independence from the period, but not the direction.

The second reading is that the alternation itself is the stimulus, and that it acts on both growth modes of the
organism at once. In *A. pullulans* a bud is not a daughter cell of the same kind but a blastoconidium released
from a mother more than twice the size of the median cell (3.2.2), and the two ways of growing, in size and by
releasing blastoconidia, traded against each other wherever the medium differed: complex medium gave cells 3.5
times larger that budded a sixth as often as in minimal medium (3.3), and within the Glc series µ_bud and
µ_area were negatively correlated (Figure 12). Cells in alternating medium did both more than cells in constant
medium: they were larger than the control mean and budded more, while growing in area more slowly. Rensink et
al. describe the cell types of *A. pullulans* as a cycle in which swollen cells, hyphae and existing
blastoconidia release new blastoconidia, which in turn differentiate into swollen cells (2026_Rensink); a
medium that alternates between feeding and withdrawal could drive this cycle faster than either constant
medium, the feast phases feeding the swelling and the famine phases the release. This is a hypothesis. It
predicts that the shift depends on the amplitude of the alternation and not on its frequency, which is what was
observed, and it predicts that an oscillation between 50 and 20 g/L would shift less than one between 50 and
0 g/L, which was not tested. The sign of the response separates *A. pullulans* from yeast. In *S. cerevisiae*
substrate oscillation reduced the budding ratio and the specific growth rate, most at 24-min intervals in one
study (2024_Blöbaum) and at 1.5 and 6 min or at 6 and 24 min, depending on the strain, in the other, in which
substrate oscillation also gave bigger and rounder cells than pH oscillation (2025_Pianale). For yeast a bud is
growth, and the oscillation cost growth. For *A. pullulans* the oscillation slowed the area growth of the
individual cell as well, but the population ended with larger cells and more blastoconidia: the perturbation
was converted into propagules.

The third reading is position. The oscillation chambers occupied the positions A3 to A12 of an array and the
constant-medium chambers the positions A1, A2, A13 and A14 (2.5.3), because the laminar layout keeps the ends of
the array in the control zones (2023_Blöbaum). The comparison of oscillation and control chambers is therefore
also a comparison between the centre and the ends of an array, and Täuber et al. report that the exchange rate
of fluorescein differs slightly with the position of a chamber in the switching zone (2020_Täuber). One
observation argues against a position effect on the flow. The immigration rate, the number of cells washed into
a chamber per cell-hour, did not differ between oscillation and control chambers in either oscillation type
(ratios 1.02 and 1.10, p 0.90 and 0.45; 3.4), so the flow that delivers cells, and with them the medium,
reached the centre and the ends alike. What differs beyond the flow is the pressure regime: the oscillation
chambers saw the switching between the two inlets at 90 mbar while the control arrays were supplied at a
passive 10 mbar (2.4.3). In the *C. glutamicum* chip a switch between two identical media at 10-s intervals
left the growth rate unchanged (2020_Täuber), which makes a pressure effect unlikely but does not exclude it
for this chip and this organism. The control that would settle both readings, constant medium in central
chambers or an oscillation between two identical media, was not run (4.8).

Under pH oscillation nothing moved, although the external pH alternated between 3.0 and 8.0 every 0.75 to
24 min. In the BioLector the wild type grew at nearly the same rate at pH 3.0 and pH 8.0 (0.36 and 0.34 h⁻¹;
3.1), the genome of *A. pullulans* carries the expansions of the HOG pathway and of alkali-metal cation
transporters that underlie its tolerance of pH extremes (2014_Gostinčar), and in the pH series glucose was
present at 20 g/L throughout the cycle. Pianale et al. draw the same conclusion for yeast: its known resistance
to low pH explains why no pH oscillation impaired growth, and in their pH-dynamic setup glucose was never
limiting, so that homeostasis and growth were maintained (2025_Pianale). The pH sensor ratio of BSpH did not
move either, which under this reading is homeostasis and under the alternative of 4.4 an unresponsive sensor.
The two oscillation types also differ in more than the oscillating variable: the Glc structures ran on oMLP at
pH 5.5 with 50 or 0 g/L glucose and budded more and grew less in area than the pH structures on 20 g/L at pH 3.0
or 8.0 (medians 0.25 against 0.16 births per cell-hour and 0.11 against 0.18 h⁻¹; Figure 12). This difference
is one between media, structures and days as well, and it is not attributed to the famine phases alone.

The oscillation result therefore meets the test that the period results failed. It is present in every strain,
independent of the interval and measured against controls on the same structure; its sign is the opposite of
what the same experiment gave in yeast; and its meaning for the process, a shift toward blastoconidia, the
morphotype that Campbell et al. found not to produce pullulan (2004_Campbell), belongs to 4.5 and 4.7.

## 4.3 Temporal and population robustness under oscillation

The robustness metrics of Trivellin et al. were computed as Blöbaum et al. defined them for the chip: R(p) from
the standard deviation and mean of a function over all cells of a chamber at one time point, as a measure of how
homogeneous the population is, and R(t) from the dispersion of a function over time at the population and the
single-cell level (2024_Blöbaum), each normalised by the mean of the function over the data set
(2022_Trivellin), each per chamber and then through the same structure logic as the readouts (2.5.5). Against
the period they behaved like the readouts: of the 90 combinations of nine metrics and ten series, 51 showed no
trend of the oscillation chambers, 23 a trend shared by a control of the same structures, 8 a trend that
vanished after the controls were subtracted and 8 met the conditions of a period effect (3.4). Five of the
eight fell in one series, BSA/Glc, in which the population became more homogeneous in area, growth rate and
budding and more heterogeneous in eccentricity with longer periods. One series carrying five of eight period
effects is a series property. Whether it is a response of this strain or a drift of these six structures that
the controls did not share is a question for a second culture of BSA under glucose oscillation.

Against the controls of the same structure the metrics did differ, all of them in the Glc series. Under glucose
oscillation the mean cell area of a chamber was less stable over time (population R(t) −0.127 against −0.096
for the feast and −0.123 for the famine controls; below the control mean on 31 of 49 structures, p 0.013), the
population more heterogeneous in cell area (R(p) −0.671 against −0.571 and −0.613; 34 of 49, p 0.003) and in
eccentricity (p 0.003), but more homogeneous in its area growth rate (R(p) of µ_area −0.333 against −0.412 and
−0.453; above the control mean on 33 of 49, p 0.008); the budding output of the mothers was as heterogeneous as
under constant medium (p 0.95), and the budding rhythm of the individual mother was less stable under glucose
oscillation (p 0.043) and not under pH oscillation (p 0.28) (Table 4). The metrics echo 4.2 rather than add an
independent finding: cells that bud more and grow more slowly in area present a population with more small new
objects and more mothers in transition, which widens the size distribution and moves the chamber mean, and a
growth rate that is reduced is also compressed.

How much of this is the metric? R is the negative Fano factor over a global mean, and Trivellin et al. note
that the mean weighs more heavily on R than on the coefficient of variation (2022_Trivellin). A function whose
mean rises at constant relative spread therefore becomes less robust by R, and one whose mean falls becomes
more robust. Both happened here: the cells under glucose oscillation were larger and grew more slowly. Measured
as the coefficient of variation instead, four of the five differences remained and one did not (3.4): the area
growth rate of the oscillation chambers was more homogeneous (ratio to the control mean 0.84, p 0.002; Glc
0.79, p < 0.001), the eccentricity more heterogeneous (1.07, p 0.001), the chamber mean area (1.05, p 0.042)
and the chamber mean eccentricity (1.23, p 0.005) less stable over time, whereas the cell-to-cell spread of the
area did not differ (1.04, p 0.13). The greater heterogeneity of cell size under oscillation in Table 4 is
therefore the larger mean size of 4.2 expressed as a Fano factor, not a wider relative distribution; the more
homogeneous growth rate is real, and so is the less stable mean area. Trade-offs of this kind are what the
metric was built to show: robustness is specific to the function, and one perturbation can raise it for one
function and lower it for another (2022_Trivellin). Here the function that became more robust under
oscillation, the area growth rate, is the one whose mean fell, and the function that became less robust, the
chamber mean area, is the one whose mean rose. The heterogeneity of cell area and that of the growth rate were
unrelated across chambers (ρ 0.06, n = 497; 3.4), so the two R(p) values measure two properties of a population
rather than one degree of order.

Two caveats of the metric limit the reading further. R is a relative quantity whose value changes when
replicates or conditions are added (2024_Blöbaum), and its normalisation by the mean of the function ties it to
the strains and conditions of the data set (2022_Trivellin); the values of Table 4 therefore compare within
this data set, oscillation chambers against the controls of their structure, and not with the yeast values of
the literature. And R(t) "is not able to differentiate between oscillating and steadily changing functions"
(2024_Blöbaum): at 10 min per frame the half-cycles of 0.75 to 6 min are not resolved, the 12-min half-cycle is
at the limit and the 24-min half-cycle is sampled 4.8 times per cycle (3.4), so the oscillation contributes to
the temporal variance of a chamber differently at every interval, as noise at the short intervals and as a
sampled cycle at the long ones. A comparison of R(t) between periods would have compared sampling artefacts.
The comparison of R(t) against the controls of the same structure, which carry no cycle, is the one that is
interpretable, and even it includes whatever part of the cycle the frames caught.

The yeast studies with the same chip and metric show what the metrics can do when the interval matters. In
*S. cerevisiae* slower feast/famine cycles made the intracellular ATP level more heterogeneous, visible as a
fall of R(p) with longer intervals, and the cells under 48-min cycles combined the highest ATP content with the
lowest temporal stability and the highest heterogeneity (2024_Blöbaum); across three strains, substrate
oscillation produced less population heterogeneity than pH oscillation, and under substrate oscillation the
industrial strain Ethanol Red combined the highest ATP level with the highest temporal stability of that level
(2025_Pianale). For *A. pullulans* the metrics show heterogeneity and stability differences between alternating
and constant medium and none along the interval axis, and the sensor ratios showed no difference in either
(3.5, Table 5). The metrics thus confirm the two answers of 4.1 and 4.2 at the level of distributions: robust
to the frequency, responsive to the alternation. What they add on their own is the growth-rate result, a
population that grows more uniformly in area when the medium alternates, which the mean values alone would not
have shown.

## 4.4 Separation of the constant-medium controls and the dynamic range of the biosensors

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

The first reading is reserves. The cells entered the chip from a preculture in oMLP with 50 g/L glucose, a
medium derived from the liamocin production medium of Haala et al.; with 50 g/L glucose and 0.46 g/L ammonium
nitrate it has the carbon-to-nitrogen ratio of 125 of that production medium (2024_Haala). Under such a ratio,
carbon in excess and nitrogen scarce, *A. pullulans* forms intracellular storage lipids, which in oleaginous
*Aureobasidium* strains reach up to 65 % of the cell dry weight (2024_Haala, citing Xue et al. 2018); the
strain of that study is the wild type used here, NRRL 62042. An inverse relation between intracellular
glycogen and extracellular pullulan has been described as well (2004_Campbell). Cells that arrive with such
reserves can bud for many hours without external carbon, and blastoconidia formation is in any case the
response that *A. pullulans* showed whenever growth in size stopped: in the static chambers, complex medium
produced large cells and few blastoconidia and minimal medium small cells and many (3.3), and under glucose
oscillation the cells budded more and grew in area more slowly than under constant feast (3.4). Under this
reading the famine medium did reach the cells, but the readout that was expected to fall, the budding rate over
the sparse phase, is the one readout that a starving *A. pullulans* population keeps up, whereas cell size,
which needs new material, responded in the expected direction. The oscillation results of 3.4 then stand, and
the missing separation of the controls is a property of the organism and of the 22-h window, not of the
measurement.

The second reading is that the famine medium was not famine. The medium without glucose still contained the
salts, the nitrogen source, the vitamins and the trace elements (2.1), and the preculture medium entered the
chip with the cells. Pianale et al. changed their own starvation condition from 0 g/L to 10 mg/L glucose
between two studies and point out that starving cells and calorie-restricted cells, which have substrate for
metabolic activity but not for reproduction, differ in the expression of more than 500 genes (2025_Pianale):
what famine means for a cell depends on residual carbon in a range that the chip experiments did not measure.
A glucose indicator in the chamber or a sampled outlet would settle this.

The third reading is transport. The control arrays were supplied at a passive 10 mbar against 90 mbar on the
oscillation line (2.4.3), and the exchange time of the chambers was not measured here. The published exchange
time of 5 s (4.1) argues against a transport limitation, and a chamber that is perfused for 22 h with
glucose-free medium, however slowly, does not keep 50 g/L glucose; the transport reading can delay the famine,
not remove it.

The sensors add an independent reason for caution. All four were expressed, and both channels of every ratio
were detected in every strain (Figure 13). The ratios were computed as in the toolbox they come from: the
excitation ratios of QUEEN-2m and sfpHluorin, which are ratiometric by themselves, and the sensor channel over a
constitutively expressed mCherry for GlyRNA and OxPro (2022_Pianale). But the ratios lay in narrow bands, 0.10
to 0.16 for QUEEN-2m and 0.63 to 0.96 for sfpHluorin over all chambers of all structures (3.5), and no sensor
separated constant feast from constant famine. In *S. cerevisiae* the same QUEEN-2m reported relative ATP
levels that fell clearly below the control under substrate oscillation and stayed at the control level under
pH oscillation (2025_Pianale), and in *E. coli* its in vivo dynamic range is about threefold (2014_Yaginuma).
A sensor that works in a population whose ATP level does not change between 50 and 0 g/L glucose, and a sensor
that does not respond in this host, give the same picture in these data. The reserves of the first reading
would produce exactly a constant ATP level and a constant intracellular pH. The decision needs a calibration
in the chip: a pH clamp of permeabilised cells in buffers of known pH, as the toolbox did for sfpHluorin in the
plate reader (2022_Pianale) and as Miesenböck et al. did by imaging cells in buffers between pH 5.3 and 7.8
(1998_Miesenböck), and a defined depletion of ATP for QUEEN-2m, as Yaginuma et al. achieved with cyanide
(2014_Yaginuma). The validation of the methods chapter was made in the BioLector and not in the chip, at a
different cell state and without the optics of the microscope, and the affinity of QUEEN depends on temperature
(2014_Yaginuma); a calibration at 30 °C under the microscope is the one that would count.

Two properties of the sensors shape what could have been seen even with a calibration. QUEEN-2m responds
within seconds, whereas the promoter-based OxPro and GlyRNA sensors respond within tens of minutes, because
transcription, translation and the maturation of the fluorescent protein lie between stimulus and signal
(2025_Pianale). At half-cycles of 0.75 to 6 min the promoter sensors integrate over many cycles and can only
move to a new steady state, which a sensor slower than the oscillation still reports (2025_Pianale). The
endpoint ratio of 3.5 is therefore the appropriate readout, and a flat endpoint ratio means that no new steady
state was reached, not that nothing happened within a cycle. And the QUEEN-2m ratio depends on pH, which is why
Yaginuma et al. expressed pHluorin beside it to account for differences in intracellular pH (2014_Yaginuma),
and why Pianale et al. found ATP concentration and intracellular pH tightly connected (2022_Pianale). Under pH
oscillation BSA and BSpH have to be read together. Both ratios were flat, which is consistent with a constant
intracellular pH between external pH 3.0 and 8.0 and with the pH tolerance that the BioLector growth rates
showed (3.1), and inconsistent only with a large change in ATP.

The growth of the sensor strains sets a last caveat. In oMLP with 50 g/L glucose the four biosensor strains grew
at 66 to 90 % of the wild-type rate, and their ranking against the wild type depended on the medium (3.1). In
*S. cerevisiae* the same sensors, integrated at a defined site, left the specific growth rate and the product
yields unchanged (2022_Pianale). In *A. pullulans* they did not, or the transformants differ from the wild type
in more than the sensor: integration site, copy number and expression level were not characterised, and no
empty-vector transformant was available as an isogenic control. The consequence for the chip data is limited,
because every comparison of 3.4 and 3.5 is made within a strain, oscillation chambers against the controls of
the same structure. A comparison between strains in the chip, however, is a comparison between cultures that
differ in growth rate, and the strain differences in cell size (3.2.2) should be read with this in mind.

Taken together, the robustness of growth and morphology against the period and against the oscillation itself is
established by the controls on the same structures. The robustness of the intracellular physiology is not
established, and it is not contradicted. The sensor data are a negative result about the measurement, not
about the cells, until a calibration in the chip shows the range in which the ratios can move; and the famine
control is a negative result about the readout, not about the medium, until the carbon state of the chamber has
been measured.

## 4.5 Growth modes and morphotypes in the microfluidic chambers

The cells in the chambers were, for the most part, the ovoid yeast-like cells of the literature. Over the
22,090 cells tracked for at least ten frames the median cell measured 22 µm² with an axis ratio of 1.5, and 90 %
of the cells had an axis ratio below 1.9 (3.2.2). Rensink et al. compile the published dimensions: yeast cells
of 9 to 11 by 3 to 6.5 µm, globular and ellipsoid swollen cells of 15 by 11 and 12 by 9 µm, and chlamydospores
of 13 by 12 µm (2026_Rensink). As projected areas these are about 20 to 55 µm² for the yeast cells, 85 to
130 µm² for the swollen cells and about 120 µm² for the chlamydospores. The median cell of this study therefore
sits in the yeast-like class, and the 30 µm² and eccentricity 0.6 that mark the "large and round" quadrant of
Figure 4 lie below the published swollen-cell size, so that quadrant collects cells on the way to swelling as
well as swollen cells. The mothers that released blastoconidia measured 71 µm² at the moment of release, the
area of a 9.5-µm sphere, between the two classes; the two cells of about 270 µm² on chip W65 (Figure 5 B)
correspond to 18.5 µm and exceed the published swollen cells. Chlamydospores and germ tubes were not seen,
pseudohyphal chains only in complex medium on W65, hyphae only on W109 in complex medium (3.2.2).

The organism grew in two ways, and the two traded against each other. In complex medium on chip W109 the cells
ended 3.5 times larger and budded a sixth as often as in minimal medium on the same chip; the population grew
into cell size rather than into blastoconidia (3.3). Within the Glc series, structures with a higher µ_bud had a
lower µ_area (Figure 12), and under glucose oscillation the cells budded more and grew in area more slowly than
under constant medium (4.2). Rensink et al. describe the same two modes as stages of one cycle, in which
blastoconidia differentiate into swollen cells and swollen cells into hyphae while all three release new
blastoconidia; in their static liquid cultures blastoconidia had ceased to be the dominant cell type after
10 h, and their model predicts 57 % blastoconidia, 32 % swollen cells and 11 % hyphae after 72 h (2026_Rensink).
The chambers show the medium deciding where on this cycle a population settles within 22 h: minimal medium
kept it at yeast-like cells that release blastoconidia from a minority of large mothers, complex medium moved
it toward swollen cells and, on W109, hyphae. Campbell et al. name the nitrogen source and the medium pH as the
most influential factors of the life cycle (2004_Campbell), and the two media differ in both: YPD supplies
peptone and yeast extract, oMLP ammonium nitrate at a carbon-to-nitrogen ratio of 125 (2.1; 2024_Haala).
Whether the large round cells of the minimal-medium chambers were the lipid-storing cells that this ratio
induces (4.4) was not examined. Pseudohyphal growth, finally, had a different trigger here than in yeast:
Blöbaum et al. saw it in *S. cerevisiae* under glucose shifts every 1.5 and 6 min and attribute it to nitrogen
starvation and stress (2024_Blöbaum); in *A. pullulans* it appeared with the complex medium and not with the
oscillation.

The two static chips disagreed, and the disagreement is instructive. W65 repeated the cell size of complex
medium but not the small cells and frequent budding of minimal medium, and one of its four minimal-medium
chambers ended at 215 µm² against 70 to 86 µm² for the other three, because it contained the two swollen cells
that released a ring of ten blastoconidia within three frames (3.3, Figure 5 B). This is the transition from
swollen cell to blastoconidia that Rensink et al. describe, caught in one chamber, and with four chambers per
medium one such cell decides the mean. The lesson is the one of 4.1 at a smaller scale: one culture per chip,
and chamber means that a single cell can move. Per-cell distributions (Figure 4) carry this information where
chamber means hide it, and more chips per medium are the remedy.

For the process the morphology matters because the morphotype decides the product. With a gold-conjugated
pullulanase probe Campbell et al. found that "only swollen cells and chlamydospores, and neither hyphae nor
unicellular blastospores, often held responsible for pullulan formation, appeared to produce pullulan-like
material" (2004_Campbell). Read with this, the glucose-oscillation shift of 4.2 is two-edged: it increases the
release of blastoconidia, which do not produce pullulan, and it enlarges the mothers, which may. Which side
wins for the titre cannot be read from morphology and needs a product measurement. The medium adds a second
caveat: oMLP was designed for liamocins (2024_Haala), whereas pullulan production is reported to favour a
carbon-to-nitrogen ratio of 10:1 (2024_Zhang), so the morphology distribution seen here is that of a liamocin
medium and not of a pullulan process. Both points return in 4.7.

## 4.6 Single-cell analysis of a blastoconidia-forming fungus

The introduction named two problems that single-cell microfluidics meets in *A. pullulans* and not in
*S. cerevisiae*: the assignment of buds to mothers when the buds are blastoconidia that move, and the secretion
of pullulan into a device that depends on micrometre-scale flow (1.4). Both have an answer.

The mother-bud assignment works in the sparse phase and only there. A new object is a bud candidate; its mother
is the cell whose mask it touched at first detection, which decided 74 % of the accepted events, or the nearest
established cell within 30 px for the 26 % that touched nothing; and a size criterion removes the washed-in
cells and split masks, because the ratio of candidate area to mother area is bimodal, with buds at 0.07 and
mother-sized objects at 0.45 and an antimode at 0.32 that rejected 15.5 % of the candidates (3.2.5, Figure 8).
The restriction to the sparse phase follows from Figure 7: above about 8 objects per frame every object gave
rise to a new touching object at a constant rate of 0.04 to 0.06 per frame, whether the chamber held 10 or 200
objects, so at high density the touching new objects are fragments of touching masks that the contact rule
cannot tell from buds (3.2.4). The sparse-phase window covered a median of 117 of 133 frames, and 245 of the 564
chambers never left it, so the restriction cost little data; it cost the dense phase, in which lineage
readouts would have needed a different method. That the problem is real and not an artefact of one tracker is
shown by 3.2.6: the budding-rate trends against the period changed when the tracker was rebuilt, the per-series
correlations of the two trackings agreeing at ρ 0.05 while the budding rates of the individual structures agreed
at 0.77. Per-series trends of a lineage readout are a property of the tracker until the controls on the same
structures confirm them, which is a second reason for the control logic of 4.1.

Pullulan is not what makes the control chambers of a structure disagree. On the knockout structure the
chamber-to-chamber variation of the cell-area level among the famine controls was 0.23 and among the feast
controls 0.37, the 67th and 96th percentile of the 49 producer structures; the temporal variation lay at the
92nd and 61st percentile; the famine-control cells were 1.5 times larger than the feast-control cells where the
producers had a median ratio of 0.98, and the famine controls budded faster than the feast controls (3.7, Figure
16). A single structure cannot show that pullulan has no effect on the flow, but it shows that removing it did
not make nominally identical chambers agree better. The remaining candidates for the disagreement, seeding
density, position within the array and flow asymmetry between arrays, are those of 4.2 and 4.8.

The segmentation was a generalist model used out of the box. Cellpose-SAM was built to generalise beyond its
training distribution, is robust to changes of channel order, cell size, noise, resolution and blur, and is
meant to run without adaptation on images acquired at other pixel sizes and under other degradations
(2025_Pachitariu). Phase-contrast images of a black yeast with cells from 8 to 270 µm² and pseudohyphal chains
are such out-of-distribution data, and the model segmented them without training on *A. pullulans*. The price
was paid in two places: in the dense phase, where touching masks merge and split (3.2.4), and in the chains of
complex medium, which were split into their single cells (3.2.2). The fine-tuning and human-in-the-loop training
that the framework offers (2025_Pachitariu) were not used here and are the obvious next step for the dense phase
and for the swollen and pseudohyphal forms. The cell filter replaced the manual curation of the earlier pipeline
with two rules, a size and persistence rule that removed 35 % of the tracks but 6.8 % of the object-frames, and
a phase-contrast rule for dead cells and debris that removed 432 tracks (3.2.1); the structure BSG/pH/6 min,
which lost 52 % of its object-frames to them, shows that the rules act where the images are bad and not
uniformly.

The tracker halved the fragmentation of its predecessor, from 0.14 to 0.070 new tracks per object-frame, and 78
% of the object-frames lay in tracks of at least ten frames, but the median track still lasted only 5 frames
(3.2.3). Together with the open chamber design, this sets what the platform can and cannot measure for this
organism. Blastoconidia leave the chamber: cells entered with the flow at 0.19 per cell-hour against 0.21 births
per cell-hour, and immigration exceeded births in 37 % of the chambers (3.4), so the object count of a chamber is
not a growth curve, µ_bud is a lower bound of the birth rate, and lineages rarely reach a second generation. The
same wash-out is what kept the chambers sparse for 22 h: where the monolayer chambers of the yeast chips hold
150 to 1,000 cells (2024_Blöbaum), the *A. pullulans* chambers started with one to three cells and reached at
most 180 objects, 245 of them never more than 20. What remains measurable are cumulative, hour-scale
quantities, the endpoint morphology, the sparse-phase budding rate, µ_bud, µ_area and the robustness metrics,
and these are the quantities the thesis reports. The cycles themselves were not observable at 10 min per frame
for intervals up to 6 min and barely for 12 min (3.4, 4.8).

## 4.7 Implications for scale-down studies and bioprocess design

The half-cycle periods tested here span the time scales of large-scale mixing. Nadal-Rey et al. give the mixing
times of bioreactors as 10 to 1,000 s depending on design and operation, with circulation times three to five
times shorter (2021_Nadal-Rey); the industrial vessels modelled by Losoi et al. had measured mixing times of
124 and 165 s, which optimally placed multiple feed points reduced in the simulation from minutes to about 10 s
(2022_Losoi); Arulrajah et al. put laboratory mixing times below 5 s and large-scale mixing times at tens to
hundreds of seconds, longer than the reaction times of the cell, which reach seconds at the level of the
transcriptome (2025_Arulrajah). A cell in such a reactor passes through zones of excess and limitation on the
scale of the circulation time, seconds to a few minutes, and sees the whole vessel on the scale of the mixing
time. The half-cycles of 0.75 to 24 min, chosen with the yeast studies to reach the mixing times of large
reactors (2025_Pianale), cover this range and extend beyond it. Within it, *A. pullulans* did not care how often
the medium switched (4.1). The cellular responses to such fluctuations span seconds to hours (2021_Nadal-Rey),
and the hour-scale readouts of this thesis integrate over all of them; whatever the cells did within a cycle, it
did not accumulate into a difference between cycles of 1.5 and of 48 min.

The amplitude, not the frequency, is the parameter to carry into process design. Two points of the results say
so. First, the one effect of the oscillation was the alternation between glucose excess and glucose absence,
present at every interval (4.2). Second, the amplitude used here is a worst case. In a 20 m³ fed-batch process
with *S. cerevisiae* the local glucose peaks measured 40 to 80 mg/L (2021_Nadal-Rey); the chambers alternated
between 50 and 0 g/L, almost three orders of magnitude more, and in a scale-down reactor for baker's yeast,
sugar pulses of only 0.45 to 1.9 g/L for 60 s on a tenth of the culture already cost 6 to 7 % of the biomass
yield (2025_Arulrajah). The chip result therefore says that even the extreme alternation did not change growth
or morphology with the frequency, and that it shifted the population toward blastoconidia and larger mothers.
For the process this shift is the two-edged one of 4.5. Fluctuating conditions raise the ATP demand for
maintenance, by 40 to 50 % in *E. coli*, and leave less carbon for the product (2021_Nadal-Rey); a population
that answers glucose alternation with more blastoconidia spends carbon on propagules that do not produce
pullulan, while the larger mothers may. Which of the two prevails is a question for a product measurement
under alternation, with a pullulan or liamocin assay, and it is the first experiment that follows from this
thesis. That the answer is specific to the organism is the general lesson of scale-down work, in which the
impact of gradients is strongly organism-dependent (2025_Arulrajah); for *A. pullulans* this thesis places the
sensitivity in the amplitude of the glucose alternation and the composition of the medium, not in the
frequency.

pH gradients are the smaller concern. Losoi et al. find the simulated pH gradients after a base pulse similar
in extent to the substrate gradients (2022_Losoi); in the chip, alternation between pH 3.0 and pH 8.0 left
growth, morphology and the pH sensor unchanged (4.2, 4.4). pH gradients in reactors are more localised and
shorter-lived than substrate or oxygen gradients (2025_Arulrajah), so the chambers' 22 h of alternation between
pH 3.0 and 8.0 exceed what a cell meets in a reactor, and unlike base pulses in a scale-down reactor, which
raise the pH and the osmolarity of the whole vessel at once (2025_Arulrajah), the chip separated the pH from the
medium. For a circular feedstock that arrives acidic or alkaline, or for a reactor in which base addition
creates local pH excursions, the organism's growth is not the limiting factor; the product may be, since the
life cycle and pullulan formation are pH-sensitive (2004_Campbell), which is again a product question.

Feedstock composition mattered more than feedstock dynamics. The largest morphological effect of the whole data
set was between complex and minimal medium in the static chambers (3.3), not between any two oscillation
regimes, and the ranking of the strains in the BioLector changed with the medium (3.1). The requirement of the
introduction, a cell factory that "maintain[s] performance while conditions fluctuate" (1.1), is met by
*A. pullulans* for the fluctuation frequencies and the pH range tested; what a process on variable feedstocks
has to control is the composition, in particular the carbon-to-nitrogen ratio that Haala et al. optimised for
liamocins (2024_Haala) and that Zhang et al. report at 10:1 for pullulan (2024_Zhang).

Two further implications concern strains and methods. The biosensor strains grew more slowly than the wild type
and their ranking depended on the medium (4.4); as production hosts they carry a cost, and as monitoring strains
they report nothing until their sensors are calibrated in the process context. The wild type remains the
production candidate. On the method side, dynamic single-cell cultivation is now named among the scale-down
tools as a complement to stirred-tank scale-down studies (2025_Arulrajah), and the thesis delivers what Olsson
et al. ask for, a quantification of robustness that can guide strain and bioprocess development (2022_Olsson),
with one addition: the control chambers on the same structure, without which the dMSCC data would have shown
dose responses that were not there (4.1). A scale-down study of *A. pullulans* that wants to rank strains or
media should therefore vary amplitude and composition, replicate per culture, keep constant-medium chambers on
every structure, and add a product readout at the single-cell level; a fluorescent pullulan reporter of the
kind Zhang et al. used for screening (2024_Zhang) could provide it in the chip.

## 4.8 Limitations and outlook

The limitations of this work follow from its design, and most of them have a direct remedy.

The experimental unit. One culture per structure and one period per structure made the structure the detection
limit (4.1). Future runs should assign periods to structures and days at random, run every period on at least
two cultures, and keep feast and famine chambers on every structure. With two cultures per period the
random-effects models of 4.1 become possible, and the between-structure variation becomes a measured quantity
instead of the limit.

Position and pressure. The oscillation chambers sat in the centre of the array and the controls at its ends,
and the control arrays were supplied passively (4.2). A run with constant medium in the central chambers, or an
oscillation between two identical media as Täuber et al. used to exclude pressure effects (2020_Täuber), would
separate position and switching from the oscillation.

The chamber as a medium. The exchange time of the chambers was taken from the *C. glutamicum* chip (2020_Täuber;
2023_Blöbaum) and not measured for the chamber height used here; a fluorescein switching test at the intervals
used, 0.75 min in particular, would close this. Whether the famine medium was famine at the cell was not
measured either (4.4); a glucose indicator in the chamber or a sampled outlet would settle it.

Sampling. At 10 min per frame the cycles of the intervals up to 6 min are not observable and the 12-min
interval is at the limit (3.4). A subset of chambers imaged at 1-min intervals would show whether cell size or a
sensor ratio follows the 24-min cycle and how much of the temporal variance of R(t) is the cycle itself (4.3).

Sensors. The four sensors are expressed and detected but not calibrated in the chip (4.4). An in-chip pH clamp
of permeabilised cells, a defined ATP depletion for QUEEN-2m and a defined oxidative stimulus for OxPro at 30 °C
under the microscope would give the range in which the ratios can move; an isogenic empty-vector transformant
would separate the cost of the sensor from clonal variation in the growth rates of 3.1.

Morphology and product. The static comparison rests on two chips with one culture each (3.3, 4.5); more chips
per medium and a per-cell product readout would turn the morphology results into process statements. The
segmentation can be fine-tuned on *A. pullulans* for the dense phase and the pseudohyphal forms
(2025_Pachitariu), and bud retention in the chamber geometry would deepen the lineages that wash-out now cuts
short (4.6).

BioLector. The scattered-light signals are uncalibrated and the five wells per condition come from one
preculture (3.1); an OD₆₀₀ calibration and biological replicates are needed before the growth parameters are
compared across strains.

Multiple testing. The 58 and 90 readout-by-series combinations of 3.6 and 3.4 were classified without a
correction for multiple testing, and so were the paired comparisons of 3.4 (2.5.6). The permutation null of 3.6
states how many period and structure effects the classification yields when the period has no effect, and it
is the safeguard on which the classification rests; the paired comparisons rest on their consistency across
strains and oscillation types rather than on a correction.

## 4.9 Conclusion

This thesis asked how robust *Aureobasidium pullulans* is under glucose and pH oscillations and what follows
for bioprocess design. In 564 microfluidic chambers on 54 structures, the wild type, four biosensor strains and
a pullulan knockout were followed for 22 h, 547 of the chambers under feast/famine half-cycle periods of 0.75
to 24 min with constant-medium control chambers on every structure and 17 under constant complex or minimal
medium, and analysed with a pipeline that tracks every cell, assigns blastoconidia to their mothers in the
sparse phase and compares every readout at the level of the structure. Growth, morphology and sensor signals
did not depend on the period; the trends that appeared were carried by the structures, as their constant-medium
controls showed, and the two remaining period effects lie below chance level. The alternation itself shifted
growth toward blastoconidia and larger cells under glucose oscillation in every strain and at every interval,
and changed nothing under pH oscillation between 3.0 and 8.0; the robustness metrics confirm both at the level
of distributions. The constant-feast and constant-famine chambers barely separated in budding and the
biosensor ratios did not separate at all, so the robustness of the intracellular physiology is neither
established nor contradicted until the sensors are calibrated in the chip and the carbon state of the chamber
is measured. Pullulan secretion does not explain the disagreement of control chambers, and the mother-bud
assignment of a blastoconidia-forming fungus works where the chambers are sparse. For process design, the
frequency of glucose and pH fluctuations is not the parameter to fear in *A. pullulans*; the amplitude of
glucose alternation and the composition of the medium are, because they decide between growth in size and
growth by blastoconidia, and with it between the morphotypes that produce pullulan and those that do not.

---

## Notes to the author (not thesis text)

**Quotations that remain verbatim** and need quotation marks and a page number in the thesis (pages are PDF
pages, see `docs/literature_4.1_4.4.md`): 2022_Trivellin p.1, the definition of robustness against tolerance
(4.1); 2010_Lazic p.2, the definition of the experimental unit (4.1); 2024_Blöbaum p.16, R(t) and oscillating
against steadily changing functions (4.3); 2004_Campbell p.1, the morphotypes that produce pullulan (4.5);
and the phrase "maintain[s] performance while conditions fluctuate" of your own introduction (4.7). Every other
passage is paraphrased; the originals stand in the literature file under the claim they support.

**Checks left to you.** The 0.46 g/L ammonium nitrate of 4.4 is the amount that gives a carbon-to-nitrogen
ratio of 125 at 50 g/L glucose; confirm it against the oMLP recipe of 2.1. The statement that the shortest
period usually occupied the same structure position (4.1) rests on the author statement in 3.6, for which no
run sheet exists. 2014_Gostinčar (4.2) and 2024_Rensink are cited from your reference list; their passages are
not in the literature file (the 2026_Rensink dimensions and cell-type cycle are, p.1 to p.2).

**Changes of 2026-10-08.** Redundancies removed: the R caveats are stated once (4.3) and referenced from 4.1;
the exchange time (95 % in 5 s) is stated once (4.1) and referenced from 4.4 and 4.8; the design lessons of 4.1
were merged into 4.8; the two-edged product reading is argued in 4.5 and referenced from 4.7; the 1-min-per-frame
proposal and the sensor-strain growth rates are referenced instead of repeated in 4.6 and 4.7. The "Statistics
and pipeline" paragraph of 4.8 became a multiple-testing paragraph, since both checks are now in the pipeline
and the methods rewrite. The conclusion lost its "separate chapter" note and now names the chamber split
(547 oscillation, 17 static) and the six strains explicitly.
