# Aureo_Analyse

Auswertungs-Pipeline für Einzelzell-Mikroskopiedaten aus Mikrofluidik-Kammern.
Sie liest die `Combined_Results`-Tabellen der Cellpose-Bildverarbeitung ein und
erzeugt daraus Tabellen und Abbildungen zu Zellmorphologie, Reproduktion
(Budding), Wachstumsraten und Robustheit — getrennt für Oszillations- und
statische Bedingungen.

**Alle Abbildungen sind auf Englisch beschriftet.** Code-Kommentare und
Log-Ausgaben sind deutsch.

---

## Installation

```bash
pip install -r requirements.txt
```

`statannotations` ist optional: fehlt es, entstehen die Violin-Plots ohne
Signifikanz-Sternchen (mit Warnung beim Start).

## Konfiguration

Alle Pfade und Parameter stehen an **einer** Stelle: [`analyse_pipeline/config.py`](analyse_pipeline/config.py).
`run_analysis.py` und `inspect_lineage.py` lesen beide von dort, damit
Kalibrierung und Auswertung garantiert dieselben Schwellen benutzen.

Den Datenpfad setzt man auf eine von drei Arten (die erste gesetzte gewinnt):

```bash
# 1. Umgebungsvariable (Cluster, Skripte, CI)
export AUREO_DATA_ROOT=/prj/microfluidic/ma_mimorde/Data
export AUREO_OUTPUT_DIR=/prj/microfluidic/ma_mimorde/analysis_output   # optional

# Windows PowerShell:
#   $env:AUREO_DATA_ROOT = "D:\Studium\...\Data"
```

2. `ACTIVE_PRESET` in `config.py` auf `"local"` oder `"cluster"` setzen.
3. Den jeweiligen Pfad in `DATA_ROOT_PRESETS` direkt editieren.

## Ausführen

```bash
cd analyse_pipeline
python run_analysis.py
```

Beim Start schreibt die Pipeline die **aktive Konfiguration** ins Log, inklusive
einer Warnung für jeden Parameter, der vom dokumentierten Standard abweicht, und
für die bekannten methodischen Eigenheiten (siehe *Bekannte Einschränkungen*).
Dieser Block gehört ins Laborbuch — er sagt, mit welchen Schwellen die
Ergebnisse entstanden sind.

Parameter kalibrieren (interaktiv, gegen die QC-Overlay-TIFFs):

```bash
ipython -i inspect_lineage.py
```

## Erwartete Ordnerstruktur

```
Data/<Biosensor>/<Oszillationstyp>/<Oszillationsfrequenz>/03_results/Combined_Results.{csv,parquet}

Beispiele:
  Data/BSG/Glc/6/03_results/Combined_Results.csv
  Data/BSA/pH/24/03_results/Combined_Results.csv
  Data/BSG/static/static_ypd/03_results/Combined_Results.csv
  Data/PKO/Glc/3/03_results/Combined_Results.csv
```

Die Metadaten `biosensor`/`osc_type`/`osc_freq` kommen aus dem **Ordnerpfad**,
nicht aus dem Dateinamen (siehe `data_loading.py`).

Daraus entstehen **drei unabhängige Zweige** mit je eigenem Output-Ordner. Sie
laufen durch dieselben Schritte, teilen aber keine Zahlen:

| Zweig | Erkennungsmerkmal | Output |
| --- | --- | --- |
| Oszillation | alles übrige | `analysis_output/` |
| statisch | `osc_type == "static"` | `analysis_output/static/` |
| PKO | `biosensor == "PKO"` (`Data/PKO/…`) | `analysis_output/pko/` |

Einzeln abschaltbar mit `--skip-static` bzw. `--skip-pko`.

Benötigte Spalten in `Combined_Results`: `track_id`, `frame`, `centroid_x`,
`centroid_y`, `area`, `condition`, `replicate`, `chamber`.
Optional, aber genutzt: `eccentricity`, `solidity`, `mean_<Kanal>`, `filename`.

## Versuchseinheiten: was ist hier ein Replikat?

Die Spalten `replicate` und `chamber` kommen aus dem Dateinamen der
Bildverarbeitung und bedeuten **nicht**, was ihre Namen nahelegen.
`experiment_units.py` leitet daraus die echten Einheiten ab und schreibt
`00_chip_overview.csv` (eine Zeile pro Chip).

| Zweig | pro (Stamm, osc_type, Periode) | `replicate` ist … | biologische Einheit | n je Bedingung |
| --- | --- | --- | --- | --- |
| Oszillation, PKO | **ein Chip**, ein Tag, **eine Vorkultur**; ~11 Kammern (5 Osc, 3 PosCtrl, 3 NegCtrl) auf mehreren Arrays | ein Array-Index, der gleiche Positionslabels (`ChamA13` ×3) auseinanderhält | der **Chip** (= der Batch) | **1** |
| statisch | — | **ein eigener Chip** mit eigener Vorkultur; `chamber` ist immer `ChamA0` | der Chip | 4–5 pro Medium und Chip-Familie |

Drei Folgen, die in die Arbeit gehören:

1. **Innerhalb einer Oszillationsbedingung gibt es keine biologische
   Replikation.** Alle Kammern einer Periode sind technische Wiederholungen
   einer Kultur auf einem Chip. Jeder Fehlerbalken dort ist ein
   Kammer-Fehlerbalken; die Tabellen sagen das in `error_unit` (`chamber`
   vs. `chip`) und zählen `n_units` entsprechend. Die Dosis-Wirkung über die
   Perioden ist **ein Chip pro Dosis** (6 bei Glc, 4 bei pH); Spearman läuft
   auf Chip-Mittelwerten mit n = Zahl der Perioden und ist damit eine
   Effektstärke, kein Test, auf den man sich stützt.
2. **Die Kontrollen liegen auf demselben Chip wie die Behandlung.** Jede
   Periode kann deshalb gegen *ihre eigenen* PosCtrl/NegCtrl normiert werden
   (Bracket-Score, Schritt 13) — das entfernt den Chip-/Tages-/Kultur-Effekt,
   ohne dass man den Chip kennen müsste.
3. **Kammerposition und Bedingung sind konfundiert**, weil der Chip die
   Medien fest verdrahtet: `A1/A2` = Feast, `A13/A14` = Famine, `A3–A12` =
   Wechsel. „Die Kontrollen verhalten sich anders“ und „die Randreihen
   verhalten sich anders“ sind innerhalb des WT dieselbe Beobachtung. PKO auf
   demselben Layout ist das Einzige im Datensatz, das beides trennt.

Die Stämme (WT, BSA, BSG, BSO, BSPH) bleiben in allen Auswertungen getrennt
(je eine n = 1-Serie); das Datum im Dateinamen (`260630_…`) ist die
Laborbuch-Referenz und steht in `00_chip_overview.csv`. Die Aufnahmen umfassen
~133 Frames à 10 min, also **~22 h**.

## Module

| Datei | Aufgabe |
| --- | --- |
| `config.py` | **Alle** Pfade & Parameter, plus `log_active_configuration()` |
| `run_analysis.py` | Hauptskript, orchestriert Schritte 00–95 |
| `data_loading.py` | Einlesen der Hierarchie, `exp_id`/`cell_uid`, Parquet-Cache |
| `qc_exclusions.py` | Nicht-destruktives manuelles QC: Track-Merges & Exclusions |
| `sensors.py` | Ratiometrische Sensoren → `ratio_*`-Spalten (`SENSOR_CONFIG`) |
| `lineage.py` | Mutter/Bud-Heuristik, Budding Ratio pro Mutter |
| `growth_rate.py` | µ_event aus Budding-Intervallen (Eq. 2) |
| `area_growth.py` | µ_area aus ln(Fläche)-Fit, plus µ_event-vs-µ_area-Scatter |
| `budding_ratio_timeseries.py` | Budding Ratio als Zeitreihe (Eq. 3) |
| `robustness.py` | R(t)/R(p) nach Eq. 1 |
| `endpoint_trends.py` | **Kumulativer Endzustand gegen die Periode + Spearman-Trendtest** |
| `control_consistency.py` | Kruskal-Wallis: driften PosCtrl/NegCtrl über Batches? |
| `queen_controls.py` | PosCtrl-vs-NegCtrl-Validierung der Sensoren selbst |
| `pko_comparison.py` | **WT gegen PKO**: Kammer-Übereinstimmung, Kontroll-Bracket, Zell-Ausbeute |
| `analysis.py`, `summary_plots.py`, `violin_plots.py`, `mother_trajectories.py` | Plots & gemeinsame Helfer |
| `inspect_lineage.py` | Interaktive Kalibrierung der Lineage-Parameter |
| `validate_lineage.py` | **Quantitative** Validierung der Mutter/Bud-Heuristik (alle Kammern) |
| `plot_qc_lineage_overlay.py` | **Visuelle** Validierung: Events ins QC-TIFF zeichnen (eine Kammer) |

## Output

Die Dateinamen sind nach Auswertungsschritt nummeriert, damit die
alphabetische Sortierung im Ordner der inhaltlichen Reihenfolge entspricht:

| Präfix | Inhalt |
| --- | --- |
| `00_` | Übersicht / Sanity-Check; `00_chip_overview.csv` = eine Zeile pro **Chip** |
| `10_`–`12_` | Zellmorphologie & Wachstum (Fläche, µ_event, µ_area) |
| `13_` | **Kumulativer Endzustand gegen die Periode** (+ Spearman) |
| `20_`–`23_` | Lineage: Budding-Events, Budding Ratio, Panel A, Stammbaum |
| `30_`–`31_` | Sensor-Intensitäten und Ratios über die Zeit |
| `40_` | Robustheit R(t)/R(p) inkl. Kontroll-Konsistenz |
| `50_` | Zusammenfassungstabelle |
| `90_`–`92_` | Anhang: Morphologie-Scatter, Einzelzell- & Mutter-Trajektorien |
| `95_` | Anhang: Sensor-Controls (PosCtrl vs. NegCtrl pro Biosensor) |
| `60_`–`61_` | **Nur in `pko/`**: Produzenten-gegen-PKO-Vergleich (siehe unten) |
| `70_` | **Nur in `qc_comparison/`**: mit QC vs. ohne QC (siehe unten) |

Statische Daten und PKO-Daten landen in denselben Präfixen unter
`analysis_output/static/` bzw. `analysis_output/pko/` — beide ohne
`30_`/`31_`/`95_`, da dort keine Fluoreszenzkanäle vorliegen.

---

## Abtastung: die Oszillation ist eine Behandlung, keine Messgröße

Die Spalte `osc_freq` enthält trotz ihres Namens die **Periode in Minuten**
(0.75 … 24), keine Frequenz. Bei `MIN_PER_FRAME = 10` liegt die kürzeste
auflösbare Periode (Nyquist) bei **20 min** — fünf der sechs Bedingungen liegen
darunter, die sechste nur knapp darüber:

| Periode | Zyklen pro Frame | Zyklen in 20 h Oszillation | auflösbar? |
| --- | --- | --- | --- |
| 0.75 min | 13.3 | 1600 | nein |
| 1.5 min | 6.7 | 800 | nein |
| 3 min | 3.3 | 400 | nein |
| 6 min | 1.7 | 200 | nein |
| 12 min | 0.8 | 100 | nein |
| 24 min | 0.4 | 50 | grenzwertig (2.4 Frames/Zyklus) |

Daraus folgen zwei Dinge, die in die Methodenbeschreibung gehören:

1. **Ein einzelner Zyklus ist nicht beobachtbar.** Eine scheinbare Periodizität
   in `31_ratio_*_over_time.pdf` wäre ein Alias-Artefakt, nicht der
   Medienwechsel — sie darf nicht als solcher interpretiert werden.
2. **Die Oszillation ist die Behandlung.** Bei gleichem Tastverhältnis erhalten
   alle Bedingungen dieselbe Gesamt-Feast- und Gesamt-Famine-Zeit und
   unterscheiden sich nur darin, wie fein sie zerhackt ist: ein **32-facher
   Dosisbereich** in der Anzahl der Wechsel. Interpretierbar ist deshalb nur
   die **kumulative** Wirkung über Stunden, nicht der Zyklusverlauf.

Die erwartete Richtung ergibt sich aus der Länge der Famine-Halbperiode: bei
0.75 min sind das ~22 s, die interne Metabolitpools mühelos überbrücken — die
Zelle sieht praktisch ein konstantes, gemitteltes Medium. Bei 24 min sind es
12 min, lang genug für echte Verarmung und eine Hungerantwort, ~50-mal in 20 h.
**Die stärkere Belastung wird bei den langsamen Zyklen erwartet**, nicht bei den
schnellen.

## Kumulativer Endzustand (Schritt 13)

Aus dem Abtast-Argument folgt, dass nur die **kumulative** Wirkung
interpretierbar ist. Schritt 13 bildet sie: pro Kammer der Mittelwert über ein
Endfenster, dann Kammer → Chip → Bedingung (`experiment_units.summarise_hierarchical`).

| Zweig | Endfenster |
| --- | --- |
| Oszillation | die letzten `ENDPOINT_LAST_FRACTION` (25 %) der Frames **jeder Kammer** |
| statisch | ein **absolutes** Fenster vor der frühesten Sättigung: ypd wächst über (bis ~4000 Tracks) und die W109-ypd-Aufnahmen enden bei 85 Frames; `detect_saturation_frame()` findet pro Kammer den Frame, ab dem die Zellzahl ≥ 90 % ihres Maximums bleibt — nur für Kammern, die überhaupt ≥ 1.5× gewachsen sind (`STATIC_SATURATION_*`). Ergebnis in `13_endpoint_saturation_per_chamber.csv`. |

| Datei | Inhalt |
| --- | --- |
| `13_endpoint_per_chamber.csv` / `_per_chip.csv` / `_summary.csv` | die drei Aggregationsstufen |
| `13_endpoint_bracket_score.csv` | pro Chip: `(osc − NegCtrl) / (PosCtrl − NegCtrl)`, 0 = wie Starvation, 1 = wie Feast, plus `bracket_degenerate` |
| `13_endpoint_spearman.csv` | Spearman ρ gegen die Periode auf Chip-Mittelwerten, `n_chips` = Perioden |
| `13_endpoint_vs_period_<spalte>_<osc_type>.pdf` | pro Periode der Chip-Wert mit **seinen** Kontrollen als Marker, gepooltes Kontrollband dahinter, Bracket-Score darunter; eine Facette pro Stamm |
| `static/13_endpoint_vs_medium_<spalte>.pdf` | statisch: Medium × Chip-Familie, Fehler über Chips |

**Bracket-Score.** Jede Periode ist ein eigener Chip; ihre Kontrollen liegen
auf demselben Chip und tragen denselben Chip-Effekt. Der Score entfernt ihn.
Wo `|PosCtrl − NegCtrl|` kleiner ist als 2× die Kammer-Streuung der Kontrollen,
trennt das Bracket nichts: der Chip wird hohl gezeichnet, der Score ist NaN und
fällt aus dem Trendtest — das ist der Befund aus Abschnitt 4, kein Fehler.

**Spearman, nicht Kruskal-Wallis**, weil die Vorhersage *monoton* in der
Periode ist. Bei n = 6 Chips braucht p < 0.05 ein |ρ| ≥ 0.83 — ρ ist eine
Effektstärke, die man berichtet. Kontrollen gehen nicht in die Korrelation
ein (sie haben keine Periode); unter drei Perioden gibt es keinen Test.


## PKO: verstopft Pullulan den Chip? (Zweig 3)

Die Kontrollkammern verhalten sich im Wildtyp nicht wie erwartet. Arbeits-
hypothese: **Pullulan** setzt die Chip-Strukturen zu. **PKO produziert kein
Pullulan** und liegt als ein Chip (Periode 3) vor.

### Die Größe, die zählt: Uneinigkeit nominell identischer Kammern

Auf jedem Chip liegen 3 PosCtrl- und 3 NegCtrl-Kammern — gleiche Kultur,
gleicher Tag, gleiches Medium, verschiedene Arrays. Ihre Uneinigkeit ist der
reinste Ausdruck dessen, was der *Chip* mit einer Kammer macht. Genau das
sollte eine Verstopfung aufblähen. `pko_comparison.py` rechnet pro Chip:

| Größe | Datei |
| --- | --- |
| **a)** CV der Kammer-Medianfläche über die 3 Kontrollkammern | `60_pko_within_chip_agreement.csv` |
| **b)** zeitlicher CV der Residuen um einen linearen Trend (Drift ohne Wachstum) | ebd. |
| **c)** Cliff's δ PosCtrl vs. NegCtrl auf µ_area (alle Zellen, alle Fits) | `60_pko_control_bracket.csv` |
| Abbildung | `61_pko_control_agreement.pdf` |

Vergleichsgruppe: **primär** die Produzenten-Chips derselben Periode (WT und
die Biosensor-Stämme bei 3 min; Entscheidung des Autors), gefüllt gezeichnet;
als **Hintergrund** alle Produzenten-Chips aller Perioden mit ihrer
10–90 %-Spanne — Kontrollkammern sind Konstant-Medium, also ist jeder Chip ein
gültiger Vergleichspunkt.

### Was nicht (mehr) geht

* **Kein Test über die Frequenz-Batches** — PKO hat eine Periode; der
  PKO-Kontext schaltet die Kontroll-Konsistenz deshalb ab.
* **Keine Zell-Ausbeute** — Blastokonidien werden zwischen Kammern gespült;
  die Zellzahl einer Kammer sagt nichts über Verstopfung. Ältere
  Ausgaben dazu werden beim nächsten Lauf gelöscht.
* **Kein Signifikanztest PKO gegen Produzenten** — PKO ist ein Chip. Die
  Aussage ist, wo er in der Verteilung der Produzenten-Chips liegt.

### Confound, der in die Diskussion gehört

PKO ist eine Mutante mit verändertem Stoffwechsel, nicht „WT ohne
Verstopfung“. Unterscheidbar sind hydraulische und metabolische Erklärung nur
über das **Muster**: hydraulisch = hohe Kammer-Streuung und Drift bei
Produzenten, metabolisch = Niveau-Verschiebung bei erhaltener Übereinstimmung.


## Mit QC vs. ohne QC

Das manuelle QC (`qc_exclusions.csv`: Track-Merges, Ausschlüsse) deckt nur die
Batches ab, die man wirklich durchgesehen hat — im echten Datensatz **einen**
(WT/pH/6, ~500 Zeilen). Die Pipeline läuft die Oszillationsschritte für genau
diese Batches zusätzlich auf den **Rohdaten** (keine Merges, keine Ausschlüsse)
nach `analysis_output/no_qc/` und stellt in `analysis_output/qc_comparison/`
Kammer für Kammer gegenüber, was das QC geändert hat:

| Datei | Inhalt |
| --- | --- |
| `70_qc_effect.pdf` | mit QC (x) gegen ohne QC (y), ein Punkt pro Kammer; auf der Diagonale = kein Effekt |
| `70_qc_effect_summary.csv` | Median der relativen Änderung je Kennzahl und Kontrollart |
| `70_qc_exclusion_inventory.csv` | Zeilen je Kammer, Aktion (merge/exclude) und Grund |
| `70_qc_exclusions_conflicts.csv` | Tracks, die mehrfach gelistet sind: `double_merge` (erste Zeile gewinnt, zweite wird verworfen), `merge_and_exclude`, `duplicate` — mit Zeilennummern |

`--skip-qc-comparison` lässt beides weg.

## Die Lineage-Heuristik validieren

`lineage.classify_mother_bud()` entscheidet über Budding Ratio, µ_event, den
Stammbaum und die Mutter/Knospe-Trennung. Es ist eine Heuristik aus
Centroid-Abstand und Tracklänge, **kein** echtes Lineage-Tracking — deshalb
gehört vor jede Aussage eine Validierung. Beide Werkzeuge lesen die Dateien,
die `run_analysis.py` erzeugt:

```
analysis_output/00_cell_positions.parquet   Zellpositionen NACH QC
analysis_output/20_budding_events.csv       erkannte Events (je eine Datei
analysis_output/static/20_budding_events.csv   pro Teil-Pipeline)
```

**1. Quantitativ, über alle Kammern:**

```bash
cd analyse_pipeline
python validate_lineage.py          # Pfade kommen aus config.py
```

Erzeugt in `analysis_output/lineage_validation/`:

| Datei | Frage, die sie beantwortet |
| --- | --- |
| `lv_01_d_over_r_distribution.pdf` | Lagen die Buds komfortabel im Suchradius, oder hat die Toleranz sie gerade noch hereingeholt? |
| `lv_02_detection_rate.pdf` + `_per_chamber.csv` | Ist die Erkennungsrate über die Bedingungen konstant? |
| `lv_03_detection_rate_kruskal.csv` | Kruskal-Wallis dazu: p < 0.05 = Erkennung mit der Bedingung konfundiert |
| `lv_03_assignment_ambiguity.pdf` | Wie oft kamen mehrere Mütter in Frage (greedy Nearest-Neighbour)? |
| `lv_04_tolerance_sweep.pdf` | Sitzt `tolerance_px` auf einem Plateau oder auf einer Flanke? |

Die Kurzfassung steht am Ende im Log — inklusive Warnung, wenn zu viele
Zuordnungen grenzwertig oder mehrdeutig sind.

**2. Visuell, eine Kammer:**

```bash
python plot_qc_lineage_overlay.py \
    --cells   ../analysis_output/00_cell_positions.parquet \
    --lineage-events ../analysis_output/20_budding_events.csv \
    --qc-tif  ".../QC/260616_Osc1.5_NegCtrl_Rep1_ChamA13_QC_overlay.tif" \
    --output  ChamA13_lineage_overlay.tif
```

`--exp-id` wird automatisch bestimmt (über `cells['filename']`, sonst über
Pfad und Dateinamen). Marker im Ausgabe-TIFF:

| Marker | Bedeutung |
| --- | --- |
| türkiser Kreis, klein | Mutter-Zentroid |
| türkiser Kreis, groß | adaptiver Suchradius dieser Mutter |
| oranger Kreis | zugeordneter Bud, `d/r` = Distanz / Suchradius |
| gelbe Linie | Mutter-Bud-Zuordnung |
| hellblauer Kreis `?` | Bud-Kandidat **ohne** Mutter — fällt aus allen Auswertungen |

Überlappen sich die Suchradien benachbarter Mütter im Bild, ist jede
Zuordnung in diesem Bereich eine Entscheidung der Heuristik, keine
Beobachtung — das ist der wichtigste Blick beim Kalibrieren.

---

## Bekannte Einschränkungen

Diese Punkte werden bei jedem Lauf als Warnung geloggt
(`config.METHOD_CAVEATS` / `config.DOCUMENTED_DEFAULTS`). Sie sind **nicht**
stillschweigend geändert worden — die Entscheidung darüber ist eine fachliche.

**Parameter weichen vom dokumentierten Standard ab**

| Parameter | aktiv | dokumentiert | Folge |
| --- | --- | --- | --- |
| `MU_MAX_THRESHOLD` | `10.0` | `0.6` h⁻¹ | Der Artefakt-Filter greift praktisch nicht mehr; `mu_is_artefact` bleibt fast überall `False`. |
| `LINEAGE_PARAMS.mother_min_frames` | `10` | `20` | Mehr, aber kürzer beobachtete Tracks gelten als Mutterzelle. |
| `LINEAGE_PARAMS.bud_max_frames` | `7` | `5` | Reine Report-Schwelle (`bud_was_washed_out`), ohne Einfluss auf die Erkennung. |

**Statistik**

* **Innerhalb einer Oszillationsbedingung gibt es keine biologische
  Replikation** (ein Chip, eine Vorkultur; siehe *Versuchseinheiten*). Alle
  Aggregate, die über `experiment_units.summarise_hierarchical()` laufen
  (`13_*`, `40_*_aggregated`, `12_*_per_chip`), beschriften den Fehlerbalken
  in `error_unit` (`chamber` = technisch, `chip` = biologisch) und zählen
  `n_units` entsprechend. Die Spalte `n_replicates` gibt es nicht mehr.
* `summarise_growth_rate()` und `summarise_area_growth()` (`11_*_summary.csv`,
  `12_*_summary_*.csv`) aggregieren weiterhin über **einzelne Zellen/Intervalle**
  — `sd_mu`/`sd_mu_area` beschreiben die Streuung über Zellen.
* Die Mann-Whitney-Sternchen in Panel A laufen über Zellen bzw. Mutterzellen
  und sind rein deskriptiv; die Abbildung sagt das in der Fußnote.
* Spearman gegen die Periode läuft auf Chip-Mittelwerten (n = Perioden) und
  ist bei n ≤ 6 eine Effektstärke, kein Test. Die Stämme werden nicht als
  Replikate gepoolt.
* Robustheit **R ist relativ**: der Normalisierungsfaktor `m` wird über den
  gesamten übergebenen Datensatz gebildet — R-Werte aus `analysis_output/`,
  `static/`, `pko/` und `no_qc/` dürfen **nicht** gegeneinander gelesen werden.
  Der PKO-Vergleich benutzt deshalb gewöhnliche Kammer-Statistiken statt R.
* R(t) ist für die **oszillierenden** Bedingungen alias-konfundiert (der
  Alias-Beitrag wächst mit der Periode, in Richtung des erwarteten Effekts)
  und bekommt dort keinen Trendtest; für Konstant-Medium-Kontrollen ist R(t)
  sauber.
* `fit_is_reliable` (R² ≥ 0.5) verzerrt flache Bedingungen (NegCtrl) nach
  oben — der Filter behält dort nur Tracks, in denen Rauschen wie ein Trend
  aussieht. `compute_control_bracket()` und `12_area_growth_rate_all.pdf`
  filtern deshalb nicht und führen den Anteil zuverlässiger Fits als Diagnose;
  `summarise_area_growth(exclude_unreliable=True)` und R(p) für µ_area
  filtern weiterhin.

**Einheiten & Proxys**

* `area` in den Zelltabellen und Plots ist die rohe Cellpose-Fläche in **px²**.
  Nur `area_growth.py` rechnet intern über `PX_TO_UM2` (1 µm = 13.63 px) in µm²
  um — das verschiebt nur den Achsenabschnitt des Fits, nicht µ_area selbst.
* `eccentricity` ist **nicht** dasselbe wie die Circularity des Referenzpapers
  (4π·A/U²) und wird in Panel A ausdrücklich als Stellvertreter beschriftet.
* Die Einheit von `osc_freq` (`0.75 … 24`) steht nirgends im Code. Die
  Kontrollkammer-Plots beschriften sie deshalb ohne Einheit.
* Die Konzentrations-Beschriftung der Kontrollen (`0 g/L` / `50 g/L`) steht in
  `config.CONTROL_CONCENTRATION_LABELS` **pro Oszillationstyp** — sie gilt nur
  für Glucose, nicht für pH.

**Heuristik**

* `lineage.py` leitet Mutter/Bud aus räumlicher Nähe und Tracklänge ab; es gibt
  **kein** echtes Lineage-Tracking aus der Bildverarbeitung. Vor jeder
  Publikation mit `inspect_lineage.py` gegen echte QC-Overlays kalibrieren.
* `tolerance_px` ist ein Pixel-Wert und damit von Kamera/Optik abhängig.
