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
| `00_` | Übersicht / Sanity-Check (Tracks pro Experiment) |
| `10_`–`12_` | Zellmorphologie & Wachstum (Fläche, µ_event, µ_area) |
| `13_` | **Kumulativer Endzustand gegen die Periode** (+ Spearman) |
| `20_`–`23_` | Lineage: Budding-Events, Budding Ratio, Panel A, Stammbaum |
| `30_`–`31_` | Sensor-Intensitäten und Ratios über die Zeit |
| `40_` | Robustheit R(t)/R(p) inkl. Kontroll-Konsistenz |
| `50_` | Zusammenfassungstabelle |
| `90_`–`92_` | Anhang: Morphologie-Scatter, Einzelzell- & Mutter-Trajektorien |
| `95_` | Anhang: Sensor-Controls (PosCtrl vs. NegCtrl pro Biosensor) |
| `60_`–`62_` | **Nur in `pko/`**: WT-gegen-PKO-Vergleich (siehe unten) |

Statische Daten und PKO-Daten landen in denselben Präfixen unter
`analysis_output/static/` bzw. `analysis_output/pko/` — beide ohne
`30_`/`31_`/`95_`, da dort keine Fluoreszenzkanäle vorliegen.

---

## Abtastung: die Oszillation ist eine Behandlung, keine Messgröße

Die Spalte `osc_freq` enthält trotz ihres Namens die **Periode in Minuten**
(0.75 … 24), keine Frequenz. Bei `MIN_PER_FRAME = 10` liegt die kürzeste
auflösbare Periode (Nyquist) bei **20 min** — fünf der sechs Bedingungen liegen
darunter, die sechste nur knapp darüber:

| Periode | Zyklen pro Frame | Zyklen in 10 h | auflösbar? |
| --- | --- | --- | --- |
| 0.75 min | 13.3 | 800 | nein |
| 1.5 min | 6.7 | 400 | nein |
| 3 min | 3.3 | 200 | nein |
| 6 min | 1.7 | 100 | nein |
| 12 min | 0.8 | 50 | nein |
| 24 min | 0.4 | 25 | grenzwertig (2.4 Frames/Zyklus) |

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
12 min, lang genug für echte Verarmung und eine Hungerantwort, 25-mal in 10 h.
**Die stärkere Belastung wird bei den langsamen Zyklen erwartet**, nicht bei den
schnellen.

## Kumulativer Endzustand (Schritt 13)

Aus dem Abtast-Argument oben folgt, dass nur die **kumulative** Wirkung
interpretierbar ist. Genau diese Auswertung fehlte: `10_`/`30_`/`31_` sind
Zeitreihen, `40_` sind Varianzmaße, und `50_summary_statistics.csv` bekommt nur
`intensity_cols` übergeben — die `ratio_*`-Spalten erscheinen dort **nie**. Für
die Sensor-Daten ist Schritt 13 die erste kumulative Auswertung überhaupt.

Der Endzustand ist der Mittelwert über die letzten `ENDPOINT_LAST_FRACTION`
(Default 25 %) der Frames **jeder Kammer** — relativ zur Kammer, nicht absolut,
weil Kammern unterschiedlich lang aufgenommen sein können. Danach wird
dreistufig aggregiert: Kammer → Replikat → Bedingung.

| Datei | Inhalt |
| --- | --- |
| `13_endpoint_per_replicate.csv` | ein Wert je Replikat und Bedingung (Grundlage des Tests) |
| `13_endpoint_summary.csv` | Mittelwert ± SEM über Replikate je Bedingung |
| `13_endpoint_spearman.csv` | Spearman ρ und p gegen die Periode, je Stamm/Oszillationstyp |
| `13_endpoint_vs_period_<spalte>.pdf` | die Abbildung, mit PosCtrl/NegCtrl als Referenzbändern |

**Spearman, nicht Kruskal-Wallis.** Die Vorhersage ist *monoton* in der Periode
(längere Famine-Halbzyklen belasten mehr, siehe oben) — geprüft wird also eine
Rangkorrelation. Kruskal-Wallis prüft „irgendeine Gruppe unterscheidet sich“ und
ist dafür das falsche Werkzeug; es war bisher der einzige verdrahtete Test, und
das nur für die Kontrollen. Getestet wird auf **Replikat**-Ebene: über Zellen
gerechnet hängt das p fast nur an der Zellzahl.

Die Kontrollen gehen **nicht** in die Korrelation ein — sie haben keine Periode,
der Ordnername ihres Batches ist keine Behandlung. Sie erscheinen als
waagerechte Bänder, was der Zweck einer Kontrolle ist.

Die x-Achse ist **logarithmisch**: die Perioden sind geometrisch gestuft (Faktor
2 von 0.75 bis 24 min, ein 32-facher Dosisbereich). Linear dargestellt drängen
sich fünf der sechs Bedingungen links zusammen.

Ein Trendtest braucht mindestens **drei** verschiedene Perioden — bei zwei
Punkten ist ρ immer ±1, das ist Arithmetik und keine Evidenz. Für die
statischen Daten und den PKO-Zweig entfällt der Test daher; das wird geloggt.

---

## PKO: verstopft Pullulan den Chip? (Zweig 3)

Die Kontrollkammern verhalten sich im Wildtyp nicht wie erwartet. Arbeits-
hypothese: **Pullulan** — das Exopolysaccharid, das der WT ausscheidet — setzt
die Chip-Strukturen zu und erzeugt unerwartete Strömungsprofile. Der
**PKO-Stamm produziert kein Pullulan**; verhalten sich *seine* Kontrollen
korrekt, stützt das die Clogging-Erklärung.

### Was nicht geht: der Test über die Frequenz-Batches

`test_control_consistency_across_freq()` vergleicht die Kontrollen **über die
`osc_freq`-Batches hinweg** und braucht dafür mindestens zwei. PKO deckt nur
**eine Periode** ab — der Test liefert dort `p = NaN`, und
`plot_control_consistency()` zeichnet einen *einzelnen Punkt* pro Kontrollart,
also eine trivial flache Linie, die wie „PKO-Kontrollen sind konsistent“
aussieht und nichts enthält.

Der PKO-Zweig schaltet die Kontroll-Konsistenz deshalb ab
(`PipelineContext.run_control_consistency=False`, automatisch gesetzt, sobald
weniger als zwei Batches vorliegen), statt eine irreführende Abbildung zu
erzeugen. Der Vergleich läuft stattdessen **innerhalb des gemeinsamen
Perioden-Batches, WT gegen PKO**. Das ist kein Notbehelf:

* PosCtrl ist durchgehend Feast, NegCtrl durchgehend Starvation — der
  Medienverlauf einer *Kontrollkammer* hängt gar nicht an der Periode des
  Batches, in dem sie mitlief.
* In den Kontrollkammern fällt am meisten Pullulan an: PosCtrl wächst 10 h
  durch. Wenn Verstopfung das Problem ist, ist das der Ort dafür.

### Drei Größen, keine davon lineage-abhängig

| Datei | Inhalt |
| --- | --- |
| `60_pko_control_chambers.csv` | ein Wert pro Kontrollkammer (zweistufiger Median) |
| `60_pko_chamber_agreement_per_replicate.csv` / `_per_strain.csv` | Kammer-zu-Kammer-CV **innerhalb** eines Replikats |
| `60_pko_control_bracket_per_replicate.csv` / `_per_strain.csv` | PosCtrl-vs-NegCtrl als Cliff's δ auf µ_area |
| `60_pko_cell_yield_timeseries.csv` / `_slopes.csv` | verfolgbare Zellen pro Kammer über die Zeit |
| `61_pko_control_agreement.pdf` | **Hauptabbildung**: Kammer-Übereinstimmung + Bracket, WT gegen PKO |
| `62_pko_cell_yield.pdf` | Zell-Ausbeute pro Kontrollkammer, WT gegen PKO |

Bewusst **ohne** Mutter/Bud-Heuristik: µ_area ist eine Regression über
ln(Fläche) eines Tracks und damit von `lineage.py` unabhängig — nur das
`cell_type`-Label hängt daran. Damit steht dieser Vergleich nicht auf einer
unvalidierten Heuristik (siehe `validate_lineage.py`).

Die Kammer-Streuung wird **innerhalb** eines Replikats gebildet: Kammern eines
Replikats sind technische Messungen desselben Chips, ihre Streuung ist der
hydraulische Anteil. Die Streuung *zwischen* Replikaten ist biologisch und
würde den gesuchten Effekt nur verwässern.

### Was dieses Design statistisch trägt

PKO ist **ein Stamm**. Ein Signifikanztest WT-gegen-PKO auf Stammebene hat
n = 1 in einer Gruppe und wird deshalb **bewusst nicht gerechnet** — auch nicht
über Replikate oder Kammern hinweg, das wäre Pseudoreplikation auf Stammebene.

Getragen wird eine **deskriptive** Aussage: die vier WT-Biosensor-Stämme
liefern vier voneinander unabhängige WT-Werte, und der PKO-Wert liegt
innerhalb oder außerhalb dieser Spanne. Die vier WT-Stämme sind damit die
interne Replikation der WT-Seite — stimmen sie untereinander nicht überein,
ist schon das Zusammenfassen zu „WT“ falsch, und `61_` zeigt genau das.

### Confound, der in die Diskussion gehört

PKO ist nicht „Wildtyp ohne Verstopfung“, sondern eine Mutante mit verändertem
Kohlenstofffluss und veränderten Oberflächeneigenschaften. **Jeder** WT-PKO-
Unterschied lässt sich auch direkt physiologisch erklären — einen Unterschied
zu finden ist noch kein Beleg für Verstopfung. Unterscheidbar sind die beiden
Erklärungen nur über das **Muster**:

| Ursache | Vorhersage |
| --- | --- |
| hydraulisch | hohe Kammer-zu-Kammer-Streuung, wegbrechende Zell-Ausbeute, kaputtes Bracket |
| metabolisch | gleichmäßige Niveau-Verschiebung, Kammer-Übereinstimmung bleibt erhalten |

Deshalb berichtet `pko_comparison.py` Streuungen und Steigungen, nicht nur
Mittelwerte.

---

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

* `summarise_growth_rate()` und `summarise_area_growth()` aggregieren über
  einzelne Zellen/Intervalle, **nicht** erst pro Replikat. `sd_mu`, `sd_mu_area`
  und `n_values` in `11_*_summary.csv` / `12_*_summary_*.csv` beschreiben damit
  die Streuung über Zellen (Pseudoreplikation), nicht über biologische
  Replikate. `analysis._aggregate_over_replicates()` (Schritte `10_`/`30_`/`31_`)
  und `queen_controls.summarise_sensor_controls()` (Schritt `95_`) gruppieren
  dagegen auf `replicate` und mitteln damit korrekt.
* `aggregate_robustness_over_replicates()` (alle `40_*_aggregated`) aggregiert
  seit der Umstellung auf `replicate_id_col="replicate"` **dreistufig**:
  Frame → Kammer → Replikat. Vorher war `exp_id` der Default, das laut
  `data_loading.py` auch `chamber` enthält — die Aggregation endete also auf
  Kammer-Ebene und zählte jede Kammer als biologisches Replikat (bei 3×2 stand
  in `n_replicates` eine 6, und `sd` mischte technische mit biologischer
  Varianz). `n_replicates` nennt jetzt die echte Replikatzahl, und `sd`
  beschreibt die Streuung zwischen Replikaten. Wer die alte Kammer-Ebene
  braucht, ruft mit `replicate_id_col="exp_id"` auf.
* Mit 3 Replikaten ist `sd` **schlecht geschätzt** — das ist der Preis dafür,
  dass sie jetzt das Richtige beschreibt. Die Fehlerbalken werden dadurch in
  der Regel breiter und ehrlicher, nicht enger.
* Die Mann-Whitney-Tests in den Violin-Plots laufen über die übergebenen Zeilen
  (eine Mutterzelle pro Zeile bei der Budding Ratio) und berücksichtigen die
  Replikat-Struktur ebenfalls nicht.
* Robustheit **R ist relativ**: der Normalisierungsfaktor `m` wird über den
  gesamten übergebenen Datensatz gebildet. R-Werte aus Läufen mit
  unterschiedlichem Datenumfang sind nicht miteinander vergleichbar. Das gilt
  auch **zwischen den drei Zweigen**: Oszillation, statisch und PKO bekommen
  je einen eigenen `PipelineContext` und damit je ein eigenes `m` — R-Werte aus
  `analysis_output/`, `static/` und `pko/` dürfen **nicht** gegeneinander
  gelesen werden. Genau deshalb benutzt `pko_comparison.py` für den
  WT-gegen-PKO-Vergleich gewöhnliche Kammer-Statistiken (Median, CV) statt R.
* `fit_is_reliable` heißt `R² >= min_r_squared` (Default 0.5, siehe
  `area_growth.compute_area_growth_rate()`). Bei einer Bedingung, die
  **tatsächlich flach ist**, erklärt die Regressionsgerade per Konstruktion
  kaum Varianz — das R² ist niedrig, *weil es nichts zu erklären gibt*. Der
  Filter wirft dann die ehrlichen flachen Fits weg und behält die, in denen
  Rauschen wie ein Trend aussieht; der überlebende Median ist nach **oben**
  verzerrt (Survivorship Bias). Das trifft genau `NegCtrl`: durchgehende
  Starvation *soll* µ_area ≈ 0 liefern. Betroffen sind
  `summarise_area_growth(exclude_unreliable=True)` (also
  `12_area_growth_rate_summary_*.csv`) und die R(p)-Rechnung für µ_area in
  Schritt 40. `pko_comparison.compute_control_bracket()` filtert deshalb
  **nicht** und berichtet stattdessen den Anteil zuverlässiger Fits pro Arm
  (`frac_reliable_*`) — ein niedriger Anteil in NegCtrl ist selbst ein Befund.
  In den übrigen Schritten ist das **nicht** korrigiert: die Entscheidung
  darüber ist eine fachliche.

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
