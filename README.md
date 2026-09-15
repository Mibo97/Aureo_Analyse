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
```

Die Metadaten `biosensor`/`osc_type`/`osc_freq` kommen aus dem **Ordnerpfad**,
nicht aus dem Dateinamen (siehe `data_loading.py`). Ordner mit
`osc_type == "static"` laufen durch dieselbe Auswertung, landen aber in
`analysis_output/static/`.

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
| `control_consistency.py` | Kruskal-Wallis: driften PosCtrl/NegCtrl über Batches? |
| `queen_controls.py` | PosCtrl-vs-NegCtrl-Validierung der Sensoren selbst |
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
| `20_`–`23_` | Lineage: Budding-Events, Budding Ratio, Panel A, Stammbaum |
| `30_`–`31_` | Sensor-Intensitäten und Ratios über die Zeit |
| `40_` | Robustheit R(t)/R(p) inkl. Kontroll-Konsistenz |
| `50_` | Zusammenfassungstabelle |
| `90_`–`92_` | Anhang: Morphologie-Scatter, Einzelzell- & Mutter-Trajektorien |
| `95_` | Anhang: Sensor-Controls (PosCtrl vs. NegCtrl pro Biosensor) |

Statische Daten landen in denselben Präfixen unter `analysis_output/static/`
(ohne `30_`/`31_`/`95_`, da dort nur der Wildtyp ohne Fluoreszenzkanäle läuft).

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
  Replikate. `aggregate_robustness_over_replicates()` und
  `analysis._aggregate_over_replicates()` mitteln dagegen korrekt zweistufig.
* Die Mann-Whitney-Tests in den Violin-Plots laufen über die übergebenen Zeilen
  (eine Mutterzelle pro Zeile bei der Budding Ratio) und berücksichtigen die
  Replikat-Struktur ebenfalls nicht.
* Robustheit **R ist relativ**: der Normalisierungsfaktor `m` wird über den
  gesamten übergebenen Datensatz gebildet. R-Werte aus Läufen mit
  unterschiedlichem Datenumfang sind nicht miteinander vergleichbar.

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
