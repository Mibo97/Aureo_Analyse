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

> **Die Argumentation der Arbeit** — welche Abbildung und Tabelle welchen Schritt
> trägt, von den Versuchseinheiten über die Tracking-Grenze bis zum Befund, dass die
> Kontrollen die Trends tragen — steht in [`docs/data_story.md`](docs/data_story.md).

> Der Befund zur Cellpose-Pipeline und zum Tracking — warum Tracks abreißen, was ein
> Re-Linker auf der Tabelle noch retten kann und welche Optionen es gibt — steht in
> [`docs/tracking_diagnosis.md`](docs/tracking_diagnosis.md); `analyse_pipeline/diagnose_tracking.py`
> reproduziert die Zahlen aus einer Combined_Results-Datei.
> Werkzeuge zum Re-Tracking aus den gespeicherten Masken, zum Segmentierungs-Sweep und fuer
> Slurm-Array-Jobs: [`imaging/README.md`](imaging/README.md); Plan: [`docs/tracking_plan.md`](docs/tracking_plan.md).

## Versuchseinheiten: was ist hier ein Replikat?

Die Spalten `replicate` und `chamber` kommen aus dem Dateinamen der
Bildverarbeitung und bedeuten **nicht**, was ihre Namen nahelegen.
`experiment_units.py` leitet daraus die echten Einheiten ab und schreibt
`00_chip_overview.csv` (eine Zeile pro Chip).

| Zweig | pro (Stamm, osc_type, Periode) | `replicate` ist … | biologische Einheit | n je Bedingung |
| --- | --- | --- | --- | --- |
| Oszillation, PKO | **eine Struktur** (im Code und in allen Tabellen `chip`); ~11 Kammern (5 Osc, 3 PosCtrl, 3 NegCtrl) auf mehreren Arrays. Ein **physischer Chip trägt 2–3 Strukturen** = 2–3 Perioden aus **einer Vorkultur an einem Tag** (Spalte `culture`, `00_chip_overview.csv`: `n_structures_in_culture`) | ein Array-Index, der gleiche Positionslabels (`ChamA13` ×3) auseinanderhält | die **Kultur** (physischer Chip); die Struktur ist die Behandlungseinheit | **1** Struktur; 2 Kulturen je Serie (Glc: 3 + 3 Perioden, pH: 2 + 2) |
| statisch W109 | — | ein Chip je `replicate` **(Annahme — die vier `replicate` eines Mediums tragen dasselbe Datum; ob sie ein Chip waren, ist offen)**; `chamber` ist immer `ChamA0` | der Chip | 4 pro Medium |
| statisch W65 | — | eine **Kammer** auf dem einen W65-Chip: ein Chip, eine Vorkultur, beide Medien (`STATIC_SINGLE_CHIP_FAMILIES`) | die Kultur (n = 1); Fehlerbalken über Kammern | 5 Kammern pro Medium |

Drei Folgen, die in die Arbeit gehören:

1. **Innerhalb einer Oszillationsbedingung gibt es keine biologische
   Replikation.** Alle Kammern einer Periode sind technische Wiederholungen
   einer Kultur auf einem Chip. Jeder Fehlerbalken dort ist ein
   Kammer-Fehlerbalken; die Tabellen sagen das in `error_unit` (`chamber`
   vs. `chip`) und zählen `n_units` entsprechend. Die Dosis-Wirkung über die
   Perioden ist **ein Chip pro Dosis** (6 bei Glc, 4 bei pH); Spearman läuft
   auf Chip-Mittelwerten mit n = Zahl der Perioden und ist damit eine
   Effektstärke, kein Test, auf den man sich stützt.
2. **Die Kontrollen liegen auf derselben Struktur wie die Behandlung.** Jede
   Periode kann deshalb gegen *ihre eigenen* PosCtrl/NegCtrl normiert werden
   (Bracket-Score, Schritt 13) — das entfernt den Struktur-/Kultur-Effekt,
   ohne dass man ihn kennen müsste. Und die Kontrollen sind der Prüfstein für
   jeden Trend gegen die Periode: sie liegen in konstantem Medium und können
   auf die Periode nicht reagieren. `13_endpoint_control_trend.csv`,
   `12_area_growth_rate_control_trend.csv` und `21_budding_rate_control_trend.csv`
   stellen dem Spearman der
   Oszillationskammern den der Kontrollen derselben Strukturen gegenüber
   (`rho_ctrl_mean`, `rho_ctrl_strongest`) und den der Differenz
   (`rho_osc_minus_ctrl`). Läuft auch nur eine Kontrollart mit, trägt die
   Struktur den Trend, nicht die Periode. Ein `period effect` verlangt drei
   Dinge zugleich: die Oszillationskammern trenden, keine Kontrollart trendet
   gleichsinnig, und die Differenz trendet ebenfalls. In den echten Daten
   (`50_control_trend_summary.csv`: 48 Readout-Serien-Kombinationen aus
   Endzustand, µ_area und Knospungsrate) erfüllen das 7, verstreut über
   Readouts und Stämme: keine kehrt für einen Stamm in beiden
   Oszillationstypen gleichsinnig wieder, bei vier läuft die stärkste
   Kontrolle gegenläufig. 18 tragen den Trend ihrer Kontrollen mit, 19 zeigen
   keinen Trend; die Sensor-Ratios (OxPro, pHluorin) laufen in den
   NegCtrl-Kammern genauso mit der Periode wie in den Oszillationskammern.
   `*_within_culture.csv` zeigt dasselbe innerhalb einer Kultur (2–3
   Perioden, eine Vorkultur), wo der Kulturvergleich entfällt.
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
| `relink.py` | Gap Closing und Sparse-Phase-Fenster gegen die Track-Fragmentierung |
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
| `00_` | Übersicht / Sanity-Check; `00_chip_overview.csv` = eine Zeile pro **Chip**; `00_chip_run_order.csv` = wurden die Perioden einer Serie in Datumsreihenfolge gefahren (dann ist ein Periodentrend nicht von einer Tagesdrift zu trennen)?; `00_track_fragmentation.csv` / `00_track_relinks.csv` = Track-Fragmentierung und automatisches Gap Closing |
| `10_`–`12_` | Zellmorphologie & Wachstum (Fläche, µ_event, µ_area) |
| `13_` | **Kumulativer Endzustand gegen die Periode** (+ Spearman) |
| `20_`–`23_` | Lineage: Budding-Events, Budding Ratio, Panel A, Stammbaum — **nur aus dem Sparse-Phase-Fenster** (`20_lineage_window.csv/.pdf`, siehe unten); `21_budding_rate_vs_period_<osc_type>.pdf` = Knospungsrate je Mutter-Stunde gegen die Periode mit eigenen Kontrollen; `20_bud_size_*` (nur direkt in `analysis_output/`) = Größenkriterium der Knospen-Heuristik, eine Schwelle für alle Zweige |
| `30_`–`31_` | Sensor-Intensitäten und Ratios über die Zeit |
| `40_` | Robustheit R(t)/R(p) inkl. Kontroll-Konsistenz |
| `50_` | Zusammenfassungstabelle; `50_control_trend_summary.pdf/.csv` = **die eine Abbildung zum Kontroll-Trend**: je Readout und Serie der Spearman der Oszillationskammern gegen den der stärksten Kontrolle derselben Strukturen (aus `12_`, `13_`, `21_`) |
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

Ausschlüsse werden **vor und nach** den Merges angewendet (ein Ausschluss auf
dem Quell-Track greift vor dem Umbenennen, einer auf dem Ziel-Track trifft
danach auch die hineingemergten Frames). `frame_from`/`frame_to` gelten auch
für Merges: zwei Merge-Zeilen derselben `cell_uid` mit disjunkten Bereichen
und verschiedenen Zielen teilen einen Track mit Tracker-Sprung auf.

| Datei | Inhalt |
| --- | --- |
| `70_qc_effect.pdf` | mit QC (x) gegen ohne QC (y), ein Punkt pro Kammer; auf der Diagonale = kein Effekt |
| `70_qc_effect_summary.csv` | Median der relativen Änderung je Kennzahl und Kontrollart |
| `70_qc_exclusion_inventory.csv` | Zeilen je Kammer, Aktion (merge/exclude) und Grund |
| `70_qc_exclusions_conflicts.csv` | Tracks, die mehrfach gelistet sind, mit Konsequenz: `split_merge` (zwei Merges mit disjunkten Frame-Bereichen — beide werden ausgeführt, der Track wird aufgeteilt), `double_merge` (zwei Merges ohne Bereiche — die erste Zeile gewinnt), `merge_and_exclude` (beides wird ausgeführt), `redundant_exclusion`, `duplicate` — mit Zeilennummern |

`--skip-qc-comparison` lässt beides weg.

Ergebnis für den QC-Batch: über den ganzen Lauf ändert das manuelle QC 2,6 %
der Budding-Events, im Sparse-Phase-Fenster 19 %. Es korrigiert einzelne
Zellen, behebt aber nicht die Track-Fragmentierung (siehe
„Sparse-Phase-Lineage“). Der Lauf ohne QC enthält dieselben automatischen
Schritte (Gap Closing, Fenster) — verglichen wird wirklich nur das manuelle QC.

## Die Lineage-Heuristik validieren

`lineage.classify_mother_bud()` entscheidet über Budding Ratio, µ_event, den
Stammbaum und die Mutter/Knospe-Trennung. Es ist eine Heuristik aus
Centroid-Abstand und Tracklänge, **kein** echtes Lineage-Tracking — deshalb
gehört vor jede Aussage eine Validierung. Beide Werkzeuge lesen die Dateien,
die `run_analysis.py` erzeugt:

```
analysis_output/00_cell_positions.parquet   Zellpositionen NACH QC und Gap Closing,
                                            Spalte in_lineage_window
analysis_output/20_budding_events.csv       erkannte Events (je eine Datei
analysis_output/static/20_budding_events.csv   pro Teil-Pipeline)
```

### Sparse-Phase-Lineage: die Heuristik läuft nur im dünn besetzten Feld

Der Tracker der Bildverarbeitung vergibt in diesen Daten mit etwa 12 % pro
Objekt und Frame eine neue Track-ID, unabhängig von der Zelldichte (QC-Batch
WT/pH/6: 9 063 Tracks in 11 Kammern, Median-Tracklänge 3 Frames, 30 %
Ein-Frame-Tracks). Jede „neu auftauchende Zelle“ ist damit zunächst ein
Fragment. Weil sich die Kammern von 2 auf bis zu 180 Objekte füllen, folgt die
Zahl der Budding-Events der Dichte: im QC-Batch 462 pro Kammer, davon 243 in
den letzten 22 Frames — Fragment-Statistik, keine Biologie. Das manuelle QC
(248 Merges, 214 Hintergrund-Ausschlüsse) entfernte davon 2,6 %.

Deshalb:

1. **Fenster je Kammer** (`20_lineage_window.csv/.pdf`, `LINEAGE_SPARSE_*` in
   `config.py`): vom Start bis zum letzten Frame, bevor der rollende Median der
   Objekte pro Frame 20 übersteigt; Kammern mit weniger als 20 solchen Frames
   fallen aus der Lineage-Auswertung. Budding Ratio, Panel A, die
   Budding-Ratio-Zeitreihe, µ_event, der Stammbaum und die Mutter/Knospe-
   Trennung von µ_area stammen nur aus diesem Fenster
   (`PipelineContext.cells_lineage`). Endzustand, µ_area und R(t)/R(p) sehen
   weiterhin den ganzen Lauf.
2. **Gap Closing im Fenster** (`00_track_relinks.csv`, `RELINK_*`): ein neu
   beginnender Track wird an einen höchstens 2 Frames vorher beendeten Track
   angehängt, wenn Abstand ≤ 150 px, Flächenverhältnis in [0.5, 2] und die
   Zuordnung in beide Richtungen eindeutig ist. Manuelle Merges aus
   `qc_exclusions.csv` laufen vorher und haben Vorrang. Kalibrierung am
   QC-Batch: von den 70 manuellen Verknüpfungen im Fenster findet die Regel
   40 %; der Rest ist mehrdeutig (Sprünge derselben Zelle von median 97 px bei
   mehreren Kandidaten) und bleibt bewusst offen — lieber ein Bruch zu viel
   als zwei Zellen vermischt. Größere Radien oder Lücken finden *weniger*,
   weil die Mehrdeutigkeit schneller wächst als die Trefferzahl.
3. **Kennzahlen** (`00_track_fragmentation.csv`): Objekte pro Frame, Tracks,
   Median-Tracklänge, neue Tracks je Objekt und Frame, Anteil der
   Objekt-Frames in Tracks ≥ 10 Frames — vor und nach dem Gap Closing, je
   Kammer. Diese Tabelle gehört in die Arbeit, sobald Lineage-Ergebnisse
   gezeigt werden.

Im QC-Batch bleiben im Fenster 136 Events in 11 Kammern (manuelles QC: 128,
beides zusammen: 116), 118 davon mit einer über ≥ 30 Frames verfolgten Mutter
— die Größenordnung, die das Zellwachstum im Fenster (2 → 20 Objekte in
etwa 11 h) erwarten lässt. Die Budding Ratio ist damit eine Aussage über die
ersten Stunden eines Laufs, nicht über den ganzen Lauf. `µ_area` fittet
außerdem nur noch Tracks mit ≥ 10 Frames (`AREA_GROWTH_MIN_FRAMES`); der alte
Default von 2 Frames fittete überwiegend Fragmente.

Die Lineage-Abbildung ist damit `21_budding_rate_vs_period_<osc_type>.pdf`:
Buds je Mutter-Stunde im Fenster (`21_budding_rate_per_chip.csv`, aus
`budding_rate_per_h` in `21_budding_ratio_per_experiment.csv`), ein Chip pro
Periode mit seinen eigenen Kontrollen, Kammer-Fehlerbalken, Bracket-Score und
Spearman über Chips — dieselbe Logik wie `13_endpoint_vs_period_*`.
Zeitnormiert, weil das Fenster je Kammer verschieden lang ist (in den echten
Daten 27 bis 133 Frames; 211 von 565 Kammern werden nie voll). Für die
statischen Chips: `21_budding_rate_vs_medium.pdf`.

### Größenkriterium: angespülte Zellen sind keine Knospen

Blastokonidien werden laufend aus anderen Kammern angespült und tauchen
„neu“ neben sitzenden Zellen auf — für die räumliche Zuordnung sind sie von
einer Knospe nicht zu unterscheiden. Eine echte Knospe beginnt aber klein,
eine angespülte Zelle ist etwa so groß wie die Zelle, neben der sie landet.
`run_analysis.py` berechnet deshalb für jeden Kandidaten das Verhältnis
*Fläche beim ersten Auftreten / Fläche der zugeordneten Mutter* und verwirft
alles über einer Schwelle (`bud_size.py`). Die Schwelle kommt aus den Daten:
Antimodus der zweigipfligen Verteilung (Kerndichte auf log10, kleinste
Bandbreite mit genau zwei Gipfeln), **eine** Schwelle für alle Zweige. Ist die
Verteilung nicht zweigipflig, wird **kein** Größenfilter angewendet
(`BUD_MAX_AREA_FRACTION_FALLBACK = None`, Schwelle `inf` in der Tabelle), und
der Grund steht in der Spalte `source`. Auf den echten Daten ist genau das der
Fall: eine breite Mode um 0.4–0.5 ohne zweiten Gipfel. Das Kriterium bleibt
dort inaktiv, bis ein tragfähiges Unterscheidungsmerkmal gefunden ist — dafür
trägt `20_bud_size_at_appearance.csv` Diagnosespalten: gerade beendeter Track
an derselben Stelle (`ended_track_*`, Tracking-Bruch statt Knospe — im
manuellen QC die häufigste Korrektur), Wachstum des Kandidaten danach,
Kontaktverhältnis, Bewegung im nächsten Frame, Flächenbilanz der Mutter (nur
Kontrolle: die Muttermaske ändert sich beim ersten Segmentieren einer Knospe
praktisch nicht).

| Datei (direkt in `analysis_output/`) | Inhalt |
| --- | --- |
| `20_bud_size_at_appearance.pdf` | Verteilung des Verhältnisses mit Dichte, beiden Gipfeln und Schwelle; rechts je Stamm/osc_type |
| `20_bud_size_threshold.csv` | angewandte Schwelle (`scope = global`), Quelle, Gipfel, Antimodus; Kontrollzeilen je Gruppe |
| `20_bud_size_at_appearance.csv` | ein Kandidat pro Zeile (`bud_area`, `mother_area`, `bud_area_fraction`) |

`validate_lineage.py` wendet dieselbe Schwelle an (liest sie aus
`20_bud_size_threshold.csv`; `--bud-size-threshold` überschreibt, `inf`
schaltet ab): zu große Kandidaten zählen nicht in den Nenner der
Erkennungsrate und stehen pro Kammer in `n_rejected_by_size`. Im QC-Overlay
(`plot_qc_lineage_overlay.py`) erscheinen sie als Kandidaten ohne Mutter.

**1. Quantitativ, über alle Kammern:**

```bash
cd analyse_pipeline
python validate_lineage.py          # Pfade kommen aus config.py
```

Erzeugt in `analysis_output/lineage_validation/`:

| Datei | Frage, die sie beantwortet |
| --- | --- |
| `lv_01_d_over_r_distribution.pdf` | Lagen die Buds komfortabel im Suchradius, oder hat die Toleranz sie gerade noch hereingeholt? |
| `lv_02_detection_rate.pdf` + `_per_chamber.csv` | Ist die Erkennungsrate über die Bedingungen konstant? Nenner: Kandidaten, die das Größenkriterium bestehen |
| `lv_03_detection_rate_kruskal.csv` | Kruskal-Wallis dazu: p < 0.05 = Erkennung mit der Bedingung konfundiert; `spearman_rho_vs_period` = läuft die Erkennung *monoton* mit der Periode (die Richtung, die einen Trend vortäuscht)? |
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

* `lineage.py` leitet Mutter/Bud aus räumlicher Nähe, Tracklänge und der Größe
  beim ersten Auftreten ab; es gibt **kein** echtes Lineage-Tracking aus der
  Bildverarbeitung, und der Tracker fragmentiert (~12 % neue IDs je Objekt und
  Frame). Die Heuristik läuft deshalb nur im Sparse-Phase-Fenster je Kammer
  (`relink.py`): alle Lineage-Ergebnisse beschreiben die ersten Stunden eines
  Laufs. Vor jeder Publikation mit `inspect_lineage.py` gegen echte
  QC-Overlays kalibrieren.
* `tolerance_px` ist ein Pixel-Wert und damit von Kamera/Optik abhängig.
