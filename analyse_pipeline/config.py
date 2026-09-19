"""
config.py
=========
EINE zentrale Stelle für alle Pfade und Analyse-Parameter.

Vorher lagen dieselben Einstellungen doppelt in run_analysis.py UND in
inspect_lineage.py - mit unterschiedlichen Werten (anderer DATA_ROOT, anderer
CACHE_PATH, mother_min_frames=10 vs. 20). Dadurch hat das Kalibrierungs-Skript
faktisch einen anderen Datenstand mit anderen Schwellen ausgewertet als die
Pipeline selbst. Beide importieren jetzt von hier.

PFADE ANPASSEN
--------------
Drei Möglichkeiten, in dieser Reihenfolge (die erste, die gesetzt ist, gewinnt):

  1. Umgebungsvariablen (praktisch auf dem Cluster / in Skripten):
         Linux/macOS:  export AUREO_DATA_ROOT=/prj/microfluidic/ma_mimorde/Data
         Windows CMD:  set AUREO_DATA_ROOT=D:\\...\\Data
         PowerShell:   $env:AUREO_DATA_ROOT = "D:\\...\\Data"
     Zusätzlich optional: AUREO_OUTPUT_DIR
  2. ACTIVE_PRESET unten auf einen Eintrag aus DATA_ROOT_PRESETS setzen.
  3. Einen Preset-Pfad direkt editieren.

PARAMETER-ABWEICHUNGEN
----------------------
Einige Werte hier weichen bewusst von dem ab, was die jeweiligen Modul-
Docstrings als Standard dokumentieren (siehe DOCUMENTED_DEFAULTS). Sie werden
NICHT stillschweigend übernommen: log_active_configuration() schreibt beim
Start jeder Auswertung eine deutliche Warnung mit aktivem Wert, dokumentiertem
Wert und Begründung. Damit bleibt die Abweichung sichtbar, ohne dass hier
ungefragt an euren Ergebnissen gedreht wird.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

from lineage import FluxChannelConfig, LineageParams

logger = logging.getLogger(__name__)


# ==============================================================================
# 1. Pfade
# ==============================================================================

DATA_ROOT_PRESETS: dict[str, Path] = {
    "local": Path(r"D:\Studium\6_SoSe2026\aureo\Daten_verarbeitung\lokal_test\Data"),
    "cluster": Path("/prj/microfluidic/ma_mimorde/Data"),
}

# Welcher Preset gilt, wenn AUREO_DATA_ROOT nicht gesetzt ist.
ACTIVE_PRESET = "local"


def _from_env_or(var: str, fallback: Path) -> Path:
    value = os.environ.get(var)
    return Path(value) if value else fallback


DATA_ROOT: Path = _from_env_or("AUREO_DATA_ROOT", DATA_ROOT_PRESETS[ACTIVE_PRESET])
OUTPUT_DIR: Path = _from_env_or("AUREO_OUTPUT_DIR", DATA_ROOT.parent / "analysis_output")

# Der Parquet-Cache liegt bewusst NEBEN analysis_output, nicht darin: so
# überlebt er ein Löschen des Output-Ordners (Neuauswertung ohne Neu-Einlesen).
CACHE_PATH: Path = OUTPUT_DIR.parent / "combined_results_cache.parquet"
QC_EXCLUSIONS_PATH: Path = OUTPUT_DIR / "qc_exclusions.csv"
OUTPUT_DIR_STATIC: Path = OUTPUT_DIR / "static"
OUTPUT_DIR_PKO: Path = OUTPUT_DIR / "pko"

FORCE_RELOAD = False  # auf True setzen, wenn neue Rohdaten dazugekommen sind


# ==============================================================================
# 2. Zeit- & Achsen-Konfiguration
# ==============================================================================

MIN_PER_FRAME = 10.0  # Minuten pro Frame

# WICHTIG ZUR SPALTE 'osc_freq': sie enthaelt trotz ihres Namens die PERIODE
# der Feast/Famine-Zyklen in MINUTEN (0.75 ... 24), keine Frequenz.
#
# Daraus folgt eine Randbedingung, die in JEDE Methodenbeschreibung gehoert:
# bei MIN_PER_FRAME=10 liegt die kuerzeste aufloesbare Periode (Nyquist) bei
# 20 min. ALLE Oszillationsbedingungen liegen also am oder unter dem
# Abtastlimit - ein einzelner Zyklus ist grundsaetzlich nicht beobachtbar,
# und eine scheinbare Periodizitaet in den Sensor-Zeitreihen waere ein
# Alias-Artefakt, nicht der Medienwechsel.
#
# Die Oszillation ist damit eine BEHANDLUNG, keine Messgroesse: bei gleicher
# Gesamtdauer erfahren die Bedingungen ~1600 (0.75 min) bis ~50 (24 min)
# Zyklen in ~20 h Oszillation (Aufnahmen: ~133 Frames a 10 min = ~22 h) - ein
# 32-facher Dosisbereich bei identischer Gesamt-Feast-
# und Gesamt-Famine-Zeit. Interpretiert werden kann entsprechend nur die
# KUMULATIVE Wirkung ueber Stunden, nicht der Verlauf innerhalb eines Zyklus.
OSC_FREQ_IS_PERIOD_IN_MINUTES = True

# Reihenfolge der Oszillationsfrequenzen auf der x-Achse.
# resolve_x_order() in run_analysis.py hängt Werte, die hier fehlen, hinten an
# (mit Warnung), statt sie stillschweigend aus den Plots zu werfen.
FREQ_ORDER = ["0.75", "1.5", "3", "6", "12", "24"]

# Statische (nicht-oszillierende) Daten liegen unter
# Data/<Biosensor>/static/static_<medium>/... und werden komplett getrennt
# ausgewertet (eigener Output-Ordner OUTPUT_DIR_STATIC). Dort steht in
# 'osc_freq' statt einer Frequenz die Vergleichsgruppe St.omlp vs. St.ypd.
# 'osc_freq' benennt bei den statischen Daten die CHIP-FAMILIE, nicht das
# Medium: W109 hat pro Medium einen Ordner (static_omlp/static_ypd), W65 einen
# gemeinsamen (static_W65) mit beiden Medien in 'condition'. Das Medium wird
# deshalb aus 'condition' abgeleitet (St.omlp -> omlp) und ist die x-Achse
# aller statischen Abbildungen; die Chip-Familie ist die Facette.
# Siehe experiment_units.add_experiment_units().
STATIC_CHIP_LABELS: dict[str, str] = {
    "static_omlp": "W109",
    "static_ypd": "W109",
    "static_W65": "W65",
}
STATIC_MEDIUM_PREFIX = "St."
STATIC_MEDIUM_ORDER = ["omlp", "ypd"]

# --- PKO: dritter, unabhaengiger Zweig ----------------------------------------
# Der PKO-Stamm produziert kein Pullulan und dient der Pruefung, ob die
# Kontrollen sich ohne Exopolysaccharid korrekt verhalten (Clogging-Hypothese,
# siehe pko_comparison.py).
#
# PKO liegt als eigener Ordner auf der BIOSENSOR-Ebene
# (Data/PKO/<osc_type>/<periode>/03_results/...) und landet damit in der Spalte
# 'biosensor' - eine echte 'strain'-Spalte liefert Combined_Results nicht
# (siehe PANEL_A_GROUP_COL unten). Ohne die Abtrennung in run_analysis.py wuerde
# PKO deshalb als zusaetzliche Farbe in JEDEN bestehenden Oszillations-Plot
# laufen (PANEL_A_GROUP_COL ist dort auch color_col) und als zusaetzliches
# Violin in Panel A - die vorhandenen Ergebnisse wuerden sich also aendern.
# Genau das soll nicht passieren: PKO bekommt einen eigenen Kontext und einen
# eigenen Output-Ordner, exakt wie die statischen Daten.
PKO_BIOSENSOR_NAME = "PKO"

# Oszillationen starten nach einer zweistündigen Kontrollphase. ANGABE IN
# MINUTEN - queen_controls._minutes_to_hours() rechnet auf die 'time_h'-Achse um.
OSCILLATION_START_MIN = 120.0


# ==============================================================================
# 3. Plot-Gruppierung
# ==============================================================================

# Panel A (Violin-Plots mit Signifikanztests): welche Spalte definiert die
# Gruppen auf der x-Achse (im Paper: der Hefe-Stamm)? Auf "strain" umstellen,
# sobald eine echte Stamm-Spalte aus der Bildverarbeitung kommt.
PANEL_A_GROUP_COL = "biosensor"
PANEL_A_FACET_COL = "osc_type"

# Panel A fuer die STATISCHEN Daten: dort ist 'biosensor' die falsche Gruppe.
# plot_panel_a() nutzt nur group_col und facet_col und sieht 'osc_freq' nie -
# bei den statischen Daten steht die Vergleichsgruppe (static_omlp/static_ypd)
# aber genau dort. Mit group_col='biosensor' landen deshalb BEIDE Medien in
# EINEM Violin, und die Abbildung kann die Frage "komplexes vs. minimales
# Medium" gar nicht beantworten. Der statische Kontext setzt daher 'osc_freq'.
PANEL_A_GROUP_COL_STATIC = "medium"

# --- Kumulativer Endzustand (endpoint_trends.py, Schritt 13) -----------------
# Anteil der Frames am ENDE jeder Kammer, der den "Endzustand" bildet. Relativ
# zur jeweiligen Kammer, nicht absolut - Kammern koennen unterschiedlich lang
# aufgenommen sein.
ENDPOINT_LAST_FRACTION = 0.25
# Spalten, fuer die der Endzustand gegen die Periode gestellt wird. Die zur
# Laufzeit erkannten ratio_*-Spalten kommen in pipeline_steps.py dazu - erst
# dadurch bekommen die Sensor-Daten ueberhaupt eine kumulative Auswertung
# (50_summary_statistics.csv sieht nur intensity_cols, nie die Ratios).
ENDPOINT_VALUE_COLS = ["area", "eccentricity"]
# Statische Daten: ypd waechst ueber (bis zu 4000 Tracks pro Kammer) und die
# W109-ypd-Aufnahmen wurden bei 85 Frames abgebrochen. Ein relatives Endfenster
# ("letzte 25 % der Frames") vergleicht dann eine ueberwachsene ypd-Kammer bei
# 14 h mit einer normalen omlp-Kammer bei 22 h. Deshalb wird fuer den statischen
# Zweig die Saettigung PRO KAMMER aus dem Zellzahl-Verlauf bestimmt, und das
# Endfenster endet an der fruehesten Saettigung ueber alle Kammern.
# endpoint_trends.detect_saturation_frame(): Saettigung = erster Frame, ab dem
# die geglaettete Zellzahl >= STATIC_SATURATION_LEVEL x ihres Maximums bleibt.
STATIC_SATURATION_LEVEL = 0.90
STATIC_SATURATION_SMOOTH_FRAMES = 5
# Nur Kammern, die ueberhaupt gewachsen sind (max/Start >= dieser Faktor),
# koennen saettigen. Eine flache Kammer (omlp) liegt sonst von Anfang an bei
# 90 % ihres Maximums und wuerde das Fenster auf die ersten Stunden ziehen.
STATIC_SATURATION_MIN_GROWTH = 1.5

# Zusatz-Zeile unter den Kontroll-Ticks in plot_sensor_control_comparison(),
# PRO OSZILLATIONSTYP. Die Konzentrationen gelten nur für den Oszillationstyp,
# mit dem sie gemessen wurden - eine Glucose-Konzentration unter einer
# pH-Kontrolle wäre schlicht falsch, deshalb steht hier kein globaler Default.
CONTROL_CONCENTRATION_LABELS: dict[str, dict[str, str]] = {
    "Glc": {"NegCtrl": "0 g/L", "PosCtrl": "50 g/L"},
}


# ==============================================================================
# 4. Lineage / Wachstum / Robustness
# ==============================================================================

# Mutter/Bud-Klassifikation - siehe lineage.py für Details & Kalibrierungshinweise.
# WICHTIG: tolerance_px ist in Pixel eurer Kamera/Optik - unbedingt an echten
# Bildern/QC-Overlays kalibrieren (inspect_lineage.py), bevor den Ergebnissen
# für eine Publikation vertraut wird.
LINEAGE_PARAMS = LineageParams(
    mother_min_frames=10,
    bud_max_frames=7,
    established_min_frames=3,
    tolerance_px=30.0,
)

# Groessenkriterium der Mutter/Bud-Heuristik (bud_size.py): eine neu
# auftauchende Zelle zaehlt nur als Knospe, wenn ihre Flaeche beim ersten
# Auftreten hoechstens diesen Anteil der Mutterflaeche hat; alles darueber
# ist eine angespuelte Blastokonidie. Die Schwelle wird aus den Daten
# abgeleitet (Antimodus der zweigipfligen Verteilung von bud_area /
# mother_area, EINE Schwelle fuer alle Zweige); der Rueckfallwert greift nur,
# wenn die Verteilung nicht zweigipflig ist, zu wenige Kandidaten hat oder der
# Antimodus ausserhalb des plausiblen Bereichs liegt. Welcher Fall eintrat,
# steht in analysis_output/20_bud_size_threshold.csv (Spalte 'source').
BUD_MAX_AREA_FRACTION_FALLBACK = 0.5
BUD_SIZE_PLAUSIBLE_RANGE = (0.15, 0.9)

# Detail-Trajektorien "stabiler" Mütter (Schritt 92, Anhang):
# coverage = n_frames / (frame_max - frame_min + 1), siehe lineage.identify_mothers().
STABLE_MOTHER_MIN_COVERAGE = 0.15
STABLE_MOTHER_GROUP_COLS = ["biosensor", "osc_type", "osc_freq", "condition"]
STABLE_MOTHER_BASE_VALUE_COLS = ["area"]

# Optional: physiologischer Flux zum Budding-Zeitpunkt (ratiometrischer Sensor).
# None = kein Flux berechnet.
FLUX_CONFIG: FluxChannelConfig | None = None
# z.B. FluxChannelConfig(channel_a="mean_mTurqouise", channel_b="mean_RFP", min_denominator=1.0)

# Spezifische Wachstumsrate (Eq. 2): µ-Werte über dieser Schwelle [h^-1] werden
# als Artefakt markiert (mu_is_artefact=True), aber NICHT gelöscht.
MU_MAX_THRESHOLD = 10.0

# Robustness R(t)/R(p) (siehe robustness.py): für welche Spalten berechnen?
# Die zur Laufzeit erkannten ratio_*-Spalten kommen in run_analysis.py dazu.
# 'budding_ratio' wird separat aus der Zeitreihe behandelt (Schritt 22).
ROBUSTNESS_VALUE_COLS = ["area", "eccentricity"]


# ==============================================================================
# 5. Abweichungen sichtbar machen
# ==============================================================================

# name -> (aktiver Wert, dokumentierter Wert, Begründung/Fundstelle)
DOCUMENTED_DEFAULTS: dict[str, tuple[object, object, str]] = {
    "MU_MAX_THRESHOLD": (
        MU_MAX_THRESHOLD, 0.6,
        "growth_rate.compute_specific_growth_rate() Default = 0.6 h^-1 (Paper Methods: "
        "maximales realistisches µ für Hefe ca. 0.5 h^-1). Beim aktiven Wert greift der "
        "Artefakt-Filter praktisch nicht mehr - mu_is_artefact bleibt fast überall False.",
    ),
    "LINEAGE_PARAMS.mother_min_frames": (
        LINEAGE_PARAMS.mother_min_frames, 20,
        "lineage.LineageParams Default = 20, dort begründet über Xiao 2025 (µ/1 = 0.29 -> "
        "20 Frames decken mindestens eine Reproduktion ab). Niedrigere Werte lassen mehr, "
        "aber kürzer beobachtete Tracks als Mutterzelle gelten.",
    ),
    "LINEAGE_PARAMS.bud_max_frames": (
        LINEAGE_PARAMS.bud_max_frames, 5,
        "lineage.LineageParams Default = 5. Reine Report-Schwelle (bud_was_washed_out), "
        "beeinflusst die Event-Erkennung nicht.",
    ),
}

# Methodische Eigenheiten, die kein einzelner Zahlenwert sind, aber beim Lesen
# der Ergebnistabellen bekannt sein müssen.
METHOD_CAVEATS: list[str] = [
    "VERSUCHSEINHEITEN: pro (Stamm, osc_type, Periode) gibt es EINEN Chip aus EINER Vorkultur; "
    "'replicate' ist ein Array-Index, kein Replikat. Innerhalb einer Oszillationsbedingung gibt "
    "es keine biologische Replikation (n = 1 Chip). Fehlerbalken dort sind Kammer-Fehlerbalken - "
    "siehe Spalte 'error_unit' in den Aggregaten. Statisch: jedes 'replicate' ist ein Chip.",
    "summarise_growth_rate() und summarise_area_growth() aggregieren ueber EINZELNE Zellen/"
    "Intervalle (11_*_summary.csv, 12_*_summary_*.csv): sd_mu/sd_mu_area sind Streuung ueber Zellen.",
    "Spearman gegen die Periode laeuft auf Chip-Mittelwerten (n = Zahl der Perioden) und ist bei "
    "n <= 6 eine Effektstaerke, kein Test. Die Staemme werden nicht als Replikate gepoolt.",
    "Robustness R ist RELATIV (Normalisierung ueber den uebergebenen Datensatz): R-Werte aus "
    "analysis_output/, static/, pko/ und no_qc/ sind nicht gegeneinander lesbar.",
    "'area' in den Zelltabellen ist die rohe Cellpose-Flaeche in px². Nur area_growth.py rechnet "
    "intern mit PX_TO_UM2 (1 µm = 13.63 px) in µm² um.",
    "Kammerposition und Bedingung sind durch die Chip-Verdrahtung konfundiert (A1/A2 Feast, "
    "A13/A14 Famine, A3-A12 Wechsel).",
    "KNOSPEN-GROESSENKRITERIUM: eine neu auftauchende Zelle zaehlt nur als Knospe, wenn ihre "
    "Flaeche beim ersten Auftreten hoechstens Schwelle x Mutterflaeche ist; groessere gelten als "
    "angespuelte Blastokonidien. EINE Schwelle fuer alle Zweige, aus den Daten (Antimodus), "
    "Rueckfall BUD_MAX_AREA_FRACTION_FALLBACK - siehe 20_bud_size_threshold.csv.",
]


def _is_number(value: object) -> bool:
    try:
        float(value)  # type: ignore[arg-type]
        return True
    except (TypeError, ValueError):
        return False


def log_active_configuration() -> None:
    """Schreibt die aktive Konfiguration ins Log, inkl. aller Abweichungen.

    Wird von run_analysis.py und inspect_lineage.py beim Start aufgerufen, damit
    in JEDEM Auswertungs-Log nachvollziehbar steht, mit welchen Pfaden und
    Schwellen die Ergebnisse entstanden sind.
    """
    logger.info("=== Aktive Konfiguration (config.py) ===")
    logger.info("  DATA_ROOT:        %s", DATA_ROOT)
    logger.info("  OUTPUT_DIR:       %s", OUTPUT_DIR)
    logger.info("  CACHE_PATH:       %s", CACHE_PATH)
    logger.info("  MIN_PER_FRAME:    %.1f min", MIN_PER_FRAME)
    if OSC_FREQ_IS_PERIOD_IN_MINUTES:
        nyquist_min = 2 * MIN_PER_FRAME
        periods = [p for p in (float(f) for f in FREQ_ORDER if _is_number(f))]
        unresolved = [p for p in periods if p <= nyquist_min]
        if unresolved:
            logger.warning(
                "ABTASTUNG: 'osc_freq' ist die Periode in Minuten. Bei %.0f min/Frame liegt die "
                "kuerzeste aufloesbare Periode bei %.0f min - %d von %d Bedingungen (%s) liegen "
                "darunter. Einzelne Zyklen sind NICHT beobachtbar; scheinbare Periodizitaet in "
                "den Sensor-Zeitreihen ist ein Alias-Artefakt. Interpretierbar ist nur die "
                "kumulative Wirkung ueber Stunden, nicht der Zyklusverlauf.",
                MIN_PER_FRAME, nyquist_min, len(unresolved), len(periods),
                ", ".join(f"{p:g}" for p in unresolved),
            )
    logger.info("  LINEAGE_PARAMS:   %s", LINEAGE_PARAMS)
    logger.info("  Knospen-Groessenkriterium: Schwelle aus den Daten, Rueckfall %.2f, plausibel %s",
                BUD_MAX_AREA_FRACTION_FALLBACK, BUD_SIZE_PLAUSIBLE_RANGE)
    logger.info("  MU_MAX_THRESHOLD: %s h^-1", MU_MAX_THRESHOLD)

    for name, (active, documented, why) in DOCUMENTED_DEFAULTS.items():
        if active != documented:
            logger.warning(
                "KONFIGURATION WEICHT AB: %s = %s (dokumentierter Standard: %s). %s",
                name, active, documented, why,
            )

    for caveat in METHOD_CAVEATS:
        logger.warning("METHODEN-HINWEIS: %s", caveat)
