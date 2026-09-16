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
# Gesamtdauer erfahren die Bedingungen ~800 (0.75 min) bis ~25 (24 min)
# Zyklen in 10 h - ein 32-facher Dosisbereich bei identischer Gesamt-Feast-
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
STATIC_ORDER = ["static_omlp", "static_ypd"]

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
    "summarise_growth_rate() und summarise_area_growth() aggregieren über EINZELNE "
    "Zellen/Intervalle, nicht erst pro Replikat. sd_mu/sd_mu_area und n_values in "
    "11_*_summary.csv und 12_*_summary_*.csv beschreiben daher die Streuung über Zellen "
    "(Pseudoreplikation), nicht über biologische Replikate. "
    "aggregate_robustness_over_replicates() und analysis._aggregate_over_replicates() "
    "mitteln dagegen korrekt zweistufig.",
    "Robustness R ist eine RELATIVE Größe: der Normalisierungsfaktor m wird über den "
    "GESAMTEN übergebenen Datensatz gebildet. R-Werte aus Läufen mit unterschiedlichem "
    "Datenumfang sind nicht miteinander vergleichbar (siehe robustness.py).",
    "'area' in den Zelltabellen ist die rohe Cellpose-Fläche in px². Nur area_growth.py "
    "rechnet intern mit PX_TO_UM2 (1 µm = 13.63 px) in µm² um.",
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
    logger.info("  MU_MAX_THRESHOLD: %s h^-1", MU_MAX_THRESHOLD)

    for name, (active, documented, why) in DOCUMENTED_DEFAULTS.items():
        if active != documented:
            logger.warning(
                "KONFIGURATION WEICHT AB: %s = %s (dokumentierter Standard: %s). %s",
                name, active, documented, why,
            )

    for caveat in METHOD_CAVEATS:
        logger.warning("METHODEN-HINWEIS: %s", caveat)
