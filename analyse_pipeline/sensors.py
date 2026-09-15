"""
sensors.py
==========
Berechnung ratiometrischer Biosensor-Werte (Verhältnis zweier Kanäle) für
Sensoren, deren Signal NICHT direkt als absolute Intensität interpretierbar
ist (z.B. FRET-basierte oder Ratio-Sensoren), im Gegensatz zu intensiometrischen
Sensoren, deren einzelner Kanal direkt aussagekräftig ist.

Da die Kanal-Namen je Biosensor unterschiedlich sind (z.B. 'mean_488'/'mean_405'
für Sensor X, 'mean_GFP'/'mean_RFP' für Sensor Y), wird die Zuordnung über
eine explizite Konfigurationstabelle gepflegt - SENSOR_CONFIG unten.

WICHTIG: Diese Konfiguration ist der einzige Ort, an dem ihr pro Biosensor
festlegt, ob und wie ein Ratio gebildet wird. Intensiometrische Sensoren
brauchen hier KEINEN Eintrag - sie werden einfach direkt verwendet
(siehe analysis.py: find_intensity_columns()).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


@dataclass
class RatioSensorConfig:
    """
    Definiert ein Ratio-Kanal-Paar für EINEN ratiometrischen Biosensor.

    channel_a / channel_b : Spaltennamen der beiden Kanäle in Combined_Results.csv
    ratio_name : Name der neu erzeugten Spalte (z.B. 'ratio_iGlucoSnFR')
    direction : 'a_over_b' (Standard) oder 'b_over_a' - je nachdem, welcher
                Kanal im Zähler stehen soll (siehe Sensor-Dokumentation/Paper,
                aus dem ihr euch orientiert, für die übliche Konvention).
    min_denominator : Nenner-Werte unterhalb dieser Schwelle werden NICHT
                durch geteilt (Ratio wird NaN statt eines extremen/instabilen
                Werts). MUSS an euer Kamera-Hintergrundrauschen angepasst
                werden - es gibt keinen sinnvollen universellen Standardwert.
                Praktischer Startpunkt: schaut euch die Verteilung eures
                Hintergrund-/Leerkanal-Signals an (z.B. außerhalb von Zellen)
                und setzt min_denominator etwas darüber.
    """
    channel_a: str
    channel_b: str
    ratio_name: str
    direction: str = "a_over_b"
    min_denominator: float = 1.0


# ==============================================================================
# HIER PRO RATIOMETRISCHEM BIOSENSOR EINTRAGEN
# Key = exakter Wert der 'biosensor' Spalte (aus dem Ordnernamen, siehe data_loading.py)
# Intensiometrische Sensoren NICHT eintragen - die werden direkt verwendet.
# ==============================================================================
SENSOR_CONFIG: dict[str, RatioSensorConfig] = {
    # Beispiel - bitte an eure echten Kanal-Spaltennamen und Biosensor-Ordnernamen anpassen:
    # "pHRed": RatioSensorConfig(channel_a="mean_488", channel_b="mean_405", ratio_name="ratio_pHRed"),
    # "iGlucoSnFR_ratio": RatioSensorConfig(channel_a="mean_GFP", channel_b="mean_RFP", ratio_name="ratio_iGlucoSnFR"),

    "BSG": RatioSensorConfig(channel_a="mean_mTurqouise", channel_b="mean_RFP", ratio_name="ratio_GlyRNA", direction="a_over_b", min_denominator=10.0),
    "BSO": RatioSensorConfig(channel_a="mean_YFP", channel_b="mean_RFP", ratio_name="ratio_OxPro", direction="a_over_b", min_denominator=20.0),
    "BSA": RatioSensorConfig(channel_a="mean_uvGFP", channel_b="mean_GFP", ratio_name="ratio_Queen-2m", direction="a_over_b", min_denominator=10.0),
    "BSPH": RatioSensorConfig(channel_a="mean_uvGFP", channel_b="mean_GFP", ratio_name="ratio_pHluorin", direction="a_over_b", min_denominator=1.5),
}


def compute_ratios(
    df: pd.DataFrame,
    config: dict[str, RatioSensorConfig] | None = None,
) -> pd.DataFrame:
    """
    Berechnet für jeden in `config` (Standard: SENSOR_CONFIG) hinterlegten
    ratiometrischen Biosensor eine neue Spalte mit dem Intensitäts-Verhältnis.

    Die Berechnung wird NUR auf die Zeilen angewendet, deren 'biosensor'
    Spalte exakt dem jeweiligen Config-Key entspricht - andere Biosensoren
    bleiben unberührt (z.B. bekommen intensiometrische Sensoren keine
    Ratio-Spalte, da sie nicht in der Config stehen).

    Division durch 0 oder sehr kleine Nenner wird zu NaN (nicht Inf), mit
    einer Warnung, wie viele Zeilen betroffen waren - das verhindert, dass
    ein paar Pixel-Rauschen-Ausreißer die ganze Skala eines Plots verzerren.

    Parameters
    ----------
    df : Datensatz mit Spalte 'biosensor' und den in config referenzierten Kanal-Spalten
    config : dict[biosensor_name -> RatioSensorConfig], Standard: SENSOR_CONFIG

    Returns
    -------
    DataFrame mit einer neuen Spalte pro konfiguriertem Sensor
    (Spaltenname = ratio_name), NaN für alle Zeilen anderer Biosensoren.
    """
    if config is None:
        config = SENSOR_CONFIG

    if not config:
        logger.info(
            "SENSOR_CONFIG ist leer - es werden keine Ratios berechnet. "
            "Falls ihr ratiometrische Biosensoren habt, in sensors.py SENSOR_CONFIG befüllen."
        )
        return df

    df = df.copy()

    if "biosensor" not in df.columns:
        logger.warning("Spalte 'biosensor' fehlt im DataFrame. Kann Ratios nicht zuordnen.")
        return df

    present_biosensors = set(df["biosensor"].dropna().unique())
    n_matched = 0

    for biosensor_name, cfg in config.items():
        # --- Check 1: Ist dieser Biosensor in den aktuellen Daten vorhanden?
        if biosensor_name not in present_biosensors:
            logger.info(
                "Sensor '%s': Überspringen (nicht in den aktuellen Daten vorhanden).", biosensor_name
            )
            continue
        # --- Check 2: Sind die benötigten Kanal-Spalten vorhanden?
        missing_cols = [col for col in (cfg.channel_a, cfg.channel_b) if col not in df.columns]
        if missing_cols:
            logger.warning(
                "Sensor '%s': Kanal-Spalten %s fehlen in den Daten - übersprungen.",
                biosensor_name, missing_cols
            )
            continue

        # --- Ab hier: Berechnung nur für den passenden Biosensor
        mask = df["biosensor"] == biosensor_name
        n_rows = mask.sum()

        if n_rows == 0:
            continue

        a = df.loc[mask, cfg.channel_a]
        b = df.loc[mask, cfg.channel_b]
        numerator, denominator = (a, b) if cfg.direction == "a_over_b" else (b, a)

        # Division durch (nahe) 0 / Hintergrundrauschen robust behandeln
        near_zero = denominator.abs() < cfg.min_denominator
        n_near_zero = int(near_zero.sum())
        if n_near_zero > 0:
            logger.warning(
                "Sensor '%s': %d von %d Zeilen haben einen Nenner unter min_denominator=%.4g "
                "- Ratio wird dort zu NaN statt eines extremen/instabilen Wertes gesetzt.",
                biosensor_name, n_near_zero, n_rows, cfg.min_denominator,
            )

        ratio = np.where(near_zero, np.nan, numerator / denominator.where(~near_zero, np.nan))

        if cfg.ratio_name not in df.columns:
            df[cfg.ratio_name] = np.nan
        df.loc[mask, cfg.ratio_name] = ratio
        n_matched += 1

        numerator_name = cfg.channel_a if cfg.direction == "a_over_b" else cfg.channel_b
        denominator_name = cfg.channel_b if cfg.direction == "a_over_b" else cfg.channel_a
        logger.info(
            "Sensor '%s': Ratio-Spalte '%s' berechnet (%s / %s, direction=%s) für %d Zeilen.",
            biosensor_name, cfg.ratio_name, numerator_name, denominator_name, cfg.direction, n_rows,
        )

    if n_matched == 0:
        logger.warning(
            "compute_ratios(): KEIN Eintrag aus SENSOR_CONFIG (%s) hat einen 'biosensor'-Wert "
            "in den Daten (%s) getroffen - es wurde keine einzige Ratio-Spalte erzeugt. "
            "Tippfehler in SENSOR_CONFIG-Keys oder im Ordnernamen?",
            sorted(config.keys()), sorted(present_biosensors),
        )

    return df
