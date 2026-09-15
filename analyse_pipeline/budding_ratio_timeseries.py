"""
budding_ratio_timeseries.py
============================
Budding Ratio als ZEITREIHE, exakt nach Eq. 3 aus Blöbaum et al. 2024
(Microbial Cell Factories):

    Budding Ratio(t) = (Anzahl neuer Buds zum Zeitpunkt t) / (Anzahl Zellen zum Zeitpunkt t-1)

WICHTIG - Unterschied zu lineage.compute_budding_ratio():
Jenes Modul liefert "Buds pro Mutterzelle, summiert über die GESAMTE
Beobachtungsdauer" (eine Zeile pro Mutter) - das ist die Paper-Panel-A-
artige Violin-Darstellung. DIESES Modul liefert die Budding Ratio als
TRAJEKTORIE pro Kammer (eine Zeile pro Kammer x Zeitpunkt) - das ist die
Größe, die im Paper in Fig. 3a (Linienplot über die Zeit) und für die
R(t)-Berechnung (siehe robustness.py) verwendet wird.

Begründung für "Zellen zum Zeitpunkt t-1" statt "Mütter": im Paper wird
explizit die GESAMTE Zellzahl der Kammer im vorherigen Zeitpunkt verwendet
(nicht nur Mutterzellen) - das ist wichtig, weil die Budding Ratio so auch
funktioniert, wenn die Kammer voll ist und Zellen herausgespült werden
(siehe Paper-Begründung für Eq. 3).

WICHTIG zu Aufbau (Abweichung vom Paper):
Im Paper ist die Kontrollbedingung 'Control' = durchgehend Feast-Medium.
Bei euch ist die Negativkontrolle dagegen durchgehend STARVATION. Dieses
Modul behandelt 'Kontrolle' nicht inhaltlich (es zählt nur Buds/Zellen pro
Kammer-Zeitpunkt, unabhängig vom Medium), aber WENN ihr eure Negativkontrolle
mit der 'Control'-Spur aus dem Paper vergleicht, denkt daran, dass die
Interpretation (z.B. 'niedrige Budding Ratio = Stress' vs. '= Normalzustand')
bei euch ggf. umgekehrt zum Paper ist.
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def compute_budding_ratio_timeseries(
    cells: pd.DataFrame,
    lineage_events: pd.DataFrame,
) -> pd.DataFrame:
    """
    Berechnet die Budding Ratio als Zeitreihe pro Kammer (exp_id), exakt
    nach Eq. 3: Budding Ratio(t) = n_buds(t) / n_cells(t-1).

    Parameters
    ----------
    cells : voller Zelldatensatz (nach QC), eine Zeile pro Zelle x Frame,
            mit Spalten exp_id, frame, track_id (und den üblichen Hierarchie-
            Metadaten biosensor/osc_type/osc_freq/condition/replicate/chamber)
    lineage_events : Ergebnis von lineage.classify_mother_bud() (eine Zeile
            pro erkanntem Budding-Event, mit Spalte budding_frame)

    Returns
    -------
    DataFrame mit einer Zeile pro exp_id x frame:
        exp_id, frame, n_cells, n_buds, budding_ratio,
        (+ Hierarchie-Metadaten falls in cells vorhanden)

    HINWEIS zu frame=0 (erster Frame jeder Kammer): hier existiert kein
    't-1', daher wird n_cells(t-1) für frame=0 auf NaN gesetzt und
    budding_ratio ist dort NaN (nicht 0) - andernfalls würde eine künstliche
    Diskontinuität am Anfang jeder Trajektorie entstehen.
    """
    required = {"exp_id", "frame", "track_id"}
    missing = required - set(cells.columns)
    if missing:
        raise ValueError(f"compute_budding_ratio_timeseries() fehlen Spalten in 'cells': {missing}")

    # n_cells(t): Anzahl distinkter Tracks pro Kammer x Frame
    n_cells = (
        cells.groupby(["exp_id", "frame"])["track_id"]
        .nunique()
        .rename("n_cells")
        .reset_index()
    )

    # n_buds(t): Anzahl Budding-Events, deren budding_frame == t, pro Kammer
    if lineage_events.empty or "budding_frame" not in lineage_events.columns:
        logger.warning(
            "lineage_events ist leer oder hat keine 'budding_frame' Spalte - "
            "n_buds wird für alle Zeitpunkte 0 gesetzt."
        )
        n_buds = pd.DataFrame(columns=["exp_id", "frame", "n_buds"])
    else:
        n_buds = (
            lineage_events.groupby(["exp_id", "budding_frame"])
            .size()
            .rename("n_buds")
            .reset_index()
            .rename(columns={"budding_frame": "frame"})
        )

    # Vollständiges Frame-Gitter pro Kammer, damit Frames ohne Buds nicht
    # fehlen (sie sollen n_buds=0 bekommen, nicht einfach aus der Tabelle
    # verschwinden - sonst entstehen Lücken in der Zeitreihe)
    out = n_cells.merge(n_buds, on=["exp_id", "frame"], how="left")
    out["n_buds"] = out["n_buds"].fillna(0).astype(int)

    out = out.sort_values(["exp_id", "frame"])
    out["n_cells_prev"] = out.groupby("exp_id")["n_cells"].shift(1)

    with np.errstate(divide="ignore", invalid="ignore"):
        out["budding_ratio"] = out["n_buds"] / out["n_cells_prev"]

    # frame=0 jeder Kammer (kein t-1 vorhanden) explizit auf NaN setzen,
    # statt eines durch Division entstandenen inf/NaN-Mischmaschs
    is_first_frame = out["n_cells_prev"].isna()
    out.loc[is_first_frame, "budding_ratio"] = np.nan

    meta_cols = [c for c in ["biosensor", "osc_type", "osc_freq", "condition", "replicate", "chamber"]
                 if c in cells.columns]
    if meta_cols:
        meta = cells[["exp_id"] + meta_cols].drop_duplicates("exp_id")
        out = out.merge(meta, on="exp_id", how="left")

    logger.info(
        "Budding Ratio Zeitreihe berechnet: %d Kammern, %d Zeitpunkte gesamt, "
        "%d Frames ohne gültige Budding Ratio (erster Frame je Kammer).",
        out["exp_id"].nunique(), len(out), int(is_first_frame.sum()),
    )

    return out.reset_index(drop=True)


def aggregate_budding_ratio_over_replicates(
    timeseries: pd.DataFrame,
    group_cols: Optional[list[str]] = None,
) -> pd.DataFrame:
    """
    Mittelt die Budding-Ratio-Zeitreihe über Replikat-Kammern hinweg, analog
    zu Fig. 3a im Paper (Linienplot mit Fehlerbalken = SD über Triplikate).

    Parameters
    ----------
    timeseries : Ergebnis von compute_budding_ratio_timeseries()
    group_cols : Spalten, die EINE Bedingung definieren (alles außer
                 'replicate'/'chamber'). Standard: biosensor, osc_type,
                 osc_freq, condition, frame - falls vorhanden.

    Returns
    -------
    DataFrame mit mean_budding_ratio, sd_budding_ratio, n_replicates pro
    Bedingung x Zeitpunkt.
    """
    if group_cols is None:
        group_cols = [c for c in ["biosensor", "osc_type", "osc_freq", "condition", "frame"]
                      if c in timeseries.columns]

    agg = (
        timeseries.dropna(subset=["budding_ratio"])
        .groupby(group_cols)["budding_ratio"]
        .agg(mean_budding_ratio="mean", sd_budding_ratio="std", n_replicates="count")
        .reset_index()
    )
    return agg
