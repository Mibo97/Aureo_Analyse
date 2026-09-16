"""
pko_comparison.py
=================
WT-gegen-PKO-Vergleich fuer die Clogging-Hypothese.

HINTERGRUND
-----------
Die Kontrollkammern (PosCtrl = durchgehend Feast, NegCtrl = durchgehend
Starvation) verhalten sich im Wildtyp nicht so, wie sie sollten. Die
Arbeitshypothese: Pullulan - das Exopolysaccharid, das der WT ausscheidet -
setzt die Chip-Strukturen zu und erzeugt dadurch unerwartete Stroemungs-
profile. Der PKO-Stamm produziert kein Pullulan. Verhalten sich SEINE
Kontrollen korrekt, stuetzt das die Clogging-Erklaerung.

WAS HIER NICHT GEHT, UND WARUM
------------------------------
Der naheliegende Test - test_control_consistency_across_freq() auf die
PKO-Daten - ist mit diesen Daten NICHT durchfuehrbar. PKO deckt nur EINE
Periode ab, der Kruskal-Wallis-Test ueber die osc_freq-Batches braucht aber
mindestens zwei; er liefert dort p = NaN. Schlimmer ist der Plot:
plot_control_consistency() zeichnet bei einer einzigen Periode einen
einzelnen Punkt pro Kontrollart, also eine trivial flache Linie, die wie
"PKO-Kontrollen sind konsistent" aussieht und keinerlei Information
enthaelt. Der PKO-Kontext schaltet die Kontroll-Konsistenz deshalb ab
(PipelineContext.run_control_consistency=False), statt eine irrefuehrende
Abbildung zu erzeugen.

Der Vergleich wandert damit von "ueber Batches innerhalb PKO" nach
"innerhalb des gemeinsamen Perioden-Batches, WT gegen PKO". Das ist kein
Notbehelf, sondern aus zwei Gruenden sogar der sauberere Schnitt:

  * PosCtrl ist durchgehend Feast, NegCtrl durchgehend Starvation - der
    Medienverlauf einer KONTROLLkammer haengt gar nicht an der Periode des
    Batches, in dem sie mitlief. Die Einschraenkung auf eine Periode kostet
    fuer den Kontrollvergleich also kaum etwas.
  * In den Kontrollkammern faellt am meisten Pullulan an: PosCtrl waechst
    10 h durch. Wenn Verstopfung das Problem ist, ist das der Ort dafuer.

DREI GROESSEN, KEINE DAVON LINEAGE-ABHAENGIG
--------------------------------------------
Bewusst ohne Mutter/Bud-Heuristik (lineage.py ist unvalidiert, siehe
validate_lineage.py), damit dieser Vergleich nicht an einer offenen Frage
haengt:

  1. Kammer-Uebereinstimmung: Streuung nominell identischer Kontrollkammern
     INNERHALB eines Replikats.
  2. Kontroll-Bracket: trennt PosCtrl von NegCtrl ueberhaupt? Cliff's Delta
     auf mu_area, pro Replikat gerechnet.
  3. Zell-Ausbeute: bricht die Zahl verfolgbarer Zellen pro Kammer ueber die
     Zeit weg?

STATISTIK: WAS DIESES DESIGN TRAEGT
-----------------------------------
PKO ist EIN Stamm. Ein Signifikanztest WT-gegen-PKO auf Stammebene hat n = 1
in einer Gruppe und ist nicht durchfuehrbar. Er wird hier bewusst NICHT
gerechnet - auch nicht ueber Replikate oder Kammern hinweg, denn das waere
Pseudoreplikation auf Stammebene und wuerde aus einem Stamm sechs machen.

Was dieses Design traegt, ist eine DESKRIPTIVE Aussage: die vier
WT-Biosensor-Staemme liefern vier voneinander unabhaengige WT-Werte, und der
PKO-Wert liegt innerhalb oder ausserhalb dieser Spanne. Die vier WT-Staemme
sind damit die interne Replikation der WT-Seite - stimmen sie untereinander
nicht ueberein, ist schon das Zusammenfassen zu "WT" falsch, und die
Abbildung zeigt genau das, statt es in einem Mittelwert zu verstecken.

CONFOUND, DER IN DIE DISKUSSION GEHOERT
---------------------------------------
PKO ist nicht "Wildtyp ohne Verstopfung", sondern eine Mutante mit
veraendertem Kohlenstofffluss und veraenderten Oberflaecheneigenschaften.
JEDER WT-PKO-Unterschied laesst sich auch direkt physiologisch erklaeren -
einen Unterschied zu finden ist daher noch kein Beleg fuer Verstopfung.

Unterscheidbar sind die beiden Erklaerungen nur ueber das MUSTER:

  hydraulisch  -> hohe Kammer-zu-Kammer-Streuung bei nominell identischen
                  Kammern, wegbrechende Zell-Ausbeute, kaputtes Bracket
  metabolisch  -> gleichmaessige Niveau-Verschiebung, Kammer-Ueberein-
                  stimmung bleibt erhalten

Deshalb berichtet dieses Modul Streuungen und Steigungen, nicht nur
Mittelwerte.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from growth_rate import classify_condition_type

logger = logging.getLogger(__name__)

CONTROL_ORDER = ["NegCtrl", "PosCtrl"]
CONTROL_COLORS = {"NegCtrl": "#4C78A8", "PosCtrl": "#E45756"}

STRAIN_GROUP_WT = "WT"
STRAIN_GROUP_PKO = "PKO"
STRAIN_GROUP_COLORS = {STRAIN_GROUP_WT: "#6E6E6E", STRAIN_GROUP_PKO: "#54A24B"}


def _minutes_to_hours(minutes: float) -> float:
    """*_min-Angaben sind MINUTEN, 'time_h' ist in STUNDEN.

    Gleiche Falle wie in queen_controls._minutes_to_hours(): ein direkter
    Vergleich liest sich als ">= 120 Stunden" und trifft stillschweigend
    nichts.
    """
    return minutes / 60.0


# ==============================================================================
# Vorbereitung
# ==============================================================================

def add_strain_group(df: pd.DataFrame, pko_biosensor: str) -> pd.DataFrame:
    """Ergaenzt 'strain_group' (WT/PKO) anhand der 'biosensor'-Spalte.

    PKO liegt als eigener Ordner auf der Biosensor-Ebene
    (Data/PKO/<osc_type>/<periode>/) und landet damit in 'biosensor' - eine
    echte 'strain'-Spalte liefert Combined_Results nicht (siehe
    config.PANEL_A_GROUP_COL). Die Zuordnung passiert deshalb hier und NUR
    hier, statt global eine Spalte anzulegen, die die bestehenden Auswertungen
    veraendern wuerde.
    """
    if "biosensor" not in df.columns:
        raise ValueError("add_strain_group() braucht eine 'biosensor'-Spalte.")
    out = df.copy()
    out["strain_group"] = np.where(
        out["biosensor"].astype(str) == str(pko_biosensor),
        STRAIN_GROUP_PKO, STRAIN_GROUP_WT,
    )
    return out


def add_condition_type(df: pd.DataFrame) -> pd.DataFrame:
    """Ergaenzt 'condition_type' (PosCtrl/NegCtrl/Oscillation/...), falls noetig."""
    if "condition_type" in df.columns:
        return df
    if "condition" not in df.columns:
        raise ValueError("add_condition_type() braucht eine 'condition'-Spalte.")
    out = df.copy()
    out["condition_type"] = out["condition"].map(classify_condition_type)
    return out


def matched_periods(cells_pko: pd.DataFrame) -> list[str]:
    """Die Perioden, die PKO abdeckt - darauf wird der WT-Teil eingeschraenkt.

    NICHT hart codiert: kommt spaeter ein PKO-Datensatz mit mehr Perioden
    dazu, verglich eine feste Liste stillschweigend gegen die falsche
    WT-Teilmenge. Der aktive Wert wird geloggt, damit im Laborbuch steht,
    welche Perioden in den Vergleich eingegangen sind.
    """
    if "osc_freq" not in cells_pko.columns:
        raise ValueError("matched_periods() braucht eine 'osc_freq'-Spalte.")
    periods = sorted(cells_pko["osc_freq"].dropna().astype(str).unique().tolist())
    if not periods:
        logger.warning("matched_periods(): PKO-Daten enthalten keine 'osc_freq'-Werte.")
    elif len(periods) == 1:
        logger.info(
            "PKO deckt genau eine Periode ab (%s). Der WT-Teil wird darauf eingeschraenkt; "
            "ein Test ueber osc_freq-Batches ist damit nicht moeglich (siehe Modul-Docstring).",
            periods[0],
        )
    else:
        logger.info("PKO deckt die Perioden %s ab; der WT-Teil wird darauf eingeschraenkt.", periods)
    return periods


def build_joint_controls(
    cells_wt: pd.DataFrame,
    cells_pko: pd.DataFrame,
    pko_biosensor: str,
) -> tuple[pd.DataFrame, list[str]]:
    """Gemeinsamer Kontroll-Datensatz beider Staemme, auf PKOs Perioden beschnitten.

    Gibt nur PosCtrl-/NegCtrl-Zeilen zurueck: die Oszillationsbedingung selbst
    ist fuer die Clogging-Frage der schlechtere Ort (dort unterscheiden sich
    WT und PKO ohnehin physiologisch), und die Kontrollen sind zwischen den
    Staemmen exakt vergleichbar, weil ihr Medienverlauf konstant ist.
    """
    periods = matched_periods(cells_pko)

    wt = add_condition_type(add_strain_group(cells_wt, pko_biosensor))
    pko = add_condition_type(add_strain_group(cells_pko, pko_biosensor))

    if periods:
        wt = wt[wt["osc_freq"].astype(str).isin(periods)]

    joint = pd.concat([wt, pko], ignore_index=True)
    joint = joint[joint["condition_type"].isin(CONTROL_ORDER)].copy()

    if joint.empty:
        logger.warning(
            "build_joint_controls(): keine PosCtrl-/NegCtrl-Zeilen im gemeinsamen Datensatz. "
            "Nutzen WT und PKO dieselben condition-Suffixe '_PosCtrl'/'_NegCtrl'? "
            "classify_condition_type() vergleicht exakt auf diese Endungen."
        )
        return joint, periods

    counts = (
        joint.groupby(["strain_group", "condition_type"])["exp_id"].nunique()
        if "exp_id" in joint.columns else pd.Series(dtype=int)
    )
    logger.info("Gemeinsame Kontrollkammern (exp_id) je Stamm/Kontrollart:\n%s", counts.to_string())
    return joint, periods


# ==============================================================================
# 1. Kammer-Uebereinstimmung
# ==============================================================================

def summarise_control_chambers(
    control_cells: pd.DataFrame,
    value_col: str = "area",
    analysis_start_min: float = 120.0,
) -> pd.DataFrame:
    """Ein Wert pro Kontrollkammer, zweistufig aggregiert.

    Stufe 1: Median ueber die Zellen EINES Frames - robust gegen einzelne
    Segmentierungs-Ausreisser. Stufe 2: Median ueber die Frames nach der
    Vorkonditionierungsphase. Damit zaehlt jede Kammer genau einmal,
    unabhaengig davon, wie viele Zellen oder Frames sie geliefert hat.

    Gleiche Konstruktion wie
    queen_controls.summarise_sensor_controls_by_chamber(), nur auf einer
    Morphologie-Spalte statt auf einem Sensor-Ratio - PKO hat keine
    Fluoreszenzkanaele.
    """
    if control_cells.empty:
        return pd.DataFrame()

    missing = [c for c in ("chamber", "time_h", value_col) if c not in control_cells.columns]
    if missing:
        logger.warning("summarise_control_chambers(): Spalten %s fehlen - uebersprungen.", missing)
        return pd.DataFrame()

    group_cols = [
        c for c in ["strain_group", "biosensor", "osc_type", "osc_freq", "condition",
                    "condition_type", "replicate", "chamber", "exp_id", "time_h"]
        if c in control_cells.columns
    ]
    per_time = (
        control_cells.dropna(subset=[value_col])
        .groupby(group_cols, dropna=False)[value_col]
        .agg(median_value="median", n_cells="count")
        .reset_index()
    )
    if per_time.empty:
        return pd.DataFrame()

    start_h = _minutes_to_hours(analysis_start_min)
    late = per_time[per_time["time_h"] >= start_h]
    if late.empty:
        logger.warning(
            "summarise_control_chambers(): keine Zeitpunkte ab %.0f min (= %.2f h); "
            "die Aufnahme deckt %.2f-%.2f h ab. OSCILLATION_START_MIN gegen die echte "
            "Laufzeit pruefen.",
            analysis_start_min, start_h, per_time["time_h"].min(), per_time["time_h"].max(),
        )
        return pd.DataFrame()

    summary_groups = [c for c in group_cols if c != "time_h"]
    return (
        late.groupby(summary_groups, dropna=False)["median_value"]
        .agg(control_value="median", n_timepoints="count")
        .reset_index()
    )


def compute_chamber_agreement(
    chamber_summary: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Streuung nominell identischer Kontrollkammern INNERHALB eines Replikats.

    Warum innerhalb und nicht ueber alle Kammern: Kammern eines Replikats sind
    technische Messungen desselben Chips (so sagt es auch
    queen_controls.summarise_sensor_controls_by_chamber()). Ihre Streuung ist
    damit der technische/hydraulische Anteil - genau der, den eine Verstopfung
    aufblaehen sollte. Die Streuung ZWISCHEN Replikaten ist biologisch; sie
    mit hineinzurechnen wuerde den gesuchten Effekt verwaessern.

    Als Streuungsmass der Variationskoeffizient (SD/Mittelwert) ueber die
    Kammern eines Replikats: dimensionslos, damit zwischen Staemmen mit
    unterschiedlichem Niveau vergleichbar - eine PKO-Zelle muss nicht
    dieselbe absolute Flaeche haben wie eine WT-Zelle.

    Rueckgabe: (pro Replikat, pro Stamm). Die Replikat-Ebene ist die
    Datengrundlage der Abbildung, die Stamm-Ebene die der Aussage.
    """
    if chamber_summary.empty:
        return pd.DataFrame(), pd.DataFrame()

    rep_groups = [
        c for c in ["strain_group", "biosensor", "osc_type", "osc_freq", "condition_type", "replicate"]
        if c in chamber_summary.columns
    ]
    per_replicate = (
        chamber_summary.groupby(rep_groups, dropna=False)["control_value"]
        .agg(
            mean_value="mean",
            sd_value=lambda s: s.std(ddof=1),
            n_chambers="count",
        )
        .reset_index()
    )

    # Ein Replikat mit nur EINER Kammer hat keine Kammer-zu-Kammer-Streuung -
    # das ist keine Null, sondern eine fehlende Messung. Explizit ausschliessen,
    # statt eine 0 in den Vergleich zu lassen, die PKO kuenstlich gut aussehen
    # liesse.
    single = per_replicate[per_replicate["n_chambers"] < 2]
    if not single.empty:
        logger.warning(
            "compute_chamber_agreement(): %d Replikat(e) haben nur eine Kontrollkammer und "
            "fallen aus dem Uebereinstimmungs-Vergleich (keine Kammer-zu-Kammer-Streuung "
            "bestimmbar).",
            len(single),
        )
    per_replicate = per_replicate[per_replicate["n_chambers"] >= 2].copy()
    if per_replicate.empty:
        logger.warning(
            "compute_chamber_agreement(): kein Replikat mit mindestens zwei Kontrollkammern - "
            "die Kammer-Uebereinstimmung ist mit diesen Daten nicht bestimmbar."
        )
        return per_replicate, pd.DataFrame()

    per_replicate["chamber_cv"] = per_replicate["sd_value"] / per_replicate["mean_value"].abs()

    strain_groups = [c for c in ["strain_group", "biosensor", "condition_type"]
                     if c in per_replicate.columns]
    per_strain = (
        per_replicate.groupby(strain_groups, dropna=False)["chamber_cv"]
        .agg(median_chamber_cv="median", min_chamber_cv="min", max_chamber_cv="max",
             n_replicates="count")
        .reset_index()
    )
    return per_replicate, per_strain


# ==============================================================================
# 2. Kontroll-Bracket
# ==============================================================================

def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """Cliff's Delta: P(a > b) - P(a < b), in [-1, 1].

    Nichtparametrisch und ohne Verteilungsannahme; 0 = keine Trennung,
    1 = a liegt vollstaendig ueber b. Aus der Mann-Whitney-U-Statistik
    berechnet (delta = 2U/(n_a*n_b) - 1), das spart die O(n^2)-Paarbildung.

    Bewusst eine EFFEKTSTAERKE und kein p-Wert: der p-Wert eines
    Mann-Whitney ueber Zellen haengt fast nur an der Zellzahl (siehe
    config.METHOD_CAVEATS zur Pseudoreplikation), die Effektstaerke nicht.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size == 0 or b.size == 0:
        return float("nan")
    u = stats.mannwhitneyu(a, b, alternative="two-sided").statistic
    return float(2.0 * u / (a.size * b.size) - 1.0)


def compute_control_bracket(
    area_table: pd.DataFrame,
    value_col: str = "mu_area",
    exclude_unreliable: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Trennt das PosCtrl/NegCtrl-Bracket ueberhaupt? Cliff's Delta auf mu_area.

    Pro Replikat gerechnet - die Effektstaerke ist eine Eigenschaft der
    Zellverteilungen DIESES Replikats -, dann ueber Replikate zusammengefasst.
    Also zweistufig und nicht ueber Zellen gepoolt.

    cell_type wird ABSICHTLICH nicht gefiltert: mu_area ist eine Regression
    ueber ln(Flaeche) eines Tracks und damit von der Mutter/Bud-Heuristik
    unabhaengig - nur das cell_type-Label haengt daran. Ein Filter auf
    "mother" wuerde die unvalidierte Heuristik ohne Gegenwert einschleppen.

    WARUM exclude_unreliable HIER PER DEFAULT *AUS* IST
    ---------------------------------------------------
    fit_is_reliable heisst R^2 >= min_r_squared (Default 0.5, siehe
    area_growth.compute_area_growth_rate()). Bei einer Bedingung, die
    tatsaechlich FLACH ist, erklaert die Regressionsgerade per Konstruktion
    kaum Varianz - das R^2 ist niedrig, WEIL es nichts zu erklaeren gibt, und
    nicht, weil der Fit schlecht waere. Der Filter wirft dann genau die
    ehrlichen flachen Fits weg und behaelt die, in denen das Rauschen zufaellig
    wie ein Trend aussieht. Der ueberlebende Median ist dadurch nach OBEN
    verzerrt (Survivorship Bias).

    Genau das trifft NegCtrl: durchgehende Starvation SOLL mu_area ~ 0 liefern.
    Im Test auf synthetischen Daten blieben von 119 NegCtrl-Fits 11 uebrig,
    und deren Median (0.080 h^-1) lag HOEHER als der des echt wachsenden
    PosCtrl (0.062 h^-1) - das Bracket kippte dadurch im Vorzeichen.

    Fuer diesen Vergleich gehen deshalb ALLE Fits ein: die Steigung eines
    flachen Tracks ist ~0, und das ist die richtige Antwort. Statt zu filtern
    wird der Anteil zuverlaessiger Fits pro Arm als Diagnose mitgeschrieben
    (frac_reliable_*) - ein niedriger Anteil in NegCtrl ist selbst ein Befund
    ("die Kontrolle ist flach"), kein Grund, Daten zu verwerfen.

    ACHTUNG fuer die uebrige Pipeline: dieselbe Verzerrung steckt in jeder
    Auswertung, die auf einer flachen Bedingung nach fit_is_reliable filtert -
    u.a. summarise_area_growth(exclude_unreliable=True) und die R(p)-Rechnung
    fuer mu_area in Schritt 40. Dort ist sie NICHT korrigiert; das ist eine
    fachliche Entscheidung und gehoert diskutiert, nicht stillschweigend
    geaendert.

    Erwartung: PosCtrl (durchgehend Feast) waechst, NegCtrl (durchgehend
    Starvation) nicht - delta deutlich > 0. Ein delta um 0 heisst, dass das
    Bracket nicht funktioniert, und genau das ist die Beobachtung, die
    erklaert werden soll.
    """
    if area_table.empty:
        return pd.DataFrame(), pd.DataFrame()

    df = add_condition_type(area_table)
    if exclude_unreliable and "fit_is_reliable" in df.columns:
        logger.warning(
            "compute_control_bracket(exclude_unreliable=True): filtert auf fit_is_reliable. "
            "Bei einer flachen Bedingung (NegCtrl) verzerrt das den Median nach oben - "
            "siehe Docstring."
        )
        df = df[df["fit_is_reliable"]]
    df = df.dropna(subset=[value_col])
    df = df[df["condition_type"].isin(CONTROL_ORDER)]
    if df.empty:
        logger.warning(
            "compute_control_bracket(): keine zuverlaessigen %s-Fits in Kontrollkammern.", value_col
        )
        return pd.DataFrame(), pd.DataFrame()

    rep_groups = [c for c in ["strain_group", "biosensor", "osc_type", "osc_freq", "replicate"]
                  if c in df.columns]
    records = []
    for keys, grp in df.groupby(rep_groups, dropna=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        pos_rows = grp[grp["condition_type"] == "PosCtrl"]
        neg_rows = grp[grp["condition_type"] == "NegCtrl"]
        pos = pos_rows[value_col].to_numpy()
        neg = neg_rows[value_col].to_numpy()
        record = dict(zip(rep_groups, keys))
        record["n_posctrl_cells"] = int(pos.size)
        record["n_negctrl_cells"] = int(neg.size)
        record["median_posctrl"] = float(np.median(pos)) if pos.size else float("nan")
        record["median_negctrl"] = float(np.median(neg)) if neg.size else float("nan")
        record["bracket_delta"] = cliffs_delta(pos, neg)
        # Diagnose statt Filter (siehe Docstring): ein niedriger Anteil
        # zuverlaessiger Fits in NegCtrl heisst "die Kontrolle ist flach" und
        # ist damit selbst Teil des Befunds.
        for label, rows in (("posctrl", pos_rows), ("negctrl", neg_rows)):
            if "fit_is_reliable" in rows.columns and len(rows):
                record[f"frac_reliable_{label}"] = float(rows["fit_is_reliable"].mean())
            else:
                record[f"frac_reliable_{label}"] = float("nan")
        records.append(record)

    per_replicate = pd.DataFrame(records)
    if per_replicate.empty:
        return per_replicate, pd.DataFrame()

    n_missing = int(per_replicate["bracket_delta"].isna().sum())
    if n_missing:
        logger.warning(
            "compute_control_bracket(): %d Replikat(e) ohne beide Kontrollarten - "
            "bracket_delta bleibt dort NaN (nicht 0).", n_missing,
        )

    strain_groups = [c for c in ["strain_group", "biosensor"] if c in per_replicate.columns]
    per_strain = (
        per_replicate.dropna(subset=["bracket_delta"])
        .groupby(strain_groups, dropna=False)["bracket_delta"]
        .agg(median_bracket_delta="median", min_bracket_delta="min",
             max_bracket_delta="max", n_replicates="count")
        .reset_index()
    )
    return per_replicate, per_strain


# ==============================================================================
# 3. Zell-Ausbeute
# ==============================================================================

def compute_cell_yield(
    control_cells: pd.DataFrame,
    cell_id_col: str = "cell_uid",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Verfolgbare Zellen pro Kontrollkammer ueber die Zeit, plus Steigung.

    Auf den ersten Zeitpunkt normiert (n(t)/n(t_0)): ohne das dominiert die
    Aussaatdichte den Vergleich zwischen Kammern, und ein Unterschied zwischen
    WT und PKO waere nicht von "unterschiedlich dicht angesetzt" zu trennen.

    Die Steigung ist eine Kleinste-Quadrate-Gerade durch die normierten Werte
    in 1/h. Negativ heisst: die Kammer verliert ueber die Zeit verfolgbare
    Zellen. Fuer PosCtrl (durchgehend Feast) ist eine POSITIVE Steigung das
    Erwartete - die Zellen wachsen -, ein Wegbrechen dort ist das eigentliche
    Signal.

    ACHTUNG bei der Interpretation: die Zellzahl haengt auch an der
    Segmentierungs- und Trackingqualitaet, nicht nur an der Stroemung. Diese
    Groesse STUETZT die Hypothese, sie beweist sie nicht - so gehoert sie auch
    in die Legende.
    """
    if control_cells.empty:
        return pd.DataFrame(), pd.DataFrame()

    missing = [c for c in ("time_h", cell_id_col) if c not in control_cells.columns]
    if missing:
        logger.warning("compute_cell_yield(): Spalten %s fehlen - uebersprungen.", missing)
        return pd.DataFrame(), pd.DataFrame()

    chamber_cols = [
        c for c in ["strain_group", "biosensor", "osc_type", "osc_freq", "condition",
                    "condition_type", "replicate", "chamber", "exp_id"]
        if c in control_cells.columns
    ]
    timeseries = (
        control_cells.groupby(chamber_cols + ["time_h"], dropna=False)[cell_id_col]
        .nunique()
        .reset_index(name="n_cells")
        .sort_values(chamber_cols + ["time_h"])
    )
    if timeseries.empty:
        return timeseries, pd.DataFrame()

    # Normierung auf den ERSTEN Zeitpunkt der jeweiligen Kammer (nicht auf das
    # globale Minimum: Kammern koennen unterschiedlich spaet starten).
    first = timeseries.groupby(chamber_cols, dropna=False)["n_cells"].transform("first")
    timeseries["n_cells_relative"] = timeseries["n_cells"] / first.replace(0, np.nan)

    records = []
    for keys, grp in timeseries.groupby(chamber_cols, dropna=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        sub = grp.dropna(subset=["n_cells_relative"])
        record = dict(zip(chamber_cols, keys))
        record["n_timepoints"] = int(len(sub))
        record["n_cells_start"] = int(grp["n_cells"].iloc[0])
        record["n_cells_end"] = int(grp["n_cells"].iloc[-1])
        if len(sub) >= 3:
            slope, _ = np.polyfit(sub["time_h"].to_numpy(float),
                                  sub["n_cells_relative"].to_numpy(float), 1)
            record["relative_slope_per_h"] = float(slope)
        else:
            # Unter drei Zeitpunkten ist eine Steigung Rauschen, kein Trend.
            record["relative_slope_per_h"] = float("nan")
        records.append(record)

    return timeseries, pd.DataFrame(records)


# ==============================================================================
# Abbildungen
# ==============================================================================

def _strain_order(df: pd.DataFrame) -> list[str]:
    """WT-Staemme alphabetisch, PKO ans Ende - PKO ist die Vergleichsgruppe."""
    if "biosensor" not in df.columns:
        return []
    wt = sorted(df.loc[df["strain_group"] == STRAIN_GROUP_WT, "biosensor"].dropna().astype(str).unique())
    pko = sorted(df.loc[df["strain_group"] == STRAIN_GROUP_PKO, "biosensor"].dropna().astype(str).unique())
    return wt + pko


def _scatter_by_strain(
    ax,
    df: pd.DataFrame,
    value_col: str,
    strain_order: list[str],
    seed: int = 0,
) -> None:
    """Punkte je Stamm, aufgeteilt nach Kontrollart, mit Median-Balken.

    Ohne 'condition_type'-Spalte (z.B. beim Bracket-Delta, das ja selbst schon
    der Vergleich BEIDER Kontrollarten ist) wird nach strain_group gefaerbt,
    NICHT in einer der Kontrollfarben: sonst bedeutete dasselbe Blau in einem
    Panel "NegCtrl" und im anderen "alle Kontrollen zusammen", was die Legende
    des ersten Panels stillschweigend falsch macht.
    """
    rng = np.random.default_rng(seed)
    has_control = "condition_type" in df.columns
    controls = CONTROL_ORDER if has_control else [None]
    for x, strain in enumerate(strain_order):
        sub_strain = df[df["biosensor"].astype(str) == strain]
        for offset, control in zip(np.linspace(-0.16, 0.16, len(controls)), controls):
            sub = sub_strain if control is None else sub_strain[sub_strain["condition_type"] == control]
            values = sub[value_col].dropna()
            if values.empty:
                continue
            if control is not None:
                color = CONTROL_COLORS[control]
            else:
                groups = sub["strain_group"].dropna().unique() if "strain_group" in sub.columns else []
                color = STRAIN_GROUP_COLORS.get(groups[0], "#6E6E6E") if len(groups) else "#6E6E6E"
            jitter = rng.uniform(-0.045, 0.045, size=len(values))
            ax.scatter(
                np.full(len(values), x + offset) + jitter, values,
                color=color, s=34, alpha=0.85,
                edgecolor="white", linewidth=0.5, zorder=3,
                label=control if x == 0 and control is not None else None,
            )
            ax.plot([x + offset - 0.07, x + offset + 0.07], [values.median()] * 2,
                    color="black", linewidth=1.6, zorder=4)
    ax.set_xticks(range(len(strain_order)))
    ax.set_xticklabels(strain_order, rotation=25, ha="right")


def plot_pko_control_agreement(
    agreement_per_replicate: pd.DataFrame,
    bracket_per_replicate: pd.DataFrame,
    out_path: Path,
    value_label: str = "cell area",
) -> None:
    """Zwei Panels: Kammer-Uebereinstimmung und Kontroll-Bracket, WT gegen PKO.

    Untereinander in EINER Abbildung, weil sie zusammen gelesen werden muessen:
    eine hydraulische Ursache sagt beides voraus (hohe Streuung UND kaputtes
    Bracket), eine metabolische keines von beiden. Die vier WT-Staemme stehen
    einzeln auf der x-Achse - sichtbar zu machen, ob sie untereinander
    uebereinstimmen, ist die Voraussetzung dafuer, sie als "WT" zusammenfassen
    zu duerfen.
    """
    panels = [
        (agreement_per_replicate, "chamber_cv",
         f"Chamber-to-chamber CV\n({value_label}, within one replicate)",
         "a) Do nominally identical control chambers agree?"),
        (bracket_per_replicate, "bracket_delta",
         "PosCtrl vs NegCtrl\nCliff's δ on µ_area",
         "b) Does the control bracket separate at all?"),
    ]
    available = [(df, col, ylab, title) for df, col, ylab, title in panels
                 if df is not None and not df.empty and col in df.columns]
    if not available:
        logger.warning("plot_pko_control_agreement(): keine Daten fuer beide Panels - uebersprungen.")
        return

    strain_order = _strain_order(pd.concat([df for df, *_ in available], ignore_index=True))
    if not strain_order:
        logger.warning("plot_pko_control_agreement(): keine Staemme bestimmbar - uebersprungen.")
        return

    fig, axes = plt.subplots(
        len(available), 1, figsize=(1.5 * len(strain_order) + 3.2, 3.6 * len(available)),
        squeeze=False, sharex=True,
    )
    axes = axes[:, 0]

    for ax, (df, col, ylabel, title) in zip(axes, available):
        _scatter_by_strain(ax, df, col, strain_order)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.set_title(title, fontsize=10, loc="left")
        ax.grid(axis="y", alpha=0.25, linewidth=0.6)
        if col == "bracket_delta":
            # delta = 0 ist die inhaltliche Nulllinie: keine Trennung der Kontrollen.
            ax.axhline(0.0, color="black", linewidth=0.9, linestyle="--", alpha=0.7)
        if ax is axes[0]:
            handles, labels = ax.get_legend_handles_labels()
            if handles:
                ax.legend(handles, labels, frameon=False, fontsize=8, loc="best")

    axes[-1].set_xlabel("Strain (WT biosensor strains, then PKO)")
    fig.suptitle("WT vs PKO: control-chamber behaviour at matched oscillation period", y=0.995)
    fig.text(
        0.5, -0.015,
        "Each point = one replicate. Black bar = median. PKO is a single strain, so no "
        "WT-vs-PKO significance test is computed;\nthe four WT strains are the internal "
        "replication of the WT side — read whether PKO falls inside or outside their spread.",
        ha="center", fontsize=8,
    )
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=180)
    plt.close(fig)
    logger.info("Plot gespeichert: %s", out_path.name)


def plot_pko_cell_yield(
    yield_timeseries: pd.DataFrame,
    yield_slopes: pd.DataFrame,
    out_path: Path,
) -> None:
    """Normierte Zell-Ausbeute pro Kontrollkammer, WT gegen PKO, plus Steigungen.

    Eine Verstopfung sollte sich als fortschreitender Verlust zeigen - und
    zwar besonders in PosCtrl, wo 10 h durchgewachsen wird und damit am
    meisten Pullulan anfaellt.
    """
    if yield_timeseries is None or yield_timeseries.empty:
        logger.warning("plot_pko_cell_yield(): keine Zeitreihe - uebersprungen.")
        return

    strain_groups = [g for g in (STRAIN_GROUP_WT, STRAIN_GROUP_PKO)
                     if g in set(yield_timeseries["strain_group"])]
    controls = [c for c in CONTROL_ORDER
                if c in set(yield_timeseries.get("condition_type", pd.Series(dtype=str)))]
    if not strain_groups or not controls:
        logger.warning("plot_pko_cell_yield(): Stamm- oder Kontrollgruppen fehlen - uebersprungen.")
        return

    n_cols = len(strain_groups) + (1 if yield_slopes is not None and not yield_slopes.empty else 0)
    fig, axes = plt.subplots(len(controls), n_cols, figsize=(3.6 * n_cols, 3.3 * len(controls)),
                             squeeze=False, sharey="row")

    chamber_cols = [c for c in ["exp_id", "replicate", "chamber"] if c in yield_timeseries.columns]
    for row, control in enumerate(controls):
        for col, group in enumerate(strain_groups):
            ax = axes[row][col]
            sub = yield_timeseries[
                (yield_timeseries["strain_group"] == group)
                & (yield_timeseries["condition_type"] == control)
            ]
            for _, chamber in sub.groupby(chamber_cols, dropna=False) if chamber_cols else []:
                ax.plot(chamber["time_h"], chamber["n_cells_relative"],
                        color=STRAIN_GROUP_COLORS.get(group, "#4C78A8"),
                        alpha=0.55, linewidth=1.1)
            ax.axhline(1.0, color="black", linewidth=0.8, linestyle=":", alpha=0.6)
            ax.set_title(f"{group} — {control}", fontsize=10)
            if col == 0:
                ax.set_ylabel("Trackable cells\n(relative to first frame)", fontsize=9)
            ax.set_xlabel("Time [h]", fontsize=9)
            ax.grid(alpha=0.22, linewidth=0.6)

        if n_cols > len(strain_groups):
            ax = axes[row][-1]
            sub = yield_slopes[yield_slopes["condition_type"] == control] \
                if "condition_type" in yield_slopes.columns else yield_slopes
            strain_order = _strain_order(sub)
            if strain_order:
                _scatter_by_strain(ax, sub.assign(condition_type=control),
                                   "relative_slope_per_h", strain_order, seed=row)
            ax.axhline(0.0, color="black", linewidth=0.9, linestyle="--", alpha=0.7)
            ax.set_ylabel("Relative yield slope [1/h]", fontsize=9)
            ax.set_title(f"Slope per chamber — {control}", fontsize=10)
            ax.grid(axis="y", alpha=0.22, linewidth=0.6)

    fig.suptitle("WT vs PKO: trackable-cell yield in control chambers", y=0.995)
    fig.text(
        0.5, -0.015,
        "One line per control chamber, normalised to its own first frame. Cell counts also "
        "reflect segmentation and seeding density,\nso this supports the occlusion "
        "explanation rather than demonstrating it.",
        ha="center", fontsize=8,
    )
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=180)
    plt.close(fig)
    logger.info("Plot gespeichert: %s", out_path.name)


# ==============================================================================
# Orchestrierung
# ==============================================================================

def run_pko_comparison(
    cells_wt: pd.DataFrame,
    cells_pko: pd.DataFrame,
    area_wt: pd.DataFrame,
    area_pko: pd.DataFrame,
    output_dir: Path,
    pko_biosensor: str = "PKO",
    value_col: str = "area",
    analysis_start_min: float = 120.0,
) -> None:
    """Der komplette WT-gegen-PKO-Vergleich: drei Groessen, zwei Abbildungen.

    Schreibt NUR in output_dir (per Default der PKO-Ordner). Die bestehenden
    Oszillations- und statischen Ergebnisse werden nicht angefasst - das ist
    der ganze Punkt eines dritten, unabhaengigen Zweigs.

    area_wt/area_pko sind die bereits berechneten area_table-Objekte der
    beiden Kontexte; sie werden hier nur zusammengefuehrt und auf die
    gemeinsamen Perioden beschnitten, statt die Fits neu zu rechnen.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    if cells_pko.empty:
        logger.warning("run_pko_comparison(): keine PKO-Daten - Vergleich uebersprungen.")
        return
    if cells_wt.empty:
        logger.warning("run_pko_comparison(): keine WT-Daten - Vergleich uebersprungen.")
        return

    joint_controls, periods = build_joint_controls(cells_wt, cells_pko, pko_biosensor)
    if joint_controls.empty:
        return

    # --- 1. Kammer-Uebereinstimmung ---------------------------------------
    chamber_summary = summarise_control_chambers(
        joint_controls, value_col=value_col, analysis_start_min=analysis_start_min,
    )
    if not chamber_summary.empty:
        chamber_summary.to_csv(output_dir / "60_pko_control_chambers.csv", index=False)
        logger.info("Tabelle gespeichert: 60_pko_control_chambers.csv")

    agreement_rep, agreement_strain = compute_chamber_agreement(chamber_summary)
    if not agreement_rep.empty:
        agreement_rep.to_csv(output_dir / "60_pko_chamber_agreement_per_replicate.csv", index=False)
    if not agreement_strain.empty:
        agreement_strain.to_csv(output_dir / "60_pko_chamber_agreement_per_strain.csv", index=False)
        logger.info(
            "Kammer-Uebereinstimmung (Median CV je Stamm/Kontrollart):\n%s",
            agreement_strain.to_string(index=False),
        )

    # --- 2. Kontroll-Bracket ----------------------------------------------
    area_frames = []
    for table, is_pko in ((area_wt, False), (area_pko, True)):
        if table is None or table.empty:
            continue
        tagged = add_strain_group(table, pko_biosensor)
        if not is_pko and periods and "osc_freq" in tagged.columns:
            tagged = tagged[tagged["osc_freq"].astype(str).isin(periods)]
        area_frames.append(tagged)

    bracket_rep, bracket_strain = (pd.DataFrame(), pd.DataFrame())
    if area_frames:
        bracket_rep, bracket_strain = compute_control_bracket(
            pd.concat(area_frames, ignore_index=True)
        )
        if not bracket_rep.empty:
            bracket_rep.to_csv(output_dir / "60_pko_control_bracket_per_replicate.csv", index=False)
        if not bracket_strain.empty:
            bracket_strain.to_csv(output_dir / "60_pko_control_bracket_per_strain.csv", index=False)
            logger.info(
                "Kontroll-Bracket (Median Cliff's delta je Stamm):\n%s",
                bracket_strain.to_string(index=False),
            )
    else:
        logger.warning("run_pko_comparison(): keine area_table-Daten - Bracket-Panel entfaellt.")

    # --- 3. Zell-Ausbeute -------------------------------------------------
    yield_ts, yield_slopes = compute_cell_yield(joint_controls)
    if not yield_ts.empty:
        yield_ts.to_csv(output_dir / "60_pko_cell_yield_timeseries.csv", index=False)
    if not yield_slopes.empty:
        yield_slopes.to_csv(output_dir / "60_pko_cell_yield_slopes.csv", index=False)
        logger.info("Tabellen gespeichert: 60_pko_cell_yield_timeseries.csv / _slopes.csv")

    # --- Abbildungen ------------------------------------------------------
    plot_pko_control_agreement(
        agreement_rep, bracket_rep,
        out_path=output_dir / "61_pko_control_agreement.pdf",
        value_label=value_col,
    )
    plot_pko_cell_yield(
        yield_ts, yield_slopes, out_path=output_dir / "62_pko_cell_yield.pdf",
    )

    logger.warning(
        "PKO-VERGLEICH: PKO ist EIN Stamm - es wird bewusst KEIN Signifikanztest "
        "WT-gegen-PKO gerechnet (n=1 auf Stammebene). Die Aussage ist deskriptiv: liegt der "
        "PKO-Wert innerhalb oder ausserhalb der Spanne der WT-Staemme? Ausserdem ist PKO eine "
        "Mutante mit veraendertem Stoffwechsel - unterscheidbar sind hydraulische und "
        "metabolische Erklaerung nur ueber das MUSTER (Streuung/Steigung), nicht ueber einen "
        "Niveau-Unterschied allein. Siehe pko_comparison.py Docstring."
    )
