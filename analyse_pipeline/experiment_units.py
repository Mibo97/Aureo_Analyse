"""
experiment_units.py
===================
Was ist in diesem Datensatz eigentlich ein Replikat? Chip, Kammer, Medium.

DIE VERSUCHSSTRUKTUR, WIE SIE WIRKLICH IST
------------------------------------------
Die Spalten 'replicate' und 'chamber' kommen aus dem Dateinamen der
Bildverarbeitung und bedeuten NICHT, was ihre Namen nahelegen:

  Oszillation / PKO
      Pro (Stamm, osc_type, Periode) gibt es EINEN Chip, an EINEM Tag, aus
      EINER Vorkultur. Auf dem Chip liegen ~11 Kammern (5 Oszillation, 3
      PosCtrl, 3 NegCtrl), verteilt auf mehrere Arrays mit wiederholten
      Positionslabels - deshalb steht 'ChamA13' dreimal in einem Batch.
      'replicate' (Rep1, Rep2, ...) ist ein Array-/Laufindex, der gleiche
      Positionen auseinanderhaelt. Es ist KEIN biologisches Replikat.

      => Biologische Einheit = der Chip = der Batch. n = 1 pro Bedingung.
         Alle Kammern eines Chips sind TECHNISCHE Replikate einer Kultur.
         Die Dosis-Wirkung ueber die Perioden ist ein Chip pro Dosis.

  Statisch
      Jedes 'replicate' ist ein eigener Chip mit eigener Vorkultur; 'chamber'
      ist durchgehend 'ChamA0' und traegt keine Information. Das Medium steht
      in 'condition' (St.omlp / St.ypd); 'osc_freq' benennt die Chip-Familie
      (W109: static_omlp/static_ypd, W65: static_W65 mit beiden Medien).

      => Biologische Einheit = das Replikat = der Chip. n = 4-5 pro Medium.

Dieses Modul leitet daraus die Spalten ab, auf die der Rest der Pipeline
aggregiert, und stellt EINE hierarchische Aggregation bereit
(Zelle -> Kammer -> Chip -> Bedingung), die den Fehlerbalken ehrlich
beschriftet: ueber Chips, wenn es mehr als einen gibt, sonst ueber Kammern -
und dann steht 'chamber' in der Spalte 'error_unit', nicht 'replicate'.

Vorher gruppierten alle 'per replicate'-Aggregationen auf den Array-Index und
nannten das Ergebnis biologisch. Das war falsch, und zwar in JEDER Abbildung
der Oszillationsdaten.
"""

from __future__ import annotations

import logging
import re
from typing import Mapping, Optional, Sequence

import numpy as np
import pandas as pd

from growth_rate import classify_condition_type

logger = logging.getLogger(__name__)

DATE_RE = re.compile(r"^(\d{6,8})[_-]")

# Dokumentiert in jeder Ausgabe, welche Regel die Einheiten erzeugt hat.
UNIT_RULE_OSC = "one chip per (biosensor, osc_type, osc_freq); 'replicate' = array index"
UNIT_RULE_STATIC = "one chip per 'replicate'; 'osc_freq' = chip family; medium from 'condition'"


def is_static_mask(df: pd.DataFrame) -> pd.Series:
    return df["osc_type"].astype(str).str.strip().str.lower() == "static"


def add_experiment_units(
    cells: pd.DataFrame,
    static_chip_labels: Mapping[str, str],
    static_medium_prefix: str = "St.",
    filename_col: str = "filename",
) -> pd.DataFrame:
    """Ergaenzt date, medium, chip_family, chip, unit_rule.

    date         : fuehrende Ziffern des Dateinamens (Laborbuch-Referenz), sonst ''
    medium       : statisch: condition ohne Prefix ('St.omlp' -> 'omlp'); sonst NaN
    chip_family  : statisch: static_chip_labels[osc_freq] (z.B. 'W65'); sonst NaN
    chip         : die biologische Einheit - siehe Modul-Docstring
    unit_rule    : welche Regel 'chip' erzeugt hat (fuer die CSVs)
    """
    required = {"biosensor", "osc_type", "osc_freq", "condition", "replicate"}
    missing = required.difference(cells.columns)
    if missing:
        raise ValueError(f"add_experiment_units(): Spalten fehlen: {sorted(missing)}")

    out = cells.copy()
    if filename_col in out.columns:
        out["date"] = out[filename_col].astype(str).map(
            lambda f: DATE_RE.match(f).group(1) if DATE_RE.match(f) else ""
        )
    else:
        out["date"] = ""
        logger.warning("add_experiment_units(): keine '%s'-Spalte - 'date' bleibt leer.", filename_col)

    static = is_static_mask(out)
    # object-dtype von Anfang an: pandas >= 3 verweigert, Strings in eine
    # float-NaN-Spalte zu schreiben.
    out["medium"] = pd.Series(None, index=out.index, dtype=object)
    out["chip_family"] = pd.Series(None, index=out.index, dtype=object)
    out["chip"] = ""
    out["unit_rule"] = ""

    if static.any():
        cond = out.loc[static, "condition"].astype(str)
        out.loc[static, "medium"] = cond.str.replace(f"^{re.escape(static_medium_prefix)}", "", regex=True)
        freq = out.loc[static, "osc_freq"].astype(str)
        unknown = sorted(set(freq.unique()) - set(static_chip_labels))
        if unknown:
            logger.warning(
                "add_experiment_units(): statische osc_freq-Werte ohne Eintrag in "
                "STATIC_CHIP_LABELS: %s - der Rohwert wird als Chip-Familie benutzt.", unknown,
            )
        out.loc[static, "chip_family"] = freq.map(lambda v: static_chip_labels.get(v, v))
        # osc_freq bleibt Teil der Chip-ID: bei W109 liegen die Medien in
        # getrennten Ordnern (static_omlp/static_ypd), und ob 'Rep1' dort
        # derselbe Chip ist, wissen die Daten nicht. Getrennte Ordner werden
        # deshalb NIE zu einem Chip zusammengelegt; bei W65 (ein Ordner, beide
        # Medien) ist Rep1 ohnehin ein Chip mit zwei Bedingungen.
        out.loc[static, "chip"] = (
            out.loc[static, "chip_family"].astype(str) + "_"
            + out.loc[static, "osc_freq"].astype(str) + "_"
            + out.loc[static, "replicate"].astype(str)
        )
        out.loc[static, "unit_rule"] = UNIT_RULE_STATIC

    osc = ~static
    if osc.any():
        out.loc[osc, "chip"] = (
            out.loc[osc, "biosensor"].astype(str) + "__"
            + out.loc[osc, "osc_type"].astype(str) + "__"
            + out.loc[osc, "osc_freq"].astype(str)
        )
        out.loc[osc, "unit_rule"] = UNIT_RULE_OSC

    n_chips = out["chip"].nunique()
    logger.info(
        "Versuchseinheiten abgeleitet: %d Chips (%d statisch, %d Oszillation/PKO). "
        "Oszillation: EIN Chip pro Bedingung - Fehlerbalken innerhalb einer Bedingung sind "
        "technisch (Kammern eines Chips).",
        n_chips, out.loc[static, "chip"].nunique(), out.loc[osc, "chip"].nunique(),
    )
    return out


def chip_overview(cells: pd.DataFrame) -> pd.DataFrame:
    """Eine Zeile pro Chip: Datum, Kammern je Kontrollart, Zellen, Frames."""
    if "chip" not in cells.columns:
        raise ValueError("chip_overview(): erst add_experiment_units() aufrufen.")
    df = cells.copy()
    df["condition_type"] = df["condition"].map(classify_condition_type)
    # 'medium' gehoert NICHT in den Chip-Schluessel: ein W65-Chip traegt beide
    # Medien, und die Medien erscheinen ohnehin als Spalten (condition_type ist
    # bei statischen Daten 'St.omlp'/'St.ypd').
    key = [c for c in ["chip", "unit_rule", "biosensor", "osc_type", "osc_freq", "chip_family"]
           if c in df.columns]
    # pivot_table verwirft Zeilen mit NaN im Index (dropna=True) - die
    # Oszillations-Chips haben keine chip_family und fielen so heraus.
    if "chip_family" in df.columns:
        df["chip_family"] = df["chip_family"].fillna("")
    per_chamber = df.groupby(key + ["condition_type", "exp_id"], dropna=False).agg(
        n_cells=("cell_uid", "nunique"), n_frames=("frame", "nunique"),
        date=("date", lambda s: ",".join(sorted(set(s.dropna().astype(str)) - {""}))),
    ).reset_index()
    wide = (
        per_chamber.pivot_table(index=key, columns="condition_type", values="exp_id",
                                aggfunc="nunique", fill_value=0)
        .add_prefix("n_chambers_").reset_index()
    )
    totals = per_chamber.groupby(key, dropna=False).agg(
        n_chambers=("exp_id", "nunique"), n_cells=("n_cells", "sum"),
        frames_min=("n_frames", "min"), frames_max=("n_frames", "max"),
        date=("date", lambda s: ",".join(sorted(set(",".join(s).split(",")) - {""}))),
    ).reset_index()
    return wide.merge(totals, on=key, how="left")


def run_order_check(overview: pd.DataFrame) -> pd.DataFrame:
    """Laufreihenfolge gegen Periode: eine Zeile je (biosensor, osc_type)-Serie.

    Jede Periode ist EIN Chip an EINEM Tag. Faellt die Reihenfolge der Tage
    mit der Reihenfolge der Perioden zusammen, kann alles, was ueber die Tage
    driftet (Vorkultur, Fokus, Medium-Charge, Beladungsdichte), einen
    monotonen Trend gegen die Periode erzeugen, der mit der Periode nichts zu
    tun hat - und mit einem Chip pro Periode ist beides nicht zu trennen.
    Spearman(Periode, Datum) ueber die Chips einer Serie: |rho| nahe 1 heisst
    'in Periodenreihenfolge gefahren', nahe 0 'gemischt gefahren'. Bei
    gemischter Reihenfolge kann eine Tagesdrift keinen monotonen
    Periodentrend vortaeuschen.
    """
    needed = {"biosensor", "osc_type", "osc_freq", "date"}
    if overview is None or overview.empty or not needed.issubset(overview.columns):
        return pd.DataFrame()
    from scipy import stats

    df = overview.copy()
    df["_period"] = pd.to_numeric(df["osc_freq"], errors="coerce")
    df["_date"] = pd.to_numeric(df["date"].astype(str).str.split(",").str[0], errors="coerce")
    df = df.dropna(subset=["_period", "_date"])
    rows = []
    for (bs, ot), g in df.groupby(["biosensor", "osc_type"]):
        g = g.sort_values("_period")
        by_date = g.sort_values(["_date", "_period"])
        rec = {
            "biosensor": bs, "osc_type": ot, "n_chips": int(len(g)),
            "periods_in_run_order": " < ".join(f"{p:g}" for p in by_date["_period"]),
            "dates_by_period": " ".join(f"{p:g}:{int(d)}" for p, d in zip(g["_period"], g["_date"])),
            "spearman_period_vs_date": np.nan, "p_value": np.nan, "verdict": "",
        }
        if len(g) >= 3 and g["_date"].nunique() >= 2:
            rho, p = stats.spearmanr(g["_period"], g["_date"])
            rec.update(spearman_period_vs_date=float(rho), p_value=float(p))
            rec["verdict"] = (
                "run in period order: a day-to-day drift would look like a period effect"
                if abs(rho) >= 0.8 else
                "mixed run order: a day-to-day drift cannot mimic a monotone period effect"
            )
        else:
            rec["verdict"] = "fewer than 3 chips or all on one date"
        rows.append(rec)
    out = pd.DataFrame(rows)
    if not out.empty:
        ordered = out[out["spearman_period_vs_date"].abs() >= 0.8]
        if not ordered.empty:
            logger.warning(
                "LAUFREIHENFOLGE: %d von %d Serien wurden in Periodenreihenfolge gefahren (|Spearman(Periode, "
                "Datum)| >= 0.8): %s. Dort ist ein monotoner Periodentrend von einer Tagesdrift nicht zu "
                "unterscheiden (ein Chip pro Periode). Siehe 00_chip_run_order.csv.",
                len(ordered), len(out), ", ".join(ordered["biosensor"] + "/" + ordered["osc_type"]),
            )
        else:
            logger.info("Laufreihenfolge: keine Serie in Periodenreihenfolge gefahren (00_chip_run_order.csv).")
    return out


def summarise_hierarchical(
    df: pd.DataFrame,
    value_col: str,
    condition_cols: Optional[Sequence[str]] = None,
    chamber_col: str = "exp_id",
    chip_col: str = "chip",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Zelle -> Kammer -> Chip -> Bedingung, mit ehrlich beschriftetem Fehler.

    Rueckgabe: (per_chamber, per_chip, per_condition).

    per_condition enthaelt:
        mean            Mittelwert ueber CHIPS (= ueber Kammern, wenn nur ein Chip)
        n_chips         Anzahl biologischer Einheiten
        n_chambers      Anzahl technischer Einheiten
        sd_chip/sem_chip   Streuung ueber Chips (NaN bei n_chips < 2)
        sd_chamber      gepoolte Streuung der Kammern innerhalb der Chips
        error_unit      'chip'  wenn n_chips >= 2, sonst 'chamber'
        sd / sem        der Fehler, der geplottet wird - aus error_unit

    So bekommt eine Oszillationsbedingung (ein Chip) automatisch einen
    Kammer-Fehlerbalken und die Beschriftung 'chamber', eine statische
    Bedingung (4-5 Chips) einen Chip-Fehlerbalken und 'chip'. Nichts wird
    stillschweigend als biologisch ausgegeben, was technisch ist.
    """
    empty = pd.DataFrame()
    if df is None or df.empty or value_col not in df.columns:
        return empty, empty, empty
    if chip_col not in df.columns:
        # Abgeleitete Tabellen (area_table, rp, mu_table ...) tragen die
        # Metadaten, aber nicht immer die Einheiten - dann hier nachziehen,
        # statt jeden Erzeuger kennen zu muessen.
        needed = {"biosensor", "osc_type", "osc_freq", "condition", "replicate"}
        if needed.issubset(df.columns):
            from config import STATIC_CHIP_LABELS, STATIC_MEDIUM_PREFIX
            df = add_experiment_units(df, STATIC_CHIP_LABELS, static_medium_prefix=STATIC_MEDIUM_PREFIX)
        else:
            raise ValueError(
                f"summarise_hierarchical(): Spalte '{chip_col}' fehlt und kann aus "
                f"{sorted(needed - set(df.columns))} nicht abgeleitet werden."
            )

    if condition_cols is None:
        condition_cols = [c for c in ["biosensor", "osc_type", "osc_freq", "condition", "medium", "chip_family"]
                          if c in df.columns]
    condition_cols = [c for c in condition_cols if c in df.columns]
    work = df.dropna(subset=[value_col])
    if work.empty:
        return empty, empty, empty

    chamber_keys = condition_cols + [chip_col] + ([chamber_col] if chamber_col in work.columns else [])
    per_chamber = work.groupby(chamber_keys, dropna=False)[value_col].mean().reset_index()
    per_chamber = per_chamber.rename(columns={value_col: "value"})

    per_chip = (
        per_chamber.groupby(condition_cols + [chip_col], dropna=False)["value"]
        .agg(value="mean", sd_chamber="std", n_chambers="count").reset_index()
    )

    def _sem(s):
        return s.std(ddof=1) / np.sqrt(len(s)) if len(s) > 1 else np.nan

    per_condition = (
        per_chip.groupby(condition_cols, dropna=False)
        .agg(mean=("value", "mean"), sd_chip=("value", "std"), sem_chip=("value", _sem),
             n_chips=("value", "count"), n_chambers=("n_chambers", "sum"),
             sd_chamber=("sd_chamber", lambda s: np.sqrt(np.nanmean(np.square(s))) if s.notna().any() else np.nan))
        .reset_index()
    )
    single = per_condition["n_chips"] < 2
    per_condition["error_unit"] = np.where(single, "chamber", "chip")
    # Bei EINEM Chip: der Kammer-Fehler dieses Chips (sd ueber seine Kammern).
    per_condition["sd"] = np.where(single, per_condition["sd_chamber"], per_condition["sd_chip"])
    per_condition["sem"] = np.where(
        single,
        per_condition["sd_chamber"] / np.sqrt(per_condition["n_chambers"].clip(lower=1)),
        per_condition["sem_chip"],
    )
    # n_units: die Zahl, die ein Leser als 'n' verstehen soll - passend zur error_unit.
    per_condition["n_units"] = np.where(single, per_condition["n_chambers"], per_condition["n_chips"])
    if "condition" in per_condition.columns:
        per_condition["condition_type"] = per_condition["condition"].map(classify_condition_type)
    if "condition" in per_chip.columns:
        per_chip["condition_type"] = per_chip["condition"].map(classify_condition_type)
    for frame in (per_chamber, per_chip, per_condition):
        frame["value_col"] = value_col
    return per_chamber, per_chip, per_condition
