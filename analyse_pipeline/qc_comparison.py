"""
qc_comparison.py
================
Was hat das manuelle QC an den Ergebnissen geaendert? Kammer fuer Kammer.

WARUM
-----
Das manuelle QC (Track-Merges, Ausschluesse von Hintergrund, Debris, Zellen am
Rand, toten Zellen) ist arbeitsintensiv - im echten Datensatz ~500 Zeilen fuer
EINEN Batch (WT/pH/6, 11 Kammern). Bevor man das fuer ~50 Batches wiederholt,
sollte man wissen, ob es die Kennzahlen, auf denen die Arbeit steht, ueberhaupt
bewegt. Dieses Modul stellt fuer jeden Chamber die Werte MIT QC (Hauptlauf)
denen OHNE QC (Lauf auf den Rohdaten, Ordner no_qc/) gegenueber.

NUR DORT, WO QC STATTGEFUNDEN HAT
---------------------------------
Der Vergleich ist auf die Batches beschraenkt, die in qc_exclusions.csv
vorkommen (qc_exclusions.qc_batches()). Ueberall sonst waere "mit QC" identisch
mit "ohne QC", und der Vergleich saehe aus wie "QC aendert nichts" - obwohl
schlicht kein QC gemacht wurde.

WAS VERGLICHEN WIRD
-------------------
Vier Kennzahlen, alle auf Kammer-Ebene, weil die Kammer die kleinste Einheit
ist, die beide Laeufe gemeinsam haben:

    endpoint area     13_endpoint_per_chamber.csv   (Flaeche im Endfenster)
    mu_area           12_area_growth_rate_per_chamber.csv (alle Zellen)
    R(p) area         40_Rp_area.csv, pro Kammer gemittelt
    budding ratio     21_budding_ratio_per_experiment.csv (Heuristik!)

Die Abbildung: mit QC auf x, ohne QC auf y, ein Punkt pro Kammer, Farbe =
Kontrollart, Diagonale = kein Effekt. Der Abstand zur Diagonalen IST der
Effekt des QC. Dazu der relative Unterschied je Kammer als Tabelle und das
Inventar der Ausschluesse nach Grund und Kammer - meist sagt allein die Zahl
der 'Background'-Zeilen schon das Meiste.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from growth_rate import classify_condition_type

logger = logging.getLogger(__name__)

CONTROL_COLORS = {"NegCtrl": "#4C78A8", "PosCtrl": "#E45756", "Oscillation": "#333333"}

# (Anzeigename, Datei, Wert-Spalte, optionaler Filter auf value_col)
READOUTS = [
    ("endpoint area [px²]", "13_endpoint_per_chamber.csv", "value", "area"),
    ("µ_area, all cells [h⁻¹]", "12_area_growth_rate_per_chamber.csv", "value", None),
    ("R(p) area", "40_Rp_area.csv", "R_p", None),
    ("budding ratio (heuristic)", "21_budding_ratio_per_experiment.csv", None, None),
]


def _load_readout(directory: Path, fname: str, value_col: Optional[str], value_filter: Optional[str]) -> pd.DataFrame:
    path = directory / fname
    if not path.exists():
        logger.warning("qc_comparison: %s fehlt in %s - Kennzahl uebersprungen.", fname, directory)
        return pd.DataFrame()
    df = pd.read_csv(path)
    if value_filter is not None and "value_col" in df.columns:
        df = df[df["value_col"] == value_filter]
    if "exp_id" not in df.columns:
        logger.warning("qc_comparison: %s hat keine exp_id-Spalte - uebersprungen.", fname)
        return pd.DataFrame()
    if value_col is None:
        # budding ratio per experiment: Spalte heuristisch finden
        cands = [c for c in df.columns if "budding_ratio" in c and df[c].dtype != object]
        if not cands:
            logger.warning("qc_comparison: keine Budding-Ratio-Spalte in %s.", fname)
            return pd.DataFrame()
        value_col = cands[0]
    out = df.groupby("exp_id")[value_col].mean().reset_index(name="value")
    return out


def compare_qc_runs(
    with_dir: Path,
    without_dir: Path,
    batches: pd.DataFrame,
    out_dir: Path,
) -> pd.DataFrame:
    """Kammerweise Gegenueberstellung mit/ohne QC fuer die QC-beruehrten Batches."""
    out_dir.mkdir(parents=True, exist_ok=True)
    if batches is None or batches.empty:
        logger.info("qc_comparison: keine QC-Batches - kein Vergleich.")
        return pd.DataFrame()
    batch_keys = set(map(tuple, batches[["biosensor", "osc_type", "osc_freq"]].astype(str).values))

    def in_qc_batch(exp_id: str) -> bool:
        parts = str(exp_id).split("__")
        return len(parts) >= 3 and tuple(parts[:3]) in batch_keys

    records = []
    for label, fname, vcol, vfilter in READOUTS:
        a = _load_readout(with_dir, fname, vcol, vfilter)
        b = _load_readout(without_dir, fname, vcol, vfilter)
        if a.empty or b.empty:
            continue
        m = a.merge(b, on="exp_id", suffixes=("_with_qc", "_without_qc"), how="outer")
        m = m[m["exp_id"].map(in_qc_batch)]
        if m.empty:
            continue
        m["readout"] = label
        m["condition"] = m["exp_id"].str.split("__").str[3]
        m["condition_type"] = m["condition"].map(classify_condition_type)
        with np.errstate(divide="ignore", invalid="ignore"):
            m["relative_change"] = (m["value_without_qc"] - m["value_with_qc"]) / m["value_with_qc"].abs()
        records.append(m)

    if not records:
        logger.warning("qc_comparison: keine vergleichbaren Kennzahlen gefunden.")
        return pd.DataFrame()
    table = pd.concat(records, ignore_index=True)
    table.to_csv(out_dir / "70_qc_effect_per_chamber.csv", index=False)

    summary = (table.dropna(subset=["relative_change"])
               .groupby(["readout", "condition_type"])["relative_change"]
               .agg(median_rel_change="median", max_abs_rel_change=lambda s: s.abs().max(), n_chambers="count")
               .reset_index())
    summary.to_csv(out_dir / "70_qc_effect_summary.csv", index=False)
    logger.info("QC-Effekt (Median der relativen Aenderung ohne vs. mit QC):\n%s",
                summary.round(3).to_string(index=False))
    plot_qc_comparison(table, out_dir / "70_qc_effect.pdf")
    return table


def qc_exclusion_inventory(exclusions: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    """Was wurde ausgeschlossen bzw. gemergt, nach Grund und Kammer."""
    if exclusions is None or exclusions.empty:
        return pd.DataFrame()
    df = exclusions.copy()
    df["exp_id"] = df["cell_uid"].astype(str).str.replace(r"__track\d+$", "", regex=True)
    df["action"] = np.where(df["merge_into_track_id"].notna(), "merge", "exclude")
    inv = (df.groupby(["exp_id", "action", "reason"]).size().reset_index(name="n_rows")
             .sort_values(["exp_id", "action", "n_rows"], ascending=[True, True, False]))
    out_dir.mkdir(parents=True, exist_ok=True)
    inv.to_csv(out_dir / "70_qc_exclusion_inventory.csv", index=False)
    by_reason = df.groupby(["action", "reason"]).size().sort_values(ascending=False)
    logger.info("QC-Inventar (Zeilen je Aktion/Grund):\n%s", by_reason.to_string())
    return inv


def plot_qc_comparison(table: pd.DataFrame, out_path: Path) -> None:
    readouts = [r for r, *_ in READOUTS if r in set(table["readout"])]
    if not readouts:
        return
    fig, axes = plt.subplots(1, len(readouts), figsize=(3.8 * len(readouts), 3.9), squeeze=False)
    for ax, readout in zip(axes[0], readouts):
        sub = table[table["readout"] == readout].dropna(subset=["value_with_qc", "value_without_qc"])
        for ct, grp in sub.groupby("condition_type"):
            ax.scatter(grp["value_with_qc"], grp["value_without_qc"], s=36, alpha=0.85,
                       color=CONTROL_COLORS.get(ct, "#777777"), edgecolor="white", linewidth=0.5,
                       label=ct, zorder=3)
        lo = np.nanmin([sub["value_with_qc"].min(), sub["value_without_qc"].min()])
        hi = np.nanmax([sub["value_with_qc"].max(), sub["value_without_qc"].max()])
        pad = 0.05 * (hi - lo if hi > lo else 1.0)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color="black", linewidth=0.9,
                linestyle="--", alpha=0.7, zorder=2)
        ax.set_xlim(lo - pad, hi + pad); ax.set_ylim(lo - pad, hi + pad)
        ax.set_aspect("equal", adjustable="box")
        med = sub["relative_change"].median()
        ax.set_title(f"{readout}\nmedian change without QC: {100 * med:+.1f}%", fontsize=9)
        ax.set_xlabel("with manual QC"); ax.set_ylabel("without QC (raw Cellpose output)")
        ax.grid(alpha=0.22, linewidth=0.6)
    handles, labels = axes[0][0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="lower center", ncol=len(labels), frameon=False,
                   bbox_to_anchor=(0.5, -0.06), fontsize=8)
    fig.suptitle("Effect of manual QC, chamber by chamber (QC-reviewed batches only)\n"
                 "point = one chamber; on the diagonal = QC changed nothing", y=1.04)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=180)
    plt.close(fig)
    logger.info("Plot gespeichert: %s", out_path.name)
