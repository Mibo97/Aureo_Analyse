"""cv_check.py - Variationskoeffizient statt R: ist ein Robustheitsunterschied zwischen Oszillations- und
Kontrollkammern eine Eigenschaft der Verteilungen oder des Masses? R = -(sigma^2/x_bar)/m traegt den
Mittelwert (Fano-Faktor: bei gleicher relativer Streuung sinkt R mit steigendem Mittel), der
Variationskoeffizient sigma/x_bar nicht. Liest die Kammer-Tabellen des Schritts 40 (40_Rp_<readout>.csv je
Kammer x Frame; 40_Rt_population_<readout>.csv je Kammer), bildet je Kammer den CV (R(p): Mittel ueber die
Frames mit >= 2 Zellen; R(t): sigma/x_bar ueber die Zeit) und den Mittelwert (x_bar), aggregiert auf Struktur x
Kontrollart und laeuft durch denselben gepaarten Vergleich wie 51_osc_vs_controls (osc_vs_controls.py).
Ausgaben (Schritt 50): 51_osc_vs_controls_cv_per_structure.csv, 51_osc_vs_controls_cv.csv.
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from endpoint_trends import add_condition_type
from osc_vs_controls import per_structure, summarise

logger = logging.getLogger(__name__)

# (Datei, Label des CV, Label des Mittelwerts, Art)
CV_SPECS = [
    ("40_Rp_area.csv", "CV(p) area", "mean area (chamber, over frames)", "p"),
    ("40_Rp_eccentricity.csv", "CV(p) eccentricity", "mean eccentricity (chamber, over frames)", "p"),
    ("40_Rp_mu_area.csv", "CV(p) mu_area", "mean mu_area (chamber)", "p"),
    ("40_Rt_population_area.csv", "CV(t) population area", "mean area (chamber, over time)", "t"),
    ("40_Rt_population_eccentricity.csv", "CV(t) population eccentricity", "mean eccentricity (chamber, over time)", "t"),
]
_META = ["biosensor", "osc_type", "osc_freq", "condition", "chip", "chip_family", "medium", "culture"]


def chamber_cv(table: pd.DataFrame, kind: str, min_cells: int = 2) -> pd.DataFrame:
    """Eine Zeile je Kammer: cv und level (x_bar), plus die Metadaten."""
    t = table.copy()
    meta = [c for c in _META if c in t.columns]
    if kind == "p" and "n_cells" in t.columns:
        t = t[t["n_cells"] >= min_cells]
    t = t[(t["x_bar"] != 0) & t["x_bar"].notna() & t["sigma"].notna()]
    t["cv"] = t["sigma"] / t["x_bar"]
    if kind == "p":
        agg = {"cv": ("cv", "mean"), "level": ("x_bar", "mean")}
        agg.update({c: (c, "first") for c in meta})
        return t.groupby("exp_id").agg(**agg).reset_index()
    out = t.rename(columns={"x_bar": "level"})
    return out[["exp_id", "cv", "level"] + meta]


def per_chip_from_chambers(chambers: pd.DataFrame, col: str) -> pd.DataFrame:
    """Kammern -> Struktur x Kontrollart (Mittel ueber die Kammern), im Format der *_per_chip.csv."""
    df = add_condition_type(chambers)
    keys = [c for c in ["biosensor", "osc_type", "osc_freq", "condition", "condition_type", "chip", "culture"]
            if c in df.columns]
    g = df.groupby(keys, dropna=False)[col].agg(value="mean", n_chambers="count").reset_index()
    return g


def cv_outputs(output_dir: Path):
    """Liest die 40_-Kammertabellen und gibt (per_structure, summary) des CV-Vergleichs zurueck."""
    output_dir = Path(output_dir)
    parts = []
    for fname, cv_label, level_label, kind in CV_SPECS:
        path = output_dir / fname
        if not path.exists():
            continue
        try:
            t = pd.read_csv(path)
        except pd.errors.EmptyDataError:
            continue
        if t.empty or not {"sigma", "x_bar", "exp_id"} <= set(t.columns):
            continue
        ch = chamber_cv(t, kind)
        if ch.empty:
            continue
        for col, label in (("cv", cv_label), ("level", level_label)):
            ps = per_structure(per_chip_from_chambers(ch, col), readout=label, mode="ratio")
            if ps is not None and not ps.empty:
                parts.append(ps)
    if not parts:
        return pd.DataFrame(), pd.DataFrame()
    per_struct = pd.concat(parts, ignore_index=True)
    return per_struct, summarise(per_struct)
