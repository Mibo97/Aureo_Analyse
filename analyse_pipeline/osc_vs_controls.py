"""osc_vs_controls.py - Oszillationskammern gegen die Kontrollen IHRER Struktur, gepaart ueber Strukturen.

Die Frage der Arbeit ist Robustheit unter Oszillation. Die Kontroll-Trend-Pruefung (endpoint_trends.py)
fragt, ob ein Readout mit der PERIODE laeuft; hier wird gefragt, ob die Oszillation SELBST etwas aendert:
je Struktur der Mittelwert der Oszillationskammern gegen den Mittelwert ihrer Feast- und Famine-Kontrolle
(dieselbe Struktur, dieselbe Kultur, derselbe Tag), ueber alle Strukturen gepaart (Wilcoxon-Vorzeichenrang-
Test). Berichtet werden das Verhaeltnis Osc / Kontrollmittel (Readouts mit positivem Niveau: Flaeche,
Knospungsrate, µ_bud, µ_area, Einwanderung, Sensor-Ratios) oder die Differenz Osc - Kontrollmittel
(Robustheitsmasse R <= 0), insgesamt, je Oszillationstyp und je Stamm; dazu, ob die Oszillationskammern
ueber oder unter BEIDEN Kontrollen liegen, und ob das Verhaeltnis mit der Periode laeuft (Spearman je
Serie, also ob der Oszillationseffekt von der Periode abhaengt).

Vorbehalt: die Oszillationskammern sitzen auf den Positionen A3-A12 eines Arrays, die Kontrollen auf
A1/A2 und A13/A14 - der Vergleich ist auch ein Positionsvergleich. Die Einwanderungsrate
(24_immigration, neue Tracks ohne Elternmaske je Zellstunde) laeuft als Kontrolle des Flusses mit: ist
sie in Oszillations- und Kontrollkammern gleich, ist der Fluss vergleichbar.

Ausgaben (Schritt 50, Praefix 51_): _per_structure.csv (eine Zeile je Struktur und Readout), .csv
(Zusammenfassung je Readout: gesamt, je osc_type, je Stamm), _vs_period.csv (Spearman des Effekts gegen
die Periode je Serie), .pdf / _robustness.pdf / _sensors.pdf (Effekt je Struktur gegen die Periode).
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from scipy import stats

from endpoint_trends import add_condition_type, period_minutes, _log_period_axis
from plot_style import INK_MUTED, finish, legend_below, marker_kwargs, ordered_strains, panel_title, strain_color, strain_handles

logger = logging.getLogger(__name__)

STRUCTURE_KEY = ["biosensor", "osc_type", "osc_freq"]


def per_structure(per_chip: pd.DataFrame, readout: str, value_col: str = "value", mode: str = "ratio") -> pd.DataFrame:
    """Eine Zeile je Struktur: Oszillationsmittel, Feast-, Famine-Kontrolle, Kontrollmittel, Differenz,
    Verhaeltnis (NaN, wenn das Kontrollmittel nicht positiv ist), ob Osc ueber/unter beiden Kontrollen liegt.
    Eingabe ist die Chip-Ebene (ein Wert je Struktur und Kontrollart, z.B. 13_endpoint_per_chip.csv).
    mode: 'ratio' (Readouts mit positivem Niveau) oder 'diff' (Robustheitsmasse R <= 0) - bestimmt nur, welcher
    Effekt in der Zusammenfassung und der Abbildung steht; beide Spalten werden immer berechnet."""
    if per_chip is None or per_chip.empty or value_col not in per_chip.columns:
        return pd.DataFrame()
    df = add_condition_type(per_chip) if "condition_type" not in per_chip.columns else per_chip.copy()
    df = df.assign(_period=period_minutes(df["osc_freq"])).dropna(subset=["_period", value_col])
    if df.empty:
        return pd.DataFrame()
    wide = df.pivot_table(index=STRUCTURE_KEY, columns="condition_type", values=value_col, aggfunc="mean")
    for ct in ("Oscillation", "PosCtrl", "NegCtrl"):
        if ct not in wide.columns:
            wide[ct] = np.nan
    wide = wide.dropna(subset=["Oscillation"])
    wide = wide[wide[["PosCtrl", "NegCtrl"]].notna().any(axis=1)]
    if wide.empty:
        return pd.DataFrame()
    out = wide[["Oscillation", "PosCtrl", "NegCtrl"]].reset_index()
    out.columns.name = None
    out = out.rename(columns={"Oscillation": "osc", "PosCtrl": "posctrl", "NegCtrl": "negctrl"})
    out["ctrl_mean"] = out[["posctrl", "negctrl"]].mean(axis=1)
    out["diff"] = out["osc"] - out["ctrl_mean"]
    with np.errstate(divide="ignore", invalid="ignore"):
        out["ratio"] = np.where(out["ctrl_mean"] > 0, out["osc"] / out["ctrl_mean"], np.nan)
    lo = out[["posctrl", "negctrl"]].min(axis=1)
    hi = out[["posctrl", "negctrl"]].max(axis=1)
    out["above_both"] = out["osc"] > hi
    out["below_both"] = out["osc"] < lo
    out["period_min"] = period_minutes(out["osc_freq"])
    extra = [c for c in ("chip", "culture") if c in df.columns]
    if extra:
        meta = df.drop_duplicates(STRUCTURE_KEY)[STRUCTURE_KEY + extra]
        out = out.merge(meta, on=STRUCTURE_KEY, how="left")
    out.insert(0, "readout", readout)
    out["mode"] = mode
    return out


def _wilcoxon(a: pd.Series, b: pd.Series) -> float:
    ok = a.notna() & b.notna()
    if ok.sum() < 3 or np.allclose(a[ok], b[ok]):
        return np.nan
    try:
        return float(stats.wilcoxon(a[ok], b[ok]).pvalue)
    except ValueError:
        return np.nan


def summarise(per_struct: pd.DataFrame) -> pd.DataFrame:
    """Je Readout und Geltungsbereich (gesamt, je osc_type, je Stamm): n Strukturen, Osc ueber dem
    Kontrollmittel auf n_osc_above_ctrl, ueber/unter beiden Kontrollen, Wilcoxon-p (Osc gegen Kontrollmittel,
    gepaart), Median und Quartile des Effekts (ratio oder diff nach 'mode')."""
    if per_struct is None or per_struct.empty:
        return pd.DataFrame()
    rows = []
    for readout, g in per_struct.groupby("readout", sort=False):
        mode = str(g["mode"].iloc[0])
        scopes = [("all", g)]
        scopes += [(f"osc_type:{ot}", s) for ot, s in g.groupby("osc_type")]
        scopes += [(f"biosensor:{b}", s) for b, s in g.groupby("biosensor")]
        for scope, s in scopes:
            eff = s["ratio"] if mode == "ratio" else s["diff"]
            rows.append({
                "readout": readout, "mode": mode, "scope": scope, "n_structures": int(len(s)),
                "n_osc_above_ctrl": int((s["osc"] > s["ctrl_mean"]).sum()),
                "n_above_both": int(s["above_both"].sum()), "n_below_both": int(s["below_both"].sum()),
                "wilcoxon_p": _wilcoxon(s["osc"], s["ctrl_mean"]),
                "effect_median": float(eff.median()) if eff.notna().any() else np.nan,
                "effect_q25": float(eff.quantile(0.25)) if eff.notna().any() else np.nan,
                "effect_q75": float(eff.quantile(0.75)) if eff.notna().any() else np.nan,
                "osc_median": float(s["osc"].median()), "posctrl_median": float(s["posctrl"].median()),
                "negctrl_median": float(s["negctrl"].median()),
            })
    return pd.DataFrame(rows)


def effect_vs_period(per_struct: pd.DataFrame, min_periods: int = 3) -> pd.DataFrame:
    """Haengt der Oszillationseffekt von der Periode ab? Spearman des Effekts (ratio oder diff) gegen die
    Periode je Readout und Serie (biosensor x osc_type), n = Perioden."""
    if per_struct is None or per_struct.empty:
        return pd.DataFrame()
    rows = []
    for (readout, b, ot), g in per_struct.groupby(["readout", "biosensor", "osc_type"], sort=False):
        mode = str(g["mode"].iloc[0])
        eff = g["ratio"] if mode == "ratio" else g["diff"]
        ok = eff.notna() & g["period_min"].notna()
        rho = p = np.nan
        if ok.sum() >= min_periods and eff[ok].nunique() > 1:
            rho, p = stats.spearmanr(g.loc[ok, "period_min"], eff[ok])
        rows.append({"readout": readout, "mode": mode, "biosensor": b, "osc_type": ot,
                     "n_periods": int(ok.sum()), "spearman_rho": float(rho), "p_value": float(p)})
    return pd.DataFrame(rows)


def plot_osc_vs_controls(per_struct: pd.DataFrame, out_path: Path, readouts: Sequence[str],
                         labels: Optional[dict] = None, title: Optional[str] = None) -> None:
    """Effekt je Struktur gegen die Periode: eine Zeile je Readout, eine Spalte je Oszillationstyp; Punkte in
    der Stammfarbe. 'ratio'-Readouts auf log2-Achse mit Linie bei 1, 'diff'-Readouts linear mit Linie bei 0.
    Paneltitel: Median des Effekts, Wilcoxon-p und n Strukturen des Panels."""
    if per_struct is None or per_struct.empty:
        return
    labels = labels or {}
    readouts = [r for r in readouts if r in set(per_struct["readout"])]
    if not readouts:
        return
    osc_types = sorted(per_struct["osc_type"].dropna().astype(str).unique())
    fig, axes = plt.subplots(len(readouts), len(osc_types),
                             figsize=(3.3 * len(osc_types) + 0.5, 1.75 * len(readouts) + 0.8),
                             squeeze=False, sharex="col")
    periods_all = sorted(per_struct["period_min"].dropna().unique())
    strains = ordered_strains(per_struct["biosensor"].dropna().unique())
    for i, readout in enumerate(readouts):
        sub_r = per_struct[per_struct["readout"] == readout]
        mode = str(sub_r["mode"].iloc[0])
        for j, ot in enumerate(osc_types):
            ax = axes[i][j]
            sub = sub_r[sub_r["osc_type"].astype(str) == ot]
            eff = sub["ratio"] if mode == "ratio" else sub["diff"]
            # log2-Achse: Verhaeltnisse <= 0 (Oszillationskammern ohne Ereignis) sind nicht zeichenbar und
            # bleiben weg; Median und p im Titel rechnen mit ihnen.
            eff_plot = eff.where(eff > 0) if mode == "ratio" else eff
            ax.axhline(1.0 if mode == "ratio" else 0.0, color=INK_MUTED, linewidth=0.8, linestyle=(0, (4, 2)), zorder=1)
            for strain in strains:
                m = (sub["biosensor"] == strain) & eff_plot.notna()
                if m.any():
                    ax.scatter(sub.loc[m, "period_min"], eff_plot[m], **marker_kwargs("Oscillation", strain_color(strain), size=26))
            periods = sorted(sub["period_min"].dropna().unique()) or periods_all
            _log_period_axis(ax, periods)
            if mode == "ratio":
                ax.set_yscale("log", base=2)
                ax.yaxis.set_major_locator(mticker.LogLocator(base=2.0))
                ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}"))
                ax.yaxis.set_minor_locator(mticker.NullLocator())
            p = _wilcoxon(sub["osc"], sub["ctrl_mean"])
            med = eff.median() if eff.notna().any() else np.nan
            p_txt = "p < 0.001" if p < 0.001 else f"p = {p:.2f}" if np.isfinite(p) else "p n.a."
            head = f"{ot}: " if i == 0 else ""
            panel_title(ax, f"{head}median {med:.2f}, {p_txt}, n = {len(sub)}")
            if j == 0:
                ax.set_ylabel(labels.get(readout, readout) + ("\n(oscillation / controls)" if mode == "ratio"
                                                             else "\n(oscillation − controls)"), fontsize=7.5)
            if i == len(readouts) - 1:
                ax.set_xlabel("cycle period [min]")
    legend_below(fig, strain_handles(strains), ncol=min(len(strains), 6), y=0.0)
    if title:
        fig.suptitle(title, fontsize=9.5, y=1.0)
    finish(fig, out_path, logger)
