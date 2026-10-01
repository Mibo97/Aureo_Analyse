"""Robustness readings for the results chapter (docs/data_story.md 5b): oscillation chambers against the
controls of their structure, R(t)/R(p) levels and control-trend verdicts. Reads the tables of a finished
analysis run; uses the pipeline functions. Usage: python docs/scratch/robust.py <analysis_output_v12 directory>
"""
import sys, warnings
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "analyse_pipeline"))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy import stats
from endpoint_trends import control_trend_check, summarise_per_replicate, add_condition_type
O = sys.argv[1] if len(sys.argv) > 1 else sys.exit("usage: python robust.py <analysis_output_v12 directory>")
pd.set_option("display.width", 250)

def per_structure(per_chip, value="value"):
    """wide: one row per (biosensor, osc_type, osc_freq/chip): Osc, PosCtrl, NegCtrl."""
    df = add_condition_type(per_chip) if "condition_type" not in per_chip.columns else per_chip.copy()
    df = df[df.osc_type.isin(["Glc", "pH"]) & (df.biosensor != "PKO")]
    w = df.pivot_table(index=["biosensor", "osc_type", "osc_freq"], columns="condition_type", values=value, aggfunc="mean")
    w["ctrl"] = w[["PosCtrl", "NegCtrl"]].mean(axis=1)
    return w.dropna(subset=["Oscillation", "ctrl"])

def osc_vs_ctrl(w, label, ratio=True):
    d = w["Oscillation"] - w["ctrl"]
    p = stats.wilcoxon(w["Oscillation"], w["ctrl"]).pvalue
    n_between = int(((w["Oscillation"] >= w[["PosCtrl", "NegCtrl"]].min(axis=1)) & (w["Oscillation"] <= w[["PosCtrl", "NegCtrl"]].max(axis=1))).sum())
    line = f"{label:22s} n={len(w):2d} Osc>ctrl on {int((d > 0).sum()):2d}, Wilcoxon p={p:.2f}, within bracket {n_between:2d}"
    if ratio:
        r = w["Oscillation"] / w["ctrl"]
        line += f", Osc/ctrl median {r.median():.2f} (q25 {r.quantile(.25):.2f}, q75 {r.quantile(.75):.2f}); per osc_type: " + ", ".join(f"{k} {v:.2f}" for k, v in r.groupby(level="osc_type").median().items())
    else:
        line += f", Osc-ctrl median {d.median():+.3f} (q25 {d.quantile(.25):+.3f}, q75 {d.quantile(.75):+.3f}); levels Osc {w['Oscillation'].median():.3f} Pos {w['PosCtrl'].median():.3f} Neg {w['NegCtrl'].median():.3f}"
    print(line)

print("##### A. readouts: oscillation chambers relative to the controls of their structure (per-chip tables)")
for f, vc in (("13_endpoint_per_chip.csv", "area"), ("13_endpoint_per_chip.csv", "eccentricity"), ("21_budding_rate_per_chip.csv", None), ("24_growth_from_budding_per_chip.csv", "mu_bud"), ("24_immigration_per_chip.csv", None), ("12_area_growth_rate_per_chip.csv", None)):
    t = pd.read_csv(f"{O}/{f}")
    if vc and "value_col" in t.columns: t = t[t.value_col == vc]
    osc_vs_ctrl(per_structure(t), vc or f.split("_per_chip")[0])

print("\n##### B. robustness R(t) population / R(p) / R(t) single cell: levels, Osc vs controls, control-trend verdicts")
rows = []
for kind, files in (("R(t) population", [("area", "R_t_population"), ("eccentricity", "R_t_population"), ("ratio_Queen-2m", "R_t_population"), ("ratio_GlyRNA", "R_t_population"), ("ratio_OxPro", "R_t_population"), ("ratio_pHluorin", "R_t_population")]),
                    ("R(p)", [("area", "R_p"), ("eccentricity", "R_p"), ("mu_area", "R_p"), ("ratio_Queen-2m", "R_p"), ("ratio_GlyRNA", "R_p"), ("ratio_OxPro", "R_p"), ("ratio_pHluorin", "R_p")]),
                    ("R(t) single cell", [("area", "R_t_single_cell"), ("eccentricity", "R_t_single_cell"), ("ratio_Queen-2m", "R_t_single_cell"), ("ratio_GlyRNA", "R_t_single_cell"), ("ratio_OxPro", "R_t_single_cell"), ("ratio_pHluorin", "R_t_single_cell")])):
    for vc, col in files:
        stem = {"R_t_population": "40_Rt_population_", "R_p": "40_Rp_", "R_t_single_cell": "40_Rt_single_cell_"}[col]
        t = pd.read_csv(f"{O}/{stem}{vc}.csv")
        if col == "R_t_single_cell" and "n_timepoints" in t.columns:
            t = t[t.n_timepoints >= 10]
        cham = t.groupby("exp_id").agg(value=(col, "mean")).reset_index().merge(t.drop_duplicates("exp_id")[["exp_id", "biosensor", "osc_type", "osc_freq", "condition", "chip"]], on="exp_id")
        per_chip, summ = summarise_per_replicate(cham, "value")
        w = per_structure(per_chip)
        print(f"--- {kind} {vc}:", end=" "); osc_vs_ctrl(w, "", ratio=False)
        ct = control_trend_check(per_chip.assign(value_col=vc))
        if not ct.empty:
            vd = ct.verdict.str.split(":").str[0].value_counts().to_dict()
            print(f"      control-trend over {len(ct)} series: {vd}")
            for _, r in ct.iterrows():
                if abs(r.rho_osc) >= 0.6:
                    print(f"        {r.biosensor}/{r.osc_type}: osc {r.rho_osc:+.2f}, strongest ctrl {r.rho_ctrl_strongest:+.2f}, diff {r.rho_osc_minus_ctrl:+.2f} -> {r.verdict.split(':')[0]}")
        rows.append(dict(kind=kind, readout=vc, n=len(w), osc=w.Oscillation.median(), pos=w.PosCtrl.median(), neg=w.NegCtrl.median(), p=stats.wilcoxon(w.Oscillation, w.ctrl).pvalue))
print("\n##### C. summary table")
print(pd.DataFrame(rows).round(3).to_string(index=False))
print("\n##### D. R(p) spearman tables of the pipeline (Osc only), for cross-check")
for vc in ("area", "eccentricity", "mu_area", "ratio_OxPro", "ratio_pHluorin", "ratio_GlyRNA", "ratio_Queen-2m"):
    s = pd.read_csv(f"{O}/40_Rp_{vc}_spearman.csv"); print(vc, [(r.biosensor, r.osc_type, round(r.spearman_rho, 2)) for _, r in s.iterrows() if abs(r.spearman_rho) >= 0.6])
