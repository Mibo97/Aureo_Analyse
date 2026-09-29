"""Robustness readings for the results chapter (docs/data_story.md 5b): oscillation chambers against the
controls of their structure, R(t)/R(p) levels and control-trend verdicts. Reads the tables of a finished
analysis run; uses the pipeline functions. Usage: python docs/scratch/robust2.py <analysis_output_v12 directory>
"""
import sys, warnings
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "analyse_pipeline")); warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy import stats
from endpoint_trends import add_condition_type, summarise_per_replicate
O = sys.argv[1] if len(sys.argv) > 1 else sys.exit("usage: python robust.py <analysis_output_v12 directory>")
def wide(per_chip, value="value"):
    df = add_condition_type(per_chip) if "condition_type" not in per_chip.columns else per_chip.copy()
    df = df[df.osc_type.isin(["Glc", "pH"]) & (df.biosensor != "PKO")]
    w = df.pivot_table(index=["biosensor", "osc_type", "osc_freq"], columns="condition_type", values=value, aggfunc="mean")
    w["ctrl"] = w[["PosCtrl", "NegCtrl"]].mean(axis=1); w["lo"] = w[["PosCtrl", "NegCtrl"]].min(axis=1); w["hi"] = w[["PosCtrl", "NegCtrl"]].max(axis=1)
    return w.dropna(subset=["Oscillation", "ctrl"])
tabs = {}
for f, vc, name in (("13_endpoint_per_chip.csv", "area", "area"), ("13_endpoint_per_chip.csv", "eccentricity", "eccentricity"), ("21_budding_rate_per_chip.csv", None, "budding rate"), ("24_growth_from_budding_per_chip.csv", "mu_bud", "mu_bud"), ("12_area_growth_rate_per_chip.csv", None, "mu_area"), ("24_immigration_per_chip.csv", None, "immigration")):
    t = pd.read_csv(f"{O}/{f}")
    if vc and "value_col" in t.columns: t = t[t.value_col == vc]
    tabs[name] = wide(t)
for name, w in tabs.items():
    above = int((w.Oscillation > w.hi).sum()); below = int((w.Oscillation < w.lo).sum())
    print(f"\n### {name}: above both controls {above}, below both {below}, between {len(w)-above-below} of {len(w)}")
    for ot, g in w.groupby(level="osc_type"):
        r = g.Oscillation / g.ctrl
        print(f"   {ot}: n={len(g)}, Osc>ctrl {int((g.Oscillation > g.ctrl).sum())}, Wilcoxon p={stats.wilcoxon(g.Oscillation, g.ctrl).pvalue:.3f}, ratio median {r.median():.2f}")
    per_strain = (w.Oscillation / w.ctrl).groupby(level="biosensor").agg(["median", "count"]).round(2)
    gt = (w.Oscillation > w.ctrl).groupby(level="biosensor").sum()
    print("   per strain ratio median:", {s: (per_strain.loc[s, "median"], f"{int(gt[s])}/{int(per_strain.loc[s, 'count'])}") for s in per_strain.index})
# does the oscillation-vs-control difference depend on the period? Spearman of the ratio vs period per series, pooled sign
print("\n### ratio Osc/ctrl vs period: share of series with |rho|>=0.6 and the sign")
for name, w in tabs.items():
    out = []
    for (b, ot), g in w.groupby(level=["biosensor", "osc_type"]):
        per = g.index.get_level_values("osc_freq").astype(float); r = (g.Oscillation / g.ctrl).values
        if len(g) >= 3:
            rho = stats.spearmanr(per, r)[0]; out.append((f"{b}/{ot}", round(rho, 2)))
    strong = [o for o in out if abs(o[1]) >= 0.6]
    print(f"   {name}: {len(strong)} of {len(out)} series |rho|>=0.6: {strong}")
# R(p) mu_area vs R(p) area per structure (Osc): more uniform growth, less uniform size?
rp_a = pd.read_csv(f"{O}/40_Rp_area.csv").groupby("exp_id").R_p.mean(); rp_m = pd.read_csv(f"{O}/40_Rp_mu_area.csv").groupby("exp_id").R_p.mean()
j = pd.concat([rp_a.rename("Rp_area"), rp_m.rename("Rp_mu")], axis=1).dropna()
print("\n### per chamber: Spearman R(p) area vs R(p) mu_area:", round(stats.spearmanr(j.Rp_area, j.Rp_mu)[0], 2), "n", len(j))
# static R values for context
for vc in ("area",):
    for kind, col in (("40_Rt_population_", "R_t_population"), ("40_Rp_", "R_p")):
        t = pd.read_csv(f"{O}/static/{kind}{vc}.csv"); print(f"static {col} {vc} by medium:", t.groupby("medium")[col].mean().round(3).to_dict())
