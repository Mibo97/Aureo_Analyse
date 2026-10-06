"""Are the robustness differences between oscillation and control chambers a property of R or of the
distributions? R = -(sigma^2/x_bar)/m carries the mean; the coefficient of variation sigma/x_bar does not.
Recomputes the paired comparison (oscillation mean of a structure against the mean of its two control types,
Wilcoxon over structures) for R and for the CV from the 40_* tables of an analysis run.
Usage: python robust_cv.py <analysis_output_dir>"""
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import wilcoxon

out = Path(sys.argv[1])


def ctype(cond):
    c = str(cond)
    return "PosCtrl" if "PosCtrl" in c else ("NegCtrl" if "NegCtrl" in c else "Oscillation")


def paired(per_chamber, col, label):
    d = per_chamber.dropna(subset=[col]).copy()
    d["ctype"] = d["condition"].map(ctype)
    d = d[d["osc_type"].isin(["Glc", "pH"])]
    per_type = d.groupby(["chip", "osc_type", "ctype"])[col].mean().unstack("ctype").dropna(subset=["Oscillation"])
    per_type["ctrl"] = per_type[["PosCtrl", "NegCtrl"]].mean(axis=1)
    per_type = per_type.dropna(subset=["ctrl"])
    osc, ctrl = per_type["Oscillation"], per_type["ctrl"]
    print(f"{label:38s} n={len(per_type)} osc>ctrl on {int((osc > ctrl).sum())}  median ratio "
          f"{(osc / ctrl).median():.3f}  p={wilcoxon(osc, ctrl).pvalue:.3f}")
    for ot in ["Glc", "pH"]:
        sub = per_type.xs(ot, level="osc_type")
        print(f"      {ot}: n={len(sub)} osc>ctrl {int((sub['Oscillation'] > sub['ctrl']).sum())} "
              f"ratio {(sub['Oscillation'] / sub['ctrl']).median():.3f} p={wilcoxon(sub['Oscillation'], sub['ctrl']).pvalue:.3f}")


def per_chamber_rp(name):
    t = pd.read_csv(out / f"40_Rp_{name}.csv")
    t = t[t["n_cells"] >= 2].copy()
    t["cv"] = t["sigma"] / t["x_bar"]
    return t.groupby("exp_id").agg(cv=("cv", "mean"), xbar=("x_bar", "mean"), R_p=("R_p", "mean"),
                                   condition=("condition", "first"), chip=("chip", "first"),
                                   osc_type=("osc_type", "first")).reset_index()


for name in ["area", "eccentricity", "mu_area"]:
    ch = per_chamber_rp(name)
    paired(ch, "R_p", f"R(p) {name}")
    paired(ch, "cv", f"CV across cells, {name}")
    paired(ch, "xbar", f"mean {name} per chamber")
rt = pd.read_csv(out / "40_Rt_population_area.csv")
rt["cv"] = rt["sigma"] / rt["x_bar"]
paired(rt, "R_t_population", "R(t) population area")
paired(rt, "cv", "CV over time of chamber mean area")
