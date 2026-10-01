"""Robustness readings of the results chapter (docs/data_story.md 5b) from the SAVED tables of a finished run,
with the functions of steps 40 and 50 (osc_vs_controls.py, _chamber_mean, control_trend_check); no cell table
needed. From the next run on the pipeline writes the same numbers itself (51_osc_vs_controls*.csv,
40_*_control_trend.csv). Usage: python docs/scratch/robust3.py <analysis_output_v12 directory>
"""
import sys, warnings, logging; warnings.filterwarnings("ignore"); logging.disable(logging.CRITICAL)
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "analyse_pipeline"))
import numpy as np, pandas as pd
from endpoint_trends import summarise_per_replicate, control_trend_check
from robustness import compute_rp
from osc_vs_controls import per_structure, summarise, effect_vs_period
from pipeline_steps import _chamber_mean
from config import ROBUSTNESS_MIN_FRAMES_SINGLE_CELL, ROBUSTNESS_MIN_INTERVALS_MU_EVENT
O = sys.argv[1] if len(sys.argv) > 1 else sys.exit("usage: python robust3.py <analysis_output_v12 directory>")
pd.set_option("display.width", 250)
parts = []
ep = pd.read_csv(f"{O}/13_endpoint_per_chip.csv")
for vc, sub in ep.groupby("value_col"):
    parts.append(per_structure(sub, readout=vc, mode="ratio"))
for name, readout in (("12_area_growth_rate_per_chip.csv", "mu_area"), ("21_budding_rate_per_chip.csv", "budding_rate_per_h"), ("24_growth_from_budding_per_chip.csv", "mu_bud"), ("24_immigration_per_chip.csv", "immigration_per_cell_h")):
    parts.append(per_structure(pd.read_csv(f"{O}/{name}"), readout=readout, mode="ratio"))
# robustness metrics from the raw per-chamber / per-cell tables, as step 40 now does it
metrics = []
for vc in ("area", "eccentricity", "ratio_Queen-2m", "ratio_GlyRNA", "ratio_OxPro", "ratio_pHluorin"):
    metrics.append((f"R(t) population {vc}", _chamber_mean(pd.read_csv(f"{O}/40_Rt_population_{vc}.csv"), "R_t_population")))
    metrics.append((f"R(t) single cell {vc}", _chamber_mean(pd.read_csv(f"{O}/40_Rt_single_cell_{vc}.csv"), "R_t_single_cell", "n_timepoints", ROBUSTNESS_MIN_FRAMES_SINGLE_CELL)))
    metrics.append((f"R(p) {vc}", _chamber_mean(pd.read_csv(f"{O}/40_Rp_{vc}.csv"), "R_p")))
metrics.append(("R(p) mu_area", _chamber_mean(pd.read_csv(f"{O}/40_Rp_mu_area.csv"), "R_p")))
metrics.append(("R(t) single cell mu_event", _chamber_mean(pd.read_csv(f"{O}/40_Rt_single_cell_mu_event.csv"), "R_t_single_cell", "n_timepoints", ROBUSTNESS_MIN_INTERVALS_MU_EVENT)))
pm = pd.read_csv(f"{O}/21_budding_ratio_per_mother.csv", low_memory=False).dropna(subset=["budding_rate_per_h"]).copy(); pm["frame"] = 0
rp_bud = compute_rp(pm, value_col="budding_rate_per_h")
metrics.append(("R(p) budding rate per mother", _chamber_mean(rp_bud, "R_p")))
verdicts = []
for label, cham in metrics:
    per_chip, _ = summarise_per_replicate(cham.dropna(subset=["value"]), "value"); per_chip = per_chip.assign(value_col=label)
    parts.append(per_structure(per_chip, readout=label, mode="diff"))
    ct = control_trend_check(per_chip)
    if not ct.empty:
        ct = ct.assign(short=ct.verdict.str.split(":").str[0])
        verdicts.append(ct)
ps = pd.concat([p for p in parts if not p.empty], ignore_index=True)
summ = summarise(ps)
print("##### 51_osc_vs_controls summary, scope all / per osc_type")
show = summ[summ.scope.isin(["all", "osc_type:Glc", "osc_type:pH"])][["readout", "scope", "n_structures", "n_osc_above_ctrl", "n_above_both", "n_below_both", "wilcoxon_p", "effect_median", "effect_q25", "effect_q75", "osc_median", "posctrl_median", "negctrl_median"]]
print(show.round(3).to_string(index=False))
vd = pd.concat(verdicts, ignore_index=True)
print("\n##### robustness control-trend verdicts per metric")
tab = vd.groupby("value_col").short.value_counts().unstack(fill_value=0)
print(tab.to_string())
print("\ntotal:", vd.short.value_counts().to_dict(), "rows", len(vd))
print("\nperiod effects:", vd[vd.short == "period effect"][["value_col", "biosensor", "osc_type", "rho_osc", "rho_ctrl_strongest", "rho_osc_minus_ctrl"]].round(2).to_string(index=False))
vp = effect_vs_period(ps)
print("\n##### effect vs period: series with |rho| >= 0.6 per readout (main readouts)")
for r in ("budding_rate_per_h", "mu_bud", "area", "mu_area", "eccentricity", "immigration_per_cell_h", "R(p) budding rate per mother", "R(t) single cell mu_event"):
    s = vp[vp.readout == r]; strong = s[s.spearman_rho.abs() >= 0.6]
    print(f"  {r}: {len(strong)} of {len(s)}: " + ", ".join(f"{b}/{o} {rho:+.2f}" for b, o, rho in strong[["biosensor", "osc_type", "spearman_rho"]].itertuples(index=False)))
