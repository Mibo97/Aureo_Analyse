"""Permutation null for the control-trend classification (discussion of 3.6).

How many 'period effects', 'structure effects' and 'not robust' verdicts does the classification of
endpoint_trends.control_trend_check() produce when the period has NO effect? The periods of a series are
shuffled among its structures (one label per structure, so a structure keeps its oscillation chambers and
its controls together, and chip effects survive), the classification is run again, and the verdict counts
are collected over many shuffles. Usage: python trend_null.py <analysis_output_dir> [n_perm]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "analyse_pipeline"))
from endpoint_trends import control_trend_check  # noqa: E402

out = Path(sys.argv[1])
n_perm = int(sys.argv[2]) if len(sys.argv) > 2 else 1000   # 1,000 took about 12 min on the final-run tables; 300 are enough for the means
files = ["13_endpoint_per_chip.csv", "12_area_growth_rate_per_chip.csv", "21_budding_rate_per_chip.csv",
         "24_growth_from_budding_per_chip.csv"]
parts = []
for f in files:
    d = pd.read_csv(out / f)
    if "value_col" not in d.columns:
        d["value_col"] = f.split("_", 1)[1].replace("_per_chip.csv", "")
    parts.append(d)
per_chip = pd.concat(parts, ignore_index=True)
per_chip = per_chip[per_chip["osc_type"].isin(["Glc", "pH"])]


def verdict_counts(table: pd.DataFrame) -> pd.Series:
    ct = control_trend_check(table)
    kind = ct["verdict"].astype(str).str.split(":").str[0].str.strip()
    return kind.value_counts()


observed = verdict_counts(per_chip)
print("observed:", observed.to_dict(), "total", int(observed.sum()))

rng = np.random.default_rng(0)
series_chips = {k: g["chip"].unique() for k, g in per_chip.groupby(["biosensor", "osc_type"])}
chip_period = per_chip.drop_duplicates("chip").set_index("chip")["osc_freq"]
records = []
for i in range(n_perm):
    mapping = {}
    for (b, ot), chips in series_chips.items():
        periods = chip_period.loc[chips].to_numpy()
        mapping.update(dict(zip(chips, rng.permutation(periods))))
    shuffled = per_chip.assign(osc_freq=per_chip["chip"].map(mapping))
    records.append(verdict_counts(shuffled))
null = pd.DataFrame(records).fillna(0)
summary = pd.DataFrame({"observed": observed, "null_mean": null.mean(), "null_q05": null.quantile(0.05),
                        "null_q95": null.quantile(0.95), "null_max": null.max()})
summary["p_at_least_observed"] = [float((null[k] >= observed.get(k, 0)).mean()) if k in null else np.nan
                                  for k in summary.index]
print(summary.round(2).to_string())
summary.to_csv(out / "50_control_trend_null.csv")
print("written", out / "50_control_trend_null.csv")
