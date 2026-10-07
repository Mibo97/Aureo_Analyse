"""trend_null.py - Permutationsnull der Kontroll-Trend-Klassifikation (Schritt 50, 50_control_trend_null.csv).

Wie viele Verdikte 'period effect', 'structure effect' und 'not robust' erzeugt die Klassifikation von
endpoint_trends.control_trend_check(), wenn die Periode KEINE Wirkung hat? Die Perioden einer Serie werden
unter ihren Strukturen gemischt - ein Label je Struktur, eine Struktur behaelt ihre Oszillations- und
Kontrollkammern, Chip-Effekte bleiben also erhalten und nur die Periodenordnung wird zerstoert -, die
Klassifikation wird mit denselben Regeln wiederholt, und die Verdikt-Zaehlungen werden ueber viele Mischungen
gesammelt. Lesart: beobachtete Periodeneffekte unter dem Nullmittel brauchen keine Periode zur Erklaerung;
Struktureffekte ueber dem Nullmittel heissen, dass die Strukturdrift in der echten Periodenordnung oefter
monoton war als in einer zufaelligen - die Ausrichtung von Periode, Position und Tag (docs/data_story.md 1.2).
"""
from __future__ import annotations

import logging
import warnings
from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats

from endpoint_trends import MIN_PERIODS_FOR_TREND, add_condition_type, control_trend_check, period_minutes

logger = logging.getLogger(__name__)

VERDICT_KINDS = [
    "no monotone trend of the oscillation chambers",
    "structure effect",
    "not robust",
    "period effect",
    "fewer than 3 periods with oscillation and control values",
]
_GROUP_COLS = ["value_col", "biosensor", "osc_type"]


def _rho(x: np.ndarray, y: np.ndarray, min_periods: int) -> float:
    ok = ~(np.isnan(x) | np.isnan(y))
    if ok.sum() < min_periods or len(np.unique(y[ok])) < 2:
        return np.nan
    return float(stats.spearmanr(x[ok], y[ok])[0])


def classify(rho_osc: float, rho_pos: float, rho_neg: float, rho_mean: float, rho_diff: float,
             strong: float = 0.6) -> str:
    """Dieselben Regeln wie control_trend_check(), auf die Verdikt-Art reduziert."""
    ctrl = [r for r in (rho_pos, rho_neg, rho_mean) if not np.isnan(r)]
    r_ctrl = max(ctrl, key=abs) if ctrl else np.nan
    same_sign_diff = (not np.isnan(rho_diff)) and abs(rho_diff) >= strong and np.sign(rho_diff) == np.sign(rho_osc)
    if np.isnan(rho_osc) or np.isnan(r_ctrl):
        return VERDICT_KINDS[4]
    if abs(rho_osc) < strong:
        return VERDICT_KINDS[0]
    if abs(r_ctrl) >= strong and np.sign(r_ctrl) == np.sign(rho_osc):
        return VERDICT_KINDS[1]
    if same_sign_diff:
        return VERDICT_KINDS[3]
    return VERDICT_KINDS[2]


def _series_arrays(per_chip: pd.DataFrame, value_col: str = "value"):
    """Je (value_col, biosensor, osc_type): Matrix Struktur x (Osc, PosCtrl, NegCtrl) und der Periodenvektor."""
    df = add_condition_type(per_chip) if "condition_type" not in per_chip.columns else per_chip.copy()
    df = df.assign(_period=period_minutes(df["osc_freq"])).dropna(subset=["_period", value_col])
    if df.empty:
        return []
    key_col = "chip" if "chip" in df.columns else "_period"
    group_cols = [c for c in _GROUP_COLS if c in df.columns]
    series = []
    for keys, grp in df.groupby(group_cols, dropna=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        wide = grp.pivot_table(index=key_col, columns="condition_type", values=value_col, aggfunc="mean")
        for ct in ("Oscillation", "PosCtrl", "NegCtrl"):
            if ct not in wide.columns:
                wide[ct] = np.nan
        period = (grp.drop_duplicates(key_col).set_index(key_col)["_period"]
                  .reindex(wide.index).to_numpy(dtype=float))
        series.append((dict(zip(group_cols, keys)), wide[["Oscillation", "PosCtrl", "NegCtrl"]].to_numpy(dtype=float),
                       period))
    return series


def _classify_series(vals: np.ndarray, period: np.ndarray, strong: float, min_periods: int) -> str:
    uniq, inv = np.unique(period, return_inverse=True)
    agg = np.full((len(uniq), 3), np.nan)
    for i in range(len(uniq)):
        block = vals[inv == i]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            agg[i] = np.nanmean(block, axis=0)
    osc, pos, neg = agg[:, 0], agg[:, 1], agg[:, 2]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        ctrl_mean = np.nanmean(np.stack([pos, neg]), axis=0)
    diff = osc - ctrl_mean
    x = uniq.astype(float)
    return classify(_rho(x, osc, min_periods), _rho(x, pos, min_periods), _rho(x, neg, min_periods),
                    _rho(x, ctrl_mean, min_periods), _rho(x, diff, min_periods), strong)


def permutation_null(per_chip: pd.DataFrame, n_perm: int = 1000, seed: int = 0, strong: float = 0.6,
                     min_periods: int = MIN_PERIODS_FOR_TREND, value_col: str = "value") -> pd.DataFrame:
    """Eine Zeile je Verdikt-Art: beobachtete Zahl, Nullmittel, 5./95. Perzentil, Maximum und der Anteil der
    Mischungen mit mindestens der beobachteten Zahl. Leer, wenn keine Serie eine Periode traegt."""
    series = _series_arrays(per_chip, value_col=value_col)
    if not series:
        return pd.DataFrame()
    rng = np.random.default_rng(seed)

    def counts(permute: bool) -> dict:
        c = dict.fromkeys(VERDICT_KINDS, 0)
        for _meta, vals, period in series:
            p = rng.permutation(period) if permute else period
            c[_classify_series(vals, p, strong, min_periods)] += 1
        return c

    observed = counts(False)
    null = pd.DataFrame([counts(True) for _ in range(int(n_perm))])
    rows = []
    for k in VERDICT_KINDS:
        rows.append({
            "verdict": k, "observed": observed[k], "null_mean": float(null[k].mean()),
            "null_q05": float(null[k].quantile(0.05)), "null_q95": float(null[k].quantile(0.95)),
            "null_max": int(null[k].max()), "p_at_least_observed": float((null[k] >= observed[k]).mean()),
            "n_combinations": len(series), "n_perm": int(n_perm), "seed": int(seed), "strong": strong,
        })
    return pd.DataFrame(rows)


def observed_counts_match(per_chip: pd.DataFrame, null_table: pd.DataFrame, value_col: str = "value") -> bool:
    """Selbstpruefung: stimmen die beobachteten Zaehlungen mit control_trend_check() ueberein?"""
    ct = control_trend_check(per_chip, value_col=value_col)
    if ct.empty or null_table.empty:
        return ct.empty and null_table.empty
    kinds = ct["verdict"].astype(str).str.split(":").str[0].str.strip()
    ref = kinds.value_counts().to_dict()
    obs = dict(zip(null_table["verdict"], null_table["observed"]))
    ok = all(int(obs.get(k, 0)) == int(ref.get(k, 0)) for k in VERDICT_KINDS)
    if not ok:
        logger.warning("trend_null: beobachtete Zaehlungen (%s) weichen von control_trend_check (%s) ab.", obs, ref)
    return ok
