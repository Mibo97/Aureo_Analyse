"""Control-validation plots for ratiometric biosensor strains."""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analysis import natural_freq_sort
from growth_rate import classify_condition_type

logger = logging.getLogger(__name__)

CONTROL_ORDER = ["NegCtrl", "PosCtrl"]
CONTROL_COLORS = {"NegCtrl": "#4C78A8", "PosCtrl": "#E45756"}


def prepare_sensor_controls(
    cells: pd.DataFrame,
    value_col: str,
    biosensor: str,
    osc_type: str,
) -> pd.DataFrame:
    """Return positive- and negative-control measurements for one sensor."""
    required = {"biosensor", "osc_type", "condition", "replicate", "time_h", value_col}
    missing = required.difference(cells.columns)
    if missing:
        raise ValueError(f"prepare_sensor_controls(): missing columns: {sorted(missing)}")

    df = cells[(cells["biosensor"] == biosensor) & (cells["osc_type"] == osc_type)].copy()
    df["condition_type"] = df["condition"].map(classify_condition_type)
    df = df[df["condition_type"].isin(CONTROL_ORDER)].dropna(subset=[value_col])
    if df.empty:
        logger.warning("No control data found for biosensor=%r and osc_type=%r.", biosensor, osc_type)
    return df


def summarise_sensor_controls(
    control_cells: pd.DataFrame,
    value_col: str,
    analysis_start_min: float = 120.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Summarise cell measurements first per replicate, then after 120 min.

    The per-timepoint median is calculated over all cells (and, if present,
    chambers) in one replicate. Thus cells are not treated as independent
    experimental replicates.
    """
    if control_cells.empty:
        return pd.DataFrame(), pd.DataFrame()

    group_cols = [
        c for c in ["biosensor", "osc_type", "osc_freq", "condition", "condition_type", "replicate", "time_h"]
        if c in control_cells.columns
    ]
    per_replicate_time = (
        control_cells.groupby(group_cols, dropna=False)[value_col]
        .agg(median_value="median", n_cell_measurements="count")
        .reset_index()
    )
    summary_groups = [c for c in group_cols if c != "time_h"]
    per_replicate_summary = (
        per_replicate_time[per_replicate_time["time_h"] >= analysis_start_min]
        .groupby(summary_groups, dropna=False)["median_value"]
        .agg(control_value="median", n_timepoints="count")
        .reset_index()
    )
    return per_replicate_time, per_replicate_summary


def plot_sensor_control_timeseries(
    per_replicate_time: pd.DataFrame,
    out_path: Path,
    sensor_label: str,
    preconditioning_end_min: float | None = None,
) -> None:
    """Plot replicate trajectories plus the median and interquartile range."""
    if per_replicate_time.empty:
        logger.warning("Sensor control time series: no data; plot skipped.")
        return

    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for control in CONTROL_ORDER:
        sub = per_replicate_time[per_replicate_time["condition_type"] == control]
        if sub.empty:
            continue
        line_groups = [c for c in ["osc_freq", "condition", "replicate"] if c in sub.columns]
        for _, rep in sub.groupby(line_groups, dropna=False):
            rep = rep.sort_values("time_h")
            ax.plot(rep["time_h"], rep["median_value"], color=CONTROL_COLORS[control], alpha=0.18, linewidth=0.9)

        aggregate = (
            sub.groupby("time_h")["median_value"]
            .agg(median="median", q25=lambda x: x.quantile(0.25), q75=lambda x: x.quantile(0.75))
            .reset_index().sort_values("time_h")
        )
        ax.plot(aggregate["time_h"], aggregate["median"], color=CONTROL_COLORS[control], linewidth=2.4, label=control)
        ax.fill_between(aggregate["time_h"], aggregate["q25"], aggregate["q75"], color=CONTROL_COLORS[control], alpha=0.18)

    if preconditioning_end_min is not None:
        ax.axvline(preconditioning_end_min, color="0.35", linestyle="--", linewidth=1.0)
        ax.text(preconditioning_end_min, 1.01, "2 h", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=8)
    ax.set_xlabel("Zeit [min]")
    ax.set_ylabel(f"{sensor_label}-Ratio")
    ax.set_title(f"{sensor_label}: Kontrollzeitreihe")
    ax.legend(title="Kontrolle", frameon=True)
    ax.text(0.01, -0.20, "Dünne Linien = Replikate; dicke Linie = Median; Band = Interquartilsabstand.",
            transform=ax.transAxes, fontsize=8, va="top")
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=180)
    plt.close(fig)
    logger.info("Plot saved: %s", out_path.name)


def plot_sensor_control_comparison(
    per_replicate_summary: pd.DataFrame,
    out_path: Path,
    sensor_label: str,
    analysis_start_min: float = 0.0,
) -> None:
    """Directly compare per-replicate control values for one sensor."""
    if per_replicate_summary.empty:
        logger.warning("Sensor control comparison: no data; plot skipped.")
        return

    fig, ax = plt.subplots(figsize=(4.8, 4.8))
    rng = np.random.default_rng(0)
    for x, control in enumerate(CONTROL_ORDER):
        sub = per_replicate_summary[per_replicate_summary["condition_type"] == control]
        if sub.empty:
            continue
        jitter = rng.uniform(-0.10, 0.10, size=len(sub))
        ax.scatter(np.full(len(sub), x) + jitter, sub["control_value"], s=38, alpha=0.8,
                   color=CONTROL_COLORS[control], edgecolor="white", linewidth=0.5, zorder=3)
        median = sub["control_value"].median()
        q25, q75 = sub["control_value"].quantile([0.25, 0.75])
        ax.errorbar(x, median, yerr=[[median - q25], [q75 - median]], fmt="_", markersize=20,
                    color="black", capsize=4, linewidth=1.4, zorder=4)

    ax.set_xticks(range(len(CONTROL_ORDER)))
    ax.set_xticklabels(["NegCtrl\n0 g/L", "PosCtrl\n50 g/L"])
    time_label = "nach 2 h" if analysis_start_min == 120 else "über gesamte Messdauer"
    ax.set_ylabel(f"{sensor_label}-Ratio {time_label}\n(Median pro Replikat)")
    ax.set_title(f"{sensor_label}: Kontrollvergleich")
    ax.text(0.01, -0.18, "Punkte = Replikate; schwarzer Balken = Median ± Interquartilsabstand.",
            transform=ax.transAxes, fontsize=8, va="top")
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=180)
    plt.close(fig)
    logger.info("Plot saved: %s", out_path.name)


def plot_sensor_raw_channel_timeseries(
    channel_a_time: pd.DataFrame,
    channel_b_time: pd.DataFrame,
    out_path: Path,
    sensor_label: str,
    channel_a_label: str,
    channel_b_label: str,
    preconditioning_end_min: float = 120.0,
) -> None:
    """Show the two raw channels separately to diagnose a flat ratio."""
    if channel_a_time.empty or channel_b_time.empty:
        logger.warning("Raw-channel diagnostic: missing data; plot skipped.")
        return

    fig, axes = plt.subplots(2, 1, figsize=(8.0, 7.0), sharex=True)
    for ax, data, channel_label in zip(axes, [channel_a_time, channel_b_time], [channel_a_label, channel_b_label]):
        for control in CONTROL_ORDER:
            sub = data[data["condition_type"] == control]
            aggregate = (
                sub.groupby("time_h")["median_value"]
                .agg(median="median", q25=lambda x: x.quantile(0.25), q75=lambda x: x.quantile(0.75))
                .reset_index().sort_values("time_h")
            )
            if aggregate.empty:
                continue
            ax.plot(aggregate["time_h"], aggregate["median"], color=CONTROL_COLORS[control], linewidth=2.0, label=control)
            ax.fill_between(aggregate["time_h"], aggregate["q25"], aggregate["q75"], color=CONTROL_COLORS[control], alpha=0.18)
        ax.axvline(preconditioning_end_min, color="0.35", linestyle="--", linewidth=1.0)
        ax.set_ylabel(channel_label)
        ax.legend(title="Kontrolle", frameon=True)
    axes[-1].set_xlabel("Zeit [min]")
    fig.suptitle(f"{sensor_label}: Rohkanäle der Kontrollen", y=0.98)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=180)
    plt.close(fig)
    logger.info("Plot saved: %s", out_path.name)


def summarise_sensor_controls_by_chamber(
    control_cells: pd.DataFrame,
    value_col: str,
    analysis_start_min: float = 120.0,
) -> pd.DataFrame:
    """Summarise every control chamber for descriptive within-chip comparisons.

    Chambers are technical measurements from one chip, not independent
    biological replicates. Keeping them separate makes PosCtrl and NegCtrl
    inspectable within a frequency batch without inventing chamber pairs.
    """
    if control_cells.empty or "chamber" not in control_cells.columns:
        return pd.DataFrame()
    group_cols = [
        c for c in ["biosensor", "osc_type", "osc_freq", "condition", "condition_type", "replicate", "chamber", "time_h"]
        if c in control_cells.columns
    ]
    per_time = (
        control_cells.groupby(group_cols, dropna=False)[value_col]
        .agg(median_value="median", n_cell_measurements="count")
        .reset_index()
    )
    summary_groups = [c for c in group_cols if c != "time_h"]
    return (
        per_time[per_time["time_h"] >= analysis_start_min]
        .groupby(summary_groups, dropna=False)["median_value"]
        .agg(control_value="median", n_timepoints="count")
        .reset_index()
    )


def plot_control_chamber_comparison(
    chamber_summary: pd.DataFrame,
    out_path: Path,
    sensor_label: str,
) -> None:
    """Compare all PosCtrl and NegCtrl chambers within every frequency batch."""
    if chamber_summary.empty:
        return
    freq_order = natural_freq_sort(chamber_summary["osc_freq"].dropna().unique())
    fig, axes = plt.subplots(1, len(freq_order), figsize=(3.0 * len(freq_order), 4.6), squeeze=False, sharey=True)
    axes = axes[0]
    rng = np.random.default_rng(0)
    for ax, freq in zip(axes, freq_order):
        sub = chamber_summary[chamber_summary["osc_freq"] == freq]
        for x, control in enumerate(CONTROL_ORDER):
            values = sub.loc[sub["condition_type"] == control, "control_value"]
            jitter = rng.uniform(-0.10, 0.10, size=len(values))
            ax.scatter(np.full(len(values), x) + jitter, values, color=CONTROL_COLORS[control], s=38,
                       alpha=0.85, edgecolor="white", linewidth=0.5, zorder=3)
            if not values.empty:
                ax.plot([x - 0.16, x + 0.16], [values.median()] * 2, color="black", linewidth=1.8)
        ax.set_xticks(range(len(CONTROL_ORDER)))
        ax.set_xticklabels(["NegCtrl", "PosCtrl"], rotation=25, ha="right")
        ax.set_title(f"{freq} min", fontsize=10)
        if ax is axes[0]:
            ax.set_ylabel(f"{sensor_label}-Ratio\n(Median pro Kontrollkammer)")
    fig.suptitle(f"{sensor_label}: Kontrollen innerhalb jedes Frequenz-Batches", y=1.02)
    fig.text(0.5, -0.02, "Punkte = technische Kontrollkammern desselben Chips; schwarzer Strich = Median.",
             ha="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=180)
    plt.close(fig)
    logger.info("Plot saved: %s", out_path.name)
