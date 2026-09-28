"""Control-validation plots for ratiometric biosensor strains."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Mapping, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analysis import natural_freq_sort
from growth_rate import classify_condition_type

logger = logging.getLogger(__name__)

from plot_style import (CONTROL_LABELS, CONTROL_LINESTYLES, INK, INK_SOFT, control_handles, finish, legend_below,
                        marker_kwargs, panel_title, strain_color, tint)
from matplotlib.lines import Line2D

CONTROL_ORDER = ["NegCtrl", "PosCtrl"]


def _minutes_to_hours(minutes: float) -> float:
    """Convert a minute-valued setting to the hour-valued 'time_h' axis.

    Every ``*_min`` argument in this module is given in MINUTES (that is how the
    experiment is described: "oscillations start after two hours" = 120 min),
    while ``time_h`` produced by ``analysis.add_time_column()`` is in HOURS.
    Comparing the two directly reads as ">= 120 hours" and silently matches
    nothing, so every comparison and every axis position goes through here.
    """
    return minutes / 60.0


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
    """Summarise cell measurements first per replicate, then after the control phase.

    The per-timepoint median is calculated over all cells (and, if present,
    chambers) in one replicate. Thus cells are not treated as independent
    experimental replicates.

    ``analysis_start_min`` is in MINUTES and is converted to hours before it is
    compared against ``time_h`` - see ``_minutes_to_hours()``.
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
    start_h = _minutes_to_hours(analysis_start_min)
    per_replicate_summary = (
        per_replicate_time[per_replicate_time["time_h"] >= start_h]
        .groupby(summary_groups, dropna=False)["median_value"]
        .agg(control_value="median", n_timepoints="count")
        .reset_index()
    )
    if per_replicate_summary.empty and not per_replicate_time.empty:
        logger.warning(
            "No control timepoints at or after %.0f min (= %.2f h); the recording only "
            "covers %.2f-%.2f h. Check OSCILLATION_START_MIN against the actual run length.",
            analysis_start_min, start_h,
            per_replicate_time["time_h"].min(), per_replicate_time["time_h"].max(),
        )
    return per_replicate_time, per_replicate_summary


def plot_sensor_control_timeseries(
    per_replicate_time: pd.DataFrame,
    out_path: Path,
    sensor_label: str,
    preconditioning_end_min: float | None = None,
    strain: Optional[str] = None,
) -> None:
    """Replicate trajectories (thin) plus the median and interquartile range per control type.
    Colour = strain, feast control dashed, famine control dotted."""
    if per_replicate_time.empty:
        logger.warning("Sensor control time series: no data; plot skipped.")
        return
    color = strain_color(strain)
    fig, ax = plt.subplots(figsize=(5.2, 3.2))
    for control in CONTROL_ORDER:
        sub = per_replicate_time[per_replicate_time["condition_type"] == control]
        if sub.empty:
            continue
        ls = CONTROL_LINESTYLES[control]
        line_groups = [c for c in ["osc_freq", "condition", "replicate"] if c in sub.columns]
        for _, rep in sub.groupby(line_groups, dropna=False):
            rep = rep.sort_values("time_h")
            ax.plot(rep["time_h"], rep["median_value"], color=tint(color, 0.5), alpha=0.6, linewidth=0.6)
        aggregate = (
            sub.groupby("time_h")["median_value"]
            .agg(median="median", q25=lambda x: x.quantile(0.25), q75=lambda x: x.quantile(0.75))
            .reset_index().sort_values("time_h")
        )
        ax.plot(aggregate["time_h"], aggregate["median"], color=color, linewidth=2.0, linestyle=ls,
                label=CONTROL_LABELS[control])
        ax.fill_between(aggregate["time_h"], aggregate["q25"], aggregate["q75"], color=color, alpha=0.12, lw=0)

    if preconditioning_end_min is not None:
        end_h = _minutes_to_hours(preconditioning_end_min)
        ax.axvline(end_h, color=INK_SOFT, linestyle="--", linewidth=0.8)
        ax.text(end_h, 1.01, f"{end_h:.0f} h", transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=7.5, color=INK_SOFT)
    ax.set_xlabel("time [h]")
    ax.set_ylabel(f"{sensor_label} ratio (median per chamber)")
    panel_title(ax, f"{sensor_label}: control time course (thin = chambers, thick = median, band = IQR)")
    legend_below(fig, ax.get_legend_handles_labels()[0], ncol=2, y=0.0)
    finish(fig, out_path, logger)


def plot_sensor_control_comparison(
    per_replicate_summary: pd.DataFrame,
    out_path: Path,
    sensor_label: str,
    analysis_start_min: float = 0.0,
    control_labels: Optional[Mapping[str, str]] = None,
    strain: Optional[str] = None,
) -> None:
    """Directly compare per-chamber control values for one sensor (points = chambers, bar = median ± IQR).

    ``control_labels`` optionally adds a second tick-label line per control, e.g.
    ``{"NegCtrl": "0 g/L", "PosCtrl": "50 g/L"}`` for a glucose sensor.
    """
    if per_replicate_summary.empty:
        logger.warning("Sensor control comparison: no data; plot skipped.")
        return
    color = strain_color(strain)
    fig, ax = plt.subplots(figsize=(3.0, 3.2))
    rng = np.random.default_rng(0)
    for x, control in enumerate(CONTROL_ORDER):
        sub = per_replicate_summary[per_replicate_summary["condition_type"] == control]
        if sub.empty:
            continue
        jitter = rng.uniform(-0.10, 0.10, size=len(sub))
        ax.scatter(np.full(len(sub), x) + jitter, sub["control_value"], alpha=0.85, **marker_kwargs(control, color, size=30))
        median = sub["control_value"].median()
        q25, q75 = sub["control_value"].quantile([0.25, 0.75])
        ax.errorbar(x, median, yerr=[[median - q25], [q75 - median]], fmt="_", markersize=18,
                    color=INK, capsize=3, linewidth=1.2, zorder=4)

    ax.set_xticks(range(len(CONTROL_ORDER)))
    ax.set_xticklabels([
        f"{CONTROL_LABELS[c].split(' (')[0]}\n{control_labels[c]}" if control_labels and c in control_labels
        else CONTROL_LABELS[c].split(" (")[0] for c in CONTROL_ORDER
    ])
    start_h = _minutes_to_hours(analysis_start_min)
    time_label = f"after {start_h:.0f} h" if start_h > 0 else "over the full recording"
    ax.set_ylabel(f"{sensor_label} ratio {time_label}\n(median per chamber)")
    panel_title(ax, f"{sensor_label}: control comparison")
    finish(fig, out_path, logger)


def plot_sensor_raw_channel_timeseries(
    channel_a_time: pd.DataFrame,
    channel_b_time: pd.DataFrame,
    out_path: Path,
    sensor_label: str,
    channel_a_label: str,
    channel_b_label: str,
    preconditioning_end_min: float = 120.0,
    strain: Optional[str] = None,
) -> None:
    """Show the two raw channels separately to diagnose a flat ratio."""
    if channel_a_time.empty or channel_b_time.empty:
        logger.warning("Raw-channel diagnostic: missing data; plot skipped.")
        return
    color = strain_color(strain)
    fig, axes = plt.subplots(2, 1, figsize=(5.2, 4.6), sharex=True)
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
            ax.plot(aggregate["time_h"], aggregate["median"], color=color, linewidth=1.8,
                    linestyle=CONTROL_LINESTYLES[control], label=CONTROL_LABELS[control])
            ax.fill_between(aggregate["time_h"], aggregate["q25"], aggregate["q75"], color=color, alpha=0.12, lw=0)
        ax.axvline(_minutes_to_hours(preconditioning_end_min), color=INK_SOFT, linestyle="--", linewidth=0.8)
        ax.set_ylabel(channel_label)
    axes[0].legend(loc="best")
    axes[-1].set_xlabel("time [h]")
    fig.suptitle(f"{sensor_label}: raw control channels (median per chamber, band = IQR)", fontsize=9.5, y=1.0)
    finish(fig, out_path, logger)


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
        per_time[per_time["time_h"] >= _minutes_to_hours(analysis_start_min)]
        .groupby(summary_groups, dropna=False)["median_value"]
        .agg(control_value="median", n_timepoints="count")
        .reset_index()
    )


def plot_control_chamber_comparison(
    chamber_summary: pd.DataFrame,
    out_path: Path,
    sensor_label: str,
    strain: Optional[str] = None,
) -> None:
    """Compare all feast and famine control chambers within every structure (period)."""
    if chamber_summary.empty:
        return
    color = strain_color(strain)
    freq_order = natural_freq_sort(chamber_summary["osc_freq"].dropna().unique())
    fig, axes = plt.subplots(1, len(freq_order), figsize=(1.8 * len(freq_order) + 0.6, 3.0), squeeze=False, sharey=True)
    axes = axes[0]
    rng = np.random.default_rng(0)
    for ax, freq in zip(axes, freq_order):
        sub = chamber_summary[chamber_summary["osc_freq"] == freq]
        for x, control in enumerate(CONTROL_ORDER):
            values = sub.loc[sub["condition_type"] == control, "control_value"]
            jitter = rng.uniform(-0.10, 0.10, size=len(values))
            ax.scatter(np.full(len(values), x) + jitter, values, alpha=0.9, **marker_kwargs(control, color, size=30))
            if not values.empty:
                ax.plot([x - 0.16, x + 0.16], [values.median()] * 2, color=INK, linewidth=1.6, zorder=4)
        ax.set_xticks(range(len(CONTROL_ORDER)))
        ax.set_xticklabels(["famine", "feast"])
        panel_title(ax, f"{freq} min")
        if ax is axes[0]:
            ax.set_ylabel(f"{sensor_label} ratio\n(median per control chamber)")
    legend_below(fig, control_handles(INK_SOFT, which=("PosCtrl", "NegCtrl")), ncol=2, y=0.0)
    fig.suptitle(f"{sensor_label}: control chambers within each structure; bar = median", fontsize=9.5, y=1.0)
    finish(fig, out_path, logger)


