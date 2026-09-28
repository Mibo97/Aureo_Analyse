"""
summary_plots.py
=================
Plots für die "höheren" Auswertungsschritte (ab Schritt 05 in run_analysis.py):
Budding-Ratio-Zeitreihe, spezifische Wachstumsrate, und Robustness R(t)/R(p).

Alle Funktionen nehmen die bereits AGGREGIERTEN Tabellen entgegen (Ergebnis
der jeweiligen aggregate_*/summarise_*-Funktion aus budding_ratio_timeseries.py,
growth_rate.py, robustness.py) - nicht die Rohtabellen pro Zelle/Frame.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

from plot_style import (CONTROL_LABELS, INK_MUTED, INK_SOFT, MEDIUM_FILLED, SURFACE, axis_label, control_handles, edge_color, errorbar_kwargs,
                        finish, legend_below, marker_kwargs, ordered_strains, panel_title, strain_color,
                        strain_handles, strain_ramp)
from matplotlib.lines import Line2D

try:
    from analysis import natural_freq_sort, pretty_label  # falls dort schon definiert
except ImportError:
    import re

    def natural_freq_sort(values: Sequence[str]) -> list[str]:
        """Sortiert Frequenz-Strings wie '2min','10min' numerisch statt alphabetisch."""
        def key(v: str):
            m = re.search(r"\d+(\.\d+)?", str(v))
            return (0, float(m.group())) if m else (1, str(v))
        return sorted(set(values), key=key)

    def pretty_label(col: str) -> str:
        """Fallback ohne analysis.py - siehe analysis.pretty_label()."""
        return col.replace("_", " ")


def _save_fig(fig: plt.Figure, out_path: Path) -> None:
    """
    Speichert eine Figure als PDF und schließt sie danach IMMER (auch bei
    Fehlern), um Speicherlecks über einen langen Pipeline-Lauf mit vielen
    Plots zu vermeiden.

    Fängt PermissionError/OSError ab (z.B. Windows verweigert den Zugriff,
    weil die PDF-Datei gerade in einem Viewer wie Adobe/Edge geöffnet ist,
    oder OneDrive/Antivirus sie kurz sperrt) und loggt das nur als Fehler,
    statt die GESAMTE restliche Auswertung abzubrechen - ein einzelner
    gesperrter Output sollte nicht verhindern, dass alle anderen Plots noch
    entstehen.
    """
    try:
        fig.savefig(out_path, bbox_inches="tight", dpi=150)
        logger.info("Plot gespeichert: %s", out_path.name)
    except (PermissionError, OSError) as exc:
        logger.error(
            "Plot konnte NICHT gespeichert werden (Datei evtl. gerade geöffnet/gesperrt?): "
            "%s - %s", out_path, exc,
        )
    finally:
        plt.close(fig)


def plot_budding_ratio_timeseries(
    timeseries: pd.DataFrame,
    out_path: Path,
    group_col: str = "biosensor",
    facet_col: str = "osc_type",
    color_col: str = "osc_freq",
    freq_order: Optional[Sequence[str]] = None,
    time_col: str = "frame",
) -> None:
    """
    Budding Ratio über die Zeit, gemittelt über Replikate (± SD), analog zu
    Paper Fig. 3a (Linienplot). Eine Spalte pro group_col-Wert (z.B. Biosensor),
    eine Zeile pro facet_col-Wert (z.B. Oszillationstyp), Farbe = color_col
    (z.B. Oszillationsfrequenz).

    Erwartet die ROHE Zeitreihe (eine Zeile pro exp_id x frame, siehe
    budding_ratio_timeseries.compute_budding_ratio_timeseries) - die Mittelung
    über Replikate passiert intern in dieser Funktion.
    """
    df = timeseries.dropna(subset=["budding_ratio"])
    if df.empty:
        logger.warning("plot_budding_ratio_timeseries(): keine gültigen budding_ratio Werte - Plot übersprungen.")
        return

    groups = sorted(df[group_col].dropna().unique()) if group_col in df.columns else [None]
    facets = sorted(df[facet_col].dropna().unique()) if facet_col in df.columns else [None]
    colors = freq_order if freq_order is not None else natural_freq_sort(df[color_col].dropna().unique())
    if group_col == "biosensor":
        groups = ordered_strains(groups)

    fig, axes = plt.subplots(
        len(facets), len(groups), figsize=(2.9 * len(groups) + 0.4, 2.6 * len(facets) + 0.4),
        squeeze=False, sharex=True,
    )

    # Per column (group), find the LOWEST row (facet) that actually carries
    # data - only that panel gets the "Frame" x-label. This has to be decided
    # from the data, not from ax.get_visible(): at this point nothing has been
    # hidden yet, so the old check always returned the last row.
    def _panel_data(facet, group):
        sub = df
        if facet_col in df.columns:
            sub = sub[sub[facet_col] == facet]
        if group_col in df.columns:
            sub = sub[sub[group_col] == group]
        return sub

    last_visible_row_per_col = {}
    for i, facet in enumerate(facets):
        for j, group in enumerate(groups):
            if not _panel_data(facet, group).empty:
                last_visible_row_per_col[j] = i

    for i, facet in enumerate(facets):
        for j, group in enumerate(groups):
            ax = axes[i][j]
            sub = _panel_data(facet, group)

            if sub.empty:
                ax.set_visible(False)
                continue

            # Hell-Dunkel-Rampe der Stammfarbe ueber die Perioden (kurz hell, lang dunkel).
            ramp = dict(zip(colors, strain_ramp(group if group_col == "biosensor" else "_", len(colors))))
            for c in colors:
                line_df = sub[sub[color_col] == c]
                if line_df.empty:
                    continue
                agg = (
                    line_df.groupby(time_col)["budding_ratio"]
                    .agg(mean="mean", sd="std")
                    .reset_index()
                    .sort_values(time_col)
                )
                ax.plot(agg[time_col], agg["mean"], color=ramp[c], label=str(c), linewidth=1.3)
                ax.fill_between(
                    agg[time_col], agg["mean"] - agg["sd"], agg["mean"] + agg["sd"],
                    color=ramp[c], alpha=0.12, lw=0,
                )

            title_parts = [str(p) for p in (facet, group) if p is not None]
            panel_title(ax, " | ".join(title_parts))
            if last_visible_row_per_col.get(j) == i:
                ax.set_xlabel("frame")
            if j == 0:
                ax.set_ylabel("budding ratio\n(buds per cell)")

    grey = dict(zip(colors, strain_ramp("_", len(colors))))
    handles = [Line2D([], [], color=grey[c], linewidth=2, label=str(c)) for c in colors]
    legend_below(fig, handles, ncol=min(len(colors), 6), y=0.0,
                 title=axis_label(color_col) + " (light = short, dark = long; hue = strain)")
    finish(fig, out_path, logger)


def plot_point_errorbar(
    summary: pd.DataFrame,
    value_col: str,
    out_path: Path,
    x_col: str = "osc_freq",
    facet_col: Optional[str] = "osc_type",
    color_col: Optional[str] = "biosensor",
    style_col: Optional[str] = None,
    sd_col: Optional[str] = None,
    x_order: Optional[Sequence[str]] = None,
    ylabel: Optional[str] = None,
    title: Optional[str] = None,
) -> None:
    """
    Generischer Punkt+Errorbar-Plot ueber eine Bedingung (z.B. die Periode), analog zu
    Paper Fig. 3b / Fig. 6 (Wachstumsrate, Robustness).

    Funktioniert mit JEDER Summary-Tabelle, die value_col und (optional) eine sd-Spalte enthaelt -
    z.B. summarise_growth_rate() (mean_mu/sd_mu) oder aggregate_robustness_over_replicates()
    (mean/sd). Farbe = color_col ('biosensor': Stammfarbe aus config.STRAIN_COLORS, sonst Grau);
    style_col (z.B. 'condition_type' mit Oscillation/PosCtrl/NegCtrl) = Marker und Fuellung
    (plot_style). style_col ist noetig, wenn (facet_col, color_col, x_col) die Zeilen nicht
    eindeutig machen, etwa weil derselbe osc_freq-Wert fuer die Oszillationsbedingung und ihre
    Kontrollkammern vorkommt. x_order = Reihenfolge der x-Achse (None = natuerliche Sortierung).
    """
    if summary is None or summary.empty or value_col not in summary.columns:
        logger.warning(
            "plot_point_errorbar(): leere Tabelle oder Spalte '%s' fehlt - Plot '%s' uebersprungen.",
            value_col, out_path.name,
        )
        return

    df = summary.dropna(subset=[value_col]).copy()
    if df.empty:
        logger.warning("plot_point_errorbar(): keine gueltigen Werte in '%s' - Plot uebersprungen.", value_col)
        return

    if sd_col is None:
        sd_col = value_col.replace("mean_", "sd_") if "mean_" in value_col else (
            "sd" if value_col == "mean" else None
        )
    if sd_col not in df.columns:
        logger.warning("plot_point_errorbar(): sd-Spalte '%s' nicht gefunden - Plot ohne Fehlerbalken.", sd_col)
        sd_col = None

    facets = sorted(df[facet_col].dropna().unique()) if facet_col and facet_col in df.columns else [None]

    # x-Kategorien nur, wenn sie in DIESER Tabelle vorkommen; x_order gibt die Reihenfolge vor.
    present_x = set(df[x_col].dropna().unique())
    x_values_configured = x_order if x_order is not None else natural_freq_sort(df[x_col].dropna().unique())
    x_values = [v for v in x_values_configured if v in present_x]
    if not x_values:
        logger.warning(
            "plot_point_errorbar(): keine der konfigurierten x_order-Kategorien %s kommt in "
            "'%s' vor - Plot uebersprungen.", list(x_values_configured), x_col,
        )
        return

    has_color = bool(color_col) and color_col in df.columns
    colors = (ordered_strains(df[color_col].dropna().unique()) if has_color and color_col == "biosensor"
              else sorted(df[color_col].dropna().unique()) if has_color else [None])
    if has_color and color_col == "biosensor":
        palette = {c: strain_color(c) for c in colors}
    else:
        palette = dict(zip(colors, strain_ramp("_", max(len(colors), 1))))

    has_style = bool(style_col) and style_col in df.columns
    styles = sorted(df[style_col].dropna().unique()) if has_style else [None]

    width_per_facet = max(2.6, 0.55 * len(x_values) + 1.2)
    fig, axes = plt.subplots(1, len(facets), figsize=(width_per_facet * len(facets) + 0.4, 3.2),
                             squeeze=False, sharey=True)
    axes = axes[0]

    n_series = len(colors) * len(styles)
    dodge_span = 0.5
    dodge_step = dodge_span / n_series if n_series > 0 else 0.0

    x_values_all = x_values
    for ax, facet in zip(axes, facets):
        sub = df[df[facet_col] == facet] if facet_col and facet_col in df.columns else df
        # x-Kategorien je Facette: nur die, die in DIESER Facette vorkommen (keine leeren Ticks).
        present_f = set(sub[x_col].dropna().unique())
        x_values = [v for v in x_values_all if v in present_f] or x_values_all
        series_idx = 0
        for c in colors:
            color_df = sub[sub[color_col] == c] if has_color else sub
            for st in styles:
                line_df = color_df[color_df[style_col] == st] if has_style else color_df
                if line_df.empty:
                    series_idx += 1
                    continue
                if line_df[x_col].duplicated().any():
                    dupes = line_df.loc[line_df[x_col].duplicated(keep=False), x_col].unique()
                    raise ValueError(
                        f"plot_point_errorbar(): mehrdeutige '{x_col}'-Werte {list(dupes)} innerhalb "
                        f"({facet_col}={facet!r}, {color_col}={c!r}, {style_col}={st!r}). "
                        f"Erwartet wird eine Zeile pro Bedingung - prueft group_cols der Summary-Tabelle."
                    )
                line_df = line_df.set_index(x_col).reindex(x_values).reset_index()
                x_pos = np.arange(len(x_values)) + (series_idx - n_series / 2 + 0.5) * dodge_step
                yerr = line_df[sd_col] if sd_col else None
                kw = errorbar_kwargs(st if st is not None else "Oscillation", palette.get(c, INK_SOFT))
                kw["linestyle"] = ""
                if x_col == "medium":
                    # statisch: komplexes Medium gefuellt, Minimalmedium hohl (plot_style.MEDIUM_FILLED)
                    for k, xv in enumerate(x_values):
                        kw_k = dict(kw, markerfacecolor=kw["color"] if MEDIUM_FILLED.get(str(xv), True) else SURFACE)
                        ax.errorbar(x_pos[k:k + 1], line_df[value_col].iloc[k:k + 1],
                                    yerr=None if yerr is None else yerr.iloc[k:k + 1], **kw_k)
                else:
                    ax.errorbar(x_pos, line_df[value_col], yerr=yerr, **kw)
                series_idx += 1

        ax.set_xticks(np.arange(len(x_values)))
        ax.set_xticklabels([str(v) for v in x_values], rotation=30 if len(x_values) > 4 else 0)
        if facet is not None:
            panel_title(ax, str(facet))
        ax.set_xlabel(axis_label(x_col))
        if ax is axes[0]:
            ax.set_ylabel(ylabel or pretty_label(value_col))

    handles = []
    if has_color and len(colors) > 1:
        handles += strain_handles(colors) if color_col == "biosensor" else [
            Line2D([], [], marker="o", linestyle="", markersize=6, markerfacecolor=palette[c],
                   markeredgecolor=edge_color(palette[c]), label=str(c)) for c in colors]
    if has_style and len(styles) > 1:
        handles += control_handles(INK_SOFT, which=[s for s in ("Oscillation", "PosCtrl", "NegCtrl") if s in styles])
    if handles:
        legend_below(fig, handles, ncol=min(len(handles), 4), y=0.0)
    if title:
        fig.suptitle(title.split("\n")[0], fontsize=9.5, y=1.0)
    finish(fig, out_path, logger)


def plot_rt_vs_rp_quadrant(
    rt_summary: pd.DataFrame,
    rp_summary: pd.DataFrame,
    out_path: Path,
    rt_value_col: str = "mean",
    rp_value_col: str = "mean",
    label_col: str = "osc_freq",
    facet_col: Optional[str] = "biosensor",
    color_col: Optional[str] = "condition_type",
    join_cols: Optional[Sequence[str]] = None,
) -> None:
    """
    R(t) vs R(p) Quadranten-Plot (Paper Fig. 9): pro Bedingung, ob die Funktion ueber die Zeit
    stabil UND die Population homogen ist. Ein Punkt = eine Bedingung; Farbe = Stamm (Facette),
    Marker = Kontrollart (color_col 'condition_type'), Beschriftung = label_col (Periode).
    Gestrichelte Linien = Mittelwert ueber ALLE Punkte (ueber Facetten vergleichbar).
    """
    if join_cols is None:
        join_cols = [c for c in rt_summary.columns if c in rp_summary.columns
                     and c not in ("mean", "sd", "n_replicates")]

    merged = rt_summary.merge(rp_summary, on=join_cols, suffixes=("_rt", "_rp"))
    rt_col = rt_value_col + "_rt" if rt_value_col + "_rt" in merged.columns else rt_value_col
    rp_col = rp_value_col + "_rp" if rp_value_col + "_rp" in merged.columns else rp_value_col

    if merged.empty:
        logger.warning("plot_rt_vs_rp_quadrant(): kein Ueberlapp zwischen rt_summary und rp_summary - Plot uebersprungen.")
        return

    facets = (ordered_strains(merged[facet_col].dropna().unique()) if facet_col == "biosensor"
              else sorted(merged[facet_col].dropna().unique())) if facet_col and facet_col in merged.columns else [None]
    has_style = bool(color_col) and color_col in merged.columns
    mean_rt, mean_rp = merged[rt_col].mean(), merged[rp_col].mean()

    fig, axes = plt.subplots(1, len(facets), figsize=(3.0 * len(facets) + 0.4, 3.2), squeeze=False,
                             sharex=True, sharey=True)
    axes = axes[0]
    styles_seen: set = set()
    for ax, facet in zip(axes, facets):
        sub = merged[merged[facet_col] == facet] if facet_col and facet_col in merged.columns else merged
        color = strain_color(facet) if facet_col == "biosensor" else INK_SOFT
        ax.axhline(mean_rt, color=INK_MUTED, linewidth=0.7, linestyle="--")
        ax.axvline(mean_rp, color=INK_MUTED, linewidth=0.7, linestyle="--")
        for st, part in (sub.groupby(color_col) if has_style else [(None, sub)]):
            styles_seen.add(st)
            ax.scatter(part[rp_col], part[rt_col], **marker_kwargs(st if st is not None else "Oscillation", color, size=34))
            for _, row in part.iterrows():
                ax.annotate(str(row[label_col]), (row[rp_col], row[rt_col]), textcoords="offset points",
                            xytext=(4, 4), fontsize=6.5, color=INK_SOFT)
        ax.set_xlabel("R(p): higher = more homogeneous population")
        if ax is axes[0]:
            ax.set_ylabel("R(t): higher = more stable over time")
        if facet is not None:
            panel_title(ax, str(facet))

    handles = control_handles(INK_SOFT, which=[s for s in ("Oscillation", "PosCtrl", "NegCtrl") if s in styles_seen])
    if handles:
        legend_below(fig, handles, ncol=3, y=0.0)
    finish(fig, out_path, logger)


def plot_rt_single_cell_distribution(
    rt_cell: pd.DataFrame,
    out_path: Path,
    value_col: str = "R_t_single_cell",
    facet_col: str = "osc_freq",
    freq_order: Optional[Sequence[str]] = None,
) -> None:
    """
    Verteilung der Einzelzell-R(t)-Werte als KDE pro Bedingung (z.B. Periode), analog zu Paper
    Fig. S9c: die VOLLE Verteilung ueber alle Zellen statt eines Mittelwerts. Perioden als
    Hell-Dunkel-Rampe (kurz hell, lang dunkel; alle Staemme gepoolt, daher grau).
    """
    df = rt_cell.dropna(subset=[value_col])
    if df.empty:
        logger.warning("plot_rt_single_cell_distribution(): keine gueltigen Werte - Plot uebersprungen.")
        return

    facets = freq_order if freq_order is not None else natural_freq_sort(df[facet_col].dropna().unique())
    palette = dict(zip(facets, strain_ramp("_", len(facets))))

    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    for f in facets:
        sub = df[df[facet_col] == f]
        if sub.empty or sub[value_col].nunique() < 2:
            continue
        sns.kdeplot(sub[value_col], ax=ax, color=palette[f], label=str(f), fill=True, alpha=0.12, linewidth=1.4)

    if not ax.get_legend_handles_labels()[0]:
        logger.warning("plot_rt_single_cell_distribution(): zu wenig Streuung in allen Gruppen - Plot uebersprungen.")
        plt.close(fig)
        return

    ax.set_xlabel(pretty_label(value_col))
    ax.set_ylabel("density")
    ax.legend(title=axis_label(facet_col), loc="upper right", fontsize=7.5)
    finish(fig, out_path, logger)


