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

sns.set_theme(style="whitegrid", context="notebook")

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
    palette = dict(zip(colors, sns.color_palette("viridis", n_colors=len(colors))))

    fig, axes = plt.subplots(
        len(facets), len(groups), figsize=(4.5 * len(groups), 3.5 * len(facets)),
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
                ax.plot(agg[time_col], agg["mean"], color=palette[c], label=str(c), linewidth=1.4)
                ax.fill_between(
                    agg[time_col], agg["mean"] - agg["sd"], agg["mean"] + agg["sd"],
                    color=palette[c], alpha=0.15,
                )

            title_parts = [str(p) for p in (facet, group) if p is not None]
            ax.set_title(" | ".join(title_parts), fontsize=10)
            if last_visible_row_per_col.get(j) == i:
                ax.set_xlabel("Frame")
            if j == 0:
                ax.set_ylabel("Budding Ratio\n(buds/cell)")

    handles, labels = axes[0][0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, title=pretty_label(color_col), loc="lower center",
                   ncol=min(len(colors), 6), bbox_to_anchor=(0.5, -0.05))
    fig.suptitle("Budding ratio over time (mean ± SD over replicates)", y=1.02)
    fig.tight_layout()
    _save_fig(fig, out_path)


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
    Generischer Punkt+Errorbar-Plot über eine Bedingung (z.B. Oszillations-
    frequenz), analog zu Paper Fig. 3b / Fig. 6 (Wachstumsrate, ATP, Robustness).

    Funktioniert mit JEDER Summary-Tabelle, die value_col und (optional)
    eine '<value_col ohne sd>_sd'-artige Spalte enthält - z.B. das Ergebnis
    von summarise_growth_rate() (Spalten mean_mu/sd_mu) oder
    aggregate_robustness_over_replicates() (Spalten mean/sd).

    Parameters
    ----------
    summary : aggregierte Tabelle (eine Zeile pro Bedingung, ggf. + style_col)
    value_col : Spaltenname des Mittelwerts (z.B. 'mean_mu' oder 'mean')
    sd_col : Spaltenname der Standardabweichung (None = automatisch aus
             value_col abgeleitet: 'mean_X' -> 'sd_X', sonst 'mean' -> 'sd')
    x_col : Spalte für die x-Achse (z.B. 'osc_freq')
    facet_col : Spalte für separate Subplots (None = ein einziger Plot)
    color_col : Spalte für die Farbkodierung innerhalb eines Subplots (None = einfarbig)
    style_col : optionale zusätzliche Spalte (z.B. 'condition_type' mit Werten
        'Oscillation'/'PosCtrl'/'NegCtrl'), die als Marker-Form kodiert wird.
        Notwendig, wenn (facet_col, color_col, x_col) die Zeilen NICHT eindeutig
        identifizieren - z.B. weil derselbe osc_freq-Wert sowohl für die
        Oszillationsbedingung als auch für ihre Kontrollkammer(n) vorkommt
        (siehe growth_rate.classify_condition_type()). Ohne style_col würden
        solche Zeilen beim Reindex auf x_order kollidieren (Duplicate-Label-Fehler).
    x_order : Reihenfolge der x-Achse (None = natürliche Frequenzsortierung)
    """
    df = summary.dropna(subset=[value_col]).copy()
    if df.empty:
        logger.warning("plot_point_errorbar(): keine gültigen Werte in '%s' - Plot übersprungen.", value_col)
        return

    if sd_col is None:
        sd_col = value_col.replace("mean_", "sd_") if "mean_" in value_col else (
            "sd" if value_col == "mean" else None
        )
    if sd_col not in df.columns:
        logger.warning("plot_point_errorbar(): sd-Spalte '%s' nicht gefunden - Plot ohne Fehlerbalken.", sd_col)
        sd_col = None

    facets = sorted(df[facet_col].dropna().unique()) if facet_col and facet_col in df.columns else [None]

    # x_values nur auf Kategorien einschränken, die in DIESER Tabelle
    # tatsächlich vorkommen - x_order gibt lediglich die WÜNSCHENSWERTE
    # Reihenfolge vor (z.B. alle möglichen Frequenzen/Bedingungen über den
    # gesamten Datensatz), nicht dass jede davon auch in "summary" vertreten
    # sein muss. Ohne diesen Schritt reserviert reindex() weiter unten für
    # fehlende Kategorien einen leeren x-Tick mit NaN-Werten - das ergibt
    # sichtbar leere Lücken auf der x-Achse.
    present_x = set(df[x_col].dropna().unique())
    x_values_configured = x_order if x_order is not None else natural_freq_sort(df[x_col].dropna().unique())
    x_values = [v for v in x_values_configured if v in present_x]
    if not x_values:
        logger.warning(
            "plot_point_errorbar(): keine der konfigurierten x_order-Kategorien %s kommt in "
            "'%s' vor - Plot übersprungen.", list(x_values_configured), x_col,
        )
        return
    dropped = [v for v in x_values_configured if v not in present_x]
    if dropped:
        logger.info(
            "plot_point_errorbar(): %d von %d konfigurierten x-Kategorien haben keine Daten "
            "in dieser Tabelle und werden nicht als leerer Tick angezeigt: %s", len(dropped), len(x_values_configured), dropped,
        )

    colors = sorted(df[color_col].dropna().unique()) if color_col and color_col in df.columns else [None]
    palette = dict(zip(colors, sns.color_palette("Set2", n_colors=max(len(colors), 1))))

    styles = sorted(df[style_col].dropna().unique()) if style_col and style_col in df.columns else [None]
    marker_map = dict(zip(styles, ["o", "^", "s", "D", "P", "X", "v"]))

    # Breite pro Facette richtet sich nach der ANZAHL DER X-KATEGORIEN, nicht
    # nur nach der Anzahl der Facetten - vorher blieb die Breite pro Facette
    # bei 4in fix, egal ob x_order 2 oder 6 Einträge hatte, wodurch die
    # Datenpunkte bei vielen Frequenzen/Kategorien eng zusammengequetscht
    # wurden und sich Marker/Errorbars sichtbar überlappt haben.
    width_per_facet = max(4.0, 0.75 * len(x_values) + 1.5)
    fig, axes = plt.subplots(1, len(facets), figsize=(width_per_facet * len(facets), 4.5), squeeze=False, sharey=True)
    axes = axes[0]

    n_series = len(colors) * len(styles)
    # Horizontaler Versatz (Dodge) zwischen den Serien (color x style) pro
    # x-Kategorie: statt eines fixen Schritts von 0.08 unabhängig von
    # n_series wird jetzt eine feste GESAMTBREITE (dodge_span) gleichmäßig
    # auf alle Serien verteilt. Dadurch rücken die Punkte bei wenigen Serien
    # nah an die Kategorie-Mitte (wie vorher), aber bei vielen Serien wird
    # der Versatz automatisch enger gehalten statt in die Nachbarkategorie
    # hineinzulaufen - und bei z.B. nur 2 Serien deutlich sichtbarer als
    # die alten pauschalen 0.08.
    dodge_span = 0.5
    dodge_step = dodge_span / n_series if n_series > 0 else 0.0
    # EINE konsolidierte Legende für die ganze Abbildung statt einer Legende
    # pro Facette - sonst verdeckt sie bei mehreren Facetten wiederholt Teile
    # des Plots. Reihenfolge nach erstem Auftreten, Duplikate (gleiches Label
    # in mehreren Facetten) werden zusammengeführt.
    legend_handles: dict[str, object] = {}

    for ax, facet in zip(axes, facets):
        sub = df[df[facet_col] == facet] if facet_col and facet_col in df.columns else df

        series_idx = 0
        for c in colors:
            color_df = sub[sub[color_col] == c] if color_col and color_col in df.columns else sub

            for st in styles:
                line_df = color_df[color_df[style_col] == st] if style_col and style_col in df.columns else color_df
                if line_df.empty:
                    series_idx += 1
                    continue

                # Pro (facet, color, style) muss x_col jetzt eindeutig sein - falls
                # nicht, ist das ein echter Datenfehler, den wir nicht verstecken wollen.
                if line_df[x_col].duplicated().any():
                    dupes = line_df.loc[line_df[x_col].duplicated(keep=False), x_col].unique()
                    raise ValueError(
                        f"plot_point_errorbar(): mehrdeutige '{x_col}'-Werte {list(dupes)} innerhalb "
                        f"({facet_col}={facet!r}, {color_col}={c!r}, {style_col}={st!r}). "
                        f"Erwartet wird eine Zeile pro Bedingung - prüft group_cols der Summary-Tabelle."
                    )
                line_df = line_df.set_index(x_col).reindex(x_values).reset_index()

                x_pos = np.arange(len(x_values)) + (series_idx - n_series / 2 + 0.5) * dodge_step
                yerr = line_df[sd_col] if sd_col else None
                label = " | ".join(str(p) for p in (c, st) if p is not None) or None
                handle = ax.errorbar(
                    x_pos, line_df[value_col], yerr=yerr,
                    fmt=marker_map.get(st, "o"), capsize=3, markersize=6,
                    color=palette[c] if c is not None else None,
                    label=label,
                )
                if label is not None and label not in legend_handles:
                    legend_handles[label] = handle
                series_idx += 1

        ax.set_xticks(np.arange(len(x_values)))
        ax.set_xticklabels(x_values, rotation=30)
        if facet is not None:
            ax.set_title(str(facet), fontsize=11, fontweight="bold")
        ax.set_xlabel(pretty_label(x_col))
        if ax is axes[0]:
            ax.set_ylabel(ylabel or pretty_label(value_col))

    if legend_handles:
        # Außerhalb rechts neben der letzten Facette - bbox_inches="tight" bei
        # savefig() sorgt dafür, dass die Legende trotzdem mit ins PDF kommt.
        fig.legend(
            legend_handles.values(), legend_handles.keys(),
            loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8, borderaxespad=0.0,
        )

    fig.suptitle(title or f"{ylabel or pretty_label(value_col)} per condition (mean ± SD over replicates)", y=1.02)
    fig.tight_layout()
    _save_fig(fig, out_path)


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
    R(t) vs R(p) Quadranten-Plot, analog zu Paper Fig. 9: zeigt für jede
    Bedingung, ob die Funktion über Zeit stabil/instabil UND die Population
    homogen/heterogen ist (vier Quadranten). Ein Punkt = eine Bedingung
    (z.B. eine Kombination aus Biosensor x Oszillationstyp x Frequenz x
    Oscillation/PosCtrl/NegCtrl - je nachdem, was `join_cols` gemeinsam hat).

    Gestrichelte Linien = Mittelwert über ALLE geplotteten Punkte (nicht pro
    Facette) - so bleiben die 4 Quadranten über Facetten hinweg vergleichbar.

    Parameters
    ----------
    rt_summary, rp_summary : aggregierte Tabellen (eine Zeile pro Bedingung),
        z.B. Ergebnis von aggregate_robustness_over_replicates() für
        R_t_population bzw. R_p
    facet_col : Spalte für separate Subplots (Standard 'biosensor', None für
        einen einzigen Plot). Ohne Facettierung landen bei mehreren
        Biosensoren/Kontrollarten schnell viele Punkte in einem einzigen,
        unübersichtlichen Panel übereinander.
    color_col : Spalte für die Punktfarbe (Standard 'condition_type', falls
        vorhanden - unterscheidet Oscillation/PosCtrl/NegCtrl farblich).
        None = einfarbig. Vorher wurden Punkte einfach in Tabellen-
        Reihenfolge eingefärbt (viridis-Verlauf ohne inhaltliche Bedeutung).
    join_cols : Spalten, über die rt_summary und rp_summary zusammengeführt
        werden (Standard: alle gemeinsamen Spalten außer den Werten selbst)
    """
    if join_cols is None:
        join_cols = [c for c in rt_summary.columns if c in rp_summary.columns
                     and c not in ("mean", "sd", "n_replicates")]

    merged = rt_summary.merge(rp_summary, on=join_cols, suffixes=("_rt", "_rp"))
    rt_col = rt_value_col + "_rt" if rt_value_col + "_rt" in merged.columns else rt_value_col
    rp_col = rp_value_col + "_rp" if rp_value_col + "_rp" in merged.columns else rp_value_col

    if merged.empty:
        logger.warning("plot_rt_vs_rp_quadrant(): kein Überlapp zwischen rt_summary und rp_summary - Plot übersprungen.")
        return

    facets = sorted(merged[facet_col].dropna().unique()) if facet_col and facet_col in merged.columns else [None]
    colors_vals = sorted(merged[color_col].dropna().unique()) if color_col and color_col in merged.columns else [None]
    palette = dict(zip(colors_vals, sns.color_palette("Set2", n_colors=max(len(colors_vals), 1))))

    mean_rt, mean_rp = merged[rt_col].mean(), merged[rp_col].mean()

    fig, axes = plt.subplots(1, len(facets), figsize=(5 * len(facets), 5), squeeze=False, sharex=True, sharey=True)
    axes = axes[0]

    legend_handles: dict[str, object] = {}
    for ax, facet in zip(axes, facets):
        sub = merged[merged[facet_col] == facet] if facet_col and facet_col in merged.columns else merged

        ax.axhline(mean_rt, color="gray", linewidth=0.8, linestyle="--")
        ax.axvline(mean_rp, color="gray", linewidth=0.8, linestyle="--")

        for cv in colors_vals:
            sub_c = sub[sub[color_col] == cv] if color_col and color_col in sub.columns else sub
            if sub_c.empty:
                continue
            handle = ax.scatter(
                sub_c[rp_col], sub_c[rt_col],
                color=palette[cv] if cv is not None else None,
                s=70, zorder=5, label=str(cv) if cv is not None else None,
            )
            if cv is not None and str(cv) not in legend_handles:
                legend_handles[str(cv)] = handle
            for _, row in sub_c.iterrows():
                ax.annotate(str(row[label_col]), (row[rp_col], row[rt_col]),
                            textcoords="offset points", xytext=(5, 5), fontsize=7)

        ax.set_xlabel("R(p) (a.u.) — higher = more homogeneous population")
        if ax is axes[0]:
            ax.set_ylabel("R(t) (a.u.) — higher = more stable over time")
        if facet is not None:
            ax.set_title(str(facet), fontsize=11, fontweight="bold")

    if legend_handles:
        fig.legend(
            legend_handles.values(), legend_handles.keys(), title=pretty_label(color_col),
            loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8, borderaxespad=0.0,
        )

    fig.suptitle("Robustness over time vs. over the population", y=1.02)
    fig.tight_layout()
    _save_fig(fig, out_path)


def plot_rt_single_cell_distribution(
    rt_cell: pd.DataFrame,
    out_path: Path,
    value_col: str = "R_t_single_cell",
    facet_col: str = "osc_freq",
    freq_order: Optional[Sequence[str]] = None,
) -> None:
    """
    Verteilung der Einzelzell-R(t)-Werte als KDE pro Bedingung (z.B.
    Oszillationsfrequenz) - analog zu Paper Fig. S9c.

    Im Unterschied zu plot_point_errorbar() (EIN aggregierter Wert pro
    Bedingung) zeigt dieser Plot die VOLLE VERTEILUNG über alle Zellen -
    sinnvoll, weil R(t) auf Einzelzell-Ebene typischerweise sehr heterogen
    ist und ein einzelner Mittelwert das verdecken würde.
    """
    df = rt_cell.dropna(subset=[value_col])
    if df.empty:
        logger.warning("plot_rt_single_cell_distribution(): keine gültigen Werte - Plot übersprungen.")
        return

    facets = freq_order if freq_order is not None else natural_freq_sort(df[facet_col].dropna().unique())
    palette = dict(zip(facets, sns.color_palette("viridis", n_colors=len(facets))))

    fig, ax = plt.subplots(figsize=(7, 5))
    for f in facets:
        sub = df[df[facet_col] == f]
        if sub.empty or sub[value_col].nunique() < 2:
            continue  # KDE braucht Streuung, sonst Fehler/leere Kurve
        sns.kdeplot(sub[value_col], ax=ax, color=palette[f], label=str(f), fill=True, alpha=0.15, linewidth=1.5)

    if not ax.get_legend_handles_labels()[0]:
        logger.warning("plot_rt_single_cell_distribution(): zu wenig Streuung in allen Gruppen - Plot übersprungen.")
        plt.close(fig)
        return

    ax.set_xlabel(pretty_label(value_col))
    ax.set_ylabel("Density")
    ax.set_title(f"Distribution of {pretty_label(value_col)} across single cells, by {pretty_label(facet_col)}")
    ax.legend(title=pretty_label(facet_col), loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8, borderaxespad=0.0)
    fig.tight_layout()
    _save_fig(fig, out_path)
