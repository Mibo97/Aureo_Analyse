"""
violin_plots.py
================
Violin-Plots im Stil von Panel A des Referenz-Papers: pro Metrik
(Budding Ratio, Area, Circularity, ...) ein Violin je Gruppe (z.B. Stamm),
facettiert nach Oszillationstyp, mit paarweisen Signifikanztests zwischen
den Gruppen (Sternchen-Notation: ns, *, **, ***, ****) und rotem Punkt für
den Mittelwert.

Nutzt `statannotations` für die Signifikanz-Klammern/Sterne, falls
installiert - sonst wird automatisch auf einen einfachen Fallback ohne
Annotationen zurückgefallen (mit deutlicher Warnung, KEIN stiller Verzicht).

WICHTIG zur Statistik: die hier verwendeten Tests (Standard: Mann-Whitney-U)
laufen über die ÜBERGEBENEN ZEILEN als Stichprobe. Wenn ihr Mutterzellen als
Einheit verwenden (siehe lineage.py: compute_budding_ratio() -> per_mother),
ist jede Zeile eine Mutterzelle - das ist näher an einer echten Stichprobe
als rohe Frame-Zeilen, aber immer noch keine vollständige Berücksichtigung
von Replikat-Pseudoreplikation. Für eine strikt korrekte Inferenzstatistik
über Replikate hinweg (z.B. gemischte Modelle) müsstet ihr einen Statistiker
konsultieren - dieses Modul liefert die deskriptive/explorative Variante,
wie sie auch im Referenz-Paper für solche Violin-Übersichten üblich ist.
"""

from __future__ import annotations

import itertools
import logging
from pathlib import Path
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from plot_style import INK, STRAIN_COLOR_OTHER, SURFACE, finish, panel_title, strain_color, tint

logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def _group_palette(df: pd.DataFrame, group_col: str, groups: Sequence) -> dict:
    """Farbe je Gruppe: Staemme in ihrer Stammfarbe; andere Gruppen (statisch: Medium) als
    Toene der einen Stammfarbe der Tabelle, letzte Gruppe voll, davor heller."""
    if group_col == "biosensor":
        return {g: strain_color(g) for g in groups}
    strains = df["biosensor"].dropna().unique() if "biosensor" in df.columns else []
    base = strain_color(strains[0]) if len(strains) == 1 else STRAIN_COLOR_OTHER
    n = len(groups)
    return {g: (base if i == n - 1 else tint(base, 0.55 * (n - 1 - i) / max(n - 1, 1))) for i, g in enumerate(groups)}


try:
    from statannotations.Annotator import Annotator
    _HAS_STATANNOTATIONS = True
except ImportError:
    _HAS_STATANNOTATIONS = False
    logger.warning(
        "Paket 'statannotations' nicht installiert - Violin-Plots werden OHNE "
        "Signifikanz-Sternchen erzeugt. Installieren mit: pip install statannotations"
    )


def plot_violin_with_significance(
    df: pd.DataFrame,
    value_col: str,
    group_col: str,
    out_path: Path,
    facet_col: Optional[str] = None,
    group_order: Optional[Sequence[str]] = None,
    test: str = "Mann-Whitney",
    pairs: Optional[list[tuple]] = None,
    ylabel: Optional[str] = None,
    show_mean: bool = True,
) -> None:
    """
    Violin-Plot mit optionalen paarweisen Signifikanztests, im Stil von
    Panel A des Referenz-Papers (z.B. Budding Ratio / Area / Circularity
    pro Stamm, facettiert nach Oszillationstyp).

    Parameters
    ----------
    df : Datensatz - JEDE ZEILE ist EIN Datenpunkt im Violin (z.B. eine
         Mutterzelle für Budding Ratio, oder eine Einzelzell-Messung für
         Area/Circularity - je nachdem was inhaltlich sinnvoll ist).
    value_col : Spalte mit dem numerischen Wert (y-Achse)
    group_col : Spalte mit der Gruppierung innerhalb einer Facette (x-Achse), z.B. 'strain'
    facet_col : Optional. Spalte für separate Subplots nebeneinander, z.B. 'osc_type'
    group_order : Reihenfolge der Gruppen auf der x-Achse (None = wie in den Daten)
    test : statistischer Test für statannotations (Standard: Mann-Whitney,
           nichtparametrisch - sinnvoller Standard für Verhältniszahlen wie
           Budding Ratio, die nicht normalverteilt sind)
    pairs : Liste von Gruppen-Paaren, die verglichen werden sollen.
            None = ALLE paarweisen Kombinationen der Gruppen.
    ylabel : Achsenbeschriftung (None = value_col)
    show_mean : roter Punkt für den Mittelwert je Gruppe (wie im Paper)
    """
    df = df.dropna(subset=[value_col, group_col])
    if df.empty:
        logger.warning("plot_violin_with_significance(): keine Daten nach dropna für '%s' - Plot übersprungen.", value_col)
        return

    groups = group_order if group_order is not None else list(df[group_col].unique())
    facets = sorted(df[facet_col].dropna().unique()) if facet_col else [None]

    if pairs is None:
        pairs = list(itertools.combinations(groups, 2))

    fig, axes = plt.subplots(1, len(facets), figsize=(2.8 * len(facets) + 0.4, 3.4), squeeze=False, sharey=True)
    axes = axes[0]

    for ax, facet in zip(axes, facets):
        sub = df[df[facet_col] == facet] if facet_col else df

        sns.violinplot(
            data=sub, x=group_col, y=value_col, order=groups,
            hue=group_col, hue_order=groups, legend=False,
            ax=ax, palette=_group_palette(sub, group_col, groups), inner="quartile", cut=0, linewidth=0.8,
        )

        if show_mean:
            means = sub.groupby(group_col)[value_col].mean().reindex(groups)
            ax.scatter(range(len(groups)), means.values, color=INK, edgecolor=SURFACE, linewidth=0.8, zorder=5, s=26)

        # Nur Paare annotieren, für die in DIESER Facette tatsächlich beide
        # Gruppen mit >= 2 Werten vorhanden sind (sonst bricht der Test ab)
        valid_pairs = [
            (g1, g2) for g1, g2 in pairs
            if sub[sub[group_col] == g1][value_col].count() >= 2
            and sub[sub[group_col] == g2][value_col].count() >= 2
        ]

        if _HAS_STATANNOTATIONS and valid_pairs:
            try:
                annot = Annotator(ax, valid_pairs, data=sub, x=group_col, y=value_col, order=groups)
                annot.configure(test=test, text_format="star", verbose=0)
                annot.apply_and_annotate()
            except Exception as e:
                logger.warning(
                    "Signifikanz-Annotation für Facette '%s' fehlgeschlagen (%s: %s) - "
                    "Plot wird ohne Sternchen für diese Facette erzeugt.",
                    facet, type(e).__name__, e,
                )
        elif not _HAS_STATANNOTATIONS:
            pass  # bereits beim Modul-Import gewarnt

        if facet is not None:
            panel_title(ax, str(facet))
        ax.set_xlabel("")
        ax.set_ylabel(ylabel or value_col)
        ax.tick_params(axis="x", rotation=30)

    # Die Sternchen stammen aus einem Mann-Whitney-U ueber die uebergebenen
    # ZEILEN (Zellen bzw. Mutterzellen), nicht ueber biologische Replikate -
    # ihr p haengt damit fast nur an der Zellzahl und ist keine Inferenz
    # ueber Replikate (siehe Modul-Docstring und config.METHOD_CAVEATS). Das
    # gehoert AUF die Abbildung: sie sieht sonst aus wie ein Test, der sie
    # nicht ist.
    # Die Sternchen sind ein Mann-Whitney-U ueber Zellen bzw. Muetter, nicht ueber Replikate -
    # das steht in der Bildunterschrift der Arbeit und in config.METHOD_CAVEATS.
    finish(fig, out_path, logger)


def plot_panel_a(
    cells: pd.DataFrame,
    per_mother: pd.DataFrame,
    out_path: Path,
    group_col: str = "biosensor",
    facet_col: str = "osc_type",
    group_order: Optional[Sequence[str]] = None,
) -> None:
    """
    Erzeugt die drei Panel-A-artigen Violin-Plots (Budding Ratio, Area,
    Circularity) UNTEREINANDER in einer Abbildung, analog zur Referenz-
    Grafik - Budding Ratio aus `per_mother` (eine Zeile pro Mutterzelle,
    siehe lineage.compute_budding_ratio), Area/Circularity aus `cells`
    (eine Zeile pro Zellbeobachtung).

    HINWEIS: 'eccentricity' aus der Cellpose-Pipeline ist NICHT identisch
    mit 'Circularity' aus dem Paper (unterschiedliche math. Definition -
    Exzentrizität misst Ellipsen-Länglichkeit, Circularity typischerweise
    4*pi*Area/Perimeter^2). Hier wird eccentricity als Stellvertreter
    verwendet und auch so benannt, bis ggf. eine echte Circularity-Spalte
    aus der Bildverarbeitung verfügbar ist.
    """
    n_panels = sum([
        not per_mother.empty,
        "area" in cells.columns,
        "eccentricity" in cells.columns,
    ])
    if n_panels == 0:
        logger.warning("plot_panel_a(): weder Budding Ratio noch area/eccentricity verfügbar - Plot übersprungen.")
        return

    facets = sorted(cells[facet_col].dropna().unique())
    groups = group_order if group_order is not None else sorted(cells[group_col].dropna().unique())
    if group_col == "biosensor":
        from plot_style import ordered_strains
        groups = ordered_strains(groups)
    # Nur die Paare gegen den Wildtyp (bzw. Nachbarn, wenn es keinen WT gibt): zehn Klammern ueber fuenf
    # Violinen drueckten die Verteilungen ins untere Drittel der Abbildung.
    if "WT" in groups and len(groups) > 2:
        pairs = [(g, "WT") for g in groups if g != "WT"]
    else:
        pairs = list(zip(groups[:-1], groups[1:])) if len(groups) > 2 else list(itertools.combinations(groups, 2))

    fig, axes = plt.subplots(n_panels, len(facets), figsize=(2.8 * len(facets) + 0.4, 3.0 * n_panels), squeeze=False)

    row = 0
    panel_specs = []
    if not per_mother.empty:
        panel_specs.append(("budding_ratio", per_mother, "Budding Ratio\n(buds/cell)"))
    if "area" in cells.columns:
        panel_specs.append(("area", cells, "Area [px²]"))
    if "eccentricity" in cells.columns:
        panel_specs.append(("eccentricity", cells, "Eccentricity (a.u.)\n[proxy for circularity]"))

    for value_col, data_source, label in panel_specs:
        sub_all = data_source.dropna(subset=[value_col, group_col, facet_col])
        for col_idx, facet in enumerate(facets):
            ax = axes[row][col_idx]
            sub = sub_all[sub_all[facet_col] == facet]

            if sub.empty:
                ax.set_visible(False)
                continue

            sns.violinplot(
                data=sub, x=group_col, y=value_col, order=groups,
                hue=group_col, hue_order=groups, legend=False,
                ax=ax, palette=_group_palette(sub, group_col, groups), inner="quartile", cut=0, linewidth=0.8,
            )
            means = sub.groupby(group_col)[value_col].mean().reindex(groups)
            ax.scatter(range(len(groups)), means.values, color=INK, edgecolor=SURFACE, linewidth=0.8, zorder=5, s=26)

            valid_pairs = [
                (g1, g2) for g1, g2 in pairs
                if sub[sub[group_col] == g1][value_col].count() >= 2
                and sub[sub[group_col] == g2][value_col].count() >= 2
            ]
            if _HAS_STATANNOTATIONS and valid_pairs:
                try:
                    annot = Annotator(ax, valid_pairs, data=sub, x=group_col, y=value_col, order=groups)
                    annot.configure(test="Mann-Whitney", text_format="star", verbose=0)
                    annot.apply_and_annotate()
                except Exception as e:
                    logger.warning("Signifikanz-Annotation fehlgeschlagen für %s/%s (%s): %s", value_col, facet, type(e).__name__, e)

            if row == 0:
                panel_title(ax, str(facet))
            ax.set_xlabel("")
            ax.set_ylabel(label if col_idx == 0 else "")
            ax.tick_params(axis="x", rotation=30)
        row += 1

    # Die Sternchen stammen aus einem Mann-Whitney-U ueber die uebergebenen
    # ZEILEN (Zellen bzw. Mutterzellen), nicht ueber biologische Replikate -
    # ihr p haengt damit fast nur an der Zellzahl und ist keine Inferenz
    # ueber Replikate (siehe Modul-Docstring und config.METHOD_CAVEATS). Das
    # gehoert AUF die Abbildung: sie sieht sonst aus wie ein Test, der sie
    # nicht ist.
    # Die Sternchen sind ein Mann-Whitney-U ueber Zellen bzw. Muetter, nicht ueber Replikate -
    # das steht in der Bildunterschrift der Arbeit und in config.METHOD_CAVEATS.
    finish(fig, out_path, logger)
