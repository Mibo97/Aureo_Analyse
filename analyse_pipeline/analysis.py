"""
analysis.py
===========
Auswertung über die volle Experiment-Hierarchie:
    Biosensor x Oszillationstyp (Glucose/pH) x Oszillationsfrequenz (6 Stufen)
    x Replikat x Kammer x Zelle x Frame

Statistik wird korrekt zuerst PRO REPLIKAT gemittelt und erst dann über
Replikate aggregiert (Mittelwert ± SEM) - alles andere wäre Pseudo-
replikation (einzelne Zellen sind keine unabhängigen biologischen Replikate).
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

from plot_style import (CONTROL_LABELS, CONTROL_LINESTYLES, INK, INK_MUTED, INK_SOFT, SURFACE, edge_color, finish,
                        legend_below, ordered_strains, panel_title, strain_color, strain_ramp, tint)
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


# Anzeigenamen für Spalten, die in generischen Plots direkt als Achsen-/
# Legendentitel landen (plot_point_errorbar(), plot_rt_single_cell_distribution(),
# ...). Ohne diese Tabelle steht dort der rohe Spaltenname ("osc_freq").
# Nur Anzeige - die Spaltennamen in den CSV-Exporten ändern sich NICHT.
COLUMN_LABELS: dict[str, str] = {
    "osc_freq": "Oscillation frequency",
    "osc_type": "Oscillation type",
    "biosensor": "Biosensor",
    "condition": "Condition",
    "condition_type": "Condition type",
    "replicate": "Replicate",
    "chamber": "Chamber",
    "frame": "Frame",
    "time_h": "Time [h]",
    "area": "Cell area [px²]",
    "eccentricity": "Eccentricity (a.u.)",
    "solidity": "Solidity (a.u.)",
    "budding_ratio": "Budding ratio (buds/cell)",
    "n_tracks": "Number of tracks",
    "mu": "Specific growth rate µ [h⁻¹]",
    "mu_area": "µ_area [h⁻¹]",
    "R_t_population": "R(t), population level",
    "R_t_single_cell": "R(t), single cell",
    "R_p": "R(p)",
}


def pretty_label(col: str) -> str:
    """Anzeigename für einen Spaltennamen (siehe COLUMN_LABELS).

    Unbekannte Spalten fallen auf eine harmlose Aufhübschung zurück
    (Unterstriche zu Leerzeichen), ratio_*-Spalten auf "<Sensor> ratio" -
    so bleibt der Plot lesbar, auch wenn eine Spalte hier noch nicht
    eingetragen ist.
    """
    if col in COLUMN_LABELS:
        return COLUMN_LABELS[col]
    if col.startswith("ratio_"):
        return f"{col[len('ratio_'):]} ratio"
    return col.replace("_", " ")


def natural_freq_sort(values: Sequence[str]) -> list[str]:
    """
    Sortiert Frequenz-Strings wie '2min','5min','10min','60min' NUMERISCH
    statt alphabetisch (alphabetisch würde '10min' vor '2min' einordnen).
    Extrahiert die erste Zahl im String als Sortierschlüssel; Strings ohne
    erkennbare Zahl werden ans Ende gehängt (alphabetisch unter sich).
    """
    def key(v: str):
        m = re.search(r"\d+(\.\d+)?", str(v))
        return (0, float(m.group())) if m else (1, str(v))

    return sorted(set(values), key=key)


def add_time_column(df: pd.DataFrame, min_per_frame: float) -> pd.DataFrame:
    """Fügt `time_h` basierend auf `frame` und der Zeitkalibrierung hinzu."""
    df = df.copy()
    df["time_h"] = df["frame"] * min_per_frame / 60
    return df


def find_intensity_columns(df: pd.DataFrame, prefix: str = "mean_") -> list[str]:
    """Findet alle Intensitäts-Spalten (Standard: mean_<kanal>)."""
    return [c for c in df.columns if c.startswith(prefix)]


def data_overview(df: pd.DataFrame) -> pd.DataFrame:
    """Übersicht: Anzahl Tracks/Frames pro Hierarchie-Ebene.
    Erster Sanity-Check nach QC - deckt auf, wenn z.B. ein Replikat fast
    komplett rausgeflogen ist."""
    overview = (
        df.groupby(["biosensor", "osc_type", "osc_freq", "condition", "replicate", "chamber"])
        .agg(
            n_tracks=("track_id", "nunique"),
            n_frames=("frame", "nunique"),
            n_rows=("track_id", "size"),
        )
        .reset_index()
    )
    return overview


def plot_n_tracks_overview(
    overview: pd.DataFrame,
    out_path: Path,
    freq_order=None,
    ctrl_pattern: str = r"(NegCtrl|PosCtrl)",
) -> None:
    """
    Erstellt ZWEI PDFs:
      1. <out_path>          – ein Balken pro Replikat/Kammer, kein CI (QC-Detail)
      2. <stem>_summary.pdf  – Replikate aggregiert (Mittelwert ± STD),
                               NegCtrl/PosCtrl als schraffierte Extra-Balken
    
    NegCtrl/PosCtrl werden anhand von ctrl_pattern in der 'condition'-Spalte
    erkannt (default: erkennt '...NegCtrl...' und '...PosCtrl...').
    """
    order = freq_order if freq_order is not None else natural_freq_sort(
        overview["osc_freq"].dropna().unique()
    )

    # ------------------------------------------------------------------
    # Plot 1: Detail – ein Balken pro Replikat/Kammer (unveraendert)
    # ------------------------------------------------------------------
    ov = overview.copy()

    g = sns.catplot(
        data=ov,
        x="osc_freq", y="n_tracks", hue="replicate",
        row="osc_type", col="biosensor",
        kind="bar", order=order, errorbar=None,
        height=3.0, aspect=1.2,
        palette=dict(zip(sorted(ov["replicate"].dropna().unique()), strain_ramp("_", ov["replicate"].nunique()))),
    )
    g.set_axis_labels("Oscillation frequency", "Number of tracks")
    g.set_titles("{row_name} | {col_name}")
    for ax in g.axes.flat:
        ax.tick_params(axis="x", rotation=45)
    g.fig.suptitle(
        "Tracked cells per experiment (after QC)\nColour = replicate",
        y=1.02,
    )
    g.savefig(out_path, bbox_inches="tight")
    plt.close(g.fig)
    logger.info("Plot gespeichert: %s", out_path.name)

    # ------------------------------------------------------------------
    # Plot 2: Summary – Replikate aggregiert, STD, Ctrls separat
    # ------------------------------------------------------------------
    out_summary = out_path.parent / (out_path.stem + "_summary.pdf")

    ov2 = overview.copy()

    # Kontroll-Flag aus condition-Spalte extrahieren
    # z.B. "Osc12_NegCtrl" -> ctrl_type="NegCtrl", "Osc12" -> ctrl_type="Exp"
    def _ctrl_label(cond: str) -> str:
        m = re.search(ctrl_pattern, str(cond), re.IGNORECASE)
        return m.group(0) if m else "Experiment"

    ov2["ctrl_type"] = ov2["condition"].apply(_ctrl_label)

    # Separate DataFrames
    exp_df  = ov2[ov2["ctrl_type"] == "Experiment"]
    ctrl_df = ov2[ov2["ctrl_type"] != "Experiment"]

    # Aggregation Experimente: Mittelwert + STD ueber Replikate
    # (erst pro Replikat summieren ueber Kammern, dann ueber Replikate)
    exp_agg = (
        exp_df.groupby(["biosensor", "osc_type", "osc_freq", "replicate"])["n_tracks"]
        .sum()                          # Kammern eines Replikats zusammenzaehlen
        .reset_index()
        .groupby(["biosensor", "osc_type", "osc_freq"])["n_tracks"]
        .agg(mean="mean", std="std", n_reps="count")
        .reset_index()
    )
    exp_agg["std"] = exp_agg["std"].fillna(0)
    exp_agg["ctrl_type"] = "Experiment"

    # Aggregation Kontrollen: Mittelwert + STD getrennt nach ctrl_type
    if not ctrl_df.empty:
        ctrl_agg = (
            ctrl_df.groupby(["biosensor", "osc_type", "osc_freq", "ctrl_type", "replicate"])["n_tracks"]
            .sum()
            .reset_index()
            .groupby(["biosensor", "osc_type", "osc_freq", "ctrl_type"])["n_tracks"]
            .agg(mean="mean", std="std", n_reps="count")
            .reset_index()
        )
        ctrl_agg["std"] = ctrl_agg["std"].fillna(0)
        plot_df = pd.concat([exp_agg, ctrl_agg], ignore_index=True)
    else:
        logger.info("Keine NegCtrl/PosCtrl in 'condition' gefunden - nur Experiment-Balken.")
        plot_df = exp_agg

    # Farbe = Stamm; Kontrollart ueber die Fuellung: Oszillation voll, Feast-Kontrolle heller Ton,
    # Famine-Kontrolle hohl (weiss mit Rand in der Stammfarbe).
    ctrl_types  = ["Experiment"] + sorted(plot_df["ctrl_type"].unique().tolist())
    ctrl_types  = list(dict.fromkeys(ctrl_types))  # Reihenfolge: Exp zuerst, dedup
    n_types     = len(ctrl_types)
    ctrl_labels = {"Experiment": "oscillation chambers", "PosCtrl": CONTROL_LABELS["PosCtrl"],
                   "NegCtrl": CONTROL_LABELS["NegCtrl"]}

    def _fill(ct: str, color: str) -> dict:
        if ct == "PosCtrl":
            return dict(color=tint(color, 0.5), edgecolor=edge_color(color), linewidth=0.5)
        if ct == "NegCtrl":
            return dict(color=SURFACE, edgecolor=edge_color(color), linewidth=0.8)
        return dict(color=color, edgecolor=edge_color(color), linewidth=0.5)

    biosensors = ordered_strains(plot_df["biosensor"].dropna().unique())
    osc_types  = sorted(plot_df["osc_type"].dropna().unique())

    fig, axes = plt.subplots(
        len(osc_types), len(biosensors),
        figsize=(3.0 * len(biosensors) + 0.4, 2.6 * len(osc_types) + 0.4),
        squeeze=False,
    )

    bar_width   = 0.8 / n_types
    x_positions = {freq: i for i, freq in enumerate(order)}

    for i, osc_type in enumerate(osc_types):
        for j, biosensor in enumerate(biosensors):
            ax = axes[i][j]
            sub = plot_df[
                (plot_df["osc_type"]  == osc_type) &
                (plot_df["biosensor"] == biosensor)
            ]
            if sub.empty:
                ax.set_visible(False)
                continue

            for k, ct in enumerate(ctrl_types):
                ct_sub = sub[sub["ctrl_type"] == ct]
                if ct_sub.empty:
                    continue

                xs = [
                    x_positions[f] + (k - n_types / 2 + 0.5) * bar_width
                    for f in ct_sub["osc_freq"]
                    if f in x_positions
                ]
                ys    = ct_sub.loc[ct_sub["osc_freq"].isin(x_positions), "mean"].tolist()
                yerrs = ct_sub.loc[ct_sub["osc_freq"].isin(x_positions), "std"].tolist()

                ax.bar(xs, ys, width=bar_width * 0.9, **_fill(ct, strain_color(biosensor)))
                ax.errorbar(xs, ys, yerr=yerrs, fmt="none", color=INK_SOFT, capsize=2, linewidth=0.8)

            ax.set_xticks(range(len(order)))
            ax.set_xticklabels(order, rotation=45, ha="right")
            panel_title(ax, f"{osc_type} | {biosensor}")
            if j == 0:
                ax.set_ylabel("number of tracks")
            if i == len(osc_types) - 1:
                ax.set_xlabel("cycle period [min]")

    handles = [Patch(label=ctrl_labels.get(ct, ct), **_fill(ct, INK_SOFT)) for ct in ctrl_types]
    legend_below(fig, handles, ncol=len(handles), y=0.0)
    fig.suptitle("tracks per chamber, mean ± SD over the chambers of each structure", fontsize=9.5, y=1.0)
    finish(fig, out_summary, logger)


def _aggregate_over_replicates(
    df: pd.DataFrame, value_col: str, group_cols: list[str]
) -> pd.DataFrame:
    """
    Korrekte zweistufige Aggregation:
      1) pro Replikat (über Zellen & Frames innerhalb group_cols+replicate) mitteln
      2) über Replikate: Mittelwert + SEM
    """
    per_rep = (
        df.groupby(group_cols + ["replicate"])[value_col]
        .mean()
        .reset_index()
    )
    agg = (
        per_rep.groupby(group_cols)[value_col]
        .agg(mean="mean", sem=lambda s: s.std(ddof=1) / np.sqrt(len(s)) if len(s) > 1 else 0.0, n_reps="count")
        .reset_index()
    )
    return agg


# Kontrollen als Referenzkurven in Tintengrau: Feast gestrichelt, Famine gepunktet (plot_style).
CONTROL_REFERENCE_STYLE = {
    "NegCtrl": {"color": INK_SOFT, "linestyle": CONTROL_LINESTYLES["NegCtrl"]},
    "PosCtrl": {"color": INK_SOFT, "linestyle": CONTROL_LINESTYLES["PosCtrl"]},
}


def plot_metric_over_time_by_frequency(
    df: pd.DataFrame,
    value_col: str,
    out_path: Path,
    freq_order: Optional[Sequence[str]] = None,
    ylabel: Optional[str] = None,
    reference_cells: Optional[pd.DataFrame] = None,
    x_col: str = "osc_freq",
    facet_col: str = "osc_type",
) -> None:
    """
    Zeitverlauf einer Metrik (z.B. Sensor-Intensität, Zellfläche), facettiert
    nach Biosensor x facet_col, Farbe = x_col (Periode; statisch: Medium).
    Linie = Mittelwert über 'replicate', Band = SEM.

    WAS DAS BAND BEDEUTET, haengt vom Zweig ab und steht deshalb im Titel:
    bei den Oszillationsdaten ist 'replicate' ein Array-Index auf EINEM Chip -
    das Band ist die Streuung der Kammern eines Chips (technisch). Bei den
    statischen Daten ist jedes 'replicate' ein eigener Chip (biologisch).
    Entschieden wird das aus der Spalte 'chip', nicht aus einer Annahme.

    reference_cells : optional die PosCtrl-/NegCtrl-Zeilen (Spalte 'condition'
        genügt). Sie werden in JEDER Facette als gestrichelte Referenzkurven
        mitgezeichnet - ohne sie steht der Plot ohne Bezugsrahmen da.

        Hintergrund: die Aufrufer übergeben hier `cells_plot`, aus dem
        exclude_controls() PosCtrl/NegCtrl bereits entfernt hat. Die
        Oszillationskurven waren damit gegen NICHTS zu lesen, obwohl genau der
        Vergleich "wo liegt eine Periode relativ zu durchgehend Feast bzw.
        durchgehend Starvation" die Aussage der Abbildung ist. Die Kontrollen
        bleiben absichtlich eine SEPARATE Eingabe und keine weitere Kategorie
        auf der Frequenzachse: sie haben keine Periode.
    """
    group_cols = ["biosensor", facet_col, x_col, "time_h"]
    agg = _aggregate_over_replicates(df, value_col, group_cols)

    if "chip" in df.columns:
        chips_per_condition = df.groupby(["biosensor", facet_col, x_col], dropna=False)["chip"].nunique()
        band_unit = "chips" if chips_per_condition.max() > 1 else "chambers of ONE chip (technical)"
    else:
        band_unit = "'replicate' rows"

    ref_agg = None
    if reference_cells is not None and not reference_cells.empty:
        ref = reference_cells.copy()
        if "condition_type" not in ref.columns:
            if "condition" in ref.columns:
                from growth_rate import classify_condition_type
                ref["condition_type"] = ref["condition"].map(classify_condition_type)
            else:
                logger.warning(
                    "plot_metric_over_time_by_frequency(): reference_cells ohne 'condition'/"
                    "'condition_type' - Referenzkurven werden übersprungen."
                )
                ref = ref.iloc[0:0]
        ref = ref[ref.get("condition_type", pd.Series(dtype=str)).isin(CONTROL_REFERENCE_STYLE)]
        if not ref.empty and value_col in ref.columns:
            ref_agg = _aggregate_over_replicates(
                ref, value_col, ["biosensor", facet_col, "condition_type", "time_h"]
            ).dropna(subset=["mean"])
            if ref_agg.empty:
                ref_agg = None

    # Nur Biosensoren/Oszillationstypen mit tatsächlichen (nicht-NaN) Werten
    # für DIESE value_col anzeigen - sonst entstehen leere Panels für
    # Biosensoren, bei denen der Kanal/das Ratio schlicht nicht existiert
    # (z.B. kein RFP-Kanal bei BSA).
    has_data = agg.dropna(subset=["mean"])
    if has_data.empty:
        logger.warning("plot_metric_over_time_by_frequency(): keine gültigen Werte für '%s' - Plot übersprungen.", value_col)
        return

    biosensors = ordered_strains(has_data["biosensor"].unique())
    osc_types = sorted(has_data[facet_col].dropna().unique())
    freqs = freq_order if freq_order is not None else natural_freq_sort(agg[x_col].dropna().unique())
    # Periode = Hell-Dunkel-Rampe der Stammfarbe (kurz hell, lang dunkel); kategoriale x
    # (statisch: Medium) = Linienstil in der Stammfarbe.
    by_period = x_col == "osc_freq"
    linestyles = ["-", "--", ":", "-."]

    fig, axes = plt.subplots(
        len(osc_types), len(biosensors),
        figsize=(2.9 * len(biosensors) + 0.4, 2.6 * len(osc_types) + 0.4),
        squeeze=False, sharex=True,
    )

    for i, osc_type in enumerate(osc_types):
        for j, biosensor in enumerate(biosensors):
            ax = axes[i][j]
            ramp = strain_ramp(biosensor, len(freqs)) if by_period else [strain_color(biosensor)] * len(freqs)
            palette = dict(zip(freqs, ramp))
            styles = dict(zip(freqs, ["-"] * len(freqs) if by_period else [linestyles[k % 4] for k in range(len(freqs))]))

            # Kontrollen zuerst und im Hintergrund (zorder), damit die
            # Oszillationskurven darüber liegen und lesbar bleiben.
            if ref_agg is not None:
                rsub = ref_agg[(ref_agg[facet_col] == osc_type)
                               & (ref_agg["biosensor"] == biosensor)]
                for ctype, style in CONTROL_REFERENCE_STYLE.items():
                    rline = rsub[rsub["condition_type"] == ctype].sort_values("time_h")
                    if rline.empty:
                        continue
                    ax.plot(rline["time_h"], rline["mean"], linewidth=1.3, zorder=1, **style)
                    ax.fill_between(
                        rline["time_h"], rline["mean"] - rline["sem"], rline["mean"] + rline["sem"],
                        color=style["color"], alpha=0.08, lw=0, zorder=0,
                    )

            sub = agg[(agg[facet_col] == osc_type) & (agg["biosensor"] == biosensor)]
            for freq in freqs:
                line = sub[sub[x_col] == freq].sort_values("time_h")
                if line.empty:
                    continue
                ax.plot(line["time_h"], line["mean"], color=palette[freq], linestyle=styles[freq], linewidth=1.4)
                ax.fill_between(
                    line["time_h"], line["mean"] - line["sem"], line["mean"] + line["sem"],
                    color=palette[freq], alpha=0.12, lw=0,
                )
            panel_title(ax, f"{osc_type} | {biosensor}")
            if i == len(osc_types) - 1:
                ax.set_xlabel("time [h]")
            if j == 0:
                ax.set_ylabel(ylabel or value_col)

    if by_period:
        grey = dict(zip(freqs, strain_ramp("_", len(freqs))))
        handles = [Line2D([], [], color=grey[f], linewidth=2, label=f"{f} min") for f in freqs]
    else:
        handles = [Line2D([], [], color=INK_SOFT, linewidth=1.6, linestyle=linestyles[k % 4], label=str(f))
                   for k, f in enumerate(freqs)]
    if ref_agg is not None:
        handles += [Line2D([], [], color=INK_SOFT, linewidth=1.3, linestyle=CONTROL_LINESTYLES[ct],
                           label=CONTROL_LABELS[ct]) for ct in ("PosCtrl", "NegCtrl")]
    title = ("cycle period [min] (light = short, dark = long; hue = strain)" if by_period else x_col.replace("_", " "))
    legend_below(fig, handles, ncol=min(len(handles), 8), y=0.0, title=title)
    fig.suptitle(f"{ylabel or value_col} over time; band = ± SEM over {band_unit}", fontsize=9.5, y=1.0)
    finish(fig, out_path, logger)


def plot_morphology_scatter(
    df: pd.DataFrame, out_path: Path, facet_cols: Sequence[str] = ("biosensor",), um_per_px: Optional[float] = None,
    min_frames: int = 10, large_um2: float = 30.0, round_ecc: float = 0.6, max_points: int = 4000,
    random_state: int = 0,
) -> None:
    """Morphologie je Zelle: mittlere Flaeche (log) gegen mittlere Exzentrizitaet, eine Zelle = ein Punkt
    (Tracks mit >= min_frames Frames), eine Facette je facet_cols-Kombination (Oszillation: Stamm; statisch:
    Chip-Familie x Medium), Punkte in der Stammfarbe. Hilfslinien bei large_um2 und round_ecc teilen die Ebene
    nach den Zelltypen von Rensink et al. 2026 (klein/laenglich = hefeartige Zellen, gross/rund = geschwollene
    Zellen); der Anteil gross-runder Zellen steht im Paneltitel. Flaeche in um2, wenn um_per_px gegeben, sonst
    in px2 (dann ohne Hilfslinien). Hoechstens max_points Punkte je Facette (Stichprobe), die Anteile im Titel
    rechnen mit allen Zellen."""
    need = {"area", "eccentricity", "cell_uid", "frame"}
    if not need.issubset(df.columns):
        logger.warning("Spalten %s fehlen - Morphologie-Plot wird übersprungen.", sorted(need - set(df.columns)))
        return
    facet_cols = [c for c in facet_cols if c in df.columns]
    per_cell = (df.dropna(subset=["area", "eccentricity"])
                .groupby("cell_uid").agg(area=("area", "mean"), ecc=("eccentricity", "mean"), n=("frame", "nunique")))
    per_cell = per_cell[per_cell["n"] >= min_frames]
    if per_cell.empty:
        logger.warning("Morphologie-Plot: keine Zelle mit >= %d Frames - übersprungen.", min_frames)
        return
    meta_cols = list(dict.fromkeys(list(facet_cols) + (["biosensor"] if "biosensor" in df.columns else [])))
    per_cell = per_cell.join(df.drop_duplicates("cell_uid").set_index("cell_uid")[meta_cols])
    scale = um_per_px ** 2 if um_per_px else 1.0
    per_cell["area_u"] = per_cell["area"] * scale
    unit = "µm²" if um_per_px else "px²"

    if facet_cols:
        keys = per_cell.drop_duplicates(facet_cols)[facet_cols].astype(str).apply(tuple, axis=1).tolist()
        if facet_cols == ["biosensor"]:
            keys = [(s,) for s in ordered_strains([k[0] for k in keys])]
        else:
            keys = sorted(keys)
    else:
        keys = [()]
    n = len(keys)
    fig, axes = plt.subplots(1, n, figsize=(2.6 * n + 0.6, 3.1), squeeze=False, sharex=True, sharey=True)
    rng = np.random.default_rng(random_state)
    for ax, key in zip(axes[0], keys):
        sub = per_cell
        for col, val in zip(facet_cols, key):
            sub = sub[sub[col].astype(str) == val]
        strain = str(sub["biosensor"].iloc[0]) if "biosensor" in sub.columns and len(sub) else "_"
        color = strain_color(strain)
        n_all = len(sub)
        if um_per_px:
            large_round = float(((sub["area_u"] >= large_um2) & (sub["ecc"] < round_ecc)).mean()) if n_all else float("nan")
        show = sub.sample(max_points, random_state=rng.integers(1 << 30)) if n_all > max_points else sub
        ax.scatter(show["area_u"], show["ecc"], s=6, facecolor=color, edgecolor="none", alpha=0.35, zorder=2)
        if um_per_px:
            ax.axvline(large_um2, color=INK_MUTED, linewidth=0.8, linestyle=(0, (4, 2)), zorder=1)
            ax.axhline(round_ecc, color=INK_MUTED, linewidth=0.8, linestyle=(0, (4, 2)), zorder=1)
        ax.set_xscale("log")
        ax.set_ylim(0, 1)
        label = " ".join(key) if key else "all cells"
        title = f"{label}  (n = {n_all:,} cells)"
        if um_per_px and n_all:
            title += f"\nlarge and round: {100 * large_round:.0f} %"
        panel_title(ax, title)
        ax.set_xlabel(f"mean cell area [{unit}]")
    axes[0][0].set_ylabel("mean eccentricity\n(0 = circle, 1 = line)")
    if um_per_px:
        handles = [Line2D([], [], color=INK_MUTED, linewidth=0.8, linestyle=(0, (4, 2)),
                          label=f"guides: {large_um2:g} {unit} and eccentricity {round_ecc:g}")]
        legend_below(fig, handles, ncol=1, y=0.0)
    finish(fig, out_path, logger)


def plot_single_cell_trajectories(
    df: pd.DataFrame, value_col: str, out_path: Path, min_coverage: float = 0.8
) -> None:
    """
    Einzelzell-Trajektorien für je eine Beispielkammer pro
    Biosensor x Oszillationstyp x Frequenz - visueller QC-Check, ob die
    Tracks nach Exclusion plausibel aussehen. Nur Tracks die >= min_coverage
    der maximalen Frame-Anzahl in dieser Kammer abdecken.
    """
    example_keys = (
        df[["biosensor", "osc_type", "osc_freq", "replicate", "chamber"]]
        .drop_duplicates()
        .groupby(["biosensor", "osc_type", "osc_freq"])
        .first()
        .reset_index()
    )

    traj = df.merge(example_keys, on=["biosensor", "osc_type", "osc_freq", "replicate", "chamber"])
    if traj.empty:
        logger.warning("Keine Daten für Einzelzell-Trajektorien gefunden.")
        return

    max_frames = traj.groupby(["biosensor", "osc_type", "osc_freq"])["frame"].transform("max")
    coverage = traj.groupby("cell_uid")["frame"].transform("count")
    stable = traj[coverage > min_coverage * max_frames]

    if stable.empty:
        logger.warning("Keine stabilen Tracks (>= %.0f%% Frame-Abdeckung) gefunden.", min_coverage * 100)
        return

    osc_types = sorted(stable["osc_type"].dropna().unique())
    freqs = natural_freq_sort(stable["osc_freq"].dropna().unique())
    biosensors = sorted(stable["biosensor"].dropna().unique())

    # NUR tatsächlich vorkommende (Frequenz, Biosensor)-Kombinationen als
    # Spalten - der volle kartesische Produkt aus allen Frequenzen x allen
    # Biosensoren (wie vorher) erzeugt eine Spalte pro Kombination, auch wenn
    # ein Biosensor z.B. nur 2 von 6 Frequenzen hat -> viele leere Panels.
    freq_rank = {f: i for i, f in enumerate(freqs)}
    col_keys = sorted(
        stable[["osc_freq", "biosensor"]].drop_duplicates().itertuples(index=False, name=None),
        key=lambda fb: (freq_rank.get(fb[0], len(freqs)), fb[1]),
    )
    n_cols = len(col_keys)

    fig, axes = plt.subplots(
        len(osc_types), max(n_cols, 1),
        figsize=(3.2 * max(n_cols, 1), 3.0 * len(osc_types)),
        squeeze=False, sharey=False,
    )
    for i, osc_type in enumerate(osc_types):
        for j, (freq, biosensor) in enumerate(col_keys):
            ax = axes[i][j]
            sub = stable[
                (stable["osc_type"] == osc_type)
                & (stable["osc_freq"] == freq)
                & (stable["biosensor"] == biosensor)
            ]
            if sub.empty:
                ax.set_visible(False)
                continue
            for _, track in sub.groupby("cell_uid"):
                track = track.sort_values("time_h")
                ax.plot(track["time_h"], track[value_col], alpha=0.45, linewidth=0.6, color=strain_color(biosensor))
            panel_title(ax, f"{osc_type} | {freq} | {biosensor}")
            if i == len(osc_types) - 1:
                ax.set_xlabel("Time [h]")
            if j == 0:
                ax.set_ylabel(value_col)

    fig.suptitle(f"{value_col}: single-cell trajectories (one example chamber per condition)", y=1.02)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    logger.info("Plot gespeichert: %s", out_path.name)


def summary_statistics(df: pd.DataFrame, intensity_cols: list[str]) -> pd.DataFrame:
    """Zusammenfassungstabelle: pro Biosensor x Oszillationstyp x Frequenz x Replikat."""
    agg_dict = {"track_id": "nunique"}
    if "area" in df.columns:
        agg_dict["area"] = "mean"
    for c in intensity_cols:
        agg_dict[c] = "mean"

    # Intensity columns are already named 'mean_<channel>' (see
    # find_intensity_columns()), so prefixing them again produced
    # 'mean_mean_GFP'. Only add the prefix where it is not there yet.
    intensity_renames = {
        c: c if c.startswith("mean_") else f"mean_{c}" for c in intensity_cols
    }
    summary = (
        df.groupby(["biosensor", "osc_type", "osc_freq", "replicate"])
        .agg(agg_dict)
        .rename(columns={"track_id": "n_tracks", "area": "mean_area", **intensity_renames})
        .reset_index()
        .sort_values(["biosensor", "osc_type", "osc_freq", "replicate"])
    )
    return summary
