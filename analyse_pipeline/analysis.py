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

sns.set_theme(style="whitegrid", context="notebook")


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
        height=3.5, aspect=1.3, palette="Set2",
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

    # Farb- & Schraffur-Schema
    ctrl_types  = ["Experiment"] + sorted(plot_df["ctrl_type"].unique().tolist())
    ctrl_types  = list(dict.fromkeys(ctrl_types))  # Reihenfolge: Exp zuerst, dedup
    n_types     = len(ctrl_types)
    palette     = dict(zip(ctrl_types, sns.color_palette("Set2", n_colors=n_types)))
    hatch_map   = {ct: ("" if ct == "Experiment" else "///") for ct in ctrl_types}

    biosensors = sorted(plot_df["biosensor"].dropna().unique())
    osc_types  = sorted(plot_df["osc_type"].dropna().unique())

    fig, axes = plt.subplots(
        len(osc_types), len(biosensors),
        figsize=(4.5 * len(biosensors), 3.5 * len(osc_types)),
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

                bars = ax.bar(
                    xs, ys, width=bar_width * 0.9,
                    color=palette[ct], hatch=hatch_map[ct],
                    label=ct, edgecolor="white", linewidth=0.5,
                )
                ax.errorbar(
                    xs, ys, yerr=yerrs,
                    fmt="none", color="black", capsize=3, linewidth=1,
                )

            ax.set_xticks(range(len(order)))
            ax.set_xticklabels(order, rotation=45, ha="right")
            ax.set_title(f"{osc_type} | {biosensor}", fontsize=10)
            if j == 0:
                ax.set_ylabel("Number of tracks")
            if i == len(osc_types) - 1:
                ax.set_xlabel("Oscillation frequency")

            # Legende nur einmal pro Axes, Duplikate entfernen
            handles, labels = ax.get_legend_handles_labels()
            seen = {}
            for h, l in zip(handles, labels):
                seen.setdefault(l, h)
            ax.legend(seen.values(), seen.keys(), fontsize=7, title="Group")

    fig.suptitle(
        "Number of tracks per frequency – aggregated over replicates\n"
        "(bars = mean over replicates; error bars = ± SD; hatching = controls)",
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(out_summary, bbox_inches="tight")
    plt.close(fig)
    logger.info("Plot gespeichert: %s", out_summary.name)


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


CONTROL_REFERENCE_STYLE = {
    "NegCtrl": {"color": "#4C78A8", "linestyle": (0, (5, 2))},
    "PosCtrl": {"color": "#E45756", "linestyle": (0, (3, 1, 1, 1))},
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

    biosensors = sorted(has_data["biosensor"].unique())
    osc_types = sorted(has_data[facet_col].dropna().unique())
    freqs = freq_order if freq_order is not None else natural_freq_sort(agg[x_col].dropna().unique())
    palette = dict(zip(freqs, sns.color_palette("viridis", n_colors=len(freqs))))

    fig, axes = plt.subplots(
        len(osc_types), len(biosensors),
        figsize=(4.5 * len(biosensors), 3.5 * len(osc_types)),
        squeeze=False, sharex=True,
    )

    for i, osc_type in enumerate(osc_types):
        for j, biosensor in enumerate(biosensors):
            ax = axes[i][j]

            # Kontrollen zuerst und im Hintergrund (zorder), damit die
            # Oszillationskurven darüber liegen und lesbar bleiben.
            if ref_agg is not None:
                rsub = ref_agg[(ref_agg[facet_col] == osc_type)
                               & (ref_agg["biosensor"] == biosensor)]
                for ctype, style in CONTROL_REFERENCE_STYLE.items():
                    rline = rsub[rsub["condition_type"] == ctype].sort_values("time_h")
                    if rline.empty:
                        continue
                    ax.plot(rline["time_h"], rline["mean"], label=ctype, linewidth=1.9,
                            zorder=1, **style)
                    ax.fill_between(
                        rline["time_h"], rline["mean"] - rline["sem"], rline["mean"] + rline["sem"],
                        color=style["color"], alpha=0.10, zorder=0,
                    )

            sub = agg[(agg[facet_col] == osc_type) & (agg["biosensor"] == biosensor)]
            for freq in freqs:
                line = sub[sub[x_col] == freq].sort_values("time_h")
                if line.empty:
                    continue
                ax.plot(line["time_h"], line["mean"], label=str(freq), color=palette[freq], linewidth=1.6)
                ax.fill_between(
                    line["time_h"], line["mean"] - line["sem"], line["mean"] + line["sem"],
                    color=palette[freq], alpha=0.15,
                )
            ax.set_title(f"{osc_type} | {biosensor}", fontsize=10)
            if i == len(osc_types) - 1:
                ax.set_xlabel("Time [h]")
            if j == 0:
                ax.set_ylabel(ylabel or value_col)

    handles, labels = axes[0][0].get_legend_handles_labels()
    # Perioden zuerst, Kontrollen hinten - die Kontrollen sind der Bezugsrahmen,
    # nicht eine weitere Stufe der Dosis.
    order = ([k for k in range(len(labels)) if labels[k] not in CONTROL_REFERENCE_STYLE]
             + [k for k in range(len(labels)) if labels[k] in CONTROL_REFERENCE_STYLE])
    handles = [handles[k] for k in order]
    labels = [labels[k] for k in order]
    legend_title = ("Cycle period [min]" if x_col == "osc_freq" else x_col) + \
        ("  /  controls" if ref_agg is not None else "")
    # -0.12 statt -0.05: mit den Kontrollen sind es bis zu 8 Legendeneinträge,
    # und bei -0.05 lag die Legende auf den "Time [h]"-Achsenbeschriftungen.
    if handles:
        fig.legend(handles, labels, title=legend_title, loc="lower center",
                   ncol=max(1, min(len(labels), 8)), bbox_to_anchor=(0.5, -0.12))
    subtitle = f"(line = mean; band = ± SEM over {band_unit}"
    subtitle += "; dashed = constant-medium controls)" if ref_agg is not None else ")"
    by = "feast/famine cycle period" if x_col == "osc_freq" else x_col
    fig.suptitle(f"{value_col} over time, by {by}\n{subtitle}", y=1.02)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    logger.info("Plot gespeichert: %s", out_path.name)


def plot_morphology_scatter(
    df: pd.DataFrame, out_path: Path, sample_n: int = 1000, random_state: int = 0
) -> None:
    """Scatter Fläche vs. Exzentrizität, facettiert nach Biosensor x Oszillationstyp,
    Farbe = Oszillationsfrequenz. Stichprobe pro Gruppe, damit der Plot nicht zu schwer wird."""
    if not {"area", "eccentricity"}.issubset(df.columns):
        logger.warning("Spalten 'area'/'eccentricity' fehlen - Morphologie-Plot wird übersprungen.")
        return

    sample = (
        df.groupby(["biosensor", "osc_type", "osc_freq"], group_keys=False)[df.columns]
        .apply(lambda g: g.sample(min(len(g), sample_n), random_state=random_state))
        .reset_index(drop=True)
    )

    g = sns.relplot(
        data=sample, x="area", y="eccentricity", hue="osc_freq",
        hue_order=natural_freq_sort(sample["osc_freq"].dropna().unique()),
        row="osc_type", col="biosensor",
        kind="scatter", alpha=0.3, s=12, palette="viridis",
        height=3.5, aspect=1.3,
    )
    for ax in g.axes.flat:
        ax.set_xscale("log")
    g.set_axis_labels("Cell area [px²] (log)", "Eccentricity (0 = circle, 1 = line)")
    g.set_titles("{row_name} | {col_name}")
    g.fig.suptitle("Morphology: area vs. eccentricity (sample after QC)", y=1.02)
    g.savefig(out_path, bbox_inches="tight")
    plt.close(g.fig)
    logger.info("Plot gespeichert: %s", out_path.name)


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
                ax.plot(track["time_h"], track[value_col], alpha=0.5, linewidth=0.6)
            ax.set_title(f"{osc_type} | {freq} | {biosensor}", fontsize=8)
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
