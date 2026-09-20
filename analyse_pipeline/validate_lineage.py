#!/usr/bin/env python3
"""
validate_lineage.py
===================
QUANTITATIVE Validierung der Mutter/Bud-Heuristik aus lineage.py - das
Gegenstück zu plot_qc_lineage_overlay.py, das einzelne Kammern per Auge prüft.

WARUM
-----
classify_mother_bud() entscheidet über Budding Ratio, µ_event, den Stammbaum
und die Mutter/Knospe-Trennung in area_growth.py. Die Funktion loggt zwar, wie
viele Bud-Kandidaten KEINER Mutter zugeordnet werden konnten - aber nur als
EINE Gesamtzahl über alle Kammern. Damit bleibt die entscheidende Frage offen:

    Unterscheidet sich die Erkennungsrate systematisch zwischen den
    Bedingungen?

Wenn ja, ist jeder Vergleich zwischen Oszillationsfrequenzen teilweise ein
Vergleich von Trackingqualität statt von Biologie. Dieses Skript macht das
messbar.

WAS ES PRÜFT
------------
1. ERKENNUNGSRATE pro Kammer, plus Kruskal-Wallis-Test über die
   Frequenz-Batches (analog control_consistency.py): ist die Rate über die
   Bedingungen hinweg konstant?
2. d/r-VERTEILUNG (Distanz / adaptiver Suchradius). Sitzen die Zuordnungen
   komfortabel innerhalb des Radius (klein), oder holt der Radius sie gerade
   noch herein (nahe 1.0)? Letzteres heißt: der Toleranzwert macht die
   Arbeit, nicht die tatsächliche Nachbarschaft.
3. KONKURRENZ: wie viele ANDERE etablierte Zellen hätten denselben Bud
   ebenfalls beanspruchen können? classify_mother_bud() nimmt greedy die
   nächstgelegene. Bei n_competing_mothers > 1 ist die Zuordnung eine
   Entscheidung, keine Beobachtung.
4. TOLERANZ-KURVE: welcher Anteil der Kandidaten würde bei welchem
   tolerance_px zugeordnet? Zeigt, ob der aktuelle Wert auf einem Plateau
   liegt (robust) oder auf einer steilen Flanke (jede kleine Änderung
   verschiebt die Ergebnisse).
5. GROESSENKRITERIUM (bud_size.py): angespuelte Blastokonidien tauchen "neu"
   neben sitzenden Zellen auf und bestehen die raeumliche Zuordnung wie eine
   Knospe - sind aber beim ersten Auftreten etwa so gross wie die Mutter.
   Kandidaten ueber der Schwelle (bud_area / mother_area) zaehlen hier NICHT
   in den Nenner der Erkennungsrate; sie stehen pro Kammer in
   n_rejected_by_size. Die Schwelle ist dieselbe, mit der run_analysis.py
   die Events erzeugt hat (20_bud_size_threshold.csv), damit die Validierung
   die Heuristik prueft, die tatsaechlich gelaufen ist.

AUSFÜHREN
---------
    python validate_lineage.py \\
        --cells   ../analysis_output/00_cell_positions.parquet \\
        --lineage-events ../analysis_output/20_budding_events.csv \\
        --out-dir ../analysis_output/lineage_validation

Ohne Argumente werden die Pfade aus config.py verwendet.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Optional

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import seaborn as sns  # noqa: E402
from scipy import stats  # noqa: E402

from analysis import natural_freq_sort, pretty_label  # noqa: E402
from lineage import LineageParams  # noqa: E402
from bud_size import load_bud_size_threshold, resolve_bud_size_threshold  # noqa: E402

logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

sns.set_theme(style="whitegrid", context="notebook")

META_COLS = ["biosensor", "osc_type", "osc_freq", "condition", "replicate", "chamber"]


# ==============================================================================
# 1. Diagnose pro Bud-Kandidat
# ==============================================================================

def compute_candidate_diagnostics(
    cells: pd.DataFrame,
    lineage_events: pd.DataFrame,
    params: LineageParams,
) -> pd.DataFrame:
    """
    Eine Zeile pro Bud-KANDIDAT (zugeordnet oder nicht), mit allem, was zur
    Beurteilung der Zuordnung nötig ist.

    Die Kandidaten-Definition und die Etabliertheits-Prüfung spiegeln
    lineage.classify_mother_bud() exakt wider: Kandidat = Track, der nach dem
    ersten Frame der Kammer neu auftaucht; Mutter-Kandidat = Zelle, die zum
    Zeitpunkt des Auftauchens bereits >= established_min_frames Frames
    gesehen wurde. Wird die Heuristik dort geändert, muss diese Funktion
    mitgezogen werden - sie ist bewusst eigenständig, damit die Validierung
    nicht einfach die Implementierung gegen sich selbst prüft.

    Groessenkriterium (lineage.py, bud_size.py): pro Kandidat wird zusaetzlich
    bud_area_fraction = Flaeche beim ersten Auftreten / Flaeche der
    Referenzmutter berechnet - der naechsten etablierten Zelle im adaptiven
    Radius (exakt die Zuordnung der Heuristik), ersatzweise der naechsten
    etablierten Zelle ueberhaupt. Ob ein Kandidat damit als Knospe in Frage
    kommt, entscheidet apply_size_eligibility() mit derselben Schwelle, die
    run_analysis.py verwendet hat.

    Returns
    -------
    exp_id, bud_cell_uid, budding_frame, assigned, mother_cell_uid,
    distance_px, adaptive_radius_px, d_over_r, n_competing_mothers,
    nearest_distance_px, nearest_radius_px, required_tolerance_px,
    bud_area, reference_mother_area, bud_area_fraction (+ Metadaten)

    required_tolerance_px ist der tolerance_px-Wert, der MINDESTENS nötig
    gewesen wäre, damit dieser Kandidat der nächstgelegenen etablierten Zelle
    zugeordnet worden wäre (negativ = lag ohnehin komfortabel im Radius).
    """
    required = {"exp_id", "cell_uid", "frame", "centroid_x", "centroid_y", "area"}
    missing = required - set(cells.columns)
    if missing:
        raise ValueError(f"compute_candidate_diagnostics() fehlen Spalten in 'cells': {missing}")

    assigned_lookup: dict[tuple[str, str], pd.Series] = {}
    if not lineage_events.empty and "bud_cell_uid" in lineage_events.columns:
        for _, ev in lineage_events.iterrows():
            assigned_lookup[(str(ev["exp_id"]), str(ev["bud_cell_uid"]))] = ev

    records = []
    for exp_id, group in cells.groupby("exp_id"):
        group = group.sort_values("frame")
        chamber_start = group["frame"].min()

        frames_by_track = group.groupby("cell_uid")["frame"].apply(lambda s: np.sort(s.unique()))
        first_appearance = group.groupby("cell_uid")["frame"].min()
        candidates = first_appearance[first_appearance > chamber_start]
        if candidates.empty:
            continue

        # Pro Frame einmal gruppieren statt pro Kandidat den ganzen
        # Kammer-DataFrame zu scannen.
        by_frame = {f: sub for f, sub in group.groupby("frame")}

        for bud_uid, bud_frame in candidates.items():
            bud_rows = by_frame.get(bud_frame)
            if bud_rows is None:
                continue
            bud_row = bud_rows[bud_rows["cell_uid"] == bud_uid]
            if bud_row.empty:
                continue
            bud_x = float(bud_row["centroid_x"].iloc[0])
            bud_y = float(bud_row["centroid_y"].iloc[0])
            bud_area = float(bud_row["area"].iloc[0])

            others = bud_rows[bud_rows["cell_uid"] != bud_uid]
            record = {
                "exp_id": exp_id,
                "bud_cell_uid": bud_uid,
                "budding_frame": int(bud_frame),
                "n_competing_mothers": 0,
                "nearest_distance_px": np.nan,
                "nearest_radius_px": np.nan,
                "required_tolerance_px": np.nan,
                "bud_area": bud_area,
                "reference_mother_area": np.nan,
                "bud_area_fraction": np.nan,
            }

            if not others.empty:
                est_counts = others["cell_uid"].map(
                    lambda uid: int(np.searchsorted(frames_by_track[uid], bud_frame))
                ).to_numpy()
                established = others[est_counts >= params.established_min_frames]

                if not established.empty:
                    dists = np.hypot(
                        established["centroid_x"].to_numpy(dtype=float) - bud_x,
                        established["centroid_y"].to_numpy(dtype=float) - bud_y,
                    )
                    cell_radii = np.sqrt(established["area"].to_numpy(dtype=float) / np.pi)
                    radii = cell_radii + params.tolerance_px

                    within = dists <= radii
                    record["n_competing_mothers"] = int(np.nansum(within))

                    nearest = int(np.nanargmin(dists))
                    record["nearest_distance_px"] = float(dists[nearest])
                    record["nearest_radius_px"] = float(radii[nearest])
                    # Welcher Toleranzwert haette diesen Kandidaten gerettet?
                    record["required_tolerance_px"] = float(dists[nearest] - cell_radii[nearest])

                    # Referenzmutter fuer das Groessenkriterium: die naechste
                    # Zelle IM Radius (= Wahl der Heuristik), sonst die
                    # naechste etablierte Zelle ueberhaupt.
                    if within.any():
                        reference = int(np.nanargmin(np.where(within, dists, np.inf)))
                    else:
                        reference = nearest
                    mother_area = float(established["area"].to_numpy(dtype=float)[reference])
                    record["reference_mother_area"] = mother_area
                    record["bud_area_fraction"] = bud_area / mother_area if mother_area > 0 else np.nan

            ev = assigned_lookup.get((str(exp_id), str(bud_uid)))
            record["assigned"] = ev is not None
            if ev is not None:
                record["mother_cell_uid"] = ev.get("mother_cell_uid")
                record["distance_px"] = ev.get("distance_px")
                record["adaptive_radius_px"] = ev.get("adaptive_radius_px")
            else:
                record["mother_cell_uid"] = pd.NA
                record["distance_px"] = np.nan
                record["adaptive_radius_px"] = np.nan

            records.append(record)

    out = pd.DataFrame(records)
    if out.empty:
        logger.warning("Keine Bud-Kandidaten gefunden - Validierung ergibt eine leere Tabelle.")
        return out

    with np.errstate(divide="ignore", invalid="ignore"):
        out["d_over_r"] = out["distance_px"] / out["adaptive_radius_px"]

    meta_cols = [c for c in META_COLS if c in cells.columns]
    if meta_cols:
        meta = cells[["exp_id"] + meta_cols].drop_duplicates("exp_id")
        out = out.merge(meta, on="exp_id", how="left")

    n_assigned = int(out["assigned"].sum())
    logger.info(
        "Kandidaten-Diagnose: %d Bud-Kandidaten, %d zugeordnet (%.1f%%), "
        "%d mit mehr als einer moeglichen Mutter (mehrdeutig).",
        len(out), n_assigned, 100 * n_assigned / len(out),
        int((out["n_competing_mothers"] > 1).sum()),
    )
    return out


def apply_size_eligibility(diagnostics: pd.DataFrame, bud_size_threshold: Optional[float]) -> pd.DataFrame:
    """Spalten size_eligible / bud_size_threshold: darf der Kandidat nach dem
    Groessenkriterium ueberhaupt eine Knospe sein?

    Kandidaten ohne berechenbares Verhaeltnis (keine etablierte Zelle in der
    Kammer) bleiben 'eligible' - sie fallen ohnehin an der raeumlichen
    Zuordnung, nicht an der Groesse. None oder eine unendliche Schwelle
    schaltet das Kriterium ab (alle eligible).
    """
    out = diagnostics.copy()
    if out.empty:
        return out
    if bud_size_threshold is not None and np.isfinite(bud_size_threshold):
        threshold = float(bud_size_threshold)
        out["size_eligible"] = ~(out["bud_area_fraction"] > threshold)
        out["bud_size_threshold"] = threshold
        n_rejected = int((~out["size_eligible"]).sum())
        logger.info(
            "Groessenkriterium (Schwelle %.2f x Mutterflaeche): %d von %d Kandidaten sind zu "
            "gross fuer eine Knospe und zaehlen nicht in die Erkennungsrate.",
            threshold, n_rejected, len(out),
        )
        inconsistent = int((out["assigned"].astype(bool) & ~out["size_eligible"]).sum())
        if inconsistent:
            logger.warning(
                "%d ZUGEORDNETE Events liegen UEBER der Groessenschwelle. Die Budding-Events "
                "stammen vermutlich aus einem Lauf ohne Groessenkriterium oder mit einer "
                "anderen Schwelle - run_analysis.py neu laufen lassen, sonst vergleicht die "
                "Validierung zwei verschiedene Heuristiken.", inconsistent,
            )
    else:
        out["size_eligible"] = True
        out["bud_size_threshold"] = np.nan
    return out


def summarise_per_chamber(diagnostics: pd.DataFrame) -> pd.DataFrame:
    """Erkennungsrate pro Kammer - die Zahl, die zwischen Bedingungen
    vergleichbar sein MUSS, damit Bedingungsvergleiche Biologie messen.

    Nenner sind nur die Kandidaten, die das Groessenkriterium bestehen
    (size_eligible); die zu grossen (angespuelte Zellen) stehen separat in
    n_rejected_by_size, damit sichtbar bleibt, wie viel das Kriterium pro
    Kammer wegnimmt. Ohne die Spalte size_eligible zaehlen alle Kandidaten.
    """
    if diagnostics.empty:
        return pd.DataFrame()

    meta_cols = [c for c in META_COLS if c in diagnostics.columns]
    d = diagnostics.copy()
    eligible = (d["size_eligible"].astype(bool) if "size_eligible" in d.columns
                else pd.Series(True, index=d.index))
    d["_eligible"] = eligible
    d["_assigned_eligible"] = d["assigned"].astype(bool) & eligible
    d["_competing_eligible"] = d["n_competing_mothers"].where(eligible)
    summary = (
        d.groupby(["exp_id"] + meta_cols, dropna=False)
        .agg(
            n_candidates_all=("assigned", "size"),
            n_rejected_by_size=("_eligible", lambda s: int((~s).sum())),
            n_candidates=("_eligible", "sum"),
            n_assigned=("_assigned_eligible", "sum"),
            median_d_over_r=("d_over_r", "median"),
            frac_ambiguous=("_competing_eligible",
                            lambda s: float((s.dropna() > 1).mean()) if s.notna().any() else np.nan),
        )
        .reset_index()
    )
    summary["n_candidates"] = summary["n_candidates"].astype(int)
    summary["n_assigned"] = summary["n_assigned"].astype(int)
    # Kammern, in denen kein Kandidat klein genug war, haben keine Rate (NaN),
    # nicht 0 - sonst zoegen sie den Batch-Test nach unten.
    summary["assignment_rate"] = summary["n_assigned"] / summary["n_candidates"].where(summary["n_candidates"] > 0)
    return summary


def test_detection_rate_across_conditions(
    per_chamber: pd.DataFrame,
    group_cols: Optional[list[str]] = None,
    freq_col: str = "osc_freq",
) -> pd.DataFrame:
    """
    Kruskal-Wallis auf die Erkennungsrate über die Frequenz-Batches, pro
    (biosensor, osc_type) - dasselbe Muster wie
    control_consistency.test_control_consistency_across_freq(), aber auf der
    Erkennungsrate statt auf einer biologischen Groesse.

    p < 0.05 heisst: die Heuristik erkennt in manchen Bedingungen
    systematisch mehr Budding-Events als in anderen. Jeder anschliessende
    Vergleich der Budding Ratio zwischen diesen Bedingungen ist dann
    konfundiert.
    """
    if per_chamber.empty or freq_col not in per_chamber.columns:
        return pd.DataFrame()

    if group_cols is None:
        group_cols = [c for c in ["biosensor", "osc_type"] if c in per_chamber.columns]
    if not group_cols:
        per_chamber = per_chamber.assign(_all="all")
        group_cols = ["_all"]

    results = []
    for keys, grp in per_chamber.groupby(group_cols, dropna=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        rates = grp.dropna(subset=["assignment_rate"])
        batches = [g["assignment_rate"].to_numpy() for _, g in rates.groupby(freq_col) if len(g) >= 2]
        record = dict(zip(group_cols, keys))
        record["n_batches"] = len(batches)
        record["n_chambers"] = int(len(rates))
        pooled = np.concatenate(batches) if batches else np.array([])
        if len(batches) < 2:
            record["h_statistic"], record["p_value"] = np.nan, np.nan
            record["note"] = "weniger als 2 Batches mit >= 2 Kammern"
        elif pooled.size and np.allclose(pooled, pooled[0]):
            # Alle Raten identisch (z.B. ueberall 100% erkannt): Kruskal-Wallis
            # dividiert dann durch eine Ties-Korrektur von 0 und gibt NaN mit
            # RuntimeWarning zurueck. Das ist kein Testergebnis, sondern der
            # bestmoegliche Fall - entsprechend benennen statt NaN zu melden.
            record["h_statistic"], record["p_value"] = np.nan, np.nan
            record["note"] = f"alle Erkennungsraten identisch ({pooled[0]:.3f}) - kein Test noetig"
        else:
            h, p = stats.kruskal(*batches)
            record["h_statistic"], record["p_value"] = h, p
            record["note"] = ""
        results.append(record)

    out = pd.DataFrame(results)
    if not out.empty and out["p_value"].notna().any():
        n_sig = int((out["p_value"] < 0.05).sum())
        if n_sig:
            logger.warning(
                "%d von %d Gruppen zeigen eine SIGNIFIKANT unterschiedliche Erkennungsrate "
                "zwischen den Frequenz-Batches (p < 0.05). Bedingungsvergleiche der Budding "
                "Ratio / µ_event sind dort mit der Trackingqualitaet konfundiert.",
                n_sig, len(out),
            )
        else:
            logger.info(
                "Erkennungsrate ist ueber die Frequenz-Batches hinweg konsistent "
                "(kein p < 0.05 in %d Gruppen).", len(out),
            )
    return out


def tolerance_sweep(
    diagnostics: pd.DataFrame,
    tolerances: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """
    Anteil zuordenbarer Kandidaten als Funktion von tolerance_px.

    Nutzt required_tolerance_px: ein Kandidat waere bei Toleranz t zugeordnet
    worden, wenn required_tolerance_px <= t. Das ist eine Untergrenze (es
    zaehlt nur die NAECHSTGELEGENE etablierte Zelle), reicht aber, um zu
    sehen, ob der gewaehlte Wert auf einem Plateau oder auf einer Flanke sitzt.
    """
    if diagnostics.empty or "required_tolerance_px" not in diagnostics.columns:
        return pd.DataFrame()
    if "size_eligible" in diagnostics.columns:
        # Zu grosse Kandidaten sind keine Knospen - sie sollen auch nicht als
        # "mit mehr Toleranz zuordenbar" in den Sweep eingehen.
        diagnostics = diagnostics[diagnostics["size_eligible"].astype(bool)]
        if diagnostics.empty:
            return pd.DataFrame()
    if tolerances is None:
        tolerances = np.arange(0, 101, 2.5)

    req = diagnostics["required_tolerance_px"]
    n_total = len(diagnostics)
    rows = [
        {
            "tolerance_px": float(t),
            "n_assignable": int((req <= t).sum()),
            "frac_assignable": float((req <= t).sum() / n_total),
        }
        for t in tolerances
    ]
    return pd.DataFrame(rows)


# ==============================================================================
# 2. Plots
# ==============================================================================

def plot_d_over_r_distribution(diagnostics: pd.DataFrame, out_path: Path) -> None:
    """Verteilung von Distanz/Suchradius, pro Oszillationstyp."""
    df = diagnostics.dropna(subset=["d_over_r"])
    if df.empty:
        logger.warning("plot_d_over_r_distribution(): keine zugeordneten Events - Plot uebersprungen.")
        return

    facet_col = "osc_type" if "osc_type" in df.columns else None
    facets = sorted(df[facet_col].dropna().unique()) if facet_col else [None]

    fig, axes = plt.subplots(1, len(facets), figsize=(5 * len(facets), 4.2), squeeze=False, sharey=True)
    axes = axes[0]
    for ax, facet in zip(axes, facets):
        sub = df[df[facet_col] == facet] if facet_col else df
        ax.hist(sub["d_over_r"], bins=40, range=(0, 1), color="#4C78A8", edgecolor="white", linewidth=0.4)
        ax.axvline(1.0, color="#E45756", linestyle="--", linewidth=1.2)
        ax.set_xlabel("Distance / search radius (d/r)")
        if ax is axes[0]:
            ax.set_ylabel("Budding events")
        ax.set_title(str(facet) if facet else "All data", fontsize=11, fontweight="bold")

    fig.suptitle(
        "How comfortably did each bud fall inside its mother's search radius?\n"
        "Mass near 1.0 means the tolerance, not real adjacency, made the assignment",
        y=1.04,
    )
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    logger.info("Plot gespeichert: %s", out_path.name)


def plot_assignment_rate(per_chamber: pd.DataFrame, out_path: Path,
                         freq_order: Optional[list[str]] = None) -> None:
    """Erkennungsrate pro Kammer ueber die Bedingungen - ein Punkt pro Kammer."""
    if per_chamber.empty or "osc_freq" not in per_chamber.columns:
        logger.warning("plot_assignment_rate(): keine osc_freq-Spalte - Plot uebersprungen.")
        return

    # Werte, die nicht in freq_order stehen (z.B. 'static_omlp'), hinten
    # anhaengen statt sie aus dem Plot fallen zu lassen - siehe
    # run_analysis.resolve_x_order() fuer dasselbe Problem.
    present = sorted(per_chamber["osc_freq"].dropna().unique().tolist(), key=str)
    if freq_order:
        order = [v for v in freq_order if v in present] + [v for v in present if v not in freq_order]
    else:
        order = natural_freq_sort(present)
    facet_col = "osc_type" if "osc_type" in per_chamber.columns else None
    facets = sorted(per_chamber[facet_col].dropna().unique()) if facet_col else [None]

    fig, axes = plt.subplots(1, len(facets), figsize=(5.2 * len(facets), 4.4), squeeze=False, sharey=True)
    axes = axes[0]
    rng = np.random.default_rng(0)
    for ax, facet in zip(axes, facets):
        sub = per_chamber[per_chamber[facet_col] == facet] if facet_col else per_chamber
        for x, freq in enumerate(order):
            vals = sub.loc[sub["osc_freq"] == freq, "assignment_rate"]
            if vals.empty:
                continue
            ax.scatter(np.full(len(vals), x) + rng.uniform(-.12, .12, len(vals)), vals,
                       s=30, alpha=.75, color="#4C78A8", edgecolor="white", linewidth=.4, zorder=3)
            ax.plot([x - .2, x + .2], [vals.median()] * 2, color="black", linewidth=1.8, zorder=4)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order, rotation=30)
        ax.set_ylim(0, 1.02)
        ax.set_xlabel(pretty_label("osc_freq"))
        if ax is axes[0]:
            ax.set_ylabel("Assigned bud candidates (fraction)")
        ax.set_title(str(facet) if facet else "All data", fontsize=11, fontweight="bold")

    fig.suptitle(
        "Detection rate of the mother/bud heuristic, per chamber\n"
        "A trend across frequencies would confound every condition comparison downstream",
        y=1.04,
    )
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    logger.info("Plot gespeichert: %s", out_path.name)


def plot_tolerance_sweep(sweep: pd.DataFrame, current_tolerance: float, out_path: Path) -> None:
    """Zuordenbare Kandidaten als Funktion von tolerance_px, mit aktuellem Wert."""
    if sweep.empty:
        return
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.plot(sweep["tolerance_px"], sweep["frac_assignable"], color="#4C78A8", linewidth=2)
    ax.axvline(current_tolerance, color="#E45756", linestyle="--", linewidth=1.4)
    ax.text(current_tolerance, 0.02, f" current: {current_tolerance:g} px",
            color="#E45756", fontsize=9, ha="left")
    ax.set_xlabel("tolerance_px")
    ax.set_ylabel("Assignable bud candidates (fraction)")
    ax.set_ylim(0, 1.02)
    ax.set_title(
        "Sensitivity of detection to the tolerance parameter\n"
        "A steep slope at the current value means results move with the setting",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    logger.info("Plot gespeichert: %s", out_path.name)


def plot_ambiguity(diagnostics: pd.DataFrame, out_path: Path) -> None:
    """Wie viele Muetter haetten denselben Bud beanspruchen koennen?"""
    df = diagnostics[diagnostics["assigned"]]
    if df.empty:
        return
    counts = df["n_competing_mothers"].value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    colors = ["#4C78A8" if k <= 1 else "#E45756" for k in counts.index]
    ax.bar(counts.index.astype(str), counts.values, color=colors, edgecolor="white", linewidth=.5)
    ax.set_xlabel("Established cells whose search radius also covered this bud")
    ax.set_ylabel("Budding events")
    ax.set_title(
        "Assignment ambiguity\n"
        "Red = more than one candidate mother; the nearest one was picked greedily",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    logger.info("Plot gespeichert: %s", out_path.name)


# ==============================================================================
# 3. CLI
# ==============================================================================

def _read_table(path: Path) -> pd.DataFrame:
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        logger.warning("'%s' ist leer - wird als Tabelle ohne Zeilen behandelt.", path)
        return pd.DataFrame()


def run_validation(
    cells: pd.DataFrame,
    lineage_events: pd.DataFrame,
    params: LineageParams,
    out_dir: Path,
    freq_order: Optional[list[str]] = None,
    bud_size_threshold: Optional[float] = None,
) -> pd.DataFrame:
    """bud_size_threshold: die Schwelle des Groessenkriteriums, mit der die
    Budding-Events erzeugt wurden (20_bud_size_threshold.csv). None = aus den
    Verhaeltnissen der Spiegel-Implementierung hier ableiten (dieselbe
    Statistik wie bud_size.py, aber nicht dieselben Zahlen)."""
    out_dir.mkdir(parents=True, exist_ok=True)

    diagnostics = compute_candidate_diagnostics(cells, lineage_events, params)
    if diagnostics.empty:
        logger.warning("Keine Kandidaten - Validierung abgebrochen.")
        return diagnostics
    if bud_size_threshold is None:
        from config import BUD_MAX_AREA_FRACTION_FALLBACK, BUD_SIZE_PLAUSIBLE_RANGE
        resolved = resolve_bud_size_threshold(
            diagnostics["bud_area_fraction"],
            fallback=BUD_MAX_AREA_FRACTION_FALLBACK, plausible=BUD_SIZE_PLAUSIBLE_RANGE,
        )
        bud_size_threshold = resolved.threshold
        logger.info(
            "Keine Groessenschwelle uebergeben - aus den eigenen Kandidaten abgeleitet: "
            "%.2f (%s, n = %d).", resolved.threshold, resolved.source, resolved.n_candidates,
        )
    diagnostics = apply_size_eligibility(diagnostics, bud_size_threshold)
    diagnostics.to_csv(out_dir / "lv_01_candidate_diagnostics.csv", index=False)

    per_chamber = summarise_per_chamber(diagnostics)
    per_chamber.to_csv(out_dir / "lv_02_detection_rate_per_chamber.csv", index=False)

    kruskal = test_detection_rate_across_conditions(per_chamber)
    if not kruskal.empty:
        kruskal.to_csv(out_dir / "lv_03_detection_rate_kruskal.csv", index=False)

    sweep = tolerance_sweep(diagnostics)
    if not sweep.empty:
        sweep.to_csv(out_dir / "lv_04_tolerance_sweep.csv", index=False)

    plot_d_over_r_distribution(diagnostics, out_dir / "lv_01_d_over_r_distribution.pdf")
    plot_assignment_rate(per_chamber, out_dir / "lv_02_detection_rate.pdf", freq_order=freq_order)
    plot_ambiguity(diagnostics, out_dir / "lv_03_assignment_ambiguity.pdf")
    plot_tolerance_sweep(sweep, params.tolerance_px, out_dir / "lv_04_tolerance_sweep.pdf")

    _log_verdict(diagnostics, per_chamber, kruskal, params)
    return diagnostics


def _log_verdict(
    diagnostics: pd.DataFrame,
    per_chamber: pd.DataFrame,
    kruskal: pd.DataFrame,
    params: LineageParams,
) -> None:
    """Kurzfassung als Log - damit man nicht erst vier CSVs lesen muss."""
    assigned = diagnostics["assigned"]
    d_over_r = diagnostics["d_over_r"].dropna()

    logger.info("=== Validierung der Lineage-Heuristik: Kurzfassung ===")
    logger.info("  Bud-Kandidaten gesamt:        %d", len(diagnostics))
    if "size_eligible" in diagnostics.columns:
        eligible = diagnostics["size_eligible"].astype(bool)
        thr = diagnostics["bud_size_threshold"].dropna()
        logger.info(
            "  zu gross fuer eine Knospe:    %d (%.1f%%; Schwelle %s x Mutterflaeche)",
            int((~eligible).sum()), 100 * float((~eligible).mean()),
            f"{thr.iloc[0]:.2f}" if not thr.empty else "keine",
        )
        logger.info("  Kandidaten nach Groesse:      %d", int(eligible.sum()))
        if eligible.any():
            logger.info("  davon zugeordnet:             %d (%.1f%%)",
                        int(assigned[eligible].sum()), 100 * float(assigned[eligible].mean()))
    else:
        logger.info("  davon zugeordnet:             %d (%.1f%%)", int(assigned.sum()), 100 * assigned.mean())
    if not d_over_r.empty:
        logger.info("  d/r Median:                   %.2f", d_over_r.median())
        frac_marginal = float((d_over_r > 0.8).mean())
        logger.info("  d/r > 0.8 (grenzwertig):      %.1f%%", 100 * frac_marginal)
        if frac_marginal > 0.25:
            logger.warning(
                "  -> Mehr als ein Viertel der Zuordnungen liegt im aeusseren Fuenftel des "
                "Suchradius. Das spricht dafuer, dass tolerance_px=%.1f die Zuordnungen "
                "erzeugt, statt echte Nachbarschaft zu bestaetigen.", params.tolerance_px,
            )
    frac_ambiguous = float((diagnostics.loc[assigned, "n_competing_mothers"] > 1).mean()) if assigned.any() else 0.0
    logger.info("  mehrdeutige Zuordnungen:      %.1f%%", 100 * frac_ambiguous)
    if frac_ambiguous > 0.2:
        logger.warning(
            "  -> Bei mehr als einem Fuenftel der Events kamen mehrere Muetter in Frage. "
            "Die greedy Nearest-Neighbour-Wahl in classify_mother_bud() ist dort eine "
            "Entscheidung, keine Messung."
        )
    if not per_chamber.empty:
        logger.info(
            "  Erkennungsrate pro Kammer:    Median %.2f (min %.2f, max %.2f ueber %d Kammern)",
            per_chamber["assignment_rate"].median(), per_chamber["assignment_rate"].min(),
            per_chamber["assignment_rate"].max(), len(per_chamber),
        )
    if not kruskal.empty and kruskal["p_value"].notna().any():
        worst = kruskal.loc[kruskal["p_value"].idxmin()]
        logger.info("  kleinster p-Wert (Batch-Test): %.4g", worst["p_value"])


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Quantitative Validierung der Mutter/Bud-Heuristik aus lineage.py.",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__,
    )
    parser.add_argument("--cells", type=Path, default=None,
                        help="Zellpositionen (Default: OUTPUT_DIR/00_cell_positions.parquet aus config.py)")
    parser.add_argument("--lineage-events", type=Path, nargs="*", default=None,
                        help="Budding-Events, auch mehrere Dateien. Default: ALLE "
                             "'20_budding_events.csv' unterhalb von OUTPUT_DIR - run_analysis.py "
                             "schreibt je eine fuer die Oszillations- und die statischen Daten, "
                             "waehrend 00_cell_positions.parquet beide enthaelt. Nur die eine "
                             "Datei zu uebergeben laesst die andere Haelfte faelschlich als "
                             "'nicht zugeordnet' erscheinen.")
    parser.add_argument("--out-dir", type=Path, default=None,
                        help="Ausgabeordner (Default: OUTPUT_DIR/lineage_validation)")
    parser.add_argument("--bud-size-threshold", type=float, default=None,
                        help="Schwelle des Groessenkriteriums (bud_area / mother_area beim ersten "
                             "Auftreten). Default: aus OUTPUT_DIR/20_bud_size_threshold.csv, also "
                             "genau die Schwelle, mit der run_analysis.py die Events erzeugt hat; "
                             "fehlt die Datei, wird sie aus den Kandidaten hier abgeleitet. "
                             "'inf' schaltet das Kriterium ab.")
    args = parser.parse_args()

    from config import LINEAGE_PARAMS, OUTPUT_DIR, FREQ_ORDER, log_active_configuration

    log_active_configuration()

    cells_path = args.cells or OUTPUT_DIR / "00_cell_positions.parquet"
    out_dir = args.out_dir or OUTPUT_DIR / "lineage_validation"

    if not cells_path.exists():
        raise FileNotFoundError(
            f"'{cells_path}' nicht gefunden. Erst run_analysis.py laufen lassen - "
            f"die Datei entsteht dort in Schritt 2c."
        )

    event_paths = args.lineage_events
    if not event_paths:
        event_paths = sorted(OUTPUT_DIR.rglob("20_budding_events.csv"))
        if not event_paths:
            logger.warning("Keine '20_budding_events.csv' unterhalb von %s gefunden.", OUTPUT_DIR)

    cells = _read_table(cells_path)
    if "in_lineage_window" in cells.columns:
        # Die Heuristik lief nur im Sparse-Phase-Fenster (relink.py); Kandidaten
        # ausserhalb waeren hier lauter "nicht zugeordnete" Fragmente.
        n_all = len(cells)
        cells = cells[cells["in_lineage_window"].astype(bool)]
        logger.info(
            "Sparse-Phase-Fenster (Spalte in_lineage_window): %d von %d Zellzeilen in %d Kammern - "
            "nur dort lief die Heuristik, nur dort wird validiert.",
            len(cells), n_all, cells["exp_id"].nunique(),
        )
    frames = [_read_table(p) for p in event_paths if Path(p).exists()]
    frames = [f for f in frames if not f.empty]
    lineage_events = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    logger.info(
        "Geladen: %d Zellzeilen, %d Budding-Events aus %d Datei(en): %s",
        len(cells), len(lineage_events), len(frames), [str(p) for p in event_paths],
    )
    if lineage_events.empty:
        logger.warning("Keine Budding-Events geladen - alle Kandidaten gelten als nicht zugeordnet.")

    bud_size_threshold = args.bud_size_threshold
    if bud_size_threshold is None:
        bud_size_threshold = load_bud_size_threshold(OUTPUT_DIR)

    run_validation(cells, lineage_events, LINEAGE_PARAMS, out_dir, freq_order=FREQ_ORDER,
                   bud_size_threshold=bud_size_threshold)
    logger.info("=== Fertig. Ergebnisse in: %s ===", out_dir)


if __name__ == "__main__":
    main()
