"""
area_growth.py

Flächenbasierte spezifische Wachstumsrate µ_area pro Einzelzell-Track,
GETRENNT nach Zelltyp (Mutter / Knospe) und mit Erkennung von
Teilungs-/Cluster-Events.

Prinzip
-------
Exponentielles Wachstum der projizierten Zellfläche:
    A(t) = A0 · exp(µ · t)      →      ln A(t) = ln A0 + µ · t
Ein Fit von ln(area) gegen die Zeit (in Stunden) liefert µ als Steigung.

Zelltyp-Trennung
----------------
- Mütter (aus lineage.identify_mothers): wachsen kontinuierlich, produzieren
  Knospen. µ_area misst Biomasse-Zunahme unabhängig von der Reproduktion.
- Knospen/Töchter: wachsen von klein → groß. µ_area misst Biomasse-Aufbau.
- Der Vergleich µ_event vs. µ_area zeigt die Ressourcen-Allokation.

Cluster-/Teilungs-Erkennung
----------------------------
Swollen Cells teilen sich irgendwann und bilden Cluster. Das verfälscht den
Fit, weil die Area dann mehrere Zellen repräsentiert. Wir detektieren:
1. Plötzlicher Flächen-Einbruch (> drop_threshold): Zelle hat sich geteilt,
   Töchter wurden separat segmentiert → Track endet hier.
2. Plötzlicher Flächen-Sprung (> jump_threshold): Cluster-Bildung oder
   Segmentierungsfehler → Frames ab hier ausschließen.
3. Morphologie-Bruch: Solidity fällt unter min_solidity ODER Eccentricity
   steigt über max_eccentricity → Cluster-Indikator.

Der Fit wird NUR auf Frames vor dem ersten Breakpoint berechnet.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Pixel-Kalibrierung: 1 µm = 13.63 px  →  1 px² = (1/13.63)² µm² ≈ 0.00538
PX_TO_UM2 = (1.0 / 13.63) ** 2


def _detect_breakpoint(
    frames: np.ndarray,
    areas: np.ndarray,
    solidities: Optional[np.ndarray] = None,
    eccentricities: Optional[np.ndarray] = None,
    drop_threshold: float = 0.6,
    jump_threshold: float = 1.8,
    min_solidity: float = 0.85,
    max_eccentricity: float = 0.95,
) -> int:
    """
    Findet den ersten Frame-Index, ab dem der Track nicht mehr eine einzelne
    wachsende Zelle repräsentiert.

    Returns
    -------
    Index des ersten „bad" Frames (0-basiert, relativ zum Input-Array).
    len(areas) wenn kein Breakpoint gefunden wurde (ganzer Track ist ok).
    """
    n = len(areas)
    if n < 3:
        return n  # zu kurz für Breakpoint-Detektion

    for i in range(1, n):
        ratio = areas[i] / areas[i - 1] if areas[i - 1] > 0 else 1.0

        # Szenario 1: plötzlicher Einbruch (Teilung, Töchter separiert)
        if ratio < drop_threshold:
            return i

        # Szenario 2: plötzlicher Sprung (Cluster-Bildung oder Segmentierungsfehler)
        if ratio > jump_threshold:
            return i

        # Szenario 3: Morphologie-Bruch (Cluster-Indikator)
        if solidities is not None and solidities[i] < min_solidity:
            return i
        if eccentricities is not None and eccentricities[i] > max_eccentricity:
            return i

    return n  # kein Breakpoint


def _fit_ln_area(t_h: np.ndarray, ln_area: np.ndarray, method: str = "theilsen") -> dict:
    """
    Fittet ln_area = intercept + mu · t_h.
    method: "theilsen" (robust gegen einzelne Ausreißer) oder "ols".
    """
    n = len(t_h)
    if n < 2:
        return {"mu": np.nan, "intercept": np.nan, "r_squared": np.nan, "n_points": n}

    if method == "theilsen":
        from scipy.stats import theilslopes
        slope, intercept, _, _ = theilslopes(ln_area, t_h)
    else:
        slope, intercept = np.polyfit(t_h, ln_area, 1)

    pred = intercept + slope * t_h
    ss_res = np.sum((ln_area - pred) ** 2)
    ss_tot = np.sum((ln_area - np.mean(ln_area)) ** 2)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    return {"mu": float(slope), "intercept": float(intercept),
            "r_squared": float(r_squared), "n_points": int(n)}


def compute_area_growth_rate(
    cells: pd.DataFrame,
    min_per_frame: float,
    lineage_events: Optional[pd.DataFrame] = None,
    mothers: Optional[pd.DataFrame] = None,
    min_frames: int = 2,
    min_area_px: float = 250.0,
    max_area_px: float = 30000.0,
    min_r_squared: float = 0.5,
    fit_method: str = "theilsen",
    drop_threshold: float = 0.6,
    jump_threshold: float = 1.8,
    min_solidity: float = 0.85,
    max_eccentricity: float = 0.95,
) -> pd.DataFrame:
    """
    Berechnet µ_area für jeden Track, getrennt nach Zelltyp (Mutter/Knospe).

    Parameters
    ----------
    cells : DataFrame mit 'frame', 'area', 'cell_uid'/'exp_id', optional
            'solidity', 'eccentricity'.
    min_per_frame : Minuten pro Frame.
    lineage_events : Ergebnis von classify_mother_bud() – wird genutzt, um
                     Tracks als „Mutter" oder „Knospe" zu klassifizieren.
                     None = alle Tracks als „unknown".
    mothers : Ergebnis von identify_mothers() – Liste von Mutter-cell_uids.
              Wird bevorzugt genutzt, wenn lineage_events None ist.
    min_frames : Mindest-Frames pro Track für einen Fit.
    min_area_px / max_area_px : Flächenfilter (px²).
    min_r_squared : Mindest-R²; darunter → fit_is_reliable=False.
    fit_method : "theilsen" oder "ols".
    drop_threshold : Flächenverhältnis < dieser Wert → Teilungs-Event.
    jump_threshold : Flächenverhältnis > dieser Wert → Cluster/Sprung.
    min_solidity : Solidity unter diesem Wert → Cluster-Indikator.
    max_eccentricity : Eccentricity über diesem Wert → Cluster-Indikator.

    Returns
    -------
    DataFrame, eine Zeile pro Track:
        exp_id, cell_uid, cell_type (mother/bud/unknown),
        mu_area [h^-1], intercept, r_squared, n_points_total, n_points_used,
        fit_is_reliable, breakpoint_frame, mean_area_um2, + Metadaten
    """
    required = {"frame", "area"}
    missing = required - set(cells.columns)
    if missing:
        raise ValueError(f"compute_area_growth_rate() fehlen Spalten: {missing}")

    df = cells.copy()
    if "cell_uid" not in df.columns:
        if {"track_id", "exp_id"} <= set(df.columns):
            df["cell_uid"] = df["exp_id"].astype(str) + "_" + df["track_id"].astype(str)
        else:
            raise ValueError("Weder 'cell_uid' noch ('track_id'+'exp_id') vorhanden.")

    # --- Zelltyp-Klassifikation ---
    mother_uids = set()
    bud_uids = set()
    if mothers is not None and "cell_uid" in mothers.columns:
        mother_uids = set(mothers["cell_uid"].dropna().unique())
    elif lineage_events is not None and "mother_cell_uid" in lineage_events.columns:
        mother_uids = set(lineage_events["mother_cell_uid"].dropna().unique())
    if lineage_events is not None and "bud_cell_uid" in lineage_events.columns:
        bud_uids = set(lineage_events["bud_cell_uid"].dropna().unique())

    t_h_factor = min_per_frame / 60.0

    results, n_too_short, n_unreliable, n_cluster = [], 0, 0, 0
    for uid, grp in df.groupby("cell_uid"):
        exp_id = grp["exp_id"].iloc[0] if "exp_id" in grp.columns else np.nan

        # Zelltyp bestimmen
        if uid in mother_uids:
            cell_type = "mother"
        elif uid in bud_uids:
            cell_type = "bud"
        else:
            cell_type = "unknown"

        # 'solidity'/'eccentricity' are optional (see docstring) - select only what
        # is actually present, otherwise this raises KeyError on datasets that
        # carry neither, and the `if "solidity" in g.columns` guards below would
        # never be reachable.
        morphology_cols = [c for c in ("solidity", "eccentricity") if c in grp.columns]
        g = grp[["frame", "area"] + morphology_cols].dropna(subset=["frame", "area"])
        g = g.sort_values("frame")
        g = g[(g["area"] >= min_area_px) & (g["area"] <= max_area_px)]

        if len(g) < min_frames:
            n_too_short += 1
            continue

        frames_arr = g["frame"].to_numpy(dtype=float)
        areas_arr = g["area"].to_numpy(dtype=float)
        sol_arr = g["solidity"].to_numpy(dtype=float) if "solidity" in g.columns else None
        ecc_arr = g["eccentricity"].to_numpy(dtype=float) if "eccentricity" in g.columns else None

        # Breakpoint detektieren (Cluster/Teilung)
        bp_idx = _detect_breakpoint(
            frames_arr, areas_arr, sol_arr, ecc_arr,
            drop_threshold=drop_threshold,
            jump_threshold=jump_threshold,
            min_solidity=min_solidity,
            max_eccentricity=max_eccentricity,
        )
        breakpoint_frame = int(frames_arr[bp_idx]) if bp_idx < len(frames_arr) else None
        if bp_idx < len(frames_arr):
            n_cluster += 1

        # Nur Frames VOR dem Breakpoint fitten
        areas_fit = areas_arr[:bp_idx]
        frames_fit = frames_arr[:bp_idx]

        if len(areas_fit) < min_frames:
            n_too_short += 1
            continue

        t_h = frames_fit * t_h_factor
        area_um2 = areas_fit * PX_TO_UM2
        ln_area = np.log(area_um2)

        fit = _fit_ln_area(t_h, ln_area, method=fit_method)
        reliable = (not np.isnan(fit["r_squared"])) and (fit["r_squared"] >= min_r_squared)
        if not reliable:
            n_unreliable += 1

        results.append({
            "exp_id": exp_id,
            "cell_uid": uid,
            "cell_type": cell_type,
            "mu_area": fit["mu"],
            "intercept": fit["intercept"],
            "r_squared": fit["r_squared"],
            "n_points_total": int(len(areas_arr)),
            "n_points_used": fit["n_points"],
            "fit_is_reliable": reliable,
            "breakpoint_frame": breakpoint_frame,
            "mean_area_um2": float(np.mean(area_um2)),
        })

    if n_too_short:
        logger.info("%d Tracks hatten zu wenige Frames (<%d) – übersprungen.", n_too_short, min_frames)
    if n_unreliable:
        logger.info("%d Fits haben R² < %.2f → fit_is_reliable=False.", n_unreliable, min_r_squared)
    if n_cluster:
        logger.info("%d Tracks hatten einen Breakpoint (Cluster/Teilung) – Fit nur auf Pre-Breakpoint-Frames.", n_cluster)

    out = pd.DataFrame(results)
    if out.empty:
        logger.warning("Keine µ_area-Werte berechnet.")
        return out

    # Metadaten anhängen
    meta_cols = [c for c in ["biosensor", "osc_type", "osc_freq", "condition", "replicate", "chamber"]
                 if c in cells.columns]
    if meta_cols and "exp_id" in cells.columns:
        meta = cells[["exp_id"] + meta_cols].drop_duplicates("exp_id")
        out = out.merge(meta, on="exp_id", how="left")

    logger.info(
        "µ_area berechnet: %d Tracks (%d Mütter, %d Knospen, %d unbekannt; %d zuverlässig).",
        len(out),
        int((out["cell_type"] == "mother").sum()),
        int((out["cell_type"] == "bud").sum()),
        int((out["cell_type"] == "unknown").sum()),
        int(out["fit_is_reliable"].sum()),
    )
    return out


def summarise_area_growth(
    area_table: pd.DataFrame,
    group_cols: Optional[list] = None,
    exclude_unreliable: bool = True,
    cell_type: Optional[str] = None,
) -> pd.DataFrame:
    """
    Aggregiert µ_area pro Bedingung.

    Parameters
    ----------
    cell_type : None = alle, "mother" = nur Mütter, "bud" = nur Knospen.
    """
    if area_table.empty:
        return pd.DataFrame()
    df = area_table
    if exclude_unreliable:
        df = df[df["fit_is_reliable"]]
    if cell_type is not None:
        df = df[df["cell_type"] == cell_type]
    if df.empty:
        return pd.DataFrame()

    if group_cols is None:
        group_cols = [c for c in ["biosensor", "osc_type", "osc_freq", "condition"] if c in df.columns]

    summary = (
        df.groupby(group_cols)
        .agg(
            mean_mu_area=("mu_area", "mean"),
            sd_mu_area=("mu_area", "std"),
            n_values=("mu_area", "count"),
            mean_r_squared=("r_squared", "mean"),
            frac_with_breakpoint=("breakpoint_frame", lambda x: x.notna().mean()),
        )
        .reset_index()
    )
    if "condition" in summary.columns:
        from growth_rate import classify_condition_type
        summary["condition_type"] = summary["condition"].map(classify_condition_type)
    return summary

# In summary_plots.py oder als neue Datei growth_scatter.py

def _paired_single_cell_stats(
    mu_event_table: Optional[pd.DataFrame],
    mu_area_table: Optional[pd.DataFrame],
    facet_col: Optional[str],
    exclude_artefacts: bool = True,
) -> dict:
    """
    Pairt µ_event und µ_area auf Einzelzell-Ebene (statt Bedingungsmittel)
    für Pearson-Korrelation und gepaarten Wilcoxon-Test.

    `mu_event_table` (Output von compute_specific_growth_rate()) hat eine
    Zeile pro REPRODUKTIONS-INTERVALL, d.h. eine Mutterzelle mit z.B. 3
    Ereignissen liefert 2 Zeilen. Für die Paarung "eine Zelle = ein
    Wertepaar" wird pro Mutterzelle (exp_id, mother_cell_uid) zuerst der
    Mittelwert über ihre Intervalle gebildet (mu_is_artefact==True Zeilen
    werden vorher ausgeschlossen, sofern exclude_artefacts=True - siehe
    Docstring von compute_specific_growth_rate()).

    Danach Join über (exp_id, cell_uid) mit `mu_area_table` (eine Zeile pro
    Track, wird auf cell_type=="mother" gefiltert, da µ_event nur für
    Mütter definiert ist).

    Returns
    -------
    dict: {facet_value: {"n": int, "r": float, "p_corr": float,
                          "p_wilcoxon": float}}, Schlüssel fehlen, wenn n
          für den jeweiligen Test zu klein ist. Leeres dict, wenn Pairing
          nicht möglich ist (fehlende Spalten/Tabellen).
    """
    if mu_event_table is None or mu_area_table is None:
        return {}
    if mu_event_table.empty or mu_area_table.empty:
        return {}
    from scipy.stats import pearsonr, wilcoxon

    et = mu_event_table
    required_event_cols = {"mu", "mother_cell_uid", "exp_id"}
    missing = required_event_cols - set(et.columns)
    if missing:
        logger.warning(
            "plot_mu_event_vs_mu_area(): mu_event_table fehlen Spalten %s "
            "(erwartet Output von compute_specific_growth_rate()) - "
            "Einzelzell-Statistik übersprungen.", missing,
        )
        return {}
    if "mu_area" not in mu_area_table.columns:
        logger.warning(
            "plot_mu_event_vs_mu_area(): mu_area_table fehlt Spalte 'mu_area' "
            "(erwartet Output von compute_area_growth_rate()) - "
            "Einzelzell-Statistik übersprungen."
        )
        return {}

    if exclude_artefacts and "mu_is_artefact" in et.columns:
        n_before = len(et)
        et = et[~et["mu_is_artefact"]]
        logger.info(
            "plot_mu_event_vs_mu_area(): %d von %d µ_event-Intervallen als Artefakt "
            "ausgeschlossen (mu_is_artefact) vor Einzelzell-Statistik.",
            n_before - len(et), n_before,
        )
    if et.empty:
        return {}

    # Eine Mutterzelle kann mehrere Intervalle (Zeilen) haben - für "eine
    # Zelle = ein Wertepaar" wird pro Mutter über ihre Intervalle gemittelt.
    cell_mu = (
        et.groupby(["exp_id", "mother_cell_uid"], as_index=False)
        .agg(mu=("mu", "mean"))
        .rename(columns={"mother_cell_uid": "cell_uid"})
    )

    at = mu_area_table
    if "cell_type" in at.columns:
        at = at[at["cell_type"] == "mother"]

    join_keys = [k for k in ("exp_id", "cell_uid") if k in cell_mu.columns and k in at.columns]
    if "cell_uid" not in join_keys:
        logger.warning(
            "plot_mu_event_vs_mu_area(): keine gemeinsame Spalte 'cell_uid' zwischen "
            "aggregierter mu_event_table und mu_area_table - Einzelzell-Statistik übersprungen."
        )
        return {}

    at_cols = join_keys + ["mu_area"] + ([facet_col] if facet_col and facet_col in at.columns else [])
    cells = cell_mu.merge(at[at_cols], on=join_keys, how="inner", validate="1:1")
    cells = cells.dropna(subset=["mu", "mu_area"])
    if cells.empty:
        return {}

    groups = cells.groupby(facet_col) if facet_col and facet_col in cells.columns else [(None, cells)]
    stats = {}
    for fval, g in groups:
        entry = {"n": len(g)}
        if len(g) >= 3:
            r, p_corr = pearsonr(g["mu"], g["mu_area"])
            entry["r"], entry["p_corr"] = r, p_corr
        diffs = g["mu"] - g["mu_area"]
        if len(g) >= 6 and np.any(diffs != 0):
            try:
                _, p_wilcoxon = wilcoxon(diffs)
                entry["p_wilcoxon"] = p_wilcoxon
            except ValueError:
                pass  # z.B. alle Differenzen 0
        stats[fval] = entry
    return stats


def plot_mu_event_vs_mu_area(
    mu_event_summary: pd.DataFrame,
    mu_area_summary: pd.DataFrame,
    out_path,
    label_col: str = "osc_freq",
    facet_col: Optional[str] = "biosensor",
    color_col: Optional[str] = "condition_type",
    join_cols: Optional[list] = None,
    mu_event_table: Optional[pd.DataFrame] = None,
    mu_area_table: Optional[pd.DataFrame] = None,
    exclude_artefacts: bool = True,
) -> None:
    """
    Scatter: µ_event (x) vs. µ_area (y), ein Punkt pro Bedingung.
    Zeigt die Entkopplung von Reproduktion und Biomasse-Wachstum.

    Punkte werden mit Fehlerbalken (sd_mu / sd_mu_area, falls vorhanden)
    dargestellt. Die gestrichelte Diagonale markiert µ_event == µ_area
    (gekoppeltes Wachstum); Achsen sind pro Facette gleich skaliert
    (aspect="equal"), damit die Diagonale visuell tatsächlich 45° entspricht.

    Pro Facette wird zusätzlich im Titel eine Statistik angegeben, wahlweise
    auf zwei Ebenen:

    - EINZELZELL-Ebene (bevorzugt, wenn `mu_event_table` UND `mu_area_table`
      übergeben werden): µ_event und µ_area werden pro Mutterzelle
      (cell_uid, exp_id) gepaart - Pearson-r und gepaarter Wilcoxon-Test
      laufen dann über alle einzelnen Zellen statt über Bedingungsmittel.
      Das erhöht die Teststärke deutlich, da n = Anzahl Zellen statt Anzahl
      Bedingungen ist. µ_area wird dafür auf cell_type=="mother" gefiltert,
      da µ_event nur für Mutterzellen definiert ist (Reproduktion).
    - BEDINGUNGS-Ebene (Fallback, wenn die Roh-Tabellen fehlen): wie bisher
      auf den Bedingungsmitteln (mean_mu/mean_mu_area) aus den Summaries.

    In beiden Fällen: Pearson-r+p ab n>=3 Punkten, Wilcoxon-p ab n>=6 Punkten
    (sonst als "n.v." markiert, da nicht aussagekräftig).

    Parameters
    ----------
    mu_event_summary : Output von summarise_growth_rate() (Spalten u.a.
                        mean_mu, sd_mu, ...). Bestimmt weiterhin die
                        geplotteten Punkte (ein Punkt pro Bedingung).
    mu_area_summary : Output von summarise_area_growth() (Spalten u.a.
                       mean_mu_area, sd_mu_area, ...).
    out_path : Pfad der Ausgabedatei (Elternverzeichnis wird bei Bedarf
               angelegt).
    label_col : Spalte, deren Werte neben die Punkte annotiert werden.
    facet_col : Spalte für Subplots nebeneinander (None = ein Panel).
    color_col : Spalte für Punktfarbe/Legende (None = einfarbig).
    join_cols : Spalten, auf denen event- und area-Summary gemerged
                werden. None = automatisch alle gemeinsamen Spalten außer
                den µ-Statistik-Spalten selbst.
    mu_event_table : Optional. Roh-Output von compute_specific_growth_rate()
                      (eine Zeile pro Reproduktions-INTERVALL einer
                      Mutterzelle, Spalten u.a. 'mu', 'mother_cell_uid',
                      'exp_id', 'mu_is_artefact'). Zusammen mit
                      `mu_area_table` aktiviert dies die Einzelzell-Statistik;
                      Mutterzellen mit mehreren Intervallen werden vor dem
                      Pairing über ihre Intervalle gemittelt.
    mu_area_table : Optional. Roh-Output von compute_area_growth_rate()
                     (eine Zeile pro Track, Spalten u.a. 'mu_area',
                     'cell_type', 'cell_uid', 'exp_id').
    exclude_artefacts : ob mu_is_artefact==True Intervalle aus
                     `mu_event_table` vor der Einzelzell-Statistik
                     ausgeschlossen werden (Standard True, analog zu
                     summarise_growth_rate() - siehe growth_rate.py).

    Returns
    -------
    None. Schreibt die Abbildung nach `out_path`.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy.stats import pearsonr, wilcoxon

    if join_cols is None:
        join_cols = [c for c in mu_event_summary.columns
                     if c in mu_area_summary.columns
                     and c not in ("mean_mu", "sd_mu", "mean_mu_area", "sd_mu_area", "n_values")]

    merged = mu_event_summary.merge(
        mu_area_summary, on=join_cols, suffixes=("_event", "_area"),
        how="inner", validate="1:1", indicator=True,
    )
    n_event_only = (mu_event_summary.merge(mu_area_summary, on=join_cols, how="left", indicator=True)["_merge"] == "left_only").sum()
    n_area_only = (mu_area_summary.merge(mu_event_summary, on=join_cols, how="left", indicator=True)["_merge"] == "left_only").sum()
    if n_event_only or n_area_only:
        logger.warning(
            "plot_mu_event_vs_mu_area(): %d Bedingungen nur in µ_event, %d nur in µ_area "
            "vorhanden (join_cols=%s) - werden im Scatter nicht angezeigt.",
            n_event_only, n_area_only, join_cols,
        )
    merged = merged.drop(columns=["_merge"])
    if merged.empty:
        logger.warning("plot_mu_event_vs_mu_area(): kein Overlap zwischen µ_event und µ_area - kein Plot erzeugt.")
        return

    facets = sorted(merged[facet_col].dropna().unique()) if facet_col and facet_col in merged.columns else [None]
    has_color = bool(color_col) and color_col in merged.columns and merged[color_col].notna().any()
    colors_vals = sorted(merged[color_col].dropna().unique()) if has_color else [None]
    palette = dict(zip(colors_vals, sns.color_palette("Set2", n_colors=max(len(colors_vals), 1))))

    fig, axes = plt.subplots(1, len(facets), figsize=(5 * len(facets), 5), squeeze=False)
    axes = axes[0]

    cell_stats = _paired_single_cell_stats(mu_event_table, mu_area_table, facet_col, exclude_artefacts=exclude_artefacts)
    using_cell_level = bool(cell_stats)
    if mu_event_table is not None and mu_area_table is not None and not using_cell_level:
        logger.warning(
            "plot_mu_event_vs_mu_area(): mu_event_table/mu_area_table übergeben, aber "
            "Pairing fehlgeschlagen - falle zurück auf Bedingungs-Statistik."
        )

    for ax, facet in zip(axes, facets):
        sub = merged[merged[facet_col] == facet] if facet_col and facet_col in merged.columns else merged
        for cv in colors_vals:
            sub_c = sub[sub[color_col] == cv] if has_color else sub
            if sub_c.empty:
                continue
            xerr = sub_c["sd_mu"] if "sd_mu" in sub_c.columns else None
            yerr = sub_c["sd_mu_area"] if "sd_mu_area" in sub_c.columns else None
            ax.errorbar(
                sub_c["mean_mu"], sub_c["mean_mu_area"], xerr=xerr, yerr=yerr,
                fmt="o", color=palette.get(cv), label=str(cv) if cv is not None else None,
                markersize=6, alpha=0.8, capsize=2, linewidth=1, elinewidth=1,
            )
            # Labels (osc_freq) neben die Punkte
            if label_col in sub_c.columns:
                for _, row in sub_c.iterrows():
                    ax.annotate(str(row[label_col]), (row["mean_mu"], row["mean_mu_area"]),
                                textcoords="offset points", xytext=(4, 4), fontsize=6)

        # Diagonale: µ_event == µ_area (gekoppeltes Wachstum).
        # Limits VOR dem Zeichnen fixieren und danach explizit erneut setzen,
        # da Matplotlib beim Plotten sonst per Auto-Margin leicht nachskaliert
        # und die Linie nicht mehr exakt in den Ecken sitzt.
        lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
                max(ax.get_xlim()[1], ax.get_ylim()[1])]
        ax.plot(lims, lims, "k--", alpha=0.3, linewidth=0.8)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_aspect("equal", adjustable="box")

        ax.set_xlabel("µ_event [h⁻¹] (reproduction)")
        ax.set_ylabel("µ_area [h⁻¹] (biomass)")

        # Statistik pro Facette: Pearson-Korrelation (Zusammenhang µ_event/µ_area)
        # und gepaarter Wilcoxon-Vorzeichen-Rang-Test auf die Differenz, als
        # quantitativer Beleg für/gegen Entkopplung - ergänzt die rein visuelle
        # Abweichung von der Diagonale. Bevorzugt auf Einzelzell-Ebene (mehr
        # Teststärke), sonst Fallback auf Bedingungsmittel.
        title = str(facet) if facet else ""
        if using_cell_level:
            entry = cell_stats.get(facet, {"n": 0})
            n_stat, level_tag = entry["n"], "cells"
        else:
            entry = {}
            n_stat, level_tag = len(sub), "conditions"
            if n_stat >= 3:
                r, p_corr = pearsonr(sub["mean_mu"], sub["mean_mu_area"])
                entry["r"], entry["p_corr"] = r, p_corr
            diffs = sub["mean_mu"] - sub["mean_mu_area"]
            if n_stat >= 6 and np.any(diffs != 0):
                try:
                    _, p_wilcoxon = wilcoxon(diffs)
                    entry["p_wilcoxon"] = p_wilcoxon
                except ValueError:
                    pass  # z.B. alle Differenzen 0

        if "r" in entry:
            title += f"\nr={entry['r']:.2f} (p={entry['p_corr']:.3f}, n={n_stat} {level_tag})"
        else:
            title += f"\nr: n/a (n={n_stat} {level_tag} < 3)"
        if "p_wilcoxon" in entry:
            title += f", Wilcoxon p={entry['p_wilcoxon']:.3f}"
        ax.set_title(title, fontsize=9)

    if has_color:
        # Handles/Labels über ALLE Facetten sammeln und deduplizieren, da eine
        # einzelne Facette nicht zwingend alle color_col-Werte enthält.
        handles_by_label = {}
        for ax in axes:
            h, l = ax.get_legend_handles_labels()
            for hi, li in zip(h, l):
                handles_by_label.setdefault(li, hi)
        fig.legend(handles_by_label.values(), handles_by_label.keys(),
                   loc="lower center", ncol=len(handles_by_label))
    fig.suptitle("Reproduction vs. biomass growth", y=1.02)
    fig.tight_layout()

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight", dpi=150)
    plt.close(fig)