"""growth_from_budding.py - spezifische Wachstumsrate der POPULATION aus den Knospungen.

Die Tabelle 11_specific_growth_rate.csv (growth_rate.py) rechnet µ = ln(1+k)/Δt aus dem Abstand
zweier Knospungen DERSELBEN Mutter (Bloebaum et al. 2024, Eq. 2). Auf den v12-Tabellen liegt der
Median dieses Abstands bei 1.3 h und jeder zehnte bei einem Frame; die Bedingungsmittel kommen auf
1.6-2.9 h^-1. Das ist eine Knospen-Abgabefrequenz in log-Einheiten, keine Wachstumsrate: die
Generationszeit der Toechter geht nirgends ein.

HIER: Geburten je Zellstunde im Sparse-Phase-Fenster einer Kammer.

    mu_bud = n_births / cell_hours

    n_births    = akzeptierte Knospungs-Events (20_budding_events.csv: gemessene Mutter aus der
                  Maskenberuehrung oder Radius-Heuristik, Persistenz, Groessenkriterium)
    cell_hours  = Summe ueber die Frames des Fensters von (Zellen im Frame) x MIN_PER_FRAME / 60

In balanciertem exponentiellem Wachstum ohne Tod gilt dN/dt = mu N, und jede Geburt erhoeht N um
eins - Geburten je Zellstunde SIND mu. Ein Mehrfach-Knospungs-Ereignis zaehlt jede Knospe, ein
Intervall wird nicht gebraucht, und der Nenner sind ALLE Zellen (Muetter und Toechter), nicht nur
die Muetter wie bei budding_rate_per_h (buds per mother-hour, lineage.compute_budding_ratio).

DANEBEN: Einwanderung. Ein neuer Track ohne Elternmaske (link_type 'new' der Pipeline v12) ist
eine angespuelte Zelle, keine Geburt. immigration_per_cell_h = solche Tracks je Zellstunde -
damit sieht man, wie viel vom Zuwachs der Objektzahl Zufluss ist. Ohne link_type-Spalte (v11-
Tabellen) bleibt die Spalte NaN. Auswanderung und Tod sind nicht zaehlbar (eine tote Zelle bleibt
als Truemmer liegen und wird vom Zellfilter entfernt, siehe cell_filter.py); mu_bud ist deshalb
die spezifische GEBURTENrate, die obere Schranke der Nettowachstumsrate.

Nur das Sparse-Phase-Fenster: im dichten Feld beruehrt jedes neue Objekt eine Zelle, und die
Events sind Fragmentstatistik (docs/data_story.md 2.2). Die Fensterlaenge steckt in cell_hours,
ein reiner Zaehler waere mit ihr konfundiert.
"""
from __future__ import annotations

import logging
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

META_COLS = ["biosensor", "osc_type", "osc_freq", "condition", "replicate", "chamber",
             "chip", "chip_family", "medium", "date", "culture"]


def compute_growth_from_budding(
    cells_window: pd.DataFrame,
    lineage_events: pd.DataFrame,
    min_per_frame: float,
) -> pd.DataFrame:
    """Eine Zeile je Kammer (exp_id): n_births, n_immigrants, cell_hours, mu_bud,
    immigration_per_cell_h, doubling_time_h, n_frames_window, cells_per_frame_mean.

    cells_window : Zellzeilen IM Sparse-Phase-Fenster (PipelineContext.cells_lineage), nach
                   Zellfilter und QC; braucht exp_id, cell_uid, frame.
    lineage_events : Ergebnis von classify_mother_bud() auf denselben Zellen.
    """
    if cells_window is None or cells_window.empty or "exp_id" not in cells_window.columns:
        return pd.DataFrame()
    hours_per_frame = float(min_per_frame) / 60.0

    per_frame = cells_window.groupby(["exp_id", "frame"])["cell_uid"].nunique().rename("n_cells").reset_index()
    per_exp = per_frame.groupby("exp_id").agg(
        n_frames_window=("frame", "nunique"),
        cells_per_frame_mean=("n_cells", "mean"),
        cell_frames=("n_cells", "sum"),
    ).reset_index()
    per_exp["cell_hours"] = per_exp["cell_frames"] * hours_per_frame

    if lineage_events is not None and not lineage_events.empty and "exp_id" in lineage_events.columns:
        births = lineage_events.groupby("exp_id").size().rename("n_births")
    else:
        births = pd.Series(dtype=int, name="n_births")
    per_exp["n_births"] = per_exp["exp_id"].map(births).fillna(0).astype(int)

    # Einwanderung: erster Frame eines Tracks im Fenster, aber nicht der erste Frame des Fensters,
    # und ohne Elternmaske (link_type 'new'). Nur mit link_type (Pipeline v12) bestimmbar.
    if "link_type" in cells_window.columns:
        first = (cells_window.sort_values("frame")
                 .groupby("cell_uid", sort=False)
                 .agg(exp_id=("exp_id", "first"), first_frame=("frame", "first"), link_type=("link_type", "first"))
                 .reset_index())
        window_start = cells_window.groupby("exp_id")["frame"].min()
        first["at_window_start"] = first["first_frame"] <= first["exp_id"].map(window_start)
        immigrants = first[(first["link_type"] == "new") & ~first["at_window_start"]].groupby("exp_id").size()
        per_exp["n_immigrants"] = per_exp["exp_id"].map(immigrants).fillna(0).astype(float)
    else:
        per_exp["n_immigrants"] = np.nan

    ok = per_exp["cell_hours"] > 0
    per_exp["mu_bud"] = np.where(ok, per_exp["n_births"] / per_exp["cell_hours"].where(ok), np.nan)
    per_exp["immigration_per_cell_h"] = per_exp["n_immigrants"] / per_exp["cell_hours"].where(ok)
    with np.errstate(divide="ignore"):
        per_exp["doubling_time_h"] = np.where(per_exp["mu_bud"] > 0, np.log(2) / per_exp["mu_bud"], np.nan)

    meta_cols = [c for c in META_COLS if c in cells_window.columns]
    meta = cells_window[["exp_id"] + meta_cols].drop_duplicates("exp_id")
    out = per_exp.merge(meta, on="exp_id", how="left")

    logger.info(
        "Wachstumsrate aus Knospungen: %d Kammern, %d Geburten in %.0f Zellstunden - mu_bud Median %.3f h^-1 "
        "(q10 %.3f, q90 %.3f), Verdopplungszeit Median %.1f h%s.",
        len(out), int(out["n_births"].sum()), float(out["cell_hours"].sum()),
        float(out["mu_bud"].median()), float(out["mu_bud"].quantile(0.1)), float(out["mu_bud"].quantile(0.9)),
        float(out["doubling_time_h"].median()),
        "" if out["n_immigrants"].isna().all() else
        f"; Einwanderung Median {float(out['immigration_per_cell_h'].median()):.3f} je Zellstunde",
    )
    return out


def growth_from_budding_table_summary(per_chamber: pd.DataFrame) -> Optional[str]:
    """Kurze Textzusammenfassung fuer das Log (None bei leerer Tabelle)."""
    if per_chamber is None or per_chamber.empty:
        return None
    cols = [c for c in ["biosensor", "osc_type"] if c in per_chamber.columns]
    if not cols:
        return None
    g = per_chamber.groupby(cols)["mu_bud"].median().round(3)
    return g.to_string()


def plot_mu_bud_vs_mu_area(bud_per_chip: pd.DataFrame, area_per_chip: pd.DataFrame, out_path) -> None:
    """Scatter je Chip und Kontrollart: µ_bud (Population, Geburten je Zellstunde) gegen µ_area
    (Einzelzelle, Flaechenwachstum). Farbe = Stamm, Marker = Kontrollart, eine Facette je
    Oszillationstyp. Die Diagonale markiert gleiche Raten."""
    import matplotlib.pyplot as plt
    from scipy import stats
    from plot_style import (INK_MUTED, INK_SOFT, control_handles, finish, legend_below, marker_kwargs,
                            ordered_strains, strain_color, strain_handles)
    from endpoint_trends import add_condition_type

    if bud_per_chip is None or bud_per_chip.empty or area_per_chip is None or area_per_chip.empty:
        return
    keys = [c for c in ["biosensor", "osc_type", "osc_freq", "condition", "chip"] if c in bud_per_chip.columns
            and c in area_per_chip.columns]
    m = bud_per_chip.merge(area_per_chip[keys + ["value"]], on=keys, suffixes=("", "_area"))
    m = add_condition_type(m).dropna(subset=["value", "value_area"])
    if m.empty:
        logger.warning("plot_mu_bud_vs_mu_area(): kein Ueberlapp zwischen µ_bud und µ_area - uebersprungen.")
        return
    facets = sorted(m["osc_type"].dropna().astype(str).unique()) if "osc_type" in m.columns else ["all"]
    fig, axes = plt.subplots(1, len(facets), figsize=(3.4 * len(facets), 3.5), squeeze=False, sharex=True, sharey=True)
    lo = float(min(m["value"].min(), m["value_area"].min()))
    hi = float(max(m["value"].max(), m["value_area"].max()))
    pad = 0.05 * (hi - lo if hi > lo else 1.0)
    for ax, ot in zip(axes[0], facets):
        sub = m if ot == "all" else m[m["osc_type"].astype(str) == ot]
        for strain in ordered_strains(sub["biosensor"].unique()):
            s = sub[sub["biosensor"].astype(str) == strain]
            for ct, part in s.groupby("condition_type"):
                ax.scatter(part["value"], part["value_area"], **marker_kwargs(ct, strain_color(strain), size=26))
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color=INK_MUTED, lw=0.8, ls="--", zorder=1)
        ok = sub[["value", "value_area"]].dropna()
        if len(ok) >= 3:
            rho = stats.spearmanr(ok["value"], ok["value_area"])[0]
            ax.text(0.97, 0.03, f"Spearman ρ = {rho:+.2f} (n = {len(ok)} chips × control types)",
                    transform=ax.transAxes, va="bottom", ha="right", fontsize=7.5, color=INK_SOFT, zorder=6,
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.85))
        ax.set_title(ot if ot != "all" else "", loc="left")
        ax.set_xlabel("µ_bud [h⁻¹] (births per cell-hour, population)")
        ax.set_xlim(lo - pad, hi + pad); ax.set_ylim(lo - pad, hi + pad)
        ax.set_aspect("equal", adjustable="box")
    axes[0][0].set_ylabel("µ_area [h⁻¹] (area growth, single cells)")
    handles = strain_handles(m["biosensor"].unique()) + control_handles()
    legend_below(fig, handles, ncol=min(len(handles), 4), y=0.02)
    fig.suptitle("population growth from budding vs single-cell area growth, one point per chip and control type",
                 fontsize=9)
    finish(fig, out_path, logger)
