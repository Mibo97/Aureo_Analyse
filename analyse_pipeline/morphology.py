"""
morphology.py
=============
KUMULATIVE Wirkung der Feast/Famine-Zyklen auf die Zellmorphologie.

FRAGESTELLUNG
-------------
Die Oszillationsperioden (0.75 - 24 min) liegen ALLE unter bzw. am Nyquist-
Limit der Zeitraffer-Aufnahme (10 min/Frame -> kuerzeste aufloesbare Periode
20 min). Ein einzelner Zyklus ist damit grundsaetzlich nicht beobachtbar - und
das ist hier auch nicht das Ziel. Die Oszillation ist eine BEHANDLUNG, keine
Messgroesse:

    Alle Bedingungen erhalten (bei gleichem Tastverhaeltnis) dieselbe Gesamt-
    Feast- und Gesamt-Famine-Zeit. Sie unterscheiden sich NUR darin, wie fein
    diese Zeit zerhackt wird - von ~800 Zyklen (0.75 min) bis ~25 Zyklen
    (24 min) in 10 h, also ein 32-facher Bereich in der Anzahl der Wechsel.

Die Frage lautet damit nicht "folgt die Zelle der Oszillation?", sondern:
kostet es die Zelle etwas, dieselbe Naehrstoffmenge feiner zerhackt zu
bekommen - und summiert sich dieser Preis ueber Stunden auf?

Morphologie ist dafuer der direkteste Messwert, den die Bildverarbeitung
ohnehin liefert: 'area', 'eccentricity' und 'solidity' liegen in jeder Zeile
von Combined_Results, wurden bisher aber nirgends als Ergebnis ausgewertet
(solidity ausschliesslich intern als Breakpoint-Trigger in area_growth.py).

WELCHE RICHTUNG ZU ERWARTEN IST
-------------------------------
Entscheidend ist die Laenge der FAMINE-HALBPERIODE im Verhaeltnis zu den
internen Metabolitpools der Zelle:

    Periode  0.75 min -> Famine-Halbperiode ~22 s  -> Pools ueberbruecken das
                         muehelos; die Zelle "sieht" praktisch ein konstantes,
                         gemitteltes Medium.
    Periode    24 min -> Famine-Halbperiode 12 min -> lang genug fuer echte
                         Verarmung und eine Hungerantwort, danach
                         Reakklimatisierung. 25x in 10 h.

Die STARKE Belastung wird also bei den LANGSAMEN Zyklen erwartet, nicht bei
den schnellen - entgegen der naheliegenden Intuition "schneller = stressiger".
Wo der Uebergang liegt, ist eine Messgroesse dieser Arbeit.

ANSATZ
------
Statt eine absolute Morphologie-Klassifikation zu erfinden, wird die
"normale" Morphologie AUS DEN DATEN definiert: aus der Referenzbedingung
(per Default PosCtrl = durchgehend Feast) wird ein Morphospace-Fenster
abgeleitet (Perzentile von eccentricity/solidity/area). Jede Zelle ausserhalb
dieses Fensters gilt als 'aberrant'. Damit ist die Kennzahl

    aberranter Anteil = Anteil Zellen ausserhalb des Normalfensters

pro Kammer und Zeitpunkt definiert, in Einheiten der eigenen Kontrolle statt
in willkuerlichen Pixel-Schwellen - und ueber Bedingungen hinweg vergleichbar.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

from analysis import natural_freq_sort, pretty_label
from growth_rate import classify_condition_type

logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

sns.set_theme(style="whitegrid", context="notebook")

MORPHOLOGY_FEATURES = ["eccentricity", "solidity", "area"]

# Reihenfolge = Prioritaet der Klassifikation (siehe classify_morphotype()).
MORPHOTYPE_ORDER = ["yeast_like", "swollen", "elongated", "clustered"]
MORPHOTYPE_COLORS = {
    "yeast_like": "#4C78A8",
    "swollen": "#F2B701",
    "elongated": "#7B5EA7",
    "clustered": "#E45756",
}


@dataclass
class MorphotypeThresholds:
    """Grenzen des 'normalen' Morphospace, abgeleitet aus einer Referenzbedingung.

    max_eccentricity : darueber -> 'elongated' (laengliche/pseudohyphale Form)
    min_solidity     : darunter -> 'clustered' (konkave Umrisse: Zellcluster,
                       Verzweigung, Knospenketten)
    max_area_px      : darueber -> 'swollen' (vergroesserte Zelle)
    source           : Herkunft der Werte, wird in die CSV geschrieben und ins
                       Log geschrieben - ohne diese Angabe ist keine
                       Morphotyp-Zahl reproduzierbar.
    """
    max_eccentricity: float
    min_solidity: float
    max_area_px: float
    source: str = "manual"

    def as_row(self) -> dict:
        return {
            "max_eccentricity": self.max_eccentricity,
            "min_solidity": self.min_solidity,
            "max_area_px": self.max_area_px,
            "source": self.source,
        }


def derive_thresholds_from_reference(
    cells: pd.DataFrame,
    reference_condition_types: Sequence[str] = ("PosCtrl",),
    percentile: float = 95.0,
) -> Optional[MorphotypeThresholds]:
    """
    Leitet das Normalfenster aus der Referenzbedingung ab (Default: PosCtrl,
    durchgehend Feast - also Zellen, die KEINE Oszillation erfahren haben).

    Bewusst NICHT aus dem Gesamtdatensatz: sonst definierten die
    Oszillationszellen mit, was 'normal' heisst, und der Effekt, den man
    messen will, verschwaende teilweise in der eigenen Referenz.

    percentile=95 heisst: die oberen 5% der Referenz-Exzentrizitaet und
    -Flaeche sowie die unteren 5% der Referenz-Solidity gelten bereits als
    aberrant. Bei perfekt gleicher Morphologie erwartet man in jeder
    Bedingung also einen aberranten Anteil nahe 5-15% (drei Kriterien,
    teils korreliert) - Abweichungen davon sind das Signal.
    """
    available = [f for f in MORPHOLOGY_FEATURES if f in cells.columns]
    if len(available) < len(MORPHOLOGY_FEATURES):
        logger.warning(
            "derive_thresholds_from_reference(): Spalten %s fehlen - Morphologie-Auswertung "
            "nicht moeglich.", sorted(set(MORPHOLOGY_FEATURES) - set(available)),
        )
        return None

    if "condition" not in cells.columns:
        logger.warning("derive_thresholds_from_reference(): Spalte 'condition' fehlt.")
        return None

    ctype = cells["condition"].apply(classify_condition_type)
    reference = cells[ctype.isin(reference_condition_types)].dropna(subset=MORPHOLOGY_FEATURES)

    if reference.empty:
        # Fallback statt Abbruch: relative Vergleiche zwischen den Bedingungen
        # bleiben moeglich, nur der absolute Bezug "gegenueber ungestoerten
        # Zellen" geht verloren. Muss in der Abbildungslegende stehen.
        logger.warning(
            "Keine Referenzzellen (%s) gefunden - Normalfenster wird aus dem GESAMTEN "
            "Datensatz abgeleitet. Die aberranten Anteile sind dann nur noch relativ "
            "zwischen den Bedingungen interpretierbar, nicht als Abweichung von "
            "ungestoerten Zellen. Bitte in der Abbildungslegende vermerken.",
            list(reference_condition_types),
        )
        reference = cells.dropna(subset=MORPHOLOGY_FEATURES)
        source = f"pooled_all_conditions_p{percentile:g}"
        if reference.empty:
            return None
    else:
        source = f"{'+'.join(reference_condition_types)}_p{percentile:g}_n{len(reference)}"

    thresholds = MorphotypeThresholds(
        max_eccentricity=float(np.percentile(reference["eccentricity"], percentile)),
        min_solidity=float(np.percentile(reference["solidity"], 100.0 - percentile)),
        max_area_px=float(np.percentile(reference["area"], percentile)),
        source=source,
    )
    logger.info(
        "Morphotyp-Schwellen aus '%s': eccentricity > %.3f = elongated, solidity < %.3f = "
        "clustered, area > %.0f px² = swollen.",
        source, thresholds.max_eccentricity, thresholds.min_solidity, thresholds.max_area_px,
    )
    return thresholds


def classify_morphotype(cells: pd.DataFrame, thresholds: MorphotypeThresholds) -> pd.DataFrame:
    """
    Ergaenzt die Spalten 'morphotype' und 'is_aberrant'.

    PRIORITAET (eine Zelle bekommt genau einen Typ):
        1. clustered  - solidity unter der Schwelle. Zuerst geprueft, weil eine
           konkave Kontur meist bedeutet, dass das Objekt gar keine EINZELNE
           Zelle mehr ist; Flaeche und Exzentrizitaet sind dann ohnehin nicht
           mehr als Einzelzellmass lesbar.
        2. elongated  - eccentricity ueber der Schwelle (laengliche Form).
        3. swollen    - area ueber der Schwelle, bei ansonsten normaler Form.
        4. yeast_like - innerhalb aller Grenzen.

    Zeilen mit fehlenden Messwerten bekommen morphotype=NA und gehen NICHT als
    'nicht aberrant' in die Anteile ein (siehe summarise_*).
    """
    df = cells.copy()
    valid = df[MORPHOLOGY_FEATURES].notna().all(axis=1)

    morphotype = pd.Series(pd.NA, index=df.index, dtype="object")
    is_clustered = valid & (df["solidity"] < thresholds.min_solidity)
    is_elongated = valid & ~is_clustered & (df["eccentricity"] > thresholds.max_eccentricity)
    is_swollen = valid & ~is_clustered & ~is_elongated & (df["area"] > thresholds.max_area_px)

    morphotype[valid] = "yeast_like"
    morphotype[is_swollen] = "swollen"
    morphotype[is_elongated] = "elongated"
    morphotype[is_clustered] = "clustered"

    df["morphotype"] = morphotype
    df["is_aberrant"] = np.where(valid, morphotype != "yeast_like", np.nan)

    n_valid = int(valid.sum())
    if n_valid == 0:
        logger.warning("classify_morphotype(): keine Zeile mit vollstaendigen Morphologie-Werten.")
        return df

    counts = df.loc[valid, "morphotype"].value_counts()
    logger.info(
        "Morphotypen klassifiziert (%d von %d Zeilen mit vollstaendigen Werten): %s",
        n_valid, len(df),
        ", ".join(f"{k}={counts.get(k, 0)} ({100 * counts.get(k, 0) / n_valid:.1f}%)"
                  for k in MORPHOTYPE_ORDER),
    )
    return df


def add_switching_dose(
    cells: pd.DataFrame,
    min_per_frame: float,
    period_col: str = "osc_freq",
) -> pd.DataFrame:
    """
    Ergaenzt 'period_min' und 'n_cycles_elapsed' - die eigentliche DOSIS-Achse.

    'osc_freq' enthaelt trotz seines Namens die PERIODE in Minuten (0.75 ...
    24), nicht eine Frequenz. Hier wird sie als solche gelesen und in die
    Anzahl bereits durchlaufener Feast/Famine-Zyklen umgerechnet:

        n_cycles_elapsed = (frame * min_per_frame) / period_min

    Damit laesst sich die Morphologie gegen die tatsaechlich erfahrene Anzahl
    Wechsel auftragen statt nur gegen die Zeit - genau der kumulative Effekt,
    um den es geht. Bedingungen ohne numerische Periode (statisch, 'static_*')
    bekommen NaN.
    """
    df = cells.copy()
    if period_col not in df.columns:
        logger.warning("add_switching_dose(): Spalte '%s' fehlt.", period_col)
        return df

    df["period_min"] = pd.to_numeric(df[period_col], errors="coerce")
    n_non_numeric = int(df["period_min"].isna().sum())
    if n_non_numeric:
        logger.info(
            "add_switching_dose(): %d Zeilen ohne numerische Periode in '%s' (z.B. statische "
            "Bedingungen) - n_cycles_elapsed bleibt dort NaN.", n_non_numeric, period_col,
        )

    elapsed_min = df["frame"] * min_per_frame
    with np.errstate(divide="ignore", invalid="ignore"):
        df["n_cycles_elapsed"] = elapsed_min / df["period_min"]
    return df


def summarise_morphotype_over_time(
    cells: pd.DataFrame,
    group_cols: Optional[list[str]] = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Morphotyp-Anteile pro Kammer x Zeitpunkt, und daraus ueber Replikate
    gemittelt - zweistufig, wie in analysis._aggregate_over_replicates()
    beschrieben (Zellen sind keine unabhaengigen Replikate).

    Returns
    -------
    per_chamber : eine Zeile pro exp_id x frame, mit n_cells, frac_aberrant
                  und einer frac_<morphotype>-Spalte je Typ
    aggregated  : ueber die Kammern einer Bedingung gemittelt, mit SEM
    """
    if "morphotype" not in cells.columns:
        raise ValueError("summarise_morphotype_over_time() braucht die Spalte 'morphotype' "
                         "(siehe classify_morphotype()).")

    valid = cells[cells["morphotype"].notna()]
    if valid.empty:
        return pd.DataFrame(), pd.DataFrame()

    meta_cols = [c for c in ["biosensor", "osc_type", "osc_freq", "condition", "replicate", "chamber"]
                 if c in valid.columns]

    counts = (
        valid.groupby(["exp_id", "frame"])["morphotype"]
        .value_counts()
        .unstack(fill_value=0)
        .reset_index()
    )
    for m in MORPHOTYPE_ORDER:
        if m not in counts.columns:
            counts[m] = 0

    counts["n_cells"] = counts[MORPHOTYPE_ORDER].sum(axis=1)
    for m in MORPHOTYPE_ORDER:
        counts[f"frac_{m}"] = counts[m] / counts["n_cells"]
    counts["frac_aberrant"] = 1.0 - counts["frac_yeast_like"]

    extra = [c for c in ["time_h", "n_cycles_elapsed", "period_min"] if c in valid.columns]
    keys = valid[["exp_id", "frame"] + meta_cols + extra].drop_duplicates(["exp_id", "frame"])
    per_chamber = counts.merge(keys, on=["exp_id", "frame"], how="left")

    if group_cols is None:
        group_cols = [c for c in ["biosensor", "osc_type", "osc_freq", "condition", "frame"]
                      if c in per_chamber.columns]

    frac_cols = [f"frac_{m}" for m in MORPHOTYPE_ORDER] + ["frac_aberrant"]
    agg_spec = {c: [("mean", "mean"), ("sem", _sem)] for c in frac_cols}
    aggregated = per_chamber.groupby(group_cols, dropna=False).agg(
        n_chambers=("exp_id", "nunique"),
        n_cells=("n_cells", "sum"),
        **{f"{c}_{stat}": (c, fn) for c in frac_cols for stat, fn in agg_spec[c]},
    ).reset_index()

    for c in ["time_h", "n_cycles_elapsed", "period_min"]:
        if c in per_chamber.columns:
            lookup = per_chamber.groupby(group_cols, dropna=False)[c].mean().reset_index()
            aggregated = aggregated.merge(lookup, on=group_cols, how="left")

    return per_chamber, aggregated


def _sem(s: pd.Series) -> float:
    s = s.dropna()
    return float(s.std(ddof=1) / np.sqrt(len(s))) if len(s) > 1 else 0.0


def summarise_morphotype_endpoint(
    per_chamber: pd.DataFrame,
    last_fraction: float = 0.25,
) -> pd.DataFrame:
    """
    Kumulatives Endergebnis pro Kammer: Morphotyp-Anteile gemittelt ueber das
    LETZTE Zeitfenster (Default: letztes Viertel der Frames dieser Kammer).

    Ein einzelner Endframe waere zu verrauscht (wenige Zellen), der Mittelwert
    ueber die gesamte Messdauer wuerde dagegen den fruehen, noch unbelasteten
    Zustand mit hineinmitteln und den kumulativen Effekt verwaessern. Das
    Spaetfenster ist der Kompromiss - und die Kennzahl, die spaeter zwischen
    den Perioden verglichen wird.
    """
    if per_chamber.empty:
        return pd.DataFrame()

    df = per_chamber.copy()
    frame_max = df.groupby("exp_id")["frame"].transform("max")
    frame_min = df.groupby("exp_id")["frame"].transform("min")
    cutoff = frame_max - (frame_max - frame_min) * last_fraction
    late = df[df["frame"] >= cutoff]

    meta_cols = [c for c in ["biosensor", "osc_type", "osc_freq", "condition", "replicate", "chamber"]
                 if c in late.columns]
    frac_cols = [f"frac_{m}" for m in MORPHOTYPE_ORDER] + ["frac_aberrant"]

    endpoint = (
        late.groupby(["exp_id"] + meta_cols, dropna=False)
        .agg(n_frames=("frame", "nunique"), n_cells=("n_cells", "sum"),
             **{c: (c, "mean") for c in frac_cols})
        .reset_index()
    )
    if "condition" in endpoint.columns:
        endpoint["condition_type"] = endpoint["condition"].apply(classify_condition_type)
    logger.info(
        "Endpunkt-Morphologie: %d Kammern, Spaetfenster = letzte %.0f%% der Frames.",
        len(endpoint), 100 * last_fraction,
    )
    return endpoint


def test_aberrant_across_periods(
    endpoint: pd.DataFrame,
    value_col: str = "frac_aberrant",
    group_cols: Optional[list[str]] = None,
    period_col: str = "osc_freq",
) -> pd.DataFrame:
    """
    Kruskal-Wallis ueber die Oszillationsperioden, plus Spearman-Korrelation
    gegen die numerische Periode.

    Kruskal beantwortet "unterscheiden sich die Perioden ueberhaupt?", Spearman
    zusaetzlich "gibt es einen MONOTONEN Trend mit der Periodenlaenge?" - und
    genau die Monotonie ist die eigentliche Vorhersage (laengere Famine-
    Halbperiode = staerkere kumulative Wirkung), nicht bloss irgendein
    Unterschied.

    Nur Oszillationsbedingungen gehen ein; PosCtrl/NegCtrl sind der Bezug, kein
    Punkt auf der Perioden-Achse.
    """
    if endpoint.empty or period_col not in endpoint.columns:
        return pd.DataFrame()

    df = endpoint.copy()
    if "condition_type" in df.columns:
        df = df[df["condition_type"] == "Oscillation"]
    df = df.dropna(subset=[value_col])
    df["period_min"] = pd.to_numeric(df[period_col], errors="coerce")
    df = df.dropna(subset=["period_min"])
    if df.empty:
        logger.warning("test_aberrant_across_periods(): keine Oszillationskammern mit '%s'.", value_col)
        return pd.DataFrame()

    if group_cols is None:
        group_cols = [c for c in ["biosensor", "osc_type"] if c in df.columns]
    if not group_cols:
        df = df.assign(_all="all")
        group_cols = ["_all"]

    results = []
    for keys, grp in df.groupby(group_cols, dropna=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        record = dict(zip(group_cols, keys))
        record["n_chambers"] = int(len(grp))
        record["n_periods"] = int(grp["period_min"].nunique())

        batches = [g[value_col].to_numpy() for _, g in grp.groupby("period_min") if len(g) >= 2]
        pooled = np.concatenate(batches) if batches else np.array([])
        if len(batches) < 2 or (pooled.size and np.allclose(pooled, pooled[0])):
            record["kruskal_h"], record["kruskal_p"] = np.nan, np.nan
        else:
            h, p = stats.kruskal(*batches)
            record["kruskal_h"], record["kruskal_p"] = h, p

        if grp["period_min"].nunique() >= 3:
            rho, p_rho = stats.spearmanr(grp["period_min"], grp[value_col])
            record["spearman_rho"], record["spearman_p"] = rho, p_rho
        else:
            record["spearman_rho"], record["spearman_p"] = np.nan, np.nan

        results.append(record)

    out = pd.DataFrame(results)
    for _, r in out.iterrows():
        if pd.notna(r.get("spearman_p")) and r["spearman_p"] < 0.05:
            direction = "steigt" if r["spearman_rho"] > 0 else "faellt"
            logger.info(
                "Monotoner Trend: aberranter Anteil %s mit laengerer Periode "
                "(rho=%.2f, p=%.4g, %d Kammern).",
                direction, r["spearman_rho"], r["spearman_p"], r["n_chambers"],
            )
    return out


# ==============================================================================
# Plots
# ==============================================================================

def plot_morphospace(
    cells: pd.DataFrame,
    thresholds: MorphotypeThresholds,
    out_path: Path,
    freq_order: Optional[Sequence[str]] = None,
    sample_n: int = 4000,
    random_state: int = 0,
) -> None:
    """Wo im Morphospace liegen die Zellen jeder Periode, relativ zum Normalfenster?"""
    df = cells.dropna(subset=MORPHOLOGY_FEATURES)
    if df.empty or "osc_freq" not in df.columns:
        logger.warning("plot_morphospace(): keine Daten - Plot uebersprungen.")
        return

    present = df["osc_freq"].dropna().unique().tolist()
    order = ([f for f in (freq_order or []) if f in present]
             + [f for f in natural_freq_sort(present) if f not in (freq_order or [])])

    fig, axes = plt.subplots(1, len(order), figsize=(2.9 * len(order), 3.4),
                             squeeze=False, sharex=True, sharey=True)
    axes = axes[0]
    rng = np.random.default_rng(random_state)

    for ax, freq in zip(axes, order):
        sub = df[df["osc_freq"] == freq]
        if len(sub) > sample_n:
            sub = sub.iloc[rng.choice(len(sub), sample_n, replace=False)]
        ax.scatter(sub["eccentricity"], sub["solidity"], s=3, alpha=0.18,
                   color="#28414F", linewidth=0, rasterized=True)
        ax.axvline(thresholds.max_eccentricity, color="#7B5EA7", linestyle="--", linewidth=1)
        ax.axhline(thresholds.min_solidity, color="#E45756", linestyle="--", linewidth=1)
        ax.set_title(f"{freq} min", fontsize=10)
        ax.set_xlabel("Eccentricity")
        if ax is axes[0]:
            ax.set_ylabel("Solidity")

    fig.suptitle(
        "Morphospace by oscillation period — dashed lines bound the normal window\n"
        "(right of purple = elongated; below red = clustered)",
        y=1.06,
    )
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    logger.info("Plot gespeichert: %s", out_path.name)


def plot_aberrant_over_time(
    aggregated: pd.DataFrame,
    out_path: Path,
    freq_order: Optional[Sequence[str]] = None,
    x_col: str = "time_h",
) -> None:
    """Die Kernabbildung: waechst der aberrante Anteil ueber die Stunden, und
    haengt dieses Wachstum von der Periode ab?"""
    df = aggregated.dropna(subset=["frac_aberrant_mean"])
    if df.empty or x_col not in df.columns:
        logger.warning("plot_aberrant_over_time(): keine Daten - Plot uebersprungen.")
        return

    facet_col = "osc_type" if "osc_type" in df.columns else None
    facets = sorted(df[facet_col].dropna().unique()) if facet_col else [None]
    present = df["osc_freq"].dropna().unique().tolist()
    order = ([f for f in (freq_order or []) if f in present]
             + [f for f in natural_freq_sort(present) if f not in (freq_order or [])])
    palette = dict(zip(order, sns.color_palette("viridis", n_colors=len(order))))

    fig, axes = plt.subplots(1, len(facets), figsize=(5.4 * len(facets), 4.2),
                             squeeze=False, sharey=True)
    axes = axes[0]
    for ax, facet in zip(axes, facets):
        sub = df[df[facet_col] == facet] if facet_col else df
        for freq in order:
            line = sub[sub["osc_freq"] == freq].sort_values(x_col)
            if line.empty:
                continue
            ax.plot(line[x_col], line["frac_aberrant_mean"], color=palette[freq],
                    linewidth=1.8, label=f"{freq} min")
            ax.fill_between(
                line[x_col],
                line["frac_aberrant_mean"] - line["frac_aberrant_sem"],
                line["frac_aberrant_mean"] + line["frac_aberrant_sem"],
                color=palette[freq], alpha=0.15, linewidth=0,
            )
        ax.set_xlabel("Time [h]" if x_col == "time_h" else pretty_label(x_col))
        if ax is axes[0]:
            ax.set_ylabel("Aberrant cells (fraction)")
        ax.set_title(str(facet) if facet else "All data", fontsize=11, fontweight="bold")

    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, title="Oscillation period", loc="center left",
                   bbox_to_anchor=(1.0, 0.5), fontsize=8, borderaxespad=0.0)
    fig.suptitle(
        "Cumulative morphological impact over the time course\n"
        "(line = mean over replicate chambers; band = ± SEM)",
        y=1.04,
    )
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    logger.info("Plot gespeichert: %s", out_path.name)


def plot_aberrant_vs_period(
    endpoint: pd.DataFrame,
    out_path: Path,
    freq_order: Optional[Sequence[str]] = None,
    value_col: str = "frac_aberrant",
) -> None:
    """Dosis-Wirkung: Endzustand gegen die Periode, mit den Kontrollen als Bezugsbaender."""
    if endpoint.empty or "osc_freq" not in endpoint.columns:
        logger.warning("plot_aberrant_vs_period(): keine Daten - Plot uebersprungen.")
        return

    df = endpoint.copy()
    if "condition_type" not in df.columns:
        df["condition_type"] = "Oscillation"
    osc = df[df["condition_type"] == "Oscillation"]
    if osc.empty:
        logger.warning("plot_aberrant_vs_period(): keine Oszillationskammern - Plot uebersprungen.")
        return

    present = osc["osc_freq"].dropna().unique().tolist()
    order = ([f for f in (freq_order or []) if f in present]
             + [f for f in natural_freq_sort(present) if f not in (freq_order or [])])
    facet_col = "osc_type" if "osc_type" in df.columns else None
    facets = sorted(osc[facet_col].dropna().unique()) if facet_col else [None]

    fig, axes = plt.subplots(1, len(facets), figsize=(5.0 * len(facets), 4.4),
                             squeeze=False, sharey=True)
    axes = axes[0]
    rng = np.random.default_rng(0)
    ctrl_style = {"PosCtrl": ("#2E7D4F", "PosCtrl (constant feast)"),
                  "NegCtrl": ("#A8442F", "NegCtrl (constant starvation)")}

    for ax, facet in zip(axes, facets):
        sub = osc[osc[facet_col] == facet] if facet_col else osc
        ctrl_sub = df[(df[facet_col] == facet)] if facet_col else df

        for ctype, (color, label) in ctrl_style.items():
            vals = ctrl_sub.loc[ctrl_sub["condition_type"] == ctype, value_col].dropna()
            if vals.empty:
                continue
            ax.axhspan(vals.quantile(0.25), vals.quantile(0.75), color=color, alpha=0.12, linewidth=0)
            ax.axhline(vals.median(), color=color, linestyle="--", linewidth=1.2, label=label)

        for x, freq in enumerate(order):
            vals = sub.loc[sub["osc_freq"] == freq, value_col].dropna()
            if vals.empty:
                continue
            ax.scatter(np.full(len(vals), x) + rng.uniform(-.11, .11, len(vals)), vals,
                       s=34, alpha=.8, color="#28414F", edgecolor="white", linewidth=.4, zorder=3)
            ax.plot([x - .2, x + .2], [vals.median()] * 2, color="black", linewidth=1.9, zorder=4)

        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order)
        ax.set_xlabel("Oscillation period [min]")
        if ax is axes[0]:
            ax.set_ylabel("Aberrant cells at end of run (fraction)")
        ax.set_title(str(facet) if facet else "All data", fontsize=11, fontweight="bold")

    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="center left", bbox_to_anchor=(1.0, 0.5),
                   fontsize=8, borderaxespad=0.0)
    fig.suptitle(
        "Dose–response: cumulative morphological damage vs. switching period\n"
        "Longer period = longer famine half-cycle. Points = chambers; bar = median",
        y=1.05,
    )
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    logger.info("Plot gespeichert: %s", out_path.name)


def plot_morphotype_composition(
    aggregated: pd.DataFrame,
    out_path: Path,
    freq_order: Optional[Sequence[str]] = None,
    x_col: str = "time_h",
) -> None:
    """Welcher Typ waechst? Gestapelte Anteile ueber die Zeit, eine Spalte pro Periode."""
    df = aggregated.dropna(subset=["frac_yeast_like_mean"])
    if df.empty or x_col not in df.columns:
        return

    present = df["osc_freq"].dropna().unique().tolist()
    order = ([f for f in (freq_order or []) if f in present]
             + [f for f in natural_freq_sort(present) if f not in (freq_order or [])])
    facet_col = "osc_type" if "osc_type" in df.columns else None
    facets = sorted(df[facet_col].dropna().unique()) if facet_col else [None]

    fig, axes = plt.subplots(len(facets), len(order),
                             figsize=(2.6 * len(order), 2.9 * len(facets)),
                             squeeze=False, sharex=True, sharey=True)
    aberrant_types = [m for m in MORPHOTYPE_ORDER if m != "yeast_like"]

    for i, facet in enumerate(facets):
        for j, freq in enumerate(order):
            ax = axes[i][j]
            sub = df[df["osc_freq"] == freq]
            if facet_col:
                sub = sub[sub[facet_col] == facet]
            sub = sub.sort_values(x_col)
            if sub.empty:
                ax.set_visible(False)
                continue
            ax.stackplot(
                sub[x_col],
                *[sub[f"frac_{m}_mean"] for m in aberrant_types],
                colors=[MORPHOTYPE_COLORS[m] for m in aberrant_types],
                labels=[m.replace("_", " ") for m in aberrant_types],
            )
            ax.set_ylim(0, max(0.05, float(sub[[f"frac_{m}_mean" for m in aberrant_types]].sum(axis=1).max()) * 1.15))
            if i == 0:
                ax.set_title(f"{freq} min", fontsize=10)
            if j == 0:
                ax.set_ylabel(f"{facet}\nfraction" if facet else "Fraction", fontsize=9)
            if i == len(facets) - 1:
                ax.set_xlabel("Time [h]" if x_col == "time_h" else pretty_label(x_col), fontsize=9)

    handles, labels = axes[0][0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, title="Morphotype", loc="center left",
                   bbox_to_anchor=(1.0, 0.5), fontsize=8, borderaxespad=0.0)
    fig.suptitle("Which aberrant morphotype accumulates, and when", y=1.03)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    logger.info("Plot gespeichert: %s", out_path.name)
