"""
mother_trajectories.py
=======================
Detaillierte Darstellung STABIL GETRACKTER Mütter - ergänzt die aggregierten
Auswertungen (Panel A Violin, Fig. 3a Zeitreihe, µ-Punkt+Errorbar) um eine
Einzelzell-Ebene, auf der sich die Lineage-Klassifikation visuell gegen die
rohen Trajektorien (Fläche, Sensor-Ratio, ...) prüfen lässt.

Baut auf lineage.py auf:
    - identify_mothers() + select_stable_mothers() liefern die Auswahl
      "stabil getrackter" Mütter (lang UND lückenarm, siehe coverage-Spalte).
    - classify_mother_bud() liefert lineage_events, aus denen hier zwei
      Dinge gebaut werden, die vorher nicht existierten:
        1. plot_mother_trajectories(): Budding-Events als Marker direkt in
           die Trajektorien-Plots einzeichnen, statt Tabelle und Plot
           unverbunden nebeneinander zu haben.
        2. build_lineage_tree() / summarise_lineage_depth(): seit dem Fix in
           lineage.py (Entkopplung von Erkennung und Qualitätsfilter,
           established_min_frames statt globaler Tracklänge) verketten sich
           mother_cell_uid -> bud_cell_uid jetzt korrekt über mehrere
           Generationen. Das war vorher kaum sinnvoll nutzbar, weil ein Bud,
           der selbst zur Mutter heranwuchs, als unverbundene zweite Mutter
           geführt wurde (siehe lineage.py Modul-Docstring "FIX").
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

try:
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "mother_trajectories.py benötigt matplotlib (pip install matplotlib)."
    ) from exc

from lineage import identify_mothers, select_stable_mothers, LineageParams

logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


# ==============================================================================
# 1. Trajektorien-Plot mit Budding-Markern für stabile Mütter
# ==============================================================================

def plot_mother_trajectories(
    cells: pd.DataFrame,
    lineage_events: pd.DataFrame,
    mothers: pd.DataFrame,
    value_cols: Sequence[str],
    out_path: Path,
    min_coverage: float = 0.7,
    min_n_frames: Optional[int] = None,
    max_mothers: int = 24,
    mothers_per_page: int = 6,
    min_per_frame: Optional[float] = None,
) -> pd.DataFrame:
    """
    Plottet die Roh-Trajektorien (z.B. Fläche + Sensor-Ratio) EINZELNER,
    stabil getrackter Mütter über die Zeit, mit vertikalen Linien an jedem
    erkannten Budding-Event dieser Mutter (aus lineage_events). Damit lässt
    sich direkt sehen, ob ein erkanntes Event optisch mit einem Flächen-
    sprung/Signalwechsel zusammenfällt - eine Einzelzell-Version der
    µ_event-vs-µ_area-Prüfung aus growth_rate.py/area_growth.py.

    Auswahl der Mütter: select_stable_mothers(mothers, min_coverage,
    min_n_frames) - siehe lineage.py für die Coverage-Definition. Bei mehr
    Treffern als max_mothers werden die max_mothers mit der höchsten
    coverage genommen (kein Zufalls-Sample, damit die Auswahl reproduzierbar
    ist und bevorzugt die saubersten Tracks zeigt).

    Parameters
    ----------
    cells : voller Zelldatensatz (nach QC), braucht exp_id, cell_uid, frame
            + die Spalten in value_cols
    lineage_events : Ergebnis von classify_mother_bud() - braucht
            mother_cell_uid, budding_frame (siehe lineage.py)
    mothers : Ergebnis von identify_mothers() - braucht coverage-Spalte
    value_cols : welche Spalten geplottet werden (je eine Subplot-Zeile pro
            Mutter x Spalte), z.B. ["area", "ratio_iGlucoSnFR"]
    out_path : Pfad der Ausgabe-PDF (mehrseitig, mothers_per_page Mütter pro Seite)
    min_coverage, min_n_frames : siehe select_stable_mothers()
    max_mothers : harte Obergrenze, wie viele Mütter insgesamt geplottet werden
            (Default 24 - mehr wird unübersichtlich; ggf. cell_uids gezielt
            vorher in 'mothers' filtern statt max_mothers hochzusetzen)
    mothers_per_page : wie viele Mütter pro PDF-Seite (Default 6)
    min_per_frame : falls gesetzt, wird die x-Achse in Stunden statt Frames
            beschriftet (frame * min_per_frame / 60)

    Returns
    -------
    DataFrame der tatsächlich geplotteten Mütter (Teilmenge von 'mothers',
    inkl. coverage), zur Nachvollziehbarkeit was warum ausgewählt wurde.
    """
    missing = {"exp_id", "cell_uid", "frame"} - set(cells.columns)
    if missing:
        raise ValueError(f"plot_mother_trajectories() fehlen Spalten in 'cells': {missing}")
    missing_vals = [c for c in value_cols if c not in cells.columns]
    if missing_vals:
        raise ValueError(f"plot_mother_trajectories(): value_cols nicht in 'cells' gefunden: {missing_vals}")

    selected = select_stable_mothers(mothers, min_coverage=min_coverage, min_n_frames=min_n_frames)
    if selected.empty:
        logger.warning(
            "plot_mother_trajectories(): keine Mutter erfüllt min_coverage=%.2f%s - "
            "kein Plot erzeugt. Ggf. min_coverage senken oder Tracking-Qualität prüfen.",
            min_coverage, f"/min_n_frames={min_n_frames}" if min_n_frames else "",
        )
        return selected

    if len(selected) > max_mothers:
        logger.info(
            "plot_mother_trajectories(): %d stabile Mütter gefunden, zeige die %d mit "
            "höchster coverage (max_mothers=%d).", len(selected), max_mothers, max_mothers,
        )
        selected = selected.head(max_mothers)

    # events_by_mother: schneller Lookup mother_cell_uid -> Liste budding_frame
    has_events = not lineage_events.empty and "mother_cell_uid" in lineage_events.columns
    events_by_mother = (
        lineage_events.groupby("mother_cell_uid")["budding_frame"].apply(list)
        if has_events else {}
    )

    x_label = "Time [h]" if min_per_frame is not None else "Frame"

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_pages = int(np.ceil(len(selected) / mothers_per_page))

    with PdfPages(out_path) as pdf:
        for page in range(n_pages):
            page_rows = selected.iloc[page * mothers_per_page:(page + 1) * mothers_per_page]
            n_rows = len(page_rows)
            n_cols = len(value_cols)
            fig, axes = plt.subplots(
                n_rows, n_cols, figsize=(5 * n_cols, 2.4 * n_rows), squeeze=False,
            )

            for row_i, (_, mom_row) in enumerate(page_rows.iterrows()):
                exp_id, cell_uid = mom_row["exp_id"], mom_row["cell_uid"]
                traj = cells[(cells["exp_id"] == exp_id) & (cells["cell_uid"] == cell_uid)].sort_values("frame")
                bud_frames = events_by_mother.get(cell_uid, []) if has_events else []

                for col_i, value_col in enumerate(value_cols):
                    ax = axes[row_i][col_i]
                    x = traj["frame"] * min_per_frame / 60.0 if min_per_frame is not None else traj["frame"]
                    ax.plot(x, traj[value_col], color="#0F6E56", linewidth=1.2, marker="o", markersize=2)

                    for bf in bud_frames:
                        bx = bf * min_per_frame / 60.0 if min_per_frame is not None else bf
                        ax.axvline(bx, color="#D85A30", linewidth=1.0, linestyle="--", alpha=0.8)

                    if row_i == 0:
                        ax.set_title(value_col, fontsize=10)
                    if col_i == 0:
                        ax.set_ylabel(
                            f"{cell_uid}\ncoverage={mom_row['coverage']:.2f}", fontsize=8,
                        )
                    if row_i == n_rows - 1:
                        ax.set_xlabel(x_label, fontsize=9)
                    ax.tick_params(labelsize=8)

            fig.suptitle(
                f"Stably tracked mother cells (coverage ≥ {min_coverage:.2f}) - "
                f"page {page + 1}/{n_pages} - dashed red lines = detected budding events",
                fontsize=10,
            )
            fig.tight_layout(rect=(0, 0, 1, 0.96))
            pdf.savefig(fig)
            plt.close(fig)

    logger.info(
        "plot_mother_trajectories(): %d Mütter x %d Metriken auf %d Seiten geplottet -> %s",
        len(selected), len(value_cols), n_pages, out_path,
    )
    return selected


# ==============================================================================
# 1b. Eine stabile Mutter PRO GRUPPE (Strain x Oszillation x Condition), mit
#     nur dem für diese Mutter tatsächlich relevanten Ratio-Kanal
# ==============================================================================

def select_one_stable_mother_per_group(
    mothers: pd.DataFrame,
    cells: pd.DataFrame,
    group_cols: Sequence[str] = ("biosensor", "osc_type", "osc_freq", "condition"),
    min_coverage: float = 0.7,
    min_n_frames: Optional[int] = None,
    n_per_group: int = 1,
) -> pd.DataFrame:
    """
    Wählt statt eines globalen Top-N (plot_mother_trajectories()) für JEDE
    Kombination aus group_cols (Default: Biosensor/Strain x Oszillationstyp x
    -frequenz x Condition) die n_per_group stabilste(n) Mutter(zellen) aus
    select_stable_mothers(). So bekommt ihr eine repräsentative Mutter pro
    Bedingung statt z.B. 24 Mütter, die zufällig alle aus derselben, am
    saubersten getrackten Bedingung stammen.

    group_cols müssen in 'cells' vorhanden sein (fehlende werden mit Warnung
    ignoriert) - die Metadaten werden aus 'cells' pro exp_id/cell_uid
    ermittelt (eine Zelle gehört immer zu genau einer Bedingung, siehe
    exp_id-Struktur in run_analysis.py).

    Returns
    -------
    Teilmenge von identify_mothers()-Output, ergänzt um group_cols, sortiert
    nach group_cols dann absteigend nach coverage. Gruppen, für die KEINE
    Mutter min_coverage/min_n_frames erfüllt, fehlen in der Ausgabe - siehe
    Log-Warnung dafür (n_missing_groups), damit das nicht stillschweigend
    passiert.
    """
    group_cols = [c for c in group_cols if c in cells.columns]
    missing_group_cols = set(group_cols) - set(cells.columns)
    if missing_group_cols:
        logger.warning(
            "select_one_stable_mother_per_group(): group_cols nicht in 'cells' "
            "gefunden, werden ignoriert: %s", missing_group_cols,
        )
    if not group_cols:
        raise ValueError("select_one_stable_mother_per_group(): keine gültigen group_cols übrig.")

    candidates = select_stable_mothers(mothers, min_coverage=min_coverage, min_n_frames=min_n_frames)
    if candidates.empty:
        return candidates

    meta = cells[["exp_id", "cell_uid"] + group_cols].drop_duplicates(subset=["exp_id", "cell_uid"])
    merged = candidates.merge(meta, on=["exp_id", "cell_uid"], how="left")

    n_unmatched = int(merged[group_cols].isna().any(axis=1).sum())
    if n_unmatched > 0:
        logger.warning(
            "select_one_stable_mother_per_group(): %d Mütter ohne vollständige "
            "Metadaten (%s) - werden ausgeschlossen.", n_unmatched, group_cols,
        )
    merged = merged.dropna(subset=group_cols)

    selected = (
        merged.sort_values("coverage", ascending=False)
        .groupby(group_cols, dropna=False, as_index=False)
        .head(n_per_group)
        .sort_values(group_cols)
        .reset_index(drop=True)
    )

    all_combos = cells[group_cols].dropna().drop_duplicates()
    covered_combos = selected[group_cols].drop_duplicates()
    n_missing = len(all_combos) - len(
        all_combos.merge(covered_combos, on=group_cols, how="inner")
    )
    if n_missing > 0:
        logger.warning(
            "select_one_stable_mother_per_group(): %d von %d Bedingungs-Kombinationen "
            "(%s) haben KEINE Mutter mit coverage >= %.2f%s - fehlen im Plot. "
            "Ggf. min_coverage senken oder Tracking dieser Bedingungen prüfen.",
            n_missing, len(all_combos), group_cols, min_coverage,
            f"/min_n_frames={min_n_frames}" if min_n_frames else "",
        )
    logger.info(
        "select_one_stable_mother_per_group(): %d Mütter für %d von %d Bedingungs-"
        "Kombinationen ausgewählt (n_per_group=%d).",
        len(selected), len(covered_combos), len(all_combos), n_per_group,
    )
    return selected


def _own_ratio_columns(cells_sub: pd.DataFrame, ratio_prefix: str = "ratio_") -> list[str]:
    """
    Ermittelt, welche ratio_*-Spalten für DIESE Zelle tatsächlich befüllt
    sind (nicht durchgängig NaN) - robust gegenüber der genauen Namens-
    konvention aus sensors.compute_ratios(). Ein Biosensor hat i.d.R. genau
    einen befüllten Ratio-Kanal, alle anderen ratio_*-Spalten sind für seine
    Zellen komplett NaN (anderer Kanalsatz).
    """
    ratio_cols = [c for c in cells_sub.columns if c.startswith(ratio_prefix)]
    return [c for c in ratio_cols if cells_sub[c].notna().any()]


def plot_stable_mother_per_group(
    cells: pd.DataFrame,
    lineage_events: pd.DataFrame,
    mothers: pd.DataFrame,
    out_path: Path,
    group_cols: Sequence[str] = ("biosensor", "osc_type", "osc_freq", "condition"),
    base_value_cols: Sequence[str] = ("area",),
    ratio_prefix: str = "ratio_",
    min_coverage: float = 0.7,
    min_n_frames: Optional[int] = None,
    n_per_group: int = 1,
    mothers_per_page: int = 6,
    min_per_frame: Optional[float] = None,
) -> pd.DataFrame:
    """
    Wie plot_mother_trajectories(), aber zwei Unterschiede:
      1. Auswahl über select_one_stable_mother_per_group() statt globalem
         Top-N - eine stabile Mutter PRO Bedingungs-Kombination (Strain x
         Oszillationstyp x -frequenz x Condition), nicht die insgesamt
         stabilsten Mütter unabhängig von der Bedingung.
      2. Pro Mutter wird nur ihr EIGENER Ratio-Kanal geplottet (ermittelt aus
         den tatsächlich befüllten ratio_*-Spalten dieser Zelle), statt für
         jede Mutter alle ratio_*-Spalten aus dem gesamten Datensatz
         durchzuplotten - vermeidet leere Panels für Sensor-Kanäle, die zum
         Biosensor dieser Mutter gar nicht gehören.

    base_value_cols gilt weiterhin für ALLE Mütter identisch (z.B. "area").
    Falls einzelne Mütter keinen befüllten Ratio-Kanal haben (sollte nicht
    vorkommen, aber Datenlücken sind möglich), bleibt das entsprechende Panel
    leer und wird als solches beschriftet, statt die Spaltenzahl der Seite zu
    verzerren.

    Returns
    -------
    Ausgewählte Mütter (Ergebnis von select_one_stable_mother_per_group()),
    zur Nachvollziehbarkeit welche Mutter welche Bedingung repräsentiert.
    """
    selected = select_one_stable_mother_per_group(
        mothers, cells, group_cols=group_cols, min_coverage=min_coverage,
        min_n_frames=min_n_frames, n_per_group=n_per_group,
    )
    if selected.empty:
        logger.warning("plot_stable_mother_per_group(): keine Mutter ausgewählt - kein Plot erzeugt.")
        return selected

    group_cols = [c for c in group_cols if c in selected.columns]
    base_value_cols = list(base_value_cols)

    has_events = not lineage_events.empty and "mother_cell_uid" in lineage_events.columns
    events_by_mother = (
        lineage_events.groupby("mother_cell_uid")["budding_frame"].apply(list)
        if has_events else {}
    )

    x_label = "Time [h]" if min_per_frame is not None else "Frame"

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_pages = int(np.ceil(len(selected) / mothers_per_page))

    with PdfPages(out_path) as pdf:
        for page in range(n_pages):
            page_rows = selected.iloc[page * mothers_per_page:(page + 1) * mothers_per_page]
            n_rows = len(page_rows)

            # Pro Zeile (Mutter) eigene value_cols ermitteln (base + ihr
            # eigener Ratio-Kanal), damit KEINE unpassenden Sensor-Spalten
            # anderer Biosensoren mitgeplottet werden.
            row_trajs, row_value_cols = [], []
            for _, mom_row in page_rows.iterrows():
                traj = cells[
                    (cells["exp_id"] == mom_row["exp_id"]) & (cells["cell_uid"] == mom_row["cell_uid"])
                ].sort_values("frame")
                own_ratio = _own_ratio_columns(traj, ratio_prefix)
                row_trajs.append(traj)
                row_value_cols.append(base_value_cols + own_ratio)

            n_cols = max(len(v) for v in row_value_cols)
            fig, axes = plt.subplots(
                n_rows, n_cols, figsize=(5 * n_cols, 2.4 * n_rows), squeeze=False,
            )

            for row_i, (_, mom_row) in enumerate(page_rows.iterrows()):
                traj = row_trajs[row_i]
                vcols = row_value_cols[row_i]
                bud_frames = events_by_mother.get(mom_row["cell_uid"], []) if has_events else []
                x = traj["frame"] * min_per_frame / 60.0 if min_per_frame is not None else traj["frame"]

                group_label = " / ".join(str(mom_row[c]) for c in group_cols)

                for col_i in range(n_cols):
                    ax = axes[row_i][col_i]
                    if col_i >= len(vcols):
                        ax.axis("off")
                        continue

                    value_col = vcols[col_i]
                    ax.plot(x, traj[value_col], color="#0F6E56", linewidth=1.2, marker="o", markersize=2)
                    for bf in bud_frames:
                        bx = bf * min_per_frame / 60.0 if min_per_frame is not None else bf
                        ax.axvline(bx, color="#D85A30", linewidth=1.0, linestyle="--", alpha=0.8)

                    ax.set_title(value_col, fontsize=9)
                    if col_i == 0:
                        ax.set_ylabel(
                            f"{group_label}\n{mom_row['cell_uid']} (cov={mom_row['coverage']:.2f})",
                            fontsize=7.5,
                        )
                    if row_i == n_rows - 1:
                        ax.set_xlabel(x_label, fontsize=9)
                    ax.tick_params(labelsize=8)

            fig.suptitle(
                f"One stably tracked mother per condition ({' x '.join(group_cols)}) - "
                f"page {page + 1}/{n_pages} - dashed red lines = budding events",
                fontsize=9.5,
            )
            fig.tight_layout(rect=(0, 0, 1, 0.95))
            pdf.savefig(fig)
            plt.close(fig)

    logger.info(
        "plot_stable_mother_per_group(): %d Mütter (eine pro Bedingung) auf %d Seiten geplottet -> %s",
        len(selected), n_pages, out_path,
    )
    return selected


# ==============================================================================
# 2. Mehrgenerationaler Stammbaum aus lineage_events
# ==============================================================================

def build_lineage_tree(lineage_events: pd.DataFrame) -> pd.DataFrame:
    """
    Rekonstruiert pro Kammer (exp_id) den Stammbaum aus lineage_events, indem
    mother_cell_uid -> bud_cell_uid Kanten zu Ketten verknüpft werden.

    Erst seit dem Entkopplungs-Fix in lineage.py (classify_mother_bud() filtert
    Bud-Kandidaten nicht mehr über ihre finale Tracklänge, sondern prüft
    Mutter-Kandidaten kausal über established_min_frames) sind diese Ketten
    verlässlich mehrgenerational: ein Bud, der selbst zur Mutter heranwächst,
    bleibt über sein eigenes mother_cell_uid mit der vorigen Generation
    verbunden, statt als unabhängige, unverbundene Mutter zu erscheinen.

    "Wurzel"-Zellen (generation=0) sind Zellen, die selbst NIE als bud_cell_uid
    in lineage_events auftauchen - ihre tatsächliche Herkunft liegt entweder
    vor Beobachtungsbeginn der Kammer oder ihr Budding-Event wurde nicht
    erkannt (siehe classify_mother_bud()-Logging für Erkennungsraten).

    Returns
    -------
    DataFrame mit einer Zeile pro Zelle, die in lineage_events als Mutter
    ODER Bud auftaucht:
        exp_id, cell_uid, parent_cell_uid (NaN für Wurzeln), generation,
        budding_frame (Frame, an dem DIESE Zelle als Bud entstand; NaN für
        Wurzeln), n_own_buds (wie viele eigene Budding-Events diese Zelle hat)
    """
    required = {"exp_id", "mother_cell_uid", "bud_cell_uid", "budding_frame"}
    missing = required - set(lineage_events.columns)
    if missing:
        raise ValueError(
            f"build_lineage_tree() fehlen Spalten in 'lineage_events': {missing}. "
            f"Braucht cell_uid-basierte Events aus classify_mother_bud() "
            f"(Datensatz mit 'cell_uid'-Spalte)."
        )

    if lineage_events.empty:
        logger.warning("build_lineage_tree(): lineage_events ist leer - leerer Stammbaum.")
        return pd.DataFrame()

    n_own_buds = lineage_events.groupby(["exp_id", "mother_cell_uid"]).size()

    results = []
    for exp_id, group in lineage_events.groupby("exp_id"):
        # parent_of[cell_uid] = (mother_cell_uid, budding_frame) - jede Zelle
        # kann in EINER Kammer nur aus EINEM Event entstanden sein
        parent_of = {}
        for _, ev in group.iterrows():
            bud = ev["bud_cell_uid"]
            if bud in parent_of:
                logger.warning(
                    "build_lineage_tree(): '%s' in Kammer '%s' hat mehr als ein "
                    "Mutter-Event - erstes wird verwendet, weitere ignoriert "
                    "(möglicher Tracking-Konflikt, mit inspect_classification() prüfen).",
                    bud, exp_id,
                )
                continue
            parent_of[bud] = (ev["mother_cell_uid"], ev["budding_frame"])

        # sorted(): die Iterationsreihenfolge eines Python-set haengt vom
        # PYTHONHASHSEED ab und wechselt damit zwischen zwei Prozessen. Ohne
        # diese Sortierung hat 23_lineage_tree.csv bei jedem Lauf eine andere
        # Zeilenreihenfolge - identischer Inhalt, aber nicht reproduzierbar
        # und bei jedem diff scheinbar veraendert.
        all_cells = sorted(set(group["mother_cell_uid"]).union(group["bud_cell_uid"]))

        # Generation via BFS von den Wurzeln aus (Zellen ohne Eintrag in parent_of)
        generation = {}

        def resolve_generation(cell_uid: str, _visiting: Optional[set] = None) -> int:
            if cell_uid in generation:
                return generation[cell_uid]
            if cell_uid not in parent_of:
                generation[cell_uid] = 0
                return 0
            _visiting = _visiting or set()
            if cell_uid in _visiting:
                logger.warning(
                    "build_lineage_tree(): Zyklus in Kammer '%s' bei '%s' erkannt - "
                    "als Wurzel behandelt (sollte bei echten Daten nicht vorkommen).",
                    exp_id, cell_uid,
                )
                generation[cell_uid] = 0
                return 0
            _visiting.add(cell_uid)
            parent_uid, _ = parent_of[cell_uid]
            generation[cell_uid] = resolve_generation(parent_uid, _visiting) + 1
            return generation[cell_uid]

        for cell_uid in all_cells:
            resolve_generation(cell_uid)

        for cell_uid in all_cells:
            parent_uid, budding_frame = parent_of.get(cell_uid, (np.nan, np.nan))
            results.append({
                "exp_id": exp_id,
                "cell_uid": cell_uid,
                "parent_cell_uid": parent_uid,
                "generation": generation[cell_uid],
                "budding_frame": budding_frame,
                "n_own_buds": int(n_own_buds.get((exp_id, cell_uid), 0)),
            })

    out = pd.DataFrame(results)
    logger.info(
        "build_lineage_tree(): %d Zellen über %d Kammern rekonstruiert, "
        "maximale Generationstiefe=%d.",
        len(out), out["exp_id"].nunique(), int(out["generation"].max()) if not out.empty else 0,
    )
    return out


def summarise_lineage_depth(tree: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregiert build_lineage_tree() pro Kammer zu einer Kennzahl, wie viele
    Generationen tatsächlich zuverlässig verfolgt wurden - eine Kennzahl, die
    es vor dem Entkopplungs-Fix in lineage.py praktisch nicht gab, weil
    Stammbäume vorher kaum über eine Generation hinausreichten.

    Returns
    -------
    DataFrame mit einer Zeile pro exp_id:
        exp_id, n_cells_in_tree, n_root_mothers, max_generation,
        n_cells_generation_2plus (wie viele Zellen sind Enkel-oder-weiter,
        also generation >= 2 - ein direktes Maß für "geht der Stammbaum
        über eine reine Mutter->Bud-Kante hinaus")
    """
    if tree.empty:
        return pd.DataFrame()

    def _agg(group: pd.DataFrame) -> pd.Series:
        return pd.Series({
            "n_cells_in_tree": len(group),
            "n_root_mothers": int(group["parent_cell_uid"].isna().sum()),
            "max_generation": int(group["generation"].max()),
            "n_cells_generation_2plus": int((group["generation"] >= 2).sum()),
        })

    out = tree.groupby("exp_id").apply(_agg, include_groups=False).reset_index()
    logger.info(
        "summarise_lineage_depth(): %d Kammern, mittlere max. Generationstiefe=%.1f.",
        len(out), out["max_generation"].mean(),
    )
    return out
