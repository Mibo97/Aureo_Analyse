"""
relink.py
=========
Zwei Werkzeuge gegen die Track-Fragmentierung, kalibriert am manuellen QC
von WT/pH/6 (11 Kammern, 9 063 Tracks, 503 QC-Zeilen):

1. gap_close_tracks(): automatisches Schliessen von Tracking-Luecken
   (Gap Closing). Ein Track, der bei Frame t neu beginnt, wird an einen Track
   angehaengt, der hoechstens max_gap Frames vorher GEENDET hat, wenn beide
   nah beieinander liegen (<= max_distance_px), aehnlich gross sind
   (Flaechenverhaeltnis in [1/r, r]) und die Zuordnung in BEIDE Richtungen
   eindeutig ist (kein zweiter Kandidat im Radius, weder ein zweiter beendeter
   noch ein zweiter neuer Track). Mehrdeutige Faelle werden NICHT verknuepft:
   lieber ein Bruch zu viel als zwei Zellen vermischt. Manuelle Merges aus
   qc_exclusions.csv laufen VORHER und haben Vorrang.

2. detect_sparse_window(): pro Kammer das Frame-Fenster, in dem das Bildfeld
   noch duenn besetzt ist (<= max_objects Objekte pro Frame, geglaettet). Nur
   dort ist die Mutter/Bud-Heuristik auswertbar. Der Tracker vergibt in diesen
   Daten mit ~12 % pro Objekt und Frame eine neue ID, unabhaengig von der
   Dichte; jede neu auftauchende Zelle ist also zunaechst ein Fragment, keine
   Knospe. Im duennen Feld (Abstaende von Hunderten Pixeln) laesst sich das
   per Gap Closing eindeutig korrigieren, im vollen Feld (50-180 Objekte am
   Ende der Laeufe, Spruenge derselben Zelle von median 83 px) nicht mehr:
   dort sind die "Budding-Events" Fragment-Statistik (im QC-Batch 462 pro
   Kammer, davon 243 in den letzten 22 Frames, und das manuelle QC entfernte
   2.6 % davon).

   track_fragmentation() liefert die Kennzahlen dazu pro Kammer - vor und
   nach dem Gap Closing - fuer 00_track_fragmentation.csv.

AUSGABEN (run_analysis.py)
--------------------------
    00_track_fragmentation.csv  Kennzahlen je Kammer und Stufe (nach QC / nach Gap Closing)
    00_new_objects_vs_density.csv/.pdf  neue Tracks gegen die Objektdichte, je Kammer und Block
                                (new_objects_vs_density(), plot_new_objects_vs_density())
    00_track_relinks.csv        eine Zeile pro automatischer Verknuepfung
    20_lineage_window.csv       Fenster je Kammer (Ende, Laenge, ob auswertbar)
    20_lineage_window.pdf       Objekte pro Frame ueber die Zeit, Fensterlaengen
"""

from __future__ import annotations

import logging
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

REQUIRED = {"exp_id", "cell_uid", "track_id", "frame", "centroid_x", "centroid_y", "area"}


def _check(cells: pd.DataFrame, who: str) -> None:
    missing = REQUIRED - set(cells.columns)
    if missing:
        raise ValueError(f"{who} fehlen Spalten: {sorted(missing)}")


# ==============================================================================
# 1. Fragmentierung messen
# ==============================================================================

def track_fragmentation(cells: pd.DataFrame, long_track_frames: int = 10, stage: str = "") -> pd.DataFrame:
    """Kennzahlen der Track-Fragmentierung je Kammer.

    new_tracks_per_object_frame: neu beginnende Tracks je Objekt und Frame
    (ohne den ersten Frame der Kammer) - die Rate, mit der der Tracker eine
    Zelle "verliert". share_objectframes_in_long_tracks: Anteil der
    Objekt-Frames, die zu Tracks mit >= long_track_frames Frames gehoeren.
    """
    _check(cells, "track_fragmentation()")
    rows = []
    for exp_id, g in cells.groupby("exp_id"):
        per_frame = g.groupby("frame").size().sort_index()
        start = per_frame.index.min()
        nfr = g.groupby("cell_uid")["frame"].nunique()
        t0 = g.groupby("cell_uid")["frame"].min()
        n_new = int((t0 > start).sum())
        obj_frames_after_start = int(per_frame[per_frame.index > start].sum())
        long_uids = nfr.index[nfr >= long_track_frames]
        rows.append({
            "exp_id": exp_id, "stage": stage,
            "n_frames": int(len(per_frame)),
            "objects_first_frame": int(per_frame.iloc[0]),
            "objects_last_frame": int(per_frame.iloc[-1]),
            "objects_per_frame_median": float(per_frame.median()),
            "n_tracks": int(len(nfr)),
            "median_track_frames": float(nfr.median()),
            "share_single_frame_tracks": float((nfr == 1).mean()),
            "new_tracks_per_object_frame": (n_new / obj_frames_after_start) if obj_frames_after_start else np.nan,
            "share_objectframes_in_long_tracks": float(g["cell_uid"].isin(long_uids).mean()),
        })
    return pd.DataFrame(rows)


def new_objects_vs_density(cells: pd.DataFrame, block_frames: int = 22) -> pd.DataFrame:
    """Neu beginnende Tracks gegen die Objektdichte, je Kammer und Block von block_frames
    Frames (die Block-Tabelle von docs/data_story.md 2.2, fuer alle Kammern).

    Je Block: objects_per_frame (Median der Objekte je Frame), new_tracks_total (Tracks, die
    im Block beginnen; der erste Frame der Kammer zaehlt nicht, dort beginnt jeder Track),
    davon new_tracks_touching (die Maske beruehrt beim ersten Auftreten eine getrackte Maske
    oder faellt aus einer heraus: link_type 'new_touching' / 'split') und new_tracks_free
    (link_type 'new', beruehrt nichts). Dazu dieselben Zahlen je Frame nach dem Kammerstart
    (*_per_frame) und je Objekt und Frame (*_per_object_frame). Ohne link_type
    (v11-Tabellen) bleiben die getrennten Spalten NaN.

    Im dichten Feld beruehrt jeder neue Track eine Maske: Fragmente aus Beruehren, Teilen und
    Wiedervereinigen, keine Knospen - der Grund fuer das Sparse-Phase-Fenster
    (detect_sparse_window()).
    """
    _check(cells, "new_objects_vs_density()")
    if cells.empty:
        return pd.DataFrame()
    has_lt = "link_type" in cells.columns
    cols = ["exp_id", "cell_uid", "frame"] + (["link_type"] if has_lt else [])
    df = cells[cols].copy()
    start_of = df.groupby("exp_id")["frame"].min()
    df["start"] = df["exp_id"].map(start_of)
    df["block"] = ((df["frame"] - df["start"]) // block_frames).astype(int)

    per_frame = df.groupby(["exp_id", "block", "frame"]).size().rename("n").reset_index()
    per_frame["after_start"] = per_frame["frame"] > per_frame["exp_id"].map(start_of)
    keys = ["exp_id", "block"]
    dens = per_frame.groupby(keys).agg(
        frame_lo=("frame", "min"), frame_hi=("frame", "max"), n_frames=("frame", "size"),
        objects_per_frame=("n", "median"),
    ).reset_index()
    after = per_frame[per_frame["after_start"]].groupby(keys).agg(
        frames_after_start=("frame", "size"), object_frames_after_start=("n", "sum"),
    ).reset_index()
    dens = dens.merge(after, on=keys, how="left")
    dens[["frames_after_start", "object_frames_after_start"]] = (
        dens[["frames_after_start", "object_frames_after_start"]].fillna(0).astype(int))

    agg = {"exp_id": ("exp_id", "first"), "first_frame": ("frame", "first"), "start": ("start", "first")}
    if has_lt:
        agg["link_type"] = ("link_type", "first")
    first = df.sort_values("frame").groupby("cell_uid", sort=False).agg(**agg).reset_index()
    first = first[first["first_frame"] > first["start"]].copy()
    first["block"] = ((first["first_frame"] - first["start"]) // block_frames).astype(int)
    out = dens.merge(first.groupby(keys).size().rename("new_tracks_total").reset_index(), on=keys, how="left")
    out["new_tracks_total"] = out["new_tracks_total"].fillna(0).astype(int)
    if has_lt:
        lt = first["link_type"].astype(str)
        touching = first[lt.isin(["new_touching", "split"])].groupby(keys).size().rename("new_tracks_touching")
        free = first[lt == "new"].groupby(keys).size().rename("new_tracks_free")
        out = (out.merge(touching.reset_index(), on=keys, how="left")
                  .merge(free.reset_index(), on=keys, how="left"))
        out["new_tracks_touching"] = out["new_tracks_touching"].fillna(0).astype(int)
        out["new_tracks_free"] = out["new_tracks_free"].fillna(0).astype(int)
    else:
        out["new_tracks_touching"] = np.nan
        out["new_tracks_free"] = np.nan
    frames = out["frames_after_start"].where(out["frames_after_start"] > 0)
    obj_frames = out["object_frames_after_start"].where(out["object_frames_after_start"] > 0)
    for kind in ("total", "touching", "free"):
        out[f"new_tracks_{kind}_per_frame"] = out[f"new_tracks_{kind}"] / frames
        out[f"new_tracks_{kind}_per_object_frame"] = out[f"new_tracks_{kind}"] / obj_frames
    meta = [c for c in ["biosensor", "osc_type", "osc_freq", "condition", "chip", "chip_family", "medium"]
            if c in cells.columns]
    if meta:
        out = out.merge(cells[["exp_id"] + meta].drop_duplicates("exp_id"), on="exp_id", how="left")
    return out


def plot_new_objects_vs_density(table: pd.DataFrame, out_path, max_objects: int, block_frames: int) -> None:
    """Links: neue Tracks je Frame gegen die Objekte je Frame, ein Punkt je Kammer und Block -
    mit beruehrender Maske dunkel, ohne hell; Linie und Band = Median und Quartile je
    Dichteklasse (Zweierpotenzen, ab 5 Bloecken). Rechts: neue beruehrende Tracks je Objekt
    und Frame gegen die Dichte, also was der Tracker je Objekt an Fragmenten erzeugt.
    Senkrechte Linie = Sparse-Grenze (max_objects)."""
    from plot_style import INK, INK_MUTED, SURFACE, dedupe_handles, finish, legend_below, panel_title
    if table is None or table.empty or "objects_per_frame" not in table.columns:
        return
    t = table[table["objects_per_frame"] > 0].copy()
    if t.empty:
        return
    has_lt = t["new_tracks_touching"].notna().any()
    left = ([("new_tracks_touching_per_frame", "new track touching a tracked mask, or split from one", INK, True),
             ("new_tracks_free_per_frame", "new track touching nothing", INK_MUTED, False)] if has_lt
            else [("new_tracks_total_per_frame", "new tracks (tables without link_type)", INK, True)])
    right_col = "new_tracks_touching_per_object_frame" if has_lt else "new_tracks_total_per_object_frame"
    n_edges = int(np.ceil(np.log2(max(float(t["objects_per_frame"].max()), 2.0)))) + 2
    edges = 2.0 ** np.arange(0, n_edges)
    t["dbin"] = pd.cut(t["objects_per_frame"], bins=edges, right=False)
    rng = np.random.default_rng(0)

    def binned(ax, col, color, filled, label):
        ok = t[col].notna()
        x_jit = t.loc[ok, "objects_per_frame"] * 2.0 ** rng.uniform(-0.06, 0.06, int(ok.sum()))
        ax.scatter(x_jit, t.loc[ok, col], marker="o", s=5, facecolor=color if filled else SURFACE,
                   edgecolor=color, linewidth=0.35, alpha=0.2, zorder=2)
        g = t.loc[ok].groupby("dbin", observed=True)[col]
        med, q25, q75, n = g.median(), g.quantile(0.25), g.quantile(0.75), g.size()
        keep = (n >= 5).to_numpy()
        if not keep.any():
            return
        x = np.array([np.sqrt(iv.left * iv.right) for iv in med.index[keep]])
        ax.fill_between(x, q25.to_numpy()[keep], q75.to_numpy()[keep], color=color, alpha=0.15,
                        linewidth=0, zorder=3)
        ax.plot(x, med.to_numpy()[keep], color=color, linewidth=1.4, marker="o", markersize=4,
                markerfacecolor=color if filled else SURFACE, markeredgecolor=color, zorder=4, label=label)

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.0))
    for col, label, color, filled in left:
        binned(ax, col, color, filled, label)
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=0.1, linscale=0.5)
    ax.set_ylim(bottom=0)
    ax.axvline(max_objects, color=INK, linestyle="--", linewidth=1.0, label=f"sparse limit ({max_objects} objects)")
    ax.set_xlabel("objects per frame (median of the block)")
    ax.set_ylabel("new tracks per frame")
    panel_title(ax, f"{t['exp_id'].nunique()} chambers, blocks of {block_frames} frames")

    binned(ax2, right_col, INK, True, None)
    ax2.set_xscale("log")
    ax2.set_yscale("symlog", linthresh=0.01, linscale=0.5)
    ax2.set_ylim(bottom=0)
    ax2.axvline(max_objects, color=INK, linestyle="--", linewidth=1.0)
    ax2.set_xlabel("objects per frame (median of the block)")
    ax2.set_ylabel("new touching tracks\nper object and frame" if has_lt else "new tracks\nper object and frame")
    panel_title(ax2, "new touching tracks per object and frame")
    handles, labels = dedupe_handles([ax])
    legend_below(fig, handles, labels, ncol=2 if len(handles) > 2 else len(handles), y=0.0)
    finish(fig, out_path, logger)


# ==============================================================================
# 2. Gap Closing
# ==============================================================================

def gap_close_tracks(
    cells: pd.DataFrame,
    max_gap: int = 2,
    max_distance_px: float = 150.0,
    max_area_ratio: float = 2.0,
    window_col: Optional[str] = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Eindeutige Tracking-Luecken schliessen (siehe Modul-Docstring).

    window_col : Name einer booleschen Spalte (z.B. in_lineage_window). Falls
        gesetzt, werden nur neue Tracks verknuepft, deren ERSTE Zeile dort True
        ist - Gap Closing bleibt damit auf das duenn besetzte Feld beschraenkt,
        wo eine eindeutige Zuordnung auch eine richtige ist. Im vollen Feld
        bleiben die Tracks unangetastet (dort laeuft auch keine Lineage).

    Returns
    -------
    cells  : Kopie mit umbenannten track_id/cell_uid fuer alle verknuepften Tracks
             (der frueheste Track der Kette gibt den Namen)
    links  : eine Zeile pro Verknuepfung: exp_id, new_cell_uid, ended_cell_uid,
             root_cell_uid, frame, gap_frames, distance_px, area_ratio,
             objects_in_frame
    stats  : je Kammer: n_new_tracks, n_linked, n_ambiguous, n_no_candidate
    """
    _check(cells, "gap_close_tracks()")
    if max_area_ratio < 1:
        raise ValueError("max_area_ratio muss >= 1 sein (Verhaeltnis-Fenster [1/r, r]).")
    r_lo, r_hi = 1.0 / max_area_ratio, max_area_ratio

    cells = cells.sort_values(["exp_id", "cell_uid", "frame"], kind="stable")
    first = cells.groupby("cell_uid").head(1).set_index("cell_uid")
    last = cells.groupby("cell_uid").tail(1).set_index("cell_uid")
    exp_of_uid = first["exp_id"]

    link_rows, stat_rows = [], []
    parent: dict[str, str] = {}

    def find(u: str) -> str:
        parent.setdefault(u, u)
        while parent[u] != u:
            parent[u] = parent[parent[u]]
            u = parent[u]
        return u

    for exp_id, uids in exp_of_uid.groupby(exp_of_uid).groups.items():
        uids = pd.Index(uids)
        f = first.loc[uids]
        l = last.loc[uids]
        start = int(f["frame"].min())
        objects_in_frame = cells.loc[cells["exp_id"] == exp_id].groupby("frame").size()
        f_frame = f["frame"].to_numpy(); f_x = f["centroid_x"].to_numpy(float); f_y = f["centroid_y"].to_numpy(float); f_a = f["area"].to_numpy(float)
        f_ok = f[window_col].to_numpy(bool) if window_col else np.ones(len(f), dtype=bool)
        l_frame = l["frame"].to_numpy(); l_x = l["centroid_x"].to_numpy(float); l_y = l["centroid_y"].to_numpy(float); l_a = l["area"].to_numpy(float)
        uid_arr = uids.to_numpy()

        n_new = n_linked = n_ambiguous = n_none = 0
        order = np.argsort(f_frame, kind="stable")
        for i in order:
            t = int(f_frame[i])
            if t <= start or not f_ok[i]:
                continue
            n_new += 1
            # beendete Tracks kurz vor t
            ended = (l_frame >= t - max_gap) & (l_frame <= t - 1)
            if not ended.any():
                n_none += 1
                continue
            d = np.hypot(l_x[ended] - f_x[i], l_y[ended] - f_y[i])
            ratio = f_a[i] / np.where(l_a[ended] > 0, l_a[ended], np.nan)
            ok = (d <= max_distance_px) & (ratio >= r_lo) & (ratio <= r_hi)
            if ok.sum() == 0:
                n_none += 1
                continue
            if ok.sum() > 1:
                n_ambiguous += 1
                continue
            j = np.flatnonzero(ended)[np.flatnonzero(ok)[0]]
            # Gegenrichtung: ist der neue Track der EINZIGE, der zu diesem Ende passt?
            t_end = int(l_frame[j])
            starters = (f_frame > t_end) & (f_frame <= t_end + max_gap)
            d_rev = np.hypot(f_x[starters] - l_x[j], f_y[starters] - l_y[j])
            ratio_rev = f_a[starters] / (l_a[j] if l_a[j] > 0 else np.nan)
            ok_rev = (d_rev <= max_distance_px) & (ratio_rev >= r_lo) & (ratio_rev <= r_hi)
            if ok_rev.sum() != 1:
                n_ambiguous += 1
                continue
            new_uid, ended_uid = uid_arr[i], uid_arr[j]
            parent[find(new_uid)] = find(ended_uid)
            n_linked += 1
            link_rows.append({
                "exp_id": exp_id, "new_cell_uid": new_uid, "ended_cell_uid": ended_uid,
                "frame": t, "gap_frames": t - t_end, "distance_px": float(d[np.flatnonzero(ok)[0]]),
                "area_ratio": float(ratio[np.flatnonzero(ok)[0]]),
                "objects_in_frame": int(objects_in_frame.get(t, 0)),
            })
        stat_rows.append({"exp_id": exp_id, "n_new_tracks": n_new, "n_linked": n_linked,
                          "n_ambiguous": n_ambiguous, "n_no_candidate": n_none})

    links = pd.DataFrame(link_rows)
    stats = pd.DataFrame(stat_rows)
    if links.empty:
        logger.info("Gap Closing: keine eindeutigen Verknuepfungen gefunden.")
        return cells.sort_index(), links, stats

    # Ketten aufloesen: der frueheste Track der Kette gibt den Namen.
    links["root_cell_uid"] = links["new_cell_uid"].map(find)
    mapping = {u: find(u) for u in set(links["new_cell_uid"]) | set(links["ended_cell_uid"])}
    mapping = {u: root for u, root in mapping.items() if u != root}
    root_track = first["track_id"]
    out = cells.copy()
    is_mapped = out["cell_uid"].isin(mapping)
    out.loc[is_mapped, "cell_uid"] = out.loc[is_mapped, "cell_uid"].map(mapping)
    out["track_id"] = out["cell_uid"].map(root_track).astype(out["track_id"].dtype)
    logger.info(
        "Gap Closing: %d Verknuepfungen in %d Kammern (%d neue Tracks: %d verknuepft, %d mehrdeutig "
        "abgelehnt, %d ohne Kandidat); Tracks %d -> %d.",
        len(links), stats["exp_id"].nunique(), int(stats["n_new_tracks"].sum()), int(stats["n_linked"].sum()),
        int(stats["n_ambiguous"].sum()), int(stats["n_no_candidate"].sum()),
        cells["cell_uid"].nunique(), out["cell_uid"].nunique(),
    )
    return out.sort_index(), links, stats


# ==============================================================================
# 3. Sparse-Phase-Fenster
# ==============================================================================

def detect_sparse_window(
    cells: pd.DataFrame,
    max_objects: int = 20,
    smooth_frames: int = 5,
    min_frames: int = 20,
) -> pd.DataFrame:
    """Je Kammer: Frames vom Start bis zum letzten Frame, bevor die geglaettete
    Objektzahl (rollender Median ueber smooth_frames) max_objects uebersteigt.
    window_ok = False, wenn das Fenster kuerzer als min_frames ist - die
    Kammer faellt dann aus der Lineage-Auswertung."""
    rows = []
    for exp_id, g in cells.groupby("exp_id"):
        counts = g.groupby("frame").size().sort_index()
        smooth = counts.rolling(smooth_frames, center=True, min_periods=1).median()
        over = smooth.index[smooth > max_objects]
        if len(over) == 0:
            end = int(counts.index.max()); reason = "never crowded"
        else:
            end = int(over.min()) - 1; reason = f"> {max_objects} objects from frame {int(over.min())}"
        in_window = counts.index <= end
        n_frames = int(in_window.sum())
        rows.append({
            "exp_id": exp_id, "start_frame": int(counts.index.min()), "end_frame": end, "n_frames": n_frames,
            "n_frames_total": int(len(counts)),
            "objects_at_end": int(counts.loc[end]) if end in counts.index else 0,
            "objects_max": int(counts.max()),
            "window_ok": n_frames >= min_frames, "reason": reason,
        })
    window = pd.DataFrame(rows)
    if not window.empty:
        logger.info(
            "Sparse-Phase-Fenster: %d von %d Kammern auswertbar (>= %d Frames mit <= %d Objekten); "
            "Fensterlaenge Median %.0f Frames, komplett duenn: %d Kammern.",
            int(window["window_ok"].sum()), len(window), min_frames, max_objects,
            window["n_frames"].median(), int((window["reason"] == "never crowded").sum()),
        )
    return window


def flag_lineage_window(cells: pd.DataFrame, window: pd.DataFrame) -> pd.DataFrame:
    """Spalte in_lineage_window: liegt die Zeile im auswertbaren Fenster ihrer Kammer?"""
    out = cells.copy()
    if window.empty:
        out["in_lineage_window"] = True
        return out
    w = window.set_index("exp_id")
    end = out["exp_id"].map(w["end_frame"])
    ok = out["exp_id"].map(w["window_ok"]).fillna(False).astype(bool)
    out["in_lineage_window"] = ok & (out["frame"] <= end.fillna(-1))
    return out


def plot_lineage_window(cells: pd.DataFrame, window: pd.DataFrame, out_path, max_objects: int,
                        min_frames: int, max_curves: int = 80, seed: int = 0) -> None:
    """Links: Objekte pro Frame ueber die Zeit (Stichprobe von Kammern), Grenze
    max_objects. Rechts: Verteilung der Fensterlaengen ueber alle Kammern."""
    from plot_style import INK, INK_MUTED, finish, panel_title, strain_color
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.9))
    counts = cells.groupby(["exp_id", "frame"]).size().rename("n").reset_index()
    exp_ids = counts["exp_id"].unique()
    strain_of = (cells.drop_duplicates("exp_id").set_index("exp_id")["biosensor"].astype(str).to_dict()
                 if "biosensor" in cells.columns else {})
    rng = np.random.default_rng(seed)
    shown = rng.choice(exp_ids, size=min(max_curves, len(exp_ids)), replace=False) if len(exp_ids) else []
    for e in shown:
        sub = counts[counts["exp_id"] == e]
        ax.plot(sub["frame"], sub["n"], lw=0.6, alpha=0.45, color=strain_color(strain_of.get(e, "_")))
    ax.axhline(max_objects, color=INK, ls="--", lw=1.0, label=f"sparse limit ({max_objects} objects)")
    ax.set_xlabel("frame"); ax.set_ylabel("objects per frame")
    panel_title(ax, f"objects per frame, {len(shown)} of {len(exp_ids)} chambers")
    ax.set_yscale("symlog", linthresh=10); ax.legend(loc="upper left")
    if not window.empty:
        ax2.hist(window["n_frames"], bins=30, color="#d9d9d9", edgecolor="white")
        ax2.axvline(min_frames, color=INK, ls="--", lw=1.0, label=f"minimum ({min_frames} frames)")
        n_ok = int(window["window_ok"].sum())
        panel_title(ax2, f"window length, {n_ok} of {len(window)} chambers usable")
        ax2.legend()
    ax2.set_xlabel("frames in the sparse window"); ax2.set_ylabel("chambers")
    if strain_of:
        from plot_style import legend_below, strain_handles
        legend_below(fig, strain_handles(set(strain_of.values())), ncol=6, y=0.0)
    finish(fig, out_path)
