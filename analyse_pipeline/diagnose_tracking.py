"""Diagnose, warum Tracks der Cellpose-Pipeline abreissen - ohne Bilder, nur aus Combined_Results.

    python diagnose_tracking.py <Combined_Results.csv> [--qc qc_exclusions.csv] [--out DIR]

Fuer jede neue Track-ID (Frame >= 1) wird geprueft, was einen Frame vorher an derselben Stelle war:

    same object: failed IoU gate        - dasselbe Objekt (Abstand <= 1.5 r + 10 px, Flaeche <= 2.5x) stand einen
                                          Frame vorher da; die geschaetzte IoU zweier Kreise liegt < 0.3, also hat
                                          der IoU-Schwellwert des Trackers die Verbindung verworfen
    same object: failed distance gate   - dito, aber Abstand > 50 px (max_dist des Trackers)
    same object: failed area gate       - dito, Flaechenverhaeltnis > 4 (max_area_growth)
    same object: lost anyway            - dito, IoU-Schaetzung >= 0.3: Zuordnungskonkurrenz oder Maskenform
    same object: dropout                - Partner endete zwei Frames vorher (ein Frame lang nicht segmentiert/gefiltert)
    new: small object touching larger   - kein Partner; beruehrt ein groesseres Objekt (Knospe oder abgespaltenes Stueck)
    new: touching similar-sized         - kein Partner; beruehrt ein aehnlich grosses Objekt (Aufspaltung, Flackern)
    new: isolated                       - kein Partner, kein Nachbar (angespuelt, oder > 1 Frame nicht segmentiert)
    after empty frame                   - der Tracker setzt bei einem leeren Frame alle Tracks zurueck

Dazu ein Prototyp-Re-Linker (Centroid + Flaeche, Gedaechtnis ueber max_gap Frames, ungarische Zuordnung, KEINE
IoU-Schwelle; in duennen Frames (<= sparse_n Objekte) ein zweiter Durchgang mit grossem Suchradius, der nur
eindeutige Paare verbindet). Ausgabe: Fragmentierung vorher/nachher je Kammer und, wenn --qc gegeben, welcher
Anteil der manuellen merge_into_track_id-Verbindungen reproduziert wird.

Gedacht als Werkzeug fuer die Entscheidung Re-Tracking vs. Re-Segmentierung; aendert keine Pipeline-Ausgabe.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

IOU_GATE, MAX_DIST_ORIG, MAX_AREA_GROWTH = 0.3, 50.0, 4.0


def circle_iou(a1: float, a2: float, d: float) -> float:
    """IoU zweier Kreise mit Flaechen a1, a2 und Mittelpunktsabstand d (Schaetzung ohne Masken)."""
    r1, r2 = np.sqrt(a1 / np.pi), np.sqrt(a2 / np.pi)
    if d >= r1 + r2:
        inter = 0.0
    elif d <= abs(r1 - r2):
        inter = np.pi * min(r1, r2) ** 2
    else:
        a = np.arccos(np.clip((d**2 + r1**2 - r2**2) / (2 * d * r1), -1, 1))
        b = np.arccos(np.clip((d**2 + r2**2 - r1**2) / (2 * d * r2), -1, 1))
        k = 0.5 * np.sqrt(max((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2), 0.0))
        inter = r1**2 * a + r2**2 * b - k
    return float(inter / (a1 + a2 - inter))


def _pair_gate(new: pd.DataFrame, old: pd.DataFrame, dist_frac: float, dist_add: float, max_ratio: float):
    D = np.hypot(new.centroid_y.values[:, None] - old.centroid_y.values[None, :],
                 new.centroid_x.values[:, None] - old.centroid_x.values[None, :])
    rmax = np.maximum(new.r.values[:, None], old.r.values[None, :])
    ratio = (np.maximum(new.area.values[:, None], old.area.values[None, :])
             / np.minimum(new.area.values[:, None], old.area.values[None, :]))
    ok = (D <= dist_frac * rmax + dist_add) & (ratio <= max_ratio)
    return D, rmax, ratio, np.where(ok, D / rmax, 1e6)


def categorise_new_ids(df: pd.DataFrame, key: str) -> pd.DataFrame:
    rows = []
    for fn, g in df.groupby(key):
        frames = {t: sub.reset_index(drop=True) for t, sub in g.groupby("frame")}
        T = int(g.frame.max()) + 1
        last = g.groupby("track_id")["frame"].max()
        for t in range(1, T):
            C, P = frames.get(t), frames.get(t - 1)
            if C is None:
                continue
            if P is None:
                rows += [dict(unit=fn, frame=t, n_obj=len(C), category="after empty frame", area=a) for a in C.area]
                continue
            prev_ids, cur_ids = set(P.track_id), set(C.track_id)
            new = C[~C.track_id.isin(prev_ids)].reset_index(drop=True)
            if new.empty:
                continue
            cat = np.array(["?"] * len(new), dtype=object)
            info = [dict() for _ in range(len(new))]
            ended1 = P[~P.track_id.isin(cur_ids)].reset_index(drop=True)
            if len(ended1):
                D, rmax, ratio, cost = _pair_gate(new, ended1, 1.5, 10.0, 2.5)
                for i, j in zip(*linear_sum_assignment(cost)):
                    if cost[i, j] >= 1e5:
                        continue
                    iou = circle_iou(ended1.area[j], new.area[i], D[i, j])
                    if D[i, j] > max(1.5 * min(new.r[i], ended1.r[j]), MAX_DIST_ORIG):
                        cat[i] = "same object: failed distance gate (>50 px)"
                    elif ratio[i, j] > MAX_AREA_GROWTH:
                        cat[i] = "same object: failed area gate (>4x)"
                    elif iou < IOU_GATE:
                        cat[i] = "same object: failed IoU gate (est. IoU < 0.3)"
                    else:
                        cat[i] = "same object: est. IoU >= 0.3, lost anyway"
                    info[i] = dict(d_px=D[i, j], d_over_r=D[i, j] / rmax[i, j], area_ratio=ratio[i, j], est_iou=iou)
            rest = np.where(cat == "?")[0]
            P2 = frames.get(t - 2)
            ended2 = P2[P2.track_id.map(last) == t - 2].reset_index(drop=True) if P2 is not None else None
            if len(rest) and ended2 is not None and len(ended2):
                D, rmax, ratio, cost = _pair_gate(new.iloc[rest].reset_index(drop=True), ended2, 2.0, 10.0, 2.5)
                for i, j in zip(*linear_sum_assignment(cost)):
                    if cost[i, j] < 1e5:
                        cat[rest[i]] = "same object: missing for one frame (dropout)"
                        info[rest[i]] = dict(d_px=D[i, j], d_over_r=D[i, j] / rmax[i, j], area_ratio=ratio[i, j])
            for i in np.where(cat == "?")[0]:
                others = C[C.track_id != new.track_id[i]]
                if len(others):
                    D = np.hypot(others.centroid_y.values - new.centroid_y[i], others.centroid_x.values - new.centroid_x[i])
                    j = int(np.argmin(D))
                    touching = D[j] <= new.r[i] + others.r.values[j] + 8
                    if touching and new.area[i] / others.area.values[j] < 0.5:
                        cat[i] = "new: small object touching a larger one (bud or split-off piece)"
                    elif touching:
                        cat[i] = "new: touching a similar-sized object (split, flicker)"
                    else:
                        cat[i] = "new: isolated (flushed in, or absent > 1 frame)"
                else:
                    cat[i] = "new: isolated (flushed in, or absent > 1 frame)"
            rows += [dict(unit=fn, frame=t, n_obj=len(C), category=cat[i], area=new.area[i], track_id=new.track_id[i], **info[i])
                     for i in range(len(new))]
    return pd.DataFrame(rows)


def relink(g: pd.DataFrame, max_gap: int = 2, dist_frac: float = 1.5, dist_add: float = 10.0, max_ratio: float = 2.5,
           sparse_n: int = 20, sparse_frac: float = 6.0, sparse_add: float = 20.0) -> pd.Series:
    """Prototyp: neue Track-IDs je Zeile (Index von g). Siehe Modul-Docstring."""
    frames = {t: sub for t, sub in g.sort_values("frame").groupby("frame")}
    T = int(g.frame.max()) + 1
    tracks: dict[int, tuple] = {}
    nid = 0
    assign: dict = {}
    for t in range(T):
        C = frames.get(t)
        if C is None:
            continue
        cand = [k for k, v in tracks.items() if t - v[0] <= max_gap + 1]
        taken: set[int] = set()
        if cand:
            lastf = np.array([tracks[k][0] for k in cand]); cy = np.array([tracks[k][1] for k in cand])
            cx = np.array([tracks[k][2] for k in cand]); ar = np.array([tracks[k][3] for k in cand])
            gap = (t - lastf)[:, None].astype(float)
            D = np.hypot(cy[:, None] - C.centroid_y.values[None, :], cx[:, None] - C.centroid_x.values[None, :])
            rmax = np.maximum(np.sqrt(ar / np.pi)[:, None], C.r.values[None, :])
            ratio = np.maximum(ar[:, None], C.area.values[None, :]) / np.minimum(ar[:, None], C.area.values[None, :])
            ok = (D <= dist_frac * rmax * np.sqrt(gap) + dist_add) & (ratio <= max_ratio)
            cost = np.where(ok, D / rmax + 0.5 * np.log(ratio) + 0.3 * (gap - 1), 1e6)
            used_i: set[int] = set()
            for i, j in zip(*linear_sum_assignment(cost)):
                if cost[i, j] < 1e5:
                    k = cand[i]; assign[C.index[j]] = k; taken.add(j); used_i.add(i)
                    tracks[k] = (t, C.centroid_y.values[j], C.centroid_x.values[j], C.area.values[j])
            if sparse_n and len(C) <= sparse_n:
                free_i = [i for i in range(len(cand)) if i not in used_i]
                free_j = [j for j in range(len(C)) if j not in taken]
                if free_i and free_j:
                    big = (D <= sparse_frac * rmax * np.sqrt(gap) + sparse_add) & (ratio <= max_ratio)
                    sub = big[np.ix_(free_i, free_j)]
                    for a, i in enumerate(free_i):
                        js = np.where(sub[a])[0]
                        if len(js) == 1 and sub[:, js[0]].sum() == 1:
                            j = free_j[js[0]]; k = cand[i]; assign[C.index[j]] = k; taken.add(j)
                            tracks[k] = (t, C.centroid_y.values[j], C.centroid_x.values[j], C.area.values[j])
        for j, idx in enumerate(C.index):
            if j not in taken:
                nid += 1; assign[idx] = nid
                tracks[nid] = (t, C.centroid_y.values[j], C.centroid_x.values[j], C.area.values[j])
    return pd.Series(assign)


def frag_stats(g: pd.DataFrame, col: str) -> dict:
    g = g.sort_values("frame")
    first = g.groupby(col)["frame"].min()
    n_new = int((g[col].map(first) == g.frame).sum() - (g.frame == 0).sum())
    n_obj = int((g.frame > 0).sum())
    length = g.groupby(col)["frame"].agg(lambda s: s.max() - s.min() + 1)
    return dict(new_ids_per_object_frame=n_new / max(n_obj, 1), median_track_frames=float(length.median()),
                share_tracks_ge_10_frames=float((length >= 10).mean()), n_tracks=int(g[col].nunique()))


def manual_link_recall(df: pd.DataFrame, key: str, qc_path: Path, col: str) -> pd.DataFrame:
    """Welche manuellen Verbindungen (merge_into_track_id) reproduziert der Prototyp? Zuordnung ueber
    (condition, replicate, chamber, track_id) - die letzten Teile der cell_uid."""
    q = pd.read_csv(qc_path)
    m = q[q.merge_into_track_id.notna()].copy()
    parts = m.cell_uid.str.split("__")
    m["condition"], m["replicate"], m["chamber"] = parts.str[-4], parts.str[-3], parts.str[-2]
    m["track_id"] = parts.str[-1].str.replace("track", "").astype(int)
    m["target_id"] = m.merge_into_track_id.astype(int)
    idx = df.drop_duplicates([key, "track_id"]).set_index(["condition", "replicate", "chamber", "track_id"])[key]
    newid = df.groupby([key, "track_id"])[col].agg(lambda s: s.mode().iloc[0])
    lastf = df.groupby([key, "track_id"])["frame"].max(); firstf = df.groupby([key, "track_id"])["frame"].min()
    rows = []
    for r in m.itertuples():
        ka, kb = (r.condition, r.replicate, r.chamber, r.track_id), (r.condition, r.replicate, r.chamber, r.target_id)
        if ka not in idx.index or kb not in idx.index:
            continue
        a, b = (idx.loc[ka], r.track_id), (idx.loc[kb], r.target_id)
        lo, hi = (a, b) if lastf.loc[a] <= lastf.loc[b] else (b, a)
        rows.append(dict(cell_uid=r.cell_uid, gap_frames=int(firstf.loc[hi] - lastf.loc[lo]), reproduced=bool(newid.loc[a] == newid.loc[b])))
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results", type=Path, help="Combined_Results.csv/.parquet der Cellpose-Pipeline")
    ap.add_argument("--qc", type=Path, default=None, help="qc_exclusions.csv mit merge_into_track_id")
    ap.add_argument("--out", type=Path, default=None, help="Ausgabeordner (Default: neben der Eingabe)")
    ap.add_argument("--max-gap", type=int, default=2)
    args = ap.parse_args()
    df = pd.read_parquet(args.results) if args.results.suffix == ".parquet" else pd.read_csv(args.results)
    key = "filename" if "filename" in df.columns else "exp_id"
    n0 = len(df)
    df = df.dropna(subset=["area", "centroid_x", "centroid_y", "track_id"]).copy()
    if len(df) < n0:
        print(f"{n0 - len(df)} rows without track_id/geometry dropped")
    df["track_id"] = df["track_id"].astype(int)
    df["r"] = np.sqrt(df.area / np.pi)
    out = args.out or args.results.parent
    out.mkdir(parents=True, exist_ok=True)

    cats = categorise_new_ids(df, key)
    cats.to_csv(out / "tracking_new_id_categories.csv", index=False)
    n_obj = int((df.frame > 0).shape[0])
    print(f"{df[key].nunique()} units, {n_obj} object-frames (frame >= 1), {len(cats)} new track IDs "
          f"({len(cats) / n_obj:.3f} per object-frame)")
    s = cats.category.value_counts()
    print(pd.DataFrame({"n": s, "share": (s / s.sum()).round(3)}).to_string())
    length = df.groupby([key, "track_id"])["frame"].agg(lambda s: s.max() - s.min() + 1)
    print(f"tracks lasting 1 frame: {(length == 1).mean():.2f}, <= 2 frames: {(length <= 2).mean():.2f}")

    rows = []
    df["relinked_id"] = np.nan
    for fn, g in df.groupby(key):
        s = relink(g, max_gap=args.max_gap)
        df.loc[s.index, "relinked_id"] = s.values
        b, a = frag_stats(g, "track_id"), frag_stats(df.loc[g.index], "relinked_id")
        rows.append(dict(unit=fn, **{f"{k}_before": v for k, v in b.items()}, **{f"{k}_after": v for k, v in a.items()}))
    res = pd.DataFrame(rows)
    res.to_csv(out / "tracking_relink_prototype.csv", index=False)
    print("\nprototype re-linker (centroid + area, memory, no IoU gate), mean over units:")
    print(res.drop(columns="unit").mean().round(3).to_string())
    if args.qc is not None and args.qc.exists():
        rec = manual_link_recall(df, key, args.qc, "relinked_id")
        rec.to_csv(out / "tracking_manual_link_recall.csv", index=False)
        if len(rec):
            print(f"\nmanual merge links usable: {len(rec)}, reproduced: {rec.reproduced.mean():.2f}")
            print(rec.groupby(pd.cut(rec.gap_frames, [-1000, 0, 1, 2, 5, 1000], labels=["overlap", "1", "2", "3-5", ">5"]),
                              observed=True)["reproduced"].agg(["size", "mean"]).round(2).to_string())


if __name__ == "__main__":
    main()
