"""Synthetischer Ground-Truth-Test fuer track_labels.py gegen einen v11-artigen Tracker.

Erzeugt einen Label-Stack mit bekannten Identitaeten: wandernde, wachsende Ellipsen mit gelegentlichen
Spruengen, Knospen (beruehren die Mutter, wachsen, loesen sich), ausgelassene Frames (nicht segmentiert),
falsche Aufspaltungen, verschmolzene Masken, kurz auftauchende Zellen und ein leerer Frame. Beide Tracker
bekommen dieselben Label-Bilder.

    python synthetic_tracking_test.py [--seed 0] [--frames 80] [--out DIR]
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import linear_sum_assignment
sys.path.insert(0, str(Path(__file__).resolve().parent))
from track_labels import LabelTracker, TrackParams, frame_features, overlap_counts, tracking_summary  # noqa: E402


def ellipse_radius(c, ang):
    """Abstand Mittelpunkt -> Rand der Ellipse c in Bildrichtung ang (sin -> y, cos -> x)."""
    th = ang - c["ang"]
    return c["a"] * c["b"] / np.sqrt((c["b"] * np.cos(th)) ** 2 + (c["a"] * np.sin(th)) ** 2)


def make_movie(T=80, H=700, W=700, n0=30, seed=0):
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:H, 0:W]
    cells = {}   # gt_id -> dict(cy, cx, a, b, ang, alive, bud_of, is_bud, detach_r, transient_until)
    nid = [0]
    def new_cell(cy, cx, a, b, is_bud=False, bud_of=None, transient_until=None):
        nid[0] += 1
        cells[nid[0]] = dict(cy=cy, cx=cx, a=a, b=b, ang=rng.uniform(0, np.pi), alive=True, is_bud=is_bud,
                             bud_of=bud_of, transient_until=transient_until, born=None)
        return nid[0]
    def overlaps(cy, cx, r, ignore=()):
        for k, c in cells.items():
            if k in ignore or not c["alive"]: continue
            if np.hypot(c["cy"] - cy, c["cx"] - cx) < r + max(c["a"], c["b"]) + 3: return True
        return False
    for _ in range(n0):
        for _try in range(200):
            a = rng.uniform(18, 38); b = a * rng.uniform(0.6, 1.0)
            cy, cx = rng.uniform(60, H - 60), rng.uniform(60, W - 60)
            if not overlaps(cy, cx, max(a, b)): new_cell(cy, cx, a, b); break
    stack = np.zeros((T, H, W), dtype=np.int32); gt_rows = []; split_counter = [200000]
    for t in range(T):
        # bewegen, wachsen, Knospen
        for k in list(cells):
            c = cells[k]
            if not c["alive"]: continue
            if c["transient_until"] is not None and t > c["transient_until"]: c["alive"] = False; continue
            if c["is_bud"] and c["bud_of"] is not None:
                m = cells[c["bud_of"]]; c["a"] *= 1.06; c["b"] *= 1.06
                ang = c["ang_to_mother"]; rr = ellipse_radius(m, ang) + max(c["a"], c["b"]) - 2
                c["cy"], c["cx"] = m["cy"] + rr * np.sin(ang), m["cx"] + rr * np.cos(ang)
                if max(c["a"], c["b"]) >= 16: c["bud_of"] = None   # abgeloest
            else:
                s = 1.002; c["a"] *= s; c["b"] *= s
                dy, dx = rng.normal(0, 3, 2)
                if rng.random() < 0.03:  # Sprung
                    ang = rng.uniform(0, 2 * np.pi); d = rng.uniform(40, 90); dy, dx = d * np.sin(ang), d * np.cos(ang)
                ny, nx = np.clip(c["cy"] + dy, 40, H - 40), np.clip(c["cx"] + dx, 40, W - 40)
                if not overlaps(ny, nx, max(c["a"], c["b"]), ignore=(k,)): c["cy"], c["cx"] = ny, nx
            if not c["is_bud"] and max(c["a"], c["b"]) >= 22 and rng.random() < 0.015 and not any(
                    v["bud_of"] == k for v in cells.values() if v["alive"]):
                ang = rng.uniform(0, 2 * np.pi); rr = ellipse_radius(c, ang) + 7 - 2
                bid = new_cell(c["cy"] + rr * np.sin(ang), c["cx"] + rr * np.cos(ang), 7, 6, is_bud=True, bud_of=k)
                cells[bid]["ang_to_mother"] = ang; cells[bid]["born"] = t
        if rng.random() < 0.15:  # kurz auftauchende Zelle
            for _try in range(50):
                a = rng.uniform(15, 30); cy, cx = rng.uniform(40, H - 40), rng.uniform(40, W - 40)
                if not overlaps(cy, cx, a): new_cell(cy, cx, a, a * 0.8, transient_until=t + rng.integers(0, 3)); break
        # zeichnen
        img = np.zeros((H, W), dtype=np.int32)
        alive = [k for k, c in cells.items() if c["alive"]]
        empty_frame = (t == 40)
        for k in alive:
            c = cells[k]; drawn = True; label = k; kind = "normal"
            if empty_frame or (rng.random() < 0.05 and cells[k]["born"] != t): drawn = False; kind = "dropout"
            if drawn:
                ca, sa = np.cos(c["ang"]), np.sin(c["ang"])
                X = (xx - c["cx"]) * ca + (yy - c["cy"]) * sa; Y = -(xx - c["cx"]) * sa + (yy - c["cy"]) * ca
                mask = (X / c["a"]) ** 2 + (Y / c["b"]) ** 2 <= 1
                if rng.random() < 0.02 and not c["is_bud"]:   # falsche Aufspaltung
                    split_counter[0] += 1; kind = "split"
                    img[mask & (X < 0)] = k; img[mask & (X >= 0)] = split_counter[0]
                else:
                    if rng.random() < 0.02:  # mit einem beruehrenden Nachbarn verschmolzen
                        for k2 in alive:
                            c2 = cells[k2]
                            if k2 != k and np.hypot(c2["cy"] - c["cy"], c2["cx"] - c["cx"]) < max(c["a"], c["b"]) + max(c2["a"], c2["b"]) + 4:
                                label = k2; kind = "merged_into_neighbour"; break
                    img[mask] = label
            gt_rows.append(dict(frame=t, gt_id=k, label=label if drawn else -1, drawn=drawn and kind != "merged_into_neighbour",
                                kind=kind, is_bud=c["is_bud"], bud_of=c["bud_of"] if c["is_bud"] else -1,
                                transient=c["transient_until"] is not None, born=c["born"]))
        stack[t] = img
    return stack, pd.DataFrame(gt_rows)


def v11_like_track(stack, iou_thr=0.3, floor=0.05, small_r=35.0, max_dist=50.0, dist_frac=1.5, max_ratio=7.0):
    """Frame-zu-Frame, IoU-Tor, kein Gedaechtnis, leerer Frame setzt zurueck (Pipeline v11 mit der Config)."""
    rows = []; prev = None; prev_objs = None; next_id = 1; ids_prev = {}
    for t in range(stack.shape[0]):
        lab = stack[t].astype(np.int64); objs = frame_features(lab)
        ids = {}
        if len(objs) == 0:
            prev = None; prev_objs = None; ids_prev = {}
            continue
        if prev is not None and len(prev_objs):
            oc = overlap_counts(prev, lab)
            cost = np.ones((len(prev_objs), len(objs)))
            for i, p in enumerate(prev_objs.itertuples()):
                for j, c in enumerate(objs.itertuples()):
                    d = np.hypot(p.centroid_y - c.centroid_y, p.centroid_x - c.centroid_x)
                    rmin = min(p.r, c.r)
                    if d > max(dist_frac * rmin, max_dist): continue
                    ratio = max(p.area, c.area) / min(p.area, c.area)
                    if ratio > max_ratio: continue
                    inter = oc.get((int(p.label), int(c.label)), 0); iou = inter / (p.area + c.area - inter)
                    if iou < (floor if rmin <= small_r else iou_thr): continue
                    cost[i, j] = 1 - iou
            for i, j in zip(*linear_sum_assignment(cost)):
                if cost[i, j] < 1.0: ids[int(objs.label.values[j])] = ids_prev[int(prev_objs.label.values[i])]
        for l in objs.label.values:
            if int(l) not in ids: ids[int(l)] = next_id; next_id += 1
        for l in objs.label.values:
            rows.append(dict(frame=t, label=int(l), track_id=ids[int(l)], parent_track_id=-1, link_type=""))
        prev, prev_objs, ids_prev = lab, objs, ids
    return pd.DataFrame(rows)


def score(table: pd.DataFrame, gt: pd.DataFrame, memory=3) -> dict:
    g = gt[gt.drawn].merge(table[["frame", "label", "track_id", "parent_track_id"]], on=["frame", "label"], how="left")
    g = g.dropna(subset=["track_id"]).sort_values(["gt_id", "frame"])
    # Verbindungen: aufeinanderfolgende Auftritte desselben gt-Objekts mit Luecke <= memory
    links = []
    for k, sub in g.groupby("gt_id"):
        f = sub.frame.values; tid = sub.track_id.values
        for i in range(1, len(f)):
            gap = f[i] - f[i - 1]
            if gap <= memory:
                links.append(dict(gt_id=k, gap=gap, ok=tid[i] == tid[i - 1], transient=bool(sub.transient.iloc[0])))
    L = pd.DataFrame(links)
    frag = g.groupby("gt_id")["track_id"].nunique(); n_app = g.groupby("gt_id").size()
    long_objs = frag[n_app >= 5]
    # Identitaetswechsel: eine track_id auf zwei gt-Objekten (Aufspaltungshaelften ausgenommen)
    merged_targets = set(zip(gt.loc[gt.kind == "merged_into_neighbour", "frame"], gt.loc[gt.kind == "merged_into_neighbour", "label"]))
    clean = g[(g.kind != "split") & ~pd.Series(list(zip(g.frame, g.label)), index=g.index).isin(merged_targets)]
    per_track = clean.groupby("track_id")["gt_id"].nunique()
    # Knospen: Elternteil im ersten gezeichneten Frame
    buds = g[g.is_bud & (g.bud_of > 0)]
    first_bud = buds.groupby("gt_id").head(1)
    mother_tid = g.set_index(["gt_id", "frame"])["track_id"]
    ok_parent = []
    for r in first_bud.itertuples():
        mt = mother_tid.get((r.bud_of, r.frame), np.nan)
        ok_parent.append(r.parent_track_id == mt)
    return dict(
        link_recall_gap1=float(L[L.gap == 1].ok.mean()), link_recall_gap2_3=float(L[L.gap > 1].ok.mean()) if (L.gap > 1).any() else np.nan,
        n_links=int(len(L)), tracks_per_object_ge5=float(long_objs.mean()), share_objects_single_track=float((long_objs == 1).mean()),
        identity_switch_tracks=int((per_track > 1).sum()), n_tracks=int(table.track_id.nunique()),
        buds=int(len(first_bud)), bud_parent_assigned=float((first_bud.parent_track_id > 0).mean()) if len(first_bud) else np.nan,
        bud_parent_correct=float(np.mean(ok_parent)) if ok_parent else np.nan,
    )


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--seed", type=int, default=0); ap.add_argument("--frames", type=int, default=80)
    ap.add_argument("--out", type=Path, default=None); a = ap.parse_args()
    stack, gt = make_movie(T=a.frames, seed=a.seed)
    n_obj = int((stack > 0).any(axis=(1, 2)).sum())
    print(f"movie: {a.frames} frames, {gt.gt_id.nunique()} objects incl. {gt[gt.is_bud].gt_id.nunique()} buds and "
          f"{gt[gt.transient].gt_id.nunique()} transient cells; kinds: {gt.kind.value_counts().to_dict()}")
    old = v11_like_track(stack)
    tr = LabelTracker(TrackParams())
    for t in range(stack.shape[0]):
        tr.update(t, stack[t].astype(np.int64))
    new = tr.table()
    res = pd.DataFrame({"v11-like (IoU gate, no memory)": score(old, gt), "track_labels (memory, cost, overlap)": score(new, gt)})
    print(res.round(3).to_string())
    print("\nnew tracker summary:", json.dumps(tracking_summary(new, pd.DataFrame(tr.events)), indent=None))
    if a.out:
        a.out.mkdir(parents=True, exist_ok=True); res.to_csv(a.out / "synthetic_tracking_scores.csv")
        gt.to_csv(a.out / "synthetic_gt.csv", index=False); new.to_csv(a.out / "synthetic_tracks.csv", index=False)


if __name__ == "__main__":
    main()
