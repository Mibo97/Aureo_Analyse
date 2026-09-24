"""track_labels.py - Tracking auf gespeicherten Label-Stacks (zarr, T x H x W), ohne Cellpose.

Eingabe ist ein Label-Stack pro Film: entweder die masks_<stem>.zarr der Pipeline v11 (Labels = alte
Track-IDs; je Frame ist jeder Wert genau ein Objekt) oder ein roher Cellpose-Label-Stack. Ausgabe sind
neue Track-IDs mit Elternbeziehung, eine Objekttabelle mit Formmerkmalen und ein Ereignisprotokoll.

Warum ein zweiter Durchgang statt Tracking im Bild-Loop: das Tracking ist der billige Teil (CPU, Minuten
pro Film) und kann so beliebig oft mit anderen Parametern laufen, ohne die Segmentierung (GPU, Stunden)
zu wiederholen.

Regeln (Begruendung in docs/tracking_diagnosis.md):
  * Gedaechtnis: eine Spur bleibt `memory` Frames lang verknuepfbar (ein nicht segmentierter Frame
    beendet sie nicht; ein leerer Frame setzt nichts zurueck).
  * Kosten statt Schwelle: Abstand in Radien, Flaechenverhaeltnis und Maskenueberlappung gehen als
    Kosten in eine ungarische Zuordnung ein; Tore gibt es nur fuer Abstand (dist_frac * r * sqrt(gap)
    + dist_add) und Flaechenverhaeltnis (max_ratio). KEIN IoU-Tor: eine Zelle, die eine Radiuslaenge
    weiterrutscht, bleibt dieselbe Zelle.
  * Duenne Frames (<= sparse_n Objekte): ein zweiter Durchgang verknuepft uebrig gebliebene Paare ueber
    einen grossen Radius, aber nur, wenn das Paar innerhalb dieses Radius eindeutig ist.
  * Merge/Split aus der Ueberlappung: deckt ein Objekt zwei Vorgaenger ab, laeuft die groessere Spur
    weiter und die kleinere endet mit Vermerk; zerfaellt ein Vorgaenger in zwei Objekte, bekommt das
    kleinere eine neue ID mit parent_track_id = Vorgaenger (link_type 'split').
  * Elternregel: ein neues Objekt, dessen Maske (um touch_px erweitert) eine getrackte Maske beruehrt,
    bekommt parent_track_id = beruehrte Spur (link_type 'new_touching'), zusammen mit parent_area_ratio.
    Ob das eine Knospe ist, entscheidet die Analyse (Groessenverhaeltnis, Persistenz), nicht der Tracker.

Aufrufe:
  python track_labels.py masks_<stem>.zarr --out DIR
  python track_labels.py --batch <02_processed> --combined <03_results>/Combined_Results.csv --out DIR
      verarbeitet alle masks_*.zarr eines Experiments und schreibt neben die alte Tabelle eine
      Combined_Results_retracked.csv: alle alten Spalten (inkl. Fluoreszenz), track_id = NEUE ID,
      track_id_v11 = alte ID, dazu parent_track_id, link_type, gap_frames und die Formmerkmale.
"""
from __future__ import annotations

import argparse
import json
import logging
import re
import time
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from scipy.optimize import linear_sum_assignment
from skimage.measure import regionprops_table

logger = logging.getLogger("track_labels")
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


@dataclass
class TrackParams:
    memory: int = 3              # Frames, die eine verlorene Spur verknuepfbar bleibt (gap <= memory)
    dist_frac: float = 1.5       # Tor: d <= dist_frac * r_max * sqrt(gap) + dist_add
    dist_add: float = 10.0
    max_ratio: float = 2.5       # Tor: max(area)/min(area) <= max_ratio
    w_dist: float = 1.0          # Kosten: Abstand in Einheiten des groesseren Radius
    w_area: float = 0.5          # Kosten: |ln Flaechenverhaeltnis|
    w_gap: float = 0.3           # Kosten je uebersprungenem Frame
    w_overlap: float = 1.0       # Gutschrift: Ueberlappung / kleinere Flaeche (0..1)
    sparse_n: int = 20           # Frames mit <= sparse_n Objekten gelten als duenn
    sparse_frac: float = 6.0     # grosser Radius im duennen Durchgang (nur eindeutige Paare)
    sparse_add: float = 20.0
    merge_min_frac: float = 0.5  # Objekt deckt >= diesen Anteil ZWEIER Vorgaenger ab -> merge
    split_min_frac: float = 0.3  # Vorgaenger deckt >= diesen Anteil zweier Nachfolger ab -> split
    touch_px: int = 2            # Erweiterung der Maske fuer die Beruehrungsregel
    skeleton: bool = False       # Skelettlaenge je Objekt (langsam, fuer Zelltypen nach Rensink 2026)


FEATURE_COLS = ["area", "centroid_y", "centroid_x", "bbox_r0", "bbox_c0", "bbox_r1", "bbox_c1",
                "eccentricity", "solidity", "axis_major", "axis_minor", "perimeter", "circularity"]


# ----------------------------------------------------------------------------------------------
# Merkmale und Ueberlappungen
# ----------------------------------------------------------------------------------------------
def frame_features(labels: np.ndarray, skeleton: bool = False) -> pd.DataFrame:
    """Eine Zeile je Objekt (= je Labelwert > 0) mit Formmerkmalen."""
    if labels.max() == 0:
        return pd.DataFrame(columns=["label"] + FEATURE_COLS)
    props = regionprops_table(
        labels, properties=("label", "area", "centroid", "bbox", "eccentricity", "solidity",
                            "axis_major_length", "axis_minor_length", "perimeter"))
    df = pd.DataFrame(props).rename(columns={
        "centroid-0": "centroid_y", "centroid-1": "centroid_x", "bbox-0": "bbox_r0", "bbox-1": "bbox_c0",
        "bbox-2": "bbox_r1", "bbox-3": "bbox_c1", "axis_major_length": "axis_major",
        "axis_minor_length": "axis_minor"})
    df["area"] = df["area"].astype(float)
    df["circularity"] = np.where(df["perimeter"] > 0, 4 * np.pi * df["area"] / np.maximum(df["perimeter"], 1e-9) ** 2, np.nan)
    df["circularity"] = df["circularity"].clip(upper=1.5)
    if skeleton:
        from skimage.morphology import skeletonize
        sk = []
        for r in df.itertuples():
            crop = labels[int(r.bbox_r0):int(r.bbox_r1), int(r.bbox_c0):int(r.bbox_c1)] == r.label
            sk.append(int(skeletonize(crop).sum()))
        df["skeleton_px"] = sk
    df["r"] = np.sqrt(df["area"] / np.pi)
    return df.reset_index(drop=True)


def overlap_counts(prev: np.ndarray, curr: np.ndarray) -> dict[tuple[int, int], int]:
    """Pixelzahl je (prev_label, curr_label)-Paar, nur wo beide > 0."""
    m = (prev > 0) & (curr > 0)
    if not m.any():
        return {}
    a = prev[m].astype(np.int64); b = curr[m].astype(np.int64)
    k = int(curr.max()) + 1
    pair, cnt = np.unique(a * k + b, return_counts=True)
    return {(int(p // k), int(p % k)): int(c) for p, c in zip(pair, cnt)}


def _touching_labels(labels: np.ndarray, row, touch_px: int) -> np.ndarray:
    """Labels (> 0, != eigenes), die die um touch_px erweiterte Maske des Objekts beruehren."""
    r0, c0 = max(int(row.bbox_r0) - touch_px, 0), max(int(row.bbox_c0) - touch_px, 0)
    r1, c1 = min(int(row.bbox_r1) + touch_px, labels.shape[0]), min(int(row.bbox_c1) + touch_px, labels.shape[1])
    crop = labels[r0:r1, c0:c1]
    own = crop == row.label
    ring = ndi.binary_dilation(own, iterations=touch_px) & ~own
    vals = crop[ring]
    return np.unique(vals[vals > 0])


# ----------------------------------------------------------------------------------------------
# Tracker
# ----------------------------------------------------------------------------------------------
class LabelTracker:
    def __init__(self, params: Optional[TrackParams] = None):
        self.p = params or TrackParams()
        self.tracks: dict[int, dict] = {}        # track_id -> dict(last_frame, label, area, cy, cx, r)
        self.next_id = 1
        self.recent: dict[int, np.ndarray] = {}  # frame -> label image (nur die letzten memory Frames)
        self.rows: list[dict] = []
        self.events: list[dict] = []
        self.n_empty_frames = 0

    # -- Hilfen
    def _candidates(self, t: int) -> list[int]:
        return [k for k, v in self.tracks.items() if 0 < t - v["last_frame"] <= self.p.memory]

    def _cost_matrix(self, t: int, cand: list[int], objs: pd.DataFrame, labels: np.ndarray):
        p = self.p
        n_c, n_o = len(cand), len(objs)
        D = np.empty((n_c, n_o)); rmax = np.empty((n_c, n_o)); ratio = np.empty((n_c, n_o))
        ov = np.zeros((n_c, n_o)); gap = np.empty((n_c, 1))
        ov_cache: dict[int, dict] = {}
        oy, ox, oa, orad = objs.centroid_y.values, objs.centroid_x.values, objs.area.values, objs.r.values
        for i, k in enumerate(cand):
            v = self.tracks[k]; g = t - v["last_frame"]; gap[i, 0] = g
            D[i] = np.hypot(v["cy"] - oy, v["cx"] - ox)
            rmax[i] = np.maximum(v["r"], orad)
            ratio[i] = np.maximum(v["area"], oa) / np.minimum(v["area"], oa)
            prev_frame = v["last_frame"]
            if prev_frame in self.recent:
                if prev_frame not in ov_cache:
                    ov_cache[prev_frame] = overlap_counts(self.recent[prev_frame], labels)
                oc = ov_cache[prev_frame]
                for j, lab in enumerate(objs.label.values):
                    c = oc.get((v["label"], int(lab)))
                    if c:
                        ov[i, j] = c / min(v["area"], oa[j])
        ok = (D <= p.dist_frac * rmax * np.sqrt(gap) + p.dist_add) & (ratio <= p.max_ratio)
        cost = np.where(ok, p.w_dist * D / rmax + p.w_area * np.log(ratio) + p.w_gap * (gap - 1) - p.w_overlap * ov, np.inf)
        return cost, ok, D, rmax, ratio, ov, gap

    # -- ein Frame
    def update(self, t: int, labels: np.ndarray, objs: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        p = self.p
        if objs is None:
            objs = frame_features(labels, skeleton=p.skeleton)
        n_o = len(objs)
        assign = np.full(n_o, -1, dtype=int)       # Index -> track_id
        link_type = np.array(["new"] * n_o, dtype=object)
        gap_frames = np.zeros(n_o, dtype=int)
        ov_prev = np.zeros(n_o); n_cand = np.zeros(n_o, dtype=int)
        parent = np.full(n_o, -1, dtype=int); parent_ratio = np.full(n_o, np.nan)
        if n_o == 0:
            self.n_empty_frames += 1
        cand = self._candidates(t)
        if cand and n_o:
            cost, ok, D, rmax, ratio, ov, gap = self._cost_matrix(t, cand, objs, labels)
            n_cand[:] = ok.sum(axis=0)
            finite = np.where(np.isfinite(cost), cost, 1e6)
            used = set()
            for i, j in zip(*linear_sum_assignment(finite)):
                if ok[i, j]:
                    k = cand[i]; assign[j] = k; used.add(i)
                    g = int(gap[i, 0]); gap_frames[j] = g; ov_prev[j] = ov[i, j]
                    link_type[j] = "continued" if g == 1 else "gap"
                    if g > 1:
                        self.events.append(dict(frame=t, type="gap", track_id=k, other_track_id=-1, gap_frames=g))
            # duenner Frame: uebrig gebliebene Paare ueber grossen Radius, nur wenn eindeutig
            if p.sparse_n and n_o <= p.sparse_n:
                free_i = [i for i in range(len(cand)) if i not in used]
                free_j = [j for j in range(n_o) if assign[j] < 0]
                if free_i and free_j:
                    big = (D <= p.sparse_frac * rmax * np.sqrt(gap) + p.sparse_add) & (ratio <= p.max_ratio)
                    sub = big[np.ix_(free_i, free_j)]
                    for a, i in enumerate(free_i):
                        js = np.where(sub[a])[0]
                        if len(js) == 1 and sub[:, js[0]].sum() == 1:
                            j = free_j[js[0]]; k = cand[i]; assign[j] = k; used.add(i)
                            g = int(gap[i, 0]); gap_frames[j] = g; link_type[j] = "long_range"
                            self.events.append(dict(frame=t, type="long_range", track_id=k, other_track_id=-1, gap_frames=g))
            # merge: ein zugeordnetes Objekt deckt zusaetzlich einen NICHT zugeordneten Vorgaenger (gap 1) ab
            prev_frame = t - 1
            if prev_frame in self.recent:
                oc = overlap_counts(self.recent[prev_frame], labels)
                lab_to_j = {int(l): j for j, l in enumerate(objs.label.values)}
                for i, k in enumerate(cand):
                    v = self.tracks[k]
                    if i in used or v["last_frame"] != prev_frame:
                        continue
                    for j, lab in enumerate(objs.label.values):
                        c = oc.get((v["label"], int(lab)), 0)
                        if assign[j] >= 0 and assign[j] != k and c >= p.merge_min_frac * v["area"]:
                            self.events.append(dict(frame=t, type="merge", track_id=k, other_track_id=int(assign[j]), gap_frames=1))
                            v["merged_into"] = int(assign[j]); v["last_frame"] = -10**9  # Spur beendet
                            break
                # split: ein Vorgaenger (gap 1, zugeordnet) deckt ein weiteres, freies Objekt ab
                for i, k in enumerate(cand):
                    v = self.tracks[k]
                    if i not in used or v["last_frame"] != prev_frame:
                        continue
                    for j, lab in enumerate(objs.label.values):
                        if assign[j] >= 0:
                            continue
                        c = oc.get((v["label"], int(lab)), 0)
                        if c >= p.split_min_frac * objs.area.values[j]:
                            nid = self.next_id; self.next_id += 1
                            assign[j] = nid; link_type[j] = "split"; parent[j] = k
                            parent_ratio[j] = objs.area.values[j] / v["area"]
                            self.events.append(dict(frame=t, type="split", track_id=nid, other_track_id=k, gap_frames=1))
        # neue Objekte: ID vergeben, Beruehrungsregel fuer die Elternschaft
        lab_to_track = {int(l): int(assign[j]) for j, l in enumerate(objs.label.values) if assign[j] >= 0}
        for j in np.where(assign < 0)[0]:
            nid = self.next_id; self.next_id += 1; assign[j] = nid
            touched = _touching_labels(labels, objs.iloc[j], p.touch_px) if n_o > 1 else np.array([])
            touched_tracks = [(lab_to_track[int(l)], int(l)) for l in touched if int(l) in lab_to_track and lab_to_track[int(l)] != nid]
            if touched_tracks:
                # Elternteil = beruehrtes, bereits getracktes Objekt mit der groessten Flaeche
                areas = {int(l): a for l, a in zip(objs.label.values, objs.area.values)}
                k, lab = max(touched_tracks, key=lambda x: areas[x[1]])
                parent[j] = k; parent_ratio[j] = objs.area.values[j] / areas[lab]; link_type[j] = "new_touching"
            lab_to_track[int(objs.label.values[j])] = nid
        # Spuren aktualisieren
        for j in range(n_o):
            k = int(assign[j]); r = objs.iloc[j]
            self.tracks[k] = dict(last_frame=t, label=int(r.label), area=float(r.area), cy=float(r.centroid_y),
                                  cx=float(r.centroid_x), r=float(r.r))
        # Label-Bilder fuer die Ueberlappung der naechsten memory Frames behalten
        self.recent[t] = labels
        for f in [f for f in self.recent if f < t - p.memory]:
            del self.recent[f]
        out = objs.copy()
        out.insert(0, "frame", t)
        out["track_id"] = assign; out["parent_track_id"] = parent; out["parent_area_ratio"] = parent_ratio
        out["link_type"] = link_type; out["gap_frames"] = gap_frames; out["overlap_prev_frac"] = ov_prev
        out["n_candidates"] = n_cand
        self.rows.append(out)
        return out

    def table(self) -> pd.DataFrame:
        if not self.rows:
            return pd.DataFrame()
        return pd.concat(self.rows, ignore_index=True)


def track_stack(stack, params: Optional[TrackParams] = None, frames: Optional[range] = None,
                progress_every: int = 25) -> tuple[pd.DataFrame, pd.DataFrame, "LabelTracker"]:
    """Trackt einen (T, H, W)-Stack (zarr-Array oder numpy). Gibt Objekttabelle, Ereignisse, Tracker zurueck."""
    tr = LabelTracker(params)
    T = stack.shape[0]
    frames = frames if frames is not None else range(T)
    t0 = time.perf_counter()
    for n, t in enumerate(frames):
        labels = np.asarray(stack[t]).astype(np.int64)
        tr.update(t, labels)
        if progress_every and (n + 1) % progress_every == 0:
            logger.info("  Frame %d/%d, %d Spuren bisher, %.1f s", n + 1, len(frames), tr.next_id - 1, time.perf_counter() - t0)
    return tr.table(), pd.DataFrame(tr.events), tr


def relabel_stack(stack, table: pd.DataFrame) -> np.ndarray:
    """Neuer Stack mit track_id als Label (gleiches Format wie masks_<stem>.zarr der Pipeline v11)."""
    T, H, W = stack.shape
    out = np.zeros((T, H, W), dtype=np.uint32)
    for t, sub in table.groupby("frame"):
        lab = np.asarray(stack[t]).astype(np.int64)
        lut = np.zeros(int(lab.max()) + 1, dtype=np.uint32)
        lut[sub.label.values.astype(int)] = sub.track_id.values.astype(np.uint32)
        out[t] = lut[lab]
    return out


def tracking_summary(table: pd.DataFrame, events: pd.DataFrame) -> dict:
    if table.empty:
        return dict(n_objects=0)
    t = table.sort_values("frame")
    first = t.groupby("track_id")["frame"].min()
    is_new = (t.track_id.map(first) == t.frame) & (t.frame > t.frame.min())
    length = t.groupby("track_id")["frame"].agg(lambda s: s.max() - s.min() + 1)
    ev = events["type"].value_counts().to_dict() if not events.empty else {}
    return dict(
        n_objects=int(len(t)), n_tracks=int(t.track_id.nunique()),
        new_ids_per_object_frame=float(is_new.sum() / max((t.frame > t.frame.min()).sum(), 1)),
        median_track_frames=float(length.median()), share_tracks_ge_10_frames=float((length >= 10).mean()),
        share_tracks_1_frame=float((length == 1).mean()),
        n_new_touching=int((t.link_type == "new_touching").sum()),
        n_gap_links=int(ev.get("gap", 0)), n_long_range=int(ev.get("long_range", 0)),
        n_merge=int(ev.get("merge", 0)), n_split=int(ev.get("split", 0)),
    )


# ----------------------------------------------------------------------------------------------
# Dateien
# ----------------------------------------------------------------------------------------------
def open_stack(path: Path):
    import zarr
    return zarr.open(str(path), mode="r")


def save_stack(path: Path, arr: np.ndarray) -> None:
    import zarr
    kw = dict(mode="w", shape=arr.shape, chunks=(1,) + arr.shape[1:], dtype=arr.dtype)
    try:
        z = zarr.open(str(path), zarr_format=2, **kw)
    except TypeError:
        z = zarr.open(str(path), **kw)
    for t in range(arr.shape[0]):
        z[t] = arr[t]


def stem_of(path: Path) -> str:
    name = path.name
    name = re.sub(r"\.zarr$", "", name)
    return re.sub(r"^(masks|labels|tracks)_", "", name)


def process_stack(path: Path, out_dir: Path, params: TrackParams, write_zarr: bool = True) -> tuple[pd.DataFrame, dict]:
    stem = stem_of(path)
    stack = open_stack(path)
    logger.info("%s: %s Frames, %sx%s px", stem, stack.shape[0], stack.shape[1], stack.shape[2])
    table, events, _ = track_stack(stack, params)
    summary = tracking_summary(table, events); summary["stem"] = stem
    out_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(out_dir / f"{stem}_tracks.csv", index=False)
    (events if not events.empty else pd.DataFrame(columns=["frame", "type", "track_id", "other_track_id", "gap_frames"])
     ).to_csv(out_dir / f"{stem}_events.csv", index=False)
    with open(out_dir / f"{stem}_tracking_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    if write_zarr:
        save_stack(out_dir / f"tracks_{stem}.zarr", relabel_stack(stack, table))
    logger.info("%s: %d Objekte, %d Spuren, %.3f neue IDs je Objekt-Frame, Median %s Frames",
                stem, summary["n_objects"], summary["n_tracks"], summary["new_ids_per_object_frame"], summary["median_track_frames"])
    return table, summary


def merge_into_results(results: pd.DataFrame, tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Alte Ergebnistabelle + neue IDs. Join ueber (filename, frame, alte track_id = Label im Stack)."""
    new_cols = ["track_id", "parent_track_id", "parent_area_ratio", "link_type", "gap_frames", "overlap_prev_frac",
                "n_candidates", "axis_major", "axis_minor", "perimeter", "circularity"]
    parts = []
    for stem, tab in tables.items():
        t = tab[["frame", "label"] + [c for c in new_cols if c in tab.columns] + (["skeleton_px"] if "skeleton_px" in tab.columns else [])]
        t = t.rename(columns={"label": "track_id_v11"}); t.insert(0, "filename", stem)
        parts.append(t)
    new = pd.concat(parts, ignore_index=True)
    old = results.rename(columns={"track_id": "track_id_v11"})
    merged = old.merge(new, on=["filename", "frame", "track_id_v11"], how="left")
    n_missing = int(merged["track_id"].isna().sum())
    if n_missing:
        logger.warning("%d Zeilen der alten Tabelle haben kein Objekt im Stack (bleiben mit track_id NaN).", n_missing)
    return merged


def batch(processed_dir: Path, combined: Optional[Path], out_dir: Path, params: TrackParams, write_zarr: bool) -> None:
    stacks = sorted(processed_dir.glob("masks_*.zarr")) + sorted(processed_dir.glob("labels_*.zarr"))
    if not stacks:
        raise SystemExit(f"Keine masks_*.zarr / labels_*.zarr in {processed_dir}")
    tables, summaries = {}, []
    for s in stacks:
        table, summary = process_stack(s, out_dir, params, write_zarr=write_zarr)
        tables[stem_of(s)] = table; summaries.append(summary)
    pd.DataFrame(summaries).to_csv(out_dir / "tracking_summary.csv", index=False)
    if combined is not None:
        old = pd.read_parquet(combined) if combined.suffix == ".parquet" else pd.read_csv(combined)
        merged = merge_into_results(old, tables)
        target = combined.with_name(combined.stem + "_retracked" + combined.suffix)
        if combined.suffix == ".parquet":
            merged.to_parquet(target, index=False)
        else:
            merged.to_csv(target, index=False)
        logger.info("Geschrieben: %s (%d Zeilen, %d Spuren)", target, len(merged), merged["track_id"].nunique())


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stack", nargs="?", type=Path, help="ein masks_<stem>.zarr / labels_<stem>.zarr")
    ap.add_argument("--batch", type=Path, help="Ordner mit masks_*.zarr (02_processed eines Experiments)")
    ap.add_argument("--combined", type=Path, help="Combined_Results.csv des Experiments (fuer --batch)")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--no-zarr", action="store_true", help="keinen tracks_<stem>.zarr schreiben")
    ap.add_argument("--skeleton", action="store_true", help="Skelettlaenge je Objekt berechnen")
    for name, default in (("memory", 3), ("dist-frac", 1.5), ("dist-add", 10.0), ("max-ratio", 2.5),
                          ("sparse-n", 20), ("sparse-frac", 6.0), ("sparse-add", 20.0), ("touch-px", 2)):
        ap.add_argument(f"--{name}", type=type(default), default=default)
    a = ap.parse_args(argv)
    params = TrackParams(memory=a.memory, dist_frac=a.dist_frac, dist_add=a.dist_add, max_ratio=a.max_ratio,
                         sparse_n=a.sparse_n, sparse_frac=a.sparse_frac, sparse_add=a.sparse_add,
                         touch_px=a.touch_px, skeleton=a.skeleton)
    logger.info("Parameter: %s", asdict(params))
    if a.batch:
        batch(a.batch, a.combined, a.out, params, write_zarr=not a.no_zarr)
    elif a.stack:
        process_stack(a.stack, a.out, params, write_zarr=not a.no_zarr)
    else:
        ap.error("stack oder --batch angeben")


if __name__ == "__main__":
    main()
