# Why the tracks break, and what would fix them

Review of `cellpose_pipeline_v11.py` (segmentation + tracking) and a diagnosis on the manually QC'd batch
(WT, pH, 6 min: 11 chambers, 132 frames at 10 min, 66,912 object-frames). Reproduce the numbers with

```
python analyse_pipeline/diagnose_tracking.py <Combined_Results.csv> --qc <qc_exclusions.csv> --out <dir>
```

which writes `tracking_new_id_categories.csv`, `tracking_relink_prototype.csv` and, with `--qc`,
`tracking_manual_link_recall.csv`. The script only reads; it changes no pipeline output.

## 1. Findings in the Cellpose pipeline

| where | what | effect |
| --- | --- | --- |
| `CellTracker.update` (l. 117–158) | Only the previous frame is remembered. An object missed or filtered in one frame ends its track; when it reappears it gets a new ID. | 12 % of all new IDs are the same object one frame later (`dropout`), more after longer gaps. |
| `CellTracker.update` (l. 118–120) | A frame without ROIs resets the tracker; every cell gets a new ID afterwards. | 6 empty frames in one chamber of the batch, each cutting all tracks. |
| `_iou_threshold_for` (l. 99–106) | IoU ≥ 0.3 is a hard gate. Two equal circles reach IoU 0.3 at a shift of 0.85 radii, so a cell pushed by its neighbours by one radius in 10 min is lost. The `small_object_iou_floor` only applies below r = 12 px (area ≤ 452 px). | 10.8 % of new IDs are the same object, same size (area ratio 1.09), moved 1.2 radii, estimated IoU 0.16. |
| `_allowed_dist` (l. 87–97) | `max(adaptive, max_dist)` makes the size-adaptive distance a floor of 50 px, never a tighter bound; for normal cells (r ≈ 34 px) the gate is 50–70 px, about 1.5 radii. | 4.8 % of new IDs are the same object beyond 50 px. |
| `segment_and_filter_frame` (l. 386–398) | The ROI filters (border, min/max area, eccentricity, solidity) run before tracking. An object that fails one filter in one frame vanishes from the tracker. A mother with an attached bud, segmented as one object, has low solidity exactly at budding. | Objects appear without any change in their neighbours' masks: 28 % of new IDs are a similar-sized object next to an unchanged neighbour, 23 % a small one. Cannot be separated from missed segmentations or arriving cells without the raw masks. |
| `process_file` (l. 505–517, 562–579) | With `rotation_angle_deg` set: the chamber is detected on the rotated first frame, the crop is cut from the unrotated frames, the phase crop is then rotated around its own centre, and the fluorescence crops are never rotated. Masks and fluorescence would be misaligned. | Only when the config sets a rotation angle. Harmless when it is null. |
| `process_file` (l. 594–598) | The zarr stack stores track IDs of the filtered ROIs, not Cellpose's raw labels. | Re-tracking from the masks is possible (unique values per frame are objects), but the filtered-out objects are gone. |
| `extract_metadata_from_filename` (l. 280) | Comment says YYYYMMDD, regex takes six digits. | None; the filenames have six digits. |

The tracker's other parts are sound: the Hungarian assignment with cost 1 − IoU, the area-growth gate of 4×,
one tracker per file, the channel index mapping between `frame_data` and `measure_features_frame`.

## 2. What a new track ID was one frame earlier

| category | share of 9,042 new IDs |
| --- | --- |
| new: touching a similar-sized object (split or flicker) | 28.0 % |
| new: small object touching a larger one (bud or split-off piece) | 22.6 % |
| new: isolated (flushed in, or absent for more than one frame) | 18.4 % |
| same object, missing for one frame (dropout) | 12.3 % |
| same object, failed the IoU gate | 10.8 % |
| same object, moved more than 50 px | 4.8 % |
| same object, estimated IoU ≥ 0.3 but lost (competition, mask shape) | 3.0 % |
| after an empty frame | 0.1 % |

- 31 % of new IDs are an object that is still in the table within two frames: the tracker alone loses them.
- Of the "touching" appearances, the neighbour's mask lost the new object's area in only 4–10 %: these are not
  splits of one mask but objects that were absent from the table a frame earlier (filtered, missed, or arrived).
- 30 % of all tracks last one frame, 43 % at most two. Blips of this kind next to a mother count as budding
  events in `lineage.py`, which does not require the bud to persist (`bud_max_frames` is a report flag only).
- In sparse frames (≤ 20 objects) new IDs are as frequent as in dense ones (0.114 vs 0.139 per object-frame),
  but 57 % of them are isolated: cells arriving, or cells that jumped. The manual merges confirm the jumps: for
  links across one frame the centroid moved a median of 2.2 radii, the 90th percentile 8 radii, at equal area.
- Continued links sit close to the gate: median shift 0.3 radii, estimated IoU 0.65, and 23 % below IoU 0.4.

## 3. Prototype re-linker on the table alone

Centroid and area only, Hungarian assignment, memory of two frames, gate 1.5 r·√gap + 10 px and area ratio
≤ 2.5, no IoU; in frames with ≤ 20 objects a second pass links leftover pairs up to 6 r + 20 px when the pair is
mutually unique.

| | original | prototype |
| --- | --- | --- |
| new IDs per object-frame | 0.137 | 0.085 |
| dense frames (> 60 objects) | 0.139 | 0.078 |
| sparse frames (≤ 20 objects) | 0.114 | 0.091 |
| median track length (frames) | 3 | 5 |
| tracks ≥ 10 frames | 21 % | 37 % |
| tracks per chamber | 824 | 503 |
| manual merges reproduced (244) | | 23 %; 41 % of the one-frame gaps |

The remaining manual merges are jumps over several radii and more than two frames (92 of 244 have a gap above
five frames; 20 overlap in time, i.e. one cell segmented as two). A linker without images cannot make those
safely, and nothing on the table can repair a merged or missing mask.

## 4. Ways to fix it, by cost

1. **Analysis side, now.** Replace the conservative gap closing in `relink.py` with the prototype (whole movie,
   memory, ambiguity-aware long-range links in sparse frames) and require a bud to persist two frames before it
   counts as an event. Gains: fragmentation −37 %, blips out of the event count. Limits: filtered and merged
   masks stay lost; jumps over several frames stay unlinked.
2. **Re-track from the saved masks** (`masks_<stem>.zarr` in `io.scratch_dir`), if they still exist. Per frame
   the unique label values are the objects, so a new tracker can run without Cellpose: overlap as cost, not as
   gate; memory of three frames; merge and split handled by overlap with two predecessors or two successors;
   a bud gets a parent link from its first mask touching the mother. This replaces the heuristic mother-bud
   matching in `lineage.py` by a measured link. Limits: the filtered-out objects are not in the masks.
3. **Re-run the Cellpose pipeline** with: tracking before filtering (filters become flags in the table), raw
   label masks saved next to the track masks, the tracker of option 2, an empty frame that does not reset, the
   rotation applied to the full frame before chamber detection and to every channel. This is the only route
   that removes the filter flicker. Costs GPU time for ~565 movies.

Open questions that decide between 2 and 3: does the config set `rotation_angle_deg`; which ROI filters are
active (`min_solidity`, `max_eccentricity`, `min_area_px`); do the scratch zarr stacks still exist on the
cluster.
