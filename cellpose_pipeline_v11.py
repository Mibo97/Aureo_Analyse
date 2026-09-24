import yaml
import logging
import time
from PIL import Image, ImageDraw, ImageFont
import re
import numpy as np
import pandas as pd
from pathlib import Path
import nd2
from cellpose import models
from skimage import exposure, filters, morphology, measure, registration, restoration, transform
from skimage.color import label2rgb
from skimage.io import imsave
from scipy.optimize import linear_sum_assignment
import zarr
import warnings
from typing import Dict, List, Optional, Tuple

warnings.filterwarnings('ignore')

# Configure logging for cluster compatibility
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger(__name__)

CONFIG_VERSION = "1.4"

# ==============================================================================
# TRACKING LOGIC
# ==============================================================================
def compute_iou(mask_a: np.ndarray, mask_b: np.ndarray) -> float:
    """IoU zwischen zwei booleschen Masken."""
    intersection = np.logical_and(mask_a, mask_b).sum()
    union = np.logical_or(mask_a, mask_b).sum()
    return intersection / union if union > 0 else 0.0

def centroid_distance(roi_a: Dict, roi_b: Dict) -> float:
    """Euklidische Distanz zwischen zwei ROI-Centroids."""
    cy_a, cx_a = roi_a['centroid']
    cy_b, cx_b = roi_b['centroid']
    return np.sqrt((cy_a - cy_b) ** 2 + (cx_a - cx_b) ** 2)

class CellTracker:
    """
    Frame-by-frame Tracker mit IoU + Centroid-Distanz und Hungarian Algorithm.

    Die Distanz- und Wachstumstoleranz werden GRÖSSENADAPTIV berechnet:
    kleine Objekte (Knospen) dürfen sich relativ zu ihrer eigenen Größe
    stärker bewegen und schneller wachsen als große, stabile Mutterzellen.
    Das verhindert, dass schnell wachsende/leicht driftende Knospen vom
    Tracker als "verloren" behandelt und fälschlich neu gezählt werden.
    """
    def __init__(
        self,
        iou_threshold: float = 0.3,
        max_dist: float = 50.0,
        size_adaptive: bool = True,
        max_dist_frac: float = 1.5,
        max_area_growth: float = 4.0,
        small_object_radius_px: float = 12.0,
        small_object_iou_floor: float = 0.05,
    ):
        self.iou_threshold = iou_threshold
        self.max_dist = max_dist
        # Größenadaptiver Modus: max_dist wird zusätzlich relativ zum
        # äquivalenten Radius des kleineren der beiden ROIs skaliert.
        self.size_adaptive = size_adaptive
        self.max_dist_frac = max_dist_frac    # Vielfaches des Radius, das ein Objekt sich bewegen darf
        self.max_area_growth = max_area_growth  # max. erlaubtes Flächenverhältnis (groß/klein) zw. Frames
        # Unterhalb dieses Radius gilt ein Objekt als "klein" (z.B. Knospe):
        # IoU reagiert bei kleinen Masken extrem empfindlich auf 1-2px
        # Subpixel-Verschiebung, daher wird die IoU-Schwelle dort abgesenkt
        # und stattdessen primär über Distanz + Wachstumsplausibilität geprüft.
        self.small_object_radius_px = small_object_radius_px
        self.small_object_iou_floor = small_object_iou_floor
        self.next_id = 1
        self.prev_rois: List[Dict] = []

    @staticmethod
    def _equiv_radius(area: float) -> float:
        """Radius eines flächengleichen Kreises - robuster Größenproxy als bbox."""
        return np.sqrt(area / np.pi)

    def _allowed_dist(self, p_roi: Dict, c_roi: Dict) -> float:
        """Erlaubte Centroid-Distanz für dieses ROI-Paar."""
        if not self.size_adaptive:
            return self.max_dist
        # Kleineres der beiden Objekte ist der limitierende Faktor:
        # eine winzige Knospe neben einer riesigen Mutterzelle soll an
        # ihrer eigenen Größe gemessen werden, nicht an der Mutterzelle.
        r_min = min(self._equiv_radius(p_roi['area']), self._equiv_radius(c_roi['area']))
        # Untere Grenze, damit sehr kleine ROIs nicht eine winzige Toleranz von <1px bekommen
        adaptive = max(r_min * self.max_dist_frac, 3.0)
        return max(adaptive, self.max_dist)

    def _iou_threshold_for(self, p_roi: Dict, c_roi: Dict) -> float:
        """Senkt die IoU-Schwelle für kleine Objekte ab."""
        if not self.size_adaptive:
            return self.iou_threshold
        r_min = min(self._equiv_radius(p_roi['area']), self._equiv_radius(c_roi['area']))
        if r_min <= self.small_object_radius_px:
            return self.small_object_iou_floor
        return self.iou_threshold

    def _area_growth_ok(self, p_roi: Dict, c_roi: Dict) -> bool:
        """Erlaubt schnelles Wachstum (z.B. Knospen), begrenzt aber Sprünge,
        die eher auf einen Identity-Switch hindeuten als auf echtes Wachstum."""
        a_prev, a_curr = p_roi['area'], c_roi['area']
        if a_prev <= 0 or a_curr <= 0:
            return False
        ratio = max(a_prev, a_curr) / min(a_prev, a_curr)
        return ratio <= self.max_area_growth

    def update(self, curr_rois: List[Dict]) -> List[Dict]:
        if not curr_rois:
            self.prev_rois = []
            return []

        if not self.prev_rois:
            for roi in curr_rois:
                roi['track_id'] = self.next_id
                self.next_id += 1
            self.prev_rois = curr_rois
            return curr_rois

        N, M = len(self.prev_rois), len(curr_rois)
        costs = np.ones((N, M))

        for i, p_roi in enumerate(self.prev_rois):
            for j, c_roi in enumerate(curr_rois):
                dist = centroid_distance(p_roi, c_roi)
                if dist > self._allowed_dist(p_roi, c_roi):
                    continue
                if not self._area_growth_ok(p_roi, c_roi):
                    continue
                iou = compute_iou(p_roi['mask'], c_roi['mask'])
                if iou < self._iou_threshold_for(p_roi, c_roi):
                    continue
                costs[i, j] = 1.0 - iou

        row_ind, col_ind = linear_sum_assignment(costs)
        assigned_curr = set()

        for r, c in zip(row_ind, col_ind):
            if costs[r, c] < 1.0:
                curr_rois[c]['track_id'] = self.prev_rois[r]['track_id']
                assigned_curr.add(c)

        for j, roi in enumerate(curr_rois):
            if j not in assigned_curr:
                roi['track_id'] = self.next_id
                self.next_id += 1

        self.prev_rois = curr_rois
        return curr_rois

    def reset(self):
        self.next_id = 1
        self.prev_rois = []

# ==============================================================================
# MAIN PIPELINE CLASS
# ==============================================================================
class MicroscopyPipeline:
    def __init__(self, config_path: str):
        self.config = self._load_config(config_path)
        self._setup_directories()
        self.rotation_angle = self.config['preprocessing'].get('rotation_angle_deg')
        self.channel_map = self._build_channel_map()

        # Rolling-Ball-Hintergrundsubtraktion für die Fluoreszenz-Messkanäle -
        # OPT-IN über measurements.rolling_ball_radius_px in der YAML-Config,
        # damit bestehende Configs sich nicht stillschweigend im Verhalten
        # ändern. None = deaktiviert, Kanäle werden roh gemessen (bisheriges
        # Verhalten). Betrifft NUR die Messkanäle, NICHT phase_contrast (das
        # dient ausschließlich der Segmentierung, nicht der Intensitätsmessung).
        self.rolling_ball_radius_px = self.config['measurements'].get('rolling_ball_radius_px')
        if self.rolling_ball_radius_px is not None:
            logger.info(
                f"Rolling-Ball-Hintergrundsubtraktion AKTIV (radius={self.rolling_ball_radius_px}px). "
                f"Radius sollte deutlich größer als der Zellradius sein, sonst wird echtes "
                f"Zellsignal mit als Hintergrund geschätzt."
            )
        else:
            logger.info(
                "Keine Rolling-Ball-Hintergrundsubtraktion konfiguriert "
                "(measurements.rolling_ball_radius_px fehlt in der YAML) - "
                "Fluoreszenzkanäle werden roh (unsubtrahiert) gemessen."
            )

        cfg_version = self.config.get('config_version', 'unset')
        logger.info(f"Pipeline init | config_version={cfg_version} | expected={CONFIG_VERSION}")

    def _subtract_background(self, frame: np.ndarray) -> np.ndarray:
        """
        Rolling-Ball-Hintergrundsubtraktion für EINEN Fluoreszenz-Frame.

        Schätzt den (räumlich variierenden) Hintergrund, indem eine "Kugel"
        mit Radius `self.rolling_ball_radius_px` unter das Intensitätsprofil
        gerollt wird, und zieht ihn ab. Im Unterschied zu einer pauschalen
        Offset-Korrektur in der Analyse (sensors.py) erfasst das auch
        ungleichmäßige Beleuchtung/Autofluoreszenz-Gradienten über die
        Kammer, nicht nur einen globalen Kamera-Offset.

        Ergebnis wird bei 0 geclippt (negative Intensitäten sind physikalisch
        nicht sinnvoll).

        PERFORMANCE: Die Hintergrundschätzung erfolgt auf einem auf
        `rolling_ball_downscale_size` px (längste Kante) verkleinerten Bild
        und wird anschließend hochskaliert. Der Hintergrund ist per Definition
        niederfrequent/glatt - Downscaling verliert dort kaum Information,
        spart aber ~80x Rechenzeit (1024x1024, r=50: ~16s -> ~0.2s).
        Ueber measurements.rolling_ball_downscale_size in der YAML konfigurierbar
        (Default 256). Auf 0 setzen, um Downscaling zu deaktivieren.

        WICHTIG: Nach dieser Subtraktion ist der Rauschboden der gemessenen
        Werte viel niedriger als bei Rohwerten - sensors.RatioSensorConfig.
        min_denominator (siehe dortiger Docstring) muss entsprechend
        NEU kalibriert werden (deutlich kleinerer Wert als bisher), sonst
        werden nach der Subtraktion zu viele echte Messwerte faelschlich als
        Rauschen zu NaN.
        """
        frame_f = frame.astype(np.float32)  # float32 reicht, float64 verdoppelt RAM + Zeit
        radius = self.rolling_ball_radius_px
        downscale_size = self.config['measurements'].get('rolling_ball_downscale_size', 256)

        h, w = frame_f.shape
        longest_edge = max(h, w)

        if downscale_size <= 0 or longest_edge <= downscale_size:
            # Bild bereits klein genug, oder Downscaling explizit deaktiviert.
            background = restoration.rolling_ball(frame_f, radius=radius)
        else:
            scale = downscale_size / longest_edge
            small_shape = (max(1, int(round(h * scale))), max(1, int(round(w * scale))))

            small = transform.resize(
                frame_f, small_shape, preserve_range=True, anti_aliasing=True
            ).astype(np.float32)

            # Radius mitskalieren: sonst deckt die Kugel im verkleinerten Bild
            # relativ zur Bildgroesse einen falschen Bereich ab.
            radius_small = max(radius * scale, 1.0)

            bg_small = restoration.rolling_ball(small, radius=radius_small)

            background = transform.resize(
                bg_small, frame_f.shape, preserve_range=True
            ).astype(np.float32)

        return np.clip(frame_f - background, 0, None)

    def _load_config(self, path: str) -> dict:
        with open(path, 'r') as f:
            cfg = yaml.safe_load(f)
        required = ['experiment', 'io', 'preprocessing', 'segmentation', 'roi_filter', 'tracking', 'channels', 'measurements']
        for key in required:
            if key not in cfg:
                raise ValueError(f"Missing required YAML section: {key}")
        return cfg

    def _setup_directories(self):
        io = self.config['io']
        self.input_dir = Path(io['input_dir'])
        self.output_dir = Path(io['output_dir'])
        self.scratch_dir = Path(io['scratch_dir'])
        for d in [self.output_dir / 'QC', self.scratch_dir]:
            d.mkdir(parents=True, exist_ok=True)

    def _build_channel_map(self) -> dict:
        return {ch['name']: {'index': ch['index'], 'measure': ch['measure_intensity']} for ch in self.config['channels']}

    def extract_metadata_from_filename(self, filepath: str) -> Dict[str, str]:
        """Extrahiert Datum, Bedingung, Replikat und Kammer aus dem Dateinamen."""
        stem = Path(filepath).stem
        # Erwartet: YYYYMMDD_Condition_Rep_Cham (z.B. 240612_Osc10min_Rep1_ChamA)
        match = re.match(r"^(\d{6})_(.+?)_(Rep\d+)_(Cham[A-Za-z0-9]+)", stem)
        if match:
            return {'date': match.group(1),
		    'condition': match.group(2),
                    'replicate': match.group(3),
                    'chamber': match.group(4)}
        logger.warning(f"Konnte Metadaten nicht aus '{stem}' extrahieren.")
        return {'date': 'Unknown', 'condition': 'Unknown', 'replicate': 'Unknown', 'chamber': 'Unknown'}

    def detect_chamber(self, frame: np.ndarray) -> Optional[Tuple[int, int, int, int]]:
        """Robuste iterative Kammer-Erkennung (aus dem Notebook)."""
        cfg = self.config['preprocessing']['chamber_detection']
        img = exposure.rescale_intensity(frame.astype(float), out_range=(0.0, 1.0))
        img = filters.gaussian(img, sigma=1)
        img = 1.0 - img  # Invertieren wie im Fiji-Makro

        base_thresh = filters.threshold_otsu(img)
        step = 25 / 255.0
        min_area = cfg['min_area_px']
        shrink = cfg.get('shrink_px', 10)

        for iteration in range(cfg.get('max_iterations', 100)):
            thresh = base_thresh + iteration * step
            if thresh >= 1.0:
                break
            binary = img > thresh
            labeled = measure.label(binary)
            props = [p for p in measure.regionprops(labeled) if p.area >= min_area]

            if props:
                best = max(props, key=lambda p: p.area)
                r0, c0, r1, c1 = best.bbox
                # Enlarge nach innen
                r0 += shrink; c0 += shrink
                r1 -= shrink; c1 -= shrink
                logger.info(f"Kammer gefunden (Iter {iteration}): y={r0}:{r1}, x={c0}:{c1}")
                return (r0, c0, r1, c1)
        return None

    # Modelle die CellposeModel statt Cellpose brauchen und keinen diameter-Parameter kennen
    _NEW_MODELS = {'cpsam', 'cpsam_v2', 'cpdino', 'cpdino-vitb'}

    def _load_model(self, seg_cfg: dict):
        """
        Laedt das richtige Cellpose-Modell je nach Config.

        Alte Modelle (cyto2, cyto3, bact_phase_omni, ...):
            models.Cellpose(model_type=...) – benoetigt diameter + channels + flow_threshold

        Neue SAM/DINO-Modelle (cpsam, cpsam_v2, cpdino, cpdino-vitb):
            models.CellposeModel(pretrained_model=...) – kein diameter, kein flow_threshold
            Gewichte werden beim ersten Aufruf automatisch von HuggingFace geladen.
        """
        model_type = seg_cfg['model_type']
        use_gpu    = seg_cfg.get('use_gpu', False)
        is_new     = model_type in self._NEW_MODELS

        if is_new:
            logger.info(f"Lade neues SAM/DINO-Modell: {model_type} (Gewichte ggf. von HuggingFace)")
            model = models.CellposeModel(pretrained_model=model_type, gpu=use_gpu)
        else:
            logger.info(f"Lade klassisches Cellpose-Modell: {model_type}")
            model = models.Cellpose(model_type=model_type, gpu=use_gpu)

        # Flag fuer segment_and_filter_frame damit der richtige eval()-Aufruf gewaehlt wird
        model._is_new_model = is_new
        return model

    def segment_and_filter_frame(self, phase_frame: np.ndarray, model) -> List[Dict]:
        """Segmentiert einen Frame und filtert ROIs sofort."""
        seg_cfg = self.config['segmentation']
        filt_cfg = self.config['roi_filter']

        norm = exposure.rescale_intensity(phase_frame.astype(float), out_range=(0.0, 1.0))

        # HINWEIS: Auch die SAM/DINO-Modelle (cpsam, cpsam_v2, cpdino, cpdino-vitb)
        # nutzen weiterhin Flow-Felder zur Maskenberechnung (CellposeModel.eval()
        # akzeptiert flow_threshold ebenso wie cellprob_threshold). Nur `diameter`
        # ist bei diesen Modellen ohne Wirkung, da die Groesse automatisch bestimmt
        # wird. flow_threshold daher bei BEIDEN Modell-Typen uebergeben.
        is_new_model = getattr(model, '_is_new_model', False)
        niter = seg_cfg.get('niter')  # optional: mehr Dynamics-Iterationen helfen bei eingeschnuerten/laenglichen Formen
        if is_new_model:
            masks, _, _ = model.eval(
                norm,
                cellprob_threshold=seg_cfg['cellprob_threshold'],
                flow_threshold=seg_cfg.get('flow_threshold', 0.4),
                min_size=seg_cfg['min_size_px'],
                niter=niter,
            )
        else:
            masks, _, _, _ = model.eval(
                norm, diameter=seg_cfg['diameter_px'], channels=[0, 0],
                cellprob_threshold=seg_cfg['cellprob_threshold'],
                flow_threshold=seg_cfg['flow_threshold'],
                min_size=seg_cfg['min_size_px'],
                niter=niter,
            )

        if masks.max() == 0:
            return []

        H, W = masks.shape
        bnd = filt_cfg['exclude_boundary_px']
        valid_rois = []

        for p in measure.regionprops(masks):
            r0, c0, r1, c1 = p.bbox
            # Rand-Check
            if r0 <= bnd or c0 <= bnd or (H - r1) <= bnd or (W - c1) <= bnd:
                continue
            # Größen-Check
            if not (filt_cfg['min_area_px'] <= p.area <= filt_cfg['max_area_px']):
                continue
            # Form-Check (optional)
            if filt_cfg.get('max_eccentricity') and p.eccentricity > filt_cfg['max_eccentricity']:
                continue
            if filt_cfg.get('min_solidity') and p.solidity < filt_cfg['min_solidity']:
                continue

            valid_rois.append({
                'roi_id':      p.label,
                'area':        p.area,
                'centroid':    p.centroid,
                'bbox':        p.bbox,
                'mask':        (masks == p.label),
                'solidity':    p.solidity,
                'eccentricity': p.eccentricity,
            })
        return valid_rois

    def measure_features_frame(self, valid_rois: List[Dict], hyperstack_frame: np.ndarray, t: int) -> List[Dict]:
        """Misst Intensitäten für die getrackten ROIs eines einzelnen Frames.

        hyperstack_frame Layout: [0]=phase_contrast, [1..N]=measure-Kanäle in Config-Reihenfolge.
        """
        records = []
        meas_cfg = self.config['measurements']

        # Baue Mapping: ch_name -> Index in hyperstack_frame
        # frame_array wurde aufgebaut als: [phase] + [ch für ch mit measure=True, ohne phase]
        arr_idx = 1  # 0 ist phase_contrast
        frame_idx_map = {}
        for ch_name, ch_info in self.channel_map.items():
            if ch_name == 'phase_contrast':
                continue
            if ch_info['measure']:
                frame_idx_map[ch_name] = arr_idx
                arr_idx += 1

        for roi in valid_rois:
            entry = {
                'frame':        t,
                'track_id':     roi.get('track_id', -1),
                'roi_id':       roi['roi_id'],
                'area':         roi['area'],
                'centroid_y':   roi['centroid'][0],
                'centroid_x':   roi['centroid'][1],
                'solidity':     roi.get('solidity', np.nan),
                'eccentricity': roi.get('eccentricity', np.nan),
            }
            mask = roi['mask']

            for ch_name, arr_position in frame_idx_map.items():
                if arr_position < hyperstack_frame.shape[0]:
                    vals = hyperstack_frame[arr_position][mask]
                    if len(vals) > 0:
                        for stat in meas_cfg.get('intensity_stats', ['mean', 'std']):
                            entry[f'{stat}_{ch_name}'] = float(getattr(np, stat)(vals))
                    else:
                        for stat in meas_cfg.get('intensity_stats', ['mean', 'std']):
                            entry[f'{stat}_{ch_name}'] = np.nan
            records.append(entry)
        return records

    def process_file(self, filepath: str) -> Optional[pd.DataFrame]:
        t0 = time.time()
        stem = Path(filepath).stem
        metadata = self.extract_metadata_from_filename(filepath)
        logger.info(f"=== Processing: {Path(filepath).name} | Cond: {metadata['condition']} ===")

        # 1. Bild öffnen mit nd2 (Lazy Loading via Dask)
        with nd2.ND2File(filepath) as nd2_file:
            # Dimensionen extrahieren
            dims = nd2_file.sizes  # dict z.B. {'T': 100, 'C': 2, 'Y': 1024, 'X': 1024, 'Z': 1}
            T = dims.get('T', 1)
            C = dims.get('C', 1)
            Z = dims.get('Z', 1)
            H = dims.get('Y', 1)
            W = dims.get('X', 1)

            phase_idx = self.channel_map['phase_contrast']['index']

            # Lazy Dask-Array laden (speicherschonend!)
            dask_stack = nd2_file.to_dask()

            # Sicherheitscheck: Shape des Dask-Arrays muss zur Reihenfolge
            # von dims.values() passen, sonst werden Achsen in get_frame
            # stillschweigend falsch zugeordnet (z.B. T/C vertauscht).
            assert dask_stack.shape == tuple(dims.values()), (
                f"Shape mismatch: dask_stack.shape={dask_stack.shape} "
                f"vs. dims={dims}"
            )

            # Achsenreihenfolge bestimmen (nd2 gibt sie in der Reihenfolge der dims)
            # Typisch: ('T', 'C', 'Z', 'Y', 'X') oder ähnlich
            axis_order = tuple(dims.keys())

            # Helper-Funktion um Frame zu extrahieren
            def get_frame(t: int, c: int, z: int = 0) -> np.ndarray:
                """Extrahiert einen einzelnen Frame als numpy array."""
                # Baue Index-Tuple basierend auf Achsenreihenfolge
                idx = []
                for ax in axis_order:
                    if ax == 'T':
                        idx.append(t)
                    elif ax == 'C':
                        idx.append(c)
                    elif ax == 'Z':
                        idx.append(z)
                    elif ax in ('Y', 'X'):
                        idx.append(slice(None))  # alle Pixel
                return dask_stack[tuple(idx)].compute()

            # 2. Preprocessing (Rotation) auf den ersten Frame anwenden
            ref_frame_phase = get_frame(0, phase_idx, 0)
            if self.rotation_angle is not None:
                ref_frame_phase = transform.rotate(
                    ref_frame_phase, self.rotation_angle,
                    resize=False, preserve_range=True
                )

            # 3. Kammer erkennen
            chamber_roi = self.detect_chamber(ref_frame_phase)
            if chamber_roi is None:
                if self.config['preprocessing']['chamber_detection'].get('fallback_to_full_image', True):
                    logger.warning(f"{stem}: Kammer nicht gefunden. Nutze gesamtes Bild.")
                    chamber_roi = (0, 0, H, W)
                else:
                    logger.error(f"{stem}: Kammer-Erkennung fehlgeschlagen. Überspringe Datei.")
                    return None

            minr, minc, maxr, maxc = chamber_roi
            crop_H, crop_W = maxr - minr, maxc - minc

            # 4. Cellpose Modell laden
            seg_cfg = self.config['segmentation']
            model = self._load_model(seg_cfg)

            # 5. Frame-by-Frame Verarbeitung
            tracking_cfg = self.config['tracking']
            tracker = CellTracker(
                iou_threshold=tracking_cfg['iou_threshold'],
                max_dist=tracking_cfg['max_dist'],
                size_adaptive=tracking_cfg.get('size_adaptive', True),
                max_dist_frac=tracking_cfg.get('max_dist_frac', 1.5),
                max_area_growth=tracking_cfg.get('max_area_growth', 4.0),
                small_object_radius_px=tracking_cfg.get('small_object_radius_px', 12.0),
                small_object_iou_floor=tracking_cfg.get('small_object_iou_floor', 0.05),
            )

            all_records = []

            # NEU: Ein einziges 3D-Array für den ganzen Stack (T, Y, X)
            masks_zarr_path = self.scratch_dir / f"masks_{stem}.zarr"
            masks_zarr = zarr.open(
                str(masks_zarr_path),
                mode='w',
                shape=(T, crop_H, crop_W),
                chunks=(1, crop_H, crop_W),  # Ein Frame pro Chunk
                dtype=np.uint16
            )

            qc_overlays = []
            logger.info(f"Starte Frame-Verarbeitung: {T} Frames | Modell: {seg_cfg['model_type']}")
            t_start_batch = time.perf_counter()

            for t in range(T):
                t_frame_start = time.perf_counter()

                # Phase-Kanal laden
                phase_idx_ch = self.channel_map['phase_contrast']['index']
                phase_frame = get_frame(t, phase_idx_ch, 0)[minr:maxr, minc:maxc]
                if self.rotation_angle is not None:
                    phase_frame = transform.rotate(
                        phase_frame, self.rotation_angle,
                        resize=False, preserve_range=True
                    )

                # Nur Kanäle mit measure_intensity laden
                frame_data = [phase_frame]
                for ch_name, ch_info in self.channel_map.items():
                    if ch_name == 'phase_contrast':
                        continue
                    if ch_info['measure']:
                        ch_idx = ch_info['index']
                        crop_data = get_frame(t, ch_idx, 0)[minr:maxr, minc:maxc]
                        if self.rolling_ball_radius_px is not None:
                            crop_data = self._subtract_background(crop_data)
                        frame_data.append(crop_data)

                frame_array = np.array(frame_data)

                # Segmentierung & Filterung
                valid_rois = self.segment_and_filter_frame(phase_frame, model)

                # Tracking
                tracked_rois = tracker.update(valid_rois)

                # Messung
                if tracked_rois:
                    frame_records = self.measure_features_frame(tracked_rois, frame_array, t)
                    all_records.extend(frame_records)

                # Masken speichern
                mask_img = np.zeros((crop_H, crop_W), dtype=np.uint16)
                for roi in tracked_rois:
                    mask_img[roi['mask']] = roi['track_id']
                masks_zarr[t, :, :] = mask_img

                # Fortschritts-Log
                t_frame_end = time.perf_counter()
                elapsed = t_frame_end - t_start_batch
                per_frame = elapsed / (t + 1)
                remaining = per_frame * (T - t - 1)
                n_rois = len(tracked_rois)
                logger.info(
                    f"Frame {t+1:>3}/{T} | ROIs: {n_rois:>3} | "
                    f"Frame-Zeit: {t_frame_end - t_frame_start:.1f}s | "
                    f"Verbleibend: ~{remaining/60:.1f} min"
                )

                # QC Overlay
                try:
                    phase_norm = exposure.rescale_intensity(
                        phase_frame.astype(float), out_range=(0, 255)
                    ).astype(np.uint8)
                    phase_rgb = np.stack([phase_norm] * 3, axis=-1)

                    overlay = label2rgb(
                        mask_img, image=phase_rgb, bg_label=0,
                        alpha=0.4, image_alpha=1.0, bg_color=(0, 0, 0)
                    )
                    overlay_uint8 = (overlay * 255).astype(np.uint8)

                    pil_img = Image.fromarray(overlay_uint8)
                    draw = ImageDraw.Draw(pil_img)
                    for roi in tracked_rois:
                        if 'track_id' in roi:
                            cy, cx = roi['centroid']
                            tid = roi['track_id']
                            draw.text((int(cx) + 1, int(cy) + 1), str(tid), fill=(0, 0, 0))
                            draw.text((int(cx), int(cy)), str(tid), fill=(255, 255, 255))
                    overlay_uint8 = np.array(pil_img)
                    qc_overlays.append(overlay_uint8)

                except Exception as e:
                    logger.warning(f"QC Overlay für Frame {t} fehlgeschlagen: {e}")

        # QC-Stack speichern (außerhalb des `with`-Blocks!)
        if qc_overlays:
            try:
                import tifffile as tifffile_mod
                qc_stack = np.stack(qc_overlays, axis=0)
                qc_path = self.output_dir / 'QC' / f"{stem}_QC_overlay.tif"
                tifffile_mod.imwrite(
                    str(qc_path), qc_stack,
                    compression='lzw',
                    imagej=True,
                    metadata={'axes': 'TYXS'},
                    photometric='rgb',
                )
                logger.info(
                    f"QC-Stack gespeichert: {qc_path.name} | {qc_stack.shape[0]} Frames, "
                    f"{qc_stack.nbytes // 1024 // 1024} MB"
                )
            except Exception as e:
                logger.error(f"QC-Stack Speicherung fehlgeschlagen: {e}")

        # 6. Ergebnisse zusammenbauen
        if not all_records:
            logger.warning(f"{stem}: Keine Zellen gefunden.")
            return None

        df = pd.DataFrame(all_records)

        for key, val in metadata.items():
            df[key] = val
        df['filename'] = stem

        # 7. Speichern
        out_fmt = self.config['io']['export_results_as']
        res_path = self.output_dir / f"Single-Cell-Results_{stem}.{out_fmt}"

        if out_fmt == 'parquet':
            df.to_parquet(res_path, index=False)
        else:
            df.to_csv(res_path, index=False)

        elapsed = time.time() - t0
        logger.info(
            f"{stem}: Fertig in {elapsed:.1f}s | {len(df)} Zeilen, "
            f"{df['track_id'].nunique()} Tracks gespeichert."
        )
        return df

    def run_batch(self) -> None:
        files = sorted(list(self.input_dir.glob(self.config['io']['file_pattern'])))
        # Filtere Overview-Dateien heraus, falls sie im selben Ordner liegen
        files = [f for f in files if 'overview' not in f.name.lower()]

        logger.info(f"Gefundene Dateien zur Verarbeitung: {len(files)}")
        all_results = []

        for f in files:
            try:
                df = self.process_file(str(f))
                if df is not None:
                    all_results.append(df)
            except Exception as e:
                logger.error(f"Fehler bei {f.name}: {e}", exc_info=True)
                continue

        if all_results:
            combined = pd.concat(all_results, ignore_index=True)
            out_fmt = self.config['io']['export_results_as']
            out_path = self.output_dir / f"Combined_Results.{out_fmt}"

            if out_fmt == 'parquet':
                combined.to_parquet(out_path, index=False)
            else:
                combined.to_csv(out_path, index=False)

            logger.info(f"Batch abgeschlossen: {len(combined)} Gesamtzeilen, {combined['track_id'].nunique()} eindeutige Tracks.")
        else:
            logger.warning("Keine Ergebnisse produziert.")

# ==============================================================================
# Budding-Analyse: siehe separates Modul budding_analysis.py
# Dieses ließt die hier erzeugten Combined_Results/Single-Cell-Results und
# kann unabhängig von dieser Bildverarbeitungspipeline mit unterschiedlichen
# Parametern (mehrfach) durchlaufen werden.
# ==============================================================================


if __name__ == "__main__":
    import sys
    config_path = sys.argv[1] if len(sys.argv) > 1 else "pipeline_config_1.2_cpsam.yaml"
    pipeline = MicroscopyPipeline(config_path)
    pipeline.run_batch()

    out_fmt = pipeline.config['io']['export_results_as']
    combined_path = pipeline.output_dir / f"Combined_Results.{out_fmt}"
    if combined_path.exists():
        logger.info(
            f"Bildverarbeitung abgeschlossen. Fuer die Budding-Analyse "
            f"separat ausfuehren:\n"
            f"  python budding_analysis.py {combined_path}"
        )
    else:
        logger.warning(f"Keine Combined_Results gefunden unter {combined_path}")
