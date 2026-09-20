"""
lineage.py
==========
Mutter/Bud-Klassifikation und Budding-Ratio-Berechnung OHNE explizites
Lineage-Tracking aus der Bildverarbeitung.

HINTERGRUND
-----------
Die Cellpose-Pipeline liefert pro Zelle nur eine zeitlich stabile track_id,
aber KEINE Mutter-Tochter-Beziehung (kein echtes Lineage-Tracking). Da reife
Buds in der Mikrofluidik-Kammer typischerweise weggespült werden statt zu
bleiben, nutzen wir das als Signal:

  MUTTERZELLE: ein Track, der über eine MINDESTANZAHL Frames hinweg sichtbar
               bleibt (Standard: >= 10 Frames - konfigurierbar über
               mother_min_frames).

  BUD-KANDIDAT: ein Track, der NEU auftaucht (nicht im allerersten Frame der
               Kammer vorhanden) und deutlich KÜRZER bleibt als eine
               Mutterzelle (Standard: <= 5 Frames - konfigurierbar über
               bud_max_frames).

  ZUORDNUNG Bud -> Mutter: räumliche Nähe (Centroid-Distanz) zum Zeitpunkt
               des ERSTEN Auftretens des Bud-Kandidaten, ABER nur unter
               Tracks, die selbst als Mutterzelle qualifizieren. Das
               vermeidet den Fehlerfall "zwei kurzlebige Zellen berühren
               sich zufällig" (beide wären gar keine Mutter-Kandidaten).

               Der Suchradius ist GRÖSSENADAPTIV statt eines festen Pixel-
               werts: radius = sqrt(mother_area / pi) + tolerance_px, d.h.
               der geschätzte Mutterzell-Radius plus ein Toleranz-Puffer.
               Eine große Mutterzelle darf also einen weiter entfernt
               auftauchenden Bud "claimen" als eine kleine - biologisch
               sinnvoller als ein für alle Zellgrößen fixer Radius.

WICHTIG - das ist eine HEURISTIK, kein Ground-Truth-Lineage-Tracking:
- Falsch-positive sind möglich, wenn zwei Mutterzellen nah beieinander
  liegen und ein Bud fälschlich der falschen Mutter zugeordnet wird
  (bei mehreren Kandidaten im Radius wird die NÄCHSTGELEGENE gewählt).
- Falsch-negative sind möglich, wenn ein Bud außerhalb des adaptiven
  Radius auftaucht (z.B. Tracking-Sprung).
- Mütter mit area=NaN (z.B. Segmentierungsfehler) bekommen dadurch einen
  NaN-Radius und werden für KEINEN Bud als Kandidat berücksichtigt - das
  passiert still, ohne Warnung (laut Rücksprache: kommt in der Praxis
  praktisch nicht vor; falls sich das ändert, hier nachrüsten).
- Alle Schwellenwerte sind Parameter, NICHT fest verdrahtet - bitte an
  echten QC-Overlays kalibrieren (siehe `inspect_classification()` unten),
  bevor die Ergebnisse für eine Publikation verwendet werden.

FIX (Entkopplung von Erkennung und Qualitätsfilter):
Frühere Version filterte Bud-KANDIDATEN bereits über `bud_max_frames`
(finale, GESAMTE Tracklänge über die ganze Beobachtungsdauer). Das führte
dazu, dass Buds, die NICHT weggespült wurden (und z.B. selbst zu einer
Mutter heranwuchsen), strukturell NIE als Budding-Event erkannt wurden -
ein Bud mit finaler Tracklänge zwischen bud_max_frames und
mother_min_frames fiel komplett durchs Raster, und ein Bud, der selbst
>= mother_min_frames erreichte, wurde als eigenständige, scheinbar
unverwandte Mutter geführt (das reale Reproduktionsereignis, das sie mit
ihrer echten Mutter verbindet, ging verloren). Das erzeugte einen
systematischen (nicht zufälligen) Bias in Budding Ratio und µ, der mit der
Auswaschrate korreliert - und die kann sich zwischen Bedingungen
unterscheiden (Konfund-Risiko).

Jetzt wird JEDE neu auftauchende Zelle als Bud-Kandidat behandelt
(kein Längenfilter mehr bei der Erkennung). Ob eine Zelle als potenzielle
MUTTER für das räumliche Matching infrage kommt, wird KAUSAL entschieden
(war sie zum Zeitpunkt des Bud-Auftauchens schon >= established_min_frames
Frames sichtbar?) statt über ihre finale, rückblickende Gesamttracklänge.
`bud_max_frames` ist dadurch kein Ausschlusskriterium mehr, sondern nur
noch eine Report-Markierung (`bud_was_washed_out` im Ergebnis). Dadurch
werden auch Buds erkannt, die selbst zu Müttern heranwachsen - die
Mutter-Tochter-Kette bricht nicht mehr künstlich nach einer Generation ab.

GROESSENKRITERIUM (angespuelte Zellen sind keine Knospen):
In den Kammern werden laufend Blastokonidien aus anderen Kammern angespuelt.
Sie tauchen "neu" auf, oft direkt neben einer sitzenden Zelle, und bestehen
damit die raeumliche Zuordnung wie eine Knospe - sind aber beim ersten
Auftreten etwa so gross wie die vermeintliche Mutter, waehrend eine echte
Knospe deutlich kleiner beginnt. Deshalb wird pro Kandidat das Verhaeltnis
bud_area / mother_area beim ersten Auftreten berechnet (Spalte
bud_area_fraction), und ein Kandidat oberhalb der Schwelle wird verworfen.
Die Schwelle kommt aus den Daten (bud_size.py: Antimodus der zweigipfligen
Verteilung) und wird von run_analysis.py als bud_size_threshold uebergeben.
Ist die Verteilung nicht zweigipflig - auf den echten Daten der Fall -, ist
die Schwelle unendlich und es greift KEIN Filter; ebenso ohne Schwelle
(None). Zusaetzlich traegt jedes Event Diagnose-Spalten (Flaechenabnahme der
Mutter beim Auftauchen, Wachstum des Kandidaten danach, Kontaktverhaeltnis),
mit denen bud_size.py nach einem tragfaehigen Unterscheidungsmerkmal sucht.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


@dataclass
class LineageParams:
    """Alle Schwellenwerte der Heuristik an einem Ort, mit Begründung im Docstring."""
    mother_min_frames: int = 20      # Mindestanzahl Frames, um als KANONISCHE (offizielle) Mutter zu gelten 2025_Xiao µ/1 =0.29 -> 20 Frames decken min eine reproduktion ab. Wird von identify_mothers()/compute_budding_ratio() als Downstream-Qualitaetsfilter genutzt, NICHT mehr zur Event-Erkennung selbst.
    bud_max_frames: int = 5          # NUR NOCH Report-Schwelle: kennzeichnet im Ergebnis, ob der Bud vermutlich weggespuelt wurde (bud_was_washed_out). Filtert NICHT mehr, ob ein Event erkannt wird - siehe Modul-Docstring "FIX".
    established_min_frames: int = 3  # NEU, kausal: wie viele Frames muss eine Zelle VOR einem moeglichen Bud-Auftauchen schon sichtbar gewesen sein, um als Mutter-KANDIDAT fuer das raeumliche Matching zu gelten. Bewusst klein gehalten (nicht mother_min_frames!), sonst waere das exakt derselbe Bug nur mit anderem Namen. Sollte gross genug sein, um Tracking-Rauschen im allerersten Frame einer Zelle abzufangen, aber klein genug, um schnelle Folge-Ereignisse nicht zu verpassen.
    tolerance_px: float = 30.0       # Zusätzlicher Puffer zum adaptiven Mutter-Radius (in Pixeln)
    bud_max_area_fraction: Optional[float] = None  # Groessenkriterium: ein Kandidat zaehlt nur als Knospe, wenn bud_area / mother_area beim ersten Auftreten <= dieser Wert ist. None = kein Filter, solange classify_mother_bud() keine Schwelle uebergeben bekommt (run_analysis.py leitet sie aus den Daten ab, siehe bud_size.py). Ein fester Wert hier ist der Weg fuer inspect_lineage.py, wenn keine abgeleitete Schwelle vorliegt.


@dataclass
class FluxChannelConfig:
    """
    Definiert das Kanal-Paar, aus dem der 'physiologische Flux' zum
    Budding-Zeitpunkt der Mutter berechnet wird (channel_a / channel_b).

    Analog zu sensors.RatioSensorConfig - bewusst getrennt gehalten, da der
    Flux hier NICHT pro Zeile im ganzen Datensatz berechnet wird, sondern
    NUR zum Zeitpunkt eines Budding-Events (siehe classify_mother_bud).

    channel_a / channel_b : Spaltennamen der beiden Kanäle in Combined_Results.csv
    min_denominator : Nenner-Werte unterhalb dieser Schwelle werden NICHT
                durch geteilt (Flux wird NaN statt eines extremen/instabilen
                Werts) - siehe sensors.py für denselben Mechanismus und
                Kalibrierungshinweise.
    """
    channel_a: str
    channel_b: str
    min_denominator: float = 1.0


def classify_mother_bud(
    df: pd.DataFrame,
    params: Optional[LineageParams] = None,
    flux_config: Optional[FluxChannelConfig] = None,
    bud_size_threshold: Optional[float] = None,
) -> pd.DataFrame:
    """
    Identifiziert Budding-Events basierend auf räumlicher Nähe (adaptiver,
    größenabhängiger Radius) und Track-Länge, und liefert eine Tabelle AUF
    EVENT-EBENE zurück (eine Zeile pro erkanntem Bud, der einer Mutter
    zugeordnet werden konnte).

    Parameters
    ----------
    df : Datensatz mit Spalten: exp_id, track_id, frame, centroid_x,
         centroid_y, area (sowie eccentricity falls vorhanden, sowie die in
         flux_config referenzierten Kanal-Spalten falls flux_config gesetzt ist)
    params : LineageParams, Standard wird verwendet falls None
    flux_config : FluxChannelConfig. Falls None, wird KEIN Flux berechnet
         (Spalte 'pre_budding_flux' fehlt dann in der Ausgabe) - praktisch
         für intensiometrische Sensoren oder wenn der Flux (noch) nicht
         interessiert.
    bud_size_threshold : Groessenkriterium (siehe Modul-Docstring). Falls
         None, gilt params.bud_max_area_fraction; ist auch das None (oder
         die Schwelle unendlich), wird KEIN Kandidat wegen seiner Groesse
         verworfen. Verworfene Kandidaten erscheinen nicht in der Ausgabe,
         ihre Zahl steht im Log.

    Returns
    -------
    DataFrame mit einer Zeile pro Bud-Event:
        exp_id, mother_track_id, mother_cell_uid, bud_track_id, bud_cell_uid,
        budding_frame, mother_area, mother_eccentricity (falls vorhanden),
        distance_px, adaptive_radius_px, pre_budding_flux (falls flux_config gesetzt),
        bud_area, bud_area_fraction (Flaeche des Buds beim ersten Auftreten,
            absolut und relativ zur Mutter), bud_size_threshold (die
            angewandte Schwelle, NaN = kein Groessenfilter),
        mother_area_prev, mother_area_next, mother_area_drop,
            mother_area_drop_over_bud, mother_age_frames, bud_area_plus1,
            bud_area_plus3, contact_ratio, bud_eccentricity (falls vorhanden):
            Diagnose fuer bud_size.py - verliert die Mutter beim Auftauchen
            des Kandidaten Flaeche (eine echte Knospe wird aus ihrer Maske
            herausgeloest), waechst der Kandidat danach, sitzt er an der
            Mutter an (contact_ratio ~ 1)?
        bud_final_track_length, bud_was_washed_out (finale Gesamttracklaenge
            des Buds bzw. ob sie <= bud_max_frames liegt - reine Report-Info,
            siehe Modul-Docstring "FIX"),
        mother_is_canonical_mother (ob die Mutter zusaetzlich die STRENGERE
            mother_min_frames-Schwelle erreicht, also auch in
            identify_mothers()/compute_budding_ratio() als offizielle Mutter
            gefuehrt wird - kann False sein, wenn die Mutter zum
            Budding-Zeitpunkt zwar schon >= established_min_frames sichtbar
            war, aber insgesamt kuerzer lebt als mother_min_frames)
    """
    if params is None:
        params = LineageParams()

    threshold = bud_size_threshold if bud_size_threshold is not None else params.bud_max_area_fraction
    apply_size = threshold is not None and np.isfinite(threshold)
    if apply_size:
        threshold = float(threshold)

    required_cols = {"exp_id", "cell_uid", "frame", "centroid_x", "centroid_y", "area"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"classify_mother_bud() fehlen Spalten: {missing}")

    if flux_config is not None:
        for col in (flux_config.channel_a, flux_config.channel_b):
            if col not in df.columns:
                raise ValueError(
                    f"flux_config referenziert Spalte '{col}', die nicht in den Daten vorkommt. "
                    f"Verfügbare Spalten: {sorted(df.columns)}"
                )

    has_eccentricity = "eccentricity" in df.columns
    has_cell_uid = "cell_uid" in df.columns

    results = []
    n_unassigned_buds = 0
    n_total_bud_candidates = 0
    n_rejected_by_size = 0

    for exp_id, group in df.groupby("exp_id"):
        group = group.sort_values("frame")
        n_total_frames = group["frame"].nunique()

        if n_total_frames < params.mother_min_frames:
            logger.warning(
                "exp_id '%s' hat nur %d Frames (< mother_min_frames=%d) - wird übersprungen.",
                exp_id, n_total_frames, params.mother_min_frames,
            )
            continue

        # Finale (rueckblickende) Gesamttracklaenge pro Track - wird NICHT
        # mehr zur Erkennung benutzt, sondern nur noch fuer Report-Flags
        # (bud_was_washed_out, mother_is_canonical_mother) und die
        # mother_min_frames-Kandidatenliste unten.
        track_lengths = group.groupby("cell_uid")["frame"].nunique()

        # Sortierte Frame-Liste pro Track - Grundlage der KAUSALEN
        # "etabliert zum Zeitpunkt bud_frame"-Pruefung weiter unten.
        frames_by_track = group.groupby("cell_uid")["frame"].apply(lambda s: np.sort(s.unique()))
        # Flaeche je (Zelle, Frame) fuer die Diagnose-Spalten: Mutterflaeche
        # im Frame VOR dem Auftauchen, Flaeche des Kandidaten danach.
        area_lookup = dict(zip(zip(group["cell_uid"], group["frame"]), group["area"]))

        # 1. ALLE neu auftauchenden Tracks sind Bud-KANDIDATEN - bewusst
        # KEIN Filter auf ihre eigene (finale) Tracklaenge mehr (siehe
        # Modul-Docstring "FIX"): ob ein Bud spaeter weggespuelt wird oder
        # selbst zur Mutter heranwaechst, darf nicht darueber entscheiden,
        # ob das Ereignis ueberhaupt erkannt wird.
        first_appearance = group.groupby("cell_uid")["frame"].min()
        potential_buds = first_appearance[first_appearance > group["frame"].min()].index

        if len(potential_buds) == 0:
            continue

        bud_data = group[group["cell_uid"].isin(potential_buds)]
        n_total_bud_candidates += len(potential_buds)

        # 2. Räumliche Zuordnung mit adaptivem Radius
        for bud_tid, bud_group in bud_data.groupby("cell_uid"):
            bud_frame = bud_group["frame"].iloc[0]
            bud_x = bud_group["centroid_x"].iloc[0]
            bud_y = bud_group["centroid_y"].iloc[0]

            # Kandidaten fuer die Mutterrolle: alle ANDEREN Zellen im selben
            # Frame, die zum Zeitpunkt bud_frame bereits KAUSAL etabliert
            # waren (>= established_min_frames Frames VOR bud_frame
            # gesehen) - nicht ihre globale/finale Gesamttracklaenge, die
            # zum Erkennungszeitpunkt "in der Zukunft" liegt und ausserdem
            # genau der Bias waere, den dieser Fix beheben soll.
            same_frame = group[(group["frame"] == bud_frame) & (group["cell_uid"] != bud_tid)]
            if same_frame.empty:
                n_unassigned_buds += 1
                continue

            est_counts = same_frame["cell_uid"].map(
                lambda uid: int(np.searchsorted(frames_by_track[uid], bud_frame))
            ).to_numpy()
            moms_in_frame = same_frame[est_counts >= params.established_min_frames]

            if moms_in_frame.empty:
                n_unassigned_buds += 1
                continue

            dists = np.sqrt((moms_in_frame["centroid_x"] - bud_x) ** 2 + (moms_in_frame["centroid_y"] - bud_y) ** 2)
            # Adaptiver Radius: geschätzter Mutter-Radius (aus Kreisflächen-
            # Annahme) + fester Toleranz-Puffer
            radii = np.sqrt(moms_in_frame["area"] / np.pi) + params.tolerance_px

            is_within_range = dists <= radii

            if not is_within_range.any():
                n_unassigned_buds += 1
                continue

            valid_indices = dists[is_within_range].index
            best_mom_idx = valid_indices[dists[valid_indices].argmin()]
            mom_row = moms_in_frame.loc[best_mom_idx]

            # Groessenkriterium: Flaeche des Kandidaten beim ERSTEN Auftreten
            # relativ zur Mutter. Eine Knospe beginnt klein; eine angespuelte
            # Zelle ist etwa so gross wie die Zelle, neben der sie landet.
            bud_area = float(bud_group["area"].iloc[0])
            mother_area = float(mom_row["area"])
            bud_area_fraction = bud_area / mother_area if mother_area > 0 else np.nan
            if apply_size and bud_area_fraction > threshold:
                n_rejected_by_size += 1
                continue

            bud_final_length = int(track_lengths[bud_tid])
            mother_final_length = int(track_lengths[mom_row["cell_uid"]])

            record = {
                "exp_id": exp_id,
                "mother_track_id": mom_row["track_id"],
                "bud_track_id": bud_tid,
                "budding_frame": bud_frame,
                "mother_area": mom_row["area"],
                "distance_px": dists[best_mom_idx],
                "adaptive_radius_px": radii[best_mom_idx],
                "bud_area": bud_area,
                "bud_area_fraction": bud_area_fraction,
                "bud_final_track_length": bud_final_length,
                "bud_was_washed_out": bud_final_length <= params.bud_max_frames,
                "mother_is_canonical_mother": mother_final_length >= params.mother_min_frames,
            }
            # Diagnose-Spalten (siehe Docstring): Flaechenbilanz der Mutter um
            # das Auftauchen herum, Wachstum des Kandidaten, Kontakt.
            mom_uid = mom_row["cell_uid"]
            mom_frames = frames_by_track[mom_uid]
            k_prev = int(np.searchsorted(mom_frames, bud_frame)) - 1
            mother_area_prev = (float(area_lookup.get((mom_uid, mom_frames[k_prev]), np.nan))
                                if k_prev >= 0 else np.nan)
            mother_area_drop = mother_area_prev - mother_area if np.isfinite(mother_area_prev) else np.nan
            contact = np.sqrt(mother_area / np.pi) + np.sqrt(bud_area / np.pi) if mother_area > 0 and bud_area > 0 else np.nan
            record.update({
                "mother_area_prev": mother_area_prev,
                "mother_area_next": float(area_lookup.get((mom_uid, bud_frame + 1), np.nan)),
                "mother_area_drop": mother_area_drop,
                "mother_area_drop_over_bud": (mother_area_drop / bud_area
                                              if np.isfinite(mother_area_drop) and bud_area > 0 else np.nan),
                "mother_age_frames": k_prev + 1,
                "bud_area_plus1": float(area_lookup.get((bud_tid, bud_frame + 1), np.nan)),
                "bud_area_plus3": float(area_lookup.get((bud_tid, bud_frame + 3), np.nan)),
                "contact_ratio": float(dists[best_mom_idx]) / contact if np.isfinite(contact) and contact > 0 else np.nan,
            })
            if has_eccentricity:
                record["mother_eccentricity"] = mom_row["eccentricity"]
                record["bud_eccentricity"] = bud_group["eccentricity"].iloc[0]
            if has_cell_uid:
                record["mother_cell_uid"] = mom_row["cell_uid"]
                record["bud_cell_uid"] = bud_group["cell_uid"].iloc[0]

            if flux_config is not None:
                denom = mom_row[flux_config.channel_b]
                if abs(denom) < flux_config.min_denominator:
                    record["pre_budding_flux"] = np.nan
                else:
                    record["pre_budding_flux"] = mom_row[flux_config.channel_a] / denom

            results.append(record)

    out = pd.DataFrame(results)
    if not out.empty:
        out["bud_size_threshold"] = threshold if apply_size else np.nan

    if apply_size and n_total_bud_candidates > 0:
        logger.info(
            "Groessenkriterium: %d von %d Bud-Kandidaten verworfen (Flaeche beim ersten "
            "Auftreten > %.2f x Mutterflaeche) - mutmasslich angespuelte Zellen, keine Knospen.",
            n_rejected_by_size, n_total_bud_candidates, threshold,
        )

    if n_total_bud_candidates > 0 and n_unassigned_buds > 0:
        logger.warning(
            "%d von %d Bud-Kandidaten konnten KEINER Mutter zugeordnet werden "
            "(keine etablierte Zelle im adaptiven Radius zum Zeitpunkt des "
            "Erstauftretens). Diese fließen NICHT in die Budding Ratio ein - "
            "ggf. tolerance_px oder established_min_frames anpassen.",
            n_unassigned_buds, n_total_bud_candidates,
        )
    if not out.empty:
        n_not_washed = int((~out["bud_was_washed_out"]).sum())
        n_non_canonical_mothers = int((~out["mother_is_canonical_mother"]).sum())
        logger.info(
            "Mutter/Bud-Klassifikation: %d Budding-Events erkannt (aus %d Bud-Kandidaten), "
            "davon %d von Buds, die NICHT weggespült wurden (bud_was_washed_out=False - "
            "diese Events wären vor dem Fix verloren gegangen), %d mit einer Mutter, die "
            "nicht die mother_min_frames-Schwelle erreicht (mother_is_canonical_mother=False, "
            "fließt daher NICHT in identify_mothers()/compute_budding_ratio() ein, wohl aber "
            "in growth_rate.compute_specific_growth_rate()).",
            len(out), n_total_bud_candidates, n_not_washed, n_non_canonical_mothers,
        )
    else:
        logger.info(
            "Mutter/Bud-Klassifikation: %d Budding-Events erkannt (aus %d Bud-Kandidaten).",
            len(out), n_total_bud_candidates,
        )

    return out


def identify_mothers(df: pd.DataFrame, params: Optional[LineageParams] = None) -> pd.DataFrame:
    """
    Identifiziert ALLE Mutterzellen (Tracks mit >= mother_min_frames Frames)
    je Kammer (exp_id) - UNABHÄNGIG davon, ob sie ein Budding-Event hatten.

    Verwendet exakt dieselben Kriterien wie Schritt 1 in classify_mother_bud()
    (inkl. der Kammer-Mindestlänge n_total_frames >= mother_min_frames), damit
    die Mütterliste konsistent zu den dort erkannten Events ist. Wird von
    compute_budding_ratio() gebraucht, damit RUHENDE (nicht-knospende) Mütter
    mit n_buds=0 in der Budding-Ratio-Verteilung erscheinen, statt
    stillschweigend zu fehlen (siehe compute_budding_ratio()-Docstring).

    COVERAGE: n_frames (>= mother_min_frames) ist eine ABSOLUTE Framezahl und
    unterscheidet nicht zwischen einer durchgehend beobachteten Zelle und
    einer mit vielen Tracking-Lücken. Eine Mutter mit 20 beobachteten Frames
    verteilt über 200 Kammer-Frames (coverage=0.1) ist trackingtechnisch
    etwas anderes als eine mit 20 von 22 möglichen Frames (coverage≈0.9) -
    frame_span/coverage machen das sichtbar, ohne n_frames selbst zu
    verändern. Nutzt select_stable_mothers() unten, um darauf zu filtern.

    Returns
    -------
    DataFrame mit einer Zeile pro Mutter-cell_uid:
        exp_id, cell_uid, n_frames, frame_min, frame_max, frame_span,
        coverage (= n_frames / frame_span, in (0, 1])
    """
    if params is None:
        params = LineageParams()

    required_cols = {"exp_id", "cell_uid", "frame"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"identify_mothers() fehlen Spalten: {missing}")

    results = []
    for exp_id, group in df.groupby("exp_id"):
        n_total_frames = group["frame"].nunique()
        if n_total_frames < params.mother_min_frames:
            continue
        per_track = group.groupby("cell_uid")["frame"].agg(n_frames="nunique", frame_min="min", frame_max="max")
        mothers = per_track[per_track["n_frames"] >= params.mother_min_frames]
        for cell_uid, row in mothers.iterrows():
            frame_span = int(row["frame_max"] - row["frame_min"] + 1)
            results.append({
                "exp_id": exp_id,
                "cell_uid": cell_uid,
                "n_frames": int(row["n_frames"]),
                "frame_min": int(row["frame_min"]),
                "frame_max": int(row["frame_max"]),
                "frame_span": frame_span,
                "coverage": row["n_frames"] / frame_span if frame_span > 0 else np.nan,
            })

    out = pd.DataFrame(results)
    if not out.empty:
        logger.info(
            "identify_mothers(): %d Mutterzellen identifiziert (>= %d Frames). "
            "Coverage-Median=%.2f, Coverage<0.5 bei %d Müttern (lückenhaft getrackt).",
            len(out), params.mother_min_frames, out["coverage"].median(),
            int((out["coverage"] < 0.5).sum()),
        )
    else:
        logger.info("identify_mothers(): %d Mutterzellen identifiziert (>= %d Frames).", 0, params.mother_min_frames)
    return out


def select_stable_mothers(
    mothers: pd.DataFrame,
    min_coverage: float = 0.7,
    min_n_frames: Optional[int] = None,
) -> pd.DataFrame:
    """
    Filtert identify_mothers()-Output auf "stabil getrackte" Mütter: sowohl
    lang genug (n_frames, per Default die mother_min_frames-Schwelle, die
    bereits in identify_mothers() steckt) ALS AUCH lückenarm beobachtet
    (coverage >= min_coverage). Reine Downstream-Auswahl für Darstellung/
    Detailanalyse (z.B. mother_trajectories.plot_mother_trajectories()) -
    verändert nichts an identify_mothers()/compute_budding_ratio() selbst.

    Parameters
    ----------
    mothers : Ergebnis von identify_mothers() (braucht Spalte 'coverage')
    min_coverage : Mindestanteil beobachteter Frames im eigenen Frame-Fenster
                   (Default 0.7 - großzügig, ggf. an eure Trackingqualität
                   anpassen; mit inspect_classification() an echten
                   QC-Overlays kalibrieren)
    min_n_frames : optionaler zusätzlicher, strengerer Mindest-Framezahl-Filter
                   (z.B. 40, wenn "stabil" für euch mehr heißt als die
                   generelle mother_min_frames-Schwelle aus LineageParams)

    Returns
    -------
    Teilmenge von 'mothers', sortiert absteigend nach coverage.
    """
    if mothers.empty:
        return mothers

    if "coverage" not in mothers.columns:
        raise ValueError(
            "select_stable_mothers() erwartet eine 'coverage'-Spalte - "
            "stelle sicher, dass 'mothers' von identify_mothers() stammt."
        )

    out = mothers[mothers["coverage"] >= min_coverage]
    if min_n_frames is not None:
        out = out[out["n_frames"] >= min_n_frames]

    logger.info(
        "select_stable_mothers(): %d von %d Müttern erfüllen coverage >= %.2f%s.",
        len(out), len(mothers), min_coverage,
        f" und n_frames >= {min_n_frames}" if min_n_frames is not None else "",
    )
    return out.sort_values("coverage", ascending=False).reset_index(drop=True)




def compute_budding_ratio(
    lineage_events: pd.DataFrame,
    cells: pd.DataFrame,
    mothers: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Berechnet die Budding Ratio (Buds pro Mutterzelle) aus den von
    classify_mother_bud() erkannten Events - für ALLE Mutterzellen aus
    `mothers` (siehe identify_mothers()), NICHT nur für die, die tatsächlich
    ein Budding-Event hatten. Mütter ohne zugeordneten Bud bekommen
    n_buds=0 / budding_ratio=0.

    Parameters
    ----------
    lineage_events : Ergebnis von classify_mother_bud() (eine Zeile pro
            Bud-Event). Darf leer sein - dann bekommen alle Mütter n_buds=0.
    cells : der volle Zelldatensatz (für Hierarchie-Metadaten)
    mothers : Ergebnis von identify_mothers() - ALLE Mutterzellen dieser
            Kammern, unabhängig davon ob sie geknospt haben.
            WICHTIG: eine frühere Version dieser Funktion leitete die
            Mütterliste ausschließlich aus lineage_events ab - d.h. NUR
            knospende Mütter wurden gezählt, ruhende Mütter (budding_ratio=0)
            fehlten komplett. Das hat die Budding-Ratio-Verteilung (Panel A)
            UND die per_experiment-Ratio (falscher, zu kleiner n_mothers-
            Nenner) systematisch nach oben verzerrt.

    Returns
    -------
    per_mother : eine Zeile pro Mutter-cell_uid mit n_buds (>= 0) und
                 budding_ratio (== n_buds, da 1 Mutter pro Zeile)
    per_experiment : aggregierte Budding Ratio pro exp_id
    """
    if "cell_uid" not in cells.columns:
        raise ValueError("compute_budding_ratio() braucht eine 'cell_uid' Spalte in 'cells'.")

    if mothers.empty:
        logger.warning(
            "Keine Mutterzellen (>= mother_min_frames) vorhanden - compute_budding_ratio() "
            "gibt leere Tabellen zurück."
        )
        return pd.DataFrame(), pd.DataFrame()

    if lineage_events.empty or "mother_cell_uid" not in lineage_events.columns:
        buds_per_mother = pd.DataFrame(columns=["cell_uid", "n_buds"])
    else:
        buds_per_mother = (
            lineage_events.groupby("mother_cell_uid")
            .size()
            .rename("n_buds")
            .reset_index()
            .rename(columns={"mother_cell_uid": "cell_uid"})
        )

    mother_uids = set(mothers["cell_uid"])
    mothers_meta = cells.drop_duplicates("cell_uid")
    mothers_meta = mothers_meta[mothers_meta["cell_uid"].isin(mother_uids)]

    n_missing_meta = len(mother_uids) - len(mothers_meta)
    if n_missing_meta > 0:
        logger.warning(
            "%d von %d Mutterzellen aus 'mothers' kommen nicht in 'cells' vor "
            "(Metadaten fehlen) - werden übersprungen.",
            n_missing_meta, len(mother_uids),
        )

    per_mother = mothers_meta.merge(buds_per_mother, on="cell_uid", how="left")
    per_mother["n_buds"] = per_mother["n_buds"].fillna(0).astype(int)
    per_mother["budding_ratio"] = per_mother["n_buds"]

    per_experiment = (
        per_mother.groupby("exp_id")
        .agg(n_mothers=("cell_uid", "nunique"), n_buds_total=("n_buds", "sum"))
        .reset_index()
    )
    per_experiment["budding_ratio"] = per_experiment["n_buds_total"] / per_experiment["n_mothers"]

    meta_cols = [c for c in ["biosensor", "osc_type", "osc_freq", "condition", "replicate", "chamber"]
                 if c in cells.columns]
    meta = cells[["exp_id"] + meta_cols].drop_duplicates("exp_id")
    per_experiment = per_experiment.merge(meta, on="exp_id", how="left")

    logger.info(
        "Budding Ratio berechnet: %d Mutterzellen (%d mit >=1 Bud, %d ruhend), mittlere Ratio=%.3f.",
        len(per_mother), int((per_mother["n_buds"] > 0).sum()), int((per_mother["n_buds"] == 0).sum()),
        per_mother["budding_ratio"].mean(),
    )

    return per_mother, per_experiment


def inspect_classification(lineage_events: pd.DataFrame, cells: pd.DataFrame, exp_id: str) -> pd.DataFrame:
    """
    Hilfsfunktion zur Kalibrierung: zeigt für EIN Experiment (Kammer) alle
    erkannten Budding-Events mit Distanz/Radius-Details, damit man die
    Schwellenwerte (mother_min_frames, bud_max_frames, tolerance_px) gegen
    echte QC-Overlay-Bilder validieren kann, bevor man den Ergebnissen
    vertraut.
    """
    sub = lineage_events[lineage_events["exp_id"] == exp_id].copy()
    n_frames = cells.loc[cells["exp_id"] == exp_id, "frame"].nunique()
    print(f"exp_id = {exp_id}  (insgesamt {n_frames} Frames in dieser Kammer, {len(sub)} Budding-Events erkannt)")
    cols = [c for c in [
        "mother_track_id", "bud_track_id", "budding_frame",
        "distance_px", "adaptive_radius_px", "mother_area", "bud_area", "bud_area_fraction",
        "pre_budding_flux",
        "bud_final_track_length", "bud_was_washed_out", "mother_is_canonical_mother",
    ] if c in sub.columns]
    return sub[cols].sort_values("budding_frame")
