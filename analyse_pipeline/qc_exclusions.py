"""
qc_exclusions.py
=================
Manuelles QC-Editieren: einzelne Tracks (oder Teile davon) ausschließen,
ohne die Rohdaten der Bildverarbeitungs-Pipeline zu verändern. Zusätzlich:
TRACK-MERGES für den Fall, dass die Bildverarbeitung eine durchgängige Zelle
fälschlich in zwei Tracks aufgespalten hat (z.B. weil sich die Zelle kurz so
stark bewegt hat, dass der Tracker sie als neue Zelle erkennt - das sieht
sonst aus wie ein neuer Bud, ist aber dieselbe Zelle).

IDEE
----
Es gibt EINE editierbare Tabelle (qc_exclusions.csv), die du von Hand in
Excel/Editor oder per Hilfsfunktionen add_qc_exclusion() / add_track_merge()
pflegst. Jede Zeile sagt entweder:
  (a) "dieser Track in diesem Experiment ist KEINE valide Zelle (oder aus
      anderem Grund raus), aus folgendem Grund." (reason gesetzt, merge_
      into_track_id leer), oder
  (b) "dieser Track ist eigentlich Teil von Track X - bitte zusammenführen,
      nicht ausschließen." (merge_into_track_id gesetzt)

Beim Laden der Analyse-Daten wird diese Tabelle gegen cell_uid gejoint.
TRACK-MERGES MÜSSEN VOR EXCLUSIONS angewendet werden (siehe
apply_track_merges() vs. apply_qc_exclusions()) - die empfohlene Reihenfolge
ist im Hauptskript (run_analysis.py) bereits so vorgesehen.

Spalten der qc_exclusions.csv
------------------------------
cell_uid            - eindeutige ID des BETROFFENEN Tracks: biosensor__osc_type__osc_freq__condition__replicate__chamber__trackN
reason              - Freitext, z.B. "Debris, keine echte Zelle", "Tracking-Fehler ab Frame 40"
excluded_by          - wer hat den Eintrag vorgenommen (Kürzel/Name), für Nachvollziehbarkeit
excluded_on          - Datum (ISO, YYYY-MM-DD)
frame_from           - optional: nur ab diesem Frame ausschließen (leer = alle Frames). NICHT für Merges relevant.
frame_to             - optional: nur bis zu diesem Frame ausschließen (leer = bis Ende). NICHT für Merges relevant.
merge_into_track_id  - optional: NUR für Track-Merges gesetzt. Die track_id
                       (innerhalb derselben exp_id!), in die cell_uid umbenannt/
                       verschmolzen werden soll. Wenn gesetzt, wird reason
                       als Begründung des Merges interpretiert, NICHT als Exclusion-Grund.

Die frame_from/frame_to Spalten erlauben PARTIELLE Exclusions: z.B. wenn ein
Track ab Frame 40 mit einem Nachbarn verschmilzt, aber bis dahin valide war,
muss nicht der ganze Track raus. (Für den Fall, dass ein Track-BRUCH - nicht
Verschmelzung mit einer ANDEREN Zelle - korrigiert werden soll, ist
merge_into_track_id der richtige Mechanismus, nicht frame_from/frame_to.)
"""

from __future__ import annotations

import logging
from datetime import date
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

QC_EXCLUSION_COLUMNS = [
    "cell_uid", "reason", "excluded_by", "excluded_on", "frame_from", "frame_to", "merge_into_track_id",
]


def init_qc_exclusions(path: str | Path) -> None:
    """Erstellt eine leere QC-Exclusion-Tabelle, falls noch keine existiert."""
    path = Path(path)
    if path.exists():
        logger.info("QC-Exclusion-Datei existiert bereits: %s (wird nicht überschrieben)", path)
        return
    empty = pd.DataFrame(columns=QC_EXCLUSION_COLUMNS)
    path.parent.mkdir(parents=True, exist_ok=True)
    empty.to_csv(path, index=False)
    logger.info("Leere QC-Exclusion-Datei angelegt: %s", path)


def read_qc_exclusions(path: str | Path) -> pd.DataFrame:
    """Liest die QC-Exclusion-Tabelle ein (leeres DataFrame falls nicht vorhanden).

    Rückwärtskompatibel: falls eine ältere qc_exclusions.csv ohne die Spalte
    'merge_into_track_id' existiert (von vor Einführung der Track-Merge-
    Funktion), wird die Spalte automatisch mit leeren Werten ergänzt, statt
    einen Fehler zu werfen.
    """
    path = Path(path)
    if not path.exists():
        logger.info("Keine QC-Exclusion-Datei unter '%s' gefunden - es wird nichts ausgeschlossen.", path)
        return pd.DataFrame(columns=QC_EXCLUSION_COLUMNS)

    df = pd.read_csv(path)

    if "merge_into_track_id" not in df.columns:
        logger.info(
            "Spalte 'merge_into_track_id' fehlt in '%s' (ältere Datei ohne Track-Merge-Funktion) "
            "- wird automatisch ergänzt (leer).",
            path,
        )
        df["merge_into_track_id"] = np.nan

    missing = [c for c in QC_EXCLUSION_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(
            f"qc_exclusions.csv fehlen Spalten: {missing}. "
            f"Erwartet werden: {QC_EXCLUSION_COLUMNS}"
        )
    return df


def add_qc_exclusion(
    path: str | Path,
    biosensor: str,
    osc_type: str,
    osc_freq: str,
    condition: str,
    replicate: str,
    chamber: str,
    track_id: int,
    reason: str,
    excluded_by: str,
    frame_from: Optional[int] = None,
    frame_to: Optional[int] = None,
) -> pd.DataFrame:
    """
    Fügt der QC-Exclusion-Tabelle EINEN neuen Eintrag hinzu und schreibt sie
    zurück auf Platte. Praktisch direkt nachdem man sich ein QC-Overlay-Tif
    angeschaut hat, z.B.:

        add_qc_exclusion(
            path="qc_exclusions.csv",
            biosensor="iGlucoSnFR", osc_type="Glucose", osc_freq="5min",
            condition="Osc5min", replicate="Rep1", chamber="ChamA",
            track_id=14,
            reason="Debris, keine Zelle (siehe QC-Overlay Frame 0-10)",
            excluded_by="MK",
        )

    frame_from/frame_to weglassen = ganzer Track wird ausgeschlossen.
    """
    path = Path(path)
    if not path.exists():
        init_qc_exclusions(path)

    exp_id = "__".join(str(x) for x in [biosensor, osc_type, osc_freq, condition, replicate, chamber])
    cell_uid = f"{exp_id}__track{track_id}"

    existing = read_qc_exclusions(path)

    new_row = pd.DataFrame([{
        "cell_uid": cell_uid,
        "reason": reason,
        "excluded_by": excluded_by,
        "excluded_on": date.today().isoformat(),
        "frame_from": frame_from if frame_from is not None else np.nan,
        "frame_to": frame_to if frame_to is not None else np.nan,
        "merge_into_track_id": np.nan,
    }])

    out = pd.concat([existing, new_row], ignore_index=True)
    out.to_csv(path, index=False)
    logger.info("Exclusion hinzugefügt: %s (Grund: %s)", cell_uid, reason)
    return out


def add_track_merge(
    path: str | Path,
    biosensor: str,
    osc_type: str,
    osc_freq: str,
    condition: str,
    replicate: str,
    chamber: str,
    track_id: int,
    merge_into_track_id: int,
    reason: str,
    merged_by: str,
) -> pd.DataFrame:
    """
    Markiert track_id als Teil von merge_into_track_id (z.B. weil die Zelle
    sich kurz so stark bewegt hat, dass der Tracker einen neuen Track
    begonnen hat - das sähe sonst wie ein neuer Bud aus, ist aber dieselbe
    Zelle). NICHT-DESTRUKTIV: die Rohdaten werden nicht verändert, die
    Korrektur lebt ausschließlich in qc_exclusions.csv. Anwendung beim
    Laden über apply_track_merges() - VOR apply_qc_exclusions() und vor
    jeglicher Lineage-Klassifikation aufrufen.

    WICHTIG: track_id und merge_into_track_id müssen INNERHALB DERSELBEN
    Kammer (exp_id) liegen - ein Merge über Kammern hinweg ist biologisch
    nicht sinnvoll und wird von apply_track_merges() abgelehnt.

    Beispiel: Track 23 verschwindet bei Frame 15, Track 31 erscheint bei
    Frame 16 an fast derselben Position und wird fälschlich als neuer Bud
    erkannt - tatsächlich ist es dieselbe Zelle, die der Tracker verloren hat:

        add_track_merge(
            path="qc_exclusions.csv",
            biosensor="iGlucoSnFR", osc_type="Glucose", osc_freq="5min",
            condition="Osc5min", replicate="Rep1", chamber="ChamA",
            track_id=31, merge_into_track_id=23,
            reason="Tracking-Bruch durch starke Zellbewegung, Frame 15->16, gleiche Zelle",
            merged_by="MK",
        )
    """
    path = Path(path)
    if not path.exists():
        init_qc_exclusions(path)

    if track_id == merge_into_track_id:
        raise ValueError("track_id und merge_into_track_id dürfen nicht identisch sein.")

    exp_id = "__".join(str(x) for x in [biosensor, osc_type, osc_freq, condition, replicate, chamber])
    cell_uid = f"{exp_id}__track{track_id}"

    existing = read_qc_exclusions(path)

    new_row = pd.DataFrame([{
        "cell_uid": cell_uid,
        "reason": reason,
        "excluded_by": merged_by,
        "excluded_on": date.today().isoformat(),
        "frame_from": np.nan,
        "frame_to": np.nan,
        "merge_into_track_id": merge_into_track_id,
    }])

    out = pd.concat([existing, new_row], ignore_index=True)
    out.to_csv(path, index=False)
    logger.info(
        "Track-Merge hinzugefügt: %s -> track%d (Grund: %s)",
        cell_uid, merge_into_track_id, reason,
    )
    return out


def apply_track_merges(df: pd.DataFrame, exclusions: pd.DataFrame) -> pd.DataFrame:
    """
    Wendet alle Track-Merges (Zeilen mit gesetztem merge_into_track_id) an:
    der Quell-Track wird umbenannt (track_id und cell_uid), sodass er danach
    als Teil des Ziel-Tracks erscheint - so, als hätte die Bildverarbeitung
    von Anfang an EINEN durchgängigen Track erkannt.

    MUSS vor apply_qc_exclusions() UND vor jeglicher Lineage-Klassifikation
    (lineage.classify_mother_bud) aufgerufen werden, sonst sieht die
    Heuristik weiterhin zwei getrennte Tracks und der Tracking-Bruch wird
    fälschlich als neuer Bud interpretiert.

    SICHERHEITSPRÜFUNG: wenn Quell- und Ziel-Track sich in mindestens einem
    Frame ÜBERLAPPEN (beide gleichzeitig vorhanden), ist der Merge
    biologisch unsinnig (eine Zelle kann nicht zwei Positionen gleichzeitig
    haben) - so ein Merge wird ABGELEHNT (Warnung, Zeile wird ignoriert),
    nicht stillschweigend durchgeführt. Das ist ein Hinweis, dass die
    Track-IDs in der Korrektur vertauscht oder falsch sind.

    Parameters
    ----------
    df : voller Zelldatensatz mit exp_id, track_id, frame, cell_uid
    exclusions : Ergebnis von read_qc_exclusions() (Zeilen ohne gesetztes
                 merge_into_track_id werden ignoriert - die behandelt
                 apply_qc_exclusions())

    Returns
    -------
    df mit umbenannten track_id/cell_uid für alle erfolgreich gemergten Tracks.
    """
    merges = exclusions[exclusions["merge_into_track_id"].notna()].copy()
    if merges.empty:
        return df

    if not {"exp_id", "track_id", "frame", "cell_uid"}.issubset(df.columns):
        raise ValueError("apply_track_merges() braucht exp_id, track_id, frame, cell_uid in 'df'.")

    df = df.copy()
    merges["merge_into_track_id"] = merges["merge_into_track_id"].astype(int)
    # exp_id und track_id des QUELL-Tracks aus cell_uid zurückgewinnen
    merges["exp_id"] = merges["cell_uid"].str.replace(r"__track\d+$", "", regex=True)
    merges["src_track_id"] = merges["cell_uid"].str.extract(r"__track(\d+)$").astype(int)

    n_applied, n_rejected = 0, 0

    for _, row in merges.iterrows():
        exp_id = row["exp_id"]
        src_tid = row["src_track_id"]
        dst_tid = row["merge_into_track_id"]

        src_mask = (df["exp_id"] == exp_id) & (df["track_id"] == src_tid)
        dst_mask = (df["exp_id"] == exp_id) & (df["track_id"] == dst_tid)

        if not src_mask.any():
            logger.warning(
                "Track-Merge übersprungen: cell_uid '%s' (Quelle) kommt im Datensatz nicht vor.",
                row["cell_uid"],
            )
            n_rejected += 1
            continue
        if not dst_mask.any():
            logger.warning(
                "Track-Merge übersprungen: Ziel-Track %d in exp_id '%s' kommt im Datensatz nicht vor.",
                dst_tid, exp_id,
            )
            n_rejected += 1
            continue

        src_frames = set(df.loc[src_mask, "frame"])
        dst_frames = set(df.loc[dst_mask, "frame"])
        overlap = src_frames & dst_frames
        if overlap:
            logger.warning(
                "Track-Merge ABGELEHNT: '%s' -> track%d in exp_id '%s' überlappen sich in %d Frame(s) "
                "(%s) - eine Zelle kann nicht zwei Positionen gleichzeitig haben. "
                "Bitte Track-IDs in der Korrektur prüfen.",
                row["cell_uid"], dst_tid, exp_id, len(overlap), sorted(overlap)[:5],
            )
            n_rejected += 1
            continue

        dst_cell_uid = f"{exp_id}__track{dst_tid}"
        df.loc[src_mask, "track_id"] = dst_tid
        df.loc[src_mask, "cell_uid"] = dst_cell_uid
        n_applied += 1
        logger.info("Track-Merge angewendet: '%s' -> '%s'", row["cell_uid"], dst_cell_uid)

    logger.info("Track-Merges: %d angewendet, %d abgelehnt/übersprungen.", n_applied, n_rejected)
    return df


def apply_qc_exclusions(
    df: pd.DataFrame,
    exclusions: pd.DataFrame,
    mode: str = "remove",
) -> pd.DataFrame:
    """
    Wendet die QC-Exclusion-Tabelle auf einen Datensatz an.

    mode="remove" (Standard): betroffene Zeilen werden komplett entfernt.
        -> für finale Analyse-Plots/Statistik.
    mode="flag": betroffene Zeilen bleiben erhalten, bekommen aber eine neue
        Spalte `qc_excluded` (bool) + `qc_reason`.
        -> nützlich zur Kontrolle ("was würde eigentlich rausfliegen?") oder
           für Vorher/Nachher-Vergleiche.

    Partielle Exclusions (frame_from/frame_to gesetzt) werden in beiden
    Modi korrekt nur für den angegebenen Frame-Bereich angewendet; ein
    leerer Wert bedeutet "von Anfang" bzw. "bis Ende".

    Parameters
    ----------
    df : Datensatz mit Spalten 'cell_uid' und 'frame' (siehe data_loading.load_all_results)
    exclusions : Ergebnis von read_qc_exclusions()
    """
    if mode not in ("remove", "flag"):
        raise ValueError("mode muss 'remove' oder 'flag' sein")

    if "cell_uid" not in df.columns:
        raise ValueError("apply_qc_exclusions() braucht eine 'cell_uid' Spalte (siehe load_all_results()).")

    # Merge-Zeilen (merge_into_track_id gesetzt) gehören zu apply_track_merges(),
    # nicht hierher - sonst würde ein gemergter Track hier fälschlich auch noch
    # komplett ausgeschlossen.
    exclusions = exclusions[exclusions["merge_into_track_id"].isna()] if "merge_into_track_id" in exclusions.columns else exclusions

    if exclusions.empty:
        if mode == "flag":
            df = df.copy()
            df["qc_excluded"] = False
            df["qc_reason"] = pd.NA
        return df

    # Tippfehler/falsche cell_uid in der Exclusion-Tabelle nicht stillschweigend ignorieren
    known_uids = set(df["cell_uid"])
    unmatched = exclusions[~exclusions["cell_uid"].isin(known_uids)]
    if not unmatched.empty:
        logger.warning(
            "%d QC-Exclusion(en) passen zu KEINER Zeile im aktuellen Datensatz "
            "(falsche cell_uid? falsches Experiment ausgewählt?):\n%s",
            len(unmatched), "\n".join(f"  - {u}" for u in unmatched["cell_uid"]),
        )

    exc = exclusions.copy()
    exc["frame_from"] = exc["frame_from"].fillna(-np.inf)
    exc["frame_to"] = exc["frame_to"].fillna(np.inf)

    df = df.copy()
    df["_row_id"] = np.arange(len(df))

    merged = df[["_row_id", "cell_uid", "frame"]].merge(exc, on="cell_uid", how="inner")
    hit_mask = (merged["frame"] >= merged["frame_from"]) & (merged["frame"] <= merged["frame_to"])
    hits = merged[hit_mask]

    hit_summary = (
        hits.groupby("_row_id")["reason"]
        .apply(lambda s: "; ".join(sorted(set(s))))
        .rename("qc_reason")
        .reset_index()
    )
    hit_summary["qc_excluded"] = True

    df = df.merge(hit_summary, on="_row_id", how="left")
    df["qc_excluded"] = df["qc_excluded"].fillna(False).astype(bool)
    df = df.drop(columns=["_row_id"])

    n_excluded_rows = int(df["qc_excluded"].sum())
    n_excluded_tracks = df.loc[df["qc_excluded"], "cell_uid"].nunique()
    logger.info(
        "QC-Exclusions angewendet (mode=%s): %d Zeilen / %d Tracks betroffen.",
        mode, n_excluded_rows, n_excluded_tracks,
    )

    if mode == "remove":
        df = df.loc[~df["qc_excluded"]].drop(columns=["qc_excluded", "qc_reason"]).reset_index(drop=True)

    return df


def summarise_qc_exclusions(exclusions: pd.DataFrame) -> pd.DataFrame:
    """Zeigt, wie viele Tracks pro Experiment aktuell ausgeschlossen sind
    (Track-Merges werden hier NICHT mitgezählt - siehe summarise_track_merges).
    Guter Sanity-Check nach jeder QC-Runde."""
    only_exclusions = (
        exclusions[exclusions["merge_into_track_id"].isna()]
        if "merge_into_track_id" in exclusions.columns else exclusions
    )
    if only_exclusions.empty:
        logger.info("Keine QC-Exclusions vorhanden.")
        return pd.DataFrame(columns=["exp_id", "n_excluded_tracks"])

    exp_id = only_exclusions["cell_uid"].str.replace(r"__track\d+$", "", regex=True)
    summary = (
        exp_id.value_counts()
        .rename_axis("exp_id")
        .reset_index(name="n_excluded_tracks")
        .sort_values("n_excluded_tracks", ascending=False)
        .reset_index(drop=True)
    )
    return summary


def summarise_track_merges(exclusions: pd.DataFrame) -> pd.DataFrame:
    """Zeigt alle aktuell konfigurierten Track-Merges (Quelle -> Ziel) pro
    Experiment. Guter Sanity-Check, um versehentliche Doppel-Merges oder
    Tippfehler in Track-IDs frühzeitig zu erkennen."""
    if "merge_into_track_id" not in exclusions.columns:
        return pd.DataFrame(columns=["exp_id", "src_track_id", "merge_into_track_id", "reason"])

    merges = exclusions[exclusions["merge_into_track_id"].notna()].copy()
    if merges.empty:
        logger.info("Keine Track-Merges vorhanden.")
        return pd.DataFrame(columns=["exp_id", "src_track_id", "merge_into_track_id", "reason"])

    merges["exp_id"] = merges["cell_uid"].str.replace(r"__track\d+$", "", regex=True)
    merges["src_track_id"] = merges["cell_uid"].str.extract(r"__track(\d+)$").astype(int)
    merges["merge_into_track_id"] = merges["merge_into_track_id"].astype(int)
    return merges[["exp_id", "src_track_id", "merge_into_track_id", "reason"]].reset_index(drop=True)
