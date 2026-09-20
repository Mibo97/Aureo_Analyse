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
REIHENFOLGE (run_analysis.py): Exclusions, dann Track-Merges, dann noch
einmal Exclusions. Vorher, weil ein Ausschluss auf dem QUELL-Track nur
greift, solange der Track noch seinen alten Namen hat; nachher, weil ein
Ausschluss auf dem ZIEL-Track auch die hineingemergten Frames treffen soll.
Ausschluesse sind idempotent, der zweite Durchlauf aendert sonst nichts.
Merges muessen in jedem Fall VOR der Lineage-Klassifikation liegen.

Spalten der qc_exclusions.csv
------------------------------
cell_uid            - eindeutige ID des BETROFFENEN Tracks: biosensor__osc_type__osc_freq__condition__replicate__chamber__trackN
reason              - Freitext, z.B. "Debris, keine echte Zelle", "Tracking-Fehler ab Frame 40"
excluded_by          - wer hat den Eintrag vorgenommen (Kürzel/Name), für Nachvollziehbarkeit
excluded_on          - Datum (ISO, YYYY-MM-DD)
frame_from           - optional: nur ab diesem Frame ausschließen (leer = alle Frames). Bei Merges: nur ab diesem Frame der Quelle ins Ziel uebernehmen.
frame_to             - optional: nur bis zu diesem Frame ausschließen (leer = bis Ende). Bei Merges: nur bis zu diesem Frame der Quelle ins Ziel uebernehmen.
merge_into_track_id  - optional: NUR für Track-Merges gesetzt. Die track_id
                       (innerhalb derselben exp_id!), in die cell_uid umbenannt/
                       verschmolzen werden soll. Wenn gesetzt, wird reason
                       als Begründung des Merges interpretiert, NICHT als Exclusion-Grund.

Die frame_from/frame_to Spalten erlauben PARTIELLE Exclusions: z.B. wenn ein
Track ab Frame 40 mit einem Nachbarn verschmilzt, aber bis dahin valide war,
muss nicht der ganze Track raus. (Für den Fall, dass ein Track-BRUCH - nicht
Verschmelzung mit einer ANDEREN Zelle - korrigiert werden soll, ist
merge_into_track_id der richtige Mechanismus, nicht frame_from/frame_to.)

Mit frame_from/frame_to UND merge_into_track_id wird ein Track AUFGETEILT:
hat der Tracker nacheinander zwei verschiedene Zellen unter einer track_id
gefuehrt (Sprung), gehoeren die Frames vor dem Sprung zur einen und die
danach zur anderen Zelle. Zwei Merge-Zeilen fuer dieselbe cell_uid mit
disjunkten Frame-Bereichen und verschiedenen Zielen fuehren beide Teile
ihrem richtigen Track zu (find_qc_conflicts() nennt das 'split_merge' und
markiert es als ausgefuehrt). Zwei Merge-Zeilen OHNE Bereiche sind dagegen
ein Widerspruch ('double_merge'): die erste gewinnt, die zweite wird
verworfen.
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
    frame_from: Optional[int] = None,
    frame_to: Optional[int] = None,
) -> pd.DataFrame:
    """
    Markiert track_id als Teil von merge_into_track_id (z.B. weil die Zelle
    sich kurz so stark bewegt hat, dass der Tracker einen neuen Track
    begonnen hat - das sähe sonst wie ein neuer Bud aus, ist aber dieselbe
    Zelle). NICHT-DESTRUKTIV: die Rohdaten werden nicht verändert, die
    Korrektur lebt ausschließlich in qc_exclusions.csv. Anwendung beim
    Laden über apply_track_merges() - vor jeglicher Lineage-Klassifikation
    (run_analysis.py: Exclusions, Merges, Exclusions).

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

    frame_from/frame_to (optional) beschraenken den Merge auf diese Frames
    der Quelle. Zwei Aufrufe mit disjunkten Bereichen und verschiedenen
    Zielen teilen einen Track mit Tracker-Sprung auf (siehe Modul-Docstring).
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
        "frame_from": frame_from if frame_from is not None else np.nan,
        "frame_to": frame_to if frame_to is not None else np.nan,
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

    MUSS vor jeglicher Lineage-Klassifikation (lineage.classify_mother_bud)
    aufgerufen werden, sonst sieht die Heuristik weiterhin zwei getrennte
    Tracks und der Tracking-Bruch wird fälschlich als neuer Bud
    interpretiert. run_analysis.py wendet die Exclusions davor UND danach an
    (siehe Modul-Docstring).

    FRAME-BEREICH: frame_from/frame_to auf einer Merge-Zeile beschraenken den
    Merge auf diese Frames der Quelle. So wird ein Track, der nacheinander
    zwei verschiedene Zellen verfolgt hat (Tracker-Sprung), mit zwei Zeilen
    und disjunkten Bereichen auf seine zwei richtigen Ziele aufgeteilt.

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

        # Frame-Bereich auf der Quelle: nur diese Frames wandern ins Ziel.
        # So laesst sich ein Track, der zwei echte Zellen nacheinander
        # verfolgt hat, AUFTEILEN (zwei Merge-Zeilen mit disjunkten Bereichen).
        lo = row["frame_from"] if "frame_from" in row and pd.notna(row["frame_from"]) else -np.inf
        hi = row["frame_to"] if "frame_to" in row and pd.notna(row["frame_to"]) else np.inf
        if lo != -np.inf or hi != np.inf:
            src_mask = src_mask & (df["frame"] >= lo) & (df["frame"] <= hi)
            if not src_mask.any():
                logger.warning(
                    "Track-Merge übersprungen: '%s' hat keine Frames im Bereich %s-%s.",
                    row["cell_uid"], lo, hi,
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
        logger.info("Track-Merge angewendet: '%s' -> '%s'%s", row["cell_uid"], dst_cell_uid,
                    f" (Frames {lo:g}-{hi:g})" if (lo != -np.inf or hi != np.inf) else "")

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


# ==============================================================================
# Konsistenz der QC-Datei und Reichweite des QC
# ==============================================================================

def qc_batches(exclusions: pd.DataFrame) -> pd.DataFrame:
    """Welche Batches (biosensor, osc_type, osc_freq) hat das manuelle QC ueberhaupt beruehrt?

    Das QC ist arbeitsintensiv und deckt in der Praxis nur einen Teil der
    Daten ab (im echten Datensatz: einen einzigen Batch, WT/pH/6). Ein
    Vergleich "mit QC vs. ohne QC" ist nur dort sinnvoll, wo QC stattgefunden
    hat - sonst misst man "kein QC vs. kein QC" und nennt es "QC hat keinen
    Effekt".
    """
    if exclusions is None or exclusions.empty or "cell_uid" not in exclusions.columns:
        return pd.DataFrame(columns=["biosensor", "osc_type", "osc_freq", "n_rows", "n_chambers"])
    parts = exclusions["cell_uid"].astype(str).str.split("__", expand=True)
    if parts.shape[1] < 6:
        return pd.DataFrame(columns=["biosensor", "osc_type", "osc_freq", "n_rows", "n_chambers"])
    df = pd.DataFrame({
        "biosensor": parts[0], "osc_type": parts[1], "osc_freq": parts[2],
        "exp_id": exclusions["cell_uid"].astype(str).str.replace(r"__track\d+$", "", regex=True),
    })
    return (df.groupby(["biosensor", "osc_type", "osc_freq"])
              .agg(n_rows=("exp_id", "size"), n_chambers=("exp_id", "nunique")).reset_index())


def find_qc_conflicts(exclusions: pd.DataFrame) -> pd.DataFrame:
    """Zeilen der QC-Datei, die mehrfach denselben Track betreffen - mit der FOLGE, die
    apply_track_merges()/apply_qc_exclusions() daraus machen.

      split_merge        derselbe Track wird in zwei Ziele gemergt, mit frame_from/frame_to:
                         eine Aufteilung. apply_track_merges() verschiebt je Zeile nur den
                         angegebenen Frame-Bereich - das wird ausgefuehrt.
      double_merge       zwei Ziele OHNE Frame-Bereiche: echter Widerspruch. Erste Zeile
                         gewinnt, zweite wird verworfen.
      merge_and_exclude  gemergt UND ausgeschlossen. run_analysis.py wendet Ausschluesse
                         VOR den Merges (auf die urspruengliche cell_uid) und danach noch
                         einmal an - beides wird ausgefuehrt.
      redundant_exclusion zwei Ausschluss-Zeilen mit verschiedenen Gruenden: harmlos,
                         beide loeschen dieselben Zeilen.
      duplicate          identische Zeile zweimal: harmlos.

    Nichts davon wird hier repariert - die Entscheidung gehoert der Person, die
    die Bilder gesehen hat. Die Tabelle nennt die Zeilen mit Nummer.
    """
    if exclusions is None or exclusions.empty or "cell_uid" not in exclusions.columns:
        return pd.DataFrame()
    df = exclusions.copy()
    df["_row"] = range(len(df))
    has_merge = df["merge_into_track_id"].notna() if "merge_into_track_id" in df.columns else pd.Series(False, index=df.index)
    has_range = pd.Series(False, index=df.index)
    for c in ("frame_from", "frame_to"):
        if c in df.columns:
            has_range |= df[c].notna()
    consequences = {
        "split_merge": "intended split: each row moves only its frame range into its target (executed)",
        "double_merge": "contradiction: first row wins, second row rejected",
        "merge_and_exclude": "exclusion applied BEFORE the merge on the original cell_uid, then the rest merged (executed)",
        "redundant_exclusion": "harmless: both rows exclude the same track",
        "duplicate": "harmless: identical rows",
    }
    records = []
    for uid, grp in df.groupby("cell_uid"):
        if len(grp) < 2:
            continue
        m = has_merge.loc[grp.index]
        merges = grp[m]
        targets = merges["merge_into_track_id"].dropna().unique()
        if len(targets) > 1:
            kind = "split_merge" if has_range.loc[merges.index].any() else "double_merge"
        elif len(merges) and len(merges) < len(grp):
            kind = "merge_and_exclude"
        elif grp.drop(columns="_row").astype(str).drop_duplicates().shape[0] == 1:
            kind = "duplicate"
        else:
            kind = "redundant_exclusion"
        for _, row in grp.iterrows():
            records.append({"kind": kind, "consequence": consequences[kind], "cell_uid": uid,
                            "row_in_file": int(row["_row"]) + 2, "reason": row.get("reason", ""),
                            "merge_into_track_id": row.get("merge_into_track_id", None),
                            "frame_from": row.get("frame_from", None), "frame_to": row.get("frame_to", None)})
    out = pd.DataFrame(records)
    if not out.empty:
        counts = out.groupby("kind")["cell_uid"].nunique().to_dict()
        serious = {k: v for k, v in counts.items() if k in ("double_merge",)}
        logger.warning(
            "QC-DATEI: %d Tracks sind mehrfach gelistet (%s). Davon NICHT wie beabsichtigt ausgefuehrt: %s. "
            "Details und Zeilennummern in 70_qc_exclusions_conflicts.csv.",
            out["cell_uid"].nunique(), ", ".join(f"{k}: {v}" for k, v in counts.items()),
            ", ".join(f"{k}: {v}" for k, v in serious.items()) or "keine",
        )
    return out
