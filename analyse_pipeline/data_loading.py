"""
data_loading.py
================
Einlesen & Zusammenführen aller Combined_Results über die volle
Experiment-Hierarchie:

    Data/<Biosensor>/<Oszillationstyp>/<Oszillationsfrequenz>/03_results/Combined_Results.{csv,parquet}

Beispiel:
    Data/BSG/Glc/6/03_results/Combined_Results.csv
    Data/BSA/pH/24/03_results/Combined_Results.csv

WICHTIG: Die Metadaten (biosensor, osc_type, osc_freq) werden NICHT aus dem
Dateinamen geparst, sondern aus dem ORDNERPFAD relativ zu data_root. Das ist
robuster, weil der Dateiname (YYYYMMDD_Condition_Rep_Cham), den die
Bildverarbeitungs-Pipeline erzeugt, eine 'condition' enthalten kann, die
etwas völlig anderes bedeutet (z.B. Well-Position) als die hier relevante
biologische Hierarchie.

Die Funktion ist bewusst tolerant: fehlt eine Ebene oder ist die Struktur an
einer Stelle anders, gibt es eine klare Warnung statt eines stillen Fehlers;
die Datei wird trotzdem (mit fehlenden Metadaten) eingelesen, damit nichts
unbemerkt verloren geht.
"""

from __future__ import annotations

import hashlib
import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pandas as pd

logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Erwartete Tiefe relativ zu data_root, VOR dem Dateinamen:
#   <biosensor>/<osc_type>/<osc_freq>/<results_subdir>/<dateiname>
EXPECTED_DEPTH = 4


@dataclass
class DiscoveredFile:
    path: Path
    biosensor: Optional[str]
    osc_type: Optional[str]
    osc_freq: Optional[str]
    results_subdir: Optional[str]
    path_ok: bool


def _parse_hierarchy(rel_path: Path) -> DiscoveredFile:
    """Zerlegt einen relativen Pfad in die Hierarchie-Ebenen."""
    parts = rel_path.parts  # z.B. ('BSG', 'Glc', '6', '03_results', 'Combined_Results.csv')
    n = len(parts)
    if n < EXPECTED_DEPTH + 1:
        return DiscoveredFile(rel_path, None, None, None, None, path_ok=False)

    biosensor, osc_type, osc_freq, results_subdir = parts[n - EXPECTED_DEPTH - 1: n - 1]
    return DiscoveredFile(rel_path, biosensor, osc_type, osc_freq, results_subdir, path_ok=True)


def discover_result_files(
    data_root: str | Path,
    filename_pattern: str = "Combined_Results.*",
    results_subdir: str = "",
) -> list[DiscoveredFile]:
    """
    Findet rekursiv alle Ergebnisdateien unterhalb von data_root und
    extrahiert Biosensor / Oszillationstyp / Oszillationsfrequenz aus dem
    Pfad relativ zu data_root.
    """
    data_root = Path(data_root)
    if not data_root.exists():
        raise FileNotFoundError(f"data_root existiert nicht: {data_root}")

    all_files = sorted(
        p for p in data_root.rglob(filename_pattern)
        if p.suffix.lower() in (".csv", ".parquet")
    )

    if not all_files:
        raise FileNotFoundError(
            f"Keine Dateien passend zu '{filename_pattern}' unter '{data_root}' gefunden. "
            f"Pfad und Muster korrekt?"
        )

    discovered = [_parse_hierarchy(p.relative_to(data_root)) for p in all_files]
    # path im DiscoveredFile soll der ABSOLUTE Pfad sein, nicht relativ
    discovered = [
        DiscoveredFile(abs_p, d.biosensor, d.osc_type, d.osc_freq, d.results_subdir, d.path_ok)
        for abs_p, d in zip(all_files, discovered)
    ]
    if results_subdir:
        n_all = len(discovered)
        discovered = [d for d in discovered if d.results_subdir == results_subdir]
        logger.info("Ergebnisordner '%s': %d von %d Dateien.", results_subdir, len(discovered), n_all)
        if not discovered:
            raise FileNotFoundError(
                f"Keine '{filename_pattern}' in einem Ordner '{results_subdir}' unter '{data_root}' "
                f"(AUREO_RESULTS_SUBDIR)."
            )

    bad = [d for d in discovered if not d.path_ok]
    if bad:
        bad_list = "\n".join(f"  - {d.path}" for d in bad)
        logger.warning(
            "%d von %d Dateien hatten nicht die erwartete Ordnertiefe "
            "(biosensor/osc_type/osc_freq/<results_subdir>/Datei). "
            "Metadaten dafür sind None, bitte Pfade prüfen:\n%s",
            len(bad), len(discovered), bad_list,
        )

    return discovered


def _fingerprint(discovered: list[DiscoveredFile]) -> str:
    """Kennung der Quelldateien: Pfad, Groesse, Aenderungszeit - aendert sich, sobald eine Datei dazukommt,
    fehlt oder neu geschrieben wurde."""
    parts = []
    for d in sorted(discovered, key=lambda d: str(d.path)):
        try:
            st = Path(d.path).stat(); parts.append(f"{d.path}|{st.st_size}|{int(st.st_mtime)}")
        except OSError:
            parts.append(f"{d.path}|missing")
    return hashlib.sha1("\n".join(parts).encode()).hexdigest()


def _read_one(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path)


def load_all_results(
    data_root: str | Path,
    cache_path: Optional[str | Path] = None,
    force_reload: bool = False,
    filename_pattern: str = "Combined_Results.*",
    results_subdir: str = "",
) -> pd.DataFrame:
    """
    Lädt und kombiniert ALLE Combined_Results unter data_root zu einem
    einzigen DataFrame, inkl. Hierarchie-Metadaten (biosensor, osc_type,
    osc_freq).

    Zusätzlich erzeugte Spalten:
        source_file : Pfad der Quelldatei (Rückverfolgung bei QC-Fragen)
        exp_id      : eindeutige Experiment-ID über die volle Hierarchie
                      (biosensor, osc_type, osc_freq, condition, replicate, chamber)
        cell_uid    : global eindeutige Zell-Track-ID (exp_id + track_id)

    Parameters
    ----------
    data_root : Pfad zur Data-Wurzel (z.B. ".../Data")
    cache_path : Optional. Wenn gesetzt, wird das kombinierte Ergebnis als
                 .parquet gecacht und bei folgenden Aufrufen (falls vorhanden
                 und force_reload=False) direkt davon geladen - deutlich
                 schneller bei großen Datenmengen.
    force_reload : Cache ignorieren und neu einlesen.
    """
    cache_path = Path(cache_path) if cache_path else None

    discovered = discover_result_files(data_root, filename_pattern, results_subdir=results_subdir)
    # Der Cache gilt nur fuer genau diese Dateien in genau diesem Zustand: Pfade, Groessen und
    # Aenderungszeiten stehen in <cache>.meta.json. Kommt eine Tabelle dazu oder wird eine neu
    # geschrieben, liest der Loader von selbst neu - FORCE_RELOAD ist dann nicht noetig.
    fingerprint = _fingerprint(discovered)
    if cache_path and cache_path.exists() and not force_reload:
        meta_path = cache_path.with_suffix(cache_path.suffix + ".meta.json")
        stored = None
        if meta_path.exists():
            try:
                stored = json.loads(meta_path.read_text()).get("fingerprint")
            except Exception:  # noqa: BLE001
                stored = None
        if stored == fingerprint:
            logger.info("Lade aus Cache (Quelldateien unveraendert): %s", cache_path)
            return pd.read_parquet(cache_path)
        logger.info("Cache %s passt nicht zu den %d gefundenen Dateien (neu, geaendert oder ohne Kennung) - "
                    "Tabellen werden neu eingelesen.", cache_path.name, len(discovered))

    logger.info("Gefundene Ergebnisdateien: %d", len(discovered))
    logger.info("  Biosensoren:            %s", sorted({d.biosensor for d in discovered if d.biosensor}))
    logger.info("  Oszillationstypen:      %s", sorted({d.osc_type for d in discovered if d.osc_type}))
    logger.info("  Oszillationsfrequenzen: %s", sorted({d.osc_freq for d in discovered if d.osc_freq}))

    frames = []
    for d in discovered:
        df = _read_one(d.path)
        df["biosensor"] = d.biosensor
        df["osc_type"] = d.osc_type
        df["osc_freq"] = d.osc_freq
        df["source_file"] = str(d.path)
        frames.append(df)

    combined = pd.concat(frames, ignore_index=True)

    required_for_id = ["biosensor", "osc_type", "osc_freq", "condition", "replicate", "chamber", "track_id"]
    missing = [c for c in required_for_id if c not in combined.columns]
    if missing:
        raise ValueError(
            f"Für exp_id/cell_uid fehlen Spalten: {missing}. "
            f"Kommen 'condition', 'replicate', 'chamber', 'track_id' aus der "
            f"Bildverarbeitungs-Pipeline? (siehe Combined_Results.csv)"
        )

    combined["exp_id"] = (
        combined[["biosensor", "osc_type", "osc_freq", "condition", "replicate", "chamber"]]
        .astype(str)
        .agg("__".join, axis=1)
    )
    combined["cell_uid"] = combined["exp_id"] + "__track" + combined["track_id"].astype(str)

    logger.info(
        "Kombiniert: %d Zeilen, %d eindeutige Experimente (exp_id), %d eindeutige Zell-Tracks (cell_uid)",
        len(combined), combined["exp_id"].nunique(), combined["cell_uid"].nunique(),
    )

    combined = _normalize_mixed_type_columns(combined)

    if cache_path:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        combined.to_parquet(cache_path, index=False)
        meta_path = cache_path.with_suffix(cache_path.suffix + ".meta.json")
        meta_path.write_text(json.dumps({"fingerprint": fingerprint, "n_files": len(discovered),
                                         "files": [str(d.path) for d in discovered]}, indent=1))
        logger.info("Cache geschrieben: %s (%d Quelldateien, Kennung in %s)", cache_path, len(discovered), meta_path.name)

    return combined


def _normalize_mixed_type_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Erzwingt einen konsistenten dtype für object-Spalten, BEVOR nach Parquet
    geschrieben wird.

    Hintergrund: wenn z.B. Datei A in der Spalte 'date' nur Zahlen hat
    (-> beim Einlesen int64) und Datei B in derselben Spalte irgendwo den
    String 'Unknown' enthält (-> beim Einlesen object/str), wird die nach
    pd.concat() entstehende Spalte zu einer gemischten object-Spalte
    (int UND str gemischt). pyarrow kann so eine Spalte nicht in einen
    einheitlichen Parquet-Typ konvertieren und bricht mit ArrowInvalid ab.

    Diese Funktion erkennt solche Spalten (object-dtype mit mehr als einem
    nicht-NA Python-Typ) und wandelt sie konsistent in Strings um - die
    Information bleibt vollständig erhalten (z.B. 20260101 -> "20260101"),
    nur der Spaltentyp wird einheitlich. Reine Zahlen- oder reine
    String-Spalten bleiben unverändert.
    """
    for col in df.columns:
        if df[col].dtype != object:
            continue
        non_na = df[col].dropna()
        if non_na.empty:
            continue
        types_seen = non_na.map(type).unique()
        if len(types_seen) > 1:
            logger.warning(
                "Spalte '%s' hat gemischte Typen (%s) über die eingelesenen Dateien hinweg "
                "- wird einheitlich zu String konvertiert, damit der Parquet-Cache geschrieben "
                "werden kann. Werte bleiben inhaltlich erhalten (z.B. 20260101 -> '20260101').",
                col, [t.__name__ for t in types_seen],
            )
            df[col] = df[col].astype(str).where(df[col].notna(), other=pd.NA)
    return df
