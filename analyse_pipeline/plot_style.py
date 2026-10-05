"""plot_style.py - EIN Aussehen fuer alle Abbildungen der Pipeline.

Regeln (config.STRAIN_COLORS, docs/data_story.md 8):
  * Farbe = Stamm, ueberall. Facetten je Stamm behalten ihre Stammfarbe.
  * Kontrollarten ueber Marker und Fuellung in der Stammfarbe: Oszillation Kreis gefuellt mit
    durchgezogener Linie, Feast-Kontrolle (PosCtrl) Dreieck hoch gefuellt, Famine-Kontrolle
    (NegCtrl) Dreieck runter hohl; das gepoolte Kontrollband ein heller Ton der Stammfarbe.
  * Perioden in Zeitreihen: Hell-Dunkel-Rampe der Stammfarbe (kurz hell, lang dunkel).
  * Statisch: komplexes Medium (ypd) gefuellt, Minimalmedium (omlp) hohl.
  * PKO und Unbekanntes: neutrales Grau.
  * Text in Tintenfarbe, nie in Serienfarbe; Legende unter den Panels; keine Erklaertexte in
    der Abbildung (die gehoeren in die Bildunterschrift der Arbeit).
  * Weisser Hintergrund, nur horizontales, blasses Gitter, linke und untere Achse, Schrift 8-9 pt,
    Marker >= 5 pt mit dunklem Rand (Gelb, Gruen, Rosa liegen unter 3:1 Kontrast auf Weiss),
    Linien 1.5 pt, PDF mit eingebetteten, editierbaren Schriften (fonttype 42).

Aufruf: apply_style() EINMAL am Anfang eines Laufs (run_analysis.py, validate_lineage.py).
"""
from __future__ import annotations

import colorsys
from typing import Iterable, Optional, Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex, to_rgb
from matplotlib.lines import Line2D

from config import STRAIN_COLORS, STRAIN_COLOR_OTHER, STRAIN_ORDER

# Tintenfarben fuer Text, Achsen, Referenzlinien.
INK = "#222222"
INK_SOFT = "#555555"
INK_MUTED = "#8a8a8a"
GRID = "#e4e4e4"
SURFACE = "#ffffff"

# Breiten fuer die Seite der Arbeit (Zoll): eine Spalte / volle Breite.
FIG_W_SINGLE = 3.35
FIG_W_FULL = 6.7

# Kontrollarten: Marker und Fuellung, Farbe kommt vom Stamm.
CONTROL_MARKERS = {"Oscillation": "o", "PosCtrl": "^", "NegCtrl": "v"}
CONTROL_FILLED = {"Oscillation": True, "PosCtrl": True, "NegCtrl": False}
CONTROL_LABELS = {"Oscillation": "oscillation chambers", "PosCtrl": "feast control (PosCtrl)",
                  "NegCtrl": "famine control (NegCtrl)"}
CONTROL_LINESTYLES = {"Oscillation": "-", "PosCtrl": (0, (4, 2)), "NegCtrl": (0, (1.5, 1.5))}
# Statisch: Medium ueber die Fuellung.
MEDIUM_FILLED = {"ypd": True, "omlp": False}

_RC = {
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "figure.dpi": 110, "savefig.dpi": 200, "savefig.bbox": "tight", "savefig.pad_inches": 0.04,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9, "xtick.labelsize": 8,
    "ytick.labelsize": 8, "legend.fontsize": 8, "legend.title_fontsize": 8, "figure.titlesize": 10,
    "axes.titleweight": "normal", "axes.titlepad": 5,
    "axes.edgecolor": INK_SOFT, "axes.linewidth": 0.8, "axes.labelcolor": INK,
    "xtick.color": INK_SOFT, "ytick.color": INK_SOFT, "xtick.labelcolor": INK, "ytick.labelcolor": INK,
    "xtick.major.size": 3, "ytick.major.size": 3, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.minor.visible": False, "ytick.minor.visible": False,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "axes.grid.axis": "y", "grid.color": GRID, "grid.linewidth": 0.6, "grid.linestyle": "-",
    "axes.axisbelow": True,
    "lines.linewidth": 1.5, "lines.markersize": 5, "lines.markeredgewidth": 0.7,
    "errorbar.capsize": 2.5,
    "legend.frameon": False, "legend.handlelength": 1.6, "legend.handletextpad": 0.5,
    "legend.columnspacing": 1.2, "legend.borderaxespad": 0.3,
    "text.color": INK, "axes.unicode_minus": True,
    "hatch.linewidth": 0.6,
}


def apply_style() -> None:
    """rcParams der Pipeline setzen (ueberschreibt seaborn-Themes)."""
    import logging
    logging.getLogger("fontTools").setLevel(logging.WARNING)   # PDF-Schrifteinbettung (fonttype 42) ist sonst geschwaetzig
    mpl.rcParams.update(_RC)
    strains = [s for s in STRAIN_ORDER if s in STRAIN_COLORS]
    mpl.rcParams["axes.prop_cycle"] = mpl.cycler(color=[STRAIN_COLORS[s] for s in strains] + [STRAIN_COLOR_OTHER])


# --------------------------------------------------------------------------------------
# Farben
# --------------------------------------------------------------------------------------
def strain_color(strain: object) -> str:
    """Stammfarbe; PKO und Unbekanntes grau."""
    return STRAIN_COLORS.get(str(strain), STRAIN_COLOR_OTHER)


def mix(color: str, other: str, t: float) -> str:
    """Linear zwischen color (t = 0) und other (t = 1) mischen."""
    a, b = to_rgb(color), to_rgb(other)
    t = min(max(float(t), 0.0), 1.0)
    return to_hex(tuple(a[i] * (1 - t) + b[i] * t for i in range(3)))


def tint(color: str, amount: float = 0.6) -> str:
    """Hellerer Ton (Richtung Weiss); amount 0 = unveraendert, 1 = weiss."""
    return mix(color, "#ffffff", amount)


def shade(color: str, amount: float = 0.35) -> str:
    """Dunklerer Ton (Richtung Schwarz)."""
    return mix(color, "#000000", amount)


def edge_color(color: str) -> str:
    """Randfarbe eines Markers: die Stammfarbe deutlich abgedunkelt (Kontrast-Relief)."""
    return shade(color, 0.45)


def strain_ramp(strain: object, n: int) -> list[str]:
    """n Stufen der Stammfarbe von hell nach dunkel (kurze -> lange Periode).

    Die Stufen laufen ueber die Helligkeit (HLS) bei fester Farbtoenung, damit die Reihenfolge
    auch fuer Farbsehschwache lesbar bleibt; die Endpunkte bleiben von Weiss und Schwarz weg.
    """
    base = strain_color(strain)
    if n <= 1:
        return [base]
    h, l, s = colorsys.rgb_to_hls(*to_rgb(base))
    light, dark = min(0.82, l + 0.30), max(0.22, l - 0.25)
    out = []
    for i in range(n):
        li = light + (dark - light) * i / (n - 1)
        out.append(to_hex(colorsys.hls_to_rgb(h, li, s)))
    return out


def ordered_strains(strains: Iterable[object]) -> list[str]:
    """Staemme in der festen Reihenfolge config.STRAIN_ORDER, Unbekanntes alphabetisch dahinter."""
    present = {str(s) for s in strains}
    known = [s for s in STRAIN_ORDER if s in present]
    return known + sorted(present - set(known))


# --------------------------------------------------------------------------------------
# Marker
# --------------------------------------------------------------------------------------
def marker_kwargs(condition_type: str, color: str, size: float = 30.0, filled: Optional[bool] = None) -> dict:
    """kwargs fuer ax.scatter(): Marker/Fuellung nach Kontrollart, Farbe vom Stamm."""
    ct = str(condition_type)
    if filled is None:
        filled = CONTROL_FILLED.get(ct, True)
    return dict(
        marker=CONTROL_MARKERS.get(ct, "o"), s=size,
        facecolor=color if filled else SURFACE,
        edgecolor=edge_color(color), linewidth=0.8, zorder=3,
    )


def errorbar_kwargs(condition_type: str, color: str, filled: Optional[bool] = None) -> dict:
    """kwargs fuer ax.errorbar() mit Marker und Linie in der Stammfarbe."""
    ct = str(condition_type)
    if filled is None:
        filled = CONTROL_FILLED.get(ct, True)
    return dict(
        marker=CONTROL_MARKERS.get(ct, "o"), markersize=5.5, linewidth=1.4, capsize=2.5,
        color=color, markerfacecolor=color if filled else SURFACE, markeredgecolor=edge_color(color),
        markeredgewidth=0.8, ecolor=color, linestyle=CONTROL_LINESTYLES.get(ct, "-"), zorder=4,
    )


def control_handles(color: str = INK_SOFT, which: Sequence[str] = ("Oscillation", "PosCtrl", "NegCtrl")) -> list:
    """Legendeneintraege fuer die Kontrollarten (Marker/Fuellung), in neutraler Farbe."""
    handles = []
    for ct in which:
        filled = CONTROL_FILLED.get(ct, True)
        handles.append(Line2D([], [], marker=CONTROL_MARKERS[ct], linestyle="", markersize=6,
                              markerfacecolor=color if filled else SURFACE, markeredgecolor=color,
                              markeredgewidth=0.9, label=CONTROL_LABELS[ct]))
    return handles


def strain_handles(strains: Iterable[object], marker: str = "s") -> list:
    """Legendeneintraege fuer die Staemme (Farbfelder)."""
    return [Line2D([], [], marker=marker, linestyle="", markersize=7, markerfacecolor=strain_color(s),
                   markeredgecolor=edge_color(strain_color(s)), markeredgewidth=0.6, label=str(s))
            for s in ordered_strains(strains)]


# --------------------------------------------------------------------------------------
# Layout
# --------------------------------------------------------------------------------------
def legend_below(fig, handles, labels=None, ncol: Optional[int] = None, y: float = 0.0, title: Optional[str] = None):
    """Eine Legende fuer die ganze Abbildung, zentriert unter den Panels."""
    if labels is None:
        labels = [h.get_label() for h in handles]
    if not handles:
        return None
    ncol = ncol or min(len(handles), 4)
    return fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, y), ncol=ncol,
                      frameon=False, title=title)


def dedupe_handles(axes) -> tuple[list, list]:
    """Handles/Labels aller Achsen einsammeln, Duplikate (gleiches Label) entfernen."""
    seen: dict[str, object] = {}
    for ax in (axes.flat if hasattr(axes, "flat") else axes):
        h, l = ax.get_legend_handles_labels()
        for hh, ll in zip(h, l):
            if ll and ll not in seen:
                seen[ll] = hh
    return list(seen.values()), list(seen.keys())


AXIS_LABELS = {"osc_freq": "half-cycle period [min]", "medium": "medium", "chip_family": "chip family",
               "biosensor": "strain", "condition_type": "chamber type"}


def axis_label(col: str) -> str:
    """Achsen-/Legendentitel fuer eine Spalte (config-Spaltennamen -> lesbarer Text)."""
    return AXIS_LABELS.get(str(col), str(col).replace("_", " "))


def panel_title(ax, text: str, loc: str = "left") -> None:
    ax.set_title(text, loc=loc, fontsize=mpl.rcParams["axes.titlesize"], color=INK)


def finish(fig, out_path, logger=None) -> None:
    """tight_layout, speichern, schliessen - und den Dateinamen loggen."""
    try:
        fig.tight_layout()
    except Exception:  # pragma: no cover - tight_layout kann bei leeren Achsen warnen
        pass
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    if logger is not None:
        logger.info("Plot gespeichert: %s", getattr(out_path, "name", out_path))
