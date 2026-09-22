"""
pipeline_steps.py
=================
Der Analyse-Kern als EINZELNE, benannte Schritte statt als eine Funktion.

WARUM
-----
Vorher lag der gesamte Ablauf (Schritte 00-95) in einer 555-Zeilen-Funktion
run_pipeline(). Das hatte drei praktische Folgen: man konnte keinen einzelnen
Schritt erneut laufen lassen, ohne alles neu zu rechnen; man konnte keinen
Schritt testen; und man musste die Funktion von oben lesen, um irgendetwas zu
finden. Hier ist jeder Schritt eine eigene Funktion mit einem Schluessel
("40_robustness"), und run_analysis.py kann eine Teilmenge davon ausfuehren:

    python run_analysis.py --list-steps
    python run_analysis.py --steps 20 40

BERECHNUNG vs. AUSGABE
----------------------
Die Zwischenergebnisse, die mehrere Schritte brauchen (Lineage-Events,
Mutterliste, µ-Tabellen), haengen am PipelineContext und werden LAZY berechnet:
beim ersten Zugriff einmal, danach aus dem Cache. Dadurch funktioniert

    --steps 40

auch allein - der Kontext rechnet die noetigen Vorstufen selbst nach, ohne dass
Schritt 10 seine Dateien noch einmal schreiben muss. Bei einem vollstaendigen
Lauf wird trotzdem jede Vorstufe genau einmal berechnet.

Die Schrittnummern entsprechen den Datei-Praefixen im Output-Ordner, damit
Ablauf und Ordnerinhalt dieselbe Reihenfolge haben.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

import pandas as pd

from config import (
    MIN_PER_FRAME,
    OSCILLATION_START_MIN,
    PANEL_A_GROUP_COL,
    PANEL_A_FACET_COL,
    CONTROL_CONCENTRATION_LABELS,
    LINEAGE_PARAMS,
    STABLE_MOTHER_MIN_COVERAGE,
    STABLE_MOTHER_GROUP_COLS,
    STABLE_MOTHER_BASE_VALUE_COLS,
    FLUX_CONFIG,
    MU_MAX_THRESHOLD,
    ROBUSTNESS_VALUE_COLS,
    ENDPOINT_LAST_FRACTION,
    ENDPOINT_VALUE_COLS,
    STATIC_SATURATION_LEVEL,
    STATIC_SATURATION_SMOOTH_FRAMES,
    STATIC_SATURATION_MIN_GROWTH,
    AREA_GROWTH_MIN_FRAMES,
)
from experiment_units import summarise_hierarchical
from endpoint_trends import (
    compute_endpoint_per_replicate,
    summarise_per_replicate,
    spearman_against_period,
    plot_endpoint_vs_period,
    bracket_normalise,
    control_trend_check,
    within_culture_trend,
    plot_control_trend_summary,
    detect_saturation_frame,
    static_endpoint_window,
)
from sensors import SENSOR_CONFIG
from lineage import classify_mother_bud, compute_budding_ratio, identify_mothers
from mother_trajectories import plot_stable_mother_per_group, build_lineage_tree, summarise_lineage_depth
from budding_ratio_timeseries import compute_budding_ratio_timeseries, aggregate_budding_ratio_over_replicates
from growth_rate import compute_specific_growth_rate, summarise_growth_rate, classify_condition_type
from area_growth import compute_area_growth_rate, summarise_area_growth, plot_mu_event_vs_mu_area
from robustness import (
    compute_rt_population,
    compute_rt_single_cell,
    compute_rp,
    aggregate_robustness_over_replicates,
)
from control_consistency import test_control_consistency_across_freq, plot_control_consistency
from queen_controls import (
    prepare_sensor_controls,
    summarise_sensor_controls,
    plot_sensor_control_timeseries,
    plot_sensor_control_comparison,
    plot_sensor_raw_channel_timeseries,
    summarise_sensor_controls_by_chamber,
    plot_control_chamber_comparison,
)
from violin_plots import plot_panel_a
from summary_plots import (
    plot_budding_ratio_timeseries,
    plot_point_errorbar,
    plot_rt_vs_rp_quadrant,
    plot_rt_single_cell_distribution,
)
from analysis import (
    pretty_label,
    data_overview,
    plot_n_tracks_overview,
    plot_metric_over_time_by_frequency,
    plot_morphology_scatter,
    plot_single_cell_trajectories,
    summary_statistics,
)

logger = logging.getLogger(__name__)


def exclude_controls(df: pd.DataFrame) -> pd.DataFrame:
    """
    Entfernt PosCtrl-/NegCtrl-Zeilen aus einer Tabelle für die 'normalen'
    Plots - die Kontrollen selbst werden bereits eigenständig über
    test_control_consistency_across_freq()/plot_control_consistency()
    ausgewertet (Suffix '_control_consistency_*'), daher sollen sie in den
    übrigen Plots nicht mehr doppelt auftauchen.

    Nutzt bevorzugt eine vorhandene 'condition_type'-Spalte (z.B. aus
    summarise_growth_rate() oder den Robustness-Aggregaten); ist die nicht
    vorhanden, wird sie aus 'condition' on-the-fly berechnet (siehe
    growth_rate.classify_condition_type()). Fehlen beide Spalten, wird die
    Tabelle unverändert zurückgegeben (mit Warnung), statt fälschlich alles
    zu verwerfen. Wirkt NUR auf das, was geplottet wird - CSV-Exporte und
    die Kontroll-Konsistenz-Berechnungen selbst bleiben unverändert und
    behalten weiterhin PosCtrl/NegCtrl.
    """
    if df is None or df.empty:
        return df
    if "condition_type" in df.columns:
        ctype = df["condition_type"]
    elif "condition" in df.columns:
        ctype = df["condition"].apply(classify_condition_type)
    else:
        logger.warning(
            "exclude_controls(): weder 'condition' noch 'condition_type' in der Tabelle "
            "vorhanden - wird ungefiltert (inkl. evtl. Kontrollen) zurückgegeben."
        )
        return df
    return df[~ctype.isin(["PosCtrl", "NegCtrl"])].copy()


def _per_facet(ctx, cells_plot, controls):
    """Zerlegt die Zeitreihen-Plots des Oszillationszweigs in eine Datei pro
    osc_type (Glc / pH): fuenf Staemme in einer Reihe statt eines 5 x 2-Gitters.
    Bei den statischen Daten (facet_col = chip_family, zwei Werte) bleibt es bei
    einer Datei - dort ist das Gitter klein."""
    if ctx.facet_col != "osc_type" or cells_plot is None or cells_plot.empty:
        yield "", cells_plot, controls
        return
    for value in sorted(cells_plot[ctx.facet_col].dropna().unique()):
        sel = cells_plot[cells_plot[ctx.facet_col] == value]
        csel = controls[controls[ctx.facet_col] == value] if controls is not None else None
        yield f"_{value}", sel, csel


@dataclass
class PipelineContext:
    """Alles, was die Schritte gemeinsam brauchen - Eingaben und Zwischenergebnisse.

    Die Zwischenergebnisse sind Properties und werden erst beim ersten Zugriff
    berechnet (und dann gecacht). Deshalb kann jeder Schritt einzeln laufen,
    ohne dass der Aufrufer wissen muss, was er voraussetzt.
    """

    cells: pd.DataFrame
    output_dir: Path
    freq_order: Optional[list[str]]
    intensity_cols: list[str]
    ratio_cols: list[str]
    run_sensor_controls: bool
    # Kontroll-Konsistenz ueber die osc_freq-Batches (control_consistency.py).
    # Braucht MINDESTENS ZWEI Batches: test_control_consistency_across_freq()
    # liefert bei einem einzigen p = NaN, und plot_control_consistency() zeichnet
    # dann einen einzelnen Punkt pro Kontrollart - eine trivial flache Linie, die
    # wie "Kontrollen sind konsistent" aussieht und nichts enthaelt. Fuer solche
    # Teilmengen (z.B. der PKO-Zweig mit nur einer Periode) deshalb False, statt
    # eine irrefuehrende Abbildung zu erzeugen. run_analysis.py setzt das anhand
    # der tatsaechlich vorhandenen Batches.
    run_control_consistency: bool = True
    # Gruppenspalte fuer Panel A. Per Default config.PANEL_A_GROUP_COL
    # ('biosensor'); der statische Kontext setzt 'osc_freq', weil dort die
    # Vergleichsgruppe (static_omlp/static_ypd) in osc_freq steht und
    # plot_panel_a() diese Spalte sonst nie sieht - beide Medien landeten dann
    # in EINEM Violin. Siehe config.PANEL_A_GROUP_COL_STATIC.
    panel_a_group_col: str = PANEL_A_GROUP_COL
    # x-Achse und Facette der Punkt/Errorbar-Plots. Oszillation: Periode x
    # osc_type. Statisch: Medium x Chip-Familie - dort ist 'osc_freq' die
    # Chip-Familie und traegt bei W65 BEIDE Medien, weshalb plot_point_errorbar()
    # mit x_col='osc_freq' auf "mehrdeutige osc_freq-Werte" lief und die
    # Schritte 10 und 40 fuer die statischen Daten komplett ausfielen.
    x_col: str = "osc_freq"
    facet_col: str = "osc_type"
    # Facette fuer Panel A und die Budding-Ratio-Zeitreihe. Oszillation:
    # osc_type. Statisch: Chip-Familie (W109/W65), damit die beiden Chips
    # nebeneinander stehen statt in einer Facette 'Static' zu verschwinden.
    panel_a_facet_col: str = PANEL_A_FACET_COL
    # Groessenkriterium der Mutter/Bud-Heuristik (bud_size.py): ein Kandidat
    # zaehlt nur als Knospe, wenn seine Flaeche beim ersten Auftreten hoechstens
    # diesen Anteil der Mutterflaeche hat - alles darueber ist eine angespuelte
    # Zelle. run_analysis.py leitet den Wert EINMAL aus allen Daten nach QC ab
    # und gibt ihn an alle Kontexte weiter (auch no_qc), damit die Zweige
    # vergleichbar bleiben. None = kein Groessenfilter.
    bud_size_threshold: Optional[float] = None
    # Von run_steps() gefuellt: Schluessel der Schritte, die eine Exception
    # geworfen haben. Ein einzelner fehlgeschlagener Plot soll den Rest des
    # Laufs nicht mitreissen, aber auch nicht unbemerkt bleiben.
    failed_steps: list[str] = field(default_factory=list)
    _cache: dict = field(default_factory=dict, repr=False)

    def _lazy(self, key: str, compute: Callable[[], object]):
        if key not in self._cache:
            logger.info("Zwischenergebnis '%s' wird berechnet ...", key)
            self._cache[key] = compute()
        return self._cache[key]

    @property
    def cells_plot(self) -> pd.DataFrame:
        """Zelldaten OHNE PosCtrl/NegCtrl - nur fuer die 'normalen' Plots.

        Lineage-Klassifikation, µ-/Robustness-Tabellen, Kontroll-Konsistenz und
        alle CSV-Exporte laufen weiterhin auf dem VOLLEN Datensatz; nur die
        Sichtbarkeit in den Plots aendert sich.
        """
        return self._lazy("cells_plot", lambda: exclude_controls(self.cells))

    @property
    def cells_lineage(self) -> pd.DataFrame:
        """Zellen im Sparse-Phase-Fenster (relink.py, Spalte in_lineage_window) -
        die EINZIGE Eingabe der Mutter/Bud-Heuristik und alles, was daran
        haengt (Budding Ratio, µ_event, Stammbaum, Panel A). Im vollen
        Bildfeld sind neu auftauchende Tracks Fragmente, keine Knospen.
        Ohne die Spalte: alle Zellen."""
        def select() -> pd.DataFrame:
            if "in_lineage_window" not in self.cells.columns:
                return self.cells
            sub = self.cells[self.cells["in_lineage_window"].astype(bool)]
            logger.info(
                "Lineage-Fenster: %d von %d Zellzeilen, %d von %d Kammern.",
                len(sub), len(self.cells), sub["exp_id"].nunique(), self.cells["exp_id"].nunique(),
            )
            return sub
        return self._lazy("cells_lineage", select)

    @property
    def lineage_events(self) -> pd.DataFrame:
        """Budding-Events. HEURISTIK ohne echtes Lineage-Tracking - siehe
        lineage.py und validate_lineage.py."""
        return self._lazy(
            "lineage_events",
            lambda: classify_mother_bud(
                self.cells_lineage, LINEAGE_PARAMS, flux_config=FLUX_CONFIG,
                bud_size_threshold=self.bud_size_threshold,
            ),
        )

    @property
    def mothers(self) -> pd.DataFrame:
        """ALLE Mutterzellen, auch ohne Budding-Event - notwendig, damit ruhende
        Muetter mit budding_ratio=0 in Panel A erscheinen."""
        return self._lazy("mothers", lambda: identify_mothers(self.cells_lineage, LINEAGE_PARAMS))

    @property
    def mu_table(self) -> pd.DataFrame:
        return self._lazy("mu_table", lambda: compute_specific_growth_rate(
            self.lineage_events, self.cells_lineage,
            min_per_frame=MIN_PER_FRAME, mu_max_threshold=MU_MAX_THRESHOLD,
        ))

    @property
    def mu_summary(self) -> pd.DataFrame:
        return self._lazy("mu_summary", lambda: (
            summarise_growth_rate(self.mu_table) if not self.mu_table.empty else pd.DataFrame()
        ))

    @property
    def area_table(self) -> pd.DataFrame:
        return self._lazy("area_table", lambda: compute_area_growth_rate(
            self.cells, min_per_frame=MIN_PER_FRAME,
            lineage_events=self.lineage_events, mothers=self.mothers,
            min_frames=AREA_GROWTH_MIN_FRAMES,
        ))



def step_00_overview(ctx: PipelineContext) -> None:
    """Übersicht / Sanity-Check."""
    cells = ctx.cells
    output_dir = ctx.output_dir
    freq_order = ctx.freq_order

    # ==================================================================
    # 00. Übersicht / Sanity-Check
    # ==================================================================
    overview = data_overview(cells)
    overview.to_csv(output_dir / "00_data_overview.csv", index=False)
    logger.info("Übersicht gespeichert: 00_data_overview.csv")
    plot_n_tracks_overview(overview, output_dir / "00_n_tracks_overview.pdf", freq_order=freq_order)



def step_10_growth(ctx: PipelineContext) -> None:
    """Zellfläche, µ_event und µ_area."""
    cells = ctx.cells
    cells_plot = ctx.cells_plot
    output_dir = ctx.output_dir
    freq_order = ctx.freq_order
    mu_table = ctx.mu_table
    mu_summary = ctx.mu_summary
    area_table = ctx.area_table

    # ==================================================================
    # 10. Zellmorphologie & Wachstum (Fläche, µ_event, µ_area)
    # ==================================================================
    if "area" in cells.columns:
        # reference_cells = die von cells_plot entfernten Kontrollen. Ohne sie
        # stünden die Oszillationskurven ohne Bezugsrahmen da, obwohl genau der
        # Vergleich gegen PosCtrl/NegCtrl die Aussage der Abbildung ist.
        controls_only = cells[~cells.index.isin(cells_plot.index)] if not cells_plot.empty else cells
        for facet_value, plot_sel, ctrl_sel in _per_facet(ctx, cells_plot, controls_only):
            plot_metric_over_time_by_frequency(
                plot_sel, "area", output_dir / f"10_cell_area_over_time{facet_value}.pdf",
                freq_order=freq_order, ylabel="Cell area [px²]",
                reference_cells=ctrl_sel, x_col=ctx.x_col, facet_col=ctx.facet_col,
            )

    # -- µ_event: spezifische Wachstumsrate nach Eq. 2 (Blöbaum et al. 2024):
    #    µ = ln(2)/t, aus der Zeit zwischen aufeinanderfolgenden Budding-
    #    Events EINER Mutterzelle. Läuft identisch für alle Bedingungen
    #    inkl. PosCtrl (durchgehend Feast) und NegCtrl (durchgehend
    #    Starvation) - µ ist pro Zelle definiert, nicht pro Bedingung.
    #
    #    style_col="condition_type" wird hier NICHT mehr gesetzt: seit
    #    exclude_controls() PosCtrl/NegCtrl vorher rausfiltert, bleibt in
    #    diesem Plot ohnehin nur noch "Oscillation" übrig (bzw. bei den
    #    statischen Daten St.omlp/St.ypd, was schon über die x-Achse
    #    sichtbar ist) - eine zusätzliche Marker-Form dafür wäre redundant.
    if not mu_table.empty:
        mu_table.to_csv(output_dir / "11_specific_growth_rate.csv", index=False)
        logger.info("Tabelle gespeichert: 11_specific_growth_rate.csv")

        mu_summary.to_csv(output_dir / "11_specific_growth_rate_summary.csv", index=False)
        logger.info("Tabelle gespeichert: 11_specific_growth_rate_summary.csv")

        plot_point_errorbar(
            exclude_controls(mu_summary), value_col="mean_mu", out_path=output_dir / "11_specific_growth_rate.pdf",
            x_col=ctx.x_col, facet_col=ctx.facet_col, color_col=PANEL_A_GROUP_COL,
            x_order=freq_order,
            ylabel="Specific growth rate µ [h⁻¹]",
        )

        # Kontroll-Konsistenz: driften PosCtrl/NegCtrl über die Frequenz-Batches
        # hinweg, obwohl sie eigentlich unabhängig von der Oszillation sein
        # sollten? Test auf Replikat-Ebene (mu_table), Plot auf aggregierter
        # Ebene (mu_summary) - siehe control_consistency.py Docstring.
        if ctx.run_control_consistency:
            mu_control_test = test_control_consistency_across_freq(mu_table, value_col="mu")
            if not mu_control_test.empty:
                mu_control_test.to_csv(output_dir / "11_control_consistency_mu_kruskal.csv", index=False)
                logger.info("Tabelle gespeichert: 11_control_consistency_mu_kruskal.csv")

            plot_control_consistency(
                mu_summary, value_col="mean_mu", out_path=output_dir / "11_control_consistency_mu.pdf",
                x_order=freq_order, ylabel="Specific growth rate µ [h⁻¹] (controls)",
            )
        else:
            logger.info(
                "Kontroll-Konsistenz fuer µ uebersprungen (run_control_consistency=False) - "
                "siehe PipelineContext."
            )
    else:
        logger.warning("Keine µ-Werte berechnet - benötigt Mutterzellen mit >= 2 Budding-Events.")

    # -- µ_area: flächenbasierte Wachstumsrate
    if not area_table.empty:
        area_table.to_csv(output_dir / "12_area_growth_rate.csv", index=False)

        # Separate Summaries für Mütter und Knospen
        area_summary_all = summarise_area_growth(area_table)
        area_summary_mother = summarise_area_growth(area_table, cell_type="mother")
        area_summary_bud = summarise_area_growth(area_table, cell_type="bud")

        area_summary_all.to_csv(output_dir / "12_area_growth_rate_summary_all.csv", index=False)
        area_summary_mother.to_csv(output_dir / "12_area_growth_rate_summary_mother.csv", index=False)
        area_summary_bud.to_csv(output_dir / "12_area_growth_rate_summary_bud.csv", index=False)

        # Plot: µ_area für Mütter (primär interessant) - ohne PosCtrl/NegCtrl,
        # siehe exclude_controls().
        plot_point_errorbar(
            exclude_controls(area_summary_mother), value_col="mean_mu_area",
            out_path=output_dir / "12_area_growth_rate_mother.pdf",
            x_col=ctx.x_col, facet_col=ctx.facet_col, color_col=PANEL_A_GROUP_COL,
            x_order=freq_order,
            ylabel="µ_area, mother cells [h⁻¹]",
        )

        # -- µ_area über ALLE Zellen, Fehlerbalken über REPLIKATE.
        #
        # Zwei Unterschiede zum Mutter-Plot darüber, und beide sind der Grund
        # für diese zusätzliche Abbildung:
        #
        #  1. KEIN cell_type-Filter. Der µ_area-Fit ist eine Regression über
        #     ln(Fläche) EINES Tracks und damit von lineage.py unabhängig -
        #     nur das cell_type-Label hängt an der unvalidierten Mutter/Bud-
        #     Heuristik. Ein Filter auf "mother" schleppt sie ohne Gegenwert
        #     ein; ohne Filter ist das hier die einzige Wachstumsabbildung der
        #     Pipeline, die ganz ohne die Heuristik auskommt.
        #  2. Fehlerbalken über biologische Replikate statt über Zellen.
        #     summarise_area_growth() poolt einzelne Tracks, sd_mu_area
        #     beschreibt also die Streuung über Zellen (Pseudoreplikation,
        #     siehe config.METHOD_CAVEATS). summarise_per_replicate() mittelt
        #     dreistufig: Kammer -> Replikat -> Bedingung.
        #
        # Bewusst OHNE fit_is_reliable-Filter (anders als
        # summarise_area_growth(exclude_unreliable=True)): bei einer flachen
        # Bedingung ist das R² niedrig, WEIL es nichts zu erklären gibt - der
        # Filter behielte dort nur die Tracks, in denen Rauschen wie ein Trend
        # aussieht, und verzerrte den Median nach oben. Statt zu filtern wird
        # der Anteil zuverlässiger Fits als eigene Spalte mitgeschrieben.
        area_rep, area_rep_summary = summarise_per_replicate(area_table, "mu_area")
        if not area_rep_summary.empty:
            # Kammer-Ebene zusaetzlich: Grundlage fuer den QC-Vergleich (qc_comparison.py),
            # der Kammern paarweise mit und ohne QC gegenueberstellt.
            area_chamber, _, _ = summarise_hierarchical(area_table, "mu_area")
            area_chamber.to_csv(output_dir / "12_area_growth_rate_per_chamber.csv", index=False)
            if "fit_is_reliable" in area_table.columns:
                frac = (area_table.groupby(
                            [c for c in ["biosensor", "osc_type", "osc_freq", "condition"]
                             if c in area_table.columns], dropna=False)["fit_is_reliable"]
                        .mean().reset_index(name="frac_fit_reliable"))
                area_rep_summary = area_rep_summary.merge(frac, how="left")
            area_rep.to_csv(output_dir / "12_area_growth_rate_per_chip.csv", index=False)
            area_rep_summary.to_csv(output_dir / "12_area_growth_rate_summary_per_chip.csv",
                                    index=False)
            logger.info("Tabellen gespeichert: 12_area_growth_rate_per_chip.csv / "
                        "_summary_per_chip.csv")

            if ctx.x_col == "osc_freq":
                # Laufen die Kontrollen der Strukturen mit der Periode? (wie Schritt 13)
                ctrl_trend = control_trend_check(area_rep)
                if not ctrl_trend.empty:
                    ctrl_trend.to_csv(output_dir / "12_area_growth_rate_control_trend.csv", index=False)
                    logger.info("Tabelle gespeichert: 12_area_growth_rate_control_trend.csv\n%s",
                                ctrl_trend[["biosensor", "osc_type", "rho_osc", "rho_ctrl_strongest",
                                            "rho_osc_minus_ctrl", "verdict"]].round(2).to_string(index=False))
                within = within_culture_trend(area_rep)
                if not within.empty:
                    within.to_csv(output_dir / "12_area_growth_rate_within_culture.csv", index=False)

            # Kontrollen bleiben IM Plot (Marker-Form = Kontrollart): ohne sie
            # ist ein Periodentrend nicht von einem Struktureffekt zu unterscheiden.
            plot_point_errorbar(
                area_rep_summary, value_col="mean", sd_col="sem",
                out_path=output_dir / "12_area_growth_rate_all.pdf",
                x_col=ctx.x_col, facet_col=ctx.facet_col, color_col=PANEL_A_GROUP_COL,
                style_col="condition_type" if "condition_type" in area_rep_summary.columns else None,
                x_order=freq_order,
                ylabel="µ_area, all cells [h⁻¹]",
                title="µ_area over all cells, oscillation and control chambers — mean ± SEM "
                      "(error unit per table: error_unit)\n"
                      "(no mother/bud filter, so independent of the lineage heuristic)",
            )

            area_trend = spearman_against_period(area_rep)
            if not area_trend.empty:
                area_trend.to_csv(output_dir / "12_area_growth_rate_spearman.csv", index=False)
                logger.info("Tabelle gespeichert: 12_area_growth_rate_spearman.csv\n%s",
                            area_trend.to_string(index=False))

        # Scatter: µ_event vs. µ_area (ebenfalls ohne Kontrollen)
        if not mu_summary.empty:
            plot_mu_event_vs_mu_area(
                exclude_controls(mu_summary), exclude_controls(area_summary_mother),
                output_dir / "12_mu_event_vs_mu_area.pdf",
                label_col=ctx.x_col, facet_col=PANEL_A_GROUP_COL,
                mu_event_table=exclude_controls(mu_table), mu_area_table=exclude_controls(area_table),
            )



def step_13_endpoint(ctx: PipelineContext) -> None:
    """Kumulativer Endzustand gegen die Zyklusperiode (+ Spearman-Trendtest)."""
    cells = ctx.cells
    output_dir = ctx.output_dir

    # ==================================================================
    # 13. Kumulativer Endzustand gegen die Periode
    #
    # Die einzige Auswertungsform, die bei dieser Abtastung interpretierbar
    # ist: einzelne Zyklen liegen unter dem Nyquist-Limit (siehe config.py),
    # ein Zyklusverlauf ist also nicht beobachtbar - die Wirkung ueber
    # Stunden dagegen schon. Vorher gab es dafuer keine Abbildung: 10_/30_/31_
    # sind Zeitreihen, 40_ sind Varianzmasse, und 50_summary_statistics.csv
    # bekommt nur intensity_cols uebergeben, sieht die ratio_*-Spalten also
    # nie. Fuer die Sensor-Daten ist das hier die erste kumulative Auswertung
    # ueberhaupt.
    #
    # Getestet wird mit SPEARMAN gegen die Periode, weil die Vorhersage
    # monoton ist (laengere Famine-Halbzyklen belasten mehr) - nicht mit
    # Kruskal-Wallis, das nur "irgendeine Gruppe unterscheidet sich" prueft.
    # ==================================================================
    value_cols = [c for c in list(dict.fromkeys(list(ENDPOINT_VALUE_COLS) + ctx.ratio_cols))
                  if c in cells.columns]
    if not value_cols:
        logger.warning("Schritt 13: keine der Spalten %s in den Daten - uebersprungen.",
                       list(ENDPOINT_VALUE_COLS) + ctx.ratio_cols)
        return

    # Statischer Zweig: absolutes Fenster VOR der Saettigung (config.py,
    # STATIC_SATURATION_*), sonst relatives Fenster (letzte 25 % je Kammer).
    frame_window = None
    if ctx.x_col != "osc_freq":
        saturation = detect_saturation_frame(
            cells, level=STATIC_SATURATION_LEVEL, smooth_frames=STATIC_SATURATION_SMOOTH_FRAMES,
            min_growth=STATIC_SATURATION_MIN_GROWTH,
        )
        if not saturation.empty:
            saturation.to_csv(output_dir / "13_endpoint_saturation_per_chamber.csv", index=False)
            frame_window = static_endpoint_window(
                saturation, last_fraction=ENDPOINT_LAST_FRACTION,
                min_frame=int(round(OSCILLATION_START_MIN / MIN_PER_FRAME)),
            )
            logger.info(
                "Schritt 13 (statisch): Endfenster = Frames %d-%d (vor der fruehesten "
                "Saettigung; Details in 13_endpoint_saturation_per_chamber.csv).", *frame_window,
            )
    if frame_window is None:
        logger.info("Endzustand (letzte %.0f%% der Frames je Kammer) wird gebildet fuer: %s",
                    100 * ENDPOINT_LAST_FRACTION, value_cols)

    all_per_replicate, all_summary, all_trend, all_scores, all_chambers = [], [], [], [], []
    all_ctrl_trend, all_within = [], []
    for value_col in value_cols:
        per_replicate, summary = compute_endpoint_per_replicate(
            cells, value_col, last_fraction=ENDPOINT_LAST_FRACTION, frame_window=frame_window,
        )
        if summary.empty:
            continue
        # Kammer-Ebene fuer den QC-Vergleich (paarweise Kammern mit/ohne QC).
        if frame_window is not None:
            lo, hi = frame_window
            ep_rows = cells[(cells["frame"] >= lo) & (cells["frame"] <= hi)]
        else:
            fmin = cells.groupby("exp_id")["frame"].transform("min")
            fmax = cells.groupby("exp_id")["frame"].transform("max")
            ep_rows = cells[cells["frame"] > fmax - (fmax - fmin + 1) * ENDPOINT_LAST_FRACTION]
        per_chamber, _, _ = summarise_hierarchical(ep_rows, value_col)
        all_chambers.append(per_chamber)
        trend = spearman_against_period(per_replicate)
        all_per_replicate.append(per_replicate)
        all_summary.append(summary)
        if not trend.empty:
            all_trend.append(trend)

        if ctx.x_col == "osc_freq":
            # Laufen die Kontrollen der Strukturen mit der Periode? Und was
            # passiert innerhalb einer Kultur (2-3 Perioden, eine Vorkultur)?
            ctrl_trend = control_trend_check(per_replicate)
            if not ctrl_trend.empty:
                all_ctrl_trend.append(ctrl_trend)
            within = within_culture_trend(per_replicate)
            if not within.empty:
                all_within.append(within)
            # Bracket-Score: jede Periode relativ zu den Kontrollen IHRES Chips.
            score = bracket_normalise(per_replicate)
            score_trend = spearman_against_period(score) if not score.empty else pd.DataFrame()
            if not score.empty:
                all_scores.append(score)
            if not score_trend.empty:
                all_trend.append(score_trend.assign(value_col=f"{value_col}__bracket_score"))
            # Eine Datei pro osc_type: die Staemme nebeneinander statt in einem
            # Staemme x osc_type-Gitter, das bei 5 Staemmen unlesbar wird.
            for osc_type in sorted(per_replicate["osc_type"].dropna().unique()):
                sel = per_replicate["osc_type"] == osc_type
                plot_endpoint_vs_period(
                    per_replicate[sel], output_dir / f"13_endpoint_vs_period_{value_col}_{osc_type}.pdf",
                    value_col=value_col, trend=trend[trend["osc_type"] == osc_type] if not trend.empty else trend,
                    score=score[score["osc_type"] == osc_type] if not score.empty else None,
                    score_trend=score_trend[score_trend["osc_type"] == osc_type] if not score_trend.empty else None,
                    ylabel=f"{pretty_label(value_col)}\n(endpoint)",
                )
        else:
            # Kategoriale x-Achse (Medium), Facette Chip-Familie, Fehler ueber Chips.
            plot_point_errorbar(
                summary, value_col="mean", sd_col="sem",
                out_path=output_dir / f"13_endpoint_vs_{ctx.x_col}_{value_col}.pdf",
                x_col=ctx.x_col, facet_col=ctx.facet_col, color_col=PANEL_A_GROUP_COL,
                x_order=ctx.freq_order,
                ylabel=f"{pretty_label(value_col)} (endpoint)",
                title=f"{value_col}: endpoint before saturation (frames {frame_window[0]}-{frame_window[1]})\n"
                      "mean ± SEM; error unit per chip family (error_unit: chips, or chambers of one chip)"
                      if frame_window else f"{value_col}: endpoint, mean ± SEM (error unit per table)",
            )

    if not all_summary:
        logger.warning("Schritt 13: kein Endzustand berechenbar - keine Dateien geschrieben.")
        return

    pd.concat(all_per_replicate, ignore_index=True).to_csv(
        output_dir / "13_endpoint_per_chip.csv", index=False)
    if all_chambers:
        pd.concat(all_chambers, ignore_index=True).to_csv(
            output_dir / "13_endpoint_per_chamber.csv", index=False)
    if all_scores:
        pd.concat(all_scores, ignore_index=True).to_csv(
            output_dir / "13_endpoint_bracket_score.csv", index=False)
        logger.info("Tabelle gespeichert: 13_endpoint_bracket_score.csv")
    pd.concat(all_summary, ignore_index=True).to_csv(
        output_dir / "13_endpoint_summary.csv", index=False)
    logger.info("Tabellen gespeichert: 13_endpoint_per_chip.csv / 13_endpoint_summary.csv")

    if all_trend:
        trend_all = pd.concat(all_trend, ignore_index=True)
        trend_all.to_csv(output_dir / "13_endpoint_spearman.csv", index=False)
        logger.info("Tabelle gespeichert: 13_endpoint_spearman.csv\n%s",
                    trend_all.to_string(index=False))
    if all_ctrl_trend:
        ctrl_all = pd.concat(all_ctrl_trend, ignore_index=True)
        ctrl_all.to_csv(output_dir / "13_endpoint_control_trend.csv", index=False)
        logger.info("Tabelle gespeichert: 13_endpoint_control_trend.csv\n%s",
                    ctrl_all[["value_col", "biosensor", "osc_type", "rho_osc", "rho_ctrl_mean",
                              "rho_osc_minus_ctrl", "verdict"]].round(2).to_string(index=False))
    if all_within:
        pd.concat(all_within, ignore_index=True).to_csv(output_dir / "13_endpoint_within_culture.csv", index=False)
        logger.info("Tabelle gespeichert: 13_endpoint_within_culture.csv")
    else:
        logger.info(
            "Schritt 13: kein Spearman-Trendtest geschrieben - dafuer braucht es mindestens "
            "drei numerische Perioden (bei statischen Daten und beim PKO-Zweig erwartet)."
        )


def _lineage_rate_outputs(ctx: PipelineContext, per_experiment: pd.DataFrame) -> None:
    """Knospungsrate (Buds je Mutter-Stunde im Sparse-Phase-Fenster) gegen die
    Periode - dieselbe Chip-Logik wie Schritt 13: ein Chip pro Periode mit
    SEINEN Kontrollen, Kammer-Fehlerbalken, Bracket-Score, Spearman ueber
    Chips. Zeitnormiert, weil das Fenster je Kammer verschieden lang ist
    (in den echten Daten 27 bis 133 Frames)."""
    value_col = "budding_rate_per_h"
    output_dir = ctx.output_dir
    if per_experiment.empty or value_col not in per_experiment.columns:
        logger.warning("Knospungsrate: Spalte '%s' fehlt - 21_budding_rate_* uebersprungen.", value_col)
        return
    per_chip, summary = summarise_per_replicate(per_experiment, value_col)
    if summary.empty:
        logger.warning("Knospungsrate: keine Werte - 21_budding_rate_* uebersprungen.")
        return
    per_chip.to_csv(output_dir / "21_budding_rate_per_chip.csv", index=False)
    summary.to_csv(output_dir / "21_budding_rate_summary.csv", index=False)

    if ctx.x_col == "osc_freq":
        trend = spearman_against_period(per_chip)
        score = bracket_normalise(per_chip)
        score_trend = spearman_against_period(score) if not score.empty else pd.DataFrame()
        trends = [trend] if not trend.empty else []
        if not score_trend.empty:
            trends.append(score_trend.assign(value_col=f"{value_col}__bracket_score"))
        if trends:
            pd.concat(trends, ignore_index=True).to_csv(output_dir / "21_budding_rate_spearman.csv", index=False)
        if not score.empty:
            score.to_csv(output_dir / "21_budding_rate_bracket_score.csv", index=False)
        ctrl_trend = control_trend_check(per_chip)
        if not ctrl_trend.empty:
            ctrl_trend.to_csv(output_dir / "21_budding_rate_control_trend.csv", index=False)
            logger.info("Tabelle gespeichert: 21_budding_rate_control_trend.csv\n%s",
                        ctrl_trend[["biosensor", "osc_type", "rho_osc", "rho_ctrl_mean",
                                    "rho_osc_minus_ctrl", "verdict"]].round(2).to_string(index=False))
        within = within_culture_trend(per_chip)
        if not within.empty:
            within.to_csv(output_dir / "21_budding_rate_within_culture.csv", index=False)
        for osc_type in sorted(per_chip["osc_type"].dropna().unique()):
            sel = per_chip["osc_type"] == osc_type
            plot_endpoint_vs_period(
                per_chip[sel], output_dir / f"21_budding_rate_vs_period_{osc_type}.pdf",
                value_col=value_col,
                trend=trend[trend["osc_type"] == osc_type] if not trend.empty else trend,
                score=score[score["osc_type"] == osc_type] if not score.empty else None,
                score_trend=score_trend[score_trend["osc_type"] == osc_type] if not score_trend.empty else None,
                ylabel="buds per mother-hour\n(sparse-phase window)",
                title="budding rate in the sparse-phase window vs cycle period — one chip per period",
            )
    else:
        plot_point_errorbar(
            summary, value_col="mean", sd_col="sem",
            out_path=output_dir / f"21_budding_rate_vs_{ctx.x_col}.pdf",
            x_col=ctx.x_col, facet_col=ctx.facet_col, color_col=PANEL_A_GROUP_COL,
            x_order=ctx.freq_order, ylabel="buds per mother-hour (sparse-phase window)",
            title="Budding rate in the sparse-phase window, mean ± SEM "
                  "(error unit per chip family: chips, or chambers of one chip)",
        )
    logger.info("Knospungsrate gespeichert: 21_budding_rate_per_chip.csv / _summary.csv (+ Plots)")


def step_20_lineage(ctx: PipelineContext) -> None:
    """Budding-Events, Budding Ratio, Panel A, Stammbaum - alles aus dem
    Sparse-Phase-Fenster (ctx.cells_lineage), siehe relink.py."""
    cells = ctx.cells_lineage
    cells_plot = exclude_controls(cells)
    output_dir = ctx.output_dir
    freq_order = ctx.freq_order
    lineage_events = ctx.lineage_events
    mothers = ctx.mothers

    # ==================================================================
    # 20. Lineage/Budding (Ausgabe der in Schritt 3 vorab berechneten
    #     Mutter/Bud-Klassifikation, plus Budding Ratio, Panel A und
    #     Lineage-Tiefe)
    # ==================================================================
    if not lineage_events.empty:
        lineage_events.to_csv(output_dir / "20_budding_events.csv", index=False)
        logger.info("Tabelle gespeichert: 20_budding_events.csv (%d Events)", len(lineage_events))
    else:
        logger.warning("Keine Budding-Events erkannt - Budding Ratio besteht dann ausschließlich aus ruhenden Müttern.")

    # -- Budding Ratio "pro Mutter über die gesamte Beobachtungsdauer"
    #    (NICHT die Paper-Eq.-3-Zeitreihe, siehe weiter unten dafür)
    per_mother, per_experiment = compute_budding_ratio(lineage_events, cells, mothers, min_per_frame=MIN_PER_FRAME)
    if not per_mother.empty:
        per_mother.to_csv(output_dir / "21_budding_ratio_per_mother.csv", index=False)
        per_experiment.to_csv(output_dir / "21_budding_ratio_per_experiment.csv", index=False)
        logger.info("Tabellen gespeichert: 21_budding_ratio_per_mother.csv / _per_experiment.csv")
        # Die Lineage-Abbildung: Knospungsrate gegen die Periode, mit eigenen Kontrollen.
        _lineage_rate_outputs(ctx, per_experiment)

        plot_panel_a(
            cells_plot, exclude_controls(per_mother), output_dir / "21_panel_a_violin.pdf",
            group_col=ctx.panel_a_group_col, facet_col=ctx.panel_a_facet_col,
        )
    else:
        logger.warning("compute_budding_ratio() lieferte keine per_mother-Tabelle - Panel A wird übersprungen.")

    # -- Budding Ratio als ZEITREIHE nach Eq. 3 (Blöbaum et al. 2024):
    #    n_buds(t) / n_cells(t-1), analog zu Fig. 3a im Paper.
    #    HINWEIS: bei euch ist die Negativkontrolle durchgehend Starvation
    #    (im Paper: durchgehend Feast) - die Berechnung selbst ist medien-
    #    unabhängig, aber die INTERPRETATION im Vergleich zum Paper ist es
    #    nicht. Siehe budding_ratio_timeseries.py Docstring.
    budding_ts = compute_budding_ratio_timeseries(cells, lineage_events)
    budding_ts.to_csv(output_dir / "22_budding_ratio_timeseries.csv", index=False)
    logger.info("Tabelle gespeichert: 22_budding_ratio_timeseries.csv")

    budding_ts_agg = aggregate_budding_ratio_over_replicates(budding_ts)
    budding_ts_agg.to_csv(output_dir / "22_budding_ratio_timeseries_aggregated.csv", index=False)
    logger.info("Tabelle gespeichert: 22_budding_ratio_timeseries_aggregated.csv")

    plot_budding_ratio_timeseries(
        exclude_controls(budding_ts), output_dir / "22_budding_ratio_timeseries.pdf",
        group_col=ctx.panel_a_group_col, facet_col=ctx.panel_a_facet_col, freq_order=freq_order,
    )

    # -- Lineage-Baum & Generationstiefe
    if not lineage_events.empty:
        lineage_tree = build_lineage_tree(lineage_events)
        if not lineage_tree.empty:
            lineage_tree.to_csv(output_dir / "23_lineage_tree.csv", index=False)
            logger.info("Tabelle gespeichert: 23_lineage_tree.csv")

            lineage_depth = summarise_lineage_depth(lineage_tree)
            lineage_depth.to_csv(output_dir / "23_lineage_depth_summary.csv", index=False)
            logger.info(
                "Tabelle gespeichert: 23_lineage_depth_summary.csv (mittlere max. "
                "Generationstiefe=%.1f über %d Kammern)",
                lineage_depth["max_generation"].mean(), len(lineage_depth),
            )



def step_30_sensors(ctx: PipelineContext) -> None:
    """Sensor-Intensitäten und Ratios über die Zeit."""
    cells_plot = ctx.cells_plot
    output_dir = ctx.output_dir
    freq_order = ctx.freq_order
    intensity_cols = ctx.intensity_cols
    ratio_cols = ctx.ratio_cols

    # ==================================================================
    # 30. Sensor-/Ratio-Zeitverläufe (unterstützend - Grundlage für die
    #     Robustness-Auswertung in Sektion 40)
    # ==================================================================
    if not intensity_cols:
        logger.warning("Keine 'mean_<kanal>' Spalten gefunden - Sensor-Intensitäts-Plots werden übersprungen.")

    for col in intensity_cols:
        for facet_value, plot_sel, _ in _per_facet(ctx, cells_plot, None):
            plot_metric_over_time_by_frequency(
                plot_sel, col, output_dir / f"30_{col}_over_time{facet_value}.pdf", freq_order=freq_order,
                x_col=ctx.x_col, facet_col=ctx.facet_col,
            )

    if not ratio_cols:
        logger.warning("Keine 'ratio_*' Spalten gefunden - Ratio-Plots werden übersprungen.")

    for col in ratio_cols:
        for facet_value, plot_sel, _ in _per_facet(ctx, cells_plot, None):
            plot_metric_over_time_by_frequency(
                plot_sel, col, output_dir / f"31_{col}_over_time{facet_value}.pdf",
                freq_order=freq_order, ylabel=pretty_label(col), x_col=ctx.x_col, facet_col=ctx.facet_col,
            )



def step_40_robustness(ctx: PipelineContext) -> None:
    """Robustness R(t)/R(p)."""
    cells = ctx.cells
    output_dir = ctx.output_dir
    freq_order = ctx.freq_order
    ratio_cols = ctx.ratio_cols
    mu_table = ctx.mu_table
    area_table = ctx.area_table

    # ==================================================================
    # 40. Robustness-Quantifizierung R(t) und R(p) (Eq. 1, Trivellin et al.
    #     2022 / Blöbaum et al. 2024) für die in ROBUSTNESS_VALUE_COLS
    #     konfigurierten Spalten PLUS die zur Laufzeit erkannten
    #     ratio_*-Spalten (Sektion 30) - dadurch fließen die Sensor-Daten
    #     hier in die zusammenfassende Auswertung mit ein.
    #
    #     Zusätzlich (siehe unten, nach der Haupt-Schleife): µ_event und
    #     µ_area aus Sektion 10, JEWEILS NUR IN DER SINNVOLLEN VARIANTE
    #     (µ_event -> R(t) Einzelzelle, µ_area -> R(p)) - die jeweils
    #     andere Variante ist mit der Datenstruktur dieser Tabellen nicht
    #     sauber definierbar, siehe Kommentare dort.
    # ==================================================================
    # dict.fromkeys() statt set(): erhält die Reihenfolge und entfernt
    # Duplikate, falls eine Spalte versehentlich in beiden Listen auftaucht.
    robustness_value_cols = list(dict.fromkeys(ROBUSTNESS_VALUE_COLS + ratio_cols))
    logger.info("Robustness R(t)/R(p) wird berechnet für: %s", robustness_value_cols)

    for value_col in robustness_value_cols:
        if value_col not in cells.columns:
            logger.warning("Robustness: Spalte '%s' nicht in den Daten - übersprungen.", value_col)
            continue

        rt_pop = compute_rt_population(cells, value_col)
        rt_pop.to_csv(output_dir / f"40_Rt_population_{value_col}.csv", index=False)

        rt_cell = compute_rt_single_cell(cells, value_col)
        rt_cell.to_csv(output_dir / f"40_Rt_single_cell_{value_col}.csv", index=False)
        plot_rt_single_cell_distribution(
            exclude_controls(rt_cell), output_dir / f"40_Rt_single_cell_{value_col}.pdf",
            value_col="R_t_single_cell", facet_col=ctx.x_col, freq_order=freq_order,
        )

        rp = compute_rp(cells, value_col)
        rp.to_csv(output_dir / f"40_Rp_{value_col}.csv", index=False)

        rt_pop_agg = aggregate_robustness_over_replicates(rt_pop, "R_t_population")
        rt_pop_agg.to_csv(output_dir / f"40_Rt_population_{value_col}_aggregated.csv", index=False)

        rp_agg = aggregate_robustness_over_replicates(rp, "R_p")
        rp_agg.to_csv(output_dir / f"40_Rp_{value_col}_aggregated.csv", index=False)

        # Kontrollen raus aus den 'normalen' R(t)/R(p)-Plots - die laufen
        # jetzt separat über die control_consistency-Plots weiter unten.
        rt_pop_agg_plot = exclude_controls(rt_pop_agg)
        rp_agg_plot = exclude_controls(rp_agg)

        plot_point_errorbar(
            rt_pop_agg_plot, value_col="mean", out_path=output_dir / f"40_Rt_population_{value_col}.pdf",
            x_col=ctx.x_col, facet_col=ctx.facet_col, color_col=PANEL_A_GROUP_COL,
            x_order=freq_order,
            ylabel=f"R(t) — {value_col}", title=f"R(t) population level — {value_col}",
        )
        plot_point_errorbar(
            rp_agg_plot, value_col="mean", out_path=output_dir / f"40_Rp_{value_col}.pdf",
            x_col=ctx.x_col, facet_col=ctx.facet_col, color_col=PANEL_A_GROUP_COL,
            x_order=freq_order,
            ylabel=f"R(p) — {value_col}", title=f"R(p) — {value_col}",
        )
        plot_rt_vs_rp_quadrant(
            rt_pop_agg_plot, rp_agg_plot, output_dir / f"40_Rt_vs_Rp_{value_col}.pdf",
            label_col=ctx.x_col, facet_col=PANEL_A_GROUP_COL,
        )

        # Trendtest gegen die Periode - NUR fuer R(p), absichtlich nicht fuer R(t).
        #
        # R(t) ist die Varianz UEBER DIE ZEIT. Ein Frame tastet eine
        # Oszillation unter dem Nyquist-Limit bei zufaelliger Phase ab, und der
        # Alias-Beitrag ist bei der LAENGSTEN Periode am groessten (24 min bei
        # 10 min/Frame = 2.4 Abtastungen pro Zyklus) - also in derselben
        # Richtung wie der erwartete biologische Effekt. Ein R(t)-Trend gegen
        # die Periode laesst sich deshalb nicht von einem Alias-Artefakt
        # unterscheiden, und ein p-Wert dazu wuerde genau diese Verwechslung
        # nur amtlich aussehen lassen.
        #
        # R(p) ist die Varianz UEBER DIE ZELLEN innerhalb eines Frames. Alle
        # Zellen einer Kammer sehen dieselbe Medienphase, die Phase traegt zur
        # Streuung zwischen ihnen also nichts bei - R(p) ist gegen das
        # Aliasing robust und damit die Groesse, bei der ein Trendtest zulaessig
        # ist.
        rp_per_replicate, _ = summarise_per_replicate(rp, "R_p")
        rp_trend = spearman_against_period(rp_per_replicate)
        if not rp_trend.empty:
            rp_trend = rp_trend.assign(value_col=value_col)
            rp_trend.to_csv(output_dir / f"40_Rp_{value_col}_spearman.csv", index=False)
            logger.info("Tabelle gespeichert: 40_Rp_%s_spearman.csv\n%s",
                        value_col, rp_trend.to_string(index=False))

        # Kontroll-Konsistenz für R(t)/R(p), analog zur Wachstumsrate oben.
        if ctx.run_control_consistency:
            rt_control_test = test_control_consistency_across_freq(rt_pop, value_col="R_t_population")
            if not rt_control_test.empty:
                rt_control_test.to_csv(output_dir / f"40_control_consistency_Rt_population_{value_col}_kruskal.csv", index=False)
            plot_control_consistency(
                rt_pop_agg, value_col="mean", out_path=output_dir / f"40_control_consistency_Rt_population_{value_col}.pdf",
                x_order=freq_order, ylabel=f"R(t) — {value_col} (controls)",
            )

            rp_control_test = test_control_consistency_across_freq(rp, value_col="R_p")
            if not rp_control_test.empty:
                rp_control_test.to_csv(output_dir / f"40_control_consistency_Rp_{value_col}_kruskal.csv", index=False)
            plot_control_consistency(
                rp_agg, value_col="mean", out_path=output_dir / f"40_control_consistency_Rp_{value_col}.pdf",
                x_order=freq_order, ylabel=f"R(p) — {value_col} (controls)",
            )

        logger.info("Robustness R(t)/R(p) für '%s' berechnet und gespeichert.", value_col)

    # -- µ_event -> R(t) Einzelzelle: wie stabil ist die Reproduktionsrate
    #    EINER Mutter über ihre eigenen Budding-Intervalle?
    #    mu_table hat KEINE 'frame'-Spalte (nur frame_start/frame_end pro
    #    Intervall) und ist damit NICHT mit compute_rt_population()/
    #    compute_rp() kompatibel - die brauchen eine gemeinsame Zeitachse
    #    über mehrere Zellen hinweg, die es bei Event-Daten mit ihren pro
    #    Mutter unterschiedlichen Intervall-Zeitpunkten nicht sinnvoll gibt
    #    (Intervall-Enden verschiedener Mütter fallen kaum je zusammen).
    #    compute_rt_single_cell() braucht dagegen KEINE 'frame'-Spalte,
    #    nur mehrere Werte PRO ZELLE - das ist hier gegeben (mehrere
    #    Intervalle pro Mutter), daher funktioniert NUR diese Variante.
    if not mu_table.empty:
        mu_valid = mu_table[~mu_table["mu_is_artefact"]]
        rt_cell_mu = compute_rt_single_cell(mu_valid, value_col="mu", cell_id_col="mother_cell_uid")
        rt_cell_mu.to_csv(output_dir / "40_Rt_single_cell_mu_event.csv", index=False)
        plot_rt_single_cell_distribution(
            exclude_controls(rt_cell_mu), output_dir / "40_Rt_single_cell_mu_event.pdf",
            value_col="R_t_single_cell", facet_col=ctx.x_col, freq_order=freq_order,
        )
        rt_cell_mu_agg = aggregate_robustness_over_replicates(rt_cell_mu, "R_t_single_cell")
        rt_cell_mu_agg.to_csv(output_dir / "40_Rt_single_cell_mu_event_aggregated.csv", index=False)
        logger.info("R(t) Einzelzelle für µ_event berechnet und gespeichert (%d Mütter).", len(rt_cell_mu))
    else:
        logger.warning("mu_table ist leer - R(t) Einzelzelle für µ_event wird übersprungen.")

    # -- µ_area -> R(p): wie homogen ist die flächenbasierte Wachstumsrate
    #    ÜBER DIE ZELLEN einer Kammer?
    #    area_table hat pro Zelle nur EINEN Wert (der ganze Track wurde
    #    bereits zu einer Steigung verdichtet) - eine Zeitachse gibt es
    #    hier nicht mehr, R(t) ist für µ_area daher NICHT definierbar
    #    (weder Populations- noch Einzelzell-Variante: beide bräuchten
    #    mehrere Zeitpunkte pro Zelle bzw. pro Kammer).
    #    compute_rp() gruppiert technisch nach exp_id x frame - da wir
    #    keine Zeit haben, setzen wir 'frame' konstant auf 0. Dadurch
    #    entspricht jede Kammer genau einer Gruppe, und R(p) misst wie
    #    beabsichtigt die Homogenität ÜBER DIE ZELLEN einer Kammer (statt
    #    über Zeitpunkte). Nur zuverlässige Fits (fit_is_reliable) gehen
    #    ein - unzuverlässige Fits würden R(p) sonst mit Rauschen aus
    #    schlecht bestimmten Steigungen aufblähen. Mutter- und Knospen-
    #    Tracks werden NICHT getrennt (analog zu 'area'/'eccentricity'
    #    oben, die ebenfalls nicht nach cell_type filtern).
    if not area_table.empty:
        area_reliable = area_table[area_table["fit_is_reliable"]].copy()
        area_reliable["frame"] = 0
        rp_mu_area = compute_rp(area_reliable, value_col="mu_area")
        rp_mu_area.to_csv(output_dir / "40_Rp_mu_area.csv", index=False)

        rp_mu_area_agg = aggregate_robustness_over_replicates(rp_mu_area, "R_p")
        rp_mu_area_agg.to_csv(output_dir / "40_Rp_mu_area_aggregated.csv", index=False)

        plot_point_errorbar(
            exclude_controls(rp_mu_area_agg), value_col="mean", out_path=output_dir / "40_Rp_mu_area.pdf",
            x_col=ctx.x_col, facet_col=ctx.facet_col, color_col=PANEL_A_GROUP_COL,
            x_order=freq_order,
            ylabel="R(p) — µ_area", title="R(p) — µ_area (homogeneity across cells)",
        )

        # Trendtest gegen die Periode, siehe Begruendung oben (R(p), nicht R(t)).
        rp_mu_area_rep, _ = summarise_per_replicate(rp_mu_area, "R_p")
        rp_mu_area_trend = spearman_against_period(rp_mu_area_rep)
        if not rp_mu_area_trend.empty:
            rp_mu_area_trend = rp_mu_area_trend.assign(value_col="mu_area")
            rp_mu_area_trend.to_csv(output_dir / "40_Rp_mu_area_spearman.csv", index=False)
            logger.info("Tabelle gespeichert: 40_Rp_mu_area_spearman.csv\n%s",
                        rp_mu_area_trend.to_string(index=False))

        if ctx.run_control_consistency:
            rp_mu_area_control_test = test_control_consistency_across_freq(rp_mu_area, value_col="R_p")
            if not rp_mu_area_control_test.empty:
                rp_mu_area_control_test.to_csv(output_dir / "40_control_consistency_Rp_mu_area_kruskal.csv", index=False)
            plot_control_consistency(
                rp_mu_area_agg, value_col="mean", out_path=output_dir / "40_control_consistency_Rp_mu_area.pdf",
                x_order=freq_order, ylabel="R(p) — µ_area (controls)",
            )
        logger.info("R(p) für µ_area berechnet und gespeichert (%d zuverlässige Tracks).", len(area_reliable))
    else:
        logger.warning("area_table ist leer - R(p) für µ_area wird übersprungen.")



def step_50_summary(ctx: PipelineContext) -> None:
    """Zusammenfassungstabelle."""
    cells = ctx.cells
    output_dir = ctx.output_dir
    intensity_cols = ctx.intensity_cols

    # ==================================================================
    # 50. Zusammenfassungstabelle
    # ==================================================================
    summary = summary_statistics(cells, intensity_cols)
    summary.to_csv(output_dir / "50_summary_statistics.csv", index=False)
    logger.info("Tabelle gespeichert: 50_summary_statistics.csv")

    # -- Kontroll-Trend-Zusammenfassung ueber alle Readouts (Schritte 12, 13, 21):
    #    die eine Abbildung fuer den Befund, dass die Struktur die Trends traegt.
    if ctx.x_col == "osc_freq":
        parts = []
        for name in ("13_endpoint_control_trend.csv", "12_area_growth_rate_control_trend.csv",
                     "21_budding_rate_control_trend.csv"):
            path = output_dir / name
            if path.exists():
                try:
                    part = pd.read_csv(path)
                except pd.errors.EmptyDataError:
                    continue
                if not part.empty:
                    parts.append(part)
        if parts:
            all_trends = pd.concat(parts, ignore_index=True)
            all_trends.to_csv(output_dir / "50_control_trend_summary.csv", index=False)
            plot_control_trend_summary(all_trends, output_dir / "50_control_trend_summary.pdf")
            counts = all_trends["verdict"].astype(str).str.split(":").str[0].value_counts().to_dict()
            logger.info("Kontroll-Trend-Zusammenfassung gespeichert: 50_control_trend_summary.csv/.pdf - %s", counts)



def step_90_appendix(ctx: PipelineContext) -> None:
    """Anhang: Morphologie-Scatter und Trajektorien."""
    cells = ctx.cells
    output_dir = ctx.output_dir
    lineage_events = ctx.lineage_events
    mothers = ctx.mothers

    # ==================================================================
    # 90. ANHANG: Morphologie-Scatter, Einzelzell-Trajektorien, stabile
    #     Mutter-Trajektorien - deskriptive Plausibilitäts-Checks ohne
    #     eigene quantitative Kennzahl, siehe Diskussion zur Restrukturierung.
    # ==================================================================
    plot_morphology_scatter(cells, output_dir / "90_morphology_scatter.pdf")

    # Dieser Plot zeigt 'area', hing aber an intensity_cols - dadurch fehlte er
    # bei den statischen Daten komplett (Wildtyp, keine Fluoreszenzkanäle).
    if "area" in cells.columns:
        plot_single_cell_trajectories(cells, "area", output_dir / "91_single_cell_trajectories.pdf")
    else:
        logger.warning("Spalte 'area' fehlt - 91_single_cell_trajectories.pdf übersprungen.")

    # Stabil getrackte Mütter im Detail: pro Bedingungs-Kombination
    # (STABLE_MOTHER_GROUP_COLS) EINE stabile Mutter, mit Budding-Markern
    # und nur ihrem eigenen Ratio-Kanal (siehe mother_trajectories.py) -
    # ergänzt die aggregierten Auswertungen oben (Panel A, Fig. 3a) um
    # eine Einzelzell-Ebene zur visuellen Plausibilitätsprüfung.
    stable_base_cols = [c for c in STABLE_MOTHER_BASE_VALUE_COLS if c in cells.columns]
    if not mothers.empty and stable_base_cols:
        plot_stable_mother_per_group(
            cells, lineage_events, mothers,
            out_path=output_dir / "92_stable_mother_trajectories.pdf",
            group_cols=STABLE_MOTHER_GROUP_COLS, base_value_cols=stable_base_cols,
            min_coverage=STABLE_MOTHER_MIN_COVERAGE, min_per_frame=MIN_PER_FRAME,
        )
        logger.info("Plot gespeichert: 92_stable_mother_trajectories.pdf")
    else:
        logger.warning("Keine stabilen Mütter oder keine base_value_cols - 92_stable_mother_trajectories.pdf übersprungen.")



def step_95_sensor_controls(ctx: PipelineContext) -> None:
    """Anhang: Sensor-Controls (PosCtrl vs. NegCtrl)."""
    cells = ctx.cells
    output_dir = ctx.output_dir
    run_sensor_controls = ctx.run_sensor_controls

    # ==================================================================
    # 95. ANHANG: Sensor-Controls (PosCtrl-vs-NegCtrl-Validierung pro
    #     Biosensor). Reine Kalibrierungs-/QC-Auswertung des Sensors selbst,
    #     nicht Teil der biologischen Kernaussage - siehe Diskussion zur
    #     Restrukturierung. Zellwerte werden vorab pro Replikat summiert,
    #     um Pseudoreplikation zu vermeiden.
    #
    #     Wird für die statischen Daten übersprungen (run_sensor_controls=
    #     False): dort wird nur der Wildtyp kultiviert, es gibt also keine
    #     Fluoreszenzkanäle/Sensor-Ratios auszuwerten.
    # ==================================================================
    if not run_sensor_controls:
        logger.info("Sensor-Control-Validierung (Schritt 95) übersprungen (run_sensor_controls=False).")
        SENSOR_CONFIG_ITEMS = {}
    else:
        SENSOR_CONFIG_ITEMS = SENSOR_CONFIG

    for biosensor, sensor_cfg in SENSOR_CONFIG_ITEMS.items():
        ratio_col = sensor_cfg.ratio_name
        if ratio_col not in cells.columns or not cells[ratio_col].notna().any():
            logger.info("Sensor-control validation: no usable '%s' values for %s; skipped.", ratio_col, biosensor)
            continue

        sensor_data = cells[cells["biosensor"] == biosensor]
        for osc_type in sorted(sensor_data["osc_type"].dropna().unique()):
            analysis_start_min = OSCILLATION_START_MIN
            controls = prepare_sensor_controls(
                cells, value_col=ratio_col, biosensor=biosensor, osc_type=osc_type,
            )
            control_time, control_summary = summarise_sensor_controls(
                controls, value_col=ratio_col, analysis_start_min=analysis_start_min,
            )
            if control_time.empty:
                continue

            output_stem = f"95_{biosensor}_{osc_type}_{ratio_col}"
            control_time.to_csv(output_dir / f"{output_stem}_timeseries_per_replicate.csv", index=False)
            control_summary.to_csv(output_dir / f"{output_stem}_summary_per_replicate.csv", index=False)
            sensor_label = ratio_col.replace("ratio_", "")
            plot_sensor_control_timeseries(
                control_time, output_dir / f"{output_stem}_timeseries.pdf",
                sensor_label=sensor_label,
                preconditioning_end_min=OSCILLATION_START_MIN,
            )
            plot_sensor_control_comparison(
                control_summary, output_dir / f"{output_stem}_comparison.pdf",
                sensor_label=sensor_label, analysis_start_min=analysis_start_min,
                control_labels=CONTROL_CONCENTRATION_LABELS.get(osc_type),
            )

            # A flat ratio can arise because both raw channels shift together
            # or because neither channel contains a detectable control effect.
            # The two raw-channel panels make this distinction inspectable.
            channel_a_controls = prepare_sensor_controls(
                cells, value_col=sensor_cfg.channel_a, biosensor=biosensor, osc_type=osc_type,
            )
            channel_b_controls = prepare_sensor_controls(
                cells, value_col=sensor_cfg.channel_b, biosensor=biosensor, osc_type=osc_type,
            )
            channel_a_time, _ = summarise_sensor_controls(
                channel_a_controls, value_col=sensor_cfg.channel_a, analysis_start_min=analysis_start_min,
            )
            channel_b_time, _ = summarise_sensor_controls(
                channel_b_controls, value_col=sensor_cfg.channel_b, analysis_start_min=analysis_start_min,
            )
            plot_sensor_raw_channel_timeseries(
                channel_a_time, channel_b_time, output_dir / f"{output_stem}_raw_channels.pdf",
                sensor_label=sensor_label, channel_a_label=sensor_cfg.channel_a,
                channel_b_label=sensor_cfg.channel_b,
                preconditioning_end_min=OSCILLATION_START_MIN,
            )

            chamber_summary = summarise_sensor_controls_by_chamber(
                controls, value_col=ratio_col, analysis_start_min=analysis_start_min,
            )
            if not chamber_summary.empty:
                chamber_summary.to_csv(output_dir / f"{output_stem}_control_chambers.csv", index=False)
                plot_control_chamber_comparison(
                    chamber_summary, output_dir / f"{output_stem}_control_chambers.pdf",
                    sensor_label=sensor_label,
                )

    logger.info("=== Fertig! Alle Plots & Tabellen in: %s ===", output_dir)




@dataclass(frozen=True)
class Step:
    """Ein benannter Auswertungsschritt."""

    key: str
    title: str
    run: Callable[[PipelineContext], None]


STEPS: list[Step] = [
    Step("00_overview", "Übersicht / Sanity-Check", step_00_overview),
    Step("10_growth", "Zellfläche, µ_event und µ_area", step_10_growth),
    Step("13_endpoint", "Kumulativer Endzustand gegen die Periode (+ Spearman)", step_13_endpoint),
    Step("20_lineage", "Budding-Events, Budding Ratio, Panel A, Stammbaum", step_20_lineage),
    Step("30_sensors", "Sensor-Intensitäten und Ratios über die Zeit", step_30_sensors),
    Step("40_robustness", "Robustness R(t)/R(p)", step_40_robustness),
    Step("50_summary", "Zusammenfassungstabelle", step_50_summary),
    Step("90_appendix", "Anhang: Morphologie-Scatter und Trajektorien", step_90_appendix),
    Step("95_sensor_controls", "Anhang: Sensor-Controls (PosCtrl vs. NegCtrl)", step_95_sensor_controls),
]
