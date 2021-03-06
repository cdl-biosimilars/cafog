import logging
import math
import sys
from typing import Optional

import pandas as pd

from PyQt5.QtChart import (QBarCategoryAxis, QBarSeries, QBarSet,
                           QChart, QChartView, QValueAxis)
from PyQt5.QtCore import (Qt, QLibraryInfo, QLocale,
                          QMargins, QRectF, QSize, QTranslator)
from PyQt5.QtGui import QBrush, QColor, QDropEvent, QMouseEvent, QPainter
from PyQt5.QtSvg import QSvgGenerator
from PyQt5.QtWidgets import (QApplication, QHeaderView, QMainWindow,
                             QMessageBox, QTableWidgetItem, QTextEdit, QWidget)

from correction import GlycationGraph, read_clean_datasets, read_library

from main_window import Ui_MainWindow
from widgets import (FileTypes, SortableTableWidgetItem,
                     get_filename, open_manual)


class TextEditHandler(logging.Handler):
    """
    A handler for Python's logging module
    which redirects logging output to a QTextEdit.

    .. automethod:: __init__
    """

    def __init__(self,
                 widget: QTextEdit) -> None:
        """
        Initialize the handler.

        :param QTextEdit widget: textedit that sould display the log
        :return: nothing
        :rtype: None
        """

        super().__init__()
        self.widget = widget

    def emit(self,
             record: logging.LogRecord) -> None:
        """
        Send the log record to the textedit.

        :param logging.LogRecord record: record to be displayed
        :return: nothing
        :rtype: None
        """

        self.widget.append(self.format(record))


class MainWindow(QMainWindow, Ui_MainWindow):
    """
    The main window.

    .. automethod:: __init__
    """

    def __init__(self,
                 parent: QWidget=None) -> None:
        """
        Initialize the main window.

        :param QWidget parent: parent widget
        :return: nothing
        :rtype: None
        """

        # initialize the GUI
        super().__init__(parent)
        self.setupUi(self)

        # instance attributes
        self.glycation = None
        self.glycation_graph = None
        self.glycoforms = None
        self.glycoforms_agg = None
        self.last_path = None
        self.library = None
        self.results = None
        self.results_agg = None

        # actions
        self.cbAggGlycoforms.clicked.connect(self.toggle_agg_glycoforms)
        self.cbAggResults.clicked.connect(self.toggle_agg_results)

        self.cvGlycation.dragEnterEvent = lambda e: e.accept()
        self.cvGlycation.dragMoveEvent = lambda e: e.accept()
        self.cvGlycation.dropEvent = lambda e: self.drop_file(
            e, self.cvGlycation)
        self.cvGlycoforms.dragEnterEvent = lambda e: e.accept()
        self.cvGlycoforms.dragMoveEvent = lambda e: e.accept()
        self.cvGlycoforms.dropEvent = lambda e: self.drop_file(
            e, self.cvGlycoforms)

        self.btCorrect.clicked.connect(self.correct_abundances)
        self.btHelp.clicked.connect(open_manual)
        self.btLoadGlycation.clicked.connect(lambda: self.load_glycation())
        self.btLoadGlycoforms.clicked.connect(lambda: self.load_glycoforms())
        self.btLoadLibrary.clicked.connect(lambda: self.load_library())
        self.btQuit.clicked.connect(QApplication.instance().quit)
        self.btSampleData.clicked.connect(self.load_sample_data)
        self.btSaveGraph.clicked.connect(self.save_graph)
        self.btSaveResults.clicked.connect(self.save_results)

        self.sbAggGlycoforms.valueChanged.connect(self.agg_glycoforms)
        self.sbAggResults.valueChanged.connect(self.agg_results)

        self.twLibrary.dragEnterEvent = lambda e: e.accept()
        self.twLibrary.dragMoveEvent = lambda e: e.accept()
        self.twLibrary.dropEvent = lambda e: self.drop_file(
            e, self.twLibrary)

        # GUI modifications
        self.cvGlycation.chart().setBackgroundRoundness(0)
        self.cvGlycation.chart().layout().setContentsMargins(0, 0, 0, 0)
        drop_hint = self.cvGlycation.scene().addText(
            self.tr("Drag and drop glycation data\nor click 'Load ...'"))
        drop_hint.setDefaultTextColor(QColor("#888888"))

        self.cvGlycoforms.chart().setBackgroundRoundness(0)
        self.cvGlycoforms.chart().layout().setContentsMargins(0, 0, 0, 0)
        self.cvGlycoforms.setRubberBand(QChartView.HorizontalRubberBand)
        self.old_gf_mouse_release_event = self.cvGlycoforms.mouseReleaseEvent
        self.cvGlycoforms.mouseReleaseEvent = self.zoom_glycoform_graph
        drop_hint = self.cvGlycoforms.scene().addText(
            self.tr("Drag and drop glycoform data\nor click 'Load ...'"))
        drop_hint.setDefaultTextColor(QColor("#888888"))

        self.cvResults.chart().setBackgroundRoundness(0)
        self.cvResults.chart().layout().setContentsMargins(0, 0, 0, 0)
        self.cvResults.setRubberBand(QChartView.HorizontalRubberBand)
        self.old_re_mouse_release_event = self.cvResults.mouseReleaseEvent
        self.cvResults.mouseReleaseEvent = self.zoom_results_graph

        self.lbGlycation.setText("")

        self.lbGlycoform.setText("")

        self.lbResults.setText("")

        self.twLibrary.horizontalHeader().setVisible(False)
        self.twLibrary.verticalHeader().setVisible(False)
        self.twLibrary.verticalHeader().setSectionResizeMode(
            QHeaderView.Fixed)
        self.twLibrary.verticalHeader().setDefaultSectionSize(22)
        self.twLibrary.setRowCount(2)
        self.twLibrary.setSpan(0, 0, 2, 2)
        item = QTableWidgetItem(
            self.tr("""Drag and drop glycan library data\n
or click 'Load ...' (optional)"""))
        item.setFlags(Qt.ItemIsEnabled)
        item.setForeground(QBrush(QColor("#888888")))
        self.twLibrary.setItem(0, 0, item)

        self.twResults.horizontalHeader().setSectionResizeMode(
            0, QHeaderView.Stretch)
        self.twResults.verticalHeader().setDefaultSectionSize(22)

        # logger
        handler = TextEditHandler(self.teLog)
        handler.setFormatter(logging.Formatter("[%(levelname)s]  %(message)s"))
        logging.getLogger().addHandler(handler)
        logging.getLogger().setLevel(logging.INFO)

    def load_sample_data(self):
        """
        Load sample input data.

        :return: nothing
        :rtype: None
        """

        self.load_glycoforms("sample_data/glycoforms.csv")
        self.load_glycation("sample_data/glycation.csv")
        self.load_library("sample_data/glycan_library.csv")

    def load_glycation(self,
                       filename: Optional[str]=None) -> None:
        """
        Load glycation data from a CSV file and display it
        in the corresponding chart view.

        :param str filename: directly load this file
        :return: nothing, sets self.se_glycation
        :rtype: None
        """

        # load and clean glycation data
        if filename is None:
            filename, self.last_path = get_filename(
                self, "open", self.tr("Load glycation data ..."),
                self.last_path, FileTypes(["csv"]))
            if filename is None:
                return

        logging.info(self.tr("Loading glycation data in '{}'")
                     .format(filename))
        try:
            self.glycation = read_clean_datasets(filename)
        except (OSError, ValueError) as e:
            logging.error(str(e))
            QMessageBox.critical(self, self.tr("Error"), str(e))
            return

        # extract x- and y-values from series
        x_values = [str(i) for i in self.glycation.index]
        y_values = [a.nominal_value for a in self.glycation]

        # assemble the chart
        bar_set = QBarSet("glycation abundance")
        bar_set.append(y_values)
        bar_set.setColor(QColor("#a1dab4"))
        bar_set.hovered.connect(self.update_glycation_label)
        bar_series = QBarSeries()
        bar_series.append(bar_set)

        x_axis = QBarCategoryAxis()
        x_axis.append(x_values)
        x_axis.setTitleText(self.tr("count"))

        y_axis = QValueAxis()
        y_axis.setRange(0, 100)
        y_axis.setTitleText(self.tr("abundance"))
        y_axis.setLabelFormat("%d")

        chart = QChart()
        chart.addSeries(bar_series)
        chart.setAxisX(x_axis, bar_series)
        chart.setAxisY(y_axis, bar_series)
        chart.legend().setVisible(False)
        chart.setBackgroundRoundness(0)
        chart.layout().setContentsMargins(0, 0, 0, 0)
        chart.setMargins(QMargins(5, 5, 5, 5))
        self.cvGlycation.setChart(chart)

    def drop_file(self,
                  e: QDropEvent,
                  source: QWidget) -> None:
        """
        Open a file dropped on any of the inpuit data widgets.

        :param QDropEvent e: drop event which prompted calling this function
        :param QWidget source: widget that received the drop
        :return: nothing
        :rtype: None
        """

        if e.mimeData().hasUrls():
            filename = e.mimeData().urls()[0].toLocalFile()
            if source == self.cvGlycoforms:
                self.load_glycoforms(filename)
            elif source == self.cvGlycation:
                self.load_glycation(filename)
            elif source == self.twLibrary:
                self.load_library(filename)

    def update_glycation_label(self,
                               hover: bool,
                               bar_index: int) -> None:
        """
        Display information on glycation abundance
        of the bar under the cursor.

        :param bool hover: True if the mouse hovers over a bar; false otherwise
        :param int bar_index: index of the bar that triggered the signal
        :return: nothing
        :rtype: None
        """

        if hover:
            if bar_index == 0:
                label = self.tr("no glycation")
            elif bar_index == 1:
                label = self.tr("1 glycation")
            else:
                label = self.tr("{} glycations").format(bar_index)
            self.lbGlycation.setText(
                "{}: <b>{:.2f}</b> ± {:.2f} %".format(
                    label,
                    self.glycation[bar_index].nominal_value,
                    self.glycation[bar_index].std_dev))
        else:
            self.lbGlycation.setText("")

    def load_glycoforms(self,
                        filename: Optional[str]=None) -> None:
        """
        Load glycoform data from a CSV file.

        :param str filename: directly load this file
        :return: nothing, sets self.se_glycoforms
        :rtype: None
        """

        # load and clean glycoform data
        if filename is None:
            filename, self.last_path = get_filename(
                self, "open", self.tr("Load glycoform data ..."),
                self.last_path, FileTypes(["csv"]))
            if filename is None:
                return

        logging.info(self.tr("Loading glycoform data in '{}'")
                     .format(filename))
        try:
            self.glycoforms = (read_clean_datasets(filename)
                               .sort_values(ascending=False))
        except (OSError, ValueError) as e:
            logging.error(str(e))
            QMessageBox.critical(self, self.tr("Error"), str(e))
            return

        for widget in (self.cbAggGlycoforms,
                       self.sbAggGlycoforms,
                       self.lbAggGlycoforms):
            widget.setEnabled(True)
        self.sbAggGlycoforms.setMaximum(len(self.glycoforms) - 2)
        self.agg_glycoforms()

    def update_glycoform_label(self,
                               hover: bool,
                               bar_index: int) -> None:
        """
        Display information on glycoform abundance
        of the bar under the cursor.

        :param bool hover: True if the mouse hovers over a bar; false otherwise
        :param int bar_index: index of the bar that triggered the signal
        :return: nothing
        :rtype: None
        """

        if hover:
            self.lbGlycoform.setText(
                "{}: <b>{:.2f}</b> ± {:.2f} %".format(
                    self.glycoforms_agg.index[bar_index],
                    self.glycoforms_agg[bar_index].nominal_value,
                    self.glycoforms_agg[bar_index].std_dev))
        else:
            self.lbGlycoform.setText("")

    def zoom_glycoform_graph(self,
                             e: QMouseEvent) -> None:
        """
        Correctly reset zoom of glycoform bar chart on richt click.

        :param QMouseEvent e: mouse event that triggered the signal
        :return: nothing
        :rtype: None
        """

        if e.button() == Qt.RightButton:
            x_axis = self.cvGlycoforms.chart().axisX()
            if x_axis is not None:
                x_axis.setRange(x_axis.categories()[0],
                                x_axis.categories()[-1])
            return
        self.old_gf_mouse_release_event(e)

    def agg_glycoforms(self) -> None:
        """
        Display glycoform data in the corresponding chart view.

        :return: nothing
        :rtype: None
        """

        # aggregate "other" abundances
        if self.cbAggGlycoforms.isChecked():
            agg_abundance = (self.glycoforms
                             .iloc[self.sbAggGlycoforms.value():]
                             .sum())
            self.glycoforms_agg = (
                self.glycoforms
                    .iloc[:self.sbAggGlycoforms.value()]
                    .append(
                        pd.Series(agg_abundance, index=[self.tr("other")])))
        else:
            self.glycoforms_agg = self.glycoforms

        # extract x- and y-values from series
        x_values = [str(i) for i in self.glycoforms_agg.index]
        y_values = [a.nominal_value for a in self.glycoforms_agg]

        # assemble the chart
        bar_set = QBarSet("glycoform abundance")
        bar_set.append(y_values)
        bar_set.setColor(QColor("#2c7fb8"))
        bar_set.hovered.connect(self.update_glycoform_label)
        bar_series = QBarSeries()
        bar_series.append(bar_set)

        x_axis = QBarCategoryAxis()
        x_axis.append(x_values)
        x_axis.setTitleVisible(False)
        x_axis.setLabelsAngle(270)

        range_max = max(self.glycoforms_agg).nominal_value
        range_max = math.ceil(range_max / 20) * 20
        tick_count = range_max // 20 + 1
        y_axis = QValueAxis()
        y_axis.setRange(0, range_max)
        y_axis.setTickCount(tick_count)
        y_axis.setTitleText(self.tr("abundance"))
        y_axis.setLabelFormat("%d")

        chart = QChart()
        chart.addSeries(bar_series)
        chart.setAxisX(x_axis, bar_series)
        chart.setAxisY(y_axis, bar_series)
        chart.legend().setVisible(False)
        chart.setBackgroundRoundness(0)
        chart.layout().setContentsMargins(0, 0, 0, 0)
        chart.setMargins(QMargins(5, 5, 5, 5))
        self.cvGlycoforms.setChart(chart)

    def toggle_agg_glycoforms(self) -> None:
        """
        Enable or disable the agg glycoforms spinbox.

        :return: nothing
        :rtype: None
        """

        if self.cbAggGlycoforms.isChecked():
            self.sbAggGlycoforms.setEnabled(True)
        else:
            self.sbAggGlycoforms.setEnabled(False)
        self.agg_glycoforms()

    def load_library(self,
                     filename: Optional[str]=None) -> None:
        """
        Load a glycan library and display in the respective table.

        :param str filename: directly load this file
        :return: nothing, changes self.se_library
        :rtype: None
        """

        if filename is None:
            filename, self.last_path = get_filename(
                self, "open", self.tr("Load glycan library ..."),
                self.last_path, FileTypes(["csv"]))
            if filename is None:
                return

        logging.info(self.tr("Loading glycan library in '{}'")
                     .format(filename))
        try:
            self.library = read_library(filename)
        except (OSError, ValueError) as e:
            logging.error(str(e))
            QMessageBox.critical(self, self.tr("Error"), str(e))
            return

        # fill the table
        self.twLibrary.clearContents()
        self.twLibrary.setRowCount(0)
        self.twLibrary.horizontalHeader().setVisible(True)
        self.twLibrary.verticalHeader().setVisible(True)
        for row_id, row in self.library.fillna("").iterrows():
            self.twLibrary.insertRow(row_id)
            for col_id in (0, 1):
                item = QTableWidgetItem(str(row.iloc[col_id]))
                item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
                self.twLibrary.setItem(row_id, col_id, item)

    def correct_abundances(self) -> None:
        """
        Calculate corrected abundances.

        :return: nothing, sets self.results
        :rtype: None
        """

        missing_input = [self.tr("Please provide the following input data:")]
        if self.glycoforms is None:
            missing_input.append(self.tr("glycoforms"))
        if self.glycation is None:
            missing_input.append(self.tr("glycation"))
        if len(missing_input) > 1:
            logging.error(" ".join(missing_input))
            QMessageBox.critical(
                self, self.tr("Error"), "\n- ".join(missing_input))
            return

        logging.info(self.tr("Correcting dataset  ..."))
        QApplication.setOverrideCursor(Qt.WaitCursor)
        QApplication.processEvents()
        try:
            self.glycation_graph = GlycationGraph(glycan_library=self.library,
                                                  glycoforms=self.glycoforms,
                                                  glycation=self.glycation)
            self.glycation_graph.correct_abundances()
            self.results = self.glycation_graph.to_dataframe()
        except ValueError as e:
            QApplication.restoreOverrideCursor()
            logging.error(str(e))
            QMessageBox.critical(self, self.tr("Error"), str(e))
            return

        logging.info(self.tr("... done!"))
        self.show_results()
        QApplication.restoreOverrideCursor()

    def show_results(self) -> None:
        """
        Show correction results in the table
        and call the chart creating method.

        :return: nothing
        :rtype: None
        """

        # fill the table
        self.twResults.clearContents()
        for row_id, row in self.results.iterrows():
            self.twResults.insertRow(row_id)

            # glycoform
            item = SortableTableWidgetItem(row.iloc[0].split(" or ", 1)[0])
            item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
            self.twResults.setItem(row_id, 0, item)

            # abundances and errors
            for col_id in range(1, 5):
                item = SortableTableWidgetItem(
                    "{:.2f}".format(row.iloc[col_id]))
                item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
                self.twResults.setItem(row_id, col_id, item)

            # change of abundance
            item = SortableTableWidgetItem(
                "{:.2f}".format(row.iloc[3] - row.iloc[1]))
            item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
            self.twResults.setItem(row_id, 5, item)

        # create chart
        for widget in (self.cbAggResults,
                       self.sbAggResults,
                       self.lbAggResults,
                       self.btSaveResults,
                       self.btSaveGraph):
            widget.setEnabled(True)
        self.sbAggResults.setMaximum(len(self.results) - 2)
        self.agg_results()
        self.update_results_label(False, 0)

    def update_results_label(self,
                             hover: bool,
                             bar_index: int) -> None:
        """
        Display information on abundance of the bars under the cursor.

        :param bool hover: True if the mouse hovers over a bar; false otherwise
        :param int bar_index: index of the bar that triggered the signal
        :return: nothing
        :rtype: None
        """

        if hover:
            self.lbResults.setText(
                "{}: "
                "{} <b><font color='#225ea8'>{:.2f}</font></b> "
                "± {:.2f} %, "
                "{} <b><font color='#41b6c4'>{:.2f}</font></b> "
                "± {:.2f} %".format(
                    self.results_agg.iloc[bar_index, 0]
                        .split(" or ", 1)[0],
                    self.tr("observed"),
                    self.results_agg.iloc[bar_index, 1],
                    self.results_agg.iloc[bar_index, 2],
                    self.tr("corrected"),
                    self.results_agg.iloc[bar_index, 3],
                    self.results_agg.iloc[bar_index, 4]))
        else:
            self.lbResults.setText(
                "<font color='#225ea8'>&#x25A0;</font> {}"
                "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"
                "<font color='#41b6c4'>&#x25A0;</font> {}"
                .format(self.tr("observed"), self.tr("corrected")))

    def zoom_results_graph(self,
                           e: QMouseEvent) -> None:
        """
        Correctly reset zoom of results bar chart on richt click.

        :param QMouseEvent e: mouse event that triggered the signal
        :return: nothing
        :rtype: None
        """

        if e.button() == Qt.RightButton:
            x_axis = self.cvResults.chart().axisX()
            if x_axis is not None:
                x_axis.setRange(x_axis.categories()[0],
                                x_axis.categories()[-1])
            return
        self.old_re_mouse_release_event(e)

    def agg_results(self) -> None:
        """
        Display results in the corresponding chart view.

        :return: nothing
        :rtype: None
        """

        # aggregate "other" abundances
        if self.cbAggResults.isChecked():
            agg_abundance = (self.results
                             .iloc[self.sbAggResults.value():]
                             .sum())
            agg_abundance["glycoform"] = self.tr("other")
            self.results_agg = (
                self.results
                    .iloc[:self.sbAggResults.value()]
                    .append(agg_abundance, ignore_index=True))
        else:
            self.results_agg = self.results

        # extract x- and y-values from series
        x_values = list(self.results_agg["glycoform"].str.split(" or").str[0])
        y_values_obs = list(self.results_agg["abundance"])
        y_values_cor = list(self.results_agg["corr_abundance"])

        # assemble the chart
        bar_set_obs = QBarSet(self.tr("observed"))
        bar_set_obs.append(y_values_obs)
        bar_set_obs.setColor(QColor("#225ea8"))
        bar_set_obs.hovered.connect(self.update_results_label)
        bar_set_cor = QBarSet(self.tr("corrected"))
        bar_set_cor.append(y_values_cor)
        bar_set_cor.setColor(QColor("#41b6c4"))
        bar_set_cor.hovered.connect(self.update_results_label)
        bar_series = QBarSeries()
        bar_series.append([bar_set_obs, bar_set_cor])

        x_axis = QBarCategoryAxis()
        x_axis.append(x_values)
        x_axis.setTitleVisible(False)
        x_axis.setLabelsAngle(270)

        range_min = min(self.results_agg[["abundance", "corr_abundance"]]
                            .min().min(),
                        0)
        range_min = math.floor(range_min / 20) * 20
        range_max = (self.results_agg[["abundance", "corr_abundance"]]
                         .max().max())
        range_max = math.ceil(range_max / 20) * 20
        tick_count = (range_max - range_min) // 20 + 1
        y_axis = QValueAxis()
        y_axis.setRange(range_min, range_max)
        y_axis.setTickCount(tick_count)
        y_axis.setTitleText(self.tr("abundance"))
        y_axis.setLabelFormat("%d")

        chart = QChart()
        chart.addSeries(bar_series)
        chart.setAxisX(x_axis, bar_series)
        chart.setAxisY(y_axis, bar_series)
        chart.legend().setVisible(False)
        chart.setBackgroundRoundness(0)
        chart.layout().setContentsMargins(0, 0, 0, 0)
        chart.setMargins(QMargins(5, 5, 5, 5))
        self.cvResults.setChart(chart)

    def toggle_agg_results(self) -> None:
        """
        Enable or disable the agg results spinbox.

        :return: nothing
        :rtype: None
        """

        if self.cbAggResults.isChecked():
            self.sbAggResults.setEnabled(True)
        else:
            self.sbAggResults.setEnabled(False)
        self.agg_results()

    def save_results(self) -> None:
        """
        Save results chart or table.

        :return: nothing
        :rtype: None
        """

        filename, self.last_path = get_filename(
            self, "save", self.tr("Save results ..."),
            self.last_path, FileTypes(["csv", "png", "svg"]))
        if filename is None:
            return

        logging.info(self.tr("Saving results to '{}'").format(filename))
        try:
            if filename.endswith("csv"):
                self.results.to_csv(filename, index=False)
            elif filename.endswith("png"):
                self.cvResults.grab().save(filename)
            elif filename.endswith("svg"):
                output_rect = QRectF(
                    self.cvResults.chart().scene().sceneRect())
                output_size = QSize(output_rect.size().toSize())
                svg_generator = QSvgGenerator()
                svg_generator.setFileName(filename)
                svg_generator.setSize(output_size)
                svg_generator.setViewBox(output_rect)
                painter = QPainter()
                painter.begin(svg_generator)
                self.cvResults.chart().scene().render(
                    painter,
                    source=output_rect,
                    target=output_rect,
                    mode=Qt.IgnoreAspectRatio)
                painter.end()
        except (OSError, ValueError) as e:
            logging.error(str(e))
            QMessageBox.critical(self, self.tr("Error"), str(e))
            return

    def save_graph(self) -> None:
        """
        Save the glycation graph.

        :return: nothing
        :rtype: None
        """

        filename, self.last_path = get_filename(
            self, "save", self.tr("Save glycation graph ..."),
            self.last_path, FileTypes(["gv", "gexf"]))
        if filename is None:
            return

        logging.info(self.tr("Saving glycation graph to '{}'")
                     .format(filename))
        try:
            if filename.endswith("gv"):
                self.glycation_graph.to_dot(filename)
            elif filename.endswith("gexf"):
                self.glycation_graph.to_gexf(filename)
        except (OSError, ValueError) as e:
            logging.error(str(e))
            QMessageBox.critical(self, self.tr("Error"), str(e))
            return


def _main() -> None:
    """
    Execute the main application loop.

    :return: nothing
    :rtype: None
    """

    app = QApplication(sys.argv)
    QLocale.setDefault(QLocale.c())

    # translation
    language = QLocale.system().name()[:2]
    translations_path = QLibraryInfo.location(QLibraryInfo.TranslationsPath)

    base_translator = QTranslator()
    base_translator.load("qtbase_{}".format(language), translations_path)
    app.installTranslator(base_translator)

    custom_translator = QTranslator()
    custom_translator.load("cafog_{}".format(language))
    app.installTranslator(custom_translator)

    # generate main window
    frame = MainWindow()
    frame.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    _main()
