import logging
import sys
from typing import Optional

import pandas as pd

from PyQt5.QtChart import (QBarCategoryAxis, QBarSeries, QBarSet,
                           QChart, QChartView, QValueAxis)
from PyQt5.QtCore import Qt, QLibraryInfo, QLocale, QMargins, QTranslator
from PyQt5.QtGui import QMouseEvent
from PyQt5.QtWidgets import (QApplication, QHeaderView, QMainWindow,
                             QMessageBox, QTableWidgetItem, QTextEdit,
                             QWhatsThis, QWidget)

from correction import GlycationGraph, read_clean_datasets, read_library

from main_window import Ui_MainWindow
from widgets import FileTypes, SortableTableWidgetItem, get_filename


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
        self.se_glycation = None
        self.se_glycoforms = None
        self.se_glycoforms_agg = None
        self.se_library = None
        self.last_path = None
        self.results = None

        # actions
        self.cbAggGlycoforms.clicked.connect(self.toggle_agg_glycoforms)

        self.btCorrect.clicked.connect(self.correct_abundances)
        self.btHelp.clicked.connect(lambda: QWhatsThis.enterWhatsThisMode())
        self.btLoadGlycation.clicked.connect(lambda: self.load_glycation())
        self.btLoadGlycoforms.clicked.connect(lambda: self.load_glycoforms())
        self.btLoadLibrary.clicked.connect(lambda: self.load_library())
        self.btQuit.clicked.connect(QApplication.instance().quit)
        self.btSampleData.clicked.connect(self.load_sample_data)

        self.sbAggGlycoforms.valueChanged.connect(self.agg_glycoforms)

        # GUI modifications
        self.cvGlycation.chart().setBackgroundRoundness(0)
        self.cvGlycation.chart().layout().setContentsMargins(0, 0, 0, 0)

        self.cvGlycoforms.chart().setBackgroundRoundness(0)
        self.cvGlycoforms.chart().layout().setContentsMargins(0, 0, 0, 0)
        self.cvGlycoforms.setRubberBand(QChartView.HorizontalRubberBand)
        self.old_gf_mouse_release_event = self.cvGlycoforms.mouseReleaseEvent
        self.cvGlycoforms.mouseReleaseEvent = self.zoom_glycoform_graph

        self.cvResults.chart().setBackgroundRoundness(0)
        self.cvResults.chart().layout().setContentsMargins(0, 0, 0, 0)

        self.lbGlycation.setText("")

        self.lbGlycoform.setText("")

        self.twLibrary.verticalHeader().setSectionResizeMode(
            QHeaderView.Fixed)
        self.twLibrary.verticalHeader().setDefaultSectionSize(22)
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

        self.load_glycoforms("data/0_day5_glycoforms.csv")
        self.load_glycation("data/0_day5_glycation.csv")
        self.load_library("data/glycan_library.csv")

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
            filename, _, self.last_path = get_filename(
                self, "open", "Load glycation data …",
                self.last_path, FileTypes(["csv"]))
            if filename is None:
                return

        logging.info("Loading glycation data in '{}'".format(filename))
        try:
            self.se_glycation = read_clean_datasets(filename)
        except (OSError, ValueError) as e:
            QMessageBox.critical(self, "Error", str(e))
            return

        # extract x- and y-values from series
        x_values = [str(i) for i in self.se_glycation.index]
        y_values = [a.nominal_value for a in self.se_glycation]

        # assemble the chart
        bar_set = QBarSet("glycation abundance")
        bar_set.append(y_values)
        bar_set.hovered.connect(self.update_glycation_label)
        bar_series = QBarSeries()
        bar_series.append(bar_set)

        x_axis = QBarCategoryAxis()
        x_axis.append(x_values)
        x_axis.setTitleText("count")

        y_axis = QValueAxis()
        y_axis.setRange(0, 100)
        y_axis.setTitleText("abundance")
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
            self.lbGlycation.setText(
                "Abundance: <b>{:.2f}</b> ± {:.2f} %".format(
                    self.se_glycation[bar_index].nominal_value,
                    self.se_glycation[bar_index].std_dev))
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
            filename, _, self.last_path = get_filename(
                self, "open", "Load glycoform data …",
                self.last_path, FileTypes(["csv"]))
            if filename is None:
                return

        logging.info("Loading glycoform data in '{}'".format(filename))
        try:
            self.se_glycoforms = (read_clean_datasets(filename)
                                  .sort_values(ascending=False))
        except (OSError, ValueError) as e:
            QMessageBox.critical(self, "Error", str(e))
            return

        for widget in (self.cbAggGlycoforms,
                       self.sbAggGlycoforms,
                       self.lbAggGlycoforms):
            widget.setEnabled(True)
        self.sbAggGlycoforms.setMaximum(len(self.se_glycoforms) - 2)
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
                    self.se_glycoforms_agg.index[bar_index],
                    self.se_glycoforms_agg[bar_index].nominal_value,
                    self.se_glycoforms_agg[bar_index].std_dev))
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
            agg_abundance = (self.se_glycoforms
                             .iloc[self.sbAggGlycoforms.value():]
                             .sum())
            self.se_glycoforms_agg = (
                self.se_glycoforms
                    .iloc[:self.sbAggGlycoforms.value()]
                    .append(pd.Series(agg_abundance, index=["other"])))

        # extract x- and y-values from series
        x_values = [str(i) for i in self.se_glycoforms_agg.index]
        y_values = [a.nominal_value for a in self.se_glycoforms_agg]

        # assemble the chart
        bar_set = QBarSet("glycoform abundance")
        bar_set.append(y_values)
        bar_set.hovered.connect(self.update_glycoform_label)
        bar_series = QBarSeries()
        bar_series.append(bar_set)

        x_axis = QBarCategoryAxis()
        x_axis.append(x_values)
        x_axis.setTitleVisible(False)
        x_axis.setLabelsAngle(270)

        y_axis = QValueAxis()
        y_axis.setRange(0, max(self.se_glycoforms_agg).nominal_value)
        y_axis.setTitleText("abundance")
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
            filename, _, self.last_path = get_filename(
                self, "open", "Load glycan library …",
                self.last_path, FileTypes(["csv"]))
            if filename is None:
                return

        logging.info("Loading glycan library in '{}'".format(filename))
        try:
            self.se_library = read_library(filename)
        except (OSError, ValueError) as e:
            QMessageBox.critical(self, "Error", str(e))
            return

        # fill the table
        self.twLibrary.clearContents()
        for row_id, row in self.se_library.iterrows():
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

        missing_input = ["Please provide the following input data:"]
        if self.se_glycoforms is None:
            missing_input.append("glycoforms")
        if self.se_glycation is None:
            missing_input.append("glycation")
        if len(missing_input) > 1:
            QMessageBox.critical(
                self, "Error", "\n- ".join(missing_input))
            return

        logging.info("Correcting dataset  …")
        QApplication.processEvents()
        try:
            G = GlycationGraph(glycan_library=self.se_library,
                               glycoforms=self.se_glycoforms,
                               glycation=self.se_glycation)
            G.correct_abundances()
            self.results = G.to_dataframe()
        except ValueError as e:
            QMessageBox.critical(self, "Error", str(e))
        logging.info("… done!")
        self.show_results()

    def show_results(self):
        # fill the table
        self.twResults.clearContents()
        for row_id, row in self.results.iterrows():
            self.twResults.insertRow(row_id)
            item = SortableTableWidgetItem(row.iloc[0].split(" or ", 1)[0])
            item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
            self.twResults.setItem(row_id, 0, item)

            for col_id in range(1, 5):
                item = SortableTableWidgetItem("{:.2f}".format(row.iloc[col_id]))
                item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
                self.twResults.setItem(row_id, col_id, item)


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
