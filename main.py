import sys

import pandas as pd
from PyQt5.QtChart import (QBarCategoryAxis, QBarSeries, QBarSet,
                           QChart, QValueAxis)
from PyQt5.QtCore import QLocale, QMargins
from PyQt5.QtWidgets import QApplication, QMainWindow, QMessageBox, QWidget

from correction import read_clean_datasets

from main_window import Ui_MainWindow
from widgets import FileTypes, get_filename


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
        self.last_path = None

        # actions
        self.btLoadGlycation.clicked.connect(self.load_glycation)
        self.btQuit.clicked.connect(QApplication.instance().quit)

        # GUI modifications
        self.lbGlycation.setText("")
        self.lbGlycoform.setText("")
        self.cvGlycation.chart().setBackgroundRoundness(0)
        self.cvGlycation.chart().layout().setContentsMargins(0, 0, 0, 0)
        self.cvGlycoforms.chart().setBackgroundRoundness(0)
        self.cvGlycoforms.chart().layout().setContentsMargins(0, 0, 0, 0)
        self.cvResults.chart().setBackgroundRoundness(0)
        self.cvResults.chart().layout().setContentsMargins(0, 0, 0, 0)

    def load_glycation(self) -> None:
        """
        Load glycation data from a CSV file and display it
        in the corresponding chart view.

        :return: nothing, sets self.se_glycation
        :rtype: None
        """

        # load and clean glycation data
        filename, _, self.last_path = get_filename(
            self, "open", "Load glycation data …",
            self.last_path, FileTypes(["csv"]))
        if filename is None:
            return

        try:
            self.se_glycation = read_clean_datasets(filename)
        except (OSError, ValueError) as e:
            QMessageBox.critical(self, "Error", str(e))
            return

        # extract x- and y-values from series
        x_values = [str(i) for i in self.se_glycation.index]
        y_values = [a.nominal_value for a in self.se_glycation]

        # assemble the chart
        bar_set = QBarSet("bar_set")
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
        chart.setAnimationOptions(QChart.SeriesAnimations)
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


def _main():
    """
    Execute the main application loop.

    :return: nothing
    """
    app = QApplication(sys.argv)
    QLocale.setDefault(QLocale.c())
    frame = MainWindow()
    frame.show()
    app.exec_()


if __name__ == "__main__":
    _main()
