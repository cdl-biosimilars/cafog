import sys


import pandas as pd
from PyQt5.QtChart import (QBarCategoryAxis, QBarSeries, QBarSet,
                           QChart, QValueAxis)
from PyQt5.QtCore import QLocale, QMargins
from PyQt5.QtWidgets import QApplication, QMainWindow, QMessageBox, QWidget


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
        self.df_glycation = None
        self.last_path = None

        # actions
        self.btLoadGlycation.clicked.connect(self.load_glycation)
        self.btQuit.clicked.connect(QApplication.instance().quit)

        self.cvGlycation.chart().setBackgroundRoundness(0)
        self.cvGlycation.chart().layout().setContentsMargins(0, 0, 0, 0)

    def load_glycation(self):
        # filename, _, self.last_path = get_filename(
        #     self, "open", "Load glycation data …",
        #     self.last_path, FileTypes(["csv"]))
        # if filename is None:
        #     return
        #
        # try:
        #     self.df_glycation = pd.read_csv(filename)
        # except OSError as e:
        #     QMessageBox.critical(
        #         self,
        #         "Error",
        #         "Error when reading {}: {}".format(filename, e))
        #     return

        bar_set = QBarSet("bar_set")
        bar_set.append([75, 15, 5, 3, 2])
        bar_series = QBarSeries()
        bar_series.append(bar_set)

        x_axis = QBarCategoryAxis()
        x_axis.append(["0", "1", "2", "3", "4"])
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
