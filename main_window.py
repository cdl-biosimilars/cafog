# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'main_window.ui'
#
# Created by: PyQt5 UI code generator 5.10.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1072, 821)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout_7 = QtWidgets.QVBoxLayout(self.centralwidget)
        self.verticalLayout_7.setObjectName("verticalLayout_7")
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.groupBox = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox.setObjectName("groupBox")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.groupBox)
        self.verticalLayout.setObjectName("verticalLayout")
        self.cvGlycoforms = QChartView(self.groupBox)
        self.cvGlycoforms.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.cvGlycoforms.setRenderHints(QtGui.QPainter.Antialiasing|QtGui.QPainter.TextAntialiasing)
        self.cvGlycoforms.setObjectName("cvGlycoforms")
        self.verticalLayout.addWidget(self.cvGlycoforms)
        self.horizontalLayout_7 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_7.setObjectName("horizontalLayout_7")
        self.cbAggGlycoforms = QtWidgets.QCheckBox(self.groupBox)
        self.cbAggGlycoforms.setEnabled(False)
        self.cbAggGlycoforms.setChecked(True)
        self.cbAggGlycoforms.setObjectName("cbAggGlycoforms")
        self.horizontalLayout_7.addWidget(self.cbAggGlycoforms)
        self.sbAggGlycoforms = QtWidgets.QSpinBox(self.groupBox)
        self.sbAggGlycoforms.setEnabled(False)
        self.sbAggGlycoforms.setProperty("value", 7)
        self.sbAggGlycoforms.setObjectName("sbAggGlycoforms")
        self.horizontalLayout_7.addWidget(self.sbAggGlycoforms)
        self.lbAggGlycoforms = QtWidgets.QLabel(self.groupBox)
        self.lbAggGlycoforms.setEnabled(False)
        self.lbAggGlycoforms.setObjectName("lbAggGlycoforms")
        self.horizontalLayout_7.addWidget(self.lbAggGlycoforms)
        spacerItem = QtWidgets.QSpacerItem(1, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_7.addItem(spacerItem)
        self.verticalLayout.addLayout(self.horizontalLayout_7)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.btLoadGlycoforms = QtWidgets.QPushButton(self.groupBox)
        self.btLoadGlycoforms.setObjectName("btLoadGlycoforms")
        self.horizontalLayout_2.addWidget(self.btLoadGlycoforms)
        spacerItem1 = QtWidgets.QSpacerItem(1, 17, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem1)
        self.lbGlycoform = QtWidgets.QLabel(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lbGlycoform.sizePolicy().hasHeightForWidth())
        self.lbGlycoform.setSizePolicy(sizePolicy)
        self.lbGlycoform.setObjectName("lbGlycoform")
        self.horizontalLayout_2.addWidget(self.lbGlycoform)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.horizontalLayout_3.addWidget(self.groupBox)
        self.groupBox_2 = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox_2.setObjectName("groupBox_2")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.groupBox_2)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.cvGlycation = QChartView(self.groupBox_2)
        self.cvGlycation.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.cvGlycation.setRenderHints(QtGui.QPainter.Antialiasing|QtGui.QPainter.TextAntialiasing)
        self.cvGlycation.setObjectName("cvGlycation")
        self.verticalLayout_2.addWidget(self.cvGlycation)
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.btLoadGlycation = QtWidgets.QPushButton(self.groupBox_2)
        self.btLoadGlycation.setObjectName("btLoadGlycation")
        self.horizontalLayout_4.addWidget(self.btLoadGlycation)
        spacerItem2 = QtWidgets.QSpacerItem(1, 17, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_4.addItem(spacerItem2)
        self.lbGlycation = QtWidgets.QLabel(self.groupBox_2)
        self.lbGlycation.setObjectName("lbGlycation")
        self.horizontalLayout_4.addWidget(self.lbGlycation)
        self.verticalLayout_2.addLayout(self.horizontalLayout_4)
        self.horizontalLayout_3.addWidget(self.groupBox_2)
        self.groupBox_3 = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox_3.setObjectName("groupBox_3")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.groupBox_3)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.twLibrary = QtWidgets.QTableWidget(self.groupBox_3)
        self.twLibrary.setAcceptDrops(True)
        self.twLibrary.setObjectName("twLibrary")
        self.twLibrary.setColumnCount(2)
        self.twLibrary.setRowCount(0)
        item = QtWidgets.QTableWidgetItem()
        self.twLibrary.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.twLibrary.setHorizontalHeaderItem(1, item)
        self.twLibrary.horizontalHeader().setStretchLastSection(True)
        self.verticalLayout_3.addWidget(self.twLibrary)
        self.horizontalLayout_5 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_5.setObjectName("horizontalLayout_5")
        self.btLoadLibrary = QtWidgets.QPushButton(self.groupBox_3)
        self.btLoadLibrary.setObjectName("btLoadLibrary")
        self.horizontalLayout_5.addWidget(self.btLoadLibrary)
        spacerItem3 = QtWidgets.QSpacerItem(1, 17, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_5.addItem(spacerItem3)
        self.verticalLayout_3.addLayout(self.horizontalLayout_5)
        self.horizontalLayout_3.addWidget(self.groupBox_3)
        self.verticalLayout_7.addLayout(self.horizontalLayout_3)
        self.groupBox_4 = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox_4.setObjectName("groupBox_4")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout(self.groupBox_4)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.horizontalLayout_9 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_9.setObjectName("horizontalLayout_9")
        self.verticalLayout_6 = QtWidgets.QVBoxLayout()
        self.verticalLayout_6.setSpacing(0)
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.cvResults = QChartView(self.groupBox_4)
        self.cvResults.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.cvResults.setRenderHints(QtGui.QPainter.Antialiasing|QtGui.QPainter.TextAntialiasing)
        self.cvResults.setObjectName("cvResults")
        self.verticalLayout_6.addWidget(self.cvResults)
        self.lbResults = QtWidgets.QLabel(self.groupBox_4)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lbResults.sizePolicy().hasHeightForWidth())
        self.lbResults.setSizePolicy(sizePolicy)
        self.lbResults.setMinimumSize(QtCore.QSize(0, 15))
        self.lbResults.setStyleSheet("QLabel {\n"
"background: rgb(255, 255, 255)\n"
"}")
        self.lbResults.setObjectName("lbResults")
        self.verticalLayout_6.addWidget(self.lbResults)
        self.horizontalLayout_9.addLayout(self.verticalLayout_6)
        self.twResults = QtWidgets.QTableWidget(self.groupBox_4)
        self.twResults.setObjectName("twResults")
        self.twResults.setColumnCount(6)
        self.twResults.setRowCount(0)
        item = QtWidgets.QTableWidgetItem()
        self.twResults.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.twResults.setHorizontalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.twResults.setHorizontalHeaderItem(2, item)
        item = QtWidgets.QTableWidgetItem()
        self.twResults.setHorizontalHeaderItem(3, item)
        item = QtWidgets.QTableWidgetItem()
        self.twResults.setHorizontalHeaderItem(4, item)
        item = QtWidgets.QTableWidgetItem()
        self.twResults.setHorizontalHeaderItem(5, item)
        self.twResults.horizontalHeader().setDefaultSectionSize(70)
        self.horizontalLayout_9.addWidget(self.twResults)
        self.verticalLayout_5.addLayout(self.horizontalLayout_9)
        self.horizontalLayout_8 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_8.setObjectName("horizontalLayout_8")
        self.cbAggResults = QtWidgets.QCheckBox(self.groupBox_4)
        self.cbAggResults.setEnabled(False)
        self.cbAggResults.setChecked(True)
        self.cbAggResults.setObjectName("cbAggResults")
        self.horizontalLayout_8.addWidget(self.cbAggResults)
        self.sbAggResults = QtWidgets.QSpinBox(self.groupBox_4)
        self.sbAggResults.setEnabled(False)
        self.sbAggResults.setProperty("value", 7)
        self.sbAggResults.setObjectName("sbAggResults")
        self.horizontalLayout_8.addWidget(self.sbAggResults)
        self.lbAggResults = QtWidgets.QLabel(self.groupBox_4)
        self.lbAggResults.setEnabled(False)
        self.lbAggResults.setObjectName("lbAggResults")
        self.horizontalLayout_8.addWidget(self.lbAggResults)
        spacerItem4 = QtWidgets.QSpacerItem(1, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_8.addItem(spacerItem4)
        self.btSaveResults = QtWidgets.QPushButton(self.groupBox_4)
        self.btSaveResults.setEnabled(False)
        self.btSaveResults.setObjectName("btSaveResults")
        self.horizontalLayout_8.addWidget(self.btSaveResults)
        self.btSaveGraph = QtWidgets.QPushButton(self.groupBox_4)
        self.btSaveGraph.setEnabled(False)
        self.btSaveGraph.setObjectName("btSaveGraph")
        self.horizontalLayout_8.addWidget(self.btSaveGraph)
        self.verticalLayout_5.addLayout(self.horizontalLayout_8)
        self.verticalLayout_7.addWidget(self.groupBox_4)
        self.groupBox_5 = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox_5.setObjectName("groupBox_5")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.groupBox_5)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.teLog = QtWidgets.QTextEdit(self.groupBox_5)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Ignored)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.teLog.sizePolicy().hasHeightForWidth())
        self.teLog.setSizePolicy(sizePolicy)
        self.teLog.setMinimumSize(QtCore.QSize(0, 80))
        self.teLog.setReadOnly(True)
        self.teLog.setObjectName("teLog")
        self.verticalLayout_4.addWidget(self.teLog)
        self.verticalLayout_7.addWidget(self.groupBox_5)
        self.horizontalLayout_6 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_6.setObjectName("horizontalLayout_6")
        self.btCorrect = QtWidgets.QPushButton(self.centralwidget)
        self.btCorrect.setObjectName("btCorrect")
        self.horizontalLayout_6.addWidget(self.btCorrect)
        spacerItem5 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_6.addItem(spacerItem5)
        self.btSampleData = QtWidgets.QPushButton(self.centralwidget)
        self.btSampleData.setObjectName("btSampleData")
        self.horizontalLayout_6.addWidget(self.btSampleData)
        self.btHelp = QtWidgets.QPushButton(self.centralwidget)
        self.btHelp.setObjectName("btHelp")
        self.horizontalLayout_6.addWidget(self.btHelp)
        self.btQuit = QtWidgets.QPushButton(self.centralwidget)
        self.btQuit.setObjectName("btQuit")
        self.horizontalLayout_6.addWidget(self.btQuit)
        self.verticalLayout_7.addLayout(self.horizontalLayout_6)
        self.groupBox_4.raise_()
        self.groupBox_5.raise_()
        MainWindow.setCentralWidget(self.centralwidget)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "cafog"))
        self.groupBox.setTitle(_translate("MainWindow", "Glycoforms"))
        self.cbAggGlycoforms.setText(_translate("MainWindow", "Only display the"))
        self.lbAggGlycoforms.setText(_translate("MainWindow", "most abundant glycoforms"))
        self.btLoadGlycoforms.setText(_translate("MainWindow", "Load ..."))
        self.lbGlycoform.setText(_translate("MainWindow", "(glycoform)"))
        self.groupBox_2.setTitle(_translate("MainWindow", "Glycation"))
        self.btLoadGlycation.setText(_translate("MainWindow", "Load ..."))
        self.lbGlycation.setText(_translate("MainWindow", "(glycation)"))
        self.groupBox_3.setTitle(_translate("MainWindow", "Glycan library"))
        item = self.twLibrary.horizontalHeaderItem(0)
        item.setText(_translate("MainWindow", "Glycan"))
        item = self.twLibrary.horizontalHeaderItem(1)
        item.setText(_translate("MainWindow", "Composition"))
        self.btLoadLibrary.setText(_translate("MainWindow", "Load ..."))
        self.groupBox_4.setTitle(_translate("MainWindow", "Results"))
        self.lbResults.setText(_translate("MainWindow", "(results)"))
        self.twResults.setSortingEnabled(True)
        item = self.twResults.horizontalHeaderItem(0)
        item.setText(_translate("MainWindow", "Glycoform"))
        item = self.twResults.horizontalHeaderItem(1)
        item.setText(_translate("MainWindow", "Observed"))
        item = self.twResults.horizontalHeaderItem(2)
        item.setText(_translate("MainWindow", "Error"))
        item = self.twResults.horizontalHeaderItem(3)
        item.setText(_translate("MainWindow", "Actual"))
        item = self.twResults.horizontalHeaderItem(4)
        item.setText(_translate("MainWindow", "Error"))
        item = self.twResults.horizontalHeaderItem(5)
        item.setText(_translate("MainWindow", "Change"))
        self.cbAggResults.setText(_translate("MainWindow", "Only display the"))
        self.lbAggResults.setText(_translate("MainWindow", "most abundant glycoforms"))
        self.btSaveResults.setText(_translate("MainWindow", "Save results ..."))
        self.btSaveGraph.setText(_translate("MainWindow", "Save glycation graph ..."))
        self.groupBox_5.setTitle(_translate("MainWindow", "Log"))
        self.btCorrect.setText(_translate("MainWindow", "Correct abundances"))
        self.btSampleData.setText(_translate("MainWindow", "Load sample data"))
        self.btHelp.setText(_translate("MainWindow", "Help"))
        self.btHelp.setShortcut(_translate("MainWindow", "F1"))
        self.btQuit.setText(_translate("MainWindow", "Quit"))

from PyQt5.QtChart import QChartView
