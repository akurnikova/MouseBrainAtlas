# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'DataManager_v4.ui'
#
# Created: Tue Feb  3 02:54:58 2015
#      by: PyQt4 UI code generator 4.10.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_DataManager(object):
    def setupUi(self, DataManager):
        DataManager.setObjectName(_fromUtf8("DataManager"))
        DataManager.resize(1024, 984)
        self.cWidget = QtGui.QWidget(DataManager)
        self.cWidget.setObjectName(_fromUtf8("cWidget"))
        self.verticalLayout = QtGui.QVBoxLayout(self.cWidget)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.vLayout = QtGui.QVBoxLayout()
        self.vLayout.setSizeConstraint(QtGui.QLayout.SetDefaultConstraint)
        self.vLayout.setObjectName(_fromUtf8("vLayout"))
        self.label_2 = QtGui.QLabel(self.cWidget)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.vLayout.addWidget(self.label_2)
        self.topLayout = QtGui.QHBoxLayout()
        self.topLayout.setObjectName(_fromUtf8("topLayout"))
        self.leftListLayout = QtGui.QVBoxLayout()
        self.leftListLayout.setObjectName(_fromUtf8("leftListLayout"))
        self.stack_list = QtGui.QListView(self.cWidget)
        self.stack_list.setObjectName(_fromUtf8("stack_list"))
        self.leftListLayout.addWidget(self.stack_list)
        self.label_3 = QtGui.QLabel(self.cWidget)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.leftListLayout.addWidget(self.label_3)
        self.section_list = QtGui.QListView(self.cWidget)
        self.section_list.setObjectName(_fromUtf8("section_list"))
        self.leftListLayout.addWidget(self.section_list)
        self.topLayout.addLayout(self.leftListLayout)
        self.vLayout.addLayout(self.topLayout)
        self.bottomLayout = QtGui.QHBoxLayout()
        self.bottomLayout.setObjectName(_fromUtf8("bottomLayout"))
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.bottomLayout.addItem(spacerItem)
        self.label_4 = QtGui.QLabel(self.cWidget)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.bottomLayout.addWidget(self.label_4)
        self.comboBoxBrowseMode = QtGui.QComboBox(self.cWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.comboBoxBrowseMode.sizePolicy().hasHeightForWidth())
        self.comboBoxBrowseMode.setSizePolicy(sizePolicy)
        self.comboBoxBrowseMode.setSizeAdjustPolicy(QtGui.QComboBox.AdjustToContents)
        self.comboBoxBrowseMode.setObjectName(_fromUtf8("comboBoxBrowseMode"))
        self.bottomLayout.addWidget(self.comboBoxBrowseMode)
        self.buttonParams = QtGui.QPushButton(self.cWidget)
        self.buttonParams.setObjectName(_fromUtf8("buttonParams"))
        self.bottomLayout.addWidget(self.buttonParams)
        self.buttonQuit = QtGui.QPushButton(self.cWidget)
        self.buttonQuit.setObjectName(_fromUtf8("buttonQuit"))
        self.bottomLayout.addWidget(self.buttonQuit)
        self.vLayout.addLayout(self.bottomLayout)
        self.verticalLayout.addLayout(self.vLayout)
        DataManager.setCentralWidget(self.cWidget)
        self.menubar = QtGui.QMenuBar(DataManager)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1024, 25))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        self.menuView = QtGui.QMenu(self.menubar)
        self.menuView.setObjectName(_fromUtf8("menuView"))
        DataManager.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(DataManager)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        DataManager.setStatusBar(self.statusbar)
        self.actionStack_Section = QtGui.QAction(DataManager)
        self.actionStack_Section.setObjectName(_fromUtf8("actionStack_Section"))
        self.actionLabeling = QtGui.QAction(DataManager)
        self.actionLabeling.setObjectName(_fromUtf8("actionLabeling"))
        self.menuView.addAction(self.actionLabeling)
        self.menubar.addAction(self.menuView.menuAction())

        self.retranslateUi(DataManager)
        QtCore.QMetaObject.connectSlotsByName(DataManager)

    def retranslateUi(self, DataManager):
        DataManager.setWindowTitle(_translate("DataManager", "Data Manager", None))
        self.label_2.setText(_translate("DataManager", "Stack", None))
        self.label_3.setText(_translate("DataManager", "Section", None))
        self.label_4.setText(_translate("DataManager", "Filter by Label Name:", None))
        self.buttonParams.setText(_translate("DataManager", "Manage Parameters", None))
        self.buttonQuit.setText(_translate("DataManager", "Quit", None))
        self.menuView.setTitle(_translate("DataManager", "View", None))
        self.actionStack_Section.setText(_translate("DataManager", "Show Original Images", None))
        self.actionLabeling.setText(_translate("DataManager", "Shows Labelings", None))


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    DataManager = QtGui.QMainWindow()
    ui = Ui_DataManager()
    ui.setupUi(DataManager)
    DataManager.show()
    sys.exit(app.exec_())

