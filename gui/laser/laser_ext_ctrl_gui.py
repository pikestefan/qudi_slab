import numpy as np
import os
import pyqtgraph as pg
import time

from core.connector import Connector
from gui.colordefs import QudiPalettePale as palette
from gui.guibase import GUIBase
from interface.simple_laser_interface import ControlMode, ShutterState, LaserState
from qtpy import QtCore
from qtpy import QtWidgets
from qtpy import uic

class LaserCtrlMainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setWindowTitle('Test')
        self.slider_widget = QtWidgets.QSlider()
        self.slider_widget.setOrientation(QtCore.Qt.Horizontal)
        self.slider_widget.minimum = 0
        self.slider_widget.setMaximum(1000)
        self.setCentralWidget(self.slider_widget)

        self.action_close = QtWidgets.QAction('Close Window')

        self.action_close.triggered.connect(self.close)


class LaserCtrlGUI(GUIBase):
    laserctrllogic = Connector(interface='LaserExtCtrlLogic')

    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)

    def on_activate(self):
        self._laserctrllogic = self.laserctrllogic()
        self._lw = LaserCtrlMainWindow()

        self._lw.slider_widget.valueChanged.connect(self.test)

    def on_deactivate(self):
        print("Fucking off")

    @QtCore.Slot(int)
    def test(self, value):
        voltage = value/self._lw.slider_widget.maximum()
        self._laserctrllogic.set_power(voltage)


    def show(self):
        self._lw.show()
        self._lw.activateWindow()
        self._lw.raise_()
