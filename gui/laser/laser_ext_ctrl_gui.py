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
        slider_widget = QtWidgets.QSlider()
        slider_widget.minimum = 0
        slider_widget.maximum = 1
        self.setCentralWidget(slider_widget)

        self.action_close.triggered.connect(self.close)


class LaserCtrlGUI(GUIBase):
    laser_ctrl_logic = Connector(interface='laser_ext_ctrl_logic')

    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)

    def on_activate(self):
        self._laser_ctrl_logic = self.laser_ctrl_logic()
        self._lw = LaserCtrlMainWindow()

        self._lw.slider_widget.sliderReleased.connect()

    def test(self):
        print("Slider moved")
