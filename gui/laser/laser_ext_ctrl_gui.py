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
    """ The main window for the ODMR measurement GUI.
    """

    def __init__(self):
        # Get the path to the *.ui file
        this_dir = os.path.dirname(__file__)
        ui_file = os.path.join(this_dir, 'ui_laser_ext_ctrl.ui')

        # Load it
        super(LaserCtrlMainWindow, self).__init__()
        uic.loadUi(ui_file, self)
        self.show()

class LaserCtrlGUI(GUIBase):
    laserctrllogic = Connector(interface='LaserExtCtrlLogic')

    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)

    def on_activate(self):
        self._laserctrllogic = self.laserctrllogic()
        self._lw = LaserCtrlMainWindow()

        self._lw.slider_widget.valueChanged.connect(self.change_extctrl_power)
        self._lw.slider_widget.sliderReleased.connect(self.get_power_atsource)

        self.get_power_atsource()
        self.show()

    def on_deactivate(self):
        print("Fucking off")

    @QtCore.Slot(int)
    def change_extctrl_power(self, value):
        power_percentage = value/self._lw.slider_widget.maximum()
        self._laserctrllogic.set_ext_ctrl_power(power_percentage)

    def get_power_atsource(self):
        power = self._laserctrllogic.get_power_atsource()
        self._lw.actual_pw.setText("{:.2f}".format(power))

    def show(self):
        self._lw.show()
        self._lw.activateWindow()
        self._lw.raise_()