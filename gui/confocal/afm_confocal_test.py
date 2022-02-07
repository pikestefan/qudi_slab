# -*- coding: utf-8 -*-

"""
This file contains the Qudi GUI for general Confocal control.

Qudi is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Qudi is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Qudi. If not, see <http://www.gnu.org/licenses/>.

Copyright (c) the Qudi Developers. See the COPYRIGHT.txt file at the
top-level directory of this distribution and at <https://github.com/Ulm-IQO/qudi/>
"""

import numpy as np
import os
import pyqtgraph as pg
import time
from math import ceil

from core.connector import Connector
from core.configoption import ConfigOption
from core.statusvariable import StatusVar
from qtwidgets.scan_plotwidget import ScanImageItem
from gui.guibase import GUIBase
from gui.guiutils import ColorBar
from gui.colordefs import ColorScaleInferno
from gui.colordefs import QudiPalettePale as palette
from gui.fitsettings import FitParametersWidget
from qtpy import QtCore
from qtpy import QtGui
from qtpy import QtWidgets
from qtpy import uic


class SnvmWindow(QtWidgets.QMainWindow):
    """ Create the Mainwindow based on the corresponding *.ui file. """

    sigPressKeyBoard = QtCore.Signal(QtCore.QEvent)
    sigDoubleClick = QtCore.Signal()

    def __init__(self):
        # Get the path to the *.ui file
        this_dir = os.path.dirname(__file__)
        ui_file = os.path.join(this_dir, 'ui_afm_cfm_v3.ui')
        self._doubleclicked = False

        # Load it
        super(SnvmWindow, self).__init__()
        uic.loadUi(ui_file, self)
        self.show()

    def keyPressEvent(self, event):
        """Pass the keyboard press event from the main window further. """
        self.sigPressKeyBoard.emit(event)

    def mouseDoubleClickEvent(self, event):
        self._doubleclicked = True
        self.sigDoubleClick.emit()

class SnvmGui(GUIBase):
    """ Main Confocal Class for xy and depth scans.
    """

    # declare connectors
    #snvm_logic = Connector(interface='SnvmLogic')

    default_meter_prefix = ConfigOption('default_meter_prefix', None)  # assume the unit prefix of position spinbox
    # FIXME: for now I fix the multipliers, put is as an option, and updated the labels accordingly in the GUI
    startstopFreq_multiplier = 1e9
    stepFreq_multiplier = 1e6

    # signals
    sigStartOptimizer = QtCore.Signal(list, str)

    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)

    def on_activate(self):
        """ Initializes all needed UI files and establishes the connectors.

        This method executes the all the inits for the differnt GUIs and passes
        the event argument from fysom to the methods.
        """

        # Getting an access to all connectors:
        #self._scanning_logic = self.snvm_logic()
        self._hardware_state = True
        self.initMainUI()      # initialize the main GUI

    def initMainUI(self):
        """ Definition, configuration and initialisation of the confocal GUI.

        This init connects all the graphic modules, which were created in the
        *.ui file and configures the event handling between the modules.
        Moreover it sets default values.
        """
        self._mainwindow = SnvmWindow()

        # All our gui elements are dockable, and so there should be no "central" widget.
        self._mainwindow.centralwidget.hide()
        self._mainwindow.setDockNestingEnabled(True)


        self.colormap = ColorScaleInferno()
        self.sample_cb = ColorBar(self.colormap.cmap_normed, width=100, cb_min=0, cb_max=1)

        self._mainwindow.afmCbarView.addItem(self.sample_cb)
        ##############
        # Connect the buttons to their signals
        ##############
        #self._mainwindow.action_snvmscan.triggered.connect()

        ########
        # AFM scanning settings
        ########


        #######
        # ODMR scanning settings
        ######
        #TODO: maybe make the freq resolution of the GUI a settable value
        self._mainwindow.mwStart.setDecimals(6)
        self._mainwindow.mwEnd.setDecimals(6)
        self._mainwindow.mwStep.setDecimals(6)

        #Connect the signals
        self._mainwindow.mwStart.valueChanged.connect(self.accept_set_frequency_ranges)
        self._mainwindow.mwEnd.valueChanged.connect(self.accept_set_frequency_ranges)
        self._mainwindow.mwStep.valueChanged.connect(self.accept_set_frequency_ranges)

        self.show()

    def on_deactivate(self):
        """ Reverse steps of activation

        @return int: error code (0:OK, -1:error)
        """
        self._mainwindow.close()
        return 0

    def show(self):
        """Make main window visible and put it above all other windows. """
        # Show the Main Confocal GUI:
        self._mainwindow.show()
        self._mainwindow.activateWindow()
        self._mainwindow.raise_()


    def accept_set_frequency_ranges(self):
        pass
    #     """
    #     Function that checks that the range is start freq + step * multiple. If it's not, update the stop frequency.
    #     """
    #     stopfreq, startfreq = self._mainwindow.mwEnd.value(), self._mainwindow.mwStart.value()
    #     freq_step = self._mainwindow.mwStep.value()
    #
    #     coeff_for_step = self.stepFreq_multiplier / self.startstopFreq_multiplier
    #
    #     freq_diff = stopfreq - startfreq
    #     freq_step = freq_step * coeff_for_step
    #     if not (freq_diff % freq_step) == 0:
    #         multiple = round(freq_diff / freq_step)
    #         stop_freq = startfreq + multiple * freq_step
    #         self._mainwindow.mwEnd.setValue(stop_freq)
    #
    #     self._scanning_logic.start_freq = startfreq * self.startstopFreq_multiplier
    #     self._scanning_logic.stop_freq = stop_freq * self.startstopFreq_multiplier
    #     self._scanning_logic.freq_resolution = freq_step * self.stepFreq_multiplier