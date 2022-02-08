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
    snvm_logic = Connector(interface='SnvmLogic')

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
        self._scanning_logic = self.snvm_logic()
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
        # Connect the actions to their slots
        ##############
        self._mainwindow.actionStart_snvmscan.triggered.connect(self.start_scanning)
        self._mainwindow.actionStart_conf_scan.triggered.connect(self.start_scanning)
        self._mainwindow.actionStop_scan.triggered.connect(self.stop_scanning)

        ########
        # AFM scanning settings
        ########
        #Put all the settings in a dictionary, for ease of access
        self._afm_widgets = dict()
        self._afm_widgets[self._mainwindow.xResolution.objectName()] = self._mainwindow.xResolution
        self._afm_widgets[self._mainwindow.yResolution.objectName()] = self._mainwindow.yResolution
        self._afm_widgets[self._mainwindow.xMinRange.objectName()] = self._mainwindow.xMinRange
        self._afm_widgets[self._mainwindow.xMaxRange.objectName()] = self._mainwindow.xMaxRange
        self._afm_widgets[self._mainwindow.yMinRange.objectName()] = self._mainwindow.yMinRange
        self._afm_widgets[self._mainwindow.yMaxRange.objectName()] = self._mainwindow.yMaxRange
        self._afm_widgets[self._mainwindow.fwpxTime.objectName()] = self._mainwindow.fwpxTime
        self._afm_widgets[self._mainwindow.bwPxTime.objectName()] = self._mainwindow.bwPxTime
        self._afm_widgets[self._mainwindow.storeRetrace.objectName()] = self._mainwindow.storeRetrace

        #######
        # ODMR scanning settings
        ######
        #Also here, store in a dictionary if the widgets need to be accessed in for loops
        self._odmr_widgets = dict()
        self._odmr_widgets[self._mainwindow.mwStart.objectName()] = self._mainwindow.mwStart
        self._odmr_widgets[self._mainwindow.mwEnd.objectName()] = self._mainwindow.mwEnd
        self._odmr_widgets[self._mainwindow.mwStep.objectName()] = self._mainwindow.mwStep
        self._odmr_widgets[self._mainwindow.mwPower.objectName()] = self._mainwindow.mwPower
        self._odmr_widgets[self._mainwindow.mwAverages.objectName()] = self._mainwindow.mwAverages

        #TODO: maybe turn the freq resolution of the GUI into a settable value
        self._mainwindow.mwStart.setDecimals(6)
        self._mainwindow.mwEnd.setDecimals(6)
        self._mainwindow.mwStep.setDecimals(6)

        #Connect the signals
        self._odmr_widgets['mwStart'].valueChanged.connect(self.accept_set_frequency_ranges)
        self._odmr_widgets['mwEnd'].valueChanged.connect(self.accept_set_frequency_ranges)
        self._odmr_widgets['mwStep'].valueChanged.connect(self.accept_set_frequency_ranges)

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

    def start_scanning(self):
        self.disable_interactions()

        #Get the scanning settings from the GUI, and set them in the logic
        #FIXME: find a way to do this more efficiently, without calling each attribute one by one
        self._scanning_logic.store_retrace = self._afm_widgets['storeRetrace'].value()

        self._scanning_logic.scanning_x_range = [self._afm_widgets['xMinRange'].value(),
                                                 self._afm_widgets['xMaxRange'].value()]
        self._scanning_logic.scanning_y_range = [self._afm_widgets['yMinRange'].value(),
                                                 self._afm_widgets['yMaxRange'].value()]
        self._scanning_logic.scanning_x_resolution = self._afm_widgets['xResolution'].value()
        self._scanning_logic.scanning_y_resolution = self._afm_widgets['yResolution'].value()


        start_name = self.sender().objectName()
        if start_name == 'actionStart_snvmscan':
            self._scanning_logic.start_snvm_scanning()
        else:
            pass

    def stop_scanning(self):
        self.activate_interactions()

    def disable_interactions(self):
        self._mainwindow.actionStart_snvmscan.setEnabled(False)
        self._mainwindow.actionStart_conf_scan.setEnabled(False)
        self._mainwindow.actionResume_snvmscan.setEnabled(False)
        self._mainwindow.actionResume_conf_scan.setEnabled(False)

        for setting in self._afm_widgets.values():
            setting.setEnabled(False)
        for setting in self._odmr_widgets.values():
            setting.setEnabled(False)

    def activate_interactions(self):
        self._mainwindow.actionStart_snvmscan.setEnabled(True)
        self._mainwindow.actionStart_conf_scan.setEnabled(True)
        self._mainwindow.actionResume_snvmscan.setEnabled(True)
        self._mainwindow.actionResume_conf_scan.setEnabled(True)

        for setting in self._afm_widgets.values():
            setting.setEnabled(True)
        for setting in self._odmr_widgets.values():
            setting.setEnabled(True)


    def accept_set_frequency_ranges(self):
        """
        Function that checks that the range is start freq + step * multiple. If it's not, update the stop frequency.
        """
        stopfreq = self._odmr_widgets['mwEnd'].value()
        startfreq = self._odmr_widgets['mwStart'].value()
        freq_step = self._odmr_widgets['mwStep'].value()

        coeff_for_step = self.stepFreq_multiplier / self.startstopFreq_multiplier

        freq_diff = stopfreq - startfreq
        freq_step = freq_step * coeff_for_step
        if not (freq_diff % freq_step) == 0:
            multiple = round(freq_diff / freq_step)
            stop_freq = startfreq + multiple * freq_step
            self._odmr_widgets['mwEnd'].setValue(stop_freq)

        self._scanning_logic.start_freq = startfreq * self.startstopFreq_multiplier
        self._scanning_logic.stop_freq = stop_freq * self.startstopFreq_multiplier
        self._scanning_logic.freq_resolution = freq_step * self.stepFreq_multiplier