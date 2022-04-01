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
from gui.colordefs import ColorScaleInferno, BlackAndWhite
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
        ui_file = os.path.join(this_dir, 'ui_snvm_gui.ui')
        self._doubleclicked = False

        # Load it
        super(SnvmWindow, self).__init__()
        uic.loadUi(ui_file, self)

        #Set the trace/retrace spinboxes names
        btn_names = ["Trace", "Retrace"]
        self.sampleTraceViewSpinBox.setStrings(btn_names)
        self.sampleTraceViewSpinBox.lineEdit().setReadOnly(True)
        self.tipTraceViewSpinBox.setStrings(btn_names)
        self.tipTraceViewSpinBox.lineEdit().setReadOnly(True)

        self.viewtracesample = True
        self.viewtracetip = True

        self.sampleTraceViewSpinBox.valueChanged.connect(self.setSampleTraceRetrace)
        self.tipTraceViewSpinBox.valueChanged.connect(self.setTipTraceRetrace)

        self.show()

    def keyPressEvent(self, event):
        """Pass the keyboard press event from the main window further. """
        self.sigPressKeyBoard.emit(event)

    def mouseDoubleClickEvent(self, event):
        self._doubleclicked = True
        self.sigDoubleClick.emit()

    def setSampleTraceRetrace(self, value):
        self.viewtracesample = False if value == 1 else True

    def setTipTraceRetrace(self, value):
        self.viewtracetip = False if value == 1 else True

class SnvmGui(GUIBase):
    """ Main Confocal Class for xy and depth scans.
    """
    #TODO: GENERAL TODO LIST
    # -Fix range of optimizer plot whenever you click
    # -Check that the optimizer actually moves to the desired spot at the end of the scan
    # -Update the crosshair position once the optimization is done
    # -Implement slow motion to the beginning of the plot when you click snvm or confocal scan
    # -Implement the selection of the area like the one in the original qudi confocal scanner

    # declare connectors
    snvm_logic = Connector(interface='SnvmLogic')
    optimizer_logic = Connector(interface='OptimizerLogicPxScan')

    default_meter_prefix = ConfigOption('default_meter_prefix', None)  # assume the unit prefix of position spinbox
    # FIXME: for now I fix the multipliers, put is as an option, and update the labels accordingly in the GUI
    startstopFreq_multiplier = 1e9
    stepFreq_multiplier = 1e6
    xy_range_multiplier = 1e-9
    px_time_multiplier = 1e-3

    # signals
    sigStartOptimizer = QtCore.Signal()

    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)

    def on_activate(self):
        """ Initializes all needed UI files and establishes the connectors.

        This method executes the all the inits for the differnt GUIs and passes
        the event argument from fysom to the methods.
        """

        # Getting an access to all connectors:
        self._scanning_logic = self.snvm_logic()
        self._optimizer_logic = self.optimizer_logic()
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

        self.photon_colormap = ColorScaleInferno()
        self.afm_cmap = BlackAndWhite()
        self._crosshair_maxrange = None

        #Set up the SNVM image and colorbar
        self.snvm_image = ScanImageItem(axisOrder='row-major')
        self.snvm_image.setLookupTable(self.photon_colormap.lut)
        self._mainwindow.multiFreqPlotView.addItem(self.snvm_image)
        self._mainwindow.multiFreqPlotView.toggle_crosshair(True, movable=True)
        self._mainwindow.multiFreqPlotView.set_crosshair_size((1,1))
        self._mainwindow.multiFreqPlotView.sigCrosshairDraggedPosChanged.connect(self.move_afm_crosshair)

        snvm_im_vb = self.get_image_viewbox(self.snvm_image)
        snvm_im_vb.setAspectLocked(True)
        snvm_im_vb.toggle_selection(True)
        snvm_im_vb.toggle_zoom_by_selection(True)

        self.multifreq_cb = ColorBar(self.photon_colormap.cmap_normed, width=100, cb_min=0, cb_max=1)
        self._mainwindow.multiFreqCbarView.addItem(self.multifreq_cb)

        #Set up the AFM image and colorbar
        self.afm_image = ScanImageItem(axisOrder='row-major')
        self._mainwindow.afmPlotView.addItem(self.afm_image)
        self._mainwindow.afmPlotView.toggle_crosshair(True, movable=True)
        self._mainwindow.afmPlotView.set_crosshair_size((1, 1))
        self._mainwindow.afmPlotView.sigCrosshairDraggedPosChanged.connect(self.move_multifreq_crosshair)

        afm_im_vb = self.get_image_viewbox(self.afm_image)
        afm_im_vb.setAspectLocked(True)
        afm_im_vb.toggle_selection(True)
        afm_im_vb.toggle_zoom_by_selection(True)

        self.afm_cb = ColorBar(self.afm_cmap.cmap_normed, width=100, cb_min=0, cb_max=1)
        self._mainwindow.afmCbarView.addItem(self.afm_cb)

        # Set up the confocal image and colorbar
        self.cfc_image = ScanImageItem(axisOrder='row-major')
        self.cfc_image.setLookupTable(self.photon_colormap.lut)
        self._mainwindow.confocalScannerView.addItem(self.cfc_image)
        self._mainwindow.confocalScannerView.toggle_crosshair(True, movable=True)
        self._mainwindow.confocalScannerView.set_crosshair_size((1, 1))

        cfc_im_vb = self.get_image_viewbox(self.cfc_image)
        cfc_im_vb.setAspectLocked(True)
        cfc_im_vb.toggle_selection(True)
        cfc_im_vb.toggle_zoom_by_selection(True)

        self.cfc_cb = ColorBar(self.photon_colormap.cmap_normed, width=100, cb_min=0, cb_max=1)
        self._mainwindow.confocalCbarView.addItem(self.cfc_cb)

        # Set up the optimizer image and colorbar
        self.optimizer_image = ScanImageItem(axisOrder='row-major')
        self.optimizer_image.setLookupTable(self.photon_colormap.lut)
        self._mainwindow.optimizerView.addItem(self.optimizer_image)

        opt_im_vb = self.get_image_viewbox(self.optimizer_image)
        opt_im_vb.setAspectLocked(True)

        self.opt_cb = ColorBar(self.photon_colormap.cmap_normed, width=100, cb_min=0, cb_max=1)
        self._mainwindow.optimizerCbarView.addItem(self.opt_cb)

        #Set up the ODMR plot
        self.curr_odmr_trace = pg.PlotDataItem(skipFiniteCheck=False, connect='finite', pen=pg.mkPen(color='w'))
        self.average_odmr_trace = pg.PlotDataItem(skipFiniteCheck=True, pen=pg.mkPen(color='r'))
        self._mainwindow.odmrPlotWidget.addItem(self.curr_odmr_trace)
        self._mainwindow.odmrPlotWidget.addItem(self.average_odmr_trace)

        #Quick settings for the spinbox to view the frequency slices
        self._mainwindow.frequencySliceSelector.lineEdit().setReadOnly(True)

        self._viewIndex = 0 #Variable used to scroll through the SNVM images.Gets updated when clicking the frequency selector

        ##############
        # Connect the actions to their slots
        ##############
        self._mainwindow.actionStart_snvmscan.triggered.connect(self.start_scanning)
        self._mainwindow.actionStart_conf_scan.triggered.connect(self.start_scanning)
        self._mainwindow.actionStop_scan.triggered.connect(self.stop_scanning_request)
        self._mainwindow.actionOptimize.triggered.connect(self.optimize_clicked)

        self._mainwindow.actionStop_scan.setEnabled(False)

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
        self._afm_widgets[self._mainwindow.bwSpeed.objectName()] = self._mainwindow.bwSpeed
        self._afm_widgets[self._mainwindow.storeRetrace.objectName()] = self._mainwindow.storeRetrace

        # TODO: maybe set the maximum and minimum limits of the xmin/xmax and ymin/ymax
        #  from the maximum ranges allowed.

        #########
        # TODO: remove these defaults after debugging ended
        self._afm_widgets['xResolution'].setValue(5)
        self._afm_widgets['yResolution'].setValue(5)
        self._afm_widgets['xMaxRange'].setValue(40e2)
        self._afm_widgets['yMaxRange'].setValue(40e2)
        self._afm_widgets['fwpxTime'].setValue(10)
        self._afm_widgets['bwSpeed'].setValue(10)
        self._afm_widgets['storeRetrace'].setChecked(True)
        #########

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

        #########
        # TODO: remove these defaults after debugging ended
        self._odmr_widgets['mwStart'].setValue(2.87)
        self._odmr_widgets['mwEnd'].setValue(2.88)
        self._odmr_widgets['mwStep'].setValue(1)
        self._odmr_widgets['mwPower'].setValue(-100)
        self._odmr_widgets['mwAverages'].setValue(1)
        #########

        #Connect the signals

        self.sigStartOptimizer.connect(self.optimize_counts, QtCore.Qt.QueuedConnection)

        self._afm_widgets['storeRetrace'].stateChanged.connect(self.deactivate_speed_box)

        self._odmr_widgets['mwStart'].valueChanged.connect(self.accept_frequency_ranges)
        self._odmr_widgets['mwEnd'].valueChanged.connect(self.accept_frequency_ranges)
        self._odmr_widgets['mwStep'].valueChanged.connect(self.accept_frequency_ranges)

        self._scanning_logic.signal_scan_finished.connect(self.snvm_confocal_finished)
        self._scanning_logic.signal_freq_px_acquired.connect(self.refresh_odmr_plot)
        self._scanning_logic.signal_snvm_image_updated.connect(self.refresh_snvm_image)
        self._scanning_logic.signal_snvm_image_updated.connect(self.refresh_afm_image)
        self._scanning_logic.signal_xy_image_updated.connect(self.refresh_confocal_image)
        self._scanning_logic.signal_snvm_initialized.connect(self.set_snvm_im_range)
        self._scanning_logic.signal_confocal_initialized.connect(self.set_confocal_im_range)

        self._optimizer_logic.sigImageUpdated.connect(self.refresh_optimizer_image)
        self._optimizer_logic.sigRefocusFinished.connect(self.activate_interactions)

        self._mainwindow.frequencySliceSelector.stepClicked.connect(self.frequency_selector_clicked)
        self._mainwindow.sampleTraceViewSpinBox.valueChanged.connect(self.refresh_snvm_image)
        self._mainwindow.sampleTraceViewSpinBox.valueChanged.connect(self.refresh_afm_image)
        self._mainwindow.tipTraceViewSpinBox.valueChanged.connect(self.refresh_confocal_image)

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
        self._scanning_logic.store_retrace = True if self._afm_widgets['storeRetrace'].checkState()==2 else False

        self._scanning_logic.scanning_x_range = [self._afm_widgets['xMinRange'].value()*self.xy_range_multiplier,
                                                 self._afm_widgets['xMaxRange'].value()*self.xy_range_multiplier]
        self._scanning_logic.scanning_y_range = [self._afm_widgets['yMinRange'].value()*self.xy_range_multiplier,
                                                 self._afm_widgets['yMaxRange'].value()*self.xy_range_multiplier]


        self._scanning_logic.scanning_x_resolution = self._afm_widgets['xResolution'].value()
        self._scanning_logic.scanning_y_resolution = self._afm_widgets['yResolution'].value()

        self._scanning_logic.backward_speed = self._afm_widgets['bwSpeed'].value() * 1e-6 #Set it back to m/s, since it is um/s

        #Set the integration time
        self._scanning_logic.px_time = self._afm_widgets['fwpxTime'].value() * self.px_time_multiplier


        start_name = self.sender().objectName()
        if start_name == 'actionStart_snvmscan':

            #First update the crosshair position
            crosshair_pos = self._mainwindow.multiFreqPlotView.crosshair_position
            if (crosshair_pos[0]*self.xy_range_multiplier not in self._scanning_logic.scanning_x_range
                or crosshair_pos[1]*self.xy_range_multiplier not in self._scanning_logic.scanning_y_range):

                newpos = (self._scanning_logic.scanning_x_range[0]/self.xy_range_multiplier,
                          self._scanning_logic.scanning_y_range[0]/self.xy_range_multiplier)
                self._mainwindow.multiFreqPlotView.set_crosshair_pos(newpos)
                self._mainwindow.afmPlotView.set_crosshair_pos(newpos)

            self.set_odmr_settings()

            #Here put the settings for the  spin box.
            self._mainwindow.frequencySliceSelector.setMinimum(self._odmr_widgets['mwStart'].value())
            self._mainwindow.frequencySliceSelector.setMaximum(self._odmr_widgets['mwEnd'].value())
            step_val_ghz = self._odmr_widgets['mwStep'].value() * self.stepFreq_multiplier / self.startstopFreq_multiplier
            self._mainwindow.frequencySliceSelector.setSingleStep(step_val_ghz)
            self._viewIndex = 0
            self._mainwindow.frequencySliceSelector.setValue(self._odmr_widgets['mwStart'].value())

            self._scanning_logic.start_snvm_scanning()

        else:
            # First update the crosshair position
            crosshair_pos = self._mainwindow.confocalScannerView.crosshair_position
            if (crosshair_pos[0] not in self._scanning_logic.scanning_x_range
                    or crosshair_pos[1] not in self._scanning_logic.scanning_y_range):
                newpos = self._scanning_logic.scanning_x_range[0], self._scanning_logic.scanning_y_range[0]
                self._mainwindow.confocalScannerView.set_crosshair_pos(newpos)

            self._scanning_logic.start_confocal_scanning()

    def snvm_confocal_finished(self, was_snvm):
        self.activate_interactions()

        xrange, yrange = self._scanning_logic.get_xy_image_range(multiplier=1 / self.xy_range_multiplier)

        if was_snvm:
            self._mainwindow.multiFreqPlotView.set_crosshair_range([xrange, yrange])
            self._mainwindow.afmPlotView.set_crosshair_range([xrange, yrange])
        else:
            self._mainwindow.confocalScannerView.set_crosshair_range([xrange, yrange])

    def disable_interactions(self):
        self._mainwindow.actionStart_snvmscan.setEnabled(False)
        self._mainwindow.actionStart_conf_scan.setEnabled(False)
        self._mainwindow.actionResume_snvmscan.setEnabled(False)
        self._mainwindow.actionResume_conf_scan.setEnabled(False)
        self._mainwindow.actionOptimize.setEnabled(False)

        self._mainwindow.actionStop_scan.setEnabled(True)

        for setting in self._afm_widgets.values():
            setting.setEnabled(False)
        for setting in self._odmr_widgets.values():
            setting.setEnabled(False)

    def activate_interactions(self):
        self._mainwindow.actionStart_snvmscan.setEnabled(True)
        self._mainwindow.actionStart_conf_scan.setEnabled(True)
        self._mainwindow.actionResume_snvmscan.setEnabled(True)
        self._mainwindow.actionResume_conf_scan.setEnabled(True)
        self._mainwindow.actionOptimize.setEnabled(True)

        self._mainwindow.actionStop_scan.setEnabled(False)

        for setting in self._afm_widgets.values():
            setting.setEnabled(True)
        for setting in self._odmr_widgets.values():
            setting.setEnabled(True)

    def accept_frequency_ranges(self):
        """
        Function that checks that the range is start freq + step * multiple. If it's not, update the stop frequency.
        """
        stopfreq = self._odmr_widgets['mwEnd'].value()
        startfreq = self._odmr_widgets['mwStart'].value()
        freqstep = self._odmr_widgets['mwStep'].value()

        coeff_for_step = self.stepFreq_multiplier / self.startstopFreq_multiplier

        freq_diff = stopfreq - startfreq
        freqstep = freqstep * coeff_for_step
        if not (freq_diff % freqstep) == 0:
            multiple = round(freq_diff / freqstep)
            stop_freq = startfreq + multiple * freqstep
            self._odmr_widgets['mwEnd'].setValue(stop_freq)

    def set_odmr_settings(self):
        stopfreq = self._odmr_widgets['mwEnd'].value()
        startfreq = self._odmr_widgets['mwStart'].value()
        freqstep = self._odmr_widgets['mwStep'].value()
        power = self._odmr_widgets['mwPower'].value()
        averages = self._odmr_widgets['mwAverages'].value()

        self._scanning_logic.start_freq = startfreq * self.startstopFreq_multiplier
        self._scanning_logic.stop_freq = stopfreq * self.startstopFreq_multiplier
        self._scanning_logic.freq_resolution = freqstep * self.stepFreq_multiplier
        self._scanning_logic.mw_power = power
        self._scanning_logic.odmr_averages = averages

    def stop_scanning_request(self):
        self._scanning_logic.stopRequested = True
        self._optimizer_logic.stop_refocus()

    def optimize_clicked(self):
        """
        This intermediated function is used because if optimize_counts is called directly it does not let the
        disable_interactions() function to finish. This causes a lag in the disabling, that hangs while the scanner is
        moving to its initial position. Instead, emitting the signal and using a queued connection to optimize_counts
        solves the problem.
        """
        self.disable_interactions()
        self.sigStartOptimizer.emit()

    def optimize_counts(self):
        self.disable_interactions()
        self._optimizer_logic.set_bw_speed(self._afm_widgets['bwSpeed'].value() * 1e-6)
        self._optimizer_logic.set_samps_per_pixel(self._afm_widgets['fwpxTime'].value() * self.px_time_multiplier)

        crosshair_pos = self._mainwindow.confocalScannerView.get_crosshair_pos()
        crosshair_pos = [pos * 1e-9 for pos in crosshair_pos]

        self._optimizer_logic.start_refocus(crosshair_pos)

    def refresh_snvm_image(self):
        if self._mainwindow.viewtracesample:
            curr_image = self._scanning_logic.snvm_matrix.mean(axis=-1)
            curr_image = curr_image[:, :, self._viewIndex]
        else:
            curr_image = self._scanning_logic.snvm_matrix_retrace.mean(axis=-1)
            curr_image = curr_image[:, :, self._viewIndex]
        minmax = [curr_image.min(), curr_image.max()]
        self.snvm_image.setImage(curr_image)

    def refresh_afm_image(self):
        if self._mainwindow.viewtracesample:
            curr_image = self._scanning_logic.xy_scan_matrix
        else:
            curr_image = self._scanning_logic.xy_scan_matrix_retrace
        minmax = [curr_image.min(), curr_image.max()]
        self.afm_image.setImage(curr_image)

    def refresh_odmr_plot(self, odmr_rep_index):
        curr_freq_matrix = self._scanning_logic.temp_freq_matrix[odmr_rep_index]
        self.curr_odmr_trace.setData(self._scanning_logic.freq_axis, curr_freq_matrix)
        if odmr_rep_index > 0:
            self.average_odmr_trace.setData(self._scanning_logic.freq_axis, self._scanning_logic.average_odmr_trace)
        else:
            self.average_odmr_trace.clear()

    def refresh_confocal_image(self):
        if self._mainwindow.viewtracetip:
            curr_image = self._scanning_logic.xy_scan_matrix
        else:
            curr_image = self._scanning_logic.xy_scan_matrix_retrace

        minmax = [curr_image.min(), curr_image.max()]
        self.cfc_image.setImage(curr_image)

    def refresh_optimizer_image(self):
        curr_image = self._optimizer_logic.xy_refocus_image[:, :, 2]

        self.optimizer_image.setImage(curr_image)

    def set_snvm_im_range(self):
        im_range = self._scanning_logic.get_xy_image_range(multiplier=1/self.xy_range_multiplier)
        xmin, xmax  = im_range[0]
        ymin, ymax = im_range[1]

        xpxsize, ypxsize = self._scanning_logic.get_xy_step_size(multiplier=1/self.xy_range_multiplier)

        for image in [self.snvm_image, self.afm_image]:
            image.set_image_extent(((xmin-xpxsize/2, xmax+xpxsize/2), (ymin-ypxsize/2, ymax+ypxsize/2)))

    def set_confocal_im_range(self):
        im_range = self._scanning_logic.get_xy_image_range(multiplier=1/self.xy_range_multiplier)
        xmin, xmax = im_range[0]
        ymin, ymax = im_range[1]

        xpxsize, ypxsize = self._scanning_logic.get_xy_step_size(multiplier=1/self.xy_range_multiplier)

        self.cfc_image.set_image_extent(((xmin-xpxsize/2, xmax+xpxsize/2), (ymin-ypxsize/2, ymax+ypxsize/2)))

    def frequency_selector_clicked(self, freq_val):
        difference = ( (freq_val - self._odmr_widgets["mwStart"].value()) *
                       self.startstopFreq_multiplier / self.stepFreq_multiplier )
        index = round( difference / self._odmr_widgets['mwStep'].value() )

        self._viewIndex = index
        self.refresh_snvm_image()

    def get_image_viewbox(self, imageitem):
        vb = imageitem.getViewBox()
        return vb

    def deactivate_speed_box(self, checkbox_state):
        self._afm_widgets['bwSpeed'].setDisabled(checkbox_state)

    def move_afm_crosshair(self, coordinates):
        self._mainwindow.afmPlotView.set_crosshair_pos(coordinates)

    def move_multifreq_crosshair(self, coordinates):
        self._mainwindow.multiFreqPlotView.set_crosshair_pos(coordinates)