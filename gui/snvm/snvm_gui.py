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

# mergeme
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

class OptimizerSettingDialog(QtWidgets.QDialog):
    """ User configurable settings for the optimizer embedded in cofocal gui"""

    def __init__(self):
        # Get the path to the *.ui file
        this_dir = os.path.dirname(__file__)
        ui_file = os.path.join(this_dir, 'ui_snvm_optim_settings.ui')

        # Load it
        super(OptimizerSettingDialog, self).__init__()
        uic.loadUi(ui_file, self)

class SnvmSettingDialog(QtWidgets.QDialog):
    """ User configurable settings for the optimizer embedded in cofocal gui"""

    def __init__(self):
        # Get the path to the *.ui file
        this_dir = os.path.dirname(__file__)
        ui_file = os.path.join(this_dir, 'ui_snvm_settings.ui')

        # Load it
        super(SnvmSettingDialog, self).__init__()
        uic.loadUi(ui_file, self)

class SnvmGui(GUIBase):
    """ Main Confocal Class for xy and depth scans.
    """
    #TODO: GENERAL TODO LIST
    # -Implement the selection of the area like the one in the original qudi confocal scanner

    # declare connectors
    snvm_logic = Connector(interface='SnvmLogic')
    optimizer_logic = Connector(interface='OptimizerLogic')

    default_meter_prefix = ConfigOption('default_meter_prefix', None)  # assume the unit prefix of position spinbox
    startstopFreq_multiplier = ConfigOption('startstopFreq_multiplier', default=1e9, missing='info')
    stepFreq_multiplier = ConfigOption('stepFreq_multiplier', default=1e6, missing='info')
    xy_range_multiplier = ConfigOption('xy_range_multiplier', default=1e-9, missing='info')
    px_time_multiplier = ConfigOption('px_time_multiplier', default=1e-3, missing='info')
    cbar_count_multiplier = ConfigOption('cbar_count_multiplier', default=1e-3, missing='info')
    cbar_afm_multiplier = ConfigOption('cbar_afm_multplier', default=20/3, missing='info')
    mwspinbox_float_resolution = ConfigOption('mwspinbox_float_resolution', default=6)

    # signals
    sigStartOptimizer = QtCore.Signal()
    sigStartScanningSnvm = QtCore.Signal()
    sigStartScanningConf = QtCore.Signal()
    sigGoTo = QtCore.Signal(str)

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
        self.initOptimizer()
        self.initSnvmSettings()

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
        self._mainwindow.multiFreqPlotView.setLabel('bottom', 'X (nm)')
        self._mainwindow.multiFreqPlotView.setLabel('left', 'Y (nm)')
        self._mainwindow.multiFreqPlotView.toggle_crosshair(True, movable=True)
        self._mainwindow.multiFreqPlotView.set_crosshair_size((1,1))
        self._mainwindow.multiFreqPlotView.sigCrosshairDraggedPosChanged.connect(self.move_afm_crosshair)
        self._mainwindow.multiFreqPlotView.toggle_selection(True)
        self._mainwindow.multiFreqPlotView.sigMouseAreaSelected.connect(self.update_scanning_range)

        snvm_im_vb = self.get_image_viewbox(self.snvm_image)
        snvm_im_vb.setAspectLocked(True)
        snvm_im_vb.toggle_selection(True)
        snvm_im_vb.toggle_zoom_by_selection(True)

        self.multifreq_cb = ColorBar(self.photon_colormap.cmap_normed, width=100, cb_min=0, cb_max=1)
        self._mainwindow.multiFreqCbarView.addItem(self.multifreq_cb)
        self._mainwindow.multiFreqCbarView.hideAxis('bottom')
        self._mainwindow.multiFreqCbarView.setMouseEnabled(x=False, y=False)

        #Set up the AFM image and colorbar
        self.afm_image = ScanImageItem(axisOrder='row-major')
        self._mainwindow.afmPlotView.addItem(self.afm_image)
        self._mainwindow.afmPlotView.setLabel('bottom', 'X (nm)')
        self._mainwindow.afmPlotView.setLabel('left', 'Y (nm)')
        self._mainwindow.afmPlotView.toggle_crosshair(True, movable=True)
        self._mainwindow.afmPlotView.set_crosshair_size((1, 1))
        self._mainwindow.afmPlotView.sigCrosshairDraggedPosChanged.connect(self.move_multifreq_crosshair)
        self._mainwindow.afmPlotView.toggle_selection(True)
        self._mainwindow.afmPlotView.sigMouseAreaSelected.connect(self.update_scanning_range)

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
        self._mainwindow.confocalScannerView.setLabel('bottom', 'X (nm)')
        self._mainwindow.confocalScannerView.setLabel('left', 'Y (nm)')
        self._mainwindow.confocalScannerView.toggle_crosshair(True, movable=True)
        self._mainwindow.confocalScannerView.set_crosshair_size((1, 1))
        self._mainwindow.confocalScannerView.toggle_selection(True)
        self._mainwindow.confocalScannerView.sigMouseAreaSelected.connect(self.update_scanning_range)

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
        self._mainwindow.optimizerView.setLabel('bottom', 'X (nm)')
        self._mainwindow.optimizerView.setLabel('left', 'Y (nm)')

        opt_im_vb = self.get_image_viewbox(self.optimizer_image)
        opt_im_vb.setAspectLocked(True)

        self.opt_cb = ColorBar(self.photon_colormap.cmap_normed, width=100, cb_min=0, cb_max=1)
        self._mainwindow.optimizerCbarView.addItem(self.opt_cb)


        #Set up the ODMR plot
        self.curr_odmr_trace = pg.PlotDataItem(skipFiniteCheck=False, connect='finite', pen=pg.mkPen(color='w'))
        self.average_odmr_trace = pg.PlotDataItem(skipFiniteCheck=True, pen=pg.mkPen(color='r'))
        self._mainwindow.odmrPlotWidget.addItem(self.curr_odmr_trace)
        self._mainwindow.odmrPlotWidget.addItem(self.average_odmr_trace)
        self._mainwindow.odmrPlotWidget.setLabel('bottom', 'Frequency (GHz)')
        self._mainwindow.odmrPlotWidget.setLabel('left', 'Counts (GHz)')

        #Quick settings for the spinbox to view the frequency slices
        self._mainwindow.frequencySliceSelector.lineEdit().setReadOnly(True)

        self._viewIndex = 0 #Variable used to scroll through the SNVM images.Gets updated when clicking the frequency selector


        # Init scanning settings
        self._init_scanning_settings()

        # Init odmr settings
        self._init_odmr_settings()


        #Connect the signals
        self.sigStartOptimizer.connect(self.optimize_counts, QtCore.Qt.QueuedConnection)
        self.sigStartScanningSnvm.connect(self.prepare_snvm_scan, QtCore.Qt.QueuedConnection)
        self.sigStartScanningConf.connect(self.prepare_conf_scan, QtCore.Qt.QueuedConnection)
        self.sigGoTo.connect(self.go_to_point, QtCore.Qt.QueuedConnection)

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
        self._scanning_logic.signal_moved_to_point.connect(self.go_to_finished)

        self._optimizer_logic.sigImageUpdated.connect(self.refresh_optimizer_image)
        self._optimizer_logic.sigRefocusStarted.connect(self.set_optimizer_im_range)
        self._optimizer_logic.sigRefocusFinished.connect(self._optimization_complete)

        self._mainwindow.frequencySliceSelector.stepClicked.connect(self.frequency_selector_clicked)
        self._mainwindow.sampleTraceViewSpinBox.valueChanged.connect(self.refresh_snvm_image)
        self._mainwindow.sampleTraceViewSpinBox.valueChanged.connect(self.refresh_afm_image)
        self._mainwindow.tipTraceViewSpinBox.valueChanged.connect(self.refresh_confocal_image)
        self._mainwindow.sampleXSlider.sliderMoved.connect(self.slider_move_sample_crosshair)
        self._mainwindow.sampleYSlider.sliderMoved.connect(self.slider_move_sample_crosshair)
        self._mainwindow.tipXSlider.sliderMoved.connect(self.slider_move_tip_crosshair)
        self._mainwindow.tipYSlider.sliderMoved.connect(self.slider_move_tip_crosshair)


        ##############
        # Connect the actions to their slots
        ##############
        self._mainwindow.actionStart_snvm_scan.triggered.connect(self.scanning_action_clicked)
        self._mainwindow.actionStart_conf_scan.triggered.connect(self.scanning_action_clicked)
        self._mainwindow.actionStop_scan.triggered.connect(self.stop_scanning_request)
        self._mainwindow.actionOptimize.triggered.connect(self.scanning_action_clicked)
        self._mainwindow.actionOptimizer_settings.triggered.connect(self.menu_optimizer_settings)
        self._mainwindow.actionSnvm_settings.triggered.connect(self.menu_snvm_settings)
        self._mainwindow.action_snvm_goToPoint.triggered.connect(self.scanning_action_clicked)
        self._mainwindow.action_cfc_goToPoint.triggered.connect(self.scanning_action_clicked)
        self._mainwindow.actionSave_snvm.triggered.connect(self.save_snvm_data)
        self._mainwindow.actionSave_confocal.triggered.connect(self.save_confocal_data)

        self._mainwindow.actionStop_scan.setEnabled(False)
        self.show()

    def _init_odmr_settings(self):
        # Also here, store in a dictionary if the widgets need to be accessed in for loops
        self._odmr_widgets = dict()
        self._odmr_widgets[self._mainwindow.mwStart.objectName()] = self._mainwindow.mwStart
        self._odmr_widgets[self._mainwindow.mwEnd.objectName()] = self._mainwindow.mwEnd
        self._odmr_widgets[self._mainwindow.mwStep.objectName()] = self._mainwindow.mwStep
        self._odmr_widgets[self._mainwindow.mwPower.objectName()] = self._mainwindow.mwPower
        self._odmr_widgets[self._mainwindow.mwAverages.objectName()] = self._mainwindow.mwAverages

        self._mainwindow.mwStart.setDecimals(self.mwspinbox_float_resolution)
        self._mainwindow.mwEnd.setDecimals(self.mwspinbox_float_resolution)
        self._mainwindow.mwStep.setDecimals(self.mwspinbox_float_resolution)

        #########
        self._odmr_widgets['mwStart'].setValue(self._scanning_logic.start_freq /
                                               self.startstopFreq_multiplier)
        self._odmr_widgets['mwEnd'].setValue(self._scanning_logic.stop_freq /
                                             self.startstopFreq_multiplier)
        self._odmr_widgets['mwStep'].setValue(self._scanning_logic.freq_resolution /
                                              self.stepFreq_multiplier)
        self._odmr_widgets['mwPower'].setValue(self._scanning_logic.mw_power)
        self._odmr_widgets['mwAverages'].setValue(self._scanning_logic.odmr_averages)
        #########

    def _init_scanning_settings(self):
        ########
        # AFM scanning settings
        ########
        # Put all the settings in a dictionary, for ease of access
        self._afm_widgets = dict()
        # Confacal settings
        self._afm_widgets[self._mainwindow.xResolutionConf.objectName()] = self._mainwindow.xResolutionConf
        self._afm_widgets[self._mainwindow.yResolutionConf.objectName()] = self._mainwindow.yResolutionConf
        self._afm_widgets[self._mainwindow.xMinRangeConf.objectName()] = self._mainwindow.xMinRangeConf
        self._afm_widgets[self._mainwindow.xMaxRangeConf.objectName()] = self._mainwindow.xMaxRangeConf
        self._afm_widgets[self._mainwindow.yMinRangeConf.objectName()] = self._mainwindow.yMinRangeConf
        self._afm_widgets[self._mainwindow.yMaxRangeConf.objectName()] = self._mainwindow.yMaxRangeConf
        self._afm_widgets[self._mainwindow.fwpxTimeConf.objectName()] = self._mainwindow.fwpxTimeConf
        self._afm_widgets[self._mainwindow.storeRetraceConf.objectName()] = self._mainwindow.storeRetraceConf

        self._afm_widgets['yResolutionConf'].setValue(self._scanning_logic.scanning_y_resolution['tip'])
        self._afm_widgets['xResolutionConf'].setValue(self._scanning_logic.scanning_x_resolution['tip'])
        self._afm_widgets['xMinRangeConf'].setValue(self._scanning_logic.scanning_x_range['tip'][0] /
                                                    self.xy_range_multiplier)
        self._afm_widgets['yMinRangeConf'].setValue(self._scanning_logic.scanning_y_range['tip'][0] /
                                                    self.xy_range_multiplier)
        self._afm_widgets['xMaxRangeConf'].setValue(self._scanning_logic.scanning_x_range['tip'][1] /
                                                    self.xy_range_multiplier)
        self._afm_widgets['yMaxRangeConf'].setValue(self._scanning_logic.scanning_y_range['tip'][1] /
                                                    self.xy_range_multiplier)
        self._afm_widgets['fwpxTimeConf'].setValue(self._scanning_logic.px_time['tip'] /
                                                   self.px_time_multiplier)
        self._afm_widgets['storeRetraceConf'].setChecked(self._scanning_logic.store_retrace['tip'])

        # Snvm settings
        self._afm_widgets[self._mainwindow.xResolutionSnvm.objectName()] = self._mainwindow.xResolutionSnvm
        self._afm_widgets[self._mainwindow.yResolutionSnvm.objectName()] = self._mainwindow.yResolutionSnvm
        self._afm_widgets[self._mainwindow.xMinRangeSnvm.objectName()] = self._mainwindow.xMinRangeSnvm
        self._afm_widgets[self._mainwindow.xMaxRangeSnvm.objectName()] = self._mainwindow.xMaxRangeSnvm
        self._afm_widgets[self._mainwindow.yMinRangeSnvm.objectName()] = self._mainwindow.yMinRangeSnvm
        self._afm_widgets[self._mainwindow.yMaxRangeSnvm.objectName()] = self._mainwindow.yMaxRangeSnvm
        self._afm_widgets[self._mainwindow.fwpxTimeSnvm.objectName()] = self._mainwindow.fwpxTimeSnvm
        self._afm_widgets[self._mainwindow.storeRetraceSnvm.objectName()] = self._mainwindow.storeRetraceSnvm

        self._afm_widgets['xResolutionSnvm'].setValue(self._scanning_logic.scanning_x_resolution['sample'])
        self._afm_widgets['yResolutionSnvm'].setValue(self._scanning_logic.scanning_y_resolution['sample'])
        self._afm_widgets['xMinRangeSnvm'].setValue(self._scanning_logic.scanning_x_range['sample'][0] /
                                                    self.xy_range_multiplier)
        self._afm_widgets['yMinRangeSnvm'].setValue(self._scanning_logic.scanning_y_range['sample'][0] /
                                                    self.xy_range_multiplier)
        self._afm_widgets['xMaxRangeSnvm'].setValue(self._scanning_logic.scanning_x_range['sample'][1] /
                                                    self.xy_range_multiplier)
        self._afm_widgets['yMaxRangeSnvm'].setValue(self._scanning_logic.scanning_y_range['sample'][1] /
                                                    self.xy_range_multiplier)
        self._afm_widgets['fwpxTimeSnvm'].setValue(self._scanning_logic.px_time['sample'] /
                                                   self.px_time_multiplier)
        self._afm_widgets['storeRetraceSnvm'].setChecked(self._scanning_logic.store_retrace['sample'])

        self.sample_ranges = [self._scanning_logic.x_maxrange['sample'][0] / self.xy_range_multiplier,
                              self._scanning_logic.x_maxrange['sample'][1] / self.xy_range_multiplier], \
                             [self._scanning_logic.y_maxrange['sample'][0] / self.xy_range_multiplier,
                              self._scanning_logic.y_maxrange['sample'][1] / self.xy_range_multiplier]
        self.tip_ranges = [self._scanning_logic.x_maxrange['tip'][0] / self.xy_range_multiplier,
                           self._scanning_logic.x_maxrange['tip'][1] / self.xy_range_multiplier], \
                          [self._scanning_logic.y_maxrange['tip'][0] / self.xy_range_multiplier,
                           self._scanning_logic.y_maxrange['tip'][1] / self.xy_range_multiplier]

        self._mainwindow.sampleXSliderSpinBox.setRange(self.sample_ranges[0][0], self.sample_ranges[0][1])
        self._mainwindow.sampleYSliderSpinBox.setRange(self.sample_ranges[1][0], self.sample_ranges[1][1])
        self._mainwindow.tipXSliderSpinBox.setRange(self.tip_ranges[0][0], self.sample_ranges[0][1])
        self._mainwindow.tipYSliderSpinBox.setRange(self.tip_ranges[1][0], self.sample_ranges[1][1])

    def initOptimizer(self):
        """ Definition, configuration and initialisation of the optimizer settings GUI.

        This init connects all the graphic modules, which were created in the
        *.ui file and configures the event handling between the modules.
        Moreover it sets default values if not existed in the logic modules.
        """
        self._optim_dialog = OptimizerSettingDialog()
        # Connect the action of the settings window with the code:
        self._optim_dialog.accepted.connect(self.update_optimizer_settings)
        self._optim_dialog.rejected.connect(self.keep_former_optimizer_settings)
        self._optim_dialog.buttonBox.button(QtWidgets.QDialogButtonBox.Apply).clicked.connect(self.update_optimizer_settings)

        # Set up and connect xy channel combobox
        stacks = self._scanning_logic.get_stack_names()
        for n, stack in enumerate(stacks):
            self._optim_dialog.optimScanner_ComboBox.addItem(stack, n)

        # write the configuration to the settings window of the GUI.
        self.keep_former_optimizer_settings()

    def initOptimizer(self):
        self._optim_dialog = OptimizerSettingDialog()
        # Connect the action of the settings window with the code:
        self._optim_dialog.accepted.connect(self.update_optimizer_settings)
        self._optim_dialog.rejected.connect(self.keep_former_optimizer_settings)
        self._optim_dialog.buttonBox.button(QtWidgets.QDialogButtonBox.Apply).clicked.connect(self.update_optimizer_settings)

        # Set up and connect xy channel combobox
        stacks = self._scanning_logic.get_stack_names()
        for n, stack in enumerate(stacks):
            self._optim_dialog.optimScanner_ComboBox.addItem(stack, n)

        # write the configuration to the settings window of the GUI.
        self.keep_former_optimizer_settings()

    def initSnvmSettings(self):
        self._snvm_dialog = SnvmSettingDialog()
        # Connect the action of the settings window with the code:
        self._snvm_dialog.accepted.connect(self.update_snvm_settings)
        self._snvm_dialog.rejected.connect(self.keep_former_snvm_settings)
        self._snvm_dialog.buttonBox.button(QtWidgets.QDialogButtonBox.Apply).clicked.connect(
            self.update_optimizer_settings)

        # write the configuration to the settings window of the GUI.
        self.keep_former_snvm_settings()

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

    def prepare_snvm_scan(self):
        #Get the scanning settings from the GUI, and set them in the logic
        #FIXME: find a way to do this more efficiently, without calling each attribute one by one
        stk = self._scanning_logic.sampleStackName
        self._scanning_logic.store_retrace[stk] = True if self._afm_widgets['storeRetraceSnvm'].checkState() == 2 else False

        self._scanning_logic.scanning_x_range[stk] = [self._afm_widgets['xMinRangeSnvm'].value()*self.xy_range_multiplier,
                                                      self._afm_widgets['xMaxRangeSnvm'].value()*self.xy_range_multiplier]
        self._scanning_logic.scanning_y_range[stk] = [self._afm_widgets['yMinRangeSnvm'].value()*self.xy_range_multiplier,
                                                      self._afm_widgets['yMaxRangeSnvm'].value()*self.xy_range_multiplier]

        self._scanning_logic.scanning_x_resolution[stk] = self._afm_widgets['xResolutionSnvm'].value()
        self._scanning_logic.scanning_y_resolution[stk] = self._afm_widgets['yResolutionSnvm'].value()
        self._scanning_logic.optimize_while_scanning = self._optim_dialog.optimizeDuringScanCheckBox.isChecked()
        self._scanning_logic.every_N_pixels = self._optim_dialog.everyNPixelsDoubleSpinBox.value()

        #Set the integration time
        self._scanning_logic.px_time[stk] = self._afm_widgets['fwpxTimeSnvm'].value() * self.px_time_multiplier

        #First update the crosshair position
        crosshair_pos = self._mainwindow.multiFreqPlotView.crosshair_position
        if (crosshair_pos[0]*self.xy_range_multiplier not in self._scanning_logic.scanning_x_range[stk]
            or crosshair_pos[1]*self.xy_range_multiplier not in self._scanning_logic.scanning_y_range[stk]):

            newpos = (self._scanning_logic.scanning_x_range[stk][0]/self.xy_range_multiplier,
                      self._scanning_logic.scanning_y_range[stk][0]/self.xy_range_multiplier)
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

    def prepare_conf_scan(self):
        # Get the scanning settings from the GUI, and set them in the logic
        # FIXME: find a way to do this more efficiently, without calling each attribute one by one
        stk = self._scanning_logic.tipStackName
        self._scanning_logic.store_retrace[stk] = True if self._afm_widgets['storeRetraceConf'].checkState() == 2 else False

        self._scanning_logic.scanning_x_range[stk] = [self._afm_widgets['xMinRangeConf'].value() * self.xy_range_multiplier,
                                                      self._afm_widgets['xMaxRangeConf'].value() * self.xy_range_multiplier]
        self._scanning_logic.scanning_y_range[stk] = [self._afm_widgets['yMinRangeConf'].value() * self.xy_range_multiplier,
                                                      self._afm_widgets['yMaxRangeConf'].value() * self.xy_range_multiplier]

        self._scanning_logic.scanning_x_resolution[stk] = self._afm_widgets['xResolutionConf'].value()
        self._scanning_logic.scanning_y_resolution[stk] = self._afm_widgets['yResolutionConf'].value()
        self._scanning_logic.optimize_while_scanning = self._optim_dialog.optimizeDuringScanCheckBox.isChecked()

        # Set the integration time
        self._scanning_logic.px_time[stk] = self._afm_widgets['fwpxTimeConf'].value() * self.px_time_multiplier

        # First update the crosshair position
        crosshair_pos = self._mainwindow.confocalScannerView.crosshair_position
        if (crosshair_pos[0] not in self._scanning_logic.scanning_x_range[stk]
                or crosshair_pos[1] not in self._scanning_logic.scanning_y_range[stk]):
            newpos = self._scanning_logic.scanning_x_range[stk][0], self._scanning_logic.scanning_y_range[stk][0]
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
        self._mainwindow.actionStart_snvm_scan.setEnabled(False)
        self._mainwindow.actionStart_conf_scan.setEnabled(False)
        self._mainwindow.actionResume_snvm_scan.setEnabled(False)
        self._mainwindow.actionResume_conf_scan.setEnabled(False)
        self._mainwindow.actionOptimize.setEnabled(False)
        self._mainwindow.action_snvm_goToPoint.setEnabled(False)
        self._mainwindow.action_cfc_goToPoint.setEnabled(False)

        self._mainwindow.actionStop_scan.setEnabled(True)

        for setting in self._afm_widgets.values():
            setting.setEnabled(False)
        for setting in self._odmr_widgets.values():
            setting.setEnabled(False)

    def activate_interactions(self):
        self._mainwindow.actionStart_snvm_scan.setEnabled(True)
        self._mainwindow.actionStart_conf_scan.setEnabled(True)
        self._mainwindow.actionResume_snvm_scan.setEnabled(True)
        self._mainwindow.actionResume_conf_scan.setEnabled(True)
        self._mainwindow.actionOptimize.setEnabled(True)
        self._mainwindow.action_snvm_goToPoint.setEnabled(True)
        self._mainwindow.action_cfc_goToPoint.setEnabled(True)

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

        self._mainwindow.odmrPlotWidget.setXRange(startfreq, stopfreq)

    def stop_scanning_request(self):
        self._scanning_logic.stopRequested = True
        self._optimizer_logic.stop_refocus()

    def scanning_action_clicked(self):
        sendername = self.sender().objectName()
        self.disable_interactions()
        if sendername == 'actionOptimize':
            self.sigStartOptimizer.emit()
        elif sendername == 'action_snvm_goToPoint':
            self._mainwindow.actionStop_scan.setEnabled(False)
            self.sigGoTo.emit('snvm')
        elif sendername == 'action_cfc_goToPoint':
            self._mainwindow.actionStop_scan.setEnabled(False)
            self.sigGoTo.emit('cfc')
        elif sendername == 'actionStart_conf_scan':
            self.sigStartScanningConf.emit()
        elif sendername == 'actionStart_snvm_scan':
            self.sigStartScanningSnvm.emit()

    def optimize_counts(self):
        self.disable_interactions()

        crosshair_pos = self._mainwindow.confocalScannerView.get_crosshair_pos()
        crosshair_pos = [pos * self.xy_range_multiplier for pos in crosshair_pos]

        self._optimizer_logic.start_refocus(crosshair_pos)

    def _optimization_complete(self, coords):
        self._mainwindow.confocalScannerView.set_crosshair_pos((coords[0]/self.xy_range_multiplier,
                                                                    coords[1]/self.xy_range_multiplier))
        if not self._scanning_logic._snvm_active:
            self._scanning_logic.go_to_point(coords, stack=self._optimizer_logic.optimizer_stack)
            self.activate_interactions()

    def refresh_snvm_image(self):
        if self._mainwindow.viewtracesample:
            curr_image = self._scanning_logic.snvm_matrix.mean(axis=-1)
            curr_image = curr_image[:, :, self._viewIndex]
        else:
            curr_image = self._scanning_logic.snvm_matrix_retrace.mean(axis=-1)
            curr_image = curr_image[:, :, self._viewIndex]

        curr_image *= self.cbar_count_multiplier

        snvm_range = self.get_cb_range(curr_image)
        self.snvm_image.setImage(curr_image, levels=(snvm_range[0], snvm_range[1]))
        self.refresh_colorbar(cbar=self.multifreq_cb, cbar_range=snvm_range)

    def refresh_afm_image(self):
        if self._mainwindow.viewtracesample:
            curr_image = self._scanning_logic.xy_scan_matrix[:]
        else:
            curr_image = self._scanning_logic.xy_scan_matrix_retrace[:]

        afm_image = curr_image * self.cbar_afm_multiplier

        afm_range = self.get_cb_range(afm_image)

        self.afm_image.setImage(afm_image, levels=(afm_range[0], afm_range[1]))
        self.refresh_colorbar(cbar=self.afm_cb, cbar_range=afm_range)

    def refresh_odmr_plot(self, odmr_rep_index):
        curr_freq_matrix = self._scanning_logic.temp_freq_matrix[odmr_rep_index]
        self.curr_odmr_trace.setData(self._scanning_logic.freq_axis/self.startstopFreq_multiplier, curr_freq_matrix)
        if odmr_rep_index > 0:
            self.average_odmr_trace.setData(self._scanning_logic.freq_axis/self.startstopFreq_multiplier,
                                            self._scanning_logic.average_odmr_trace)
        else:
            self.average_odmr_trace.clear()

    def refresh_confocal_image(self):
        if self._mainwindow.viewtracetip:
            curr_image = self._scanning_logic.xy_scan_matrix[:]
        else:
            curr_image = self._scanning_logic.xy_scan_matrix_retrace[:]

        cfc_image = curr_image * self.cbar_count_multiplier

        cfc_range = self.get_cb_range(cfc_image)
        self.cfc_image.setImage(cfc_image, levels=(cfc_range[0], cfc_range[1]))
        self.refresh_colorbar(cbar=self.cfc_cb, cbar_range=cfc_range)

    def refresh_optimizer_image(self):
        curr_image = self._optimizer_logic.xy_refocus_image[:, :, 2]

        opt_image = curr_image * self.cbar_count_multiplier
        opt_range = self.get_cb_range(opt_image)

        self.optimizer_image.setImage(opt_image, levels=(opt_range[0], opt_range[1]))
        self.refresh_colorbar(cbar=self.opt_cb, cbar_range=opt_range)

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

    def set_optimizer_im_range(self):
        xmin, xmax = self._optimizer_logic.xy_refocus_image[0, np.array([0, -1]), 0] #x vals in the rows of the first slice
        ymin, ymax = self._optimizer_logic.xy_refocus_image[np.array([0, -1]), 0, 1] #y vals in the columns of the second one

        xpxsize = self._optimizer_logic.xy_refocus_image[0, 1, 0] - self._optimizer_logic.xy_refocus_image[0, 0, 0]
        ypxsize = self._optimizer_logic.xy_refocus_image[1, 0, 1] - self._optimizer_logic.xy_refocus_image[0, 0, 1]

        xmin, xmax, ymin, ymax, xpxsize, ypxsize = [val / self.xy_range_multiplier
                                                    for val in [xmin, xmax, ymin, ymax, xpxsize, ypxsize]]

        # FIXME: the following is a dirty trick to set the image range correctly (otherwise set_image_extent throws an error
        #  if no image has been previously loaded): how did I do for the other images?
        self.optimizer_image.setImage(self._optimizer_logic.xy_refocus_image[..., -1])

        self.optimizer_image.set_image_extent(((xmin-xpxsize/2, xmax+xpxsize/2), (ymin-ypxsize/2, ymax+ypxsize/2)))

    def update_optimizer_settings(self):
        self._optimizer_logic.refocus_XY_size = self._optim_dialog.xy_optimizer_range_DoubleSpinBox.value()*1e-6
        self._optimizer_logic.optimizer_XY_res = self._optim_dialog.xy_optimizer_resolution_SpinBox.value()
        self._optimizer_logic.integration_time = self._optim_dialog.intTime_SpinBox.value()
        value = self._optim_dialog.optimScanner_ComboBox.currentText()
        self._optimizer_logic.optimizer_stack = value
        self._scanning_logic.optimize_while_scanning = self._optim_dialog.optimizeDuringScanCheckBox.isChecked()
        self._scanning_logic.every_N_pixels = self._optim_dialog.everyNPixelsDoubleSpinBox.value()

    def keep_former_optimizer_settings(self):
        self._optim_dialog.xy_optimizer_range_DoubleSpinBox.setValue(self._optimizer_logic.refocus_XY_size*1e6)
        self._optim_dialog.xy_optimizer_resolution_SpinBox.setValue(self._optimizer_logic.optimizer_XY_res)
        self._optim_dialog.intTime_SpinBox.setValue(self._optimizer_logic.integration_time)

        opt_stack = self._optimizer_logic.optimizer_stack
        index = self._optim_dialog.optimScanner_ComboBox.findText(opt_stack)
        self._optim_dialog.optimScanner_ComboBox.setCurrentIndex(index)

        self._optim_dialog.optimizeDuringScanCheckBox.setChecked(self._scanning_logic.optimize_while_scanning)
        self._optim_dialog.everyNPixelsDoubleSpinBox.setValue(self._scanning_logic.every_N_pixels)



    def update_snvm_settings(self):
        self._scanning_logic.set_motion_speed(self._snvm_dialog.slowSpeedConfSpinBox.value(), stack='tip')
        self._scanning_logic.set_motion_speed(self._snvm_dialog.slowSpeedSnvmSpinBox.value(), stack='sample')
        self._scanning_logic.set_slowmotion_clockrate(self._snvm_dialog.motionClockRate_Spinbox.value())

    def update_scanning_range(self, rect):
        x, y, width, height = rect.getRect()
        print('x, y, width, height = ', x, y, width, height)
        x_end = x + width
        y_end = y + height
        if x > x_end:
            x, x_end = x_end, x
        if y > y_end:
            y, y_end = y_end, y

        self._afm_widgets['xMinRange'].setValue(x)
        self._afm_widgets['yMinRange'].setValue(y)
        self._afm_widgets['xMaxRange'].setValue(x_end)
        self._afm_widgets['yMaxRange'].setValue(y_end)

    def get_cb_range(self, image):
        imrange = [0, 1]
        image_nonzero = image[image != 0]
        if len(image_nonzero) > 0:
            imrange = [image_nonzero.min(), image_nonzero.max()]
        if imrange[0] == imrange[1]:
            imrange[1] = imrange[0] * 1.01
        return imrange

    def refresh_colorbar(self, cbar, cbar_range):
        cbar.refresh_colorbar(*cbar_range)

    def keep_former_snvm_settings(self):
        self._snvm_dialog.slowSpeedConfSpinBox.setValue(self._scanning_logic.backward_speed_conf)
        self._snvm_dialog.slowSpeedSnvmSpinBox.setValue(self._scanning_logic.backward_speed_snvm)
        self._snvm_dialog.motionClockRate_Spinbox.setValue(self._scanning_logic.get_slowmotion_clockrate())

    def frequency_selector_clicked(self, freq_val):
        difference = ( (freq_val - self._odmr_widgets["mwStart"].value()) *
                       self.startstopFreq_multiplier / self.stepFreq_multiplier)
        index = round( difference / self._odmr_widgets['mwStep'].value())

        self._viewIndex = index
        self.refresh_snvm_image()

    def get_image_viewbox(self, imageitem):
        vb = imageitem.getViewBox()
        return vb

    def move_afm_crosshair(self, coordinates):
        self._mainwindow.afmPlotView.set_crosshair_pos(coordinates)

    def move_multifreq_crosshair(self, coordinates):
        self._mainwindow.multiFreqPlotView.set_crosshair_pos(coordinates)

    def menu_optimizer_settings(self):
        """ This method opens the settings menu. """
        self.keep_former_optimizer_settings()
        self._optim_dialog.exec_()

    def menu_snvm_settings(self):
        """ This method opens the settings menu. """
        self.keep_former_snvm_settings()
        self._snvm_dialog.exec_()

    def slider_move_sample_crosshair(self):
        x_coeff = self.sample_ranges[0][1] - self.sample_ranges[0][0]
        y_coeff = self.sample_ranges[1][1] - self.sample_ranges[1][0]

        pos_slider_x, pos_slider_y = self._mainwindow.sampleXSlider.value(), self._mainwindow.sampleYSlider.value()
        new_x = pos_slider_x * x_coeff / self._mainwindow.sampleXSlider.maximum()
        new_y = pos_slider_y * y_coeff / self._mainwindow.sampleYSlider.maximum()

        #self._mainwindow.sampleXSliderSpinBox.setValue(new_x)
        #self._mainwindow.sampleYSliderSpinBox.setValue(new_y)

        self._mainwindow.multiFreqPlotView.set_crosshair_pos((new_x, new_y))
        self._mainwindow.afmPlotView.set_crosshair_pos((new_x, new_y))

    def slider_move_tip_crosshair(self):
        x_coeff = self.tip_ranges[0][1] - self.tip_ranges[0][0]
        y_coeff = self.tip_ranges[1][1] - self.tip_ranges[1][0]


        pos_slider_x, pos_slider_y = self._mainwindow.tipXSlider.value(), self._mainwindow.tipYSlider.value()
        new_x = pos_slider_x * x_coeff / self._mainwindow.tipXSlider.maximum()
        new_y = pos_slider_y * y_coeff / self._mainwindow.tipYSlider.maximum()

        self._mainwindow.tipXSliderSpinBox.setValue(new_x)
        self._mainwindow.tipYSliderSpinBox.setValue(new_y)

        self._mainwindow.confocalScannerView.set_crosshair_pos((new_x, new_y))

    def go_to_point(self, scanner):
        if scanner == 'snvm':
            position = self._mainwindow.multiFreqPlotView.crosshair_position
            position = [pos * self.xy_range_multiplier for pos in position]
            self._scanning_logic.go_to_point(position, stack=self._scanning_logic.sampleStackName,
                                             caller="gui")
        elif scanner == 'cfc':
            position = self._mainwindow.confocalScannerView.crosshair_position
            position = [pos * self.xy_range_multiplier for pos in position]
            self._scanning_logic.go_to_point(position, stack=self._scanning_logic.tipStackName,
                                             caller="gui")
        else:
            self.log.exception("Invalid name.")

    def go_to_finished(self, callertag):
        if callertag=="gui":
            self.activate_interactions()
            self._mainwindow.actionStop_scan.setEnabled(True)
        else:
            pass

    def save_snvm_data(self):
        self._scanning_logic.save_snvm()

    def save_confocal_data(self):
        self._scanning_logic.save_confocal()

