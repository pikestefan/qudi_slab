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
from core.connector import Connector
from core.configoption import ConfigOption
from qtwidgets.scan_plotwidget import ScanImageItem
from gui.guibase import GUIBase
from gui.guiutils import InteractiveColBar
from gui.colordefs import ColorScaleInferno, BlackAndWhite
from qtpy import QtCore
from qtpy import QtWidgets
from qtpy import uic


class SnvmWindow(QtWidgets.QMainWindow):
    """Create the Mainwindow based on the corresponding *.ui file."""

    sigPressKeyBoard = QtCore.Signal(QtCore.QEvent)
    sigDoubleClick = QtCore.Signal()

    def __init__(self):
        # Get the path to the *.ui file
        this_dir = os.path.dirname(__file__)
        ui_file = os.path.join(this_dir, "ui_snvm_gui.ui")
        self._doubleclicked = False

        # Load it
        super(SnvmWindow, self).__init__()
        uic.loadUi(ui_file, self)

        # Set the trace/retrace spinboxes names
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
        """Pass the keyboard press event from the main window further."""
        self.sigPressKeyBoard.emit(event)

    def mouseDoubleClickEvent(self, event):
        self._doubleclicked = True
        self.sigDoubleClick.emit()

    def setSampleTraceRetrace(self, value):
        self.viewtracesample = False if value == 1 else True

    def setTipTraceRetrace(self, value):
        self.viewtracetip = False if value == 1 else True


class OptimizerSettingDialog(QtWidgets.QDialog):
    """User configurable settings for the optimizer embedded in cofocal gui"""

    def __init__(self):
        # Get the path to the *.ui file
        this_dir = os.path.dirname(__file__)
        ui_file = os.path.join(this_dir, "ui_snvm_optim_settings.ui")

        # Load it
        super(OptimizerSettingDialog, self).__init__()
        uic.loadUi(ui_file, self)


class SnvmSettingDialog(QtWidgets.QDialog):
    """User configurable settings for the optimizer embedded in cofocal gui"""

    def __init__(self):
        # Get the path to the *.ui file
        this_dir = os.path.dirname(__file__)
        ui_file = os.path.join(this_dir, "ui_snvm_settings.ui")

        # Load it
        super(SnvmSettingDialog, self).__init__()
        uic.loadUi(ui_file, self)


class PulsedODMRSettingDialog(QtWidgets.QDialog):
    """User configurable settings for the optimizer embedded in cofocal gui"""

    def __init__(self):
        # Get the path to the *.ui file
        this_dir = os.path.dirname(__file__)
        ui_file = os.path.join(this_dir, "ui_podmr_settings.ui")

        # Load it
        super(PulsedODMRSettingDialog, self).__init__()
        uic.loadUi(ui_file, self)

        self.podmr_clkrate.editingFinished.connect(self.check_clk_rate)

        self.spinboxes_list = [
            self.podmr_start_delay,
            self.podmr_laser_init,
            self.podmr_laser_read,
            self.podmr_init_delay,
            self.podmr_read_delay,
            self.podmr_apd_delay,
            self.podmr_apd_read,
            self.podmr_final_delay,
        ]

        for spinbox in self.spinboxes_list:
            spinbox.editingFinished.connect(self.time_spinbox_edited)

    def check_clk_rate(self):
        clk_val = self.podmr_clkrate.value()
        if clk_val % 1e6:
            clk_val = round(clk_val / 1e6)
            self.podmr_clkrate.setValue(clk_val)

        for spinbox in self.spinboxes_list:
            self.time_spinbox_edited(spinbox=spinbox)

    def time_spinbox_edited(self, spinbox=None):
        if spinbox is None:
            spinbox = self.sender()
        clkrate = self.podmr_clkrate.value()

        timeval = spinbox.value()
        if not (timeval * clkrate).is_integer():
            samples = round(timeval * clkrate)
            spinbox.setValue(samples / clkrate)


class SnvmGui(GUIBase):
    """Main Confocal Class for xy and depth scans."""

    # TODO: GENERAL TODO LIST
    # -Implement the selection of the area like the one in the original qudi confocal scanner

    # declare connectors
    snvm_logic = Connector(interface="SnvmLogic")
    optimizer_logic = Connector(interface="OptimizerLogic")

    default_meter_prefix = ConfigOption(
        "default_meter_prefix", None
    )  # assume the unit prefix of position spinbox
    startstopFreq_multiplier = ConfigOption(
        "startstopFreq_multiplier", default=1e9, missing="info"
    )
    stepFreq_multiplier = ConfigOption(
        "stepFreq_multiplier", default=1e6, missing="info"
    )
    xy_range_multiplier = ConfigOption(
        "xy_range_multiplier", default=1e-9, missing="info"
    )
    px_time_multiplier = ConfigOption(
        "px_time_multiplier", default=1e-3, missing="info"
    )
    cbar_count_multiplier = ConfigOption(
        "cbar_count_multiplier", default=1e-3, missing="info"
    )
    cbar_afm_multiplier = ConfigOption(
        "cbar_afm_multplier", default=4e3 / 10, missing="info"
    )
    mwspinbox_float_resolution = ConfigOption("mwspinbox_float_resolution", default=6)

    # signals
    sigStartOptimizer = QtCore.Signal()
    sigStartScanningSnvm = QtCore.Signal()
    sigStartScanningConf = QtCore.Signal()
    sigGoTo = QtCore.Signal(str)

    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)

    def on_activate(self):
        """Initializes all needed UI files and establishes the connectors.

        This method executes the all the inits for the differnt GUIs and passes
        the event argument from fysom to the methods.
        """

        # Getting an access to all connectors:
        self._scanning_logic = self.snvm_logic()
        self._optimizer_logic = self.optimizer_logic()
        self.initMainUI()  # initialize the main GUI
        self.initOptimizer()
        self.initSnvmSettings()
        self.initPodmrSettings()

    def on_deactivate(self):
        """Reverse steps of activation

        @return int: error code (0:OK, -1:error)
        """
        self._mainwindow.close()
        return 0

    def show(self):
        """Make main window visible and put it above all other windows."""
        # Show the Main Confocal GUI:
        self._mainwindow.show()
        self._mainwindow.activateWindow()
        self._mainwindow.raise_()

    def initMainUI(self):
        """Definition, configuration and initialisation of the confocal GUI.

        This init connects all the graphic modules, which were created in the
        *.ui file and configures the event handling between the modules.
        Moreover it sets default values.
        """
        self._mainwindow = SnvmWindow()

        # Init scanning settings
        self._init_scanning_settings()

        # Init odmr settings
        self._init_odmr_settings()

        # All our gui elements are dockable, and so there should be no "central" widget.
        self._mainwindow.centralwidget.hide()
        self._mainwindow.setDockNestingEnabled(True)

        self.photon_colormap = ColorScaleInferno()
        self.afm_cmap = BlackAndWhite()
        self._crosshair_maxrange = None

        ###################################
        # Set up the SNVM image and colorbar
        ###################################
        self.snvm_image = ScanImageItem(axisOrder="row-major")
        self.snvm_image.setLookupTable(self.photon_colormap.lut)
        self._mainwindow.multiFreqPlotView.addItem(self.snvm_image)
        self._mainwindow.multiFreqPlotView.setLabel("bottom", "X (nm)")
        self._mainwindow.multiFreqPlotView.setLabel("left", "Y (nm)")
        self._mainwindow.multiFreqPlotView.toggle_crosshair(True, movable=True)
        self._mainwindow.multiFreqPlotView.set_crosshair_size((1, 1))
        self._mainwindow.multiFreqPlotView.set_crosshair_range(self.sample_ranges)
        self._mainwindow.multiFreqPlotView.sigCrosshairDraggedPosChanged.connect(
            self.move_afm_crosshair
        )
        self._mainwindow.multiFreqPlotView.toggle_selection(True)
        self._mainwindow.multiFreqPlotView.sigMouseAreaSelected.connect(
            self.update_scanning_range_snvm
        )

        snvm_im_vb = self.get_image_viewbox(self.snvm_image)
        snvm_im_vb.setAspectLocked(True)
        snvm_im_vb.toggle_selection(True)
        snvm_im_vb.toggle_zoom_by_selection(True)

        self.multifreq_cb = InteractiveColBar(
            self.snvm_image,
            self._mainwindow.multiFreqCbarView,
            lut=self.photon_colormap.lut,
            multiplier=self.cbar_count_multiplier,
            label="Counts (kHz)",
            width=80,
        )
        ###################################
        ###################################

        ###################################
        # Set up the AFM image and colorbar
        ###################################
        self.afm_image = ScanImageItem(axisOrder="row-major")
        self._mainwindow.afmPlotView.addItem(self.afm_image)
        self._mainwindow.afmPlotView.setLabel("bottom", "X (nm)")
        self._mainwindow.afmPlotView.setLabel("left", "Y (nm)")
        self._mainwindow.afmPlotView.toggle_crosshair(True, movable=True)
        self._mainwindow.afmPlotView.set_crosshair_size((1, 1))
        self._mainwindow.afmPlotView.set_crosshair_range(self.sample_ranges)
        self._mainwindow.afmPlotView.sigCrosshairDraggedPosChanged.connect(
            self.move_multifreq_crosshair
        )
        self._mainwindow.afmPlotView.toggle_selection(True)
        self._mainwindow.afmPlotView.sigMouseAreaSelected.connect(
            self.update_scanning_range_snvm
        )

        afm_im_vb = self.get_image_viewbox(self.afm_image)
        afm_im_vb.setAspectLocked(True)
        afm_im_vb.toggle_selection(True)
        afm_im_vb.toggle_zoom_by_selection(True)

        self.afm_cb = InteractiveColBar(
            self.afm_image,
            self._mainwindow.afmCbarView,
            lut=self.afm_cmap.lut,
            multiplier=self.cbar_afm_multiplier,
            label="Height (nm)",
            width=80,
        )
        ###################################
        ###################################

        ###################################
        # Set up the confocal image and colorbar
        ###################################
        self.cfc_image = ScanImageItem(axisOrder="row-major")
        self.cfc_image.setLookupTable(self.photon_colormap.lut)
        self._mainwindow.confocalScannerView.addItem(self.cfc_image)
        self._mainwindow.confocalScannerView.setLabel("bottom", "X (nm)")
        self._mainwindow.confocalScannerView.setLabel("left", "Y (nm)")
        self._mainwindow.confocalScannerView.toggle_crosshair(True, movable=True)
        self._mainwindow.confocalScannerView.set_crosshair_size((1, 1))
        self._mainwindow.confocalScannerView.set_crosshair_range(self.tip_ranges)
        self._mainwindow.confocalScannerView.toggle_selection(True)
        self._mainwindow.confocalScannerView.sigMouseAreaSelected.connect(
            self.update_scanning_range_conf
        )

        cfc_im_vb = self.get_image_viewbox(self.cfc_image)
        cfc_im_vb.setAspectLocked(True)
        cfc_im_vb.toggle_selection(True)
        cfc_im_vb.toggle_zoom_by_selection(True)

        self.cfc_cb = InteractiveColBar(
            self.cfc_image,
            self._mainwindow.confocalCbarView,
            lut=self.photon_colormap.lut,
            multiplier=self.cbar_count_multiplier,
            label="Counts (kHz)",
            width=80,
        )
        ###################################
        ###################################

        ###################################
        # Set up the optimizer image and colorbar
        ###################################
        self.optimizer_image = ScanImageItem(axisOrder="row-major")
        self.optimizer_image.setLookupTable(self.photon_colormap.lut)
        self._mainwindow.optimizerView.addItem(self.optimizer_image)
        self._mainwindow.optimizerView.setLabel("bottom", "X (nm)")
        self._mainwindow.optimizerView.setLabel("left", "Y (nm)")

        opt_im_vb = self.get_image_viewbox(self.optimizer_image)
        opt_im_vb.setAspectLocked(True)

        self.opt_cb = InteractiveColBar(
            self.snvm_image,
            self._mainwindow.optimizerCbarView,
            lut=self.photon_colormap.lut,
            multiplier=self.cbar_count_multiplier,
            label="Counts (kHz)",
            width=80,
        )
        ###################################
        ###################################

        ###################################
        # Set up the ODMR plot
        ###################################
        self.curr_odmr_trace = pg.PlotDataItem(
            skipFiniteCheck=False, connect="finite", pen=pg.mkPen(color="w")
        )
        self.average_odmr_trace = pg.PlotDataItem(
            skipFiniteCheck=True, pen=pg.mkPen(color="r")
        )
        self._mainwindow.odmrPlotWidget.addItem(self.curr_odmr_trace)
        self._mainwindow.odmrPlotWidget.addItem(self.average_odmr_trace)
        self._mainwindow.odmrPlotWidget.setLabel("bottom", "Frequency (GHz)")
        self._mainwindow.odmrPlotWidget.setLabel("left", "Counts (GHz)")
        ###################################
        ###################################

        # Quick settings for the spinbox to view the frequency slices
        self._mainwindow.frequencySliceSelector.lineEdit().setReadOnly(True)

        self._viewIndex = 0  # Variable used to scroll through the SNVM images.Gets updated when clicking the frequency selector

        # Connect the signals
        self.sigStartOptimizer.connect(self.optimize_counts, QtCore.Qt.QueuedConnection)
        self.sigStartScanningSnvm.connect(
            self.prepare_snvm_scan, QtCore.Qt.QueuedConnection
        )
        self.sigStartScanningConf.connect(
            self.prepare_conf_scan, QtCore.Qt.QueuedConnection
        )
        self.sigGoTo.connect(self.go_to_point, QtCore.Qt.QueuedConnection)

        self._odmr_widgets["mwStart"].valueChanged.connect(self.accept_frequency_ranges)
        self._odmr_widgets["mwEnd"].valueChanged.connect(self.accept_frequency_ranges)
        self._odmr_widgets["mwStep"].valueChanged.connect(self.accept_frequency_ranges)

        self._scanning_logic.signal_scan_finished.connect(self.snvm_confocal_finished)
        self._scanning_logic.signal_odmr_trace_updated.connect(
            self.refresh_odmr_plot_pxbypx
        )
        self._scanning_logic.signal_odmr_line_acquired.connect(self.refresh_odmr_plot)
        self._scanning_logic.signal_snvm_image_updated.connect(self.refresh_snvm_image)
        self._scanning_logic.signal_snvm_image_updated.connect(self.refresh_afm_image)
        self._scanning_logic.signal_xy_image_updated.connect(
            self.refresh_confocal_image
        )
        self._scanning_logic.signal_snvm_initialized.connect(self.set_snvm_im_range)
        self._scanning_logic.signal_confocal_initialized.connect(
            self.set_confocal_im_range
        )
        self._scanning_logic.signal_moved_to_point.connect(self.go_to_finished)

        self._optimizer_logic.sigImageUpdated.connect(self.refresh_optimizer_image)
        self._optimizer_logic.sigRefocusStarted.connect(self.set_optimizer_im_range)
        self._optimizer_logic.sigRefocusFinished.connect(self._optimization_complete)

        self._mainwindow.frequencySliceSelector.stepClicked.connect(
            self.frequency_selector_clicked
        )
        self._mainwindow.sampleTraceViewSpinBox.valueChanged.connect(
            self.refresh_snvm_image
        )
        self._mainwindow.sampleTraceViewSpinBox.valueChanged.connect(
            self.refresh_afm_image
        )
        self._mainwindow.tipTraceViewSpinBox.valueChanged.connect(
            self.refresh_confocal_image
        )
        self._mainwindow.sampleXSlider.sliderMoved.connect(self.slider_move_crosshair)
        self._mainwindow.sampleYSlider.sliderMoved.connect(self.slider_move_crosshair)

        self._mainwindow.tipXSlider.sliderMoved.connect(self.slider_move_crosshair)
        self._mainwindow.tipYSlider.sliderMoved.connect(self.slider_move_crosshair)

        self._mainwindow.pxbypxodmr_plotting.stateChanged.connect(
            self.set_pxbypxodmr_plot
        )

        self.snvm_range_spinboxes = [
            [self._mainwindow.xMinRangeSnvm, self._mainwindow.xMaxRangeSnvm],
            [self._mainwindow.yMinRangeSnvm, self._mainwindow.yMaxRangeSnvm],
        ]

        self.cfc_range_spinboxes = [
            [self._mainwindow.xMinRangeConf, self._mainwindow.xMaxRangeConf],
            [self._mainwindow.yMinRangeConf, self._mainwindow.yMaxRangeConf],
        ]

        for spinbox_group, maxrange in zip(
            [self.snvm_range_spinboxes, self.cfc_range_spinboxes],
            [self.sample_ranges, self.tip_ranges],
        ):

            for ii, spinboxrow in enumerate(spinbox_group):
                for spinbox in spinboxrow:
                    spinbox.setMinimum(maxrange[ii][0])
                    spinbox.setMaximum(maxrange[ii][1])
                    spinbox.editingFinished.connect(self.scanning_ranges_edited)

        self._mainwindow.scanningSettingsTab.currentChanged.connect(
            self.scanning_tab_pressed
        )

        self._mainwindow.podmr_selected.toggled.connect(self.update_podmr_active)
        self._mainwindow.podmr_pipulse.editingFinished.connect(self.adjust_pipulse_time)
        self._mainwindow.podmr_showsettings_btn.clicked.connect(
            self.menu_podmr_settings
        )

        ##############
        # Connect the actions to their slots
        ##############
        self._mainwindow.actionStart_snvm_scan.triggered.connect(
            self.scanning_action_clicked
        )
        self._mainwindow.actionStart_conf_scan.triggered.connect(
            self.scanning_action_clicked
        )
        self._mainwindow.actionStop_scan.triggered.connect(self.stop_scanning_request)
        self._mainwindow.actionOptimize.triggered.connect(self.scanning_action_clicked)
        self._mainwindow.actionOptimizer_settings.triggered.connect(
            self.menu_optimizer_settings
        )
        self._mainwindow.actionSnvm_settings.triggered.connect(self.menu_snvm_settings)
        self._mainwindow.actionPulsing_settings.triggered.connect(
            self.menu_podmr_settings
        )
        self._mainwindow.action_snvm_goToPoint.triggered.connect(
            self.scanning_action_clicked
        )
        self._mainwindow.action_cfc_goToPoint.triggered.connect(
            self.scanning_action_clicked
        )
        self._mainwindow.actionSave_snvm.triggered.connect(self.save_snvm_data)
        self._mainwindow.actionSave_confocal.triggered.connect(self.save_confocal_data)

        self._mainwindow.actionStop_scan.setEnabled(False)

        self.snvm_interactions_enabled(enabled=False)
        self.conf_interactions_enabled()
        self.show()

    def initOptimizer(self):
        self._optim_dialog = OptimizerSettingDialog()
        # Connect the action of the settings window with the code:
        self._optim_dialog.accepted.connect(self.update_optimizer_settings)
        self._optim_dialog.rejected.connect(self.keep_former_optimizer_settings)
        self._optim_dialog.buttonBox.button(
            QtWidgets.QDialogButtonBox.Apply
        ).clicked.connect(self.update_optimizer_settings)

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
        self._snvm_dialog.buttonBox.button(
            QtWidgets.QDialogButtonBox.Apply
        ).clicked.connect(self.update_optimizer_settings)

        # write the configuration to the settings window of the GUI.
        self.keep_former_snvm_settings()

    def initPodmrSettings(self):
        self._podmr_dialog = PulsedODMRSettingDialog()
        # Connect the action of the settings window with the code:
        self._podmr_dialog.accepted.connect(self.update_podmr_settings)
        self._podmr_dialog.rejected.connect(self.keep_former_podmr_settings)

        # write the configuration to the settings window of the GUI.
        self.keep_former_podmr_settings()

    def _init_odmr_settings(self):
        # Also here, store in a dictionary if the widgets need to be accessed in for loops
        self._odmr_widgets = dict()
        self._odmr_widgets[
            self._mainwindow.mwStart.objectName()
        ] = self._mainwindow.mwStart
        self._odmr_widgets[self._mainwindow.mwEnd.objectName()] = self._mainwindow.mwEnd
        self._odmr_widgets[
            self._mainwindow.mwStep.objectName()
        ] = self._mainwindow.mwStep
        self._odmr_widgets[
            self._mainwindow.mwPower.objectName()
        ] = self._mainwindow.mwPower
        self._odmr_widgets[
            self._mainwindow.mwAverages.objectName()
        ] = self._mainwindow.mwAverages
        self._odmr_widgets[
            self._mainwindow.podmr_pipulse.objectName()
        ] = self._mainwindow.podmr_pipulse
        self._odmr_widgets[
            self._mainwindow.podmr_selected.objectName()
        ] = self._mainwindow.podmr_selected

        self._mainwindow.mwStart.setDecimals(self.mwspinbox_float_resolution)
        self._mainwindow.mwEnd.setDecimals(self.mwspinbox_float_resolution)
        self._mainwindow.mwStep.setDecimals(self.mwspinbox_float_resolution)

        #########
        self._odmr_widgets["mwStart"].setValue(
            self._scanning_logic.start_freq / self.startstopFreq_multiplier
        )
        self._odmr_widgets["mwEnd"].setValue(
            self._scanning_logic.stop_freq / self.startstopFreq_multiplier
        )
        self._odmr_widgets["mwStep"].setValue(
            self._scanning_logic.freq_resolution / self.stepFreq_multiplier
        )
        self._odmr_widgets["mwPower"].setValue(self._scanning_logic.mw_power)
        self._odmr_widgets["mwAverages"].setValue(self._scanning_logic.odmr_averages)

        podmr_selected = self._scanning_logic.podmr_active

        self._odmr_widgets["podmr_selected"].setCheckState(podmr_selected)
        self._odmr_widgets["podmr_pipulse"].setValue(self._scanning_logic.podmr_pipulse)

        self._odmr_widgets["podmr_selected"].setTristate(False)
        #########

    def _init_scanning_settings(self):
        ########
        # AFM scanning settings
        ########
        # Put all the settings in a dictionary, for ease of access
        self._conf_widgets = dict()
        # Confacal settings
        self._conf_widgets[
            self._mainwindow.xResolutionConf.objectName()
        ] = self._mainwindow.xResolutionConf
        self._conf_widgets[
            self._mainwindow.yResolutionConf.objectName()
        ] = self._mainwindow.yResolutionConf
        self._conf_widgets[
            self._mainwindow.xMinRangeConf.objectName()
        ] = self._mainwindow.xMinRangeConf
        self._conf_widgets[
            self._mainwindow.xMaxRangeConf.objectName()
        ] = self._mainwindow.xMaxRangeConf
        self._conf_widgets[
            self._mainwindow.yMinRangeConf.objectName()
        ] = self._mainwindow.yMinRangeConf
        self._conf_widgets[
            self._mainwindow.yMaxRangeConf.objectName()
        ] = self._mainwindow.yMaxRangeConf
        self._conf_widgets[
            self._mainwindow.fwpxTimeConf.objectName()
        ] = self._mainwindow.fwpxTimeConf
        self._conf_widgets[
            self._mainwindow.storeRetraceConf.objectName()
        ] = self._mainwindow.storeRetraceConf

        self._conf_widgets["yResolutionConf"].setValue(
            self._scanning_logic.scanning_y_resolution["tip"]
        )
        self._conf_widgets["xResolutionConf"].setValue(
            self._scanning_logic.scanning_x_resolution["tip"]
        )
        self._conf_widgets["xMinRangeConf"].setValue(
            self._scanning_logic.scanning_x_range["tip"][0] / self.xy_range_multiplier
        )
        self._conf_widgets["yMinRangeConf"].setValue(
            self._scanning_logic.scanning_y_range["tip"][0] / self.xy_range_multiplier
        )
        self._conf_widgets["xMaxRangeConf"].setValue(
            self._scanning_logic.scanning_x_range["tip"][1] / self.xy_range_multiplier
        )
        self._conf_widgets["yMaxRangeConf"].setValue(
            self._scanning_logic.scanning_y_range["tip"][1] / self.xy_range_multiplier
        )
        self._conf_widgets["fwpxTimeConf"].setValue(
            self._scanning_logic.px_time["tip"] / self.px_time_multiplier
        )
        self._conf_widgets["storeRetraceConf"].setChecked(
            self._scanning_logic.store_retrace["tip"]
        )

        # Snvm settings
        self._afm_widgets = dict()
        self._afm_widgets[
            self._mainwindow.xResolutionSnvm.objectName()
        ] = self._mainwindow.xResolutionSnvm
        self._afm_widgets[
            self._mainwindow.yResolutionSnvm.objectName()
        ] = self._mainwindow.yResolutionSnvm
        self._afm_widgets[
            self._mainwindow.xMinRangeSnvm.objectName()
        ] = self._mainwindow.xMinRangeSnvm
        self._afm_widgets[
            self._mainwindow.xMaxRangeSnvm.objectName()
        ] = self._mainwindow.xMaxRangeSnvm
        self._afm_widgets[
            self._mainwindow.yMinRangeSnvm.objectName()
        ] = self._mainwindow.yMinRangeSnvm
        self._afm_widgets[
            self._mainwindow.yMaxRangeSnvm.objectName()
        ] = self._mainwindow.yMaxRangeSnvm
        self._afm_widgets[
            self._mainwindow.fwpxTimeSnvm.objectName()
        ] = self._mainwindow.fwpxTimeSnvm
        self._afm_widgets[
            self._mainwindow.storeRetraceSnvm.objectName()
        ] = self._mainwindow.storeRetraceSnvm

        self._afm_widgets["xResolutionSnvm"].setValue(
            self._scanning_logic.scanning_x_resolution["sample"]
        )
        self._afm_widgets["yResolutionSnvm"].setValue(
            self._scanning_logic.scanning_y_resolution["sample"]
        )
        self._afm_widgets["xMinRangeSnvm"].setValue(
            self._scanning_logic.scanning_x_range["sample"][0]
            / self.xy_range_multiplier
        )
        self._afm_widgets["yMinRangeSnvm"].setValue(
            self._scanning_logic.scanning_y_range["sample"][0]
            / self.xy_range_multiplier
        )
        self._afm_widgets["xMaxRangeSnvm"].setValue(
            self._scanning_logic.scanning_x_range["sample"][1]
            / self.xy_range_multiplier
        )
        self._afm_widgets["yMaxRangeSnvm"].setValue(
            self._scanning_logic.scanning_y_range["sample"][1]
            / self.xy_range_multiplier
        )
        self._afm_widgets["fwpxTimeSnvm"].setValue(
            self._scanning_logic.px_time["sample"] / self.px_time_multiplier
        )
        self._afm_widgets["storeRetraceSnvm"].setChecked(
            self._scanning_logic.store_retrace["sample"]
        )

        self.sample_ranges = [
            self._scanning_logic.x_maxrange["sample"][0] / self.xy_range_multiplier,
            self._scanning_logic.x_maxrange["sample"][1] / self.xy_range_multiplier,
        ], [
            self._scanning_logic.y_maxrange["sample"][0] / self.xy_range_multiplier,
            self._scanning_logic.y_maxrange["sample"][1] / self.xy_range_multiplier,
        ]
        self.tip_ranges = [
            self._scanning_logic.x_maxrange["tip"][0] / self.xy_range_multiplier,
            self._scanning_logic.x_maxrange["tip"][1] / self.xy_range_multiplier,
        ], [
            self._scanning_logic.y_maxrange["tip"][0] / self.xy_range_multiplier,
            self._scanning_logic.y_maxrange["tip"][1] / self.xy_range_multiplier,
        ]

        self._mainwindow.sampleXSliderSpinBox.setRange(
            self.sample_ranges[0][0], self.sample_ranges[0][1]
        )
        self._mainwindow.sampleYSliderSpinBox.setRange(
            self.sample_ranges[1][0], self.sample_ranges[1][1]
        )
        self._mainwindow.tipXSliderSpinBox.setRange(
            self.tip_ranges[0][0], self.sample_ranges[0][1]
        )
        self._mainwindow.tipYSliderSpinBox.setRange(
            self.tip_ranges[1][0], self.sample_ranges[1][1]
        )

        self._sample_spinboxedsliders = dict()
        self._sample_spinboxedsliders["x"] = [
            self._mainwindow.sampleXSlider,
            self._mainwindow.sampleXSliderSpinBox,
        ]
        self._sample_spinboxedsliders["y"] = [
            self._mainwindow.sampleYSlider,
            self._mainwindow.sampleYSliderSpinBox,
        ]
        self._tip_spinboxedsliders = dict()
        self._tip_spinboxedsliders["x"] = [
            self._mainwindow.tipXSlider,
            self._mainwindow.tipXSliderSpinBox,
        ]
        self._tip_spinboxedsliders["y"] = [
            self._mainwindow.tipYSlider,
            self._mainwindow.tipYSliderSpinBox,
        ]

        # FIXME: there is a bug with the sliders: after using them, if the crosshair is grabbed it freezes the GUI.
        #  Disable the sliders and SpinBoxes for now.
        self._mainwindow.sampleXSlider.setEnabled(False)
        self._mainwindow.sampleXSliderSpinBox.setEnabled(False)
        self._mainwindow.sampleYSlider.setEnabled(False)
        self._mainwindow.sampleYSliderSpinBox.setEnabled(False)
        self._mainwindow.tipXSlider.setEnabled(False)
        self._mainwindow.tipXSliderSpinBox.setEnabled(False)
        self._mainwindow.tipYSlider.setEnabled(False)
        self._mainwindow.tipYSliderSpinBox.setEnabled(False)
        # FIXME: end of FIXME

        self.pxbypx_odmr = self._mainwindow.pxbypxodmr_plotting.isChecked()

    def initOptimizer(self):
        """Definition, configuration and initialisation of the optimizer settings GUI.

        This init connects all the graphic modules, which were created in the
        *.ui file and configures the event handling between the modules.
        Moreover it sets default values if not existed in the logic modules.
        """
        self._optim_dialog = OptimizerSettingDialog()
        # Connect the action of the settings window with the code:
        self._optim_dialog.accepted.connect(self.update_optimizer_settings)
        self._optim_dialog.rejected.connect(self.keep_former_optimizer_settings)
        self._optim_dialog.buttonBox.button(
            QtWidgets.QDialogButtonBox.Apply
        ).clicked.connect(self.update_optimizer_settings)

        # Set up and connect xy channel combobox
        stacks = self._scanning_logic.get_stack_names()
        for n, stack in enumerate(stacks):
            self._optim_dialog.optimScanner_ComboBox.addItem(stack, n)

        # write the configuration to the settings window of the GUI.
        self.keep_former_optimizer_settings()

    def prepare_snvm_scan(self):
        # Get the scanning settings from the GUI, and set them in the logic
        # FIXME: find a way to do this more efficiently, without calling each attribute one by one
        stk = self._scanning_logic.sampleStackName
        self._scanning_logic.store_retrace[stk] = (
            True if self._afm_widgets["storeRetraceSnvm"].checkState() == 2 else False
        )

        self._scanning_logic.scanning_x_range[stk] = [
            self._afm_widgets["xMinRangeSnvm"].value() * self.xy_range_multiplier,
            self._afm_widgets["xMaxRangeSnvm"].value() * self.xy_range_multiplier,
        ]
        self._scanning_logic.scanning_y_range[stk] = [
            self._afm_widgets["yMinRangeSnvm"].value() * self.xy_range_multiplier,
            self._afm_widgets["yMaxRangeSnvm"].value() * self.xy_range_multiplier,
        ]

        self._scanning_logic.scanning_x_resolution[stk] = self._afm_widgets[
            "xResolutionSnvm"
        ].value()
        self._scanning_logic.scanning_y_resolution[stk] = self._afm_widgets[
            "yResolutionSnvm"
        ].value()
        self._scanning_logic.optimize_while_scanning = (
            self._optim_dialog.optimizeDuringScanCheckBox.isChecked()
        )
        self._scanning_logic.every_N_pixels = (
            self._optim_dialog.everyNPixelsDoubleSpinBox.value()
        )

        # Set the integration time
        self._scanning_logic.px_time[stk] = (
            self._afm_widgets["fwpxTimeSnvm"].value() * self.px_time_multiplier
        )

        # First update the crosshair position
        crosshair_pos = self._mainwindow.multiFreqPlotView.crosshair_position
        if (
            crosshair_pos[0] * self.xy_range_multiplier
            not in self._scanning_logic.scanning_x_range[stk]
            or crosshair_pos[1] * self.xy_range_multiplier
            not in self._scanning_logic.scanning_y_range[stk]
        ):

            newpos = (
                self._scanning_logic.scanning_x_range[stk][0]
                / self.xy_range_multiplier,
                self._scanning_logic.scanning_y_range[stk][0]
                / self.xy_range_multiplier,
            )
            self._mainwindow.multiFreqPlotView.set_crosshair_pos(newpos)
            self._mainwindow.afmPlotView.set_crosshair_pos(newpos)

        self.set_odmr_settings()

        # Here put the settings for the  spin box.
        self._mainwindow.frequencySliceSelector.setMinimum(
            self._odmr_widgets["mwStart"].value()
        )
        self._mainwindow.frequencySliceSelector.setMaximum(
            self._odmr_widgets["mwEnd"].value()
        )
        step_val_ghz = (
            self._odmr_widgets["mwStep"].value()
            * self.stepFreq_multiplier
            / self.startstopFreq_multiplier
        )
        self._mainwindow.frequencySliceSelector.setSingleStep(step_val_ghz)
        self._viewIndex = 0
        self._mainwindow.frequencySliceSelector.setValue(
            self._odmr_widgets["mwStart"].value()
        )

        self._scanning_logic.start_snvm_scanning()

    def prepare_conf_scan(self):
        # Get the scanning settings from the GUI, and set them in the logic
        # FIXME: find a way to do this more efficiently, without calling each attribute one by one
        stk = self._scanning_logic.tipStackName
        self._scanning_logic.store_retrace[stk] = (
            True if self._conf_widgets["storeRetraceConf"].checkState() == 2 else False
        )

        self._scanning_logic.scanning_x_range[stk] = [
            self._conf_widgets["xMinRangeConf"].value() * self.xy_range_multiplier,
            self._conf_widgets["xMaxRangeConf"].value() * self.xy_range_multiplier,
        ]
        self._scanning_logic.scanning_y_range[stk] = [
            self._conf_widgets["yMinRangeConf"].value() * self.xy_range_multiplier,
            self._conf_widgets["yMaxRangeConf"].value() * self.xy_range_multiplier,
        ]

        self._scanning_logic.scanning_x_resolution[stk] = self._conf_widgets[
            "xResolutionConf"
        ].value()
        self._scanning_logic.scanning_y_resolution[stk] = self._conf_widgets[
            "yResolutionConf"
        ].value()
        self._scanning_logic.optimize_while_scanning = (
            self._optim_dialog.optimizeDuringScanCheckBox.isChecked()
        )

        # Set the integration time
        self._scanning_logic.px_time[stk] = (
            self._conf_widgets["fwpxTimeConf"].value() * self.px_time_multiplier
        )

        # First update the crosshair position
        crosshair_pos = self._mainwindow.confocalScannerView.crosshair_position
        if (
            crosshair_pos[0] not in self._scanning_logic.scanning_x_range[stk]
            or crosshair_pos[1] not in self._scanning_logic.scanning_y_range[stk]
        ):
            newpos = (
                self._scanning_logic.scanning_x_range[stk][0],
                self._scanning_logic.scanning_y_range[stk][0],
            )
            self._mainwindow.confocalScannerView.set_crosshair_pos(newpos)

        self._scanning_logic.start_confocal_scanning()

    def snvm_confocal_finished(self, was_snvm):
        self._mainwindow.actionOptimize.setEnabled(True)

        if was_snvm:
            self.snvm_interactions_enabled()
        else:
            self.conf_interactions_enabled()

    def snvm_interactions_enabled(self, enabled=True):
        self._mainwindow.actionStart_snvm_scan.setEnabled(enabled)
        self._mainwindow.actionResume_snvm_scan.setEnabled(enabled)
        self._mainwindow.action_snvm_goToPoint.setEnabled(enabled)
        for setting in self._afm_widgets.values():
            setting.setEnabled(enabled)
        for setting in self._odmr_widgets.values():
            setting.setEnabled(enabled)
        self._mainwindow.actionStop_scan.setEnabled(not enabled)

        if enabled:
            podmr_active = self._mainwindow.podmr_selected.checkState()
            self.podmr_interactions_enabled(podmr_active)
        else:
            self._mainwindow.podmr_selected.setEnabled(False)

    def conf_interactions_enabled(self, enabled=True):
        self._mainwindow.actionStart_conf_scan.setEnabled(enabled)
        self._mainwindow.actionResume_conf_scan.setEnabled(enabled)
        self._mainwindow.action_cfc_goToPoint.setEnabled(enabled)
        for setting in self._conf_widgets.values():
            setting.setEnabled(enabled)
        self._mainwindow.actionStop_scan.setEnabled(not enabled)

    def opti_interactions_enabled(self, enabled=True):
        tabval = self._mainwindow.scanningSettingsTab.currentIndex()
        if tabval == 0:
            self.snvm_interactions_enabled(enabled=enabled)
        else:
            self.conf_interactions_enabled(enabled=enabled)
        self._mainwindow.actionOptimize.setEnabled(enabled)

    def podmr_interactions_enabled(self, enabled):
        self._mainwindow.podmr_pipulse.setEnabled(enabled)
        self._mainwindow.podmr_showsettings_btn.setEnabled(enabled)

    def update_podmr_active(self, value):
        self._scanning_logic.podmr_active = value
        self.podmr_interactions_enabled(value)

    def accept_frequency_ranges(self):
        """
        Function that checks that the range is start freq + step * multiple. If it's not, update the stop frequency.
        """
        stopfreq = self._odmr_widgets["mwEnd"].value()
        startfreq = self._odmr_widgets["mwStart"].value()
        freqstep = self._odmr_widgets["mwStep"].value()

        coeff_for_step = self.stepFreq_multiplier / self.startstopFreq_multiplier

        freq_diff = stopfreq - startfreq
        freqstep = freqstep * coeff_for_step
        if not (freq_diff % freqstep) == 0:
            multiple = round(freq_diff / freqstep)
            stop_freq = startfreq + multiple * freqstep
            self._odmr_widgets["mwEnd"].setValue(stop_freq)

    def set_odmr_settings(self):
        stopfreq = self._odmr_widgets["mwEnd"].value()
        startfreq = self._odmr_widgets["mwStart"].value()
        freqstep = self._odmr_widgets["mwStep"].value()
        power = self._odmr_widgets["mwPower"].value()
        averages = self._odmr_widgets["mwAverages"].value()

        self._scanning_logic.start_freq = startfreq * self.startstopFreq_multiplier
        self._scanning_logic.stop_freq = stopfreq * self.startstopFreq_multiplier
        self._scanning_logic.freq_resolution = freqstep * self.stepFreq_multiplier
        self._scanning_logic.mw_power = power
        self._scanning_logic.odmr_averages = averages

        self._scanning_logic.podmr_active = self._mainwindow.podmr_selected.isChecked()
        self._scanning_logic.podmr_pipulse = self._mainwindow.podmr_pipulse.value()

        self._mainwindow.odmrPlotWidget.setXRange(startfreq, stopfreq)

    def stop_scanning_request(self):
        self._mainwindow.actionStop_scan.setEnabled(False)
        self._scanning_logic.stopRequested = True
        self._optimizer_logic.stop_refocus()

    def scanning_action_clicked(self):
        sendername = self.sender().objectName()
        self.opti_interactions_enabled(enabled=False)
        if sendername == "actionOptimize":
            self.sigStartOptimizer.emit()
        elif sendername == "action_snvm_goToPoint":
            self._mainwindow.actionStop_scan.setEnabled(False)
            self.sigGoTo.emit("snvm")
        elif sendername == "action_cfc_goToPoint":
            self._mainwindow.actionStop_scan.setEnabled(False)
            self.sigGoTo.emit("cfc")
        elif sendername == "actionStart_conf_scan":
            self.conf_interactions_enabled(enabled=False)
            self.sigStartScanningConf.emit()
        elif sendername == "actionStart_snvm_scan":
            self.snvm_interactions_enabled(enabled=False)
            self.sigStartScanningSnvm.emit()

    def optimize_counts(self):
        self.opti_interactions_enabled(enabled=False)

        crosshair_pos = self._mainwindow.confocalScannerView.get_crosshair_pos()
        crosshair_pos = [pos * self.xy_range_multiplier for pos in crosshair_pos]

        self._optimizer_logic.start_refocus(crosshair_pos)

    def _optimization_complete(self, coords):
        self._mainwindow.confocalScannerView.set_crosshair_pos(
            (coords[0] / self.xy_range_multiplier, coords[1] / self.xy_range_multiplier)
        )
        if not self._scanning_logic._snvm_active:
            self._scanning_logic.go_to_point(
                coords, stack=self._optimizer_logic.optimizer_stack
            )
            self.opti_interactions_enabled()

    def refresh_snvm_image(self):
        if self._mainwindow.viewtracesample:
            curr_image = self._scanning_logic.snvm_matrix.mean(axis=-1)
            mask = self._scanning_logic.completed_pixels_matrix[0]
        else:
            curr_image = self._scanning_logic.snvm_matrix_retrace.mean(axis=-1)
            mask = self._scanning_logic.completed_pixels_matrix[1]

        curr_image = curr_image[:, :, self._viewIndex] * mask

        self.update_image_levels(curr_image, self.snvm_image, self.multifreq_cb)

    def refresh_afm_image(self):
        if self._mainwindow.viewtracesample:
            curr_image = self._scanning_logic.xy_scan_matrix[:]
            mask = self._scanning_logic.completed_pixels_matrix[0]
        else:
            curr_image = self._scanning_logic.xy_scan_matrix_retrace[:]
            mask = self._scanning_logic.completed_pixels_matrix[1]

        curr_image = curr_image * mask

        self.update_image_levels(curr_image, self.afm_image, self.afm_cb)

    def refresh_confocal_image(self):
        if self._mainwindow.viewtracetip:
            curr_image = self._scanning_logic.xy_scan_matrix[:]
            mask = self._scanning_logic.completed_pixels_matrix[0]
        else:
            curr_image = self._scanning_logic.xy_scan_matrix_retrace[:]
            mask = self._scanning_logic.completed_pixels_matrix[1]

        curr_image = curr_image * mask

        self.update_image_levels(curr_image, self.cfc_image, self.cfc_cb)

    def refresh_optimizer_image(self):
        curr_image = self._optimizer_logic.xy_refocus_image[:, :, 2]

        opt_image = curr_image
        self.update_image_levels(opt_image, self.optimizer_image, self.opt_cb)

    def refresh_odmr_plot_pxbypx(self, odmr_rep_index=None):
        curr_freq_matrix = self._scanning_logic.temp_freq_matrix[odmr_rep_index]
        self.curr_odmr_trace.setData(
            self._scanning_logic.freq_axis / self.startstopFreq_multiplier,
            curr_freq_matrix,
        )
        if odmr_rep_index > 0:
            self.average_odmr_trace.setData(
                self._scanning_logic.freq_axis / self.startstopFreq_multiplier,
                self._scanning_logic.average_odmr_trace,
            )
        else:
            self.average_odmr_trace.clear()

    def refresh_odmr_plot(self):
        self.average_odmr_trace.setData(
            self._scanning_logic.freq_axis / self.startstopFreq_multiplier,
            self._scanning_logic.last_odmr_trace,
        )

    def refresh_colorbar(self, cbar, cbar_range):
        cbar.refresh_colorbar(*cbar_range)

    def get_cb_range(self, image):
        imrange = [0, 1]
        image_nonzero = image[image != 0]
        if len(image_nonzero) > 0:
            imrange = [image_nonzero.min(), image_nonzero.max()]
        if imrange[0] == imrange[1]:
            imrange[1] = imrange[0] * 1.01
        return imrange

    def update_image_levels(self, image, imageitem, cbar):
        minimage, maximage = self.get_cb_range(image)
        if not cbar.regionAutoUpdate:
            image_lvs = cbar.get_levels()
        else:
            image_lvs = [minimage, maximage]

        imageitem.setImage(image, levels=image_lvs)
        cbar.set_levels(minimage, maximage)

    def update_image_levels(self, image, imageitem, cbar):
        minimage, maximage = self.get_cb_range(image)
        if not cbar.regionAutoUpdate:
            image_lvs = cbar.get_levels()
        else:
            image_lvs = [minimage, maximage]

        imageitem.setImage(image, levels=image_lvs)
        cbar.set_levels(minimage, maximage)

    def set_snvm_im_range(self):
        im_range = self._scanning_logic.get_xy_image_range(
            multiplier=1 / self.xy_range_multiplier
        )
        xmin, xmax = im_range[0]
        ymin, ymax = im_range[1]

        xpxsize, ypxsize = self._scanning_logic.get_xy_step_size(
            multiplier=1 / self.xy_range_multiplier
        )

        for image in [self.snvm_image, self.afm_image]:
            image.set_image_extent(
                (
                    (xmin - xpxsize / 2, xmax + xpxsize / 2),
                    (ymin - ypxsize / 2, ymax + ypxsize / 2),
                )
            )

    def set_confocal_im_range(self):
        im_range = self._scanning_logic.get_xy_image_range(
            multiplier=1 / self.xy_range_multiplier
        )
        xmin, xmax = im_range[0]
        ymin, ymax = im_range[1]

        xpxsize, ypxsize = self._scanning_logic.get_xy_step_size(
            multiplier=1 / self.xy_range_multiplier
        )
        self.cfc_image.set_image_extent(
            (
                (xmin - xpxsize / 2, xmax + xpxsize / 2),
                (ymin - ypxsize / 2, ymax + ypxsize / 2),
            )
        )

    def set_optimizer_im_range(self):
        xmin, xmax = self._optimizer_logic.xy_refocus_image[
            0, np.array([0, -1]), 0
        ]  # x vals in the rows of the first slice
        ymin, ymax = self._optimizer_logic.xy_refocus_image[
            np.array([0, -1]), 0, 1
        ]  # y vals in the columns of the second one

        xpxsize = (
            self._optimizer_logic.xy_refocus_image[0, 1, 0]
            - self._optimizer_logic.xy_refocus_image[0, 0, 0]
        )
        ypxsize = (
            self._optimizer_logic.xy_refocus_image[1, 0, 1]
            - self._optimizer_logic.xy_refocus_image[0, 0, 1]
        )

        xmin, xmax, ymin, ymax, xpxsize, ypxsize = [
            val / self.xy_range_multiplier
            for val in [xmin, xmax, ymin, ymax, xpxsize, ypxsize]
        ]

        # FIXME: the following is a dirty trick to set the image range correctly (otherwise set_image_extent throws an error
        #  if no image has been previously loaded): how did I do for the other images?
        self.optimizer_image.setImage(self._optimizer_logic.xy_refocus_image[..., -1])

        self.optimizer_image.set_image_extent(
            (
                (xmin - xpxsize / 2, xmax + xpxsize / 2),
                (ymin - ypxsize / 2, ymax + ypxsize / 2),
            )
        )

    def set_pxbypxodmr_plot(self):
        self.pxbypx_odmr = self._mainwindow.pxbypxodmr_plotting.isChecked()
        self._scanning_logic.pxbypx_odmr = self.pxbypx_odmr
        self.curr_odmr_trace.clear()

    def update_optimizer_settings(self):
        self._optimizer_logic.refocus_XY_size = (
            self._optim_dialog.xy_optimizer_range_DoubleSpinBox.value() * 1e-6
        )
        self._optimizer_logic.optimizer_XY_res = (
            self._optim_dialog.xy_optimizer_resolution_SpinBox.value()
        )
        self._optimizer_logic.integration_time = (
            self._optim_dialog.intTime_SpinBox.value()
        )
        value = self._optim_dialog.optimScanner_ComboBox.currentText()
        self._optimizer_logic.optimizer_stack = value
        self._scanning_logic.optimize_while_scanning = (
            self._optim_dialog.optimizeDuringScanCheckBox.isChecked()
        )
        self._scanning_logic.every_N_pixels = (
            self._optim_dialog.everyNPixelsDoubleSpinBox.value()
        )

    def update_snvm_settings(self):
        self._scanning_logic.set_motion_speed(
            self._snvm_dialog.slowSpeedConfSpinBox.value(), stack="tip"
        )
        self._scanning_logic.set_motion_speed(
            self._snvm_dialog.slowSpeedSnvmSpinBox.value(), stack="sample"
        )
        self._scanning_logic.set_slowmotion_clockrate(
            self._snvm_dialog.motionClockRate_Spinbox.value()
        )

    def update_podmr_settings(self):
        self._scanning_logic.podmr_start_delay = (
            self._podmr_dialog.podmr_start_delay.value()
        )
        self._scanning_logic.podmr_laser_init = (
            self._podmr_dialog.podmr_laser_init.value()
        )
        self._scanning_logic.podmr_laser_read = (
            self._podmr_dialog.podmr_laser_read.value()
        )
        self._scanning_logic.podmr_init_delay = (
            self._podmr_dialog.podmr_init_delay.value()
        )
        self._scanning_logic.podmr_read_delay = (
            self._podmr_dialog.podmr_read_delay.value()
        )
        self._scanning_logic.podmr_apd_delay = (
            self._podmr_dialog.podmr_apd_delay.value()
        )
        self._scanning_logic.podmr_apd_read = self._podmr_dialog.podmr_apd_read.value()
        self._scanning_logic.podmr_final_delay = (
            self._podmr_dialog.podmr_final_delay.value()
        )
        self._scanning_logic.podmr_clkrate = self._podmr_dialog.podmr_clkrate.value()

    def update_scanning_range_snvm(self, rect):
        x, y, width, height = rect.getRect()
        x_end = x + width
        y_end = y + height
        if x > x_end:
            x, x_end = x_end, x
        if y > y_end:
            y, y_end = y_end, y

        self._afm_widgets["xMinRangeSnvm"].setValue(x)
        self._afm_widgets["yMinRangeSnvm"].setValue(y)
        self._afm_widgets["xMaxRangeSnvm"].setValue(x_end)
        self._afm_widgets["yMaxRangeSnvm"].setValue(y_end)

    def update_scanning_range_conf(self, rect):
        x, y, width, height = rect.getRect()
        x_end = x + width
        y_end = y + height
        if x > x_end:
            x, x_end = x_end, x
        if y > y_end:
            y, y_end = y_end, y

        self._afm_widgets["xMinRangeConf"].setValue(x)
        self._afm_widgets["yMinRangeConf"].setValue(y)
        self._afm_widgets["xMaxRangeConf"].setValue(x_end)
        self._afm_widgets["yMaxRangeConf"].setValue(y_end)

    def keep_former_optimizer_settings(self):
        self._optim_dialog.xy_optimizer_range_DoubleSpinBox.setValue(
            self._optimizer_logic.refocus_XY_size * 1e6
        )
        self._optim_dialog.xy_optimizer_resolution_SpinBox.setValue(
            self._optimizer_logic.optimizer_XY_res
        )
        self._optim_dialog.intTime_SpinBox.setValue(
            self._optimizer_logic.integration_time
        )

        opt_stack = self._optimizer_logic.optimizer_stack
        index = self._optim_dialog.optimScanner_ComboBox.findText(opt_stack)
        self._optim_dialog.optimScanner_ComboBox.setCurrentIndex(index)

        self._optim_dialog.optimizeDuringScanCheckBox.setChecked(
            self._scanning_logic.optimize_while_scanning
        )
        self._optim_dialog.everyNPixelsDoubleSpinBox.setValue(
            self._scanning_logic.every_N_pixels
        )

    def keep_former_podmr_settings(self):
        self._podmr_dialog.podmr_start_delay.setValue(
            self._scanning_logic.podmr_start_delay
        )
        self._podmr_dialog.podmr_laser_init.setValue(
            self._scanning_logic.podmr_laser_init
        )
        self._podmr_dialog.podmr_laser_read.setValue(
            self._scanning_logic.podmr_laser_read
        )
        self._podmr_dialog.podmr_init_delay.setValue(
            self._scanning_logic.podmr_init_delay
        )
        self._podmr_dialog.podmr_read_delay.setValue(
            self._scanning_logic.podmr_read_delay
        )
        self._podmr_dialog.podmr_apd_delay.setValue(
            self._scanning_logic.podmr_apd_delay
        )
        self._podmr_dialog.podmr_apd_read.setValue(self._scanning_logic.podmr_apd_read)
        self._podmr_dialog.podmr_final_delay.setValue(
            self._scanning_logic.podmr_final_delay
        )
        self._podmr_dialog.podmr_clkrate.setValue(self._scanning_logic.podmr_clkrate)

    def keep_former_snvm_settings(self):
        self._snvm_dialog.slowSpeedConfSpinBox.setValue(
            self._scanning_logic.backward_speed_conf
        )
        self._snvm_dialog.slowSpeedSnvmSpinBox.setValue(
            self._scanning_logic.backward_speed_snvm
        )
        self._snvm_dialog.motionClockRate_Spinbox.setValue(
            self._scanning_logic.get_slowmotion_clockrate()
        )

    def frequency_selector_clicked(self, freq_val):
        difference = (
            (freq_val - self._odmr_widgets["mwStart"].value())
            * self.startstopFreq_multiplier
            / self.stepFreq_multiplier
        )
        index = round(difference / self._odmr_widgets["mwStep"].value())

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
        """This method opens the settings menu."""
        self.keep_former_optimizer_settings()
        self._optim_dialog.exec_()

    def menu_snvm_settings(self):
        """This method opens the settings menu."""
        self.keep_former_snvm_settings()
        self._snvm_dialog.exec_()

    def menu_podmr_settings(self):
        """This method opens the settings menu."""
        self.keep_former_podmr_settings()
        self._podmr_dialog.exec_()

    def adjust_pipulse_time(self):
        clkrate = self._scanning_logic.podmr_clkrate
        timeval = self._mainwindow.podmr_pipulse.value()
        samples = timeval * clkrate
        if not samples.is_integer():
            samples = round(samples)
            self._mainwindow.podmr_pipulse.setValue(samples / clkrate)

    def slider_move_crosshair(self):
        sender = self.sender()
        sendername = sender.objectName()

        if sendername[:3] == "tip":
            range = self.tip_ranges
            sliders = self._tip_spinboxedsliders
        else:
            range = self.sample_ranges
            sliders = self._sample_spinboxedsliders

        xslider, xspinbox = sliders["x"]
        yslider, yspinbox = sliders["y"]

        x_coeff = range[0][1] - range[0][0]
        y_coeff = range[1][1] - range[1][0]

        pos_slider_x, pos_slider_y = (
            xslider.value(),
            yslider.value(),
        )

        new_x = pos_slider_x * x_coeff / xslider.maximum()
        new_y = pos_slider_y * y_coeff / yslider.maximum()

        xspinbox.setValue(new_x)
        yspinbox.setValue(new_y)

        if sendername[:3] == "tip":
            self._mainwindow.confocalScannerView.set_crosshair_pos((new_x, new_y))
        else:
            self._mainwindow.multiFreqPlotView.set_crosshair_pos((new_x, new_y))
            self._mainwindow.afmPlotView.set_crosshair_pos((new_x, new_y))

    def scanning_tab_pressed(self, tab_index):
        if not self._mainwindow.actionStop_scan.isEnabled():
            if tab_index == 0:
                self.conf_interactions_enabled(enabled=False)
                self.snvm_interactions_enabled()
            if tab_index == 1:
                self.snvm_interactions_enabled(enabled=False)
                self.conf_interactions_enabled()

    def scanning_ranges_edited(self):
        sender = self.sender()
        sendername = sender.objectName()

        if sendername[-4:] == "Snvm":
            spinbox_group = self.snvm_range_spinboxes
        else:
            spinbox_group = self.cfc_range_spinboxes

        if sendername[0] == "x":
            row_idx = 0
        else:
            row_idx = 1

        minbox_val, maxbox_val = (
            spinbox_group[row_idx][0].value(),
            spinbox_group[row_idx][1].value(),
        )
        if minbox_val > maxbox_val:
            spinbox_group[row_idx][0].setValue(maxbox_val)
            spinbox_group[row_idx][1].setValue(minbox_val)

    def go_to_point(self, scanner):
        if scanner == "snvm":
            position = self._mainwindow.multiFreqPlotView.crosshair_position
            position = [pos * self.xy_range_multiplier for pos in position]
            self._scanning_logic.go_to_point(
                position, stack=self._scanning_logic.sampleStackName, caller="gui"
            )
        elif scanner == "cfc":
            position = self._mainwindow.confocalScannerView.crosshair_position
            position = [pos * self.xy_range_multiplier for pos in position]
            self._scanning_logic.go_to_point(
                position, stack=self._scanning_logic.tipStackName, caller="gui"
            )
        else:
            self.log.exception("Invalid name.")

    def go_to_finished(self, callertag):
        if callertag == "gui":
            self.opti_interactions_enabled()
            self._mainwindow.actionStop_scan.setEnabled(True)
        else:
            pass

    def save_snvm_data(self):
        self._scanning_logic.save_snvm()

    def save_confocal_data(self):
        self._scanning_logic.save_confocal()
