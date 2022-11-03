# -*- coding: utf-8 -*-
"""
This file contains the Qudi GUI module for ODMR control.

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
from core.util import units
from core.configoption import ConfigOption
from gui.guibase import GUIBase
from gui.guiutils import ColorBar
from gui.colordefs import ColorScaleInferno
from gui.colordefs import QudiPalettePale as palette
from gui.fitsettings import FitSettingsDialog, FitSettingsComboBox
from qtpy import QtCore
from qtpy import QtCore, QtWidgets, uic
from qtwidgets.scientific_spinbox import ScienDSpinBox
from qtpy import uic
from functools import partial


class PulsedMainWindow(QtWidgets.QMainWindow):
    """ The main window for the ODMR measurement GUI.
    """

    def __init__(self):
        # Get the path to the *.ui file
        this_dir = os.path.dirname(__file__)
        ui_file = os.path.join(this_dir, 'pulsed.ui') #takes the qt designer file for layout

        # Load it
        super(PulsedMainWindow, self).__init__()
        uic.loadUi(ui_file, self)
        self.show()


# class ODMRSettingDialog(QtWidgets.QDialog):
#     """ The settings dialog for ODMR measurements.
#     """
#
#     def __init__(self):
#         # Get the path to the *.ui file
#         this_dir = os.path.dirname(__file__)
#         ui_file = os.path.join(this_dir, 'ui_odmr_settings.ui')
#
#         # Load it
#         super(ODMRSettingDialog, self).__init__()
#         uic.loadUi(ui_file, self)


class PulsedGui(GUIBase):
    """
    This is the GUI Class for pulsed measurements
    """

    # declare connectors to master pulse logic
    master_pulselogic = Connector(interface='MasterPulse')

    # some signals which are also used in the master pulse logic
    sigStartMeasurement = QtCore.Signal() # This starts the measurement
    sigStopMesurement = QtCore.Signal() # This stops the measurement
    # sigClearData = QtCore.Signal()
    # sigClearAWQ = qtCoreSignal()
    # sigMwOff = QtCore.Signal()
    # sigMwPowerChanged = QtCore.Signal(float)
    # sigMwCwParamsChanged = QtCore.Signal(float, float)
    # sigMwSweepParamsChanged = QtCore.Signal(list, list, list, float)
    # sigClockFreqChanged = QtCore.Signal(float)
    # sigOversamplingChanged = QtCore.Signal(int)
    # sigLockInChanged = QtCore.Signal(bool)
    # sigFitChanged = QtCore.Signal(str)
    # sigNumberOfLinesChanged = QtCore.Signal(int)
    # sigRuntimeChanged = QtCore.Signal(float)
    # sigDoFit = QtCore.Signal(str, object, object, int)
    # sigSaveMeasurement = QtCore.Signal(str, list, list)
    # sigAverageLinesChanged = QtCore.Signal(int)

    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)

    def on_activate(self):
        """ Definition, configuration and initialisation of the ODMR GUI.

        This init connects all the graphic modules, which were created in the
        *.ui file and configures the event handling between the modules.
        """

        self._master_pulselogic = self.master_pulselogic()
        # Use the inherited class 'Ui_ODMRGuiUI' to create now the GUI element:
        self._mw = PulsedMainWindow()
        # self._sd = ODMRSettingDialog()
        # Create a QSettings object for the mainwindow and store the actual GUI layout
        self.mwsettings = QtCore.QSettings("QUDI", "ODMR")
        self.mwsettings.setValue("geometry", self._mw.saveGeometry())
        self.mwsettings.setValue("windowState", self._mw.saveState())

        # Adjust range of scientific spinboxes above what is possible in Qt Designer
        # This makes sense for the microwave at least?
        self._mw.mw_power.setMaximum(-5)
        self._mw.mw_power.setMinimum(-50)
        self._mw.mw_freq.setMaximum(4000e6)
        self._mw.mw_freq.setMaximum(2000e6)



    # def _setup_connections(self):
        ########################################################################
        #                       Connect signals                                #
        ########################################################################
        # # Internal user input changed signals
        # # These are the boxes of the GUI?
        # # self._mw.cw_frequency_DoubleSpinBox.editingFinished.connect(self.change_cw_params)
        # self._mw.mw_power.editingFinished.connect(self.change_mw_power)
        # self._mw.mw_freq.editingFinished.connect(self.change_mw_freq)
        # self._mw.mw_freq.editingFinished.connect(self.change_mw_freq)
        # self._mw.comboBox.valueChanged.connect(self.change_method)
        #
        # # Internal trigger signals
        # self._mw.clear_awg_PushButton.clicked.connect(self.clear_odmr_data)
        # self._mw.action_run_stop.triggered.connect(self.run_stop_odmr)
        # self._mw.action_resume_odmr.triggered.connect(self.resume_odmr)
        # self._mw.action_toggle_cw.triggered.connect(self.toggle_cw_mode)
        # self._mw.action_Save.triggered.connect(self.save_data)
        # self._mw.action_RestoreDefault.triggered.connect(self.restore_defaultview)
        # self._mw.do_fit_PushButton.clicked.connect(self.do_fit)
        # self._mw.fit_range_SpinBox.editingFinished.connect(self.update_fit_range)
        #
        # # Control/values-changed signals to logic
        # # self.sigCwMwOn.connect(self._odmr_logic.mw_cw_on, QtCore.Qt.QueuedConnection)
        # # self.sigMwOff.connect(self._odmr_logic.mw_off, QtCore.Qt.QueuedConnection)
        # self.sigClearData.connect(self._master_pulselogic.clear_all, QtCore.Qt.QueuedConnection)
        # self.sigStartOdmrScan.connect(self.start_odmr,
        #                               QtCore.Qt.QueuedConnection)
        # self.sigStopOdmrScan.connect(self._odmr_logic.stop_odmr, QtCore.Qt.QueuedConnection)
        # self.sigContinueOdmrScan.connect(self._odmr_logic.continue_odmr,
        #                                  QtCore.Qt.QueuedConnection)
        # self.sigDoFit.connect(self._odmr_logic.do_fit, QtCore.Qt.QueuedConnection)
        # # self.sigMwCwParamsChanged.connect(self._odmr_logic.set_cw_parameters,
        # #                                   QtCore.Qt.QueuedConnection)
        # # self.sigMwSweepParamsChanged.connect(self._odmr_logic.set_sweep_parameters,
        # #                                      QtCore.Qt.QueuedConnection)
        # # self.sigRuntimeChanged.connect(self._odmr_logic.set_runtime, QtCore.Qt.QueuedConnection)
        # # self.sigNumberOfLinesChanged.connect(self._odmr_logic.set_matrix_line_number,
        # #                                      QtCore.Qt.QueuedConnection)
        # # self.sigClockFreqChanged.connect(self._odmr_logic.set_clock_frequency,
        # #                                  QtCore.Qt.QueuedConnection)
        # # self.sigOversamplingChanged.connect(self._odmr_logic.set_oversampling, QtCore.Qt.QueuedConnection)
        # # self.sigLockInChanged.connect(self._odmr_logic.set_lock_in, QtCore.Qt.QueuedConnection)
        # # self.sigSaveMeasurement.connect(self._odmr_logic.save_odmr_data, QtCore.Qt.QueuedConnection)
        # # self.sigAverageLinesChanged.connect(self._odmr_logic.set_average_length,
        # #                                     QtCore.Qt.QueuedConnection)
        #
        # # Update signals coming from logic:
        # self._odmr_logic.sigFreqPxAcquired.connect(self.refresh_odmr_plot,
        #                                            QtCore.Qt.QueuedConnection)
        # self._odmr_logic.sigStopOdmr.connect(self.stop_odmr,
        #                                      QtCore.Qt.QueuedConnection)
        # self._odmr_logic.sigOdmrTraceAcquired.connect(self.update_elapsed_sweeps)
        # # self._odmr_logic.sigParameterUpdated.connect(self.update_parameter,
        # #                                              QtCore.Qt.QueuedConnection)
        # # self._odmr_logic.sigOutputStateUpdated.connect(self.update_status,
        # #                                                QtCore.Qt.QueuedConnection)
        # # self._odmr_logic.sigOdmrPlotsUpdated.connect(self.update_plots, QtCore.Qt.QueuedConnection)
        # self._odmr_logic.sigOdmrFitUpdated.connect(self.update_fit, QtCore.Qt.QueuedConnection)
        # # self._odmr_logic.sigOdmrElapsedTimeUpdated.connect(self.update_elapsedtime,
        # #                                                    QtCore.Qt.QueuedConnection)
        #
        # # connect settings signals
        # self._mw.action_Settings.triggered.connect(self._menu_settings)
        # self._sd.accepted.connect(self.update_settings)
        # self._sd.rejected.connect(self.reject_settings)
        # self._sd.buttonBox.button(QtWidgets.QDialogButtonBox.Apply).clicked.connect(
        #     self.update_settings)
        # # self.reject_settings()
        #
        # # Show the Main ODMR GUI:
        # self.show()

    def on_deactivate(self):
        """ Reverse steps of activation

        @return int: error code (0:OK, -1:error)
        """
        # Disconnect signals
        # self._sd.buttonBox.button(QtWidgets.QDialogButtonBox.Apply).clicked.disconnect()
        # self._sd.accepted.disconnect()
        # self._sd.rejected.disconnect()
        # self._mw.action_Settings.triggered.disconnect()
        # self.sigCwMwOn.disconnect()
        # self.sigMwOff.disconnect()
        # self.sigClearData.disconnect()
        # self.sigStartOdmrScan.disconnect()
        # self.sigStopOdmrScan.disconnect()
        # self.sigContinueOdmrScan.disconnect()
        # self.sigDoFit.disconnect()
        # self.sigMwCwParamsChanged.disconnect()
        # self.sigMwSweepParamsChanged.disconnect()
        # self.sigRuntimeChanged.disconnect()
        # self.sigNumberOfLinesChanged.disconnect()
        # self.sigClockFreqChanged.disconnect()
        # self.sigOversamplingChanged.disconnect()
        # self.sigLockInChanged.disconnect()
        # self.sigSaveMeasurement.disconnect()
        # self.sigAverageLinesChanged.disconnect()
        # self._mw.odmr_cb_manual_RadioButton.clicked.disconnect()
        # self._mw.odmr_cb_centiles_RadioButton.clicked.disconnect()
        # self._mw.clear_odmr_PushButton.clicked.disconnect()
        # self._mw.action_run_stop.triggered.disconnect()
        # # self._mw.action_resume_odmr.triggered.disconnect()
        # self._mw.action_Save.triggered.disconnect()
        # self._mw.action_toggle_cw.triggered.disconnect()
        # self._mw.action_toggle_cw.triggered.disconnect()
        # self._mw.action_RestoreDefault.triggered.disconnect()
        # self._mw.do_fit_PushButton.clicked.disconnect()
        # # self._mw.cw_frequency_DoubleSpinBox.editingFinished.disconnect()
        # dspinbox_dict = self.get_all_dspinboxes_from_groupbox()
        # for identifier_name in dspinbox_dict:
        #     dspinbox_type_list = dspinbox_dict[identifier_name]
        #     [dspinbox_type.editingFinished.disconnect() for dspinbox_type in dspinbox_type_list]
        #
        # # self._mw.cw_power_DoubleSpinBox.editingFinished.disconnect()
        # # self._mw.sweep_power_DoubleSpinBox.editingFinished.disconnect()
        # # self._mw.runtime_DoubleSpinBox.editingFinished.disconnect()
        # # self._mw.odmr_cb_max_DoubleSpinBox.valueChanged.disconnect()
        # # self._mw.odmr_cb_min_DoubleSpinBox.valueChanged.disconnect()
        # # self._mw.odmr_cb_high_percentile_DoubleSpinBox.valueChanged.disconnect()
        # # self._mw.odmr_cb_low_percentile_DoubleSpinBox.valueChanged.disconnect()
        # self._mw.average_level_SpinBox.valueChanged.disconnect()
        # self._fsd.sigFitsUpdated.disconnect()
        # self._mw.fit_range_SpinBox.editingFinished.disconnect()
        # self._mw.action_FitSettings.triggered.disconnect()
        self._mw.close()
        return 0

    # def _set_enabled_odmr_ui(self, val):
        # """Set the enabled/disabled state for the odmr setting ui"""
        # self._mw.action_toggle_cw.setEnabled(val)
        # self._mw.cw_power_DoubleSpinBox.setEnabled(val)
        # self._mw.sweep_power_DoubleSpinBox.setEnabled(val)
        # self._mw.cw_frequency_DoubleSpinBox.setEnabled(val)
        # self._mw.average_level_SpinBox.setEnabled(val)
        # self._mw.integration_time_doubleSpinBox.setEnabled(val)
        # self._mw.do_fit_PushButton.setEnabled(val)
        # dspinbox_dict = self.get_all_dspinboxes_from_groupbox()
        # for identifier_name in dspinbox_dict:
        #     dspinbox_type_list = dspinbox_dict[identifier_name]
        #     [dspinbox_type.setEnabled(val) for dspinbox_type in dspinbox_type_list]
        # self._mw.odmr_control_DockWidget.add_range_button.setEnabled(val)
        # self._mw.odmr_control_DockWidget.remove_range_button.setEnabled(val)
        # self._mw.action_Save.setEnabled(val)
        # self._sd.clock_frequency_DoubleSpinBox.setEnabled(val)
        # self._sd.oversampling_SpinBox.setEnabled(val)
        # self._sd.lock_in_CheckBox.setEnabled(val)

    # def run_stop_measurement(self, is_checked):
    #     """ Manages what happens if measurement is started/stopped. """
        # if is_checked:
        #     # change the axes appearance according to input values:
        #     # self._mw.action_run_stop.setEnabled(False)
        #     # self._mw.action_resume_odmr.setEnabled(False)
        #     self._mw.odmr_PlotWidget.removeItem(self.odmr_fit_image)
        #
        #     # Disable the ui
        #     self._set_enabled_odmr_ui(False)
        #
        #     # Reset the sweeps counter
        #     self._mw.elapsed_sweeps_DisplayWidget.display(0)
        #
        #     self.sigStartOdmrScan.emit()
        # else:
        #     # self._mw.action_run_stop.setEnabled(False)
        #     # self._mw.action_resume_odmr.setEnabled(False)
        #     # self._mw.action_toggle_cw.setEnabled(False)
        #
        #     # Enable the ui
        #     self._set_enabled_odmr_ui(True)
        #
        #     self._odmr_logic.stopRequested = True
        #     self.sigStopOdmrScan.emit()
        # return

    def start_measurement(self): # This get enabled with the start button?
        # Grab all the parameters from the GUI
        # general tab
        mw_power = self._mw.mw_power.value()
        mw_freq = self._mw.mw_freq.value()
        method = self._mw.combobox_method.currentText()
        averages = self._mw.averages.value()
        int_time = self._mw.integration_time.value()
        seq_len = self._mw.seq_len.value()
        sampling_rate_awg = self._mw.clk_awg.value()

        apd_start = self._mw.apd_start_time.value()
        apd_len = self._mw.apd_len.value()

        apd_ref_start = self._mw.apd_ref_start_time.value()
        apd_ref_len = self._mw.apd_ref_len.value()

        # laser times tab
        laser_in = self._mw.laser_in.value()
        laser_off = self._mw.laser_off.value()
        laser_re = self._mw.laser_re.value()

        # delay_sweep tab
        apd_len_delay = self._mw.apd_len_pulse_delay.value()
        apd_min_start_delay = self._mw.apd_min_start_delay.value()
        apd_max_start_delay = self._mw.apd_max_start_delay.value()
        apd_steps_delay = self._mw.apd_steps_delay.value()

        use_mw_delay = self._mw.use_mw_delay.value()
        mw_start_delay = self._mw.mw_start_time_delay.value()
        mw_len_delay = self._mw.mw_len_delay.value()

        # rabi tab
        mw_start_time_rabi = self._mw.mw_start_time_rabi.value()
        mw_min_len_rabi = self._mw.mw_min_len_rabi.value()
        mw_max_len_rabi = self._mw.mw_max_len_rabi.value()
        mw_steps_rabi = self._mw.mw_steps_rabi.value()

        # ramsey tab
        mw_start_time_ramsey = self._mw.mw_start_time_ramsey.value()
        mw_min_len_ramsey = self._mw.mw_min_len_ramsey.value()
        mw_max_len_ramsey = self._mw.mw_max_len_ramsey.value()
        mw_steps_ramsey = self._mw.mw_steps_ramsey.value()

        # now send all the parameters from above to the masterlogic

        self._master_pulselogic.laser_times = [laser_in, laser_off, laser_re]
        self._master_pulselogic.mw_times_ramsey = [mw_start_time_ramsey, mw_min_len_ramsey, mw_max_len_ramsey, mw_steps_ramsey]
        self._master_pulselogic.mw_times_rabi = [mw_start_time_rabi, mw_min_len_rabi, mw_max_len_rabi, mw_steps_rabi]
        self._master_pulselogic.apd_times = [apd_start, apd_len]
        self._master_pulselogic.apd_ref_times = [apd_ref_start, apd_ref_len]
        self._master_pulselogic.apd_times_sweep = [apd_len_delay, apd_min_start_delay, apd_max_start_delay, apd_steps_delay]
        self._master_pulselogic.averages = averages
        self._master_pulselogic.mw_power = mw_power
        self._master_pulselogic.mw_frequency = mw_freq
        self._master_pulselogic.integration_time = int_time
        self._master_pulselogic.method = str(method)
        self._master_pulselogic.seq_len = seq_len
        self._master_pulselogic.clk_rate_awg = sampling_rate_awg
        self._master_pulselogic.mw_pulse_setting = use_mw_delay
        self._master_pulselogic.mw_times_sweep = [mw_start_delay, mw_len_delay]

        # # Setup the plot? later
        # # The connect parameter controls which points of the plot are drawn connec
        # connect = self._calculate_connect(starts, stops, steps)
        # self.curr_odmr_trace.opts['connect'] = connect
        # self.curr_odmr_trace.setProperty('styleWasChanged', True)
        # self.average_odmr_trace.opts['connect'] = connect
        # self.average_odmr_trace.setProperty('styleWasChanged', True)


        # self._mw.odmr_PlotWidget.setXRange(min(starts), max(stops))
        #This is where the measurement really starts
        self._master_pulselogic.start_measurement()

    def stop_measurement(self): # This gets enabled by the stop button?
        # self._set_enabled_pulsed_ui(True)
        # self._mw.action_run_stop.setChecked(False)
        self._master_pulselogic.stop_all()

    # def _calculate_connect(self, starts, stops, steps):
    #     if len(starts) == 1:
    #         return 'all'
    #
    #     hide_indices = []
    #     length = 0
    #     for start, step, stop in zip(starts, steps, stops):
    #         step_num = 1 + round((stop - start) / step)
    #         length += step_num
    #         hide_indices.append(length - 1)
    #
    #     connect = np.ones(length)
    #     connect[hide_indices[:-1]] = 0
    #     return connect
    #

    # def refresh_odmr_plot(self):
    #     # Draw current odmr trace
    #     self.curr_odmr_trace.setData(self._odmr_logic.freq_axis,
    #                                  self._odmr_logic.curr_odmr_trace)
    #
    #     # Draw average odmr trace
    #     if self._odmr_logic._average_index > 0:
    #         self.average_odmr_trace.setData(self._odmr_logic.freq_axis,
    #                                         self._odmr_logic.average_odmr_trace)
    #     else:
    #         self.average_odmr_trace.clear()

    # def show(self):
    #     """Make window visible and put it above all other windows. """
    #     self._mw.show()
    #     self._mw.activateWindow()
    #     self._mw.raise_()

    # def resume_odmr(self, is_checked):
    #     if is_checked:
    #         self._mw.action_run_stop.setEnabled(False)
    #         self._mw.action_resume_odmr.setEnabled(False)
    #         self._mw.action_toggle_cw.setEnabled(False)
    #         self._mw.cw_power_DoubleSpinBox.setEnabled(False)
    #         self._mw.sweep_power_DoubleSpinBox.setEnabled(False)
    #         self._mw.cw_frequency_DoubleSpinBox.setEnabled(False)
    #         dspinbox_dict = self.get_all_dspinboxes_from_groupbox()
    #         for identifier_name in dspinbox_dict:
    #             dspinbox_type_list = dspinbox_dict[identifier_name]
    #             [dspinbox_type.setEnabled(False) for dspinbox_type in dspinbox_type_list]
    #         self._mw.odmr_control_DockWidget.add_range_button.setEnabled(False)
    #         self._mw.odmr_control_DockWidget.remove_range_button.setEnabled(False)
    #         self._mw.runtime_DoubleSpinBox.setEnabled(False)
    #         self._sd.clock_frequency_DoubleSpinBox.setEnabled(False)
    #         self._sd.oversampling_SpinBox.setEnabled(False)
    #         self._sd.lock_in_CheckBox.setEnabled(False)
    #         self.sigContinueOdmrScan.emit()
    #     else:
    #         self._mw.action_run_stop.setEnabled(False)
    #         self._mw.action_resume_odmr.setEnabled(False)
    #         self._mw.action_toggle_cw.setEnabled(False)
    #         self.sigStopOdmrScan.emit()
    #     return

    # def _menu_settings(self):
    #     """ Open the settings menu """
    #     self._sd.exec_()

    # def add_ranges_gui_elements_clicked(self):
    #     """
    #     When button >>add range<< is pushed add some buttons to the gui and connect accordingly to the
    #     logic.
    #     :return:
    #     """
    #     # make sure the logic keeps track
    #     groupBox = self._mw.odmr_control_DockWidget.ranges_groupBox
    #     gridLayout = groupBox.layout()
    #     constraints = self._odmr_logic.get_hw_constraints()
    #
    #     insertion_row = self._odmr_logic.ranges
    #     # start
    #     start_label = QtWidgets.QLabel(groupBox)
    #     start_label.setText('Start:')
    #     setattr(self._mw.odmr_control_DockWidget, 'start_label_{}'.format(insertion_row), start_label)
    #     start_freq_DoubleSpinBox = ScienDSpinBox(groupBox)
    #     start_freq_DoubleSpinBox.setSuffix('Hz')
    #     start_freq_DoubleSpinBox.setMaximum(constraints.max_frequency)
    #     start_freq_DoubleSpinBox.setMinimum(constraints.min_frequency)
    #     start_freq_DoubleSpinBox.setMinimumSize(QtCore.QSize(80, 0))
    #     start_freq_DoubleSpinBox.setValue(self._odmr_logic.mw_starts[0])
    #     start_freq_DoubleSpinBox.setMinimumWidth(75)
    #     start_freq_DoubleSpinBox.setMaximumWidth(100)
    #     start_freq_DoubleSpinBox.editingFinished.connect(self.change_sweep_params)
    #     setattr(self._mw.odmr_control_DockWidget, 'start_freq_DoubleSpinBox_{}'.format(insertion_row),
    #             start_freq_DoubleSpinBox)
    #     gridLayout.addWidget(start_label, insertion_row, 1, 1, 1)
    #     gridLayout.addWidget(start_freq_DoubleSpinBox, insertion_row, 2, 1, 1)
    #
    #     # step
    #     step_label = QtWidgets.QLabel(groupBox)
    #     step_label.setText('Step:')
    #     setattr(self._mw.odmr_control_DockWidget, 'step_label_{}'.format(insertion_row), step_label)
    #     step_freq_DoubleSpinBox = ScienDSpinBox(groupBox)
    #     step_freq_DoubleSpinBox.setSuffix('Hz')
    #     step_freq_DoubleSpinBox.setMaximum(100e9)
    #     step_freq_DoubleSpinBox.setMinimumSize(QtCore.QSize(80, 0))
    #     step_freq_DoubleSpinBox.setValue(self._odmr_logic.mw_steps[0])
    #     step_freq_DoubleSpinBox.setMinimumWidth(75)
    #     step_freq_DoubleSpinBox.setMaximumWidth(100)
    #     step_freq_DoubleSpinBox.editingFinished.connect(self.change_sweep_params)
    #     setattr(self._mw.odmr_control_DockWidget, 'step_freq_DoubleSpinBox_{}'.format(insertion_row),
    #             step_freq_DoubleSpinBox)
    #     gridLayout.addWidget(step_label, insertion_row, 3, 1, 1)
    #     gridLayout.addWidget(step_freq_DoubleSpinBox, insertion_row, 4, 1, 1)
    #
    #     # stop
    #     stop_label = QtWidgets.QLabel(groupBox)
    #     stop_label.setText('Stop:')
    #     setattr(self._mw.odmr_control_DockWidget, 'stop_label_{}'.format(insertion_row), stop_label)
    #     stop_freq_DoubleSpinBox = ScienDSpinBox(groupBox)
    #     stop_freq_DoubleSpinBox.setSuffix('Hz')
    #     stop_freq_DoubleSpinBox.setMaximum(constraints.max_frequency)
    #     stop_freq_DoubleSpinBox.setMinimum(constraints.min_frequency)
    #     stop_freq_DoubleSpinBox.setMinimumSize(QtCore.QSize(80, 0))
    #     stop_freq_DoubleSpinBox.setValue(self._odmr_logic.mw_stops[0])
    #     stop_freq_DoubleSpinBox.setMinimumWidth(75)
    #     stop_freq_DoubleSpinBox.setMaximumWidth(100)
    #     stop_freq_DoubleSpinBox.editingFinished.connect(self.change_sweep_params)
    #     setattr(self._mw.odmr_control_DockWidget, 'stop_freq_DoubleSpinBox_{}'.format(insertion_row),
    #             stop_freq_DoubleSpinBox)
    #
    #     gridLayout.addWidget(stop_label, insertion_row, 5, 1, 1)
    #     gridLayout.addWidget(stop_freq_DoubleSpinBox, insertion_row, 6, 1, 1)
    #
    #     starts = self.get_frequencies_from_spinboxes('start')
    #     stops = self.get_frequencies_from_spinboxes('stop')
    #     steps = self.get_frequencies_from_spinboxes('step')
    #     power = self._mw.sweep_power_DoubleSpinBox.value()
    #
    #     self.sigMwSweepParamsChanged.emit(starts, stops, steps, power)
    #     self._mw.fit_range_SpinBox.setMaximum(self._odmr_logic.ranges)
    #     self._odmr_logic.ranges += 1

        # # remove stuff that remained from the old range that might have been in place there
        # key = 'channel: {0}, range: {1}'.format(self.display_channel, self._odmr_logic.ranges - 1)
        # if key in self._odmr_logic.fits_performed:
        #     self._odmr_logic.fits_performed.pop(key)
        # return

    # def remove_ranges_gui_elements_clicked(self):
    #     if self._odmr_logic.ranges == 1:
    #         return
    #
    #     remove_row = self._odmr_logic.ranges - 1
    #
    #     groupBox = self._mw.odmr_control_DockWidget.ranges_groupBox
    #     gridLayout = groupBox.layout()
    #
    #     object_dict = self.get_objects_from_groupbox_row(remove_row)
    #
    #     for object_name in object_dict:
    #         if 'DoubleSpinBox' in object_name:
    #             object_dict[object_name].editingFinished.disconnect()
    #         object_dict[object_name].hide()
    #         gridLayout.removeWidget(object_dict[object_name])
    #         del self._mw.odmr_control_DockWidget.__dict__[object_name]
    #
    #     starts = self.get_frequencies_from_spinboxes('start')
    #     stops = self.get_frequencies_from_spinboxes('stop')
    #     steps = self.get_frequencies_from_spinboxes('step')
    #     power = self._mw.sweep_power_DoubleSpinBox.value()
    #     self.sigMwSweepParamsChanged.emit(starts, stops, steps, power)
    #
    #     # in case the removed range is the one selected for fitting right now adjust the value
    #     self._odmr_logic.ranges -= 1
        # max_val = self._odmr_logic.ranges - 1
        # self._mw.fit_range_SpinBox.setMaximum(max_val)
        # if self._odmr_logic.range_to_fit > max_val:
        #     self._odmr_logic.range_to_fit = max_val
        #
        # self._mw.fit_range_SpinBox.setMaximum(max_val)

    # def get_objects_from_groupbox_row(self, row):
        # get elements from the row
        # first strings
        # start_label_str = 'start_label_{}'.format(row)
        # step_label_str = 'step_label_{}'.format(row)
        # stop_label_str = 'stop_label_{}'.format(row)
        #
        # # get widgets
        # start_freq_DoubleSpinBox_str = 'start_freq_DoubleSpinBox_{}'.format(row)
        # step_freq_DoubleSpinBox_str = 'step_freq_DoubleSpinBox_{}'.format(row)
        # stop_freq_DoubleSpinBox_str = 'stop_freq_DoubleSpinBox_{}'.format(row)
        #
        # # now get the objects
        # start_label = getattr(self._mw.odmr_control_DockWidget, start_label_str)
        # step_label = getattr(self._mw.odmr_control_DockWidget, step_label_str)
        # stop_label = getattr(self._mw.odmr_control_DockWidget, stop_label_str)
        #
        # start_freq_DoubleSpinBox = getattr(self._mw.odmr_control_DockWidget, start_freq_DoubleSpinBox_str)
        # step_freq_DoubleSpinBox = getattr(self._mw.odmr_control_DockWidget, step_freq_DoubleSpinBox_str)
        # stop_freq_DoubleSpinBox = getattr(self._mw.odmr_control_DockWidget, stop_freq_DoubleSpinBox_str)
        #
        # return_dict = {start_label_str: start_label, step_label_str: step_label,
        #                stop_label_str: stop_label,
        #                start_freq_DoubleSpinBox_str: start_freq_DoubleSpinBox,
        #                step_freq_DoubleSpinBox_str: step_freq_DoubleSpinBox,
        #                stop_freq_DoubleSpinBox_str: stop_freq_DoubleSpinBox
        #                }
        #
        # return return_dict

    # def get_freq_dspinboxes_from_groubpox(self, identifier):
    #     dspinboxes = []
    #     for name in self._mw.odmr_control_DockWidget.__dict__:
    #         box_name = identifier + '_freq_DoubleSpinBox'
    #         if box_name in name:
    #             freq_DoubleSpinBox = getattr(self._mw.odmr_control_DockWidget, name)
    #             dspinboxes.append(freq_DoubleSpinBox)
    #
    #     return dspinboxes

    # def get_all_dspinboxes_from_groupbox(self):
    #     identifiers = ['start', 'step', 'stop']
    #
    #     all_spinboxes = {}
    #     for identifier in identifiers:
    #         all_spinboxes[identifier] = self.get_freq_dspinboxes_from_groubpox(identifier)
    #
    #     return all_spinboxes

    # def get_frequencies_from_spinboxes(self, identifier):
    #     dspinboxes = self.get_freq_dspinboxes_from_groubpox(identifier)
    #     freqs = [dspinbox.value() for dspinbox in dspinboxes]
    #     return freqs

    def toggle_cw_mode(self, is_checked): # This is the button on the right
        """ Starts or stops CW microwave output if no measurement is running. """
        if is_checked:
            pass
            # self._mw.action_run_stop.setEnabled(False)
            # self._mw.action_resume_odmr.setEnabled(False)
            # self._mw.action_toggle_cw.setEnabled(False)
            # self._mw.cw_power_DoubleSpinBox.setEnabled(False)
            # self._mw.cw_frequency_DoubleSpinBox.setEnabled(False)
            # self.sigCwMwOn.emit()
        else:
            # self._mw.action_toggle_cw.setEnabled(False)
            # self.sigMwOff.emit()
            pass
        return

    # def update_status(self, mw_mode, is_running):
    #     """
    #     Update the display for a change in the microwave status (mode and output).
    #
    #     @param str mw_mode: is the microwave output active?
    #     @param bool is_running: is the microwave output active?
    #     """
    #     # Block signals from firing
    #     self._mw.action_run_stop.blockSignals(True)
    #     self._mw.action_resume_odmr.blockSignals(True)
    #     self._mw.action_toggle_cw.blockSignals(True)
    #
    #     # Update measurement status (activate/deactivate widgets/actions)
    #     if is_running:
    #         self._mw.action_resume_odmr.setEnabled(False)
    #         self._mw.cw_power_DoubleSpinBox.setEnabled(False)
    #         self._mw.cw_frequency_DoubleSpinBox.setEnabled(False)
    #         if mw_mode != 'cw':
    #             self._mw.clear_odmr_PushButton.setEnabled(True)
    #             self._mw.action_run_stop.setEnabled(True)
    #             self._mw.action_toggle_cw.setEnabled(False)
    #             dspinbox_dict = self.get_all_dspinboxes_from_groupbox()
    #             for identifier_name in dspinbox_dict:
    #                 dspinbox_type_list = dspinbox_dict[identifier_name]
    #                 [dspinbox_type.setEnabled(False) for dspinbox_type in dspinbox_type_list]
    #             self._mw.odmr_control_DockWidget.add_range_button.setEnabled(False)
    #             self._mw.odmr_control_DockWidget.remove_range_button.setEnabled(False)
    #             self._mw.sweep_power_DoubleSpinBox.setEnabled(False)
    #             self._mw.runtime_DoubleSpinBox.setEnabled(False)
    #             self._sd.clock_frequency_DoubleSpinBox.setEnabled(False)
    #             self._sd.oversampling_SpinBox.setEnabled(False)
    #             self._sd.lock_in_CheckBox.setEnabled(False)
    #             self._mw.action_run_stop.setChecked(True)
    #             self._mw.action_resume_odmr.setChecked(True)
    #             self._mw.action_toggle_cw.setChecked(False)
    #         else:
    #             self._mw.clear_odmr_PushButton.setEnabled(False)
    #             self._mw.action_run_stop.setEnabled(False)
    #             self._mw.action_toggle_cw.setEnabled(True)
    #             dspinbox_dict = self.get_all_dspinboxes_from_groupbox()
    #             for identifier_name in dspinbox_dict:
    #                 dspinbox_type_list = dspinbox_dict[identifier_name]
    #                 [dspinbox_type.setEnabled(True) for dspinbox_type in dspinbox_type_list]
    #             self._mw.odmr_control_DockWidget.add_range_button.setEnabled(True)
    #             self._mw.odmr_control_DockWidget.remove_range_button.setEnabled(True)
    #             self._mw.sweep_power_DoubleSpinBox.setEnabled(True)
    #             self._mw.runtime_DoubleSpinBox.setEnabled(True)
    #             self._sd.clock_frequency_DoubleSpinBox.setEnabled(True)
    #             self._sd.oversampling_SpinBox.setEnabled(True)
    #             self._sd.lock_in_CheckBox.setEnabled(True)
    #             self._mw.action_run_stop.setChecked(False)
    #             self._mw.action_resume_odmr.setChecked(False)
    #             self._mw.action_toggle_cw.setChecked(True)
    #     else:
    #         self._mw.action_resume_odmr.setEnabled(True)
    #         self._mw.cw_power_DoubleSpinBox.setEnabled(True)
    #         self._mw.sweep_power_DoubleSpinBox.setEnabled(True)
    #         self._mw.cw_frequency_DoubleSpinBox.setEnabled(True)
    #         self._mw.clear_odmr_PushButton.setEnabled(False)
    #         self._mw.action_run_stop.setEnabled(True)
    #         self._mw.action_toggle_cw.setEnabled(True)
    #         dspinbox_dict = self.get_all_dspinboxes_from_groupbox()
    #         for identifier_name in dspinbox_dict:
    #             dspinbox_type_list = dspinbox_dict[identifier_name]
    #             [dspinbox_type.setEnabled(True) for dspinbox_type in dspinbox_type_list]
    #         if self._odmr_logic.mw_scanmode.name == 'SWEEP':
    #             self._mw.odmr_control_DockWidget.add_range_button.setDisabled(True)
    #         elif self._odmr_logic.mw_scanmode.name == 'LIST':
    #             self._mw.odmr_control_DockWidget.add_range_button.setEnabled(True)
    #         self._mw.odmr_control_DockWidget.remove_range_button.setEnabled(True)
    #         self._mw.runtime_DoubleSpinBox.setEnabled(True)
    #         self._sd.clock_frequency_DoubleSpinBox.setEnabled(True)
    #         self._sd.oversampling_SpinBox.setEnabled(True)
    #         self._sd.lock_in_CheckBox.setEnabled(True)
    #         self._mw.action_run_stop.setChecked(False)
    #         self._mw.action_resume_odmr.setChecked(False)
    #         self._mw.action_toggle_cw.setChecked(False)
    #
    #     # Unblock signal firing
    #     self._mw.action_run_stop.blockSignals(False)
    #     self._mw.action_resume_odmr.blockSignals(False)
    #     self._mw.action_toggle_cw.blockSignals(False)
    #     return

    # def clear_data(self):
    #     """ Clear the data. """
    #     self.sigClearData.emit()
    #     return

    # def update_plots(self, odmr_data_x, odmr_data_y, odmr_matrix):
    #     """ Refresh the plot widgets with new data. """
    #     # Update mean signal plot
    #     self.odmr_image.setData(odmr_data_x, odmr_data_y[self.display_channel])
    #     # Update raw data matrix plot
    #     cb_range = self.get_matrix_cb_range()
    #     self.update_colorbar(cb_range)
    #     matrix_range = self._mw.odmr_control_DockWidget.matrix_range_SpinBox.value()
    #     start = self._odmr_logic.mw_starts[matrix_range]
    #     step = self._odmr_logic.mw_steps[matrix_range]
    #     stop = self._odmr_logic.mw_stops[matrix_range]
    #     selected_odmr_data_x = np.arange(start, stop, step)
    #
    #     self.odmr_matrix_image.setRect(
    #         QtCore.QRectF(
    #             selected_odmr_data_x[0],
    #             0,
    #             np.abs(selected_odmr_data_x[-1] - selected_odmr_data_x[0]),
    #             odmr_matrix.shape[0])
    #     )
    #
    #     odmr_matrix_range = self._odmr_logic.select_odmr_matrix_data(odmr_matrix, self.display_channel, matrix_range)
    #     self.odmr_matrix_image.setImage(
    #         image=odmr_matrix_range,
    #         axisOrder='row-major',
    #         levels=(cb_range[0], cb_range[1]))

    # def update_channel(self, index):
    #     self.display_channel = int(
    #         self._mw.odmr_channel_ComboBox.itemData(index, QtCore.Qt.UserRole))
    #     self.update_plots(
    #         self._odmr_logic.odmr_plot_x,
    #         self._odmr_logic.odmr_plot_y,
    #         self._odmr_logic.odmr_plot_xy)

    # def average_level_changed(self):
    #     """
    #     Sends to lines to average to the logic
    #     """
    #     self.sigAverageLinesChanged.emit(self._mw.average_level_SpinBox.value())
    #     return

    # def colorscale_changed(self):
    #     """
    #     Updates the range of the displayed colorscale in both the colorbar and the matrix plot.
    #     """
    #     cb_range = self.get_matrix_cb_range()
    #     self.update_colorbar(cb_range)
    #     matrix_image = self.odmr_matrix_image.image
    #     self.odmr_matrix_image.setImage(image=matrix_image, levels=(cb_range[0], cb_range[1]))
    #     return

    # def update_colorbar(self, cb_range):
    #     """
    #     Update the colorbar to a new range.
    #
    #     @param list cb_range: List or tuple containing the min and max values for the cb range
    #     """
    #     self.odmr_cb.refresh_colorbar(cb_range[0], cb_range[1])
    #     return

    # def get_matrix_cb_range(self):
    #     """
    #     Determines the cb_min and cb_max values for the matrix plot
    #     """
    #     matrix_image = self.odmr_matrix_image.image
    #
    #     # If "Manual" is checked or the image is empty (all zeros), then take manual cb range.
    #     # Otherwise, calculate cb range from percentiles.
    #     if self._mw.odmr_cb_manual_RadioButton.isChecked() or np.count_nonzero(matrix_image) < 1:
    #         cb_min = self._mw.odmr_cb_min_DoubleSpinBox.value()
    #         cb_max = self._mw.odmr_cb_max_DoubleSpinBox.value()
    #     else:
    #         # Exclude any zeros (which are typically due to unfinished scan)
    #         matrix_image_nonzero = matrix_image[np.nonzero(matrix_image)]
    #
    #         # Read centile range
    #         low_centile = self._mw.odmr_cb_low_percentile_DoubleSpinBox.value()
    #         high_centile = self._mw.odmr_cb_high_percentile_DoubleSpinBox.value()
    #
    #         cb_min = np.percentile(matrix_image_nonzero, low_centile)
    #         cb_max = np.percentile(matrix_image_nonzero, high_centile)
    #
    #     cb_range = [cb_min, cb_max]
    #     return cb_range

    # def restore_defaultview(self):
    #     self._mw.restoreGeometry(self.mwsettings.value("geometry", ""))
    #     self._mw.restoreState(self.mwsettings.value("windowState", ""))
    #
    # def update_elapsedtime(self):
    #     """ Updates current completed frequency sweeps """
    #     self._mw.elapsed_sweeps_DisplayWidget.display(self._odmr_logic._average_index)
    #
    # def update_elapsed_sweeps(self):
    #     self._mw.elapsed_sweeps_DisplayWidget.display(self._odmr_logic._average_index)

    # def update_settings(self):
    #     """ Write the new settings from the gui to the file. """
    #     number_of_lines = self._sd.matrix_lines_SpinBox.value()
    #     clock_frequency = self._sd.clock_frequency_DoubleSpinBox.value()
    #     oversampling = self._sd.oversampling_SpinBox.value()
    #     lock_in = self._sd.lock_in_CheckBox.isChecked()
    #     self.sigOversamplingChanged.emit(oversampling)
    #     self.sigLockInChanged.emit(lock_in)
    #     self.sigClockFreqChanged.emit(clock_frequency)
    #     self.sigNumberOfLinesChanged.emit(number_of_lines)
    #     return

    # def reject_settings(self):
    #     """ Keep the old settings and restores the old settings in the gui. """
    #     self._sd.matrix_lines_SpinBox.setValue(self._odmr_logic.number_of_lines)
    #     self._sd.clock_frequency_DoubleSpinBox.setValue(self._odmr_logic.clock_frequency)
    #     self._sd.oversampling_SpinBox.setValue(self._odmr_logic.oversampling)
    #     self._sd.lock_in_CheckBox.setChecked(self._odmr_logic.lock_in)
    #     return

    def do_fit(self):
        fit_function = self._mw.fit_methods_ComboBox.getCurrentFit()[0]
        self.sigDoFit.emit(fit_function, None, None, self._mw.fit_range_SpinBox.value())
        return

    def update_fit(self, x_data, y_data, result_str_dict, current_fit):
        """ Update the shown fit. """
        if current_fit != 'No Fit':
            # display results as formatted text
            self._mw.odmr_fit_results_DisplayWidget.clear()
            try:
                formated_results = units.create_formatted_output(result_str_dict)
            except:
                formated_results = 'this fit does not return formatted results'
            self._mw.odmr_fit_results_DisplayWidget.setPlainText(formated_results)

        self._mw.fit_methods_ComboBox.blockSignals(True)
        self._mw.fit_methods_ComboBox.setCurrentFit(current_fit)
        self._mw.fit_methods_ComboBox.blockSignals(False)

        # check which Fit method is used and remove or add again the
        # odmr_fit_image, check also whether a odmr_fit_image already exists.
        if current_fit != 'No Fit':
            self.odmr_fit_image.setData(x=x_data, y=y_data)
            if self.odmr_fit_image not in self._mw.odmr_PlotWidget.listDataItems():
                self._mw.odmr_PlotWidget.addItem(self.odmr_fit_image)
        else:
            if self.odmr_fit_image in self._mw.odmr_PlotWidget.listDataItems():
                self._mw.odmr_PlotWidget.removeItem(self.odmr_fit_image)

        self._mw.odmr_PlotWidget.getViewBox().updateAutoRange()
        return

    def update_fit_range(self):
        self._odmr_logic.range_to_fit = self._mw.fit_range_SpinBox.value()
        return

    # def update_matrix_range(self):
    #     self._odmr_logic.matrix_range = self._mw.odmr_control_DockWidget.matrix_range_SpinBox.value()
    #     # need to update the plot that is showed
    #     key = 'Matrix range: {}'.format(self._odmr_logic.matrix_range)
    #     self.odmr_matrix_image.setImage(self._odmr_logic.select_odmr_matrix_data(self._odmr_logic.odmr_plot_xy,
    #                                                                              self.display_channel,
    #                                                                              self._odmr_logic.matrix_range))
    #     return

    # def update_parameter(self, param_dict):
    #     """ Update the parameter display in the GUI.
    #
    #     @param param_dict:
    #     @return:
    #
    #     Any change event from the logic should call this update function.
    #     The update will block the GUI signals from emitting a change back to the
    #     logic.
    #     """
    #     param = param_dict.get('sweep_mw_power')
    #     if param is not None:
    #         self._mw.sweep_power_DoubleSpinBox.blockSignals(True)
    #         self._mw.sweep_power_DoubleSpinBox.setValue(param)
    #         self._mw.sweep_power_DoubleSpinBox.blockSignals(False)
    #
    #     mw_starts = param_dict.get('mw_starts')
    #     mw_steps = param_dict.get('mw_steps')
    #     mw_stops = param_dict.get('mw_stops')
    #
    #     if mw_starts is not None:
    #         start_frequency_boxes = self.get_freq_dspinboxes_from_groubpox('start')
    #         for mw_start, start_frequency_box in zip(mw_starts, start_frequency_boxes):
    #             start_frequency_box.blockSignals(True)
    #             start_frequency_box.setValue(mw_start)
    #             start_frequency_box.blockSignals(False)
    #
    #     if mw_steps is not None:
    #         step_frequency_boxes = self.get_freq_dspinboxes_from_groubpox('step')
    #         for mw_step, step_frequency_box in zip(mw_steps, step_frequency_boxes):
    #             step_frequency_box.blockSignals(True)
    #             step_frequency_box.setValue(mw_step)
    #             step_frequency_box.blockSignals(False)
    #
    #     if mw_stops is not None:
    #         stop_frequency_boxes = self.get_freq_dspinboxes_from_groubpox('stop')
    #         for mw_stop, stop_frequency_box in zip(mw_stops, stop_frequency_boxes):
    #             stop_frequency_box.blockSignals(True)
    #             stop_frequency_box.setValue(mw_stop)
    #             stop_frequency_box.blockSignals(False)
    #
    #     param = param_dict.get('run_time')
    #     if param is not None:
    #         self._mw.runtime_DoubleSpinBox.blockSignals(True)
    #         self._mw.runtime_DoubleSpinBox.setValue(param)
    #         self._mw.runtime_DoubleSpinBox.blockSignals(False)
    #
    #     param = param_dict.get('number_of_lines')
    #     if param is not None:
    #         self._sd.matrix_lines_SpinBox.blockSignals(True)
    #         self._sd.matrix_lines_SpinBox.setValue(param)
    #         self._sd.matrix_lines_SpinBox.blockSignals(False)
    #
    #     param = param_dict.get('clock_frequency')
    #     if param is not None:
    #         self._sd.clock_frequency_DoubleSpinBox.blockSignals(True)
    #         self._sd.clock_frequency_DoubleSpinBox.setValue(param)
    #         self._sd.clock_frequency_DoubleSpinBox.blockSignals(False)
    #
    #     param = param_dict.get('oversampling')
    #     if param is not None:
    #         self._sd.oversampling_SpinBox.blockSignals(True)
    #         self._sd.oversampling_SpinBox.setValue(param)
    #         self._sd.oversampling_SpinBox.blockSignals(False)
    #
    #     param = param_dict.get('lock_in')
    #     if param is not None:
    #         self._sd.lock_in_CheckBox.blockSignals(True)
    #         self._sd.lock_in_CheckBox.setChecked(param)
    #         self._sd.lock_in_CheckBox.blockSignals(False)
    #
    #     param = param_dict.get('cw_mw_frequency')
    #     if param is not None:
    #         self._mw.cw_frequency_DoubleSpinBox.blockSignals(True)
    #         self._mw.cw_frequency_DoubleSpinBox.setValue(param)
    #         self._mw.cw_frequency_DoubleSpinBox.blockSignals(False)
    #
    #     param = param_dict.get('cw_mw_power')
    #     if param is not None:
    #         self._mw.cw_power_DoubleSpinBox.blockSignals(True)
    #         self._mw.cw_power_DoubleSpinBox.setValue(param)
    #         self._mw.cw_power_DoubleSpinBox.blockSignals(False)
    #
    #     param = param_dict.get('average_length')
    #     if param is not None:
    #         self._mw.average_level_SpinBox.blockSignals(True)
    #         self._mw.average_level_SpinBox.setValue(param)
    #         self._mw.average_level_SpinBox.blockSignals(False)
    #     return

    ############################################################################
    #                           Change Methods                                 #
    ############################################################################
    # Here signals are used to change parameters

    # def change_mw_params(self):
    #     """ Change mw frequency and power """
    #     frequency = self._mw.cw_frequency_DoubleSpinBox.value()
    #     power = self._mw.cw_power_DoubleSpinBox.value()
    #     self.sigMwCwParamsChanged.emit(frequency, power)
    #     return

    def change_sweep_params(self):
        """ Change start, stop and step frequency of frequency sweep """
        starts = []
        steps = []
        stops = []

        num = self._odmr_logic.ranges

        for counter in range(num):
            # construct strings
            start, stop, step = self.get_frequencies_from_row(counter)

            starts.append(start)
            steps.append(step)
            stops.append(stop)

        power = self._mw.sweep_power_DoubleSpinBox.value()
        self.sigMwSweepParamsChanged.emit(starts, stops, steps, power)
        return

    def change_fit_range(self):
        self._odmr_logic.fit_range = self._mw.fit_range_SpinBox.value()
        return

    # def get_frequencies_from_row(self, row):
    #     object_dict = self.get_objects_from_groupbox_row(row)
    #     for object_name in object_dict:
    #         if "DoubleSpinBox" in object_name:
    #             if "start" in object_name:
    #                 start = object_dict[object_name].value()
    #             elif "step" in object_name:
    #                 step = object_dict[object_name].value()
    #             elif "stop" in object_name:
    #                 stop = object_dict[object_name].value()
    #
    #     return start, stop, step

    def change_runtime(self):
        """ Change time after which microwave sweep is stopped """
        runtime = self._mw.runtime_DoubleSpinBox.value()
        self.sigRuntimeChanged.emit(runtime)
        return

    def save_data(self):
        self._odmr_logic.save_data()
