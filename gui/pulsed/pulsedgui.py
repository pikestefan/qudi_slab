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


class PulsedGui(GUIBase):
    """
    This is the GUI Class for pulsed measurements
    """

    # declare connectors to master pulse logic
    master_pulselogic = Connector(interface='MasterPulse')

    # some signals which are also used in the master pulse logic
    sigStartMeasurement = QtCore.Signal() # This starts the measurement
    sigGetValues = QtCore.Signal() # This gets the value from the masterpulse logic when the ref button is pushed
    sigDoFit = QtCore.Signal(str)

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
        self._mw.laser_power_2.setMaximum(1)
        self._mw.laser_power_2.setMinimum(0)

        self.index = 0

        # fit settings
        self._fsd = FitSettingsDialog(self._master_pulselogic.fc)
        self._fsd.sigFitsUpdated.connect(self._mw.fit_methods_ComboBox.setFitFunctions)
        self._fsd.applySettings()
        self._mw.action_FitSettings.triggered.connect(self._fsd.show)

        # set up all the important connections:
        self._setup_connections()

        # Set up plot
        self.curr_trace = pg.PlotDataItem(
            pen=pg.mkPen(palette.c1, style=QtCore.Qt.DotLine),
            symbol='o',
            symbolPen=palette.c1,
            symbolBrush=palette.c1,
            symbolSize=7,
        )
        # self.curr_trace_ref = pg.PlotDataItem(
        #     pen=pg.mkPen(palette.c1, style=QtCore.Qt.DotLine),
        #     symbol='x',
        #     symbolPen=palette.c1,
        #     symbolBrush=palette.c1,
        #     symbolSize=7,
        # )
        self.average_trace = pg.PlotDataItem(skipFiniteCheck=True, pen=pg.mkPen(color='r'))
        self.average_trace_ref = pg.PlotDataItem(skipFiniteCheck=True, pen=pg.mkPen(color='g'))
        self.fit_image = pg.PlotDataItem(pen=pg.mkPen(palette.c3))

        # This one makes every
        self._mw.pulsed_PlotWidget.addItem(self.curr_trace)
        # This one makes the average signal visible
        self._mw.pulsed_PlotWidget.addItem(self.average_trace)
        self._mw.pulsed_PlotWidget.addItem(self.fit_image)
        self._mw.pulsed_PlotWidget.setLabel('bottom', 'Time in us')
        self._mw.pulsed_PlotWidget.setLabel('left', 'Counts')

    def _setup_connections(self):
        # this happens during on_activate
        ########################################################################
        #                       Connect signals                                #
        ########################################################################
        # # Internal user input changed signals


        #
        # # Internal trigger signals --> connect the signals from the buttons to the gui methods
        # These are the buttons in the toolbar of the GUI
        # This is the red run button on the left
        self._mw.action_run_stop.triggered.connect(self.run_stop_measurement)
        # This is the clear awg button
        self._mw.clear_awg.triggered.connect(self.clear_all)
        # This is the save button
        self._mw.action_Save.triggered.connect(self.save_data)
        # This should be the clear data in the plot button:
        #self._mw.clear_data.triggered.connect(self.clear_plot)
        # This button goes back to cw mode where we can see the laser
        self._mw.action_cw_mode.triggered.connect(self.cw_mode)
        self._mw.do_fit_PushButton.clicked.connect(self.do_fit)
        # This button enables the reference counts in the plot
        self._mw.show_reference_counts.triggered.connect(self.ref_button)

        # # Control/values-changed signals to logic
        self.sigStartMeasurement.connect(self.start_measurement, QtCore.Qt.QueuedConnection)
        self.sigGetValues.connect(self.show_ref_counts, QtCore.Qt.QueuedConnection)


        self.sigDoFit.connect(self._master_pulselogic.do_fit, QtCore.Qt.QueuedConnection)

        # # Update signals coming from logic:
        self._master_pulselogic.sigMeasurementDone.connect(self.measurement_done,
                                                   QtCore.Qt.QueuedConnection)
        self._master_pulselogic.sigAverageDone.connect(self.refresh_plot, QtCore.Qt.QueuedConnection)
        self._master_pulselogic.sigAverageDone.connect(self.show_ref_counts, QtCore.Qt.QueuedConnection)
        self._master_pulselogic.sigFitUpdated.connect(self.update_fit, QtCore.Qt.QueuedConnection)



    def on_deactivate(self):
        """ Reverse steps of activation

        @return int: error code (0:OK, -1:error)
        """

        # self.sigDoFit.disconnect()
        self._mw.action_run_stop.triggered.disconnect()
        self._mw.clear_awg.triggered.disconnect()
        self._mw.action_Save.triggered.disconnect()
        self._mw.action_cw_mode.triggered.disconnect()
        self._mw.show_reference_counts.triggered.disconnect()
        self.sigStartMeasurement.disconnect()
        self._mw.close()
        return 0

    def _set_enabled_ui(self, val):
        """Set the enabled/disabled state for the pulsed setting ui"""
        # All the spinBoxes
        self._mw.mw_power.setEnabled(val)
        self._mw.mw_freq.setEnabled(val)
        self._mw.comboBox_method.setEnabled(val)
        self._mw.averages.setEnabled(val)
        self._mw.integration_time.setEnabled(val)
        self._mw.seq_len.setEnabled(val)
        self._mw.clk_awg.setEnabled(val)
        self._mw.apd_start_time.setEnabled(val)
        self._mw.apd_ref_len.setEnabled(val)
        self._mw.apd_len.setEnabled(val)
        self._mw.apd_ref_start_time.setEnabled(val)
        self._mw.laser_in.setEnabled(val)
        self._mw.laser_off.setEnabled(val)
        self._mw.laser_re.setEnabled(val)
        self._mw.apd_len_pulse_delay.setEnabled(val)
        self._mw.apd_min_start_delay.setEnabled(val)
        self._mw.apd_max_start_delay.setEnabled(val)
        self._mw.apd_steps_delay.setEnabled(val)
        self._mw.use_mw_delay.setEnabled(val)
        self._mw.mw_start_time_delay.setEnabled(val)
        self._mw.mw_len_delay.setEnabled(val)
        self._mw.mw_start_time_rabi.setEnabled(val)
        self._mw.mw_min_len_rabi.setEnabled(val)
        self._mw.mw_max_len_rabi.setEnabled(val)
        self._mw.mw_steps_rabi.setEnabled(val)
        self._mw.use_neg_mw.setEnabled(val)
        self._mw.mw_start_time_ramsey.setEnabled(val)
        self._mw.mw_min_len_ramsey.setEnabled(val)
        self._mw.mw_max_len_ramsey.setEnabled(val)
        self._mw.mw_steps_ramsey.setEnabled(val)
        self._mw.laser_power_2.setEnabled(val)
        self._mw.mw_len_ramsey.setEnabled(val)
        # buttons
        # self._mw.do_fit_PushButton.setEnabled(val)
        # This is the save button
        self._mw.action_Save.setEnabled(val)


    def run_stop_measurement(self, is_checked):
        """ Manages what happens if measurement is started/stopped. """

        if is_checked:
            # change the axes appearance according to input values:

            # Disable the ui later
            self._set_enabled_ui(False)

            # Reset the sweeps counter
            # self._mw.elapsed_sweeps_DisplayWidget.display(0)
            # emit signal to get scan running
            self.sigStartMeasurement.emit()
        else:
            self._master_pulselogic.stopRequested = True
            self._mw.action_run_stop.setEnabled(True)

            self.index = 0
            # Enable the ui later
            self._set_enabled_ui(True)

        return

    def start_measurement(self): # This get enabled with the start button?
        # Grab all the parameters from the GUI
        # general tab
        mw_power = self._mw.mw_power.value()
        mw_freq = self._mw.mw_freq.value()  #in Hz
        method = self._mw.comboBox_method.currentText()
        averages = self._mw.averages.value()
        int_time = self._mw.integration_time.value()*1e-3  #in sec
        seq_len = self._mw.seq_len.value()*1e6 #this should be in us
        sampling_rate_awg = int(self._mw.clk_awg.value())

        apd_start = self._mw.apd_start_time.value()*1e6
        apd_len = self._mw.apd_len.value()*1e6

        apd_ref_start = self._mw.apd_ref_start_time.value()*1e6
        apd_ref_len = self._mw.apd_ref_len.value()*1e6

        # laser times tab
        laser_in = self._mw.laser_in.value()*1e6
        laser_off = self._mw.laser_off.value()*1e6
        laser_re = self._mw.laser_re.value()*1e6
        laser_power = self._mw.laser_power_2.value()

        # delay_sweep tab
        apd_len_delay = self._mw.apd_len_pulse_delay.value()*1e6
        apd_min_start_delay = self._mw.apd_min_start_delay.value()*1e6
        apd_max_start_delay = self._mw.apd_max_start_delay.value()*1e6
        apd_steps_delay = self._mw.apd_steps_delay.value()*1e6

        use_mw_delay = self._mw.use_mw_delay.isChecked()
        mw_start_delay = self._mw.mw_start_time_delay.value()*1e6
        mw_len_delay = self._mw.mw_len_delay.value()*1e6

        # rabi tab
        mw_start_time_rabi = self._mw.mw_start_time_rabi.value()*1e6
        mw_min_len_rabi = self._mw.mw_min_len_rabi.value()*1e6
        mw_max_len_rabi = self._mw.mw_max_len_rabi.value()*1e6
        mw_steps_rabi = self._mw.mw_steps_rabi.value()*1e6
        use_neg_mw = self._mw.use_neg_mw.isChecked()

        # ramsey tab
        mw_start_time_ramsey = self._mw.mw_start_time_ramsey.value()*1e6
        mw_min_len_ramsey = self._mw.mw_min_len_ramsey.value()*1e6
        mw_max_len_ramsey = self._mw.mw_max_len_ramsey.value()*1e6
        mw_steps_ramsey = self._mw.mw_steps_ramsey.value()*1e6
        mw_len_ramsey = self._mw.mw_len_ramsey.value()*1e6

        # now send all the parameters from above to the masterlogic

        self._master_pulselogic.laser_times = [laser_in, laser_off, laser_re]
        self._master_pulselogic.laser_power = laser_power
        self._master_pulselogic.mw_times_ramsey = [mw_start_time_ramsey, mw_min_len_ramsey, mw_max_len_ramsey, mw_steps_ramsey, mw_len_ramsey]
        self._master_pulselogic.mw_times_rabi = [mw_start_time_rabi, mw_min_len_rabi, mw_max_len_rabi, mw_steps_rabi]
        self._master_pulselogic.neg_mw = use_neg_mw
        self._master_pulselogic.apd_times = [apd_start, apd_len]
        self._master_pulselogic.apd_ref_times = [apd_ref_start, apd_ref_len]
        self._master_pulselogic.apd_times_sweep = [apd_len_delay, apd_min_start_delay, apd_max_start_delay, apd_steps_delay]
        self._master_pulselogic.averages = averages
        self._master_pulselogic.mw_power = mw_power
        self._master_pulselogic.mw_frequency = mw_freq
        self._master_pulselogic.integration_time = int_time
        self._master_pulselogic.method = method
        self._master_pulselogic.seq_len = seq_len
        self._master_pulselogic.clk_rate_awg = sampling_rate_awg
        self._master_pulselogic.mw_pulse_setting = use_mw_delay
        self._master_pulselogic.mw_times_sweep = [mw_start_delay, mw_len_delay]


        # This one sets the x axis right
        x_axis = self._master_pulselogic.get_x_axis()
        self._mw.pulsed_PlotWidget.setXRange(min(x_axis), max(x_axis))
        #This is where the measurement really starts
        self._master_pulselogic.start_measurement()

    def clear_all(self): # This gets enabled by the clear awg button?
        # It stops the replay and clears the memory
        self._master_pulselogic.stop_awg()

    def measurement_done(self):
        # self._set_enabled_pulsed_ui(True)
        # Turn the button back to the play icon
        # self._set_enabled_ui(True)
        # Enable the saving button
        self._set_enabled_ui(True)
        self._mw.action_run_stop.setChecked(False)
        self._mw.action_Save.setEnabled(True)
        self.index = 0
        #print('ready for saving process')

    def cw_mode(self, is_checked):
        if is_checked:
            laser_power = self._mw.laser_power_2.value()
            self._master_pulselogic.laser_power = laser_power
            self._master_pulselogic.cw(True)
            self._master_pulselogic.set_laser_power()
            self._mw.action_run_stop.setEnabled(False)
            self._mw.action_Save.setEnabled(False)
            # Laser is on and counts are in the timeseries
        else:
            self._master_pulselogic.cw(False)
            self._mw.action_Save.setEnabled(True)
            self._mw.action_run_stop.setEnabled(True)
            # Laser is off but the counts are still in the timeseries
        return

    def show(self):
        """Make window visible and put it above all other windows. """
        self._mw.show()
        self._mw.activateWindow()
        self._mw.raise_()

    def get_value(self, counts, ref_counts):
        self.index += 1
        x_axis_len = len(self._master_pulselogic.get_x_axis())
        steps = self._master_pulselogic.step_counter()[0]
        if self.index == steps:
            self.index == 0
            # Average is done
            # begin a new average
        else:

            self.index +=1
            #continue plotting

    def refresh_plot(self, av_index, current_row, current_average, ref_average):
        # Draw current odmr trace
        # print('curent row', current_row)
        # This is where the data from the logic should be
        #maybe have a method on the masterpulselogic which gives all values in one row? They give one plot in the GUI
        self.curr_trace.setData(self._master_pulselogic.get_x_axis(), current_row)

        # Draw average trace
        if av_index > 0:
            self.average_trace.setData(self._master_pulselogic.get_x_axis(), current_average)
        else:
            self.average_trace.clear()

    def ref_button(self, is_checked):
        """ Manages what happens if measurement is started/stopped. """

        if is_checked:
            self.sigGetValues.emit()
            # get the values from the get_values method
            # and plots them

        else:
            # dont do anything
            pass

    def show_ref_counts(self, av_index, current_row, current_average, ref_average):
        # This should fish the signal from sigAverageDone from masterpulse logic
        # Draw average trace
        if av_index > 0:
            self.average_trace_ref.setData(self._master_pulselogic.get_x_axis(), ref_average)
        else:
            self.average_trace.clear()

    def update_fit(self, x_fit, y_fit, result_str_dict, current_fit):
        print(x_fit, y_fit, result_str_dict, current_fit)
        """ Update the shown fit. """
        if current_fit != 'No Fit':
            # display results as formatted text
            self._mw.fit_results_DisplayWidget.clear()
            try:
                formated_results = units.create_formatted_output(result_str_dict)
            except:
                formated_results = 'this fit does not return formatted results'
            self._mw.fit_results_DisplayWidget.setPlainText(formated_results)

        self._mw.fit_methods_ComboBox.blockSignals(True)
        self._mw.fit_methods_ComboBox.setCurrentFit(current_fit)
        self._mw.fit_methods_ComboBox.blockSignals(False)

        # check which Fit method is used and remove or add again the
        # odmr_fit_image, check also whether a odmr_fit_image already exists.
        # if current_fit == 'No Fit':
        #     return

        print('x_fit')
        print(x_fit)
        print('y_fit')
        print(y_fit)

        self.fit_image.setData(x=x_fit, y=y_fit)
        # if self.fit_image not in self._mw.pulsed_PlotWidget.listDataItems():
        #     self._mw.pulsed_PlotWidget.addItem(self.fit_image)
        # else:
        #     if self.fit_image in self._mw.pulsed_PlotWidget.listDataItems():
        #         self._mw.pulsed_PlotWidget.removeItem(self.fit_image)

        self._mw.pulsed_PlotWidget.getViewBox().updateAutoRange()
        return


    def do_fit(self):
        fit_function = self._mw.fit_methods_ComboBox.getCurrentFit()[0]
        self.sigDoFit.emit(fit_function)

    # def update_fit(self, x_data, y_data, result_str_dict, current_fit):
    #     """ Update the shown fit. """
    #     if current_fit != 'No Fit':
    #         # display results as formatted text
    #         self._mw.fit_results_DisplayWidget.clear()
    #         try:
    #             formated_results = units.create_formatted_output(result_str_dict)
    #         except:
    #             formated_results = 'this fit does not return formatted results'
    #         self._mw.fit_results_DisplayWidget.setPlainText(formated_results)
    #
    #     self._mw.fit_methods_ComboBox.blockSignals(True)
    #     self._mw.fit_methods_ComboBox.setCurrentFit(current_fit)
    #     self._mw.fit_methods_ComboBox.blockSignals(False)
    #
    #     # check which Fit method is used and remove or add again the
    #     # odmr_fit_image, check also whether a odmr_fit_image already exists.
    #     if current_fit == 'No Fit':
    #         return
    #
    #     self.fit_image.setData(x=x_data, y=y_data)
    #     self._mw.pulsed_PlotWidget.getViewBox().updateAutoRange()
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


    def save_data(self):
        method = self._mw.comboBox_method.currentText()
        self._master_pulselogic.save_data(method)

