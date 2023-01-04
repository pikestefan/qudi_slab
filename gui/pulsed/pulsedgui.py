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
from gui.guibase import GUIBase
from gui.colordefs import QudiPalettePale as palette
from gui.fitsettings import FitSettingsDialog, FitSettingsComboBox
from qtpy import QtCore, QtWidgets, uic
from qtpy import uic



class PulsedMainWindow(QtWidgets.QMainWindow):
    """ The main window for the ODMR measurement GUI."""
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
    sigDoFit = QtCore.Signal(str)


    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)

    def on_activate(self):
        """ Definition, configuration and initialisation of the ODMR GUI.

        This init connects all the graphic modules, which were created in the
        *.ui file and configures the event handling between the modules.
        """
        self._master_pulselogic = self.master_pulselogic()
        self._mw = PulsedMainWindow()
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
        self._setup_plots()
        self._setup_instance_plot()
        # Grab the saved status variable parameters from the logic
        self._get_parameters_from_logic()


    def _setup_connections(self):
        # this happens during on_activate
        ########################################################################
        #                       Connect signals                                #
        ########################################################################
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
        self._mw.fill_puls_times_pushButton.clicked.connect(self.fill_pulse_times)
        # This button enables the reference counts in the plot
        self._mw.pulse_type_tabWidget.currentChanged.connect(self.tab_changed)
        self._mw.laser_power_2.valueChanged.connect(self.changed_laser_power)
        self._mw.laser_button_cw.clicked.connect(self.changed_laser_power)
        self._mw.laser_button_high.clicked.connect(self.changed_laser_power)
        # instance plot:
        self._mw.plot_instance_Button.clicked.connect(self.instance_plot)
        self._mw.get_values_Button.clicked.connect(self.get_values_plot)
        self._mw.update_rabi_steps.clicked.connect(self.update_mw_values)
        self._mw.update_mw_ramsey.clicked.connect(self.update_mw_values)
        self._mw.update_mw_delay.clicked.connect(self.update_mw_values)
        # # Control/values-changed signals to logic
        self.sigStartMeasurement.connect(self.start_measurement, QtCore.Qt.QueuedConnection)

        self.sigDoFit.connect(self._master_pulselogic.do_fit, QtCore.Qt.QueuedConnection)

        # # Update signals coming from logic:
        self._master_pulselogic.sigMeasurementDone.connect(self.measurement_done,
                                                   QtCore.Qt.QueuedConnection)
        self._master_pulselogic.sigAverageDone.connect(self.refresh_plot, QtCore.Qt.QueuedConnection)
        self._master_pulselogic.sigAverageDone.connect(self.update_curr_av, QtCore.Qt.QueuedConnection)
        self._master_pulselogic.sigFitUpdated.connect(self.update_fit, QtCore.Qt.QueuedConnection)

    def _setup_plots(self):
        # Set up the plots
        self.curr_trace = pg.PlotDataItem(
            pen=pg.mkPen(palette.c1, style=QtCore.Qt.DotLine),
            symbol='o',
            symbolPen=palette.c1,
            symbolBrush=palette.c1,
            symbolSize=7,
        )
        self.average_trace = pg.PlotDataItem(skipFiniteCheck=True, pen=pg.mkPen(color='r'))
        self.average_trace_ref = pg.PlotDataItem(skipFiniteCheck=True, pen=pg.mkPen(color='g'))
        self.fit_image = pg.PlotDataItem(pen=pg.mkPen(palette.c2))

        # This one makes every
        self._mw.pulsed_PlotWidget.addItem(self.curr_trace)
        # This one makes the average signal visible
        self._mw.pulsed_PlotWidget.addItem(self.average_trace_ref)
        self._mw.pulsed_PlotWidget.addItem(self.average_trace)
        self._mw.pulsed_PlotWidget.addItem(self.fit_image)
        self._mw.pulsed_PlotWidget.setLabel('bottom', 'Time in us')
        self._mw.pulsed_PlotWidget.setLabel('left', 'Counts')

    def _setup_instance_plot(self):
        # Set up the plots
        self.mw_trace_i = pg.PlotDataItem(skipFiniteCheck=True, pen=pg.mkPen(color='b'))
        self.mw_trace_q = pg.PlotDataItem(skipFiniteCheck=True, pen=pg.mkPen(color='m'))
        self.apd_trace = pg.PlotDataItem(skipFiniteCheck=True, pen=pg.mkPen(color='r'))
        self.apd_ref_trace = pg.PlotDataItem(skipFiniteCheck=True, pen=pg.mkPen(color='y'))
        self.laser_trace = pg.PlotDataItem(skipFiniteCheck=True, pen=pg.mkPen(color='g'))
        # colors: r, g, b, c, m, y, k, w
        # This one makes every
        self._mw.instance_PlotWidget.addItem(self.mw_trace_i)
        self._mw.instance_PlotWidget.addItem(self.mw_trace_q)
        # This one makes the average signal visible
        self._mw.instance_PlotWidget.addItem(self.apd_trace)
        self._mw.instance_PlotWidget.addItem(self.apd_ref_trace)
        self._mw.instance_PlotWidget.addItem(self.laser_trace)
        self._mw.instance_PlotWidget.setLabel('bottom', 'Time in us')
        self._mw.instance_PlotWidget.setLabel('left', 'Amplitude')

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

        self._mw.mw_stop_distance_rabi.setEnabled(val)
        self._mw.mw_min_len_rabi.setEnabled(val)
        self._mw.mw_start_distance_rabi.setEnabled(val)
        self._mw.mw_steps_rabi.setEnabled(val)

        self._mw.mw_start_distance_ramsey.setEnabled(val)
        self._mw.mw_min_len_ramsey.setEnabled(val)
        self._mw.mw_stop_distance_ramsey.setEnabled(val)
        self._mw.mw_steps_ramsey.setEnabled(val)
        self._mw.laser_power_2.setEnabled(val)
        self._mw.mw_len_ramsey.setEnabled(val)
        # This is the save button
        self._mw.action_Save.setEnabled(val)


    def run_stop_measurement(self, is_checked):
        """ Manages what happens if measurement is started/stopped. """

        if is_checked:
            self.update_mw_values()
            # change the axes appearance according to input values:
            self._set_enabled_ui(False)
            # Reset the sweeps counter
            self._mw.curr_av_DisplayWidget.display(0)

            # emit signal to get scan running
            self.sigStartMeasurement.emit()
        else:
            self._master_pulselogic.stopRequested = True
            self._mw.action_run_stop.setEnabled(True)

            self.index = 0
            # Enable the ui later
            self._set_enabled_ui(True)

        return

    def _send_parameters_to_logic(self):
        # Grabs all the parameters from the GUI and sends them to the master_pulslogic

        # general tab
        mw_power = self._mw.mw_power.value()
        mw_freq = self._mw.mw_freq.value()  # in Hz
        method = self._mw.comboBox_method.currentText()
        averages = self._mw.averages.value()
        int_time = self._mw.integration_time.value()  # in sec
        seq_len = self._mw.seq_len.value() * 1e6  # this should be in us
        sampling_rate_awg = int(self._mw.clk_awg.value() * 1e-6)  # in MHz
        apd_start = self._mw.apd_start_time.value() * 1e6
        apd_len = self._mw.apd_len.value() * 1e6
        apd_ref_start = self._mw.apd_ref_start_time.value() * 1e6
        apd_ref_len = self._mw.apd_ref_len.value() * 1e6

        # laser times tab
        laser_in = self._mw.laser_in.value() * 1e6
        laser_off = self._mw.laser_off.value() * 1e6
        laser_re = self._mw.laser_re.value() * 1e6
        laser_power = self._mw.laser_power_2.value()

        # delay_sweep tab
        apd_len_delay = self._mw.apd_len_pulse_delay.value() * 1e6
        apd_min_start_delay = self._mw.apd_min_start_delay.value() * 1e6
        apd_max_start_delay = self._mw.apd_max_start_delay.value() * 1e6
        apd_steps_delay = self._mw.apd_steps_delay.value() * 1e6
        use_mw_delay = self._mw.use_mw_delay.isChecked()
        mw_start_delay = self._mw.mw_start_time_delay.value() * 1e6
        mw_len_delay = self._mw.mw_len_delay.value() * 1e6

        # rabi tab
        mw_start_distance_rabi = self._mw.mw_start_distance_rabi.value() * 1e6
        mw_min_len_rabi = self._mw.mw_min_len_rabi.value() * 1e6
        mw_stop_distance_rabi = self._mw.mw_stop_distance_rabi.value() * 1e6
        mw_steps_rabi = self._mw.mw_steps_rabi.value() * 1e6

        # ramsey tab
        mw_start_distance_ramsey = self._mw.mw_start_distance_ramsey.value() * 1e6
        mw_min_len_ramsey = self._mw.mw_min_len_ramsey.value() * 1e6
        mw_stop_distance_ramsey = self._mw.mw_stop_distance_ramsey.value() * 1e6
        mw_steps_ramsey = self._mw.mw_steps_ramsey.value() * 1e6
        mw_len_ramsey = self._mw.mw_len_ramsey.value() * 1e6

        # Now send all the parameters from above to the masterlogic
        self._master_pulselogic.mw_power = mw_power
        self._master_pulselogic.mw_freq = mw_freq
        self._master_pulselogic.method = method
        self._master_pulselogic.int_time = int_time
        self._master_pulselogic.seq_len = seq_len
        self._master_pulselogic.clk_rate_awg = sampling_rate_awg
        self._master_pulselogic.averages = averages

        self._master_pulselogic.apd_start = apd_start
        self._master_pulselogic.apd_len = apd_len
        self._master_pulselogic.apd_ref_start = apd_ref_start
        self._master_pulselogic.apd_ref_len = apd_ref_len

        # laser times tab
        self._master_pulselogic.laser_in = laser_in
        self._master_pulselogic.laser_off = laser_off
        self._master_pulselogic.laser_re = laser_re
        self._master_pulselogic.laser_power = laser_power

        # delay_sweep tab
        self._master_pulselogic.apd_len_delay = apd_len_delay
        self._master_pulselogic.apd_min_start_delay = apd_min_start_delay
        self._master_pulselogic.apd_max_start_delay = apd_max_start_delay
        self._master_pulselogic.apd_steps_delay = apd_steps_delay
        self._master_pulselogic.mw_pulse_setting = use_mw_delay
        self._master_pulselogic.mw_start_delay = mw_start_delay
        self._master_pulselogic.mw_len_delay = mw_len_delay

        # rabi tab
        self._master_pulselogic.mw_start_distance_rabi = mw_start_distance_rabi
        self._master_pulselogic.mw_min_len_rabi = mw_min_len_rabi
        self._master_pulselogic.mw_stop_distance_rabi = mw_stop_distance_rabi
        self._master_pulselogic.mw_steps_rabi = mw_steps_rabi

        # ramsey tab
        self._master_pulselogic.mw_start_distance_ramsey = mw_start_distance_ramsey
        self._master_pulselogic.mw_min_len_ramsey = mw_min_len_ramsey
        self._master_pulselogic.mw_stop_distance_ramsey = mw_stop_distance_ramsey
        self._master_pulselogic.mw_steps_ramsey = mw_steps_ramsey
        self._master_pulselogic.mw_len_ramsey = mw_len_ramsey

    def _get_parameters_from_logic(self):
        # Grabs all the parameters from the logic and fills them into the GUI
        # The parameters from the logic are in microseconds

        # Grab the parameters
        mw_power = self._master_pulselogic.mw_power
        mw_freq = self._master_pulselogic.mw_frequency
        method = self._master_pulselogic.method
        int_time = self._master_pulselogic.integration_time
        seq_len = self._master_pulselogic.seq_len
        sampling_rate_awg = self._master_pulselogic.clk_rate_awg
        averages = self._master_pulselogic.averages

        apd_start = self._master_pulselogic.apd_start
        apd_len = self._master_pulselogic.apd_len
        apd_ref_start = self._master_pulselogic.apd_ref_start
        apd_ref_len = self._master_pulselogic.apd_ref_len

        # laser times tab
        laser_in = self._master_pulselogic.laser_in
        laser_off = self._master_pulselogic.laser_off
        laser_re = self._master_pulselogic.laser_re
        laser_power = self._master_pulselogic.laser_power

        # delay_sweep tab
        apd_len_delay = self._master_pulselogic.apd_len_delay
        apd_min_start_delay = self._master_pulselogic.apd_min_start_delay
        apd_max_start_delay = self._master_pulselogic.apd_max_start_delay
        apd_steps_delay = self._master_pulselogic.apd_steps_delay
        use_mw_delay = self._master_pulselogic.mw_pulse_setting
        mw_start_delay = self._master_pulselogic.mw_start_delay
        mw_len_delay = self._master_pulselogic.mw_len_delay

        # rabi tab
        mw_start_distance_rabi = self._master_pulselogic.mw_start_distance_rabi
        mw_min_len_rabi = self._master_pulselogic.mw_min_len_rabi
        mw_stop_distance_rabi = self._master_pulselogic.mw_stop_distance_rabi
        mw_steps_rabi = self._master_pulselogic.mw_steps_rabi

        # ramsey tab
        mw_start_distance_ramsey = self._master_pulselogic.mw_start_distance_ramsey
        mw_min_len_ramsey = self._master_pulselogic.mw_min_len_ramsey
        mw_stop_distance_ramsey = self._master_pulselogic.mw_stop_distance_ramsey
        mw_steps_ramsey = self._master_pulselogic.mw_steps_ramsey
        mw_len_ramsey = self._master_pulselogic.mw_len_ramsey

        # Fill in the values into the GUI
        # Some values need to be converted from us -> s

        # general tab
        self._mw.mw_power.setValue(mw_power)
        self._mw.mw_freq.setValue(mw_freq)  # in Hz
        self._mw.comboBox_method.setCurrentText(method)
        self._mw.averages.setValue(averages)
        self._mw.integration_time.setValue(int_time)  # in sec
        self._mw.seq_len.setValue(seq_len * 1e-6)  # this should be in us
        self._mw.clk_awg.setValue(int(sampling_rate_awg * 1e6))  # in MHz
        self._mw.apd_start_time.setValue(apd_start * 1e-6)
        self._mw.apd_len.setValue(apd_len * 1e-6)
        self._mw.apd_ref_start_time.setValue(apd_ref_start * 1e-6)
        self._mw.apd_ref_len.setValue(apd_ref_len * 1e-6)

        # laser times tab
        self._mw.laser_in.setValue(laser_in * 1e-6)
        self._mw.laser_off.setValue(laser_off * 1e-6)
        self._mw.laser_re.setValue(laser_re * 1e-6)
        self._mw.laser_power_2.setValue(laser_power)

        # delay_sweep tab
        self._mw.apd_len_pulse_delay.setValue(apd_len_delay * 1e-6)
        self._mw.apd_min_start_delay.setValue(apd_min_start_delay * 1e-6)
        self._mw.apd_max_start_delay.setValue(apd_max_start_delay * 1e-6)
        self._mw.apd_steps_delay.setValue(apd_steps_delay * 1e-6)
        self._mw.use_mw_delay.setChecked(use_mw_delay)
        self._mw.mw_start_time_delay.setValue(mw_start_delay * 1e-6)
        self._mw.mw_len_delay.setValue(mw_len_delay * 1e-6)

        # rabi tab
        self._mw.mw_start_distance_rabi.setValue(mw_start_distance_rabi * 1e-6)
        self._mw.mw_min_len_rabi.setValue(mw_min_len_rabi * 1e-6)
        self._mw.mw_stop_distance_rabi.setValue(mw_stop_distance_rabi * 1e-6)
        self._mw.mw_steps_rabi.setValue(mw_steps_rabi * 1e-6)

        # ramsey tab
        self._mw.mw_start_distance_ramsey.setValue(mw_start_distance_ramsey * 1e-6)
        self._mw.mw_min_len_ramsey.setValue(mw_min_len_ramsey * 1e-6)
        self._mw.mw_stop_distance_ramsey.setValue(mw_stop_distance_ramsey * 1e-6)
        self._mw.mw_steps_ramsey.setValue(mw_steps_ramsey * 1e-6)
        self._mw.mw_len_ramsey.setValue(mw_len_ramsey * 1e-6)

    def start_measurement(self):
        self._master_pulselogic.start_stop_timetrace(False)
        self._send_parameters_to_logic()
        self.fit_image.clear()
        # This one sets the x axis right
        x_axis = self._master_pulselogic.get_x_axis() # in microseconds
        self._mw.pulsed_PlotWidget.setXRange(min(x_axis), max(x_axis))
        # This is where the measurement really starts
        self._master_pulselogic.start_measurement()

    def clear_all(self): # This gets enabled by the clear awg button?
        # It stops the replay and clears the memory
        self._master_pulselogic.stop_awg()

    def measurement_done(self):
        # Turn the button back to the play icon
        # Enable the saving button
        self._set_enabled_ui(True)
        self._mw.action_run_stop.setChecked(False)
        self._mw.action_Save.setEnabled(True)
        self.index = 0

    def cw_mode(self, is_checked):
    #checks if the House button is checked in general
        if is_checked:
            self._master_pulselogic.start_stop_timetrace(True)
            # print(self._mw.laser_power_2.valueChanged(), self._mw.laser_power_2.valueChanged)
            laser_power = self._mw.laser_power_2.value() # Takes the value from the box
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

    def changed_laser_power(self):
        if self._mw.action_cw_mode.isChecked():
            # print('is checked')
            if self._mw.laser_button_high.isChecked():
                self._master_pulselogic.laser_power = 0.8
                self._master_pulselogic.set_laser_power()
            elif self._mw.laser_button_cw.isChecked():
                self._master_pulselogic.laser_power = 0.628
                self._master_pulselogic.set_laser_power()
            else:
                laser_power = self._mw.laser_power_2.value()
                self._master_pulselogic.laser_power = laser_power
                self._master_pulselogic.set_laser_power()
        else:
            pass
            # print('is not checked')
        return


    def show(self):
        """Make window visible and put it above all other windows. """
        self._mw.show()
        self._mw.activateWindow()
        self._mw.raise_()

    def tab_changed(self, current_tab):
        # For now, ignore the tab and just send the parameters to the logic
        # This way they can be saved as status variables
        self._send_parameters_to_logic()

    # def get_value(self, counts, ref_counts):
    #     self.index += 1
    #     x_axis_len = len(self._master_pulselogic.get_x_axis())
    #     steps = self._master_pulselogic.step_counter()[0]
    #     if self.index == steps:
    #         self.index == 0
    #         # Average is done
    #         # begin a new average
    #     else:
    #         self.index +=1
    #         #continue plotting

    def refresh_plot(self, av_index, current_row, current_average, av_counts_ref):
        # Draw current odmr trace
        # This is where the data from the logic should be
        # maybe have a method on the masterpulselogic which gives all values in one row? They give one plot in the GUI
        self.curr_trace.setData(self._master_pulselogic.get_x_axis(), current_row)

        # Draw average trace
        if av_index > 0:
            self.average_trace.setData(self._master_pulselogic.get_x_axis(), current_average)
        else:
            self.average_trace.clear()

        if self._mw.show_reference_counts.isChecked():
            self.average_trace_ref.setData(self._master_pulselogic.get_x_axis(), av_counts_ref)
        else:
            self.average_trace_ref.clear()

    def show_ref_counts(self, av_index, ref_average):
        # Get ref_count data from the masterpulse_logic
        if av_index > 0: #leave out the first measurement
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
        self.fit_image.setData(x=x_fit, y=y_fit)
        self._mw.pulsed_PlotWidget.getViewBox().updateAutoRange()

        return

    def fill_pulse_times(self):
        if 'Rabi' not in self._master_pulselogic.fits_performed:
            return

        # Extract Rabi frequency
        _, _, result = self._master_pulselogic.fits_performed['Rabi']
        rabi_freq = result.values['frequency']
        pi_puls = 1 / (2 * rabi_freq)
        pi_half_puls = pi_puls / 2

        # Fill in Puls Times
        self._mw.mw_len_delay.setValue(pi_puls)

    def do_fit(self):
        fit_function = self._mw.fit_methods_ComboBox.getCurrentFit()[0]
        self.sigDoFit.emit(fit_function)


    def change_fit_range(self):
        self._odmr_logic.fit_range = self._mw.fit_range_SpinBox.value()
        return


    def save_data(self):
        method = self._mw.comboBox_method.currentText()
        self._master_pulselogic.save_data(method)
    def get_values_plot(self):
        self._mw.get_values_Button.setEnabled(False)
        self._send_parameters_to_logic()    #### This helps!
        self.laser_array = []
        self.mw_i_array = []
        self.mw_q_array = []
        self.apd_array = []
        self.apd_ref_array = []
        # This provides the right values
        method, laser_array, mw_i_array, mw_q_array, apd_array, apd_ref_array = self._master_pulselogic.get_revant_parameters()
        self.method = method
        self.laser_array = laser_array
        self.mw_i_array = mw_i_array
        self.mw_q_array = mw_q_array
        self.apd_array = apd_array
        self.apd_ref_array = apd_ref_array
        self.t_axis = self._master_pulselogic.get_t_axis_pulselogic()  # in microseconds

        # # This tries to replace the values in the mask with real values if the get values button is clicked!
        # method, min_time, req_len, real_len = self._master_pulselogic.len_errors()
        #
        # if method == 'rabi':
        #     self._mw.mw_steps_rabi.setValue(real_len * 1e-6)
        # elif method == 'ramsey':
        #     self._mw.mw_len_ramsey.setValue(real_len * 1e-6)

        self._mw.plot_instance_Button.setEnabled(True)

    def instance_plot(self):
        self._mw.get_values_Button.setEnabled(True)
        instance_number = int(self._mw.instance_number.value())
        self._mw.pulsed_PlotWidget.setXRange(min(self.t_axis), max(self.t_axis))
        if self._mw.laser_plot.isChecked():
            self.laser_trace.setData(self.t_axis, self.laser_array[instance_number])
        else:
            self.laser_trace.clear()
        if self._mw.mw_plot_i.isChecked():
            self.mw_trace_i.setData(self.t_axis, self.mw_i_array[instance_number])
        else:
            self.mw_trace_i.clear()
        if self._mw.mw_plot_q.isChecked():
            self.mw_trace_q.setData(self.t_axis, self.mw_q_array[instance_number])
        else:
            self.mw_trace_q.clear()
        if self._mw.apd_plot.isChecked():
            self.apd_trace.setData(self.t_axis, self.apd_array[instance_number])
        else:
            self.apd_trace.clear()
        if self._mw.apd_ref_plot.isChecked():
            self.apd_ref_trace.setData(self.t_axis, self.apd_ref_array[instance_number])
        else:
            self.apd_ref_trace.clear()

    def update_curr_av(self):
        """ Updates current completed frequency sweeps """
        self._mw.curr_av_DisplayWidget.display(self._master_pulselogic.av_index)

    ### new attemt for rouding due to clock_rate
    def update_mw_values(self):
        ''' This tries to replace the values in the mask with real values
        #  for now it only works for rabi
        # call this function before the values get fished from the GUI
        # to make sure that the updated values get used in the actual plot and measurement!
        # this one gets the needed values from the masterpulselogic
        # These values have 5 digits
        # Here we need to differentiate between different methods
        '''
        clk_rate = self._mw.clk_awg.value() * 1e-6
        method = self._mw.comboBox_method.currentText()
        if method == 'rabi':
            mw_steps = self._mw.mw_steps_rabi.value() * 1e6
        elif method == 'ramsey':
            mw_steps = self._mw.mw_len_ramsey.value() * 1e6
        elif method == 'delaysweep' or method == 'delaysweep_ref':
            mw_steps = self._mw.mw_len_delay.value() * 1e6

        method, min_time, req_len, real_len = self._master_pulselogic.len_errors(mw_steps, clk_rate)

        if method == 'rabi':
            self._mw.mw_steps_rabi.setValue(real_len * 1e-6)
        elif method == 'ramsey':
            self._mw.mw_len_ramsey.setValue(real_len * 1e-6)
        elif method == 'delaysweep' or method == 'delaysweep_ref':
            self._mw.mw_len_delay.setValue(real_len * 1e-6)






