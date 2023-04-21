"""
This file contains the Qudi data object classes needed for pulse sequence generation.

Qudi is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Qudi is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Qudi. If not, see http://www.gnu.org/licenses/.

Copyright (c) the Qudi Developers. See the COPYRIGHT.txt file at the
top-level directory of this distribution and at https://github.com/Ulm-IQO/qudi/
"""

from qtpy import QtCore
from collections import OrderedDict
import numpy as np
import math
import time
import datetime
from core.connector import Connector
from core.util.mutex import Mutex
from logic.generic_logic import GenericLogic
from core.configoption import ConfigOption
from core.statusvariable import StatusVar


# Just hardcode the channels to the different hardware modules somewehere
# Microwave: analog channel 2 and 3 for I and Q parameter. For now we just need I
# Laser: digital channel X0
# Photon counter: digital channel X1 and X2 (for reference)


class MasterPulse(GenericLogic):
    pulsegenerator = Connector(interface="PulserInterface")
    laser = Connector(interface="SimpleLaserInterface")
    photon_counter = Connector(interface="SnvmScannerInterface")
    mw_source = Connector(interface="MicrowaveInterface")
    savelogic = Connector(interface="SaveLogic")
    pulselogic = Connector(interface="Pulse")
    fitlogic = Connector(interface="FitLogic")
    timeseries_logic = Connector(interface="TimeSeriesReaderLogic")

    ## Define status variables
    fc = StatusVar("fits", None)
    # status vars
    averages = StatusVar("averages", 10)  # max value for rows
    mw_frequency = StatusVar("mw_frequency", 2.8685e9)
    mw_power = StatusVar("mw_power", -30)  # in dBm
    integration_time = StatusVar("integration_time", 30e-3)  # in s
    laser_power = StatusVar("laser_power", 0.628)  # in V... Max is 1V!
    clk_rate_awg = StatusVar("clk_rate_awg", int(1000e6))  # in MHz
    seq_map_req = StatusVar("seq_map_req", False)
    plaintext_field = StatusVar("notes", "Notes: Magnet position:")

    ## Delay Sweep
    mw_pulse_setting = StatusVar(
        "mw_pulse_setting", False
    )  # this is only for delay sweep
    neg_mw = StatusVar(
        "neg_mw", False
    )  # if False it plays +1, if True it applies negative mw pulses -1 (all the pulses are squared)
    method = StatusVar(
        "method", "ramsey"
    )  # change this to 'ramsey', 'delaysweep', 'rabi' or 'delaysweep_ref' if needed

    # StatusVar for arrays:
    apd_start = StatusVar("apd_start", 3.05e-6)
    apd_len = StatusVar("apd_len", 250e-9)
    apd_ref_start = StatusVar("apd_ref_start", 4.5e-6)
    apd_ref_len = StatusVar("apd_ref_len", 100e-9)

    # laser times tab
    laser_in = StatusVar("laser_in", 1e-6)
    laser_off = StatusVar("laser_off", 2e-6)
    laser_re = StatusVar("laser_re", 2e-6)

    # delay_sweep tab
    apd_len_delay = StatusVar("apd_len_delay", 50e-9)
    apd_min_start_delay = StatusVar("apd_min_start_delay", 0e-6)
    apd_max_start_delay = StatusVar("apd_max_start_delay", 2.6e-6)
    apd_steps_delay = StatusVar("apd_steps_delay", 50e-9)

    mw_start_delay = StatusVar("mw_start_delay", 1.3e-6)
    mw_len_delay = StatusVar("mw_len_delay", 500e-9)  # What ever the pi pulse is

    # rabi tab
    mw_start_distance_rabi = StatusVar("mw_start_distance_rabi", 1.3e-6)
    mw_min_len_rabi = StatusVar("mw_min_len_rabi", 0e-6)
    mw_stop_distance_rabi = StatusVar("mw_stop_distance_rabi", 900e-9)
    mw_steps_rabi = StatusVar("mw_steps_rabi", 500e-9)

    # ramsey tab
    mw_start_distance_ramsey = StatusVar("mw_start_distance_ramsey", 1.3e-6)
    mw_min_len_ramsey = StatusVar("mw_min_len_ramsey", 0e-6)
    mw_stop_distance_ramsey = StatusVar("mw_stop_distance_ramsey", 900e-9)
    mw_steps_ramsey = StatusVar("mw_steps_ramsey", 50e-9)
    mw_len_ramsey = StatusVar("mw_len_ramsey", 16e-9)  # whatever the rabi says
    phase_shift = StatusVar("phase shift", False)

    # Internal signals
    sigContinueLoop = QtCore.Signal()
    sigStopMeasurement = QtCore.Signal()
    sigStopAll = QtCore.Signal()
    sigPxAcquired = QtCore.Signal()

    # Update signals, e.g. for GUI module
    sigMeasurementDone = QtCore.Signal()
    # Gives the current count to the GUi
    sigAverageDone = QtCore.Signal(
        int, np.ndarray, np.ndarray, np.ndarray
    )  # number of averages, counts and average counts
    sigFitUpdated = QtCore.Signal(np.ndarray, np.ndarray, dict, str)
    sigStartMeasurent = QtCore.Signal()
    sigSeqPlaying = QtCore.Signal(bool)

    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)
        self.threadlock = Mutex()

    def on_activate(self):
        # Initialisation performed during activation of the module

        # Get connectors
        self._pulser = self.pulsegenerator()
        self._laser = self.laser()
        self._photon_counter = self.photon_counter()
        self._mw_device = self.mw_source()
        self._savelogic = self.savelogic()
        self._pulselogic = self.pulselogic()
        self._fitlogic = self.fitlogic()
        self._timeseries = self.timeseries_logic()
        active_channels = self._pulser.get_active_channels()
        self._pulser.clear_all()
        self._pulser.reset()

        # These are the index for the matrix (step_index: columns and av_index: rows)
        #  Elapsed measurement time and number of sweeps
        self.elapsed_time = 0.0
        self._photon_counter.get_counter_clock_frequency()
        # Initalize the values which get stored later
        self.count_matrix = None  # This matrix is used to store the DAQ card values.
        self.av_index = 0  # value right now, row
        self.step_index = 0  # value right now, column

        # self.step_count = self._pulselogic.get_step_count([2.1, 1, 2, 0.5])  #max value for columns
        self.clk_rate_daq = (
            self._photon_counter.get_counter_clock_frequency()
        )  # which is 100000

        self.rep = 100
        self.trigger_setting = True

        ## Fit
        self.fit_x = None
        self.fit_y = None
        self.fits_performed = {}

        # Declare where the signals should lead to
        self.sigContinueLoop.connect(self.continue_loop, QtCore.Qt.QueuedConnection)
        self.sigPxAcquired.connect(self._next_pixel, QtCore.Qt.QueuedConnection)
        self.sigStopMeasurement.connect(
            self.stop_measurement, QtCore.Qt.QueuedConnection
        )
        self.sigStopAll.connect(self.stop_all, QtCore.Qt.QueuedConnection)
        self.sigStartMeasurent.connect(
            self.start_measurement, QtCore.Qt.QueuedConnection
        )
        # self._pulselogic.sigUsedArrays.connect(self.get_used_arrays, QtCore.Qt.QueuedConnection)
        self.stopRequested = False
        self._build_pulse_arrays()

        return active_channels

    def _build_pulse_arrays(self):
        # and calculate the real values from the input parameters:
        # For Rabi:
        """Input parameter:
        For Rabi
        min_len: Minimal length of the mw pulse applied (in us) ...usually that's zero (unchanged)
        start_distance: The distance between the end of the laserpulse and the start of the mw pulse (usual value: 300ns)
        stop_distance: The distance between the end of the last mw pulse and the start of the laser_re pulse (in us)
        steps: steps in us (unchanged)
        """
        self.mw_start_time_rabi = self.laser_in + self.mw_start_distance_rabi
        self.mw_max_len_rabi = round(
            ((self.laser_in + self.laser_off) - self.mw_stop_distance_rabi)
            - self.mw_start_time_rabi,
            ndigits=4,
        )
        """Input parameter: 
            For Ramsey
            min_distance: Minimal distance between the two pi/2 pulses (in us) ...usually that's zero (not changed)
            start_distance: The distance between the end of the laserpulse and the start of the mw pulse (usual value: 300ns)
            stop_distance: The distance between the end of the last mw pulse and the start of the laser_re pulse (in us)
            steps: steps in us (unchanged)
            mw_len_ramsey: length of the mw pulse itself (unchanged)
        """
        self.mw_start_time_ramsey = self.laser_in + self.mw_start_distance_ramsey
        self.mw_max_len_ramsey = (
            ((self.laser_in + self.laser_off) - self.mw_stop_distance_ramsey)
            - self.mw_start_time_ramsey
        ) - 2 * self.mw_len_ramsey

        # Build arrays from the StatusVar
        self.apd_times = [
            self.apd_start,
            self.apd_len,
        ]  # [time to start, length] all in microseconds
        self.apd_ref_times = [self.apd_ref_start, self.apd_ref_len]

        self.laser_times = [
            self.laser_in,
            self.laser_off,
            self.laser_re,
        ]  # on off on ...the rest is zeros
        # [mw_start_time, mw_len_0, mw_len_max, steps] all in microsecond
        self.mw_times_rabi = [
            self.mw_start_time_rabi,
            self.mw_min_len_rabi,
            self.mw_max_len_rabi,
            self.mw_steps_rabi,
        ]
        # [start, distance_min, distance_max, steps, pulse length] all in microseconds
        self.mw_times_ramsey = [
            self.mw_start_time_ramsey,
            self.mw_min_len_ramsey,
            self.mw_max_len_ramsey,
            self.mw_steps_ramsey,
            self.mw_len_ramsey,
        ]

        # mw_wait_time between laser off and mw on, mw_pulse_time
        self.mw_times_sweep = [self.mw_start_delay, self.mw_len_delay]
        self.apd_times_sweep = [
            self.apd_len_delay,
            self.apd_min_start_delay,
            self.apd_max_start_delay,
            self.apd_steps_delay,
        ]  # length, min_start, max_start, steps in microseconds

    def on_deactivate(self):
        self._pulser.clear_all()
        self._pulser.reset()
        self._pulser.pulser_off()
        # self.sigAverageDone.disconnect()
        self.sigContinueLoop.disconnect()
        self.sigStopMeasurement.disconnect()
        self.sigStopAll.disconnect()
        self.mw_off()

    @fc.constructor
    def sv_set_fits(self, val):
        # Setup fit container
        fc = self.fitlogic().make_fit_container("pulsed", "1d")
        fc.set_units(["us", "counts"])
        if isinstance(val, dict) and len(val) > 0:
            fc.load_from_dict(val)
        else:
            d1 = OrderedDict()
            d1["rabi"] = {
                "fit_function": "sineexponentialdecay",
                "estimator": "generic",
            }
            d1["ramsey"] = {
                "fit_function": "sinedoublewithtwoexpdecay",
                "estimator": "generic",
            }
            default_fits = OrderedDict()
            default_fits["1d"] = d1
            fc.load_from_dict(default_fits)
        return fc

    @fc.representer
    def sv_get_fits(self, val):
        """save configured fits"""
        if len(val.fit_list) > 0:
            return val.save_to_dict()
        else:
            return None

    def trigger(self):
        self._pulselogic.trigger()

    def cw(self, on=True):
        """
        If on is true: cw mode where the laser is on and the TTls direct the counts to the PFI 7 output on the daq card (both low/False)
        If on is False: all the digital outputs of the AWG are zero
        """
        self._pulser.reset()
        if on:
            self._pulselogic.play_ttl()
        else:
            self._pulselogic.play_ttl(
                seq_len=2, laser_out=False, apd_sig_out=False, apd_ref_out=False
            )

    def stop_awg(self):  # Makes everything stop
        self._pulselogic.stop_awg()

    def awg(self):
        """
        Setup and load awg
        clk_rate =  clock rate of the AWG in microseconds
        seq_len =   whole length of the sequence in microseconds
        laser_times = [laser_in, wait, laser_re] all in microseconds
                laser_in:   length of initialisation laser pulse (it starts right at the beginning) in microseconds
                wait:       time between the end of the initialisation and the reinitialisation in microseconds
                laser_re:   length of the reinitialisation pulse in microseconds
        Depending on the method:
        For Rabi
            apd_times = [time to start, length] all in microseconds
            mw_times = [mw_start_time, mw_len_0, mw_len_max, steps] all in microseconds
                mw_start_time:  time from the beginning when MW pulse should start, does not change
                mw_len_0:       minimal length of the mw pulse
                mw_len_max:     maximal length of the mw pulse
                steps:          stepsize for increasing the mw pulse duration
        For Ramsey:
            apd_times = [time to start, length of the pulse] both in microseconds
            mw_times = [start, distance_min, distance_max, steps, pulse length] all in microseconds
                start:          start position of the fist microwave pulse
                duration_min:   minimal distance between the two pulses
                duration_max:   maximal distance between the two pulses
                steps:          step size for increasing the mw pulse duration
                pulse length:   duration of pulses in microseconds (pi/2)
        For Delay sweep:
            apd_times = [length, min_start, max_start, steps]
                length:      time to wait until the pulse should start
                min_start:    minimal length of the apd pulse
                max_start:    maximal length of the apd pulse
                steps:      step size of the apd pulse sweep
            mw_times = [mw_wait_time, mw_pulse_time] all in microseconds
                mw_wait_time: time between the end of the first laser pulse and the beginning of the microwave (delay)
                mw_pulse_time: length of the mw pulse in between the two laser pulses
        apd_ref = [start of reference pulse, length of reference pulse]
        method =    can be rabi, ramsey or delaysweep
        rep =       repetition time when no trigger is used
        mw_pulse =  With or without MW pulse in the middle
        trigger =   either True: software trigger or False: no trigger at all
        """
        method = self.method
        if method == "rabi":
            # print('do rabi')
            if (
                len(self.apd_times) == 2
                and len(self.mw_times_rabi) == 4
                and len(self.laser_times) == 3
            ):
                apd_times = self.apd_times
                mw_times = self.mw_times_rabi
            else:
                self.log.warning("The input arrays don't have the right size.")
        elif method == "ramsey":
            if (
                len(self.apd_times) == 2
                and len(self.mw_times_ramsey) == 5
                and len(self.laser_times) == 3
            ):
                # print ('do ramsey')
                apd_times = self.apd_times
                mw_times = self.mw_times_ramsey
            else:
                self.log.warning("The input arrays don't have the right size.")
        elif method == "delaysweep":
            if (
                len(self.apd_times_sweep) == 4
                and len(self.mw_times_sweep) == 2
                and len(self.laser_times) == 3
            ):
                # print('do delay sweep')
                apd_times = self.apd_times_sweep
                mw_times = self.mw_times_sweep
            else:
                self.log.warning("The input arrays don't have the right size.")
        elif method == "delaysweep_ref":
            if (
                len(self.apd_times_sweep) == 4
                and len(self.mw_times_sweep) == 2
                and len(self.laser_times) == 3
            ):
                # print('do delay sweep_ref')
                apd_times = self.apd_times_sweep
                mw_times = self.mw_times_sweep
            else:
                self.log.warning("The input arrays don't have the right size.")
        else:
            self.log.warning(
                "The methode must be one of the following: delay_sweep, delay_sweep_ref, rabi, ramsey."
            )
        self.seq_len = self.calc_seq_len()
        method, laser_times, apd_times, apd_ref, mw_times = self._pulselogic.play_any(
            self.clk_rate_awg,
            self.seq_len,
            self.laser_times,
            apd_times,
            self.apd_ref_times,
            mw_times,
            self.mw_frequency,
            method=method,
            rep=self.rep,
            mw_pulse=self.mw_pulse_setting,
            trigger=self.trigger_setting,
            seq_map_req=self.seq_map_req,
            phase_shift=self.phase_shift,
        )

        self._sequence_map = self._pulselogic.get_seq_map()

        return method, laser_times, apd_times, apd_ref, mw_times

    def get_x_axis(self):  # in microseconds for the plot in the GUI
        steps, max_val, min_val = self.step_counter()
        x_axis = np.linspace(min_val, max_val, steps)
        return x_axis

    def get_t_axis(self):  # in microseconds for the instance plot in the GUI
        self.seq_len = self.calc_seq_len()
        number_samples = self.clk_rate_awg * self.seq_len  # this is MHz * microseconds
        t_axis = np.linspace(0, number_samples, number_samples)
        return t_axis

    def get_revant_parameters(self):
        method = self.method
        # laser_times = self.laser_array
        if method == "delaysweep":
            (
                laser_array,
                mw_i_array,
                mw_q_array,
                apd_array,
                apd_ref_array,
            ) = self.get_delay_sweep_arrays()
        elif method == "delaysweep_ref":
            (
                laser_array,
                mw_i_array,
                mw_q_array,
                apd_ref_array,
                apd_array,
            ) = self.get_delay_sweep_arrays()
        elif method == "rabi":
            (
                laser_array,
                mw_i_array,
                mw_q_array,
                apd_array,
                apd_ref_array,
            ) = self.get_rabi_arrays()
        elif method == "ramsey":
            (
                laser_array,
                mw_i_array,
                mw_q_array,
                apd_array,
                apd_ref_array,
            ) = self.get_ramsey_arrays()
        else:
            self.log.warning("method is weird.")
            laser_array = [0]
            apd_array = [0]
            apd_ref_array = [0]
            mw_i_array = [0]
            mw_q_array = [0]

        return method, laser_array, mw_i_array, mw_q_array, apd_array, apd_ref_array

    def get_used_arrays(self, analogs, digitals):
        if analogs != None and digitals != None:
            self.iwaveform = analogs["i_chan"]
            self.qwwaveform = analogs["q_chan"]
            self.laser_array = digitals["laser"]
            self.apd_array = digitals["apd_sig"]
            self.apd_ref_array = digitals["apd_read"]
        else:
            self.log.warning("Arrays are empty, run the sequence first.")
            t_axis = self.get_t_axis()
            self.iwaveform = np.repeat(0, len(t_axis))
            self.qwwaveform = np.repeat(0, len(t_axis))
            self.laser_array = np.repeat(0, len(t_axis))
            self.apd_array = np.repeat(0, len(t_axis))
            self.apd_ref_array = np.repeat(0, len(t_axis))
        return (
            self.iwaveform,
            self.qwwaveform,
            self.laser_array,
            self.apd_array,
            self.apd_ref_array,
        )

    def prepare_count_matrix(self):
        """
        creates a matrix with av_number of rows and step_count of columns
        """
        if self.method == "rabi":
            self.step_count = self._pulselogic.get_step_count(self.mw_times_rabi)
        elif self.method == "ramsey":
            self.step_count = self._pulselogic.get_step_count(self.mw_times_ramsey)
        else:
            self.step_count = self._pulselogic.get_step_count(self.apd_times_sweep)
        count_matrix = np.full((int(self.averages), int(self.step_count)), 0)
        count_matrix_ref = np.full((int(self.averages), int(self.step_count)), 0)
        # These are the index for the matrix (step_index: columns and av_index: rows)
        self.count_matrix = count_matrix
        self.count_matrix_ref = count_matrix_ref

        return count_matrix, count_matrix_ref

    def prepare_devices(self):
        """
        Initialize the counter and the settings of the MW device prior to starting scanning.
        """
        # Photon counter
        self._photon_samples = self._pxtime_to_samples()
        self._photon_counter.prepare_counters(
            samples_to_acquire=self._photon_samples, mode="pulsing"
        )

        # MW device: Put in the parameters and turn it on
        self._mw_device.set_power(self.mw_power)
        self._mw_device.set_frequency(self.mw_frequency)
        self.mw_on()
        self._mw_device.set_mod(True)

    def start_measurement(self):
        self._build_pulse_arrays()
        self.stop_awg()
        self.set_laser_power()  # this one talks to the awg via awg.cw
        self.prepare_count_matrix()
        self.prepare_devices()
        self.awg()
        self.sigSeqPlaying.emit(True)
        self.sigContinueLoop.emit()

    # TODO: original continue_loop in here!
    #  ###################

    # def continue_loop(self):
    #     #Play trigger
    #     # acquire the count value
    #     # write this value in the matrix
    #     # stop after 5 averages
    #     self.trigger() # sends a software trigger to the AWG
    #     if self.stopRequested == True:
    #         self.mw_off()
    #         self._mw_device.set_mod(False)
    #         self._photon_counter.close_counters()
    #         self.stop_awg()
    #         # make sure after all the inputs are low after the mesurement
    #         self.cw(False) # Pulls everything to zero
    #         # self._pulselogic.post_measurement(self.seq_len)
    #         self.av_index = 0
    #         self.step_index = 0
    #         self.stopRequested = False
    #         self.sigSeqPlaying.emit(False)
    #         self.log.info("Pulsed measurement stopped.")
    #     else:
    #         if self.step_index == self.step_count - 1 and self.av_index == self.averages -1:
    #             counts, ref_counts = self.acquire_pixel()
    #             # print(f'Current Average: {self.av_index}')
    #             self.count_matrix[self.av_index, self.step_index] = counts
    #             self.count_matrix_ref[self.av_index, self.step_index] = ref_counts
    #             self.sigStopMeasurement.emit()
    #
    #         elif self.step_index == self.step_count - 1 and self.av_index < self.averages - 1: # moves to the next average row
    #
    #             # print(f'Current Average: {self.av_index}')
    #             counts, ref_counts = self.acquire_pixel()
    #             # print(counts)
    #             self.count_matrix[self.av_index, self.step_index] = counts
    #             self.count_matrix_ref[self.av_index, self.step_index] = ref_counts
    #             self.step_index = 0
    #             idx = self.av_index + 1
    #             self.av_counts = self.count_matrix[:idx, :].mean(axis=0)
    #             self.av_counts_ref = self.count_matrix_ref[:idx, :].mean(axis=0)
    #             # ## Recalculate with the current seq_map
    #             # This part is only for the plotting in the gui
    #             # print(self.count_matrix[self.av_index])
    #             # if self.seq_map_req:
    #             #     current_row = self.recalc_from_seq_map(self.count_matrix[self.av_index])
    #             #     current_average = self.recalc_from_seq_map(self.av_counts)
    #             #     av_counts_ref = self.recalc_from_seq_map(self.av_counts_ref)
    #             # else: # Case: no special seq map
    #             current_row = self.count_matrix[self.av_index]
    #             current_average = self.av_counts
    #             av_counts_ref = self.av_counts_ref
    #             # print(current_row)
    #             self.sigAverageDone.emit(self.av_index, current_row, current_average, av_counts_ref)
    #             # print(current_row)
    #             self.av_index = self.av_index + 1
    #             self.sigContinueLoop.emit()
    #
    #         else:       # just continue in this --> direction
    #             counts, ref_counts = self.acquire_pixel()
    #             self.count_matrix[self.av_index, self.step_index] = counts
    #             self.count_matrix_ref[self.av_index, self.step_index] = ref_counts
    #             self.step_index = self.step_index + 1
    #
    #             self.sigContinueLoop.emit()
    #         print(self.count_matrix.mean(axis=0))
    #     return self.count_matrix
    # TODO: end of the original continue_loop

    def continue_loop(self):
        self.trigger()  # sends a software trigger to the AWG
        if not self.stopRequested:
            counts, ref_counts = self.acquire_pixel()

            if self._sequence_map:
                step = self._sequence_map[self.step_index]
            else:
                step = self.step_index

            self.count_matrix[self.av_index, step] = counts
            self.count_matrix_ref[self.av_index, step] = ref_counts

        self.sigPxAcquired.emit()
        return self.count_matrix

    def _next_pixel(self):
        self.step_index += 1
        if self.step_index == self.step_count:
            self.av_counts = self.count_matrix[: self.av_index + 1, :].mean(axis=0)
            self.av_counts_ref = self.count_matrix_ref[: self.av_index + 1, :].mean(
                axis=0
            )
            self.sigAverageDone.emit(
                self.av_index,
                self.count_matrix[self.av_index],
                self.av_counts,
                self.av_counts_ref,
            )

            self.av_index += 1
            self.step_index = 0

        if (self.av_index == self.averages) or self.stopRequested:

            self.stopRequested = True
            self.sigStopMeasurement.emit()
        else:
            self.sigContinueLoop.emit()

    def stop_measurement(self):
        self.mw_off()
        self._photon_counter.close_counters()
        self.stop_awg()  # awg is empty and ready for something new
        # make sure after all the inputs are low after the mesurement
        self.cw(False)
        # self._pulselogic.post_measurement(self.seq_len)
        self._mw_device.set_mod(False)
        self.av_index = 0
        self.step_index = 0
        self.stopRequested = False
        self.sigSeqPlaying.emit(False)
        self.log.info("Pulsed measurement stopped.")
        self.sigMeasurementDone.emit()

        #TODO: clear also the _sequence_map

    def save_data(self, method):
        timestamp = datetime.datetime.now()

        if method == "rabi":
            apd_times = self.apd_times
            mw_times = self.mw_times_rabi
        elif method == "ramsey":
            apd_times = self.apd_times
            mw_times = self.mw_times_ramsey
        else:
            apd_times = self.apd_times_sweep
            mw_times = self.mw_times_sweep

        data_raw = OrderedDict()
        data_raw["count data"] = self.count_matrix
        data_raw["REF count data"] = self.count_matrix_ref
        data_raw["averaged signal"] = self.av_counts

        parameters = OrderedDict()
        parameters["Microwave Power"] = self.mw_power
        parameters["Run Time"] = self.elapsed_time
        parameters["Frequency"] = self.mw_frequency
        parameters["Number of steps"] = self.step_count
        parameters["Number of averages"] = self.averages
        parameters["integration time"] = self.integration_time
        parameters["Clock rate of DAQ"] = self.clk_rate_daq
        parameters["Clock rate of AWG"] = self.clk_rate_awg
        parameters["Laser times"] = self.laser_times
        parameters["APD times"] = apd_times
        parameters["MW times"] = mw_times
        parameters["Used method"] = method
        parameters["Rep"] = self.rep
        parameters["MW pulse used?"] = self.mw_pulse_setting
        parameters["software trigger"] = self.trigger_setting
        parameters["laser power"] = self.laser_power
        parameters["extra notes"] = self.plaintext_field
        parameters["timestamp"] = timestamp.strftime("%Y-%m-%d, %H:%M:%S")

        attributes = {"count data": parameters, "REF count data": parameters}

        # Check if a fit has been performed
        if method in self.fits_performed:
            x_fit, y_fit, result = self.fits_performed[method]
            data_raw["fit values"] = np.stack((x_fit, y_fit))
            attributes["fit values"] = result.values

        self._savelogic.save_hdf5_data(
            data_raw,
            attributes=attributes,
            filename=timestamp.strftime("%Y-%m-%d-%H-%M-%S")
            + "_"
            + self.method
            + ".h5",
        )

    # TODO: this is the beginning of the original save_data
    # def save_data(self, method):
    #     timestamp = datetime.datetime.now()
    #     ## The matrix should always be saved in the right order! Use the right sequence map for every row in the matrix
    #     # This takes the seq_map into account:
    #     if self.seq_map_req:
    #         self.count_matrix = self.recalc_matrix(self.count_matrix)
    #
    #     if method == "rabi":
    #         apd_times = self.apd_times
    #         mw_times = self.mw_times_rabi
    #     elif method == "ramsey":
    #         apd_times = self.apd_times
    #         mw_times = self.mw_times_ramsey
    #     else:
    #         apd_times = self.apd_times_sweep
    #         mw_times = self.mw_times_sweep
    #
    #     data_raw = OrderedDict()
    #     data_raw["count data"] = self.count_matrix
    #     data_raw["REF count data"] = self.count_matrix_ref
    #     data_raw["averaged signal"] = self.av_counts
    #
    #     parameters = OrderedDict()
    #     parameters["Microwave Power"] = self.mw_power
    #     parameters["Run Time"] = self.elapsed_time
    #     parameters["Frequency"] = self.mw_frequency
    #     parameters["Number of steps"] = self.step_count
    #     parameters["Number of averages"] = self.averages
    #     parameters["integration time"] = self.integration_time
    #     parameters["Clock rate of DAQ"] = self.clk_rate_daq
    #     parameters["Clock rate of AWG"] = self.clk_rate_awg
    #     parameters["Laser times"] = self.laser_times
    #     parameters["APD times"] = apd_times
    #     parameters["MW times"] = mw_times
    #     parameters["Used method"] = method
    #     parameters["Rep"] = self.rep
    #     parameters["MW pulse used?"] = self.mw_pulse_setting
    #     parameters["software trigger"] = self.trigger_setting
    #     parameters["laser power"] = self.laser_power
    #     parameters["extra notes"] = self.plaintext_field
    #     parameters["timestamp"] = timestamp.strftime("%Y-%m-%d, %H:%M:%S")
    #
    #     attributes = {"count data": parameters, "REF count data": parameters}
    #
    #     # Check if a fit has been performed
    #     if method in self.fits_performed:
    #         x_fit, y_fit, result = self.fits_performed[method]
    #         data_raw["fit values"] = np.stack((x_fit, y_fit))
    #         attributes["fit values"] = result.values
    #
    #     self._savelogic.save_hdf5_data(
    #         data_raw,
    #         attributes=attributes,
    #         filename=timestamp.strftime("%Y-%m-%d-%H-%M-%S")
    #         + "_"
    #         + self.method
    #         + ".h5",
    #     )
    # TODO: this is the end of the original saving

    def acquire_pixel(self):
        # Requests the counts from the DAQ card, gives one number only
        samples = int(self.integration_time * self.clk_rate_daq)
        data, _ = self._photon_counter.read_pixel(samples)
        counts, ref_counts = data
        return counts, ref_counts

    def mw_off(self):
        """Switching off the MW source.

        @return str, bool: active mode ['cw', 'list', 'sweep'], is_running
        """
        error_code = self._mw_device.off()
        if error_code < 0:  # This means there is an error
            self.log.warning("Switching off microwave source failed.")
        else:
            self.log.info("Pulsed: MW device is off")
        return

    def mw_on(self):
        """
        Switching on the mw source.
        Parameters must be set somewhere else
        """
        error_code = self._mw_device.on()

        if error_code < 0:  # This means there is an error
            self.log.warning("Switching on microwave source failed.")
        else:
            self.log.info("Pulsed: MW device is on")
        return

    def _pxtime_to_samples(self):
        return round(
            self.integration_time * self._photon_counter.get_counter_clock_frequency()
        )

    def stop_all(self):
        self.stopRequested = True

    def step_counter(self):
        self._build_pulse_arrays()
        # in microseconds
        if self.method == "rabi":
            step_count = self._pulselogic.get_step_count(self.mw_times_rabi)
            max_val = self.mw_times_rabi[2]
            min_val = self.mw_times_rabi[1]
        elif self.method == "ramsey":
            step_count = self._pulselogic.get_step_count(self.mw_times_ramsey)
            max_val = self.mw_times_ramsey[2]
            min_val = self.mw_times_ramsey[1]
        else:
            step_count = self._pulselogic.get_step_count(self.apd_times_sweep)
            max_val = self.apd_times_sweep[2]
            min_val = self.apd_times_sweep[1]

        return step_count, max_val, min_val

    def set_laser_power(self):
        if self.laser_power > 0.95:  # for safety reasons
            self._laser.set_voltage(0.95)
        else:
            self._laser.set_voltage(self.laser_power)

    def do_fit(self, fit_function=None):
        """
        Execute the currently configured fit on the measurement data. Optionally on passed data
        """
        # load x_data and y_data
        x_data = self.get_x_axis()
        y_data = self.av_counts

        if isinstance(fit_function, str):
            if fit_function in self.get_fit_functions():
                self.fc.set_current_fit(fit_function)
            else:
                self.fc.set_current_fit("No Fit")
                if fit_function != "No Fit":
                    self.log.warning(
                        'Fit function "{0}" not available in ODMRLogic fit container.'
                        "".format(fit_function)
                    )

        self.fit_x, self.fit_y, result = self.fc.do_fit(x_data, y_data)

        if fit_function != "No Fit":
            self.fits_performed[self.fc.current_fit] = (self.fit_x, self.fit_y, result)

        if result is None:
            result_str_dict = {}
        else:
            result_str_dict = result.result_str_dict

        # print('emitting update fit')
        self.sigFitUpdated.emit(
            self.fit_x, self.fit_y, result_str_dict, self.fc.current_fit
        )
        return

    def get_fit_functions(self):
        """Return the hardware constraints/limits
        @return list(str): list of fit function names
        """
        return list(self.fc.fit_list)

    def start_stop_timetrace(self, val):
        if val == True:
            # Start timetrace
            self._timeseries.start_reading()
        else:
            self._timeseries.stop_reading()

    # just to get the right parameter for instance plot
    def get_rabi_arrays(self):
        """
        It just creates the arrays. They are used in the masterpulse logic and the GUI to build the instance plot on the right
        """
        self._build_pulse_arrays()
        mw_len = []
        card_idx = 1
        self.clk_rate = (
            self._pulser.get_sample_rate(card_idx) / 1e6
        )  # To adapt it to the microsecond timing
        step_count = int(
            ((self.mw_times_rabi[2] - self.mw_times_rabi[1]) / self.mw_times_rabi[3])
            + 1
        )
        for i in range(step_count):
            mw_len.append(
                round((self.mw_times_rabi[1] + (i * self.mw_times_rabi[3])), 3)
            )  # starting with minimal length

        t_centrewidth_list_mw = (
            []
        )  # get the size of the np.array right?np.array([[0, 0]])
        t_centrewidth_list_laser = self._pulselogic.convert_laser_val(self.laser_times)
        t_centrewidth_list_apd = self._pulselogic.convert_apd_val(self.apd_times)
        t_centrewidth_list_apd_ref = self._pulselogic.convert_apd_val(
            self.apd_ref_times
        )
        self.seq_len = self.calc_seq_len()
        for i in mw_len:
            if self.seq_len < (self.mw_times_rabi[0] + self.mw_times_rabi[2]):
                raise ValueError(
                    "The total length needs to be larger that the apd_time sum"
                )
            else:
                t0, width = self._pulselogic.convert_mw_val_rabi(
                    self.mw_times_rabi, i
                )  # in microseconds
                t_centrewidth_list_mw.append([np.array([[t0, width]])])
        self._pulser.set_frequency(self.mw_frequency)
        output_frequency = self._pulser.get_frequency()
        time_ax = self._pulselogic.time_axis(
            self.seq_len
        )  # this includes waveform_padding
        laser_array = []
        mw_i_array = []
        mw_q_array = []
        apd_array = []
        apd_ref_array = []
        for step in range(len(t_centrewidth_list_mw)):
            iwaveform, qwwaveform = self._pulser.iq_pulses(
                time_ax, t_centrewidth_list_mw[step], output_frequency
            )
            laser_waveform = self._pulser.box_envelope(
                time_ax, t_centrewidth_list_laser
            )
            apd_waveform = self._pulser.box_envelope(time_ax, t_centrewidth_list_apd)  #
            apd_ref_waveform = self._pulser.box_envelope(
                time_ax, t_centrewidth_list_apd_ref
            )  #
            mw_i_array.append(iwaveform)
            mw_q_array.append(qwwaveform)
            laser_array.append(laser_waveform)
            apd_array.append(apd_waveform)
            apd_ref_array.append(apd_ref_waveform)
        # print(laser_array[1])
        # self.sigUsedArrays.emit(analogs, digitals)
        # This does the right thing and the values are updated
        return laser_array, mw_i_array, mw_q_array, apd_array, apd_ref_array

    def get_ramsey_arrays(self):
        """
        It just creates the arrays. They are used in the masterpulse logic and the GUI to build the instance plot on the right
        """
        self._build_pulse_arrays()
        t_centrewidth_list_mw = []
        break_len = []
        laser_array = []
        mw_i_array = []
        mw_q_array = []
        apd_array = []
        apd_ref_array = []
        card_idx = 1
        step_count = int(
            (
                (self.mw_times_ramsey[2] - self.mw_times_ramsey[1])
                / self.mw_times_ramsey[3]
            )
            + 1
        )  # without unit

        for i in range(step_count):
            break_len.append(
                (self.mw_times_ramsey[1] + (i * self.mw_times_ramsey[3]))
            )  # starting with minimal length
        t_centrewidth_list_laser = self._pulselogic.convert_laser_val(self.laser_times)
        t_centrewidth_list_apd = self._pulselogic.convert_apd_val(self.apd_times)
        t_centrewidth_list_apd_ref = self._pulselogic.convert_apd_val(
            self.apd_ref_times
        )
        self.seq_len = self.calc_seq_len()
        for i in break_len:
            if self.seq_len < (
                self.laser_times[0] + self.laser_times[1] + self.laser_times[2]
            ):
                self.log.warning(
                    "The total length needs to be larger that the sequence time"
                )
            else:
                t0, t1, width = self._pulselogic.convert_mw_val_ramsey(
                    self.mw_times_ramsey, i
                )  # these values are in microseconds
                t_centrewidth_list_mw.append([np.array([[t0, width], [t1, width]])])
        time_ax = self._pulselogic.time_axis(self.seq_len)  # in microsecondes
        self._pulser.set_frequency(self.mw_frequency)  # this one does freq - if_freq
        output_frequency = self._pulser.get_frequency()
        if self.phase_shift:
            phases = np.empty((len(t_centrewidth_list_mw)))
            phases[::2] = 0
            phases[1::2] = np.pi
        else:
            phases = np.repeat(0, len(t_centrewidth_list_mw))
        for step in range(len(t_centrewidth_list_mw)):
            iwaveform, qwwaveform = self._pulser.iq_pulses(
                time_ax, t_centrewidth_list_mw[step], output_frequency, phases
            )
            laser_waveform = self._pulser.box_envelope(
                time_ax, t_centrewidth_list_laser
            )
            apd_waveform = self._pulser.box_envelope(time_ax, t_centrewidth_list_apd)  #
            apd_ref_waveform = self._pulser.box_envelope(
                time_ax, t_centrewidth_list_apd_ref
            )  #
            mw_i_array.append(iwaveform)
            mw_q_array.append(qwwaveform)
            laser_array.append(laser_waveform)
            apd_array.append(apd_waveform)
            apd_ref_array.append(apd_ref_waveform)

        # self.sigUsedArrays.emit(analogs, digitals)
        return laser_array, mw_i_array, mw_q_array, apd_array, apd_ref_array

    def get_delay_sweep_arrays(self):
        """
        It just creates the arrays. They are used in the masterpulse logic and the GUI to build the instance plot on the right
        """
        self._pulser.set_sample_rate(1, int(800000000))
        laser_array = []
        mw_i_array = []
        mw_q_array = []
        apd_array = []
        apd_ref_array = []
        # seq_len_s = self.seq_len * 1e-6
        # laser_times_s = np.multiply(self.laser_times, 1e-6)
        # mw_times_s = np.multiply(self.mw_times_sweep, 1e-6)
        apd_times_s = np.multiply(self.apd_times_sweep, 1e-6)
        # apd_ref_s = np.multiply(self.apd_ref_times, 1e-6)
        # So far we just use card1
        card_idx = 1
        # clk_rate = self._pulser.get_sample_rate(card_idx)
        # segment_map = []
        apd_start = []
        step_count = int(
            (
                (self.apd_times_sweep[2] - self.apd_times_sweep[1])
                / self.apd_times_sweep[3]
            )
            + 1
        )

        self.seq_len = self.calc_seq_len()
        time_ax = self._pulselogic.time_axis(
            self.seq_len
        )  # this includes waveform_padding
        for i in range(step_count):
            apd_start.append((apd_times_s[1] + (i * apd_times_s[3])))

        # Create laser and apd waveform as they are constant (put it in microseconds)
        t_centrewidth_list_laser = self._pulselogic.convert_laser_val(self.laser_times)
        t_centrewidth_list_apd_ref = self._pulselogic.convert_apd_val(
            self.apd_ref_times
        )
        laser_waveform = self._pulser.box_envelope(time_ax, t_centrewidth_list_laser)
        apd_ref_waveform = self._pulser.box_envelope(
            time_ax, t_centrewidth_list_apd_ref
        )
        if self.mw_pulse_setting:
            # replace this one
            t_centrewidth_list_mw = self._pulselogic.convert_mw_sweep(
                self.mw_times_sweep
            )
            mw_waveform = self._pulser.box_envelope(time_ax, t_centrewidth_list_mw)

        else:  # if no MW is needed
            zeros_mw = self._pulselogic.convert_mw_sweep([-1, 0])
            mw_waveform = self._pulser.box_envelope(time_ax, zeros_mw)

        for i in apd_start:
            current_t_centerwidth = self._pulselogic.convert_apd_val(
                [i * 1e6, self.apd_times_sweep[0]]
            )  # [time to start, length]
            # print(current_t_centerwidth)
            #!!! apd_times = [length, min_start, max_start, steps] # changing length
            apd_waveform = self._pulser.box_envelope(time_ax, current_t_centerwidth)
            mw_i_array.append(mw_waveform)
            mw_q_array.append(np.zeros(len(mw_waveform)))
            laser_array.append(laser_waveform)
            apd_array.append(apd_waveform)
            apd_ref_array.append(apd_ref_waveform)

        return laser_array, mw_i_array, mw_q_array, apd_array, apd_ref_array

    def get_t_axis_pulselogic(self):  # in microseconds for the instance plot in the GUI
        self.seq_len = self.calc_seq_len()
        t_axis = self._pulselogic.time_axis(self.seq_len)
        return t_axis

    def len_errors(self, mw_val, clk_rate):
        # this is general...mw_val is just any relevant value. The separation between the different methods happens in the GUI
        min_time = round(
            1 / clk_rate, 5
        )  # clock_rate in MHz ...min time is in ns therefore
        req_len = round(mw_val, 5)
        if (
            round(req_len % min_time, 5) == 0
            or round(req_len % min_time, 5) == min_time
        ):
            # if it works for the length of one step and we start at 0 length then it works for all the steps.
            real_len = round(req_len, 5)
        else:
            # The real value will be different.
            sample_number = round(req_len / min_time, 0)
            real_len = round(sample_number * min_time, 5)
        # print('real_len:', real_len)
        return self.method, min_time, req_len, real_len

    def get_seq_map(self):
        """
        fishes the seq_map and the step_count from the pulse_logic
        """
        if self.method == "rabi":
            self.step_count = self._pulselogic.get_step_count(self.mw_times_rabi)
        elif self.method == "ramsey":
            self.step_count = self._pulselogic.get_step_count(self.mw_times_ramsey)
        else:
            self.step_count = self._pulselogic.get_step_count(self.apd_times_sweep)
        self.seq_map = self._pulselogic.build_seq_map(self.step_count, self.seq_map_req)

        return self.step_count, self.seq_map

    def recalc_from_seq_map(self, array):
        """
        array = the data array in the order of the seq_map
        seq_map = gives the oder of the data points
        # what happens if seq_map = []

        """
        step_count, seq_map = self.get_seq_map()
        array_new = np.zeros(len(seq_map))
        for n in range(step_count):
            array_new[seq_map[n]] = array[n]
        # print('seq_map', seq_map)
        # print('array', array)
        # print('array_new', array_new)
        return array_new

    def calc_seq_len(self):
        self.seq_len = self.laser_in + self.laser_off + self.laser_re
        if self.method == "delaysweep" or self.method == "delaysweep_ref":
            self.seq_len = self.laser_in + self.laser_off + self.laser_re + 1
        return self.seq_len

    # def recalc_seq_map_test(self, array, seq_map, step_count):
    #     """
    #     array = the data array in the order of the seq_map
    #     seq_map = gives the oder of the data points
    #     # what happens if seq_map = []
    #     """
    #     # step_count, seq_map = self.get_seq_map()
    #     array_new = np.zeros(len(seq_map))
    #     for n in range(step_count):
    #         array_new[seq_map[n]] = array[n]
    #     # print('seq_map', seq_map)
    #     # print('array', array)
    #     #     print('array_new', array_new)
    #     return array_new
    def recalc_matrix(self, count_matrix):
        ## It is improtant that the sequence map starts with 0!
        rows = np.size(count_matrix, 0)
        count_maxtrix_new = []
        seq_map = [0, 4, 1, 3, 2]
        step_count = len(seq_map)
        print(rows)
        for n in range(rows):
            nrow = self.recalc_from_seq_map(count_matrix[n])
            count_maxtrix_new.append(nrow)
        return count_maxtrix_new
