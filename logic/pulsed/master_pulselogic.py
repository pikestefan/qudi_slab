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
from logic.generic_logic import GenericLogic



# Just hardcode the channels to the different hardware modules somewehere
# Microwave: analog channel 2 and 3 for I and Q parameter. For now we just need I
# Laser: digital channel X0
# Photon counter: digital channel X1 and X2 (for reference)

class Masterpulse(GenericLogic):
    pulsegenerator = Connector(interface='PulserInterface')
    savelogic = Connector(interface='SaveLogic')
    pulselogic = Connector(interface='Pulse')
    photon_counter = Connector(interface='SnvmScannerInterface')
    mw_source = Connector(interface='MicrowaveInterface')


    # integration_time = 30e-3 # In ms for timing...
    #
    # Internal signals
    sigContinueLoop = QtCore.Signal()
    sigStopMeasurement = QtCore.Signal()
    sigStopAll = QtCore.Signal()
    sigNextLine = QtCore.Signal()
    #
    # # Update signals, e.g. for GUI module
    sigParameterUpdated = QtCore.Signal(dict)
    sigOutputStateUpdated = QtCore.Signal(str, bool)
    sigElapsedTimeUpdated = QtCore.Signal(float, int)

    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)
        # self.threadlock = Mutex()

    def on_activate(self):
        # Initialisation performed during activation of the module

        # Get connectors
        self._pulser = self.pulsegenerator()
        self._photon_counter = self.photon_counter()
        self._mw_device = self.mw_source()
        self._savelogic = self.savelogic()
        self._pulselogic = self.pulselogic()
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
        self.step_index = 0 # value right now, column
        self.averages = 3 # max value for rows
        self.mw_frequency = 2870e6 # This should be the slope frequency?
        self.mw_power = -30
        self.integration_time = 30e-3  # In ms
        # self.step_count = self._pulselogic.get_step_count([2.1, 1, 2, 0.5])  #max value for columns
        self.clk_rate_daq = self._photon_counter.get_counter_clock_frequency() # which is 100000
        self.clk_rate_awg = 50  # in ms
        self.seq_len = 7       # in ms
        self.laser_times = [1.4, 2.1, 3.2] # on off on ...the rest is zeros
        self.mw_times_rabi = [1.5, 0.1, 2.0, 0.1] # [mw_start_time, mw_len_0, mw_len_max, steps] all in ms
        self.mw_times_ramsey = [1.45, 0.05, 0.8, 0.01, 0.01] #[start, distance_min, distance_max, steps, pulse length] all in ms
        self.mw_times_sweep = [0.06, 0.2] # mw_wait_time between laser off and mw on, mw_pulse_time
        self.apd_times = [3.51, 0.2] #[time to start, length] all in ms
        self.apd_times_sweep = [0.1, 3.25, 4.2, 0.1] #length, min_start, max_start, steps
        self.apd_ref_times = [5.1, 0.1]
        self.rep = 100
        self.mw_pulse_setting = True # this is only for delay sweep
        self.trigger_setting = True
        self.method = 'delaysweep'    # change this to 'ramsey' or 'delaysweep' if needed
        # Declare where the signals should lead to
        self.sigContinueLoop.connect(self.continue_loop, QtCore.Qt.QueuedConnection)
        self.sigStopMeasurement.connect(self.stop_measurement, QtCore.Qt.QueuedConnection)
        self.sigStopAll.connect(self.stop_all, QtCore.Qt.QueuedConnection)
        self.stopRequested = False

        return active_channels

    def on_deactivate(self):
        self._pulser.clear_all()
        self._pulser.reset()
        self._pulser.pulser_off()
        self.mw_off()

    def trigger(self):
        self._pulselogic.trigger()

    def stop_awg(self):  # Makes everything stop
        self._pulselogic.stop_awg()

    def awg(self):
        '''
        Setup and load awg
        clk_rate =  clock rate of the AWG in ms
        seq_len =   whole length of the sequence in ms
        laser_times = [laser_in, wait, laser_re] all in ms
                laser_in:   length of initialisation laser pulse (it starts right at the beginning) in ms
                wait:       time between the end of the initialisation and the reinitialisation in ms
                laser_re:   length of the reinitialisation pulse in ms
        Depending on the method:
        For Rabi
            apd_times = [time to start, length] all in ms
            mw_times = [mw_start_time, mw_len_0, mw_len_max, steps] all in ms
                mw_start_time:  time from the beginning when MW pulse should start, does not change
                mw_len_0:       minimal length of the mw pulse
                mw_len_max:     maximal length of the mw pulse
                steps:          stepsize for increasing the mw pulse duration
        For Ramsey:
            apd_times = [time to start, length of the pulse] both in ms
            mw_times = [start, distance_min, distance_max, steps, pulse length] all in ms
                start:          start position of the fist microwave pulse
                duration_min:   minimal distance between the two pulses
                duration_max:   maximal distance between the two pulses
                steps:          step size for increasing the mw pulse duration
                pulse length:   duration of pulses in ms (pi/2)
        For Delay sweep:
            apd_times = [length, min_start, max_start, steps]
                length:      time to wait until the pulse should start
                min_start:    minimal length of the apd pulse
                max_start:    maximal length of the apd pulse
                steps:      step size of the apd pulse sweep
            mw_times = [mw_wait_time, mw_pulse_time] all in ms
                mw_wait_time: time between the end of the first laser pulse and the beginning of the microwave (delay)
                mw_pulse_time: length of the mw pulse in between the two laser pulses
        apd_ref = [start of reference pulse, length of reference pulse]
        method =    can be rabi, ramsey or delaysweep
        rep =       repetition time when no trigger is used
        mw_pulse =  With or without MW pulse in the middle
        trigger =   either True: software trigger or False: no trigger at all
        '''
        method = self.method

        if method == 'rabi':
            # print ('do rabi')
            if len(self.apd_times) == 2 and len(self.mw_times_rabi) == 4 and len(self.laser_times) == 3:
                apd_times = self.apd_times
                mw_times = self.mw_times_rabi
            else:
                self.log.error("The input arrays don't have the right size.")
        elif method == 'ramsey':
            if len(self.apd_times) == 2 and len(self.mw_times_ramsey) == 5 and len(self.laser_times) == 3:
                # print ('do ramsey')
                apd_times = self.apd_times
                mw_times = self.mw_times_ramsey
            else:
                self.log.error("The input arrays don't have the right size.")
        else:
            if len(self.apd_times_sweep) == 4 and len(self.mw_times_sweep) == 2 and len(self.laser_times) == 3:
                print ('do delay sweep')
                apd_times = self.apd_times_sweep
                mw_times = self.mw_times_sweep
            else:
                # print (len(self.apd_times_sweep))
                # print(len(self.mw_times_sweep))
                # print (len(self.mw_times_sweep))
                self.log.error("The input arrays don't have the right size.")

        self._pulselogic.play_any(self.clk_rate_awg, self.seq_len, self.laser_times, apd_times, self.apd_ref_times, mw_times, method=method,
                                  rep=self.rep, mw_pulse=self.mw_pulse_setting, trigger=self.trigger_setting)
        return method

    def prepare_count_matrix(self):
        '''
        creates a matrix with av_number of rows and step_count of columns

        '''
        if self.method == 'rabi':
            self.step_count = self._pulselogic.get_step_count(self.mw_times_rabi)
        elif self.method == 'ramsey':
            self.step_count = self._pulselogic.get_step_count(self.mw_times_ramsey)
        else:
            self.step_count = self._pulselogic.get_step_count(self.apd_times_sweep)
        count_matrix = np.full((self.averages, self.step_count), 0)
        count_matrix_ref = np.full((self.averages, self.step_count), 0)
        # These are the index for the matrix (step_index: columns and av_index: rows)
        self.count_matrix = count_matrix
        self.count_matrix_ref = count_matrix_ref
        # print(count_matrix)
        return count_matrix, count_matrix_ref
    #
    def prepare_devices(self):
        """
        Initialize the counter and the settings of the MW device prior to starting scanning.
        """
        # Photon counter
        self._photon_samples = self._pxtime_to_samples()
        self._photon_counter.prepare_counters(samples_to_acquire=self._photon_samples)

        # MW device: Put in the parameters and turn it on
        self._mw_device.set_power(self.mw_power)
        self._mw_device.set_frequency(self.mw_frequency)
        self.mw_on()

    def start_measurement(self):
        self.stop_awg()
        self.prepare_count_matrix()
        self.prepare_devices()
        self.awg()
        self.sigContinueLoop.emit()

    def continue_loop(self):
        #Play trigger
        # acquire the count value
        # write this value in the matrix
        # stop after 5 averages
        self.trigger() # sends a software trigger to the AWG
        time.sleep(0.1) # in sec
        if self.stopRequested == True:
            self.mw_off()
            self._photon_counter.close_counters()
            self.stop_awg()
            self.av_index = 0
            self.step_index = 0
            self.stopRequested = False
            print('measurement stopped')
            # It does not save anything
        else:
            if self.step_index == self.step_count - 1 and self.av_index == self.averages - 1:
                # print('very last value')
                counts = self.acquire_pixel()
                ref_counts = self.acquire_pixel()
                print ('counts: ', counts)
                self.count_matrix[self.av_index, self.step_index] = counts
                self.count_matrix_ref[self.av_index, self.step_index] = ref_counts
                self.sigStopMeasurement.emit()
            elif self.step_index == self.step_count - 1 and self.av_index < self.averages - 1: # moves to the next average row
                # print('last value in the row')
                counts = self.acquire_pixel()
                ref_counts = self.acquire_pixel()
                print('counts: ', counts)
                self.count_matrix[self.av_index, self.step_index] = counts
                self.count_matrix_ref[self.av_index, self.step_index] = ref_counts
                self.step_index = 0
                self.av_index = self.av_index + 1
                self.sigContinueLoop.emit()
            else:       # just continue in this --> direction
                # print('switch to the next column, sequence')
                counts = self.acquire_pixel()
                ref_counts = self.acquire_pixel()
                print('counts: ', counts)
                self.count_matrix[self.av_index, self.step_index] = counts
                self.count_matrix_ref[self.av_index, self.step_index] = ref_counts
                self.step_index = self.step_index + 1
                self.sigContinueLoop.emit()

    def stop_measurement(self):
        self.mw_off()
        self._photon_counter.close_counters()
        self.stop_awg()
        self.av_index = 0
        self.step_index = 0
        self.stopRequested = False
        print('measurement done')
        #SAVE?
        self.save_data(self.method)

    def save_data(self, method):

        timestamp = datetime.datetime.now()
        if method == 'rabi':
            apd_times = self.apd_times
            mw_times = self.mw_times_rabi
        elif method == 'ramsey':
            apd_times = self.apd_times
            mw_times = self.mw_times_ramsey
        else:
            apd_times = self.apd_times_sweep
            mw_times = self.mw_times_sweep

        data_raw = OrderedDict()
        data_raw['count data'] = self.count_matrix
        data_raw['REF count data'] = self.count_matrix_ref
        parameters = OrderedDict()
        parameters['Microwave Power'] = self.mw_power
        parameters['Run Time'] = self.elapsed_time
        parameters['Frequency'] = self.mw_frequency
        parameters['Number of steps'] = self.step_count
        parameters['Number of averages'] = self.averages
        parameters['integration time'] = self.integration_time
        parameters['Clock rate of DAQ'] = self.clk_rate_daq
        parameters['Clock rate of AWG'] = self.clk_rate_awg
        parameters['Laser times'] = self.laser_times
        parameters['APD times'] = apd_times
        parameters['MW times'] = mw_times
        parameters['Used method'] = method
        parameters['Rep'] = self.rep
        parameters['MW pulse used?'] = self.mw_pulse_setting
        parameters['software trigger'] = self.trigger_setting
        parameters['timestamp'] = timestamp.strftime('%Y-%m-%d, %H:%M:%S')
        self._savelogic.save_hdf5_data(data_raw, attributes={'count data': parameters, 'REF count data': parameters}, filename= timestamp.strftime('%Y-%m-%d-%H-%M-%S') + '_' + self.method +'.h5')


    def acquire_pixel(self):
        # Requests the counts from the DAQ card, gives one number only
        x = int(self.integration_time * self.clk_rate_daq)
        counts, _ = self._photon_counter.read_pixel(x)
        ref_counts, _ = self._photon_counter.read_pixel(x)
        counts = counts[0] # to just get the value and not the array
        ref_counts = ref_counts[0]  # to just get the value and not the array
        return counts

    def mw_off(self):
        """ Switching off the MW source.

        @return str, bool: active mode ['cw', 'list', 'sweep'], is_running
        """
        error_code = self._mw_device.off()
        if error_code < 0:  # This means there is an error
            self.log.error('Switching off microwave source failed.')
        else:
            print('MW device is off')
        return

    def mw_on(self):
        """
        Switching on the mw source.
        Parameters must be set somewhere else
        """
        error_code = self._mw_device.on()
        if error_code < 0:  # This means there is an error
            self.log.error('Switching on microwave source failed.')
        else:
            print('MW device is on')
        return

    def _pxtime_to_samples(self):
        return round(self.integration_time * self._photon_counter.get_counter_clock_frequency())

    def stop_all(self):
        self.stopRequested = True
