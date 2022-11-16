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
    pulsegenerator = Connector(interface='PulserInterface')
    laser = Connector(interface='SimpleLaserInterface')
    photon_counter = Connector(interface='SnvmScannerInterface')
    mw_source = Connector(interface='MicrowaveInterface')
    savelogic = Connector(interface='SaveLogic')
    pulselogic = Connector(interface='Pulse')
    fitlogic = Connector(interface='FitLogic')
    timeseries_logic = Connector(interface='TimeSeriesReaderLogic')

    ## Define status variables
    fc = StatusVar('fits', None)
    # status vars
    averages = StatusVar('averages', 10)  # max value for rows
    mw_frequency = StatusVar('mw_frequency', 2.8685e9)
    mw_power = StatusVar('mw_power', -30)  # in dBm
    integration_time = StatusVar('integration_time', 30e-3)  # in s
    laser_power = StatusVar('laser_power', 0.628)  # in V... Max is 1V!
    clk_rate_awg = StatusVar('clk_rate_awg', int(800e6))  # in MHz

    ## Delay Sweep
    mw_pulse_setting = StatusVar('mw_pulse_setting', True)  # this is only for delay sweep
    seq_len = StatusVar('seq_len', 8.5e-6)
    neg_mw = StatusVar('neg_mw', False)  # if False it plays +1, if True it applies negative mw pulses -1 (all the pulses are squared)
    method = StatusVar('method', 'rabi')  # change this to 'ramsey', 'delaysweep', 'rabi' or 'delaysweep_ref' if needed

    # StatusVar for arrays:
    apd_start = StatusVar('apd_start', 3.05e-6)
    apd_len = StatusVar('apd_len', 250e-9)
    apd_ref_start = StatusVar('apd_ref_start', 4.5e-6)
    apd_ref_len = StatusVar('apd_ref_len', 250e-9)

    # laser times tab
    laser_in = StatusVar('laser_in', 1e-6)
    laser_off = StatusVar('laser_off', 2e-6)
    laser_re = StatusVar('laser_re', 2e-6)

    # delay_sweep tab
    apd_len_delay = StatusVar('apd_len_delay', 50e-9)
    apd_min_start_delay = StatusVar('apd_min_start_delay', 0e-6)
    apd_max_start_delay = StatusVar('apd_max_start_delay', 2.6e-6)
    apd_steps_delay = StatusVar('apd_steps_delay', 50e-9)

    mw_start_delay = StatusVar('mw_start_delay', 1.3e-6)
    mw_len_delay = StatusVar('mw_len_delay', 500e-9) #What ever the pi pulse is

    # rabi tab
    mw_start_time_rabi = StatusVar('mw_start_time_rabi', 1.3e-6)
    mw_min_len_rabi = StatusVar('mw_min_len_rabi', 0e-6)
    mw_max_len_rabi = StatusVar('mw_max_len_rabi', 900e-9)
    mw_steps_rabi = StatusVar('mw_steps_rabi', 100e-9)

    # ramsey tab
    mw_start_time_ramsey = StatusVar('mw_start_time_ramsey', 1.3e-6)
    mw_min_len_ramsey = StatusVar('mw_min_len_ramsey', 0e-6)
    mw_max_len_ramsey = StatusVar('mw_max_len_ramsey', 900e-9)
    mw_steps_ramsey = StatusVar('mw_steps_ramsey', 50e-9)
    mw_len_ramsey = StatusVar('mw_len_ramsey', 200e-9) #whatever the rabi says

    #
    # Internal signals
    sigContinueLoop = QtCore.Signal()
    sigStopMeasurement = QtCore.Signal()
    sigStopAll = QtCore.Signal()
    # sigNextLine = QtCore.Signal()
    #
    ## Update signals, e.g. for GUI module
    sigMeasurementDone = QtCore.Signal()
    # Gives the current count to the GUi
    sigAverageDone = QtCore.Signal(int, np.ndarray, np.ndarray, np.ndarray) # number of averages, counts and average counts
    # sigGetRefValues = QtCore.Signal(int, np.ndarray) # number of averages, ref_average counts
    sigFitUpdated = QtCore.Signal(np.ndarray, np.ndarray, dict, str)
    # sigOutputStateUpdated = QtCore.Signal(str, bool)
    # sigElapsedTimeUpdated = QtCore.Signal(float, int)

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
        self.step_index = 0 # value right now, column


        # self.step_count = self._pulselogic.get_step_count([2.1, 1, 2, 0.5])  #max value for columns
        self.clk_rate_daq = self._photon_counter.get_counter_clock_frequency() # which is 100000

        self.rep = 100
        self.trigger_setting = True

        ## Fit
        self.fit_x = None
        self.fit_y = None
        self.fits_performed = {}

        # Declare where the signals should lead to
        self.sigContinueLoop.connect(self.continue_loop, QtCore.Qt.QueuedConnection)
        self.sigStopMeasurement.connect(self.stop_measurement, QtCore.Qt.QueuedConnection)
        self.sigStopAll.connect(self.stop_all, QtCore.Qt.QueuedConnection)
        self.stopRequested = False

        self._build_pulse_arrays()

        return active_channels

    def _build_pulse_arrays(self):
        # Build arrays from the StatusVar
        self.apd_times = [self.apd_start, self.apd_len]  # [time to start, length] all in microseconds
        self.apd_ref_times = [self.apd_ref_start, self.apd_ref_len]

        self.laser_times = [self.laser_in, self.laser_off, self.laser_re]  # on off on ...the rest is zeros
        # [mw_start_time, mw_len_0, mw_len_max, steps] all in microsecond
        self.mw_times_rabi = [self.mw_start_time_rabi, self.mw_min_len_rabi, self.mw_max_len_rabi, self.mw_steps_rabi]
        # [start, distance_min, distance_max, steps, pulse length] all in microseconds
        self.mw_times_ramsey = [self.mw_start_time_ramsey, self.mw_min_len_ramsey, self.mw_max_len_ramsey,
                                self.mw_steps_ramsey, self.mw_len_ramsey]

        # mw_wait_time between laser off and mw on, mw_pulse_time
        self.mw_times_sweep = [self.mw_start_delay, self.mw_len_delay]
        self.apd_times_sweep = [self.apd_len_delay, self.apd_min_start_delay, self.apd_max_start_delay, self.apd_steps_delay]  # length, min_start, max_start, steps in microseconds


    def on_deactivate(self):
        self._pulser.clear_all()
        self._pulser.reset()
        self._pulser.pulser_off()
        self.sigAverageDone.disconnect()
        self.sigContinueLoop.disconnect()
        self.sigStopMeasurement.disconnect()
        # self.sigNextLine.disconnect()
        self.sigStopAll.disconnect()
        self.mw_off()

    @fc.constructor
    def sv_set_fits(self, val):
        # Setup fit container
        fc = self.fitlogic().make_fit_container('pulsed', '1d')
        fc.set_units(['us', 'counts'])
        if isinstance(val, dict) and len(val) > 0:
            fc.load_from_dict(val)
        else:
            d1 = OrderedDict()
            d1['rabi'] = {
                'fit_function': 'sineexponentialdecay',
                'estimator': 'generic'
            }
            d1['ramsey'] = {
                'fit_function': 'sinedoublewithtwoexpdecay',
                'estimator': 'generic'
            }
            default_fits = OrderedDict()
            default_fits['1d'] = d1
            fc.load_from_dict(default_fits)
        return fc

    @fc.representer
    def sv_get_fits(self, val):
        """ save configured fits """
        if len(val.fit_list) > 0:
            return val.save_to_dict()
        else:
            return None

    def trigger(self):
        self._pulselogic.trigger()

    def cw(self, on=True):
        '''
        plays a zero waveform on all channels but ones on the APD channels. Makes sure the Photon counter can read the counts without changing the cables
        '''
        #FIXME: initialize properly the waveform, remove the test method!
        # change to paly waveform
        #self._pulselogic.do_cw(self.seq_len)
        self._pulser.reset()
        # self._time_series_logic_con.start_reading()
        if on:
            first_out = 1. # The laser
            second_out = 0. # The APD
            third_out = 0.  #The ref APD
            # Counts end up at the very and of the chain in the sepaerate DAQ counter
        else:
            first_out = 0. # The laser
            second_out = 0. # the APD
            third_out = 0. # The ref APD
            # Counts end up at the very and of the chain in the sepaerate DAQ counter and laser is off
            # To do: Here should be the function: play_waveform
        self._pulser.waveform_test(msecondsplay=0.5, first_out=first_out, second_out=second_out,
                                   third_out=third_out, clk_mega=50)


    def stop_awg(self):  # Makes everything stop
        self._pulselogic.stop_awg()

    def awg(self):
        '''
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
        '''
        method = self.method

        if method == 'rabi':
            print('do rabi')
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
        elif method == 'delaysweep':
            if len(self.apd_times_sweep) == 4 and len(self.mw_times_sweep) == 2 and len(self.laser_times) == 3:
                print('do delay sweep')
                apd_times = self.apd_times_sweep
                mw_times = self.mw_times_sweep
            else:
                # print (len(self.apd_times_sweep))
                # print(len(self.mw_times_sweep))
                # print (len(self.mw_times_sweep))
                self.log.error("The input arrays don't have the right size.")
        elif method == 'delaysweep_ref':
            if len(self.apd_times_sweep) == 4 and len(self.mw_times_sweep) == 2 and len(self.laser_times) == 3:
                print('do delay sweep_ref')
                apd_times = self.apd_times_sweep
                mw_times = self.mw_times_sweep
            else:
                # print (len(self.apd_times_sweep))
                # print(len(self.mw_times_sweep))
                # print (len(self.mw_times_sweep))
                self.log.error("The input arrays don't have the right size.")
        else:
            self.log.error("The methode must be one of the following: delay_sweep, delay_sweep_ref, rabi, ramsey.")

        self._pulselogic.play_any(self.clk_rate_awg, self.seq_len, self.laser_times, apd_times, self.apd_ref_times, mw_times, method=method,
                                  rep=self.rep, mw_pulse=self.mw_pulse_setting, trigger=self.trigger_setting, neg_mw=self.neg_mw)
        return method

    def get_x_axis(self):
        steps, max_val, min_val = self.step_counter()
        x_axis = np.linspace(min_val, max_val, steps)
        return x_axis

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
        count_matrix = np.full((int(self.averages), int(self.step_count)), 0)
        count_matrix_ref = np.full((int(self.averages), int(self.step_count)), 0)
        # These are the index for the matrix (step_index: columns and av_index: rows)
        self.count_matrix = count_matrix
        self.count_matrix_ref = count_matrix_ref

        return count_matrix, count_matrix_ref
    #
    def prepare_devices(self):
        """
        Initialize the counter and the settings of the MW device prior to starting scanning.
        """
        # Photon counter
        self._photon_samples = self._pxtime_to_samples()
        self._photon_counter.prepare_counters(samples_to_acquire=self._photon_samples,
                                              mode='pulsing')

        # MW device: Put in the parameters and turn it on
        self._mw_device.set_power(self.mw_power)
        self._mw_device.set_frequency(self.mw_frequency)
        self.mw_on()
        self._mw_device.set_mod(True)

    def start_measurement(self):
        self._build_pulse_arrays()
        self.stop_awg()
        self.set_laser_power() #this one talks to the qwg via awg.cw
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
        #time.sleep(0.1) # in sec
        if self.stopRequested == True:
            self.mw_off()
            self._mw_device.set_mod(False)
            self._photon_counter.close_counters()
            self.stop_awg()
            # make sure after all the inputs are low after the mesurement
            # todo: write a real function for that
            self._pulser.waveform_test(msecondsplay=0.5, loops=0, first_out=0, second_out=0, clk_mega=50)
            # self._pulselogic.post_measurement(self.seq_len)
            self.av_index = 0
            self.step_index = 0
            self.stopRequested = False
            print('measurement stopped')
            # It does not save anything
        else:
            if self.step_index == self.step_count - 1 and self.av_index == self.averages -1:
                # print('very last value')
                counts, ref_counts = self.acquire_pixel()
                # print ('counts: ', counts)
                # print('ref counts: ', ref_counts)
                print(f'Current Average: {self.av_index}')
                self.count_matrix[self.av_index, self.step_index] = counts
                self.count_matrix_ref[self.av_index, self.step_index] = ref_counts
                self.sigStopMeasurement.emit()

            elif self.step_index == self.step_count - 1 and self.av_index < self.averages - 1: # moves to the next average row
                # print('last value in the row')
                print(f'Current Average: {self.av_index}')
                counts, ref_counts = self.acquire_pixel()
                # print('counts: ', counts)
                # print('ref counts: ', ref_counts)
                self.count_matrix[self.av_index, self.step_index] = counts
                self.count_matrix_ref[self.av_index, self.step_index] = ref_counts
                #print(self.count_matrix[self.av_index])
                self.step_index = 0
                idx = self.av_index + 1
                self.av_counts = self.count_matrix[:idx, :].mean(axis=0)
                self.av_counts_ref = self.count_matrix_ref[:idx, :].mean(axis=0)
                self.sigAverageDone.emit(self.av_index, self.count_matrix[self.av_index], self.av_counts, self.av_counts_ref)
                # self.sigGetRefValues.emit(self.av_index, self.av_counts_ref)
                # print('averageds array', self.av_counts)
                self.av_index = self.av_index + 1
                self.sigContinueLoop.emit()


            else:       # just continue in this --> direction
                # print('switch to the next column, sequence')
                counts, ref_counts = self.acquire_pixel()
                # print('counts: ', counts)
                # print('ref counts: ', ref_counts)
                self.count_matrix[self.av_index, self.step_index] = counts
                self.count_matrix_ref[self.av_index, self.step_index] = ref_counts
                self.step_index = self.step_index + 1

                self.sigContinueLoop.emit()
        return self.count_matrix

    def stop_measurement(self):
        self.mw_off()
        self._photon_counter.close_counters()
        self.stop_awg()
        #make sure after all the inputs are low after the mesurement
        # todo: write a real function for that
        self._pulser.waveform_test(msecondsplay=0.5, loops=0, first_out=0, second_out=0, clk_mega=50)
        # self._pulselogic.post_measurement(self.seq_len)
        self._mw_device.set_mod(False)
        self.av_index = 0
        self.step_index = 0
        self.stopRequested = False
        print('measurement done')
        self.sigMeasurementDone.emit()
        #SAVE?
        #self.save_data(self.method)

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
        data_raw['averaged signal'] = self.av_counts

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

        attributes = {'count data': parameters, 'REF count data': parameters}

        # Check if a fit has been performed
        if method in self.fits_performed:
            x_fit, y_fit, result = self.fits_performed[method]
            data_raw['fit values'] = np.stack((x_fit, y_fit))
            attributes['fit values'] = result.values

        self._savelogic.save_hdf5_data(data_raw,
                                       attributes=attributes,
                                       filename=timestamp.strftime('%Y-%m-%d-%H-%M-%S') + '_' + self.method + '.h5')

    def acquire_pixel(self):
        # Requests the counts from the DAQ card, gives one number only
        samples = int(self.integration_time * self.clk_rate_daq)
        data, _ = self._photon_counter.read_pixel(samples)
        # ref_counts, _ = self._photon_counter.read_pixel(x)
        counts, ref_counts = data
        return counts, ref_counts

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

    # def cw_ttls_easy(self): # remove this
    #     #todo: write a real function for that...
    #     self._pulser.waveform_test(msecondsplay=0.5, loops=0, first_out=1, second_out=0, clk_mega=50)

    def address_ttls(self, seq_len=5.0, first_out=1, second_out=0, third_out=0):
        '''
        Fix this one and put it in continue loop and stop measurement
        Play ones to the APD channel to make sure the counts get through and can be used in a cw measurement.
        Here I send zeros to all the other channels: Laser, mw, and APD ref
        fist_out = channelX0 connected to laser, this should be high for cw mode
        second_out = channelX1 is connected to the first TTL switch --> should be low for cw mode
        third_out = channelX2 is connected to 2nd TTL switch --> should be low for the cw mode
        in the end the counts get all the way through both ttls and end up it the User BNC connector
        '''
        seq_len *= 1e-6
        card_idx = 1
        loops = 0 # This lets the waveform play infinitely
        self._pulser.set_sample_rate(card_idx = card_idx, clk_rate=int(50000000)) # low sampling rate is fine for that
        clk_rate = self._pulser.get_sample_rate(card_idx)

        #self.set_sample_rate(card_idx, clk_rate)

        # # This sequence immediately starts after the sequences are loaded
        self._pulser.configure_ORmask(card_idx, 'immediate')
        self._pulser.configure_ANDmask(card_idx, None)

        samples = self._pulser.waveform_padding(int(seq_len * clk_rate))
        time_ax = np.linspace(0, samples / clk_rate, samples)

        do_chan = 1
        do0 = first_out * np.ones(time_ax.shape, dtype=np.int64) # laser
        do1 = second_out * np.ones(time_ax.shape, dtype=np.int64) # apd
        do2 = third_out * np.ones(time_ax.shape, dtype=np.int64) # apd ref
        # print('do0', do0)
        # print('do1', do1)
        # print('do2', do2)
        outchan = 0
        do_output = {1: do0, 2: do1, 3: do2}
        # digital_output_map = {1: do0, 2: do1, 3: do2}
        digital_output_map = {1: [0, 1, 2]}
        # digital_output_map = {1: [0], 2: [1], 3: [2]}
        self._pulser.load_waveform(do_waveform_dictionary=do_output,
                           digital_output_map=digital_output_map)
        self._pulser.play_waveform(self._pulser._waveform_container[-1], loops=loops) # What does this [-1] do? This referes to another class in the hardware code

    def step_counter(self):
        self._build_pulse_arrays()

        if self.method == 'rabi':
            step_count = self._pulselogic.get_step_count(self.mw_times_rabi)
            max_val = self.mw_times_rabi[2]
            min_val = self.mw_times_rabi[1]
        elif self.method == 'ramsey':
            step_count = self._pulselogic.get_step_count(self.mw_times_ramsey)
            max_val = self.mw_times_ramsey[2]
            min_val = self.mw_times_ramsey[1]
        else:
            step_count = self._pulselogic.get_step_count(self.apd_times_sweep)
            max_val = self.apd_times_sweep[2]
            min_val = self.apd_times_sweep[1]

        return step_count, max_val, min_val

    def set_laser_power(self):
        # self.cw()
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
                self.fc.set_current_fit('No Fit')
                if fit_function != 'No Fit':
                    self.log.warning('Fit function "{0}" not available in ODMRLogic fit container.'
                                     ''.format(fit_function))

        self.fit_x, self.fit_y, result = self.fc.do_fit(x_data, y_data)

        if fit_function != 'No Fit':
            self.fits_performed[self.fc.current_fit] = (self.fit_x, self.fit_y, result)
        # else:
        #     if key in self.fits_performed:
        #         self.fits_performed.pop(key)

        if result is None:
            result_str_dict = {}
        else:
            result_str_dict = result.result_str_dict

        print('emitting update fit')
        self.sigFitUpdated.emit(
            self.fit_x, self.fit_y, result_str_dict, self.fc.current_fit
        )
        return

    def get_fit_functions(self):
        """ Return the hardware constraints/limits
        @return list(str): list of fit function names
        """
        return list(self.fc.fit_list)

        # self._pulser.start_card(card_idx)
        # self._pulser.arm_trigger(card_idx)

    def start_stop_timetrace(self, val):
        if val == True:
            #Start timetrace
            self._timeseries.start_reading()
        else:
            self._timeseries.stop_reading()
