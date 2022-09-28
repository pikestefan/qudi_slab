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
from hardware.awg.spectrum_awg.py_header.regs import *
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

    integration_time = StatusVar('integration_time', 30e-3)  # In ms for timing...

    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)

    def on_activate(self):
    # Initialisation performed during activation of the module

    # Get connectors
        self._pulser = self.pulsegenerator()
        self._photon_counter = self.photon_counter()
        self._savelogic = self.savelogic()
        self._pulselogic = self.pulselogic()
        active_channels = self._pulser.get_active_channels()
        self._pulser.clear_all()
        self._pulser.reset()

        # Initalize the values which get stored later
        self.count_matrix = None  # This matrix is used to store the DAQ card values.

        return active_channels

    def on_deactivate(self):
        self._pulser.clear_all()
        self._pulser.reset()
        self._pulser.pulser_off()

    def trigger(self):
        self._pulselogic.trigger()

    def stop_awg(self):  # Makes everything stop
        self._pulselogic.stop_awg()


    def awg(self, clk_rate, seq_len, laser_times, apd_times=[], mw_times=[], method=str, rep=100, mw_pulse=bool,
                 trigger=bool):
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
        method =    can be rabi, ramsey or delaysweep
        rep =       repetition time when no trigger is used
        mw_pulse =  With or without MW pulse in the middle
        trigger =   either True: software trigger or False: no trigger at all
        '''

        self._pulselogic.play_any(clk_rate, seq_len, laser_times, apd_times=apd_times, mw_times=mw_times, method=method, rep=100, mw_pulse=mw_pulse, trigger=trigger)


    def _prepare_count_matrix(self, step_count, av_number):
        '''
        creates a matrix with av_number of rows and step_count of columns

        '''
        count_matrix = np.full((av_number, step_count), self.invalid)
        self._freq_scanning_index = 0
        self._average_index = 0
        self.count_matrix = count_matrix

    def _prepare_devices(self):
        """
        Initialize the counter and the settings of the MW device prior to starting scanning.
        """

        self._photon_counter.module_state.lock()
        self._photon_samples = self._pxtime_to_samples()
        self._photon_counter.prepare_counters(samples_to_acquire=self._photon_samples)


    # def save_data(self, method, laser_times, apd_times, mw_times):
    #      pulse = {
    #             'method': self.method,
    #             'laser_times': self.laser_times,
    #             'apd_times': self.apd_times,
    #             'mw_times': self.mw_times,
    #             'AWG_settings': self.clk_rate
    #      }
    #
    #     data = {**pulse}
    #     self._savelogic.save_hdf5_data(data)


    ### Get the actual data from card
    # this is from continue_odmr
    # sometimes its just empty
    # if len(counts) == 0:
    #     return
    #
    # self.count_matrix[self._average_index, self._freq_scanning_index] = counts / self.integration_time
    # self.curr_odmr_trace = self.count_matrix[self._average_index]
    #
    # if self._average_index > 0:
    #     self.average_odmr_trace = np.nanmean(self.count_matrix, axis=0)


    def _acquire_pixel(self):
        counts, _ = self._photon_counter.read_pixel(self._photon_samples)
        return counts