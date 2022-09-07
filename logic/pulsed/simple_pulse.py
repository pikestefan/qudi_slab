# -*- coding: utf-8 -*-

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
along with Qudi. If not, see <http://www.gnu.org/licenses/>.

Copyright (c) the Qudi Developers. See the COPYRIGHT.txt file at the
top-level directory of this distribution and at <https://github.com/Ulm-IQO/qudi/>
"""

from qtpy import QtCore
from collections import OrderedDict
import numpy as np
import copy
import time
import math
import datetime
import matplotlib.pyplot as plt
from scipy import signal
from hardware.awg.spectrum_awg.py_header.regs import *

from core.connector import Connector
from core.configoption import ConfigOption
from core.statusvariable import StatusVar
from core.util.mutex import Mutex
from core.util.network import netobtain
from core.util import units
from core.util.math import compute_ft
from logic.generic_logic import GenericLogic


# Just hardcode the channels to the different hardware modules
    # Microwave: analog channel 2 and 3
    # Laser: digital channel X0
    # Photon counter: digital channel X1 and X2 (for reference)
class Simplepulse(GenericLogic):

    pulsegenerator = Connector(interface="PulserInterface")

    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)

    def on_activate(self):
        #     """ Initialisation performed during activation of the module.
        #     """

        # Get connectors
        self._pulser = self.pulsegenerator()
        active_channels = self._pulser.get_active_channels()
        return active_channels
        # Turn off pulse generator somehow
        # awg.pulser_off()
        # pulsegenerator.on_activate()
        # # pulsegenerator.set_chan_amplitude(channels, amplitudes)
        # # pulsegenerator.set_output_filters(channels, filter_active)
        # pulsegenerator.set_sample_rate(1, clk_rate=MEGA(50))

    #     return

    def on_deactivate(self):
        self._pulser.clear_all()
        self._pulser.pulser_off()

<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
    def set_up_pulser(self,card_idx,channels, clk_rate, filter_active):
    # def set_up_pulser(self):

        clk_rate = MEGA(clk_rate)
=======
    def set_up_pulser(self, card_idx, channels, clk_rate_raw, filter_active):
        # def set_up_pulser(self):
        clk_rate = 10 * 1000 * clk_rate_raw
>>>>>>> Stashed changes
=======
    def set_up_pulser(self, card_idx, channels, clk_rate_raw, filter_active):
        # def set_up_pulser(self):
        clk_rate = 10 * 1000 * clk_rate_raw
>>>>>>> Stashed changes
=======
    def set_up_pulser(self, card_idx, channels, clk_rate_raw, filter_active):
        # def set_up_pulser(self):
        clk_rate = 10 * 1000 * clk_rate_raw
>>>>>>> Stashed changes
=======
    def set_up_pulser(self, card_idx, channels, clk_rate_raw, filter_active):
        # def set_up_pulser(self):
        clk_rate = 10 * 1000 * clk_rate_raw
>>>>>>> Stashed changes
=======
    def set_up_pulser(self, card_idx, channels, clk_rate_raw, filter_active):
        # def set_up_pulser(self):
        clk_rate = 10 * 1000 * clk_rate_raw
>>>>>>> Stashed changes
        # What needs to be done:
        # pulser_off
        # 1. clear all
        # set_chan_amplitude
        # set_output_filters
        # set_sample_rate
        self._pulser.pulser_off()
        self._pulser.clear_all()
        self._pulser.set_output_filters(channels, filter_active)
        self._pulser.set_sample_rate(card_idx, clk_rate)
        output_sample_rate = self._pulser.get_sample_rate(card_idx)
        return output_sample_rate

<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
### Create Blocks ###
    def create_d_on(self,msplay, clk_rate):
        '''
        Creates an array of ones for a digital channel
        len_ms:         is the desired duration in ms
        clk_rate:       samples/second of the AWG
        '''
        clk_rate = MEGA(clk_rate)
        samples = self._pulser.waveform_padding(msplay * clk_rate)
        d_on=np.ones(samples)
=======
    ### For digital signals ###
    def create_d_on(self, samples, clk_rate):
        """
        Creates an array of ones for a digital channel
        len_ms:         is the desired duration in ms
        clk_rate:       samples/second of the AWG
=======
    ### For digital signals ###
    def create_d_on(self, samples, clk_rate):
        """
        Creates an array of ones for a digital channel
        len_ms:         is the desired duration in ms
        clk_rate:       samples/second of the AWG
>>>>>>> Stashed changes
=======
    ### For digital signals ###
    def create_d_on(self, samples, clk_rate):
        """
        Creates an array of ones for a digital channel
        len_ms:         is the desired duration in ms
        clk_rate:       samples/second of the AWG
>>>>>>> Stashed changes
=======
    ### For digital signals ###
    def create_d_on(self, samples, clk_rate):
        """
        Creates an array of ones for a digital channel
        len_ms:         is the desired duration in ms
        clk_rate:       samples/second of the AWG
>>>>>>> Stashed changes
=======
    ### For digital signals ###
    def create_d_on(self, samples, clk_rate):
        """
        Creates an array of ones for a digital channel
        len_ms:         is the desired duration in ms
        clk_rate:       samples/second of the AWG
>>>>>>> Stashed changes
        """
        clk_rate = 10 * clk_rate
        # samples = self._pulser.waveform_padding(len_ms * clk_rate)
        d_on = np.ones(samples)
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
        # print('d_on=', d_on)
        # print('len of d_on=', len(d_on))
        return d_on

<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
    def create_d_off(self,msplay, clk_rate):
        '''
        Creates an array of zeros for a digital channel
        len_ms:         is the desired duration in ms
        clk_rate:       samples/second of the AWG
        '''
        clk_rate = MEGA(clk_rate)
        samples = self._pulser.waveform_padding(msplay * clk_rate)
=======
    def create_d_off(self, samples, clk_rate):
        """
        Creates an array of zeros for a digital channel
        len_ms:         is the desired duration in ms
        clk_rate:       samples/second of the AWG
=======
    def create_d_off(self, samples, clk_rate):
        """
        Creates an array of zeros for a digital channel
        len_ms:         is the desired duration in ms
        clk_rate:       samples/second of the AWG
>>>>>>> Stashed changes
=======
    def create_d_off(self, samples, clk_rate):
        """
        Creates an array of zeros for a digital channel
        len_ms:         is the desired duration in ms
        clk_rate:       samples/second of the AWG
>>>>>>> Stashed changes
=======
    def create_d_off(self, samples, clk_rate):
        """
        Creates an array of zeros for a digital channel
        len_ms:         is the desired duration in ms
        clk_rate:       samples/second of the AWG
>>>>>>> Stashed changes
=======
    def create_d_off(self, samples, clk_rate):
        """
        Creates an array of zeros for a digital channel
        len_ms:         is the desired duration in ms
        clk_rate:       samples/second of the AWG
>>>>>>> Stashed changes
        """
        clk_rate = 10 * clk_rate
        # samples = self._pulser.waveform_padding(len_ms * clk_rate)
>>>>>>> Stashed changes
        d_off = np.zeros(samples)
        # print('d_off=', d_off)
        # print('len of d_off=', len(d_off))
        return d_off

    ### For analog signals ###
    def sine_wave(self, msplay, clk_rate):
        msplay *= 1e-3
        samples = self._pulser.waveform_padding(msplay * clk_rate)
        time_ax = np.linspace(0, samples / clk_rate, samples)
        waveform = np.sin(2 * np.pi * time_ax / msplay)
        return waveform
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream

    def create_ao_dict(self,msplay,clk_rate,channels=[]):
        channels = [2, 3]
        # Building a dict for anaolog signal which can be read by the AWG
        d_off = self.create_d_off(msplay,clk_rate)
        d_on = self.create_d_on(msplay, clk_rate)
        aosequence = [{channels[0]: d_on, channels[1]: d_off}]
        return aosequence

    def create_do_dict(self,msplay,clk_rate,channels=[], outchan = 0):
        # Building a dict for analog signal which can be read by the AWG
        d_off = self.create_d_off(msplay,clk_rate)
        d_on = self.create_d_on(msplay, clk_rate)
        # Building a dict for anaolog signal which can be read by the AWG
        aosequence = [{channels[0]: d_on}]
        dosequence = [{outchan: d_on}]
        return dosequence

## Plays the dict
    def load_n_play(self, msplay, clk_rate, card_idx, segment_map=[], trigger=bool,rep=100):
        msplay *= 1e-3
        aosequence= self.create_ao_dict(msplay, clk_rate)
        dosequence= self.create_do_dict(msplay, clk_rate)
        # print('aos:', aosequence)
        # print('dos:', dosequence)
        digital_output_map = {1: [0]}  # directs to the right output channel
        for aos, dos in zip(aosequence, dosequence):
            # Loads waveforms to the awg
            self._pulser.load_waveform(ao_waveform_dictionary=aos, do_waveform_dictionary=dos,
                                       digital_output_map=digital_output_map)

        if trigger == True:
            # This sequence waits for a software trigger to start playing and moving to the next step.
            self._pulser.configure_ORmask(card_idx, None)
            self._pulser.configure_ANDmask(card_idx, None)
            stop_condition_list = np.array([0x40000000], dtype=np.int64)
            loops = np.ones((len(aosequence)), dtype=np.int64)

        else:
          # This sequence immediately starts after the sequences are loaded
            self._pulser.configure_ORmask(card_idx, 'immediate')
            self._pulser.configure_ANDmask(card_idx, None)
            loops = np.array([rep], dtype=np.int64)
        stop_condition_list = np.array([])

        # makes the waveforms repeat in loops
        self._pulser.load_sequence(  # digital_sequences=dosequence,
            loops_list=loops, segment_map=segment_map, stop_condition_list=stop_condition_list)
        self._pulser.start_card(card_idx)
        self._pulser.arm_trigger(card_idx)
        self._pulser.clear_all()
        # self._pulser.reset()

    def software_trigger(self):
        self._pulser.send_software_trig(1)

    # build a sequence where you can deterine the length and distances, Its only digital
    def laser_trace(self, clk_rate, laser_in, wait, laser_re, trigger=bool, rep=100, segment_map=[]):
        laser_in *= 1e-3
        laser_re *= 1e-3
        wait *= 1e-3
        msplay = laser_in + wait + laser_re
        clk_rate = MEGA(clk_rate)
        samples = self._pulser.waveform_padding(msplay * clk_rate)
        print('len of samples:', samples)
        time_ax = np.linspace(0, samples / clk_rate, samples)
        first_pulse = np.ones(int(clk_rate * laser_in))
        wait_time = np.zeros(int(clk_rate * wait))
        sec_pulse = np.ones(int(clk_rate * laser_re))
        print('len of first pulse in samples:', len(first_pulse))
        print('len of first pulse in s:', len(first_pulse)/ clk_rate)
        print('len of wait time in samples:', len(wait_time))
        print('len of wait time in s:', len(wait_time) / clk_rate)
        print('len of sec_pusle time in samples:', len(sec_pulse))
        print('len of sec_pulse in s:', len(sec_pulse) / clk_rate)
        do_waveform = np.concatenate((first_pulse, wait_time, sec_pulse))
        print('len of whole thing:', len(do_waveform))
        print('whole array(do_waveform):', do_waveform)
        outchan = 0
        # Building a dict for analog signal which can be read by the AWG
        d_off = np.zeros(int(clk_rate * msplay))
        # Building a dict for anaolog signal which can be read by the AWG
        aosequence = [{2: d_off}]
        dosequence = [{outchan: do_waveform}]
        print('aosequence:', aosequence)
        print('dosequence:', dosequence)
        card_idx = 1
        # load and play
        digital_output_map = {1: [0]}  # directs to the right output channel
        for aos, dos in zip(aosequence, dosequence):
            # Loads waveforms to the awg
            self._pulser.load_waveform(ao_waveform_dictionary=aos, do_waveform_dictionary=dos,
                                       digital_output_map=digital_output_map)

        if trigger == True:
            # This sequence waits for a software trigger to start playing and moving to the next step.
            self._pulser.configure_ORmask(card_idx, None)
            self._pulser.configure_ANDmask(card_idx, None)
            stop_condition_list = np.array([0x40000000], dtype=np.int64)
            loops = np.ones((len(aosequence)), dtype=np.int64)

        else:
            # This sequence immediately starts after the sequences are loaded
            self._pulser.configure_ORmask(card_idx, 'immediate')
            self._pulser.configure_ANDmask(card_idx, None)
            loops = np.array([rep], dtype=np.int64)
            stop_condition_list = np.array([])

        # makes the waveforms repeat in loops
        self._pulser.load_sequence(  # digital_sequences=dosequence,
            loops_list=loops, segment_map=segment_map, stop_condition_list=stop_condition_list)
        self._pulser.start_card(card_idx)
        self._pulser.arm_trigger(card_idx)
        self._pulser.clear_all()
        # self._pulser.reset()

    def delay_block(self, clk_rate, laser_in, wait, laser_re, apd_start, apd_stop, rest, trigger=bool, rep=100, segment_map=[]):
        '''
        This can play a laser sequence and the apd readout at the same time
        clk_rate:   something like 100
        laser_in:   length of initialisation laser pulse (it starts right at the beginning) in ms
        wait:       time between the end of the initialisation and the reinitialisation in ms
        laser_re:   length of the reinitialisation pulse in ms
        apd_start:  time in ms from the start until the apd should turn on (this needs to be varied to do the delay measurement)
        apd_stop:   time in ms where the apd is on (this needs to be varied to do the delay measurement)
        rest:       time in ms after the reinitialisation pulse to have a break before the next cycle starts
        trigger:    True if the software trigger should be used and False if the loop array should define the length
        rep:        If there is no software trigger, rep gives the number of repetitions the waveform is played
        segment_map:If the order of different segments should be mixed (not there yet)
        '''
        laser_in *= 1e-3
        laser_re *= 1e-3
        wait *= 1e-3
        apd_start *= 1e-3
        apd_stop *= 1e-3
        rest *= 1e-3
        msplay = laser_in + wait + laser_re + rest
        clk_rate = MEGA(clk_rate)
        samples = self._pulser.waveform_padding(msplay * clk_rate)
        #for laser
        first_pulse = np.ones(self._pulser.waveform_padding(clk_rate * laser_in))
        laser_wait_time = np.zeros(self._pulser.waveform_padding(clk_rate * wait))
        sec_pulse = np.ones(self._pulser.waveform_padding(clk_rate * laser_re))
        rest_array = np.zeros(self._pulser.waveform_padding(clk_rate * rest))
        #for apd
        apd_wait_time = np.zeros(self._pulser.waveform_padding(clk_rate * apd_start))
        apd_on = np.ones(self._pulser.waveform_padding(clk_rate * apd_stop))
        apd_off = np.zeros(self._pulser.waveform_padding(clk_rate * (msplay - (apd_start + apd_stop))))

        print('len of first pulse in samples:', len(first_pulse))
        print('len of first pulse in s:', len(first_pulse)/ clk_rate)
        print('len of wait time in samples:', len(laser_wait_time))
        print('len of wait time in s:', len(laser_wait_time) / clk_rate)
        print('len of sec_pusle time in samples:', len(sec_pulse))
        print('len of sec_pulse in s:', len(sec_pulse) / clk_rate)
        laser_do_waveform = np.concatenate((first_pulse, laser_wait_time, sec_pulse, rest_array))
        apd_do_waveform = np.concatenate((apd_wait_time, apd_on, apd_off))
        print('len of laser thing:', len(laser_do_waveform))
        print('len of apd thing:', len(apd_do_waveform))
        print('whole array laser thing:', laser_do_waveform)
        print('whole array apd thing:', apd_do_waveform)

        if len(laser_do_waveform) == len(apd_do_waveform):
            print('both waveforms have the same length')
        elif len(laser_do_waveform) < len(apd_do_waveform):
            print('laser waveform is shorter than apd waveform')
            difference = len(apd_do_waveform) - len(laser_do_waveform)
            print(difference)
            rest_array = np.zeros(int(difference))
            print(rest_array)
            laser_do_waveform = np.concatenate((laser_do_waveform, rest_array))
            print('new len of laser thing:', len(laser_do_waveform))
            print('len of apd thing:', len(apd_do_waveform))
        else:
            print('apd waveform is shorter than laser waveform')
            difference = len(laser_do_waveform) - len(apd_do_waveform)
            rest_array = np.zeros(int(difference))
            apd_do_waveform = np.concatenate((apd_do_waveform, rest_array))
            print('new len of laser thing:', len(laser_do_waveform))
            print('new len of apd thing:', len(apd_do_waveform))



        laser_outchan = 0 # Laser
        apd_outchan = 1
        # Building a dict for analog signal which can be read by the AWG
        d_off = np.zeros(len(apd_do_waveform))
        # Building a dict for anaolog signal which can be read by the AWG
        aosequence = [{2: d_off}]
        dosequence = [{laser_outchan: laser_do_waveform, apd_outchan: apd_do_waveform}]
        print('aosequence:', aosequence)
        print('dosequence:', dosequence)
        card_idx = 1

        # load and play
        digital_output_map = {0: [0,1]} # directs to the right output channel analog channel0 directs to digital x0 and x1
        for aos, dos in zip(aosequence, dosequence):
            # Loads waveforms to the awg
            self._pulser.load_waveform(ao_waveform_dictionary=aos, do_waveform_dictionary=dos,
                                       digital_output_map=digital_output_map)

        if trigger == True:
            # This sequence waits for a software trigger to start playing and moving to the next step.
            self._pulser.configure_ORmask(card_idx, None)
            self._pulser.configure_ANDmask(card_idx, None)
            stop_condition_list = np.array([0x40000000], dtype=np.int64)
            loops = np.ones((len(aosequence)), dtype=np.int64)

        else:
            # This sequence immediately starts after the sequences are loaded
            self._pulser.configure_ORmask(card_idx, 'immediate')
            self._pulser.configure_ANDmask(card_idx, None)
            loops = np.array([rep], dtype=np.int64)
            stop_condition_list = np.array([])

        # makes the waveforms repeat in loops
        self._pulser.load_sequence(  # digital_sequences=dosequence,
            loops_list=loops, segment_map=segment_map, stop_condition_list=stop_condition_list)
        self._pulser.start_card(card_idx)
        self._pulser.arm_trigger(card_idx)
        self._pulser.clear_all()
        # self._pulser.reset()


    def delay_sweep(self, clk_rate, laser_in, wait, laser_re, apd_stop, rest,start, steps, trigger=bool, rep=100, segment_map=[]):
        '''
        This can play a laser sequence and the apd readout at the same time
        clk_rate:   something like 100
        laser_in:   length of initialisation laser pulse (it starts right at the beginning) in ms
        wait:       time between the end of the initialisation and the reinitialisation in ms
        laser_re:   length of the reinitialisation pulse in ms

        apd_stop:   time in ms where the apd is on (this needs to be varied to do the delay measurement)
        rest:       time in ms after the reinitialisation pulse to have a break before the next cycle starts
        trigger:    True if the software trigger should be used and False if the loop array should define the length
        rep:        If there is no software trigger, rep gives the number of repetitions the waveform is played
        segment_map:If the order of different segments should be mixed (not there yet)
        '''
        #building the array to sweep apd_start
        msplay = laser_in + wait + laser_re + rest
        clk_rate = MEGA(clk_rate)
        stop = (msplay - rest) - apd_stop
        # print('stop value: ', stop)
        #build array with all the desired apd_start values]
        apd_start = []
        step_count = math.floor((stop-start)/steps)
        # print('step_count: ', step_count)
        for i in range(step_count):
            apd_start.append((start+(i*steps)))
        # print('adpstart array: ', apd_start)
        laser_in *= 1e-3
        laser_re *= 1e-3
        wait *= 1e-3
        apd_stop *= 1e-3
        rest *= 1e-3
        start *= 1e-3
        steps *= 1e-3
        msplay *= 1e-3
        apd_start = np.multiply(1e-3, apd_start)
        # print('apd_start: ', apd_start)

        #Laser waveform does not change the whole time
        first_pulse = np.ones(self._pulser.waveform_padding(clk_rate * laser_in))
        laser_wait_time = np.zeros(self._pulser.waveform_padding(clk_rate * wait))
        sec_pulse = np.ones(self._pulser.waveform_padding(clk_rate * laser_re))
        rest_array = np.zeros(self._pulser.waveform_padding(clk_rate * rest))
        laser_do_waveform = np.concatenate((first_pulse, laser_wait_time, sec_pulse, rest_array))
        dosequence = []
        aosequence = []
        laser_outchan = 0  # Laser
        apd_outchan = 1  # APD
        card_idx = 1
        for i in apd_start:
            apd_wait_time = np.zeros(self._pulser.waveform_padding(clk_rate * i))
            apd_on = np.ones(self._pulser.waveform_padding(clk_rate * apd_stop))
            apd_off = np.zeros(self._pulser.waveform_padding(clk_rate * (msplay - (i + apd_stop))))
            apd_do_waveform = np.concatenate((apd_wait_time, apd_on, apd_off))
            if len(laser_do_waveform) == len(apd_do_waveform):
                # print('len of apd_do_waveform: ', len(apd_do_waveform))
                print('both waveforms have the same length')
            elif len(laser_do_waveform) < len(apd_do_waveform):
                # print('laser waveform is shorter than apd waveform')
                difference = len(apd_do_waveform) - len(laser_do_waveform)
                # print(difference)
                rest_array = np.zeros(int(difference))
                # print(rest_array)
                laser_do_waveform = np.concatenate((laser_do_waveform, rest_array))
                # print('new len of laser thing:', len(laser_do_waveform))
                # print('len of apd thing:', len(apd_do_waveform))
            else:
                # print('apd waveform is shorter than laser waveform')
                difference = len(laser_do_waveform) - len(apd_do_waveform)
                rest_array = np.zeros(int(difference))
                apd_do_waveform = np.concatenate((apd_do_waveform, rest_array))
                # print('new len of laser thing:', len(laser_do_waveform))
                # print('new len of apd thing:', len(apd_do_waveform))
            d_off = np.zeros(len(apd_do_waveform))
            do_dict = {laser_outchan: laser_do_waveform, apd_outchan: apd_do_waveform}
            ao_dict = {2: d_off}
            dosequence.append((do_dict))
            aosequence.append((ao_dict))

        # print('dosequence: ', dosequence)
        # print('aosequence: ', aosequence)
        # load and play
        # print('len of dosequence: ', len(dosequence))
        # print('len of aosequence: ', len(aosequence))
        digital_output_map = {0: [0, 1]}  # directs to the right output channel analog channel0 directs to digital x0 and x1
        for aos, dos in zip(aosequence, dosequence):
            # Loads waveforms to the awg
            self._pulser.load_waveform(ao_waveform_dictionary=aos, do_waveform_dictionary=dos,
                                           digital_output_map=digital_output_map)

        if trigger == True:
            # This sequence waits for a software trigger to start playing and moving to the next step.
            self._pulser.configure_ORmask(card_idx, None)
            self._pulser.configure_ANDmask(card_idx, None)
            array=np.array([0x40000000], dtype=np.int64)
            stop_condition_list = np.repeat(array, step_count)
            # print('stop_condition_list: ', stop_condition_list)
            loops = np.ones((len(aosequence)), dtype=np.int64)

        else:
            # This sequence immediately starts after the sequences are loaded
            self._pulser.configure_ORmask(card_idx, 'immediate')
            self._pulser.configure_ANDmask(card_idx, None)
            loops = np.array([rep], dtype=np.int64)
            stop_condition_list = np.array([])

        # makes the waveforms repeat in loops
        self._pulser.load_sequence(  # digital_sequences=dosequence,
            loops_list=loops, segment_map=segment_map, stop_condition_list=stop_condition_list)
        self._pulser.start_card(card_idx)
        self._pulser.arm_trigger(card_idx)
        # self._pulser.clear_all()
        # self._pulser.reset()


    def delay_sweep_mw(self, clk_rate, laser_in, wait, laser_re, apd_stop, rest, start, steps, mw_wait_time, mw_pulse_time, rep, mw_pulse=bool, trigger=bool, segment_map=[]):

        '''
        This can play a laser sequence and the apd readout at the same time
        clk_rate:   something like 100
        laser_in:   length of initialisation laser pulse (it starts right at the beginning) in ms
        wait:       time between the end of the initialisation and the reinitialisation in ms
        laser_re:   length of the reinitialisation pulse in ms

        apd_stop:   time in ms where the apd is on (this needs to be varied to do the delay measurement)
        rest:       time in ms after the reinitialisation pulse to have a break before the next cycle starts
        start:      where in time should the apd pulse start shifting
        steps:      step size of the apd pulse sweep
        mw_wait_time: when should the mw pulse should happen?
        mw_pulse_time: length of the mw pulse in between the two laser pulses
        trigger:    True if the software trigger should be used and False if the loop array should define the length
        rep:        If there is no software trigger, rep gives the number of repetitions the waveform is played
        segment_map:If the order of different segments should be mixed (not there yet)
        '''
        #building the array to change the length of the apd pulse
        msplay = laser_in + wait + laser_re + rest
        clk_rate = MEGA(clk_rate)
        stop = (msplay - rest) - apd_stop
        # print('stop value: ', stop)

        #build array with all the desired apd_start values]
        apd_start = []
        step_count = math.floor((stop-start)/steps)
        # print('step_count: ', step_count)
        for i in range(step_count):
            apd_start.append((start+(i*steps)))
        # print('adpstart array: ', apd_start)
        laser_in *= 1e-3
        mw_wait_time *= 1e-3
        mw_pulse_time *= 1e-3
        laser_re *= 1e-3
        wait *= 1e-3
        apd_stop *= 1e-3
        rest *= 1e-3
        start *= 1e-3
        steps *= 1e-3
        msplay *= 1e-3
        apd_start = np.multiply(1e-3, apd_start)
        # print('apd_start: ', apd_start)

        #Laser waveform does not change the whole time
        first_pulse = np.ones(self._pulser.waveform_padding(clk_rate * laser_in))
        laser_wait_time = np.zeros(self._pulser.waveform_padding(clk_rate * wait))
        sec_pulse = np.ones(self._pulser.waveform_padding(clk_rate * laser_re))
        rest_array = np.zeros(self._pulser.waveform_padding(clk_rate * rest))
        laser_do_waveform = np.concatenate((first_pulse, laser_wait_time, sec_pulse, rest_array))
        dosequence = []
        aosequence = []
        laser_outchan = 0  # Laser
        apd_outchan = 1  # APD
        card_idx = 1
        for i in apd_start:
            apd_wait_time = np.zeros(self._pulser.waveform_padding(clk_rate * i))
            apd_on = np.ones(self._pulser.waveform_padding(clk_rate * apd_stop))
            apd_off = np.zeros(self._pulser.waveform_padding(clk_rate * (msplay - (i + apd_stop))))
            apd_do_waveform = np.concatenate((apd_wait_time, apd_on, apd_off))
            if len(laser_do_waveform) == len(apd_do_waveform):
                print('both waveforms have the same length')
            elif len(laser_do_waveform) < len(apd_do_waveform):
                print('laser waveform is shorter than apd waveform')
                difference = len(apd_do_waveform) - len(laser_do_waveform)
                # print(difference)
                rest_array = np.zeros(int(difference))
                # print(rest_array)
                laser_do_waveform = np.concatenate((laser_do_waveform, rest_array))
                # print('new len of laser thing:', len(laser_do_waveform))
                # print('len of apd thing:', len(apd_do_waveform))
            else:
                print('apd waveform is shorter than laser waveform')
                difference = len(laser_do_waveform) - len(apd_do_waveform)
                rest_array = np.zeros(int(difference))
                apd_do_waveform = np.concatenate((apd_do_waveform, rest_array))
                # print('new len of laser thing:', len(laser_do_waveform))
                # print('new len of apd thing:', len(apd_do_waveform))
            if mw_pulse == True:
                mw_wait = np.zeros(self._pulser.waveform_padding(clk_rate * mw_wait_time))
                pi_pulse_array = np.ones(self._pulser.waveform_padding(clk_rate * mw_pulse_time))
                mw_off = np.zeros(self._pulser.waveform_padding(clk_rate * msplay-(mw_wait_time+mw_pulse_time)))
                mw_ao_waveform = np.concatenate((mw_wait, pi_pulse_array, mw_off))
                # print('mw_ao_waveform: ', mw_ao_waveform)
                # print('len of mw_ao_waveform: ', len(mw_ao_waveform))

                # Make sure the arrays have the same length

                if len(laser_do_waveform) == len(mw_ao_waveform):
                    print('both waveforms have the same length')
                elif len(mw_ao_waveform) < len(laser_do_waveform):
                    print('mw waveform is shorter than laser waveform')
                    difference = len(laser_do_waveform) - len(mw_ao_waveform)
                    rest_array = np.zeros(int(difference))
                    mw_ao_waveform = np.concatenate((mw_ao_waveform, rest_array))
                    # print('new len of laser thing:', len(laser_do_waveform))
                    # print('len of apd thing:', len(apd_do_waveform))
                else:
                    print('mw waveform is longer than laser waveform')
                    difference = len(mw_ao_waveform) - len(laser_do_waveform)
                    mw_ao_waveform = mw_ao_waveform[:-difference]
                    # print('new len of laser thing:', len(laser_do_waveform))
                    # print('new len of apd thing:', len(apd_do_waveform))

                ao_dict= {2: mw_ao_waveform}

            else:
                d_off = np.zeros(len(apd_do_waveform))
                ao_dict = {2: d_off}
            # print('len of mw_ao_waveform: ', len(mw_ao_waveform))
            # print('len of laser_do_waveform: ', len(laser_do_waveform))
            # print('len of apd_do_waveform: ', len(apd_do_waveform))

            do_dict = {laser_outchan: laser_do_waveform, apd_outchan: apd_do_waveform}
            dosequence.append((do_dict))
            aosequence.append((ao_dict))

        # print('dosequence: ', dosequence)
        # print('aosequence: ', aosequence)
        # load and play
        # print('len of dosequence: ', len(dosequence))
        # print('len of aosequence: ', len(aosequence))
        digital_output_map = {0: [0, 1]}  # directs to the right output channel analog channel0 directs to digital x0 and x1
        for aos, dos in zip(aosequence, dosequence):
            # Loads waveforms to the awg
            self._pulser.load_waveform(ao_waveform_dictionary=aos, do_waveform_dictionary=dos,
                                           digital_output_map=digital_output_map)

        if trigger == True:
            # This sequence waits for a software trigger to start playing and moving to the next step.
            self._pulser.configure_ORmask(card_idx, None)
            self._pulser.configure_ANDmask(card_idx, None)
            array=np.array([0x40000000], dtype=np.int64)
            stop_condition_list = np.repeat(array, step_count)
            # print('stop_condition_list: ', stop_condition_list)
            loops = np.ones((len(aosequence)), dtype=np.int64)

        else:
            # This sequence immediately starts after the sequences are loaded
            self._pulser.configure_ORmask(card_idx, 'immediate')
            self._pulser.configure_ANDmask(card_idx, None)
            loops = np.array([rep], dtype=np.int64)
            stop_condition_list = np.array([])

        # makes the waveforms repeat in loops
        self._pulser.load_sequence(  # digital_sequences=dosequence,
            loops_list=loops, segment_map=segment_map, stop_condition_list=stop_condition_list)
        self._pulser.start_card(card_idx)
        self._pulser.arm_trigger(card_idx)


    def delay_sweep_mw_1(self, clk_rate, laser_in, wait, laser_re, apd_stop, rest, start, steps, mw_wait_time, mw_pulse_time, rep, mw_pulse=bool, trigger=bool, segment_map=[]):

        '''
        check for the longest trace for sure
        This can play a laser sequence and the apd readout at the same time
        clk_rate:   something like 100
        laser_in:   length of initialisation laser pulse (it starts right at the beginning) in ms
        wait:       time between the end of the initialisation and the reinitialisation in ms
        laser_re:   length of the reinitialisation pulse in ms

        apd_stop:   time in ms where the apd is on (this needs to be varied to do the delay measurement)
        rest:       time in ms after the reinitialisation pulse to have a break before the next cycle starts
        start:      where in time should the apd pulse start shifting
        steps:      step size of the apd pulse sweep
        mw_wait_time: when should the mw pulse should happen?
        mw_pulse_time: length of the mw pulse in between the two laser pulses
        trigger:    True if the software trigger should be used and False if the loop array should define the length
        rep:        If there is no software trigger, rep gives the number of repetitions the waveform is played
        segment_map:If the order of different segments should be mixed (not there yet)
        '''
        #building the array to change the length of the apd pulse
        msplay = laser_in + wait + laser_re + rest
        clk_rate = MEGA(clk_rate)
        stop = (msplay - rest) - apd_stop
        # print('stop value: ', stop)

        #build array with all the desired apd_start values]
        apd_start = []
        step_count = math.floor((stop-start)/steps)
        # print('step_count: ', step_count)
        for i in range(step_count):
            apd_start.append((start+(i*steps)))
        # print('adpstart array: ', apd_start)
        laser_in *= 1e-3
        mw_wait_time *= 1e-3
        mw_pulse_time *= 1e-3
        laser_re *= 1e-3
        wait *= 1e-3
        apd_stop *= 1e-3
        rest *= 1e-3
        start *= 1e-3
        steps *= 1e-3
        msplay *= 1e-3
        apd_start = np.multiply(1e-3, apd_start)
        # print('apd_start: ', apd_start)

        #Laser waveform does not change the whole time
        first_pulse = np.ones(self._pulser.waveform_padding(clk_rate * laser_in))
        laser_wait_time = np.zeros(self._pulser.waveform_padding(clk_rate * wait))
        sec_pulse = np.ones(self._pulser.waveform_padding(clk_rate * laser_re))
        rest_array = np.zeros(self._pulser.waveform_padding(clk_rate * rest))
        laser_do_waveform = np.concatenate((first_pulse, laser_wait_time, sec_pulse, rest_array))
        dosequence = []
        aosequence = []
        laser_outchan = 0  # Laser
        apd_outchan = 1  # APD
        card_idx = 1
        for i in apd_start:
            apd_wait_time = np.zeros(self._pulser.waveform_padding(clk_rate * i))
            apd_on = np.ones(self._pulser.waveform_padding(clk_rate * apd_stop))
            apd_off = np.zeros(self._pulser.waveform_padding(clk_rate * (msplay - (i + apd_stop))))
            apd_do_waveform = np.concatenate((apd_wait_time, apd_on, apd_off))

            if mw_pulse == True:
                # if a mw pulse is desired:
                mw_wait = np.zeros(self._pulser.waveform_padding(clk_rate * mw_wait_time))
                pi_pulse_array = np.ones(self._pulser.waveform_padding(clk_rate * mw_pulse_time))
                mw_off = np.zeros(self._pulser.waveform_padding(clk_rate * msplay-(mw_wait_time+mw_pulse_time)))
                mw_ao_waveform = np.concatenate((mw_wait, pi_pulse_array, mw_off))
                # print('mw_ao_waveform: ', mw_ao_waveform)
                # print('len of mw_ao_waveform: ', len(mw_ao_waveform))

                # Make sure the arrays have the same length:
                # Check which array is the longest and append the two others:
                max_val = np.max(np.array([len(laser_do_waveform), len(apd_do_waveform), len(mw_ao_waveform)]))
                if max_val == 0:
                    # print('laser waveforms is the longest')
                    difference_mw = len(laser_do_waveform) - len(mw_ao_waveform)
                    difference_apd = len(laser_do_waveform) - len(apd_do_waveform)
                    rest_array_mw = np.zeros(int(difference_mw))
                    rest_array_apd = np.zeros(int(difference_apd))
                    mw_ao_waveform = np.concatenate((mw_ao_waveform, rest_array_mw))
                    apd_do_waveform = np.concatenate((apd_do_waveform, rest_array_apd))
                    # print('len of mw_ao_waveform: ', len(mw_ao_waveform))
                    # print('len of laser_do_waveform: ', len(laser_do_waveform))
                    # print('len of apd_do_waveform: ', len(apd_do_waveform))
                elif max_val ==1:
                    # print('apd waveforms is the longest')
                    difference_mw = len(apd_do_waveform) - len(mw_ao_waveform)
                    difference_laser = len(apd_do_waveform) - len(laser_do_waveform)
                    rest_array_mw = np.zeros(int(difference_mw))
                    rest_array_laser = np.zeros(int(difference_laser))
                    mw_ao_waveform = np.concatenate((mw_ao_waveform, rest_array_mw))
                    laser_do_waveform = np.concatenate((laser_do_waveform, rest_array_laser))
                    # print('len of mw_ao_waveform: ', len(mw_ao_waveform))
                    # print('len of laser_do_waveform: ', len(laser_do_waveform))
                    # print('len of apd_do_waveform: ', len(apd_do_waveform))
                else:
                    # print('mw waveforms is the longest')
                    difference_apd = len(mw_ao_waveform) - len(apd_do_waveform)
                    difference_laser = len(mw_ao_waveform) - len(laser_do_waveform)
                    rest_array_apd = np.zeros(int(difference_apd))
                    rest_array_laser = np.zeros(int(difference_laser))
                    apd_do_waveform = np.concatenate((apd_do_waveform, rest_array_apd))
                    laser_do_waveform = np.concatenate((laser_do_waveform, rest_array_laser))
                    # print('len of mw_ao_waveform: ', len(mw_ao_waveform))
                    # print('len of laser_do_waveform: ', len(laser_do_waveform))
                    # print('len of apd_do_waveform: ', len(apd_do_waveform))

                ao_dict= {2: mw_ao_waveform}

            else: #if no MW is needed
                # check which array is longer from laser and apd and use that one
                if len(apd_do_waveform) < len(laser_do_waveform):
                    # print('laser waveform is the longer')
                    difference = len(laser_do_waveform) - len(apd_do_waveform)
                    rest_array = np.zeros(int(difference))
                    apd_do_waveform = np.concatenate((apd_do_waveform, rest_array))
                    #Create zero array of same size:
                    d_off = np.zeros(len(laser_do_waveform))
                    # print('len of mw_ao_waveform: ', len(d_off))
                    # print('len of laser_do_waveform: ', len(laser_do_waveform))
                    # print('len of apd_do_waveform: ', len(apd_do_waveform))
                elif len(laser_do_waveform) < len(apd_do_waveform):
                    # print('apd waveform is the longer')
                    difference = len(apd_do_waveform) - len(laser_do_waveform)
                    rest_array = np.zeros(int(difference))
                    laser_do_waveform = np.concatenate((laser_do_waveform, rest_array))
                    #Create zero array of same size:
                    d_off = np.zeros(len(apd_do_waveform))
                    # print('len of mw_ao_waveform: ', len(d_off))
                    # print('len of laser_do_waveform: ', len(laser_do_waveform))
                    # print('len of apd_do_waveform: ', len(apd_do_waveform))
                else:
                    # print('apd waveform and laser_waveform have the same size')
                    d_off = np.zeros(len(apd_do_waveform))
                    # print('len of mw_ao_waveform: ', len(d_off))
                    # print('len of laser_do_waveform: ', len(laser_do_waveform))
                    # print('len of apd_do_waveform: ', len(apd_do_waveform))
                #Creates a dict for the analog channel if no MW is needed
                ao_dict = {2: d_off}

            #Now all the waveformes are prepared and have the same length
            do_dict = {laser_outchan: laser_do_waveform, apd_outchan: apd_do_waveform}
            dosequence.append((do_dict))
            aosequence.append((ao_dict))

        # print('dosequence: ', dosequence)
        # print('aosequence: ', aosequence)
        # load and play
        # print('len of dosequence: ', len(dosequence))
        # print('len of aosequence: ', len(aosequence))
        digital_output_map = {0: [0, 1]}  # directs to the right output channel analog channel0 directs to digital x0 and x1
        for aos, dos in zip(aosequence, dosequence):
            # Loads waveforms to the awg
            self._pulser.load_waveform(ao_waveform_dictionary=aos, do_waveform_dictionary=dos,
                                           digital_output_map=digital_output_map)

        if trigger == True:
            # This sequence waits for a software trigger to start playing and moving to the next step.
            self._pulser.configure_ORmask(card_idx, None)
            self._pulser.configure_ANDmask(card_idx, None)
            array=np.array([0x40000000], dtype=np.int64)
            stop_condition_list = np.repeat(array, step_count)
            # print('stop_condition_list: ', stop_condition_list)
            loops = np.ones((len(aosequence)), dtype=np.int64)

        else:
            # This sequence immediately starts after the sequences are loaded
            self._pulser.configure_ORmask(card_idx, 'immediate')
            self._pulser.configure_ANDmask(card_idx, None)
            loops = np.array([rep], dtype=np.int64)
            stop_condition_list = np.array([])
        # makes the waveforms repeat in loops
        self._pulser.load_sequence(  # digital_sequences=dosequence,
            loops_list=loops, segment_map=segment_map, stop_condition_list=stop_condition_list)
        self._pulser.start_card(card_idx)
        self._pulser.arm_trigger(card_idx)






    ########## old and weird things ############
=======

    #
    # # Load analog waveform to the awg and play it
    #     def load_analog_waveform(self,msplay,clk_rate):
    #     self.sine_wave()
    #     aosequence = [{2: first_seq_ch0, 3: first_seq_ch1},
    #                   {2: first_seq_ch1, 3: first_seq_ch0},
    #                   {2: first_seq_ch0, 3: first_seq_ch1},
    #                   {2: first_seq_ch1, 3: first_seq_ch0}
    #                   ]
>>>>>>> Stashed changes

=======

    #
    # # Load analog waveform to the awg and play it
    #     def load_analog_waveform(self,msplay,clk_rate):
    #     self.sine_wave()
    #     aosequence = [{2: first_seq_ch0, 3: first_seq_ch1},
    #                   {2: first_seq_ch1, 3: first_seq_ch0},
    #                   {2: first_seq_ch0, 3: first_seq_ch1},
    #                   {2: first_seq_ch1, 3: first_seq_ch0}
    #                   ]

>>>>>>> Stashed changes
=======

    #
    # # Load analog waveform to the awg and play it
    #     def load_analog_waveform(self,msplay,clk_rate):
    #     self.sine_wave()
    #     aosequence = [{2: first_seq_ch0, 3: first_seq_ch1},
    #                   {2: first_seq_ch1, 3: first_seq_ch0},
    #                   {2: first_seq_ch0, 3: first_seq_ch1},
    #                   {2: first_seq_ch1, 3: first_seq_ch0}
    #                   ]

>>>>>>> Stashed changes
=======

    #
    # # Load analog waveform to the awg and play it
    #     def load_analog_waveform(self,msplay,clk_rate):
    #     self.sine_wave()
    #     aosequence = [{2: first_seq_ch0, 3: first_seq_ch1},
    #                   {2: first_seq_ch1, 3: first_seq_ch0},
    #                   {2: first_seq_ch0, 3: first_seq_ch1},
    #                   {2: first_seq_ch1, 3: first_seq_ch0}
    #                   ]

>>>>>>> Stashed changes
=======

    #
    # # Load analog waveform to the awg and play it
    #     def load_analog_waveform(self,msplay,clk_rate):
    #     self.sine_wave()
    #     aosequence = [{2: first_seq_ch0, 3: first_seq_ch1},
    #                   {2: first_seq_ch1, 3: first_seq_ch0},
    #                   {2: first_seq_ch0, 3: first_seq_ch1},
    #                   {2: first_seq_ch1, 3: first_seq_ch0}
    #                   ]

>>>>>>> Stashed changes
    def build_block_element(self, clk_rate, len0, len1):
        """
        Creates an array of zeros and ones for a digital channel --> its just one step
        len0:           is the desired duration of zeros before the pulse in ms
        len1:           is the desired duration of ones in ms
        clk_rate:       samples/second of the AWG
        Returns the block_element as an array and the length of the array in units of samples
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
        '''
        msplay *= 1e-3
        clk_rate = MEGA(clk_rate)
        block_element=[]
        low = self.create_d_off(len0,clk_rate)
        high = self.create_d_on(len1,clk_rate)
        block_element = np.concatenate((low,high))
        block_element_len= len(block_element)
=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
        """
        clk_rate = 10 * clk_rate
        block_element = []
        low = self.create_d_off(len0, clk_rate)
        high = self.create_d_on(len1, clk_rate)
        block_element = np.concatenate((low, high))
        block_element_len = len(block_element)
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
        # print(len(block_element)
        return block_element, block_element_len

    def build_sequence(self, clk_rate, len0, len1, msplay):
        """
        Creates an array of block elements with the length of msplay for a digital channel
        len0:          is the desired duration of zeros before the pulse in ms
        len1:           is the desired duration of ones in ms
        clk_rate:       samples/second of the AWG
        msplay:          is the length of the whole sequence
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
        '''
        msplay *= 1e-3
        clk_rate = MEGA(clk_rate)
=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
        """
        clk_rate = 10 * clk_rate
>>>>>>> Stashed changes
        samples = self._pulser.waveform_padding(msplay * clk_rate)
        print("totalsamples: ", samples)

        sequence = []
        block_element = self.build_block_element(clk_rate, len0, len1)[0]
        print(
            "len of block element: ", self.build_block_element(clk_rate, len0, len1)[1]
        )
        rep = round(
            np.divide(samples, self.build_block_element(clk_rate, len0, len1)[1])
        )

        print("Reps: ", rep)
        sequence = np.repeat(block_element, rep)
        # sequence_len = len(sequence)
        # return sequence

    def build_waveform(
        self, msplay, clk_rate, loops, card_idx, segment_map=[]
    ):  # This works
        """
        Makes the AWG play a analog and digital waveform
        clk_rate:       samples/second of the AWG
        loops:            number of repetitions of each block
        """
        msplay *= 1e-3
        clk_rate = MEGA(clk_rate)
        # Makes sure samples has the right length
        samples = self._pulser.waveform_padding(msplay * clk_rate)
        # #Creates a time axis of the same length
        time_ax = np.linspace(0, samples / clk_rate, samples)
        # Just a sine with amplitude 1
        first_seq_ch0 = np.sin(2 * np.pi * time_ax / msplay)
        # linear function
        first_seq_ch1 = np.linspace(0, 1, len(time_ax))
        # Building a dict for anaolog signal which can be read by the AWG
        aosequence = [{2: first_seq_ch0, 3: first_seq_ch1}]

        # for digital signal
        do_chan = 1
        do1 = np.zeros(
            first_seq_ch0.shape
        )  # does zeros array with the length of the sine function which has the same length as samples
        do1[first_seq_ch0 > 0] = 1  # Creates one step function -_ first ones then zeros

        outchan = 0
        # Building a dict for anaolog signal which can be read by the AWG
        dosequence = [{outchan: do1}]

        digital_output_map = {1: [outchan]}  # directs to the right output channel
        for aos, dos in zip(aosequence, dosequence):
            # Loads waveforms to the awg
            self._pulser.load_waveform(
                ao_waveform_dictionary=aos,
                do_waveform_dictionary=dos,
                digital_output_map=digital_output_map,
            )

        # # This sequence immediately starts after the sequences are loaded
        # self.configure_ORmask(card_idx, 'immediate')
        # self.configure_ANDmask(card_idx, None)
        # stop_condition_list = np.array([])

        # This sequence waits for a software trigger to start playing and moving to the next step.
        self._pulser.configure_ORmask(card_idx, None)
        self._pulser.configure_ANDmask(card_idx, None)
        loops = np.ones((len(aosequence)), dtype=np.int64)
        SPCSEQ_ENDLOOPONTRIG = 0x40000000
        stop_condition_list = np.array([SPCSEQ_ENDLOOPONTRIG], dtype=np.int64)
        # makes the waveforms repeat in loops
        self._pulser.load_sequence(  # digital_sequences=dosequence,
            loops_list=loops,
            segment_map=segment_map,
            stop_condition_list=stop_condition_list,
        )

        self._pulser.start_card(card_idx)
        self._pulser.arm_trigger(card_idx)
        self._pulser.clear_all()

    def build_ao_signal(self, clk_rate, msplay):
        # this one should hava a dict as output
        msplay *= 1e-3
        clk_rate = MEGA(clk_rate)
        # Makes sure samples has the right length
        samples = self._pulser.waveform_padding(msplay * clk_rate)
        # #Creates a time axis of the same length
        time_ax = np.linspace(0, samples / clk_rate, samples)
        # Just a sine with amplitude 1
        sine = np.sin(2 * np.pi * time_ax / msplay)
        # linear function
        zeros = np.zeros(len(time_ax))
        # Building a dict for anaolog signal which can be read by the AWG
        aosequence = [{2: sine, 3: zeros}]
        return aosequence

    def build_do_signal(self, clk_rate, msplay):
        # this one should have a dict as output
        # this one should hava a dict as output
        msplay *= 1e-3
        clk_rate = MEGA(clk_rate)
        # Makes sure samples has the right length
        samples = self._pulser.waveform_padding(msplay * clk_rate)
        # #Creates a time axis of the same length
        time_ax = np.linspace(0, samples / clk_rate, samples)
        # Just a sine with amplitude 1
        sine = np.sin(2 * np.pi * time_ax / msplay)
        do1 = np.zeros(
            sine.shape
        )  # does zeros array with the length of the sine function which has the same length as samples
        do1[sine > 0] = 1  # Creates one step function -_ first ones then zeros
        outchan = 0
        # Building a dict for anaolog signal which can be read by the AWG
        dosequence = [{outchan: do1}]
        return dosequence

    def play_sequence(self, clk_rate, msplay, d_outchan, card_idx, segment_map=[]):
        # This one takes an ao input and do input as two dicts and makes them play
        # Get aosequence and dosequence
        aosequence = self.build_ao_signal(clk_rate, msplay)
        dosequence = self.build_do_signal(clk_rate, msplay)
        digital_output_map = {1: [d_outchan]}  # directs to the right output channel
        for aos, dos in zip(aosequence, dosequence):
            # Loads waveforms to the awg
            self._pulser.load_waveform(
                ao_waveform_dictionary=aos,
                do_waveform_dictionary=dos,
                digital_output_map=digital_output_map,
            )

        # # This sequence immediately starts after the sequences are loaded
        # self.configure_ORmask(card_idx, 'immediate')
        # self.configure_ANDmask(card_idx, None)
        # stop_condition_list = np.array([])

        # This sequence waits for a software trigger to start playing and moving to the next step.
        self._pulser.configure_ORmask(card_idx, None)
        self._pulser.configure_ANDmask(card_idx, None)
        loops = np.ones((len(aosequence)), dtype=np.int64)
        SPCSEQ_ENDLOOPONTRIG = 0x40000000
        stop_condition_list = np.array([SPCSEQ_ENDLOOPONTRIG], dtype=np.int64)
        # makes the waveforms repeat in loops
        self._pulser.load_sequence(  # digital_sequences=dosequence,
            loops_list=loops,
            segment_map=segment_map,
            stop_condition_list=stop_condition_list,
        )

        self._pulser.start_card(card_idx)
        self._pulser.arm_trigger(card_idx)
        self._pulser.clear_all()

<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
    # def input_parameters(self, msplay, clk_rate, ao_type=str, rep=100, number):
    #
    #     '''
    #     ao_type = shape of the analog signal: sine or square so far for each channel and every number of sequence
    #                 for now all of them are the same either square or sine so just str is fine
    #     msplay = time one sequence should be played
    #     clk_rate
    #     rep = number of repetitions per sequence (just in use for the no trigger mode)
    #     number = how many different sequences should be in one dict (played after each other or changed by trigger). This must be the same for analog and digital dict
    #     '''
    # Instead of that: just hardcode the channels to the different hardware modules
    # Microwave: analog channel 2 and 3
    # Laser: digital channel X0
    # Photon counter: digital channel X1 and X2 (for reference)

    # if n_ao <= 0 and n_do <= 0:
    #     print('Number of input channels is zero or < 0')
    # elif n_do > 3 or n_ao > 4:
    #     print('Too many channels: There are 3 digital channels and 4 analog channels')
    # else: # if the number of requested channels makes sense in general
    #     # Create digital outchan array
    #     if  n_do == 0: # --- ONLY ANALOG ---
    #         # No digital channels only analog channels a= 1,2,3 or 4
    #         print('no digital channel')
    #     elif n_do == 1:
    #         # one digital channel and n analog channels a= 0,1,2,3 or 4
    #         outchan = [0]
    #         print('1 digital channel: Use: X0')
    #     elif n_do == 2:
    #         outchan = [0, 1]
    #         # 2 digital channels and n analog channels a= 0,1,2,3 or 4
    #         print('2 digital channels: Use: X0 and X1')
    #     elif n_do == 3:
    #         outchan = [0, 1, 2]
    #         # 2 digital channels and n analog channels a= 0,1,2,3 or 4
    #         print('3 digital channel: Use: X0, X1 and X2')
    #     else:
    #         print('n_do is not a number between 0 and 3')
    #
    #     # Create analog channel array
    #     if n_ao == 0:
    #         # no analog channels
    #         channels = [0]
    #         print('no analog channel: Use: (analog channel 0)')
    #     elif n_ao == 1:
    #         # Case: Just one digital channel and no analog channels
    #         channels = [0]
    #         print('one analog channel: Use: analog channel 0')
    #     elif n_ao == 2:
    #         # Case: Just one digital channel and 1 analog channels
    #         channels = [0, 1]
    #         print('2 analog channels: Use analog channel 1 and 2')
    #     elif n_ao == 3:
    #         # Case: Just one digital channel and 2 analog channels
    #         channels = [1, 2, 3]
    #         print('3 analog channels: Use analog channel 1, 2 and 3')
    #     elif n_ao == 4:
    #         # Case: Just one digital channel and 3 analog channels
    #         channels = [1, 2, 3, 4]
    #         print('4 analog channels: Use analog channel 1, 2, 3 and 4')
    #     else:
    #         print('n_ao is not a number between 0 and 4')
    #
    #     print('channels:', channels)
    #     print('digital channels:', outchan)
    #     Things work until here!

    # #Check if the shape of ao_type=[] fits the number of ao_channels
    # # later when different ensembles are defined
    #     if ao_type == 'sine':
    #         waveform_ao = self.sine_wave(msplay, clk_rate)
    #         print('waveform: sine')
    #     elif ao_type == 'square':
    #         waveform_ao = self.create_d_on(msplay, clk_rate)
    #         print('waveform: square')
    #     else:
    #         print('sofar ao_type can only be sine or square')
    # else:
    #     print('there is no analog channel so we can only use square')
    #     #Still create a empty waveform:
    #     waveform_ao = self.create_d_off(msplay, clk_rate)
    #
    # # create a digital waveform if needed:
    # if n_do > 0:
    #     waveform_do = self.create_d_on(msplay, clk_rate)
    # else:
    #     print('no digital waveform needed')
    #
    # #  Now we don't care about channels anymore
    # #  now build the dict for ao and do
    # #  first ao_dict
    # v= len(channels)
    # aosequence = []     #This is an array with {} inside
    # for i in range(channels):
    #     part = {i: first_seq_ch0, 3: first_seq_ch1}
    # aosequence = [{2: first_seq_ch0, 3: first_seq_ch1},
    #               {2: first_seq_ch1, 3: first_seq_ch0},
    #               {2: first_seq_ch0, 3: first_seq_ch1},
    #               {2: first_seq_ch1, 3: first_seq_ch0}
    #               ]
    # aosequence = {k: v for k, v in (('a', 1), ('b', 2), ('c', 3))}
    # print(aosequence)
    #
    #
    #
    # #  then do_dict
    # d_off = self.create_d_off(msplay, clk_rate)
    # d_on = self.create_d_on(msplay, clk_rate)
    # # Building a dict for anaolog signal which can be read by the AWG
    # aosequence = [{channels[0]: d_on}]
    # dosequence = [{outchan: d_on}]
    #
    #
    #
    #

    #
    #     # Just digital channel in use
    #     # create an array with zeros for the analog dict
    #     d_off= self.create_d_off(msplay, clk_rate)
    #     aosequences
    #     #use only one analog channel
    #     channels_n=
    #     self.create_ao_dict(msplay,clk_rate, channels=[channels_n])
    #     self
    #     if n_do == 1:
    #     elif n_do == 2:
    #
    # elif n_ao ==

=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
    # def sequence_test(self, msecondsplay, loops, segment_map=[]):
    #     """
    #     Function for early debugging of the awg. Remove from the final class.
    #
    #     The steps to a successful initialization of a sequence are:
    #     - Define the clock rate
    #     - create the arrays that will represent the sequences
    #     - configure the trigger masks
    #     - load the sequence
    #     - start the card
    #     - arm the trigger
    #     """
    #     msecondsplay *= 1e-3
    #
    #     clk_rate = 1000 * 1000 * 100
    #
    #     card_idx = 1
    #     self._pulser.set_sample_rate(card_idx, clk_rate)
    #
    #     samples = self._pulser.waveform_padding(msecondsplay * clk_rate)
    #     time_ax = np.linspace(0, samples/clk_rate, samples)
    #
    #     first_seq_ch0 = np.sin(2*np.pi * time_ax / msecondsplay)
    #     first_seq_ch1 = np.linspace(0, 1, len(time_ax))
    #     aosequence = [{2: first_seq_ch0, 3:first_seq_ch1}]
    #
    #
    #     do_chan = 1
    #
    #
    #
    #     outchan = 0
    #     dosequence = [{outchan: do1}]
    #
    #     digital_output_map = {1: [outchan]}
    #     for aos,dos in zip(aosequence, dosequence):
    #         self._pulser.load_waveform(ao_waveform_dictionary=aos,do_waveform_dictionary=dos, digital_output_map=digital_output_map)
    #
    #     # # This sequence immediately starts after the sequences are loaded
    #     # self.configure_ORmask(card_idx, 'immediate')
    #     # self.configure_ANDmask(card_idx, None)
    #     # stop_condition_list = np.array([])
    #
    #     #This sequence waits for a software trigger to start playing and moving to the next step.
    #     self._pulser.configure_ORmask(card_idx, None)
    #     self._pulser.configure_ANDmask(card_idx, None)
    #     loops = np.ones((len(aosequence)), dtype=np.int64)
    #     stop_condition_list = np.array([0x40000000], dtype=np.int64)
    #
    #
    #     self._pulser.load_sequence(#digital_sequences=dosequence,
    #                        loops_list=loops, segment_map=segment_map, stop_condition_list=stop_condition_list)
    #
    #     self._pulser.start_card(card_idx)
    #     self._pulser.arm_trigger(card_idx)
    #     self._pulser.clear_all()
    #
    #
    # def test_seq_logic(self, msecondsplay,l1,l2,l3):
    #     loop_array= np.array([l1, l2, l3], dtype=np.int64)
    #     segment_map= [0, 1, 2]
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
    #     self._pulser.sequence_test(msecondsplay, loop_array, segment_map)








=======
    #     self._pulser.test_seq(msecondsplay, loop_array, segment_map)
>>>>>>> Stashed changes
=======
    #     self._pulser.test_seq(msecondsplay, loop_array, segment_map)
>>>>>>> Stashed changes
=======
    #     self._pulser.test_seq(msecondsplay, loop_array, segment_map)
>>>>>>> Stashed changes
=======
    #     self._pulser.test_seq(msecondsplay, loop_array, segment_map)
>>>>>>> Stashed changes
=======
    #     self._pulser.test_seq(msecondsplay, loop_array, segment_map)
>>>>>>> Stashed changes
