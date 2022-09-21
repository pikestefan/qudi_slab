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
import decimal
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


# Just hardcode the channels to the different hardware modules somewehere
# Microwave: analog channel 2 and 3 for I and Q parameter. For now we just need I
# Laser: digital channel X0
# Photon counter: digital channel X1 and X2 (for reference)


class Pulse(GenericLogic):
    pulsegenerator = Connector(interface='PulserInterface')

    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)

    def on_activate(self):
        # Initialisation performed during activation of the module

        # Get connectors
        self._pulser = self.pulsegenerator()
        active_channels = self._pulser.get_active_channels()
        self._pulser.clear_all()
        self._pulser.reset()
        return active_channels

    def on_deactivate(self):
        self._pulser.clear_all()
        self._pulser.reset()
        self._pulser.pulser_off()

    def set_up_pulser(self, card_idx, clk_rate):
        # def set_up_pulser(self):

        # clk_rate = MEGA(clk_rate)
        # What needs to be done:
        # pulser_off
        # 1. clear all
        # set_chan_amplitude
        # set_output_filters
        # set_sample_rate
        self._pulser.pulser_off()
        self._pulser.clear_all()
        # self._pulser.set_output_filters(channels, filter_active)
        self._pulser.set_sample_rate(card_idx, clk_rate)
        output_sample_rate = self._pulser.get_sample_rate(card_idx)
        return output_sample_rate

    def trigger(self):
        self._pulser.send_software_trig(1)


    # def laser_seq_old(self,seq_len, laser_times):
    #     '''
    #     len =   whole length of the sequence in ms
    #     laser_times = [laser_in, wait, laser_re] all in ms
    #             laser_in:   length of initialisation laser pulse (it starts right at the beginning) in ms
    #             wait:       time between the end of the initialisation and the reinitialisation in ms
    #             laser_re:   length of the reinitialisation pulse in ms
    #     '''
    #     seq_len *= 1e-3
    #     laser_times = np.multiply(laser_times, 1e-3)
    #     card_idx = 1
    #     # takes sampling rate which is there already
    #     clk_rate = self._pulser.get_sample_rate(card_idx)
    #     laser_channel = 0  # digital    X0 Input of the laser is 0-1 V
    #     laser_len = laser_times[0] + laser_times[1] + laser_times[2]
    #     if seq_len < laser_len:
    #         print('The total length needs to be larger that the laser_time sum')
    #     else:
    #         # Laser waveform does not change the whole time
    #         first_pulse = np.ones(self._pulser.waveform_padding(clk_rate * laser_times[0]))
    #         laser_wait_time = np.zeros(self._pulser.waveform_padding(clk_rate * laser_times[1]))
    #         sec_pulse = np.ones(self._pulser.waveform_padding(clk_rate * laser_times[2]))
    #         rest_array = np.zeros(self._pulser.waveform_padding(clk_rate * (seq_len - laser_len)))
    #         laser_do_waveform = np.concatenate((first_pulse, laser_wait_time, sec_pulse, rest_array))
    #     return laser_do_waveform
    def laser_seq(self, seq_len, laser_times):
        '''
        len =   whole length of the sequence in ms
        laser_times = [laser_in, wait, laser_re] all in ms
                laser_in:   length of initialisation laser pulse (it starts right at the beginning) in ms
                wait:       time between the end of the initialisation and the reinitialisation in ms
                laser_re:   length of the reinitialisation pulse in ms
        '''
        seq_len *= 1e-3
        laser_times = np.multiply(laser_times, 1e-3)
        card_idx = 1
        # takes sampling rate which is there already
        clk_rate = self._pulser.get_sample_rate(card_idx)
        laser_len = laser_times[0] + laser_times[1] + laser_times[2] # len of the important part in s
        len_rest = seq_len - laser_len # len of the rest in s
        if seq_len < laser_len:
            print('The total length needs to be larger that the laser_time sum')
        else:
            # Laser waveform does not change the whole time
            first_pulse = np.ones(int(clk_rate * laser_times[0]))
            laser_wait_time = np.zeros(int(clk_rate * laser_times[1]))
            sec_pulse = np.ones(int(clk_rate * laser_times[2]))
            rest_array = np.zeros(int(clk_rate * len_rest))
            laser_do_waveform = np.concatenate((first_pulse, laser_wait_time, sec_pulse, rest_array))
            len_laser_do_waveform = self._pulser.waveform_padding(len(laser_do_waveform))
            difference = len_laser_do_waveform - (clk_rate * seq_len)
            laser_do_waveform = np.concatenate((laser_do_waveform, np.zeros(int(difference))))
            # print('new length of waveform: ', len(laser_do_waveform))

        return laser_do_waveform
    def apd_seq(self,seq_len, apd_times):
        '''
        len =   whole length of the sequence in ms
        apd_times = [time to start, length] all in ms
        '''
        seq_len *= 1e-3
        apd_times = np.multiply(apd_times, 1e-3)
        apd_channel = 1  # Laser
        card_idx = 1
        #takes sampling rate which is there already
        clk_rate = self._pulser.get_sample_rate(card_idx)
        apd_len = apd_times[0] + apd_times[1]
        len_rest = seq_len - apd_len  # len of the rest in s
        if seq_len < apd_len:
            print('The total length needs to be larger that the apd_time sum')
        else:
            # apd waveform
            apd_wait = np.zeros(int(clk_rate * apd_times[0]))
            apd_on = np.ones(int(clk_rate * apd_times[1]))
            rest_array = np.zeros(int(clk_rate * (seq_len - apd_len)))
            apd_do_waveform = np.concatenate((apd_wait, apd_on, rest_array))
            len_apd_do_waveform = self._pulser.waveform_padding(len(apd_do_waveform))
            difference = len_apd_do_waveform - (clk_rate * seq_len)
            apd_do_waveform = np.concatenate((apd_do_waveform, np.zeros(int(difference))))
            # print('new length of waveform: ', len(apd_do_waveform))
        return apd_do_waveform


    def play_rabi(self, seq_len, laser_times, apd_times, mw_times, rep=100, trigger=bool):
        '''
        laser_times = [laser_in, wait, laser_re] all in ms
                laser_in:   length of initialisation laser pulse (it starts right at the beginning) in ms
                wait:       time between the end of the initialisation and the reinitialisation in ms
                laser_re:   length of the reinitialisation pulse in ms
        seq_len =   whole length of the sequence in ms
        apd_times = [time to start, length] all in ms
        mw_times = [mw_start_time, mw_len_0, mw_len_max, steps] all in ms
                mw_start_time:  time from the beginning when MW pulse should start, does not change
                mw_len_0:       minimal length of the mw pulse
                mw_len_max:     maximal length of the mw pulse
                steps:          stepsize for increasing the mw pulse duration
        '''
        dosequence = []
        aosequence = []
        segment_map = []
        mw_len = []
        #Channels
        laser_outchan = 0  # Laser
        apd_outchan = 1
        mw_outchan = 2  # Laser

        card_idx = 1
        # convert to seconds
        seq_len_s = seq_len * 1e-3
        # laser_times_s = np.multiply(laser_times, 1e-3)
        # apd_times_s = np.multiply(apd_times, 1e-3)
        mw_times_s = np.multiply(mw_times, 1e-3)

        clk_rate = self._pulser.get_sample_rate(card_idx)
        step_count = math.floor((mw_times[2] - mw_times[1]) / mw_times[3])
        print('step_count: ', step_count)
        #Create laser and apd waveform as they are costant
        laser_do_waveform = self.laser_seq(seq_len, laser_times)
        print('laser len: ', len(laser_do_waveform))
        apd_do_waveform = self.apd_seq(seq_len, apd_times)
        print('apd len: ', len(apd_do_waveform))

        for i in range(step_count):
            mw_len.append((mw_times[1] + (i * mw_times[3])))  # starting with minimal length
        # print ('mw_len: ', mw_len)
        mw_len = np.multiply(1e-3, mw_len)
        print('mw_len in sec: ', mw_len)
        for i in mw_len:
            # build the mw waveform for every
            if seq_len < (mw_times[0] + mw_times[2]):
                print('The total length needs to be larger that the apd_time sum')
            else:
                # build mw waveform for every possible length
                mw_wait = np.zeros(int(clk_rate * mw_times_s[0]))
                mw_on = np.ones(int(clk_rate * i))  # this changes
                rest_array = np.zeros(int(clk_rate * (seq_len_s - (mw_times_s[0] + i))))
                mw_ao_waveform = np.concatenate((mw_wait, mw_on, rest_array))
                # print('length of mw waveform before: ', len(mw_ao_waveform))
                len_mw_ao_waveform = self._pulser.waveform_padding(len(mw_ao_waveform))
                print('length of mw waveform after: ', len_mw_ao_waveform)
                difference = len_mw_ao_waveform - len(mw_ao_waveform)
                print('difference: ', difference)
                mw_ao_waveform = np.concatenate((mw_ao_waveform, np.zeros(int(difference))))
                # print('new length of mw waveform: ', len(mw_ao_waveform))
                # check_out = self.check_len(laser_do_waveform, apd_do_waveform, mw_ao_waveform, mw_pulse= True)
                # mw_ao_waveform = check_out[0]
                # apd_do_waveform = check_out[0]
                # laser_do_waveform = check_out[0]
                print('laser_do_waveform: ', len(laser_do_waveform))
                print('apd_do_waveform: ', len(apd_do_waveform))
                print('mw_ao_waveform: ', len(mw_ao_waveform))
                # Now build the  do_dict and ao_dict
                do_dict = {laser_outchan: laser_do_waveform, apd_outchan: apd_do_waveform}
                ao_dict = {mw_outchan: mw_ao_waveform}
                aosequence.append(ao_dict)  # maybe put ()
                dosequence.append(do_dict)  # maybe put ()
                # print('aosequence: ', aosequence)
                # print('dosequence: ', dosequence)

        digital_output_map = {0: [laser_outchan, apd_outchan]}
        # print('digital_output_map: ', digital_output_map)
        for aos, dos in zip(aosequence, dosequence):
            # Loads waveforms to the awg
            #check_out = self.check_len(laser_do_waveform, apd_do_waveform, mw_ao_waveform, mw_pulse=True)
            self._pulser.load_waveform(ao_waveform_dictionary=aos, do_waveform_dictionary=dos,
                                   digital_output_map=digital_output_map)

        if trigger == True:
            # This sequence waits for a software trigger to start playing and moving to the next step.
            self._pulser.configure_ORmask(card_idx, None)
            self._pulser.configure_ANDmask(card_idx, None)
            array = np.array([0x40000000], dtype=np.int64)
            stop_condition_list = np.repeat(array, step_count)
            loops = np.ones((len(aosequence)), dtype=np.int64)

        else:
            # This sequence immediately starts after the sequences are loaded
            self._pulser.configure_ORmask(card_idx, 'immediate')
            self._pulser.configure_ANDmask(card_idx, None)
            loop_array = np.array([rep], dtype=np.int64)
            loops = np.repeat(loop_array, step_count)
            stop_condition_list = np.array([])

            # makes the waveforms repeat in loops
        self._pulser.load_sequence(  # digital_sequences=dosequence,
            loops_list=loops, segment_map=segment_map, stop_condition_list=stop_condition_list)
        self._pulser.start_card(card_idx)
        self._pulser.arm_trigger(card_idx)

    def setup_and_play(self, laser_times, mw_times, apd_times, rep=100, rest=100, clk_rate=100, mw_pulse=bool,
                       trigger=bool):
        '''
            laser_times, mw_times, apd_times should be arrays with the right values inside
            laser_times = [laser_in, wait, laser_re] all in ms
                laser_in:   length of initialisation laser pulse (it starts right at the beginning) in ms
                wait:       time between the end of the initialisation and the reinitialisation in ms
                laser_re:   length of the reinitialisation pulse in ms
            mw_times = [mw_wait_time, mw_pulse_time] all in ms
                mw_wait_time: time between the end of the first laser pulse and the beginning of the microwave (delay)
                mw_pulse_time: length of the mw pulse in between the two laser pulses
            apd_times = [steps, apd_width]
                steps:      step size of the apd pulse sweep
                apd_width:   time in ms where the apd is on (length of the pulse)
            rep:        If there is no software trigger, rep gives the number of repetitions the waveform is played
            rest:       Time between single sequences (just to be able to distinguish between beginning of 1 laser pulse and end of second laser lulse)
            clk_rate:  set to 100 if None
            mw_pulse:   True: plays MW pulse and False sets the MW sequence to zero
            trigger:    True if the software trigger should be used and False if the loop array should define the length

        These parameters work:
        laser_in=2,wait=1,laser_re=7,
        steps=1, apd_width=2
         rest=100,
         mw_wait_time=0.5, mw_pulse_time=0.25,
         rep=10,
         mw_pulse=True, trigger=True)
        '''

        laser_times = np.multiply(laser_times, 1e-3)
        mw_times = np.multiply(mw_times, 1e-3)
        apd_times = np.multiply(apd_times, 1e-3)

        # Sofar we just use card1
        card_idx = 1

        # Set up the right channels:
        # Microwave: analog channel 2 and 3 for I and Q parameter. For now we just need I
        # Laser: digital channel X0
        # Photon counter: digital channel X1 and X2 (for reference)
        mw_I = 2  # analog     A0 for card 1   Input at MW generator is +- 0.5 V
        mw_Q = 3  # analog     A1 for card 1   --> add later
        laser_channel = 0  # digital    X0 Input of the laser is 0-1 V
        apd_channel1 = 1  # digital    X1
        # build an array for that [0,2,1]
        channels = np.array([laser_channel, mw_I, apd_channel1])

        # if not given differently put the clk_rate to 100ms
        if clk_rate == 100:
            clk_rate = MEGA(100)
            # print('clk_rate not changed')
            print('clk_rate in 1/s: ', clk_rate)
        else:
            clk_rate = MEGA(clk_rate)
            # print('clk_rate changed')
            print('clk_rate in 1/s: ', clk_rate)

        # This one sets up the right clk_rate
        self.set_up_pulser(card_idx, clk_rate)
        print('AWG is ready')

        if rest == 100:
            rest = 100 * 1e-3
            print('break in between single sequences set to 100ms')
        else:
            rest *= 1e-3
            print('break in between single sequences set to the new value')
        if rep == 100:
            print('set rep to 100')
        else:
            print('set rep to the new value')

        segment_map = []
        # calculate the whole length of the waveform --> seq_len
        # The parameters for the laser determine how long the waveform is
        seq_len = laser_times[0] + laser_times[1] + laser_times[2] + rest  # in s

        # use all these parameters and put it to the general sweep function:
        self.delay_sweep_general(card_idx, clk_rate, seq_len, laser_times[0], laser_times[1], laser_times[2], rest,
                                 apd_times[0], apd_times[1], mw_times[0],
                                 mw_times[1], rep, channels, mw_pulse, trigger, segment_map)

    # This is the main general block...works fine with laser, mw and apd
    def delay_sweep_general(self, card_idx, clk_rate, seq_len, laser_in, wait, laser_re, rest, steps, apd_width,
                            mw_wait_time, mw_pulse_time, rep, channels, mw_pulse=bool, trigger=bool, segment_map=[]):
        '''
            This can play a laser sequence, microwaves and the apd readout at the same time
            card_idx:   sofar just 1
            clk_rate:   something like 100
            seq_len:    length of one whole sequence
            laser times:
                laser_in:   length of initialisation laser pulse (it starts right at the beginning) in ms
                wait:       time between the end of the initialisation and the reinitialisation in ms
                laser_re:   length of the reinitialisation pulse in ms
            rest:       time in ms after the reinitialisation pulse to have a break before the next sequences starts
            apd times
                start:      where in time should the apd pulse start shifting --> this should be the time as the start of the laser reinitialisation = laser_in + wait
                steps:      step size of the apd pulse sweep ( this is in ms right?)
                apd_width:   time in ms where the apd is on (length of the pulse)
            MW times
                mw_wait_time: time between the end of the first laser pulse and the beginning of the microwave (delay)
                mw_pulse_time: length of the mw pulse in between the two laser pulses
            rep:        If there is no software trigger, rep gives the number of repetitions the waveform is played
            channels:   Connetors of the AWG (hardcoded)
            mw_pulse:   True: sends MW pulse, False: sends zeros instead
            trigger:    True if the software trigger should be used and False if the loop array should define the length
            segment_map:If the order of different segments should be mixed (not there yet)
            '''

        apd_start_val = laser_in + wait
        stop = (seq_len - rest) - apd_width
        laser_outchan = channels[0]  # Laser
        apd_outchan = channels[2]  # APD
        mw_outchan = channels[1]

        # build array with all the desired apd_start values
        apd_start = []
        step_count = math.floor((stop - apd_start_val) / steps)
        for i in range(step_count):
            apd_start.append((apd_start_val + (i * steps)))

        # creating empty arrays
        dosequence = []
        aosequence = []

        # Laser waveform does not change the whole time
        first_pulse = np.ones(self._pulser.waveform_padding(clk_rate * laser_in))
        laser_wait_time = np.zeros(self._pulser.waveform_padding(clk_rate * wait))
        sec_pulse = np.ones(self._pulser.waveform_padding(clk_rate * laser_re))
        rest_array = np.zeros(self._pulser.waveform_padding(clk_rate * rest))
        laser_do_waveform = np.concatenate((first_pulse, laser_wait_time, sec_pulse, rest_array))

        for i in apd_start:
            apd_wait_time = np.zeros(self._pulser.waveform_padding(clk_rate * i))
            apd_on = np.ones(self._pulser.waveform_padding(clk_rate * apd_width))
            apd_off = np.zeros(self._pulser.waveform_padding(clk_rate * (seq_len - (i + apd_width))))
            apd_do_waveform = np.concatenate((apd_wait_time, apd_on, apd_off))

            if mw_pulse == True:
                mw_wait = np.zeros(self._pulser.waveform_padding(clk_rate * (laser_in + mw_wait_time)))
                mw_pulse_array = np.ones(self._pulser.waveform_padding(clk_rate * mw_pulse_time))
                mw_off = np.zeros(self._pulser.waveform_padding(clk_rate * seq_len - (laser_in + mw_wait_time + mw_pulse_time)))
                mw_pulse_array = np.subtract(mw_pulse_array, 0.5)  # because the MW takes only +-0.5V
                mw_ao_waveform = np.concatenate((mw_wait, mw_pulse_array, mw_off))

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

                elif max_val == 1:
                    # print('apd waveforms is the longest')
                    difference_mw = len(apd_do_waveform) - len(mw_ao_waveform)
                    difference_laser = len(apd_do_waveform) - len(laser_do_waveform)
                    rest_array_mw = np.zeros(int(difference_mw))
                    rest_array_laser = np.zeros(int(difference_laser))
                    mw_ao_waveform = np.concatenate((mw_ao_waveform, rest_array_mw))
                    laser_do_waveform = np.concatenate((laser_do_waveform, rest_array_laser))

                else:
                    # print('mw waveforms is the longest')
                    difference_apd = len(mw_ao_waveform) - len(apd_do_waveform)
                    difference_laser = len(mw_ao_waveform) - len(laser_do_waveform)
                    rest_array_apd = np.zeros(int(difference_apd))
                    rest_array_laser = np.zeros(int(difference_laser))
                    apd_do_waveform = np.concatenate((apd_do_waveform, rest_array_apd))
                    laser_do_waveform = np.concatenate((laser_do_waveform, rest_array_laser))

                ao_dict = {mw_outchan: mw_ao_waveform}

            else:  # if no MW is needed
                # check which array is longer laser or apd
                if len(apd_do_waveform) < len(laser_do_waveform):
                    # print('laser waveform is longer than apd waveform')
                    difference = len(laser_do_waveform) - len(apd_do_waveform)
                    rest_array = np.zeros(int(difference))
                    apd_do_waveform = np.concatenate((apd_do_waveform, rest_array))
                    # Create zero array of same size:
                    d_off = np.zeros(len(laser_do_waveform))
                elif len(laser_do_waveform) < len(apd_do_waveform):
                    # print('apd waveform is the longer')
                    difference = len(apd_do_waveform) - len(laser_do_waveform)
                    rest_array = np.zeros(int(difference))
                    laser_do_waveform = np.concatenate((laser_do_waveform, rest_array))
                    # Create zero array of same size:
                    d_off = np.zeros(len(apd_do_waveform))
                else:
                    print('apd waveform and laser_waveform have the same size')
                    # d_off = np.zeros(len(apd_do_waveform))
                # Creates a dict for the analog channel if no MW is needed
                ao_dict = {mw_outchan: d_off}

            # Now all the waveformes are prepared and have the same length
            do_dict = {laser_outchan: laser_do_waveform, apd_outchan: apd_do_waveform}
            dosequence.append((do_dict))
            aosequence.append((ao_dict))

        digital_output_map = {0: [laser_outchan, apd_outchan]}  # directs to the right output channel analog channel0 directs to digital x0 and x1

        for aos, dos in zip(aosequence, dosequence):
            # Loads waveforms to the awg
            self._pulser.load_waveform(ao_waveform_dictionary=aos, do_waveform_dictionary=dos,
                                       digital_output_map=digital_output_map)

        if trigger == True:
            # This sequence waits for a software trigger to start playing and moving to the next step.
            self._pulser.configure_ORmask(card_idx, None)
            self._pulser.configure_ANDmask(card_idx, None)
            array = np.array([0x40000000], dtype=np.int64)
            stop_condition_list = np.repeat(array, step_count)
            loops = np.ones((len(aosequence)), dtype=np.int64)

        else:
            # This sequence immediately starts after the sequences are loaded
            self._pulser.configure_ORmask(card_idx, 'immediate')
            self._pulser.configure_ANDmask(card_idx, None)
            loop_array = np.array([rep], dtype=np.int64)
            loops = np.repeat(loop_array, step_count)
            stop_condition_list = np.array([])

            # makes the waveforms repeat in loops
        self._pulser.load_sequence(  # digital_sequences=dosequence,
            loops_list=loops, segment_map=segment_map, stop_condition_list=stop_condition_list)
        self._pulser.start_card(card_idx)
        self._pulser.arm_trigger(card_idx)

    def test_sequences(self, seq_len, clk_rate=100, rep=100, laser=bool, mw=bool, apd=bool, trigger=bool):
        '''
        This is just for testing the different channels. Just sends ones ore zeros to the 3 different channels
        clk_rate:   something like 100
        msplay:     time of one sequence in ms
        trigger:    True if the software trigger should be used and False if the loop array should define the length
        rep:        If there is no software trigger, rep gives the number of repetitions the sequence is played
        '''

        # Sofar we just use card1
        card_idx = 1
        seq_len *= 1e-3

        if clk_rate == 100:
            clk_rate = MEGA(100)
            # print('clk_rate not changed')
            print('clk_rate in 1/s: ', clk_rate)
        else:
            clk_rate = MEGA(clk_rate)
            # print('clk_rate changed')
            print('clk_rate in 1/s: ', clk_rate)

        # This one sets up the right clk_rate
        self.set_up_pulser(card_idx, clk_rate)
        print('setup of the pulser is done')

        # Set up the right channels:
        # Microwave: analog channel 2 and 3 for I and Q parameter. For now we just need I
        # Laser: digital channel X0
        # Photon counter: digital channel X1 and X2 (for reference)
        mw_I = 2  # analog     A0 for card 1   Input at MW generator is +- 0.5 V
        mw_Q = 3  # analog     A1 for card 1   --> add later
        laser_channel = 0  # digital    X0 Input of the laser is 0-1 V
        apd_channel1 = 1  # digital    X1
        channels = np.array([laser_channel, mw_I, apd_channel1])    #[0,2,1]

        if rep == 100:
            print('set rep to 100')
        else:
            print('set rep to something else')

        segment_map = []

        # Create arrays:
        # For Laser
        if laser == True:
            laser_sequence = np.ones(self._pulser.waveform_padding(clk_rate * seq_len))
        else:
            laser_sequence = np.zeros(self._pulser.waveform_padding(clk_rate * seq_len))

        # For MW
        if mw == True:
            mw_sequence = np.multiply(np.ones(self._pulser.waveform_padding(clk_rate * seq_len)), 0.5)
        else:
            mw_sequence = np.zeros(self._pulser.waveform_padding(clk_rate * seq_len))

        # For APD
        if apd == True:
            apd_sequence = np.ones(self._pulser.waveform_padding(clk_rate * seq_len))*3
        else:
            apd_sequence = np.zeros(self._pulser.waveform_padding(clk_rate * seq_len))

        #create some zeros
        zero = np.zeros(self._pulser.waveform_padding(clk_rate * seq_len))

        aosequence = [{channels[1]: mw_sequence}, {channels[0]: zero}]
        dosequence = [{channels[0]: laser_sequence, channels[2]: apd_sequence},
                      {channels[0]: zero, channels[2]: zero}]

        # Play it 2 times because it doesn't work with only one waveform
        digital_output_map = {0: [channels[0]], 1: [channels[2]]}  # directs to the right output channel, analog channel0 directs to digital x0 and x1
        # print('digital output map: ', digital_output_map)
        for aos, dos in zip(aosequence, dosequence):
            # Loads waveforms to the awg
            self._pulser.load_waveform(ao_waveform_dictionary=aos, do_waveform_dictionary=dos,
                                       digital_output_map=digital_output_map)
        if trigger == True:
            # This sequence waits for a software trigger to start playing and moving to the next step.
            self._pulser.configure_ORmask(card_idx, None)
            self._pulser.configure_ANDmask(card_idx, None)
            array = np.array([0x40000000], dtype=np.int64)
            loops = np.ones((len(aosequence)), dtype=np.int64)
            stop_condition_list = np.repeat(array,len(aosequence))
            # print('loops: ', loops)
            # print('stop con list: ', stop_condition_list)

        else:
        # This sequence immediately starts after the sequences are loaded
            self._pulser.configure_ORmask(card_idx, 'immediate')
            self._pulser.configure_ANDmask(card_idx, None)
            array = np.array([0x40000000], dtype=np.int64)
            stop_condition_list = np.repeat(array, len(aosequence))
            loops= np.array([100,100], dtype=np.int64)
            # print ('loops: ', loops)

        # makes the waveforms repeat in loops
        self._pulser.load_sequence(  # digital_sequences=dosequence,
            loops_list=loops, segment_map=segment_map, stop_condition_list=stop_condition_list)
        self._pulser.start_card(card_idx)
        self._pulser.arm_trigger(card_idx)

    def check_len(self, laser_do_waveform, apd_do_waveform, mw_ao_waveform, mw_pulse=bool):
        # checks and adjusts the length of the different waveform because waveform_padding creates a mess
        if mw_pulse == True:
            val_list = [len(laser_do_waveform), len(apd_do_waveform), len(mw_ao_waveform)]
            max_val = max(val_list)
            # print ('max_val: ', max_value)
            max_index = val_list.index(max_val)
            print(max_index)
            if max_index == 0:
                print('laser waveforms is the longest')
                difference_mw = len(laser_do_waveform) - len(mw_ao_waveform)
                difference_apd = len(laser_do_waveform) - len(apd_do_waveform)
                rest_array_mw = np.zeros(int(difference_mw))
                rest_array_apd = np.zeros(int(difference_apd))
                mw_ao_waveform = np.concatenate((mw_ao_waveform, rest_array_mw))
                apd_do_waveform = np.concatenate((apd_do_waveform, rest_array_apd))
                d_off = []
            elif max_index == 1:
                print('apd waveforms is the longest')
                difference_mw = len(apd_do_waveform) - len(mw_ao_waveform)
                difference_laser = len(apd_do_waveform) - len(laser_do_waveform)
                rest_array_mw = np.zeros(int(difference_mw))
                rest_array_laser = np.zeros(int(difference_laser))
                mw_ao_waveform = np.concatenate((mw_ao_waveform, rest_array_mw))
                laser_do_waveform = np.concatenate((laser_do_waveform, rest_array_laser))
                d_off = []
            else:
                print('mw waveforms is the longest')
                difference_apd = len(mw_ao_waveform) - len(apd_do_waveform)
                difference_laser = len(mw_ao_waveform) - len(laser_do_waveform)
                rest_array_apd = np.zeros(int(difference_apd))
                rest_array_laser = np.zeros(int(difference_laser))
                apd_do_waveform = np.concatenate((apd_do_waveform, rest_array_apd))
                laser_do_waveform = np.concatenate((laser_do_waveform, rest_array_laser))
                d_off = []
        else:
            if len(apd_do_waveform) < len(laser_do_waveform):
                print('laser waveform is longer than apd waveform')
                difference = len(laser_do_waveform) - len(apd_do_waveform)
                rest_array = np.zeros(int(difference))
                apd_do_waveform = np.concatenate((apd_do_waveform, rest_array))
                # Create zero array of same size:
                d_off = np.zeros(len(laser_do_waveform))
            elif len(laser_do_waveform) < len(apd_do_waveform):
                print('apd waveform is the longer')
                difference = len(apd_do_waveform) - len(laser_do_waveform)
                rest_array = np.zeros(int(difference))
                laser_do_waveform = np.concatenate((laser_do_waveform, rest_array))
                # Create zero array of same size:
                d_off = np.zeros(len(apd_do_waveform))
            else:
                print('apd_waveform and laser_waveform have the same size')
                d_off = np.zeros(len(apd_do_waveform))
        return mw_ao_waveform, apd_do_waveform, laser_do_waveform, d_off

    def setup_and_play_new(self, laser_times, mw_times, apd_times, rep=100, rest=100, clk_rate=100, mw_pulse=bool,
                       trigger=bool):
        '''
            laser_times, mw_times, apd_times should be arrays with the right values inside
            laser_times = [laser_in, wait, laser_re] all in ms
                laser_in:   length of initialisation laser pulse (it starts right at the beginning) in ms
                wait:       time between the end of the initialisation and the reinitialisation in ms
                laser_re:   length of the reinitialisation pulse in ms
            mw_times = [mw_wait_time, mw_pulse_time] all in ms
                mw_wait_time: time between the end of the first laser pulse and the beginning of the microwave (delay)
                mw_pulse_time: length of the mw pulse in between the two laser pulses
            apd_times = [steps, apd_width]
                steps:      step size of the apd pulse sweep
                apd_width:   time in ms where the apd is on (length of the pulse)
            rep:        If there is no software trigger, rep gives the number of repetitions the waveform is played
            rest:       Time between single sequences (just to be able to distinguish between beginning of 1 laser pulse and end of second laser lulse)
            clk_rate:  set to 100 if None
            mw_pulse:   True: plays MW pulse and False sets the MW sequence to zero
            trigger:    True if the software trigger should be used and False if the loop array should define the length

        These parameters work:
        laser_in=2,wait=1,laser_re=7,
        steps=1, apd_width=2
         rest=100,
         mw_wait_time=0.5, mw_pulse_time=0.25,
         rep=10,
         mw_pulse=True, trigger=True)
        '''

        laser_times = np.multiply(laser_times, 1e-3)
        mw_times = np.multiply(mw_times, 1e-3)
        apd_times = np.multiply(apd_times, 1e-3)

        # Sofar we just use card1
        card_idx = 1

        # Set up the right channels:
        # Microwave: analog channel 2 and 3 for I and Q parameter. For now we just need I
        # Laser: digital channel X0
        # Photon counter: digital channel X1 and X2 (for reference)
        mw_I = 2  # analog     A0 for card 1   Input at MW generator is +- 0.5 V
        mw_Q = 3  # analog     A1 for card 1   --> add later
        laser_channel = 0  # digital    X0 Input of the laser is 0-1 V
        apd_channel1 = 1  # digital    X1
        # build an array for that [0,2,1]
        channels = np.array([laser_channel, mw_I, apd_channel1])

        # if not given differently put the clk_rate to 100ms
        if clk_rate == 100:
            clk_rate = MEGA(100)
            # print('clk_rate not changed')
            print('clk_rate in 1/s: ', clk_rate)
        else:
            clk_rate = MEGA(clk_rate)
            # print('clk_rate changed')
            print('clk_rate in 1/s: ', clk_rate)

        # This one sets up the right clk_rate
        self.set_up_pulser(card_idx, clk_rate)
        print('AWG is ready')

        if rest == 100:
            rest = 100 * 1e-3
            print('break in between single sequences set to 100ms')
        else:
            rest *= 1e-3
            print('break in between single sequences set to the new value')
        if rep == 100:
            print('set rep to 100')
        else:
            print('set rep to the new value')

        segment_map = []
        # calculate the whole length of the waveform --> seq_len
        # The parameters for the laser determine how long the waveform is
        seq_len = laser_times[0] + laser_times[1] + laser_times[2] + rest  # in s
        apd_start_val = laser_times[0] + laser_times[1]
        stop = (seq_len - rest) - apd_times[1]

        # build array with all the desired apd_start values
        apd_start = []
        step_count = math.floor((stop - apd_start_val) / apd_times[0])
        for i in range(step_count):
            apd_start.append((apd_start_val + (i * apd_times[0])))

        # use all these parameters and put it to the general sweep function:
        self.delay_sweep_general_new(card_idx, clk_rate, seq_len, laser_times[0], laser_times[1], laser_times[2], rest,
                                 apd_start, apd_times[1], mw_times[0],
                                 mw_times[1],step_count, rep, channels, mw_pulse, trigger, segment_map)


    # This is the main general block...works fine with laser, mw and apd NEW
    def delay_sweep_general_new(self, card_idx, clk_rate, seq_len, laser_in, wait, laser_re, rest, apd_start, apd_width,
                            mw_wait_time, mw_pulse_time, step_count, rep, channels, mw_pulse=bool, trigger=bool, segment_map=[]):
        '''
            This can play a laser sequence, microwaves and the apd readout at the same time
            card_idx:   sofar just 1
            clk_rate:   something like 100
            seq_len:    length of one whole sequence
            laser times:
                laser_in:   length of initialisation laser pulse (it starts right at the beginning) in ms
                wait:       time between the end of the initialisation and the reinitialisation in ms
                laser_re:   length of the reinitialisation pulse in ms
            rest:       time in ms after the reinitialisation pulse to have a break before the next sequences starts
            apd times
                start:      where in time should the apd pulse start shifting --> this should be the time as the start of the laser reinitialisation = laser_in + wait
                steps:      step size of the apd pulse sweep ( this is in ms right?)
                apd_width:   time in ms where the apd is on (length of the pulse)
            MW times
                mw_wait_time: time between the end of the first laser pulse and the beginning of the microwave (delay)
                mw_pulse_time: length of the mw pulse in between the two laser pulses
            rep:        If there is no software trigger, rep gives the number of repetitions the waveform is played
            channels:   Connetors of the AWG (hardcoded)
            mw_pulse:   True: sends MW pulse, False: sends zeros instead
            trigger:    True if the software trigger should be used and False if the loop array should define the length
            segment_map:If the order of different segments should be mixed (not there yet)
            '''

        laser_outchan = channels[0]  # Laser
        apd_outchan = channels[2]  # APD
        mw_outchan = channels[1]

        # creating empty arrays
        dosequence = []
        aosequence = []

        # Laser waveform does not change the whole time
        first_pulse = np.ones(self._pulser.waveform_padding(clk_rate * laser_in))
        laser_wait_time = np.zeros(self._pulser.waveform_padding(clk_rate * wait))
        sec_pulse = np.ones(self._pulser.waveform_padding(clk_rate * laser_re))
        rest_array = np.zeros(self._pulser.waveform_padding(clk_rate * rest))
        laser_do_waveform = np.concatenate((first_pulse, laser_wait_time, sec_pulse, rest_array))

        for i in apd_start:
            apd_wait_time = np.zeros(self._pulser.waveform_padding(clk_rate * i))
            apd_on = np.ones(self._pulser.waveform_padding(clk_rate * apd_width))
            apd_off = np.zeros(self._pulser.waveform_padding(clk_rate * (seq_len - (i + apd_width))))
            apd_do_waveform = np.concatenate((apd_wait_time, apd_on, apd_off))

            if mw_pulse == True:
                mw_wait = np.zeros(self._pulser.waveform_padding(clk_rate * (laser_in + mw_wait_time)))
                mw_pulse_array = np.ones(self._pulser.waveform_padding(clk_rate * mw_pulse_time))
                mw_off = np.zeros(
                    self._pulser.waveform_padding(clk_rate * seq_len - (laser_in + mw_wait_time + mw_pulse_time)))
                mw_pulse_array = np.subtract(mw_pulse_array, 0.5)  # because the MW takes only +-0.5V
                mw_ao_waveform = np.concatenate((mw_wait, mw_pulse_array, mw_off))

                # Make sure the arrays have the same length:
                # Check which array is the longest and append the two others:
                check_out = self.check_len(laser_do_waveform,apd_do_waveform,mw_ao_waveform, mw_pulse= True)
                # print ('check_out: ', check_out)
                # order: mw_ao_waveform, apd_do_waveform, laser_do_waveform, d_off
                mw_ao_waveform = check_out[0]
                apd_do_waveform = check_out[1]
                laser_do_waveform = check_out[2]
                d_off = check_out[3]
                ao_dict = {mw_outchan: mw_ao_waveform}

            else:  # if no MW is needed
                # check which array is longer laser or apd
                check_out = self.check_len(laser_do_waveform, apd_do_waveform, mw_ao_waveform=[], mw_pulse=False)
                # print('check_out: ', check_out)
                # order: mw_ao_waveform, apd_do_waveform, laser_do_waveform, d_off
                mw_ao_waveform = check_out[0]
                apd_do_waveform = check_out[1]
                laser_do_waveform = check_out[2]
                d_off = check_out[3]
                # Creates a dict for the analog channel if no MW is needed
                ao_dict = {mw_outchan: d_off}

            # Now all the waveformes are prepared and have the same length
            do_dict = {laser_outchan: laser_do_waveform, apd_outchan: apd_do_waveform}
            dosequence.append((do_dict))
            aosequence.append((ao_dict))

        digital_output_map = {0: [laser_outchan,
                                  apd_outchan]}  # directs to the right output channel analog channel0 directs to digital x0 and x1

        for aos, dos in zip(aosequence, dosequence):
            # Loads waveforms to the awg
            self._pulser.load_waveform(ao_waveform_dictionary=aos, do_waveform_dictionary=dos,
                                       digital_output_map=digital_output_map)

        if trigger == True:
            # This sequence waits for a software trigger to start playing and moving to the next step.
            self._pulser.configure_ORmask(card_idx, None)
            self._pulser.configure_ANDmask(card_idx, None)
            array = np.array([0x40000000], dtype=np.int64)
            stop_condition_list = np.repeat(array, step_count)
            loops = np.ones((len(aosequence)), dtype=np.int64)

        else:
            # This sequence immediately starts after the sequences are loaded
            self._pulser.configure_ORmask(card_idx, 'immediate')
            self._pulser.configure_ANDmask(card_idx, None)
            loop_array = np.array([rep], dtype=np.int64)
            loops = np.repeat(loop_array, step_count)
            stop_condition_list = np.array([])

            # makes the waveforms repeat in loops
        self._pulser.load_sequence(  # digital_sequences=dosequence,
            loops_list=loops, segment_map=segment_map, stop_condition_list=stop_condition_list)
        self._pulser.start_card(card_idx)
        self._pulser.arm_trigger(card_idx)


    def stop_awg(self): # Makes everything stop
        self._pulser.clear_all()
        self._pulser.reset()

    # This is the main block...not that general but works with the laser
    def delay_sweep_mw(self, clk_rate, laser_in, wait, laser_re, apd_width, rest, start, steps, mw_wait_time,
                       mw_pulse_time, rep, mw_pulse=bool, trigger=bool, segment_map=[]):

        '''
        check for the longest trace for sure
        This can play a laser sequence and the apd readout at the same time
        clk_rate:   something like 100
        laser_in:   length of initialisation laser pulse (it starts right at the beginning) in ms
        wait:       time between the end of the initialisation and the reinitialisation in ms
        laser_re:   length of the reinitialisation pulse in ms

        apd_width:   time in ms where the apd is on (this needs to be varied to do the delay measurement)
        rest:       time in ms after the reinitialisation pulse to have a break before the next cycle starts
        start:      where in time should the apd pulse start shifting
        steps:      step size of the apd pulse sweep
        mw_wait_time: when should the mw pulse should happen?
        mw_pulse_time: length of the mw pulse in between the two laser pulses
        trigger:    True if the software trigger should be used and False if the loop array should define the length
        rep:        If there is no software trigger, rep gives the number of repetitions the waveform is played
        segment_map:If the order of different segments should be mixed (not there yet)
        '''
        # building the array to change the length of the apd pulse
        msplay = laser_in + wait + laser_re + rest
        clk_rate = MEGA(clk_rate)
        stop = (msplay - rest) - apd_width
        # print('stop value: ', stop)

        # build array with all the desired apd_start values]
        apd_start = []
        step_count = math.floor((stop - start) / steps)
        # print('step_count: ', step_count)
        for i in range(step_count):
            apd_start.append((start + (i * steps)))
        # print('adpstart array: ', apd_start)
        laser_in *= 1e-3
        mw_wait_time *= 1e-3
        mw_pulse_time *= 1e-3
        print('mw specs: ', mw_wait_time)
        print('mw specs: ', mw_pulse_time)
        laser_re *= 1e-3
        wait *= 1e-3
        apd_width *= 1e-3
        rest *= 1e-3
        start *= 1e-3
        steps *= 1e-3
        msplay *= 1e-3
        apd_start = np.multiply(1e-3, apd_start)
        # print('apd_start: ', apd_start)

        # Laser waveform does not change the whole time
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
            apd_on = np.ones(self._pulser.waveform_padding(clk_rate * apd_width))
            apd_off = np.zeros(self._pulser.waveform_padding(clk_rate * (msplay - (i + apd_width))))
            apd_do_waveform = np.concatenate((apd_wait_time, apd_on, apd_off))

            if mw_pulse == True:
                # if a mw pulse is desired:
                mw_wait = np.zeros(self._pulser.waveform_padding(clk_rate * mw_wait_time))
                pi_pulse_array = np.ones(self._pulser.waveform_padding(clk_rate * mw_pulse_time))
                mw_off = np.zeros(self._pulser.waveform_padding(clk_rate * msplay - (mw_wait_time + mw_pulse_time)))
                mw_ao_waveform = np.concatenate((mw_wait, pi_pulse_array, mw_off))
                mw_ao_waveform = np.subtract(mw_ao_waveform, 0.5)  # because the MW takes only +-0.5V
                print('mw_ao_waveform max: ', np.max(mw_ao_waveform))
                # print('len of mw_ao_waveform: ', len(mw_ao_waveform))

                # Make sure the arrays have the same length:
                # Check which array is the longest and append the two others:
                max_val = np.max(np.array([len(laser_do_waveform), len(apd_do_waveform), len(mw_ao_waveform)]))
                if max_val == 0:
                    print('laser waveforms is the longest')
                    difference_mw = len(laser_do_waveform) - len(mw_ao_waveform)
                    difference_apd = len(laser_do_waveform) - len(apd_do_waveform)
                    rest_array_mw = np.zeros(int(difference_mw))
                    rest_array_apd = np.zeros(int(difference_apd))
                    mw_ao_waveform = np.concatenate((mw_ao_waveform, rest_array_mw))
                    apd_do_waveform = np.concatenate((apd_do_waveform, rest_array_apd))
                    print('len of mw_ao_waveform: ', len(mw_ao_waveform))
                    print('len of laser_do_waveform: ', len(laser_do_waveform))
                    print('len of apd_do_waveform: ', len(apd_do_waveform))
                elif max_val == 1:
                    print('apd waveforms is the longest')
                    difference_mw = len(apd_do_waveform) - len(mw_ao_waveform)
                    difference_laser = len(apd_do_waveform) - len(laser_do_waveform)
                    rest_array_mw = np.zeros(int(difference_mw))
                    rest_array_laser = np.zeros(int(difference_laser))
                    mw_ao_waveform = np.concatenate((mw_ao_waveform, rest_array_mw))
                    laser_do_waveform = np.concatenate((laser_do_waveform, rest_array_laser))
                    print('len of mw_ao_waveform: ', len(mw_ao_waveform))
                    print('len of laser_do_waveform: ', len(laser_do_waveform))
                    print('len of apd_do_waveform: ', len(apd_do_waveform))
                else:
                    print('mw waveforms is the longest')
                    difference_apd = len(mw_ao_waveform) - len(apd_do_waveform)
                    difference_laser = len(mw_ao_waveform) - len(laser_do_waveform)
                    rest_array_apd = np.zeros(int(difference_apd))
                    rest_array_laser = np.zeros(int(difference_laser))
                    apd_do_waveform = np.concatenate((apd_do_waveform, rest_array_apd))
                    laser_do_waveform = np.concatenate((laser_do_waveform, rest_array_laser))
                    print('len of mw_ao_waveform: ', len(mw_ao_waveform))
                    print('len of laser_do_waveform: ', len(laser_do_waveform))
                    print('len of apd_do_waveform: ', len(apd_do_waveform))

                ao_dict = {2: mw_ao_waveform}

            else:  # if no MW is needed
                # check which array is longer from laser and apd and use that one
                if len(apd_do_waveform) < len(laser_do_waveform):
                    print('laser waveform is the longer')
                    difference = len(laser_do_waveform) - len(apd_do_waveform)
                    rest_array = np.zeros(int(difference))
                    apd_do_waveform = np.concatenate((apd_do_waveform, rest_array))
                    # Create zero array of same size:
                    d_off = np.zeros(len(laser_do_waveform))
                    print('len of mw_ao_waveform: ', len(d_off))
                    print('len of laser_do_waveform: ', len(laser_do_waveform))
                    print('len of apd_do_waveform: ', len(apd_do_waveform))
                elif len(laser_do_waveform) < len(apd_do_waveform):
                    print('apd waveform is the longer')
                    difference = len(apd_do_waveform) - len(laser_do_waveform)
                    rest_array = np.zeros(int(difference))
                    laser_do_waveform = np.concatenate((laser_do_waveform, rest_array))
                    # Create zero array of same size:
                    d_off = np.zeros(len(apd_do_waveform))
                    print('len of mw_ao_waveform: ', len(d_off))
                    print('len of laser_do_waveform: ', len(laser_do_waveform))
                    print('len of apd_do_waveform: ', len(apd_do_waveform))
                else:
                    print('apd waveform and laser_waveform have the same size')
                    d_off = np.zeros(len(apd_do_waveform))
                    print('len of mw_ao_waveform: ', len(d_off))
                    print('len of laser_do_waveform: ', len(laser_do_waveform))
                    print('len of apd_do_waveform: ', len(apd_do_waveform))
                # Creates a dict for the analog channel if no MW is needed
                ao_dict = {2: d_off}

            # Now all the waveformes are prepared and have the same length
            do_dict = {laser_outchan: laser_do_waveform, apd_outchan: apd_do_waveform}
            dosequence.append((do_dict))
            aosequence.append((ao_dict))

        # print('dosequence: ', dosequence)
        # print('aosequence: ', aosequence)
        # load and play
        # print('len of dosequence: ', len(dosequence))
        # print('len of aosequence: ', len(aosequence))
        digital_output_map = {
            0: [0, 1]}  # directs to the right output channel analog channel0 directs to digital x0 and x1
        for aos, dos in zip(aosequence, dosequence):
            # Loads waveforms to the awg
            self._pulser.load_waveform(ao_waveform_dictionary=aos, do_waveform_dictionary=dos,
                                       digital_output_map=digital_output_map)

        if trigger == True:
            # This sequence waits for a software trigger to start playing and moving to the next step.
            self._pulser.configure_ORmask(card_idx, None)
            self._pulser.configure_ANDmask(card_idx, None)
            array = np.array([0x40000000], dtype=np.int64)
            stop_condition_list = np.repeat(array, step_count)
            # print('stop_condition_list: ', stop_condition_list)
            loops = np.ones((len(aosequence)), dtype=np.int64)

        else:
            # This sequence immediately starts after the sequences are loaded
            self._pulser.configure_ORmask(card_idx, 'immediate')
            self._pulser.configure_ANDmask(card_idx, None)
            loop_array = np.array([rep], dtype=np.int64)
            loops = np.repeat(loop_array, step_count)
            print('loops: ', loops)
            stop_condition_list = np.array([])

        # makes the waveforms repeat in loops
        self._pulser.load_sequence(  # digital_sequences=dosequence,
            loops_list=loops, segment_map=segment_map, stop_condition_list=stop_condition_list)
        self._pulser.start_card(card_idx)
        self._pulser.arm_trigger(card_idx)

    def mw_test(self, msecondsplay, loops, segment_map=[]):
        """
        Function for early debugging of the awg. Remove from the final class.

        The steps to a successful initialization of a sequence are:
        - Define the clock rate
        - create the arrays that will represent the sequences
        - configure the trigger masks
        - load the sequence
        - start the card
        - arm the trigger
        """
        msecondsplay *= 1e-3

        clk_rate = MEGA(100)

        card_idx = 1
        self.set_sample_rate(card_idx, clk_rate)

        samples = self.waveform_padding(msecondsplay * clk_rate)
        time_ax = np.linspace(0, samples / clk_rate, samples)
        first_seq_ch0 = np.sin(2 * np.pi * time_ax / msecondsplay)
        second_seq_ch0 = np.sin(3 * 2 * np.pi * time_ax / msecondsplay)
        first_seq_ch1 = np.linspace(0, 1, len(time_ax))
        # func_1= 0.5 * (signal.square(2 * np.pi * 20 * time_ax)) + 0.5
        #
        sigma = (time_ax[-1]) / 20
        second_seq_ch1 = np.exp(-(time_ax - time_ax[-1] / 2) ** 2 / (2 * sigma ** 2)) * np.sin(
            20 * 2 * np.pi * time_ax / msecondsplay)

        aosequence = [{2: first_seq_ch0, 3: first_seq_ch1},
                      {2: first_seq_ch1, 3: first_seq_ch0},
                      {2: first_seq_ch0, 3: first_seq_ch1},
                      {2: first_seq_ch1, 3: first_seq_ch0}
                      ]

        do_chan = 1
        do1 = np.zeros(first_seq_ch0.shape)
        do2 = np.zeros(second_seq_ch0.shape)
        do3 = np.zeros(first_seq_ch1.shape)
        do4 = np.copy(do3)
        do1[first_seq_ch0 > 0] = 1
        do2[second_seq_ch0 > 0] = 1
        do3[:] = 1
        print('len do1', len(do1))
        print('len do2', len(do2))
        print('len do3', len(do3))
        print('len do4', len(do4))
        outchan = 0
        dosequence = [{outchan: do1},
                      {outchan: do2},
                      {outchan: do3},
                      {outchan: do4}]

        digital_output_map = {1: [outchan]}
        for aos, dos in zip(aosequence, dosequence):
            self.load_waveform(ao_waveform_dictionary=aos, do_waveform_dictionary=dos,
                               digital_output_map=digital_output_map)

        # This sequence immediately starts after the sequences are loaded
        # self.configure_ORmask(card_idx, 'immediate')
        # self.configure_ANDmask(card_idx, None)
        # stop_condition_list = np.array([])

        # This sequence waits for a software trigger to start playing and moving to the next step.
        self.configure_ORmask(card_idx, None)
        self.configure_ANDmask(card_idx, None)
        loops = np.ones((len(aosequence)), dtype=np.int64)
        # print('loops: ', loops)
        stop_condition_list = np.array(
            [SPCSEQ_ENDLOOPONTRIG, SPCSEQ_ENDLOOPONTRIG, SPCSEQ_ENDLOOPONTRIG, SPCSEQ_ENDLOOPONTRIG], dtype=np.int64)
        # print('segment map: ', segment_map)
        # print('stop cond.:', stop_condition_list)

        self.load_sequence(  # digital_sequences=dosequence,
            loops_list=loops, segment_map=segment_map, stop_condition_list=stop_condition_list)

        self.start_card(card_idx)
        self.arm_trigger(card_idx)
