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


class Simplepulse(GenericLogic):

    pulsegenerator = Connector(interface='PulserInterface')


    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)




    def on_activate(self):
    #     """ Initialisation performed during activation of the module.
    #     """

        # Get connectors
        self._pulser = self.pulsegenerator()
        active_channels= self._pulser.get_active_channels()
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

    def set_up_pulser(self,card_idx,channels, clk_rate_raw, filter_active):
    # def set_up_pulser(self):
        clk_rate= 10 * 1000 * clk_rate_raw
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
        output_sample_rate= self._pulser.get_sample_rate(card_idx)
        return output_sample_rate

### For digital signals ###
    def create_d_on(self,samples, clk_rate):
        '''
        Creates an array of ones for a digital channel
        len_ms:         is the desired duration in ms
        clk_rate:       samples/second of the AWG
        '''
        clk_rate = 10 * clk_rate
        # samples = self._pulser.waveform_padding(len_ms * clk_rate)
        d_on=np.ones(samples)
        # print('d_on=', d_on)
        # print('len of d_on=', len(d_on))
        return d_on
    def create_d_off(self,samples, clk_rate):
        '''
        Creates an array of zeros for a digital channel
        len_ms:         is the desired duration in ms
        clk_rate:       samples/second of the AWG
        '''
        clk_rate = 10 * clk_rate
        # samples = self._pulser.waveform_padding(len_ms * clk_rate)
        d_off = np.zeros(samples)
        # print('d_off=', d_off)
        # print('len of d_off=', len(d_off))
        return d_off


### For analog signals ###
    def sine_wave(self, msplay, clk_rate):
        samples = self.waveform_padding(msplay * clk_rate)
        time_ax = np.linspace(0, samples / clk_rate, samples)
        waveform = np.sin(2 * np.pi * time_ax / msplay)
        return waveform
#
# # Load analog waveform to the awg and play it
#     def load_analog_waveform(self,msplay,clk_rate):
#     self.sine_wave()
#     aosequence = [{2: first_seq_ch0, 3: first_seq_ch1},
#                   {2: first_seq_ch1, 3: first_seq_ch0},
#                   {2: first_seq_ch0, 3: first_seq_ch1},
#                   {2: first_seq_ch1, 3: first_seq_ch0}
#                   ]



    def build_block_element(self, clk_rate, len0, len1):
        '''
        Creates an array of zeros and ones for a digital channel --> its just one step
        len0:           is the desired duration of zeros before the pulse in ms
        len1:           is the desired duration of ones in ms
        clk_rate:       samples/second of the AWG
        Returns the block_element as an array and the length of the array in units of samples
        '''
        clk_rate = 10 * clk_rate
        block_element=[]
        low = self.create_d_off(len0,clk_rate)
        high = self.create_d_on(len1,clk_rate)
        block_element = np.concatenate((low,high))
        block_element_len= len(block_element)
        # print(len(block_element)
        return block_element, block_element_len


    def build_sequence(self, clk_rate, len0, len1, msplay):
        '''
        Creates an array of block elements with the length of msplay for a digital channel
        len0:          is the desired duration of zeros before the pulse in ms
        len1:           is the desired duration of ones in ms
        clk_rate:       samples/second of the AWG
        msplay:          is the length of the whole sequence
        '''
        clk_rate = 10 * clk_rate
        samples = self._pulser.waveform_padding(msplay * clk_rate)
        print('totalsamples: ', samples)

        sequence=[]
        block_element = self.build_block_element(clk_rate,len0,len1)[0]
        print('len of block element: ',self.build_block_element(clk_rate,len0,len1)[1])
        rep = round(np.divide(samples, self.build_block_element(clk_rate,len0,len1)[1]))

        print('Reps: ',rep)
        sequence = np.repeat(block_element, rep)
        # sequence_len = len(sequence)
        # return sequence


    def build_waveform(self, msplay, clk_rate, loops, card_idx, segment_map=[]):  #This works
        '''
        Makes the AWG play a analog and digital waveform
        clk_rate:       samples/second of the AWG
        loops:            number of repetitions of each block
        '''
        msplay *= 1e-3
        clk_rate = MEGA(clk_rate)
        #Makes sure samples has the right length
        samples = self._pulser.waveform_padding(msplay * clk_rate)
        # #Creates a time axis of the same length
        time_ax = np.linspace(0, samples / clk_rate, samples)
        #Just a sine with amplitude 1
        first_seq_ch0 = np.sin(2 * np.pi * time_ax / msplay)
        # linear function
        first_seq_ch1 = np.linspace(0, 1, len(time_ax))
        #Building a dict for anaolog signal which can be read by the AWG
        aosequence = [{2: first_seq_ch0, 3: first_seq_ch1}]

        # for digital signal
        do_chan = 1
        do1 = np.zeros(first_seq_ch0.shape) #does zeros array with the length of the sine function which has the same length as samples
        do1[first_seq_ch0 > 0] = 1   #Creates one step function -_ first ones then zeros

        outchan = 0
        # Building a dict for anaolog signal which can be read by the AWG
        dosequence = [{outchan: do1}]

        digital_output_map = {1: [outchan]} #directs to the right output channel
        for aos, dos in zip(aosequence, dosequence):
            #Loads waveforms to the awg
            self._pulser.load_waveform(ao_waveform_dictionary=aos, do_waveform_dictionary=dos,
                               digital_output_map=digital_output_map)

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
        #makes the waveforms repeat in loops
        self._pulser.load_sequence(  # digital_sequences=dosequence,
            loops_list=loops, segment_map=segment_map, stop_condition_list=stop_condition_list)

        self._pulser.start_card(card_idx)
        self._pulser.arm_trigger(card_idx)
        self._pulser.clear_all()

    def build_ao_signal(self,clk_rate, msplay):
        #this one should hava a dict as output
        msplay *= 1e-3
        clk_rate = MEGA(clk_rate)
        # Makes sure samples has the right length
        samples = self._pulser.waveform_padding(msplay * clk_rate)
        # #Creates a time axis of the same length
        time_ax = np.linspace(0, samples / clk_rate, samples)
        # Just a sine with amplitude 1
        sine = np.sin(2 * np.pi * time_ax / msplay)
        # linear function
        zeros= np.zeros(len(time_ax))
        # Building a dict for anaolog signal which can be read by the AWG
        aosequence = [{2: sine, 3: zeros}]
        return aosequence

    def build_do_signal(self,clk_rate, msplay):
        #this one should have a dict as output
        # this one should hava a dict as output
        msplay *= 1e-3
        clk_rate = MEGA(clk_rate)
        # Makes sure samples has the right length
        samples = self._pulser.waveform_padding(msplay * clk_rate)
        # #Creates a time axis of the same length
        time_ax = np.linspace(0, samples / clk_rate, samples)
        # Just a sine with amplitude 1
        sine = np.sin(2 * np.pi * time_ax / msplay)
        do1 = np.zeros(sine.shape)  # does zeros array with the length of the sine function which has the same length as samples
        do1[sine > 0] = 1  # Creates one step function -_ first ones then zeros
        outchan = 0
        # Building a dict for anaolog signal which can be read by the AWG
        dosequence = [{outchan: do1}]
        return dosequence

    def play_sequence(self, clk_rate, msplay, d_outchan, card_idx, segment_map=[]):
        # This one takes an ao input and do input as two dicts and makes them play
        #Get aosequence and dosequence
        aosequence = self.build_ao_signal(clk_rate,msplay)
        dosequence = self.build_do_signal(clk_rate,msplay)
        digital_output_map = {1: [d_outchan]}  # directs to the right output channel
        for aos, dos in zip(aosequence, dosequence):
            # Loads waveforms to the awg
            self._pulser.load_waveform(ao_waveform_dictionary=aos, do_waveform_dictionary=dos,
                                       digital_output_map=digital_output_map)

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
            loops_list=loops, segment_map=segment_map, stop_condition_list=stop_condition_list)

        self._pulser.start_card(card_idx)
        self._pulser.arm_trigger(card_idx)
        self._pulser.clear_all()



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
    #     self._pulser.test_seq(msecondsplay, loop_array, segment_map)