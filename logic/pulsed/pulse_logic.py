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


import numpy as np
import math
from qtpy import QtCore
from hardware.awg.spectrum_awg.py_header.regs import *
from core.connector import Connector
from logic.generic_logic import GenericLogic
from core.statusvariable import StatusVar


# Just hardcode the channels to the different hardware modules somewhere
# Microwave: analog channel 2 and 3 for I and Q parameter. For now we just need I
# Laser: digital channel X0
# Photon counter: digital channel X1 and X2 (for reference)


class Pulse(GenericLogic):
    pulsegenerator = Connector(interface="PulserInterface")
    if_modulation_freq = StatusVar(default=100e6)
    # Use this signal to transfer the created Arrays to the masterpulse logic
    # sigUsedArrays = QtCore.Signal(dict, dict)

    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)

    def on_activate(self):
        # Initialisation performed during activation of the module

        # Get connectors
        self._pulser = self.pulsegenerator()
        active_channels = self._pulser.get_active_channels()
        self._pulser.clear_all()
        self._pulser.reset()
        self._if_freq = self.if_modulation_freq
        self._if_freq_mhz = self._if_freq / 1e6 # if_req in MHz
        return active_channels

    def on_deactivate(self):
        self._pulser.clear_all()
        self._pulser.reset()
        self._pulser.pulser_off()

    ###############     Gerneral Methods    ###############
    def set_up_pulser(self, clk_rate=1250):
        card_idx = 1
        clk_rate = MEGA(clk_rate)
        self._pulser.pulser_off()
        self._pulser.clear_all()
        self._pulser.set_sample_rate(card_idx, int(clk_rate))
        self.sample_rate = self.get_clk_rate()


    def get_clk_rate(self):
        self.sample_rate = self._pulser.get_sample_rate(card_idx=1)
        return self.sample_rate

    def trigger(self):
        self._pulser.send_software_trig(1)

    def stop_awg(self):  # Makes everything stop
        self._pulser.clear_all()
        self._pulser.stop_replay(1)

    def play_any(
        self,
        clk_rate,
        seq_len,
        laser_times,
        apd_times=[],
        apd_ref=[],
        mw_times=[],
        freq=2.87e9,
        method=None,
        rep=100,
        mw_pulse=False,
        trigger=True,
        neg_mw=False,
        seq_map_req=True,
        phase_shift=False
    ):
        """
        clk_rate =  clock rate of the AWG in microseconds
        seq_len =   whole length of the sequence in microseconds
        laser_times = [laser_in, wait, laser_re] all in microseconds
                laser_in:   length of initialisation laser pulse (it starts right at the beginning) in microseconds
                wait:       time between the end of the initialisation and the reinitialisation in microseconds
                laser_re:   length of the reinitialisation pulse in microseconds
        mw_sweep = [mw_wait_time, mw_pulse_time] all in microseconds
                mw_wait_time: time between the end of the first lase r pulse and the beginning of the microwave (delay)
                mw_pulse_time: length of the mw pulse in between the two laser pulses
        apd_times = [steps, apd_width]
                steps:      step size of the apd pulse sweep
                apd_width:   time in microseconds where the apd is on (length of the pulse)
        method =    can be rabi, ramsey or delaysweep
        rep =       repetition time when no trigger is used
        mw_pulse =  With or without MW pulse in the middle
        trigger =   either True: software trigger or False: no trigger at all
        """
        self.set_up_pulser(clk_rate)
        if method == "delaysweep":
            # print("It does delay sweep.")
            self.delay_sweep(
                seq_len,
                laser_times,
                apd_times,
                apd_ref,
                mw_times,
                rep,
                trigger=trigger,
                mw_pulse=mw_pulse,
            )

        elif method == "delaysweep_ref":
            # do the delaysweep with the reference counter
            self.delay_sweep_ref(
                seq_len,
                laser_times,
                apd_times,
                apd_ref,
                mw_times,
                rep,
                trigger=trigger,
                mw_pulse=mw_pulse,
            )

        elif method == "ramsey":
            # Do ramsey
            # print("It does Ramsey.")
            self.play_ramsey(
                seq_len, laser_times, apd_times, apd_ref, mw_times, freq, seq_map_req, rep,
                trigger, phase_shift
            )

        elif method == "rabi":
            # Do Rabi
            # print("It does Rabi.")
            self.play_rabi(
                seq_len, laser_times, apd_times, apd_ref, mw_times, freq, rep, seq_map_req, trigger)

        else:
            self.log.warning(
                "Method can only be: rabi, ramsey, delaysweep or delaysweep_ref"
            )
        return method, laser_times, apd_times, apd_ref, mw_times
    def get_step_count(self, array):
        step_count = int(((array[2] - array[1]) / array[3]) + 1)
        return step_count

    ###############     Laser and APD blocks    ###############

    # def laser_seq(self, seq_len, laser_times):
    #     """
    #     seq_len =   whole length of the sequence in microseconds
    #     laser_times = [laser_in, wait, laser_re] all in microseconds
    #             laser_in:   length of initialisation laser pulse (it starts right at the beginning) in microseconds
    #             wait:       time between the end of the initialisation and the reinitialisation in microseconds
    #             laser_re:   length of the reinitialisation pulse in microseconds
    #     """
    #     seq_len *= 1e-6  # in seconds
    #     laser_times = np.multiply(laser_times, 1e-6)  # in seconds
    #     card_idx = 1
    #     # takes sampling rate which is there already
    #     clk_rate = self._pulser.get_sample_rate(card_idx)
    #     laser_len = (
    #         laser_times[0] + laser_times[1] + laser_times[2]
    #     )  # len of the important part in s
    #     len_rest = seq_len - laser_len  # len of the rest in s
    #     if seq_len < laser_len:
    #         self.log.warning(
    #             "The total length needs to be larger that the laser_time sum"
    #         )
    #     else:
    #         # Laser waveform does not change the whole time
    #         first_pulse = np.ones(int(round(clk_rate * laser_times[0])))
    #         laser_wait_time = np.zeros(int(round(clk_rate * laser_times[1])))
    #         sec_pulse = np.ones(int(round(clk_rate * laser_times[2])))
    #         rest_array = np.zeros(int(round(clk_rate * len_rest)))
    #         laser_do_waveform = np.concatenate(
    #             (first_pulse, laser_wait_time, sec_pulse, rest_array)
    #         )
    #         len_laser_do_waveform = self._pulser.waveform_padding(
    #             len(laser_do_waveform)
    #         )
    #         difference = len_laser_do_waveform - (clk_rate * seq_len)
    #         laser_do_waveform = np.concatenate(
    #             (laser_do_waveform, np.zeros(int(difference)))
    #         )
    #
    #     return laser_do_waveform

    # def apd_seq(self, seq_len, apd_times):
    #     """
    #     len =   whole length of the sequence microseconds
    #             apd_times = [time to start, length] all in microseconds
    #     """
    #     seq_len *= 1e-6
    #     apd_times = np.multiply(apd_times, 1e-6)
    #     card_idx = 1
    #     # takes sampling rate which is there already
    #     clk_rate = self._pulser.get_sample_rate(card_idx)
    #     apd_len = apd_times[0] + apd_times[1]
    #     len_rest = seq_len - apd_len  # len of the rest in s
    #     if seq_len < apd_len:
    #         raise ValueError(
    #             "The total length needs to be larger that the apd_time sum"
    #         )
    #     else:
    #         # apd waveform
    #         apd_wait = np.zeros(int(round(clk_rate * apd_times[0])))
    #         apd_on = np.ones(int(round(clk_rate * apd_times[1])))
    #         rest_array = np.zeros(int(round(clk_rate * (seq_len - apd_len))))
    #         apd_do_waveform = np.concatenate((apd_wait, apd_on, rest_array))
    #         len_apd_do_waveform = self._pulser.waveform_padding(len(apd_do_waveform))
    #         difference = len_apd_do_waveform - (clk_rate * seq_len)
    #         apd_do_waveform = np.concatenate(
    #             (apd_do_waveform, np.zeros(int(difference)))
    #         )
    #
    #     return apd_do_waveform

    # def just_zeros(self, seq_len):
    #     """
    #     len =   whole length of the sequence in microseconds
    #     apd_times = [time to start, length] all in microseconds
    #     """
    #     seq_len *= 1e-6
    #     analog_channel = 1  # Laser
    #     card_idx = 1
    #     # takes sampling rate which is there already
    #     clk_rate = self._pulser.get_sample_rate(card_idx)
    #     array = np.zeros(int(round(clk_rate * seq_len)))
    #     len_waveform = self._pulser.waveform_padding(len(array))
    #     difference = len_waveform - (clk_rate * seq_len)
    #     array_zeros = np.concatenate((array, np.zeros(int(difference))))
    # #     return array_zeros
    #
    # def just_ones(self, seq_len):
    #     """
    #     len =   whole length of the sequence in microseconds
    #     apd_times = [time to start, length] all in microseconds
    #     """
    #     seq_len *= 1e-6
    #     analog_channel = 1  # Laser
    #     card_idx = 1
    #     # takes sampling rate which is there already
    #     clk_rate = self._pulser.get_sample_rate(card_idx)
    #     array = np.ones(int(round(clk_rate * seq_len)))
    #     len_waveform = self._pulser.waveform_padding(len(array))
    #     difference = len_waveform - (clk_rate * seq_len)
    #     array_ones = np.concatenate((array, np.zeros(int(difference))))
    #     return array_ones

    # def apd_ref(self, seq_len, apd_ref):
    #     """
    #     len =   whole length of the sequence microseconds
    #             apd_times = [time to start, length] all in microseconds
    #     """
    #     seq_len *= 1e-6
    #     apd_ref = np.multiply(apd_ref, 1e-6)
    #     card_idx = 1
    #     # takes sampling rate which is there already
    #     clk_rate = self._pulser.get_sample_rate(card_idx)
    #     apd_len = apd_ref[0] + apd_ref[1]
    #     if seq_len < apd_len:
    #         raise ValueError(
    #             "The total length needs to be larger that the apd_time sum"
    #         )
    #     else:
    #         # apd waveform
    #         apd_wait = np.zeros(int(round(clk_rate * apd_ref[0])))
    #         apd_on = np.ones(int(round(clk_rate * apd_ref[1])))
    #         rest_array = np.zeros(int(round(clk_rate * (seq_len - apd_len))))
    #         apd_do_waveform = np.concatenate((apd_wait, apd_on, rest_array))
    #         len_apd_do_waveform = self._pulser.waveform_padding(len(apd_do_waveform))
    #         difference = len_apd_do_waveform - (clk_rate * seq_len)
    #         apd_ref_do_waveform = np.concatenate(
    #             (apd_do_waveform, np.zeros(int(difference)))
    #         )
    #     return apd_ref_do_waveform

    ###############     Different seqeuences    ###############

    def play_rabi(
        self,
        seq_len = 5,
        laser_times = [1, 1, 2],
        apd_times = [2.05, 0.25],
        apd_ref = [3.05, 0.25],
        mw_times = [1.05, 0, 0.9, 0.1],
        freq = 2.87e9,
        rep=100,
        seq_map_req=True,
        trigger=True):
        """
        all are in microseconds
        laser_times = [laser_in, wait, laser_re] all in microseconds
                laser_in:   length of initialisation laser pulse (it starts right at the beginning) in microseconds
                wait:       time between the end of the initialisation and the reinitialisation in microseconds
                laser_re:   length of the reinitialisation pulse in microseconds
        seq_len =   whole length of the sequence in microseconds
        apd_times = [time to start, length] all in microseconds
        apd_ref = [time to start, length] all in microseconds
        mw_times = [mw_start_time, mw_len_0, mw_len_max, steps] all in microseconds
                mw_start_time:  time from the beginning when MW pulse should start, does not change
                mw_len_0:       minimal length of the mw pulse
                mw_len_max:     maximal length of the mw pulse
                steps:          stepsize for increasing the mw pulse duration
        This works fine:        pulselogic.play_rabi(seq_len=25, laser_times=[2,2,6], apd_times=[4.5,1], mw_times=[2,0.5,2,0.5], rep=100, trigger=True)
        """

        mw_len = []


        card_idx = 1

        self.clk_rate = self._pulser.get_sample_rate(card_idx) / 1e6 # To adapt it to the microsecond timing
        step_count = int(((mw_times[2] - mw_times[1]) / mw_times[3]) + 1)

        for i in range(step_count):
            mw_len.append(round((mw_times[1] + (i * mw_times[3])), 3))  # starting with minimal length

        t_centrewidth_list_mw = [] # get the size of the np.array right?np.array([[0, 0]])
        t_centrewidth_list_laser = self.convert_laser_val(laser_times)
        t_centrewidth_list_apd = self.convert_apd_val(apd_times)
        t_centrewidth_list_apd_ref = self.convert_apd_val(apd_ref)

        for i in mw_len:
            if seq_len < (mw_times[0] + mw_times[2]):
                raise ValueError(
                    "The total length needs to be larger that the apd_time sum"
                )
            else:
                t0, width = self.convert_mw_val_rabi(mw_times, i) # in microseconds
                t_centrewidth_list_mw.append([np.array([[t0, width]])])

        time_ax = self.time_axis(seq_len) # this includes waveform_padding
        self._pulser.set_frequency(freq)
        output_frequency = self._pulser.get_frequency()

        for step in range(len(t_centrewidth_list_mw)):
            iwaveform, qwwaveform = self._pulser.iq_pulses(time_ax, t_centrewidth_list_mw[step], output_frequency)
            laser_waveform = self._pulser.box_envelope(time_ax, t_centrewidth_list_laser)
            apd_waveform = self._pulser.box_envelope(time_ax, t_centrewidth_list_apd)  #
            apd_ref_waveform = self._pulser.box_envelope(time_ax, t_centrewidth_list_apd_ref)  #
            analogs = {"i_chan": iwaveform, "q_chan": qwwaveform}
            digitals = {"laser": laser_waveform, "apd_sig": apd_waveform, "apd_read": apd_ref_waveform} # read/sig?

            self._pulser.load_waveform(
                iq_dictionary=analogs,
                digital_pulses=digitals,
                digital_output_map={0: [0, 1, 2]}
            )

        if trigger:
            # This sequence waits for a software trigger to start playing and moving to the next step.
            self._pulser.configure_ORmask(card_idx, None)
            self._pulser.configure_ANDmask(card_idx, None)
            array = np.array([0x40000000], dtype=np.int64)
            stop_condition_list = np.repeat(array, step_count)
            loops = np.ones((len(mw_len)), dtype=np.int64)

        else:
            # This sequence immediately starts after the sequences are loaded
            self._pulser.configure_ORmask(card_idx, "immediate")
            self._pulser.configure_ANDmask(card_idx, None)
            loop_array = np.array([rep], dtype=np.int64)
            loops = np.repeat(loop_array, step_count)
            stop_condition_list = np.array([])
            #only build a seq_map when its requested, if not dont do anything
        if seq_map_req:
            self.seq_map = self.build_seq_map(step_count)
        else:
            self.seq_map = []

        self._pulser.load_sequence(
            loops_list=loops,
            segment_map=self.seq_map,
            stop_condition_list=stop_condition_list,)


        self._pulser.start_card(card_idx)
        self._pulser.arm_trigger(card_idx)
        # self.sigUsedArrays.emit(analogs, digitals)
        return self.seq_map

    def play_ramsey(self,
                    seq_len = 5,
                    laser_times = [1, 1, 2],
                    apd_times = [3.25, 0.25],
                    apd_ref = [4.25, 0.25],
                    mw_times = [1.1, 0, 0.8, 0.1, 0.1],
                    freq=2.87e9,
                    seq_map_req=True,
                    trigger=True,
                    rep=100,
                    phase_shift=True):
        """
        The arrays are given in microseconds
        It changes the distance between two pi/2 pusles, the position of the first pulse stays the same
        seq_len =   whole length of the sequence in microseconds
        laser_times = [laser_in, wait, laser_re] all in microseconds
                laser_in:   length of initialisation laser pulse (it starts right at the beginning) in microseconds
                wait:       time between the end of the initialisation and the reinitialisation in microseconds
                laser_re:   length of the reinitialisation pulse in microseconds
        apd_times = [time to start, length of the pulse] both in microseconds
        apd_ref = [time to start, length of the pulse] both in microseconds
        mw_times = [start, distance_min, distance_max, steps, pulse length] all in microseconds
                start:          start position of the fist microwave pulse
                duration_min:   minimal distance between the two pulses
                duration_max:   maximal distance between the two pulses
                steps:          step size for increasing the mw pulse duration
                pulse length:   duration of pulses in microseconds (pi/2)
        Amplitude of the MW pulse should be +- 0.5V
        rep =       repetition time when no trigger is used
        trigger =   either True: software trigger or False: no trigger at all
        """

        t_centrewidth_list_mw = []

        break_len = []
        card_idx = 1

        step_count = int(((mw_times[2] - mw_times[1]) / mw_times[3]) + 1) # without unit

        for i in range(step_count):
            break_len.append((mw_times[1] + (i * mw_times[3])))  # starting with minimal length
        t_centrewidth_list_laser = self.convert_laser_val(laser_times)
        t_centrewidth_list_apd = self.convert_apd_val(apd_times)
        t_centrewidth_list_apd_ref = self.convert_apd_val(apd_ref)

        if seq_len < (laser_times[0] + laser_times[1] + laser_times[2]):
            self.log.warning(
                "The total length needs to be larger that the sequence time"
            )
            return
        else:
            for i in break_len:
                t0, t1, width = self.convert_mw_val_ramsey(mw_times, i)  # these values are in microseconds
                if phase_shift:
                    mwpulse2append = [np.array([[t0, width]]), np.array([[t1, width]])]
                else:
                    mwpulse2append = [np.array([[t0, width], [t1, width]])]
                t_centrewidth_list_mw.append(mwpulse2append)

        time_ax = self.time_axis(seq_len) # in microseconds
        self._pulser.set_frequency(freq) # this one does freq - if_freq
        output_frequency = self._pulser.get_frequency()
        if phase_shift:
            phases = np.array([0, np.pi])
        else:
            phases = None

        for mw_sequence_step in t_centrewidth_list_mw:
            iwaveform, qwwaveform = self._pulser.iq_pulses(time_ax, mw_sequence_step, output_frequency, phases)
            laser_waveform = self._pulser.box_envelope(time_ax, t_centrewidth_list_laser)
            apd_waveform = self._pulser.box_envelope(time_ax, t_centrewidth_list_apd)  #
            apd_ref_waveform = self._pulser.box_envelope(time_ax, t_centrewidth_list_apd_ref)  #
            ## Add a minus sign here?
            analogs = {"i_chan": iwaveform, "q_chan": qwwaveform}
            digitals = {"laser": laser_waveform, "apd_sig": apd_waveform, "apd_read": apd_ref_waveform} # read/sig?
            self._pulser.load_waveform(
                iq_dictionary=analogs,
                digital_pulses=digitals,
                digital_output_map={0: [0, 1, 2]}
            )

        if trigger:
            # This sequence waits for a software trigger to start playing and moving to the next step.
            self._pulser.configure_ORmask(card_idx, None)
            self._pulser.configure_ANDmask(card_idx, None)
            array = np.array([0x40000000], dtype=np.int64)
            stop_condition_list = np.repeat(array, step_count)
            loops = np.ones((len(break_len)), dtype=np.int64)

        else:
            # This sequence immediately starts after the sequences are loaded
            self._pulser.configure_ORmask(card_idx, "immediate")
            self._pulser.configure_ANDmask(card_idx, None)
            loop_array = np.array([rep], dtype=np.int64)
            loops = np.repeat(loop_array, step_count)
            stop_condition_list = np.array([])

        if seq_map_req:
            self.seq_map = self.build_seq_map(step_count)
        else:
            self.seq_map = []

        self._pulser.load_sequence(
            loops_list=loops,
            segment_map=self.seq_map,
            stop_condition_list=stop_condition_list)
        self._pulser.start_card(card_idx)
        self._pulser.arm_trigger(card_idx)
        # self.sigUsedArrays.emit(analogs, digitals)
        # return self.seq_map

    def delay_sweep(
        self,
            seq_len=7,
            laser_times=[0, 1, 5],
            apd_times=[0.1, 1, 4, 0.1],
            apd_ref=[4.5, 0.1],
            mw_times=[0.2, 0.2],
            rep=100,
            trigger=True,
            mw_pulse=False,
    ):
        """
        seq_len =   whole length of the sequence in microseconds
        laser_times = [laser_in, wait, laser_re] all in microseconds
                laser_in:   length of initialisation laser pulse (it starts right at the beginning) in microseconds
                wait:       time between the end of the initialisation and the reinitialisation in microseconds
                laser_re:   length of the reinitialisation pulse in microseconds
        apd_times = [length, min_start, max_start, steps] # changing length
                length:      time to wait until the pulse should start
                min_start:    minimal length of the apd pulse
                max_start:    maximal length of the apd pulse
                steps:      step size of the apd pulse sweep

        apd_ref = [time to start, length of the pulse] both in microseconds
        mw_times = [mw_wait_time, mw_pulse_time] all in microseconds
                mw_wait_time: time between the end of the first laser pulse and the beginning of the microwave (delay)
                mw_pulse_time: length of the mw pulse in between the two laser pulses
        rep =       repetition time when no trigger is used
        mw_pulse =  With or without MW pulse in the middle
        trigger =   either True: software trigger or False: no trigger at all
        """
        apd_times_s = np.multiply(apd_times, 1e-6)


        # So far we just use card1
        card_idx = 1

        segment_map = []
        apd_start = []
        step_count = int(((apd_times[2] - apd_times[1]) / apd_times[3]) + 1)

        time_ax = self.time_axis(seq_len)  # this includes waveform_padding
        for i in range(step_count):
            apd_start.append((apd_times_s[1] + (i * apd_times_s[3])))

        # Create laser and apd waveform as they are constant (put it in microseconds)
        t_centrewidth_list_laser = self.convert_laser_val(laser_times)
        t_centrewidth_list_apd_ref = self.convert_apd_val(apd_ref)
        laser_waveform = self._pulser.box_envelope(time_ax, t_centrewidth_list_laser)
        apd_ref_waveform = self._pulser.box_envelope(time_ax, t_centrewidth_list_apd_ref)
        if mw_pulse:
            # replace this one
            t_centrewidth_list_mw = self.convert_mw_sweep(mw_times)
            mw_waveform = self._pulser.box_envelope(time_ax, t_centrewidth_list_mw)

        else:  # if no MW is needed
            zeros_mw = self.convert_mw_sweep([-1,0])
            mw_waveform = self._pulser.box_envelope(time_ax, zeros_mw)


        for i in apd_start:
            current_t_centerwidth = self.convert_apd_val([i * 1e6, apd_times[0]])
            apd_waveform = self._pulser.box_envelope(time_ax, current_t_centerwidth)
            analogs = {"i_chan": mw_waveform, "q_chan": np.zeros(len(mw_waveform))}
            digitals = {"laser": laser_waveform, "apd_sig": apd_waveform, "apd_read": apd_ref_waveform}
            self._pulser.load_waveform(
                iq_dictionary=analogs,
                digital_pulses=digitals,
                digital_output_map={0: [0, 1, 2]})

        if trigger:
            # This sequence waits for a software trigger to start playing and moving to the next step.
            self._pulser.configure_ORmask(card_idx, None)
            self._pulser.configure_ANDmask(card_idx, None)
            array = np.array([0x40000000], dtype=np.int64)
            stop_condition_list = np.repeat(array, step_count)
            loops = np.ones((len(apd_start)), dtype=np.int64)

        else:
            # This sequence immediately starts after the sequences are loaded
            self._pulser.configure_ORmask(card_idx, "immediate")
            self._pulser.configure_ANDmask(card_idx, None)
            loop_array = np.array([rep], dtype=np.int64)
            loops = np.repeat(loop_array, step_count)
            stop_condition_list = np.array([])

        self._pulser.load_sequence(
            loops_list=loops,
            segment_map=segment_map,
            stop_condition_list=stop_condition_list,
        )
        self._pulser.start_card(card_idx)
        self._pulser.arm_trigger(card_idx)

    def delay_sweep_ref(
        self,
            seq_len=7,
            laser_times=[0, 1, 5],
            apd_times=[0.1, 1, 4, 0.1],
            apd_ref=[4.5, 0.1],
            mw_times=[0.2, 0.2],
            rep=100,
            trigger=True,
            mw_pulse=False,
    ):
        """
        The same but asign them the other way around
        seq_len =   whole length of the sequence in microseconds
        laser_times = [laser_in, wait, laser_re] all in microseconds
                laser_in:   length of initialisation laser pulse (it starts right at the beginning) in microseconds
                wait:       time between the end of the initialisation and the reinitialisation in microseconds
                laser_re:   length of the reinitialisation pulse in microseconds
        apd_times = [time_to_wait, min_start, max_start, steps] this will be connected to the reference apd pulse but still doing the sweep
                length:      time to wait until the pulse should start
                min_start:    minimal length of the apd pulse
                max_start:    maximal length of the apd pulse
                steps:      step size of the apd pulse sweep
        apd_ref = [time to start, length of the pulse] both in microseconds this will be connected to the real apd put it somewhere and dont sweep it
        mw_times = [mw_wait_time, mw_pulse_time] all in microseconds
                mw_wait_time: time between the end of the first laser pulse and the beginning of the microwave (delay)
                mw_pulse_time: length of the mw pulse in between the two laser pulses
        rep =       repetition time when no trigger is used
        mw_pulse =  With or without MW pulse in the middle
        trigger =   either True: software trigger or False: no trigger at all
        """

        apd_times_s = np.multiply(apd_times, 1e-6)
        # Sofar we just use card1
        card_idx = 1

        segment_map = []
        time_ax = self.time_axis(seq_len)  # this includes waveform_padding
        apd_start = []
        step_count = int(((apd_times[2] - apd_times[1]) / apd_times[3]) + 1)

        for i in range(step_count):
            apd_start.append((apd_times_s[1] + (i * apd_times_s[3])))

        # Create laser and apd waveform as they are constant (put it in microseconds)
        t_centrewidth_list_laser = self.convert_laser_val(laser_times)
        t_centrewidth_list_apd_ref = self.convert_apd_val(apd_ref)
        laser_waveform = self._pulser.box_envelope(time_ax, t_centrewidth_list_laser)
        apd_ref_waveform = self._pulser.box_envelope(time_ax, t_centrewidth_list_apd_ref)
        if mw_pulse:
            t_centrewidth_list_mw = self.convert_mw_sweep(mw_times)
            mw_waveform = self._pulser.box_envelope(time_ax, t_centrewidth_list_mw)

        else:  # if no MW is needed
            zeros_mw = self.convert_mw_sweep([-1, 0])
            mw_waveform = self._pulser.box_envelope(time_ax, zeros_mw)

        for i in apd_start:
            current_t_centerwidth = self.convert_apd_val([i * 1e6, apd_times[0]])
            apd_waveform = self._pulser.box_envelope(time_ax, current_t_centerwidth)
            analogs = {"i_chan": mw_waveform, "q_chan": np.zeros(len(mw_waveform))}
            digitals = {"laser": laser_waveform, "apd_sig": apd_ref_waveform, "apd_read": apd_waveform}
            self._pulser.load_waveform(
                iq_dictionary=analogs,
                digital_pulses=digitals,
                digital_output_map={0: [0, 1, 2]})

        if trigger == True:
            # This sequence waits for a software trigger to start playing and moving to the next step.
            self._pulser.configure_ORmask(card_idx, None)
            self._pulser.configure_ANDmask(card_idx, None)
            array = np.array([0x40000000], dtype=np.int64)
            stop_condition_list = np.repeat(array, step_count)
            loops = np.ones((len(apd_start)), dtype=np.int64)

        else:
            # This sequence immediately starts after the sequences are loaded
            self._pulser.configure_ORmask(card_idx, "immediate")
            self._pulser.configure_ANDmask(card_idx, None)
            loop_array = np.array([rep], dtype=np.int64)
            loops = np.repeat(loop_array, step_count)
            stop_condition_list = np.array([])

            # makes the waveforms repeat in loops
        self._pulser.load_sequence(  # digital_sequences=dosequence,
            loops_list=loops,
            segment_map=segment_map,
            stop_condition_list=stop_condition_list,
        )
        self._pulser.start_card(card_idx)
        self._pulser.arm_trigger(card_idx)


    def play_ttl(
            self,
            seq_len=2,
            laser_out=True,
            apd_sig_out=False,
            apd_ref_out=False,
            rep=100,
            trigger=True):
        """
        all are in microseconds
        This thing replaces the waveform test.
        It can switch the ttls according to the booleans
        True = high and False = low
        laser= True and
        apd, apd_ref = False ...means that the laser is on and the counts get to the very last port of the ttls (PFI 7)
        which is used in the cw mode
        We choose a small clock-rate for the awg because we dont care about the resolution
        The sequence len can be short but not too short.. Careful with the combination of seq len and clock rate!
        We need at least 2 step_counts to be able to use the software trigger
        """
        card_idx = 1
        self._pulser.set_sample_rate(card_idx, int(180 * 1e6)) # if the sampling rate is too low it throws an error
        self.clk_rate = self._pulser.get_sample_rate(card_idx) / 1e6  # To adapt it to the microsecond timing

        step_count = 2
        # print("step_count: ", step_count)  # Without unit
        timeaxis = self.time_axis(seq_len)
        if laser_out:  # high
            t_centrewidth_list_laser = self.convert_laser_val([2*seq_len, 0, 0])
        else:
            t_centrewidth_list_laser = self.convert_laser_val([-1, seq_len*2, 0])
        if apd_sig_out:  # high
            t_centrewidth_list_apd = self.convert_apd_val([0, seq_len*2])
        else:
            t_centrewidth_list_apd = self.convert_apd_val([seq_len*2, 0])
        if apd_ref_out:  # high
            t_centrewidth_list_apd_ref = self.convert_apd_val([0, seq_len*2])
        else:
            t_centrewidth_list_apd_ref = self.convert_apd_val([seq_len*2, 0])

        for step in range(step_count):
            laser_waveform = self._pulser.box_envelope(timeaxis, t_centrewidth_list_laser)
            apd_waveform = self._pulser.box_envelope(timeaxis, t_centrewidth_list_apd)  #
            apd_ref_waveform = self._pulser.box_envelope(timeaxis, t_centrewidth_list_apd_ref)  #
            digitals = {"laser": laser_waveform, "apd_sig": apd_waveform, "apd_read": apd_ref_waveform}  # read/sig?

            self._pulser.load_waveform(
                digital_pulses=digitals,
                digital_output_map={0: [0, 1, 2]}
            )
        if trigger:
            self._pulser.configure_ORmask(card_idx, None)
            self._pulser.configure_ANDmask(card_idx, None)
            array = np.array([0x40000000], dtype=np.int64)
            stop_condition_list = np.repeat(array, step_count)
            loops = np.ones(step_count, dtype=np.int64)

        else:
            self._pulser.configure_ORmask(card_idx, "immediate")
            self._pulser.configure_ANDmask(card_idx, None)
            loop_array = np.array([rep], dtype=np.int64)
            loops = np.repeat(loop_array, step_count)
            stop_condition_list = np.array([])

        self._pulser.load_sequence(
            loops_list=loops,
            stop_condition_list=stop_condition_list,)

        self._pulser.start_card(card_idx)
        self._pulser.arm_trigger(card_idx)
        self._pulser.send_software_trig(card_idx)
        return 0

    def convert_mw_val_rabi(self, mw_times, current_len):
        # everything is in microseconds
        if len(mw_times) == 4: # Rabi
            t0 = mw_times[0] + current_len / 2
            width = current_len
        else:
            raise ValueError('mw_times has a weird length')
        return t0, width

    def convert_mw_val_ramsey(self, mw_times, current_wait_time): #mw_times = [start, distance_min, distance_max, steps, pulse length]
        # everything is in microseconds
        if len(mw_times) == 5: # ramsey
            t0 = mw_times[0] + mw_times[4] / 2
            t1 = mw_times[0] + current_wait_time + mw_times[4] + 0.5 * mw_times[4]
            width = mw_times[4]
        else:
            raise ValueError('mw_times has a weird length')
        return t0, t1, width


    def convert_laser_val(self, laser_times): # initialisation time, wait time, 2nd pulse time
        # everything is in microseconds
        if len(laser_times) == 3:
            t0_1 = laser_times[0] / 2
            width_1 = laser_times[0]
            t0_2 = laser_times[0] + laser_times [1] + laser_times [2] / 2
            width_2 = laser_times[2]
        else:
            raise ValueError(
                "The laser_times array has a weird length! It should have [initialisation time, wait time, 2nd pulse time]")
            t0_1, t0_2, width_1, width_2 = 0
        t_centrewidth_list_laser = np.array([[t0_1, width_1], [t0_2, width_2]])
        return t_centrewidth_list_laser
    def convert_apd_val(self, apd_times): # [time to start, length]
        # everything is in microseconds
        if len(apd_times) == 2:
            t0_1 = apd_times[0] + apd_times[1] / 2
            width_1 = apd_times[1]
        else:
            raise ValueError(
                "The apd_(ref)_times array has a weird length! It should have [time to start, length]")
            t0_1, t0_2 = 0
        t_centrewidth_list_apd = np.array([[t0_1, width_1]])
        return t_centrewidth_list_apd

    def convert_mw_sweep(self, mw_times_sweep): # [mw_wait_time, mw_pulse_time]
        # everything is in microseconds
        if len(mw_times_sweep) == 2:
            t0_1 = mw_times_sweep[0] + mw_times_sweep[1] / 2
            width_1 = mw_times_sweep[1]
        else:
            raise ValueError(
                "The apd_(ref)_times array has a weird length! It should have [time to start, length]")
            t0_1, t0_2 = 0
        t_centrewidth_list_mw_sweep = np.array([[t0_1, width_1]])
        return t_centrewidth_list_mw_sweep


    def time_axis(self, seq_len=5): # we need the axis for iq stuff, it depends on the sequence length
        # samples = len of waveform ...how can we fix that
        # this includes waveform_padding
        # seq_len = seq_len * 1e-6 # this should be in seconds
        self.samples = self.get_samples(seq_len)
        self.sample_rate = self.get_clk_rate() / 1e6 # To adapt it to the microsecond timing
        self.time_ax = np.linspace(0, self.samples/ self.sample_rate, self.samples)
        return self.time_ax # in units of microseconds

    def get_samples(self, seq_len=5):
        # seq_len = seq_len*1e-6
        self.sample_rate = self.get_clk_rate() / 1e6 # To adapt it to the microsecond timing
        samples = self._pulser.waveform_padding((seq_len * self.sample_rate))
        # this includes waveform_padding
        return samples # in units of microseconds

    def calc_sample_time(self):
        card_idx = 1
        clk_rate = self._pulser.get_sample_rate(card_idx)
        sample_time = 1 / clk_rate
        return sample_time, clk_rate

    def calc_real_time(self, microsec):
        sample_time, clk_rate = self.calc_sample_time()
        sec = microsec * 1e-6
        # This must also include the waveform padding
        samples_number = int(clk_rate * sec)
        real_duration = samples_number * sample_time
        samples = self.waveform_padding((sec * clk_rate))
        time_ax = np.linspace(0, samples / clk_rate, samples)
        return samples_number, real_duration

    def get_start_stop_arrays(self, t_centrewidths):
        """
        @param np.ndarray timeaxis: the time axis
        @param np.ndarray t_centrewidths: 2D matrix. Each row contains the central time and width of the square pulse.
        @return np.ndarray: the pulse train.
        """
        card_idx = 1
        clk_rate = self._pulser.get_sample_rate(card_idx)
        timeaxis = np.linspace(0, 20, 21)
        tstart = t_centrewidths[:, 0] - t_centrewidths[:, 1] / 2
        tend = t_centrewidths[:, 0] + t_centrewidths[:, 1] / 2
        # start_edges = timeaxis[None, :] >= tstart[:, None]
        # end_edges = timeaxis[None, :] <= tend[:, None]
        start_array = []
        stop_array = []
        start_array.append(tstart)
        stop_array.append(tend)
        # pulses = np.logical_and(start_edges, end_edges).astype(np.float64)
        # The values in the start_array give the index of the last zero before the edge.
        # The values in the stop_array give the index of the last 1 before the end of the pulse
        # both are indices referring to the time_axis
        return start_array, stop_array, timeaxis

    def calc_len_error(self, t_centerwidths, seq_len):

        real_lenghts = []
        for i in range(len(t_centerwidths)):
            if (t_centerwidths[i][1] % 2) == 0: # even
                if t_centerwidths[i][1] == 0:
                    pulse_len = t_centerwidths[i][1]
                else:
                    pulse_len = t_centerwidths[i][1] + 1
            else: # odd
                pulse_len = t_centerwidths[i][1]
            real_lenghts.append(pulse_len)
        # convert the number of samples to time:
        sample_time, clk_rate = self.calc_sample_time()
        time_axis = self.time_axis(seq_len)


        return real_lenghts

    def build_seq_map(self, step_count):
        seq_map = []
        array = np.linspace(0, step_count - 1, step_count)
        e = 0
        o = 1
        for n in range(step_count):
            if (n % 2) == 0:
                seq_map.append(int(array[e]))
                e = e + 1
            else:
                seq_map.append(int(array[-(o)]))
                o = o + 1

        return seq_map

    def get_seq_map(self):
        return self.seq_map




# ###############     OLD STUFF!  Measure delay time of the system       ###############
    # def delay_sweep_old(
    #     self,
    #     seq_len,
    #     laser_times,
    #     apd_times,
    #     apd_ref,
    #     mw_times,
    #     rep=100,
    #     trigger=True,
    #     mw_pulse=False,
    # ):
    #     """
    #     seq_len =   whole length of the sequence in microseconds
    #     laser_times = [laser_in, wait, laser_re] all in microseconds
    #             laser_in:   length of initialisation laser pulse (it starts right at the beginning) in microseconds
    #             wait:       time between the end of the initialisation and the reinitialisation in microseconds
    #             laser_re:   length of the reinitialisation pulse in microseconds
    #     apd_times = [length, min_start, max_start, steps]
    #             length:      time to wait until the pulse should start
    #             min_start:    minimal length of the apd pulse
    #             max_start:    maximal length of the apd pulse
    #             steps:      step size of the apd pulse sweep
    #
    #     apd_ref = [time to start, length of the pulse] both in microseconds
    #     mw_times = [mw_wait_time, mw_pulse_time] all in microseconds
    #             mw_wait_time: time between the end of the first laser pulse and the beginning of the microwave (delay)
    #             mw_pulse_time: length of the mw pulse in between the two laser pulses
    #     rep =       repetition time when no trigger is used
    #     mw_pulse =  With or without MW pulse in the middle
    #     trigger =   either True: software trigger or False: no trigger at all
    #     """
    #     seq_len_s = seq_len * 1e-6
    #     laser_times_s = np.multiply(laser_times, 1e-6)
    #     mw_times_s = np.multiply(mw_times, 1e-6)
    #     apd_times_s = np.multiply(apd_times, 1e-6)
    #     apd_ref_s = np.multiply(apd_ref, 1e-6)
    #
    #     # So far we just use card1
    #     card_idx = 1
    #
    #     # Set up the right channels:
    #     mw_outchan = 2  # analog     A0 for card 1   Input at MW generator is +- 0.5 V
    #     laser_outchan = 0  # digital    X0 Input of the laser is 0-1 V
    #     apd_outchan = 1  # digital    X1
    #     apd_ref_outchan = 2  # ADP reference this happens on Card1
    #
    #     clk_rate = self._pulser.get_sample_rate(card_idx)
    #     segment_map = []
    #
    #     # build array with all the desired apd_start values
    #     apd_start = []
    #     step_count = int(((apd_times[2] - apd_times[1]) / apd_times[3]) + 1)
    #     print("step_count: ", step_count)
    #
    #     for i in range(step_count):
    #         apd_start.append((apd_times_s[1] + (i * apd_times_s[3])))
    #     # creating empty arrays
    #     dosequence = []
    #     aosequence = []
    #
    #     # Create laser and apd waveform as they are constant (put it in microseconds)
    #     laser_do_waveform = self.laser_seq(seq_len, laser_times)
    #     apd_ref_do_waveform = self.apd_ref(
    #         seq_len, apd_ref
    #     )
    #
    #     for i in apd_start:
    #         apd_wait_time = np.zeros(int(round(clk_rate * i)))
    #         apd_on = np.ones(int(round(clk_rate * apd_times_s[0])))
    #         apd_off = np.zeros(
    #             int(round(clk_rate * (seq_len_s - (i + apd_times_s[0]))))
    #         )
    #         apd_do_waveform = np.concatenate((apd_wait_time, apd_on, apd_off))
    #         len_apd_do_waveform = self._pulser.waveform_padding(len(apd_do_waveform))
    #         difference = len_apd_do_waveform - (clk_rate * seq_len_s)
    #         apd_do_waveform = np.concatenate(
    #             (apd_do_waveform, np.zeros(int(difference)))
    #         )
    #         if mw_pulse:
    #             mw_wait = np.zeros(
    #                 int(round(clk_rate * (laser_times_s[0] + mw_times_s[0])))
    #             )
    #             mw_pulse_array = np.ones(int(round(clk_rate * mw_times_s[1])))
    #             mw_off = np.zeros(
    #                 int(
    #                     round(
    #                         clk_rate
    #                         * (
    #                             seq_len_s
    #                             - (laser_times_s[0] + mw_times_s[0] + mw_times_s[1])
    #                         )
    #                     )
    #                 )
    #             )
    #
    #             mw_ao_waveform = np.concatenate((mw_wait, mw_pulse_array, mw_off))
    #
    #             len_mw_ao_waveform = self._pulser.waveform_padding(len(mw_ao_waveform))
    #             difference = len_mw_ao_waveform - (clk_rate * seq_len_s)
    #             mw_ao_waveform = np.concatenate(
    #                 (mw_ao_waveform, np.zeros(int(difference)))
    #             )
    #
    #             ao_dict = {mw_outchan: mw_ao_waveform}
    #
    #         else:  # if no MW is needed
    #             d_off = np.zeros(int(round(clk_rate * seq_len_s)))
    #             len_d_off = self._pulser.waveform_padding(len(d_off))
    #             difference = len_d_off - (clk_rate * seq_len_s)
    #             d_off = np.concatenate((d_off, np.zeros(int(difference))))
    #             # Creates a dict for the analog channel if no MW is needed
    #             ao_dict = {mw_outchan: d_off}
    #
    #         # Now all the waveformes are prepared and have the same length
    #         do_dict = {
    #             laser_outchan: laser_do_waveform,
    #             apd_outchan: apd_do_waveform,
    #             apd_ref_outchan: apd_ref_do_waveform,
    #         }
    #         dosequence.append(do_dict)
    #         aosequence.append(ao_dict)
    #
    #     self.aosequence = aosequence
    #     self.dosequence = dosequence
#
#
#     # Otherwise its not possible to send 3.5 and 3V at the same time
#     digital_output_map = {0: [laser_outchan, apd_outchan], 1: [apd_ref_outchan]}
#     for aos, dos in zip(aosequence, dosequence):
#         # Loads waveforms to the awg
#         self._pulser.load_waveform(
#             iq_dictionary=aos, # Change this to the new format!
#             digital_pulses=dos, # Change this to the new format!
#             digital_output_map=digital_output_map,
#         )
#
#     if trigger:
#         # This sequence waits for a software trigger to start playing and moving to the next step.
#         self._pulser.configure_ORmask(card_idx, None)
#         self._pulser.configure_ANDmask(card_idx, None)
#         array = np.array([0x40000000], dtype=np.int64)
#         stop_condition_list = np.repeat(array, step_count)
#         loops = np.ones((len(aosequence)), dtype=np.int64)
#
#     else:
#         # This sequence immediately starts after the sequences are loaded
#         self._pulser.configure_ORmask(card_idx, "immediate")
#         self._pulser.configure_ANDmask(card_idx, None)
#         loop_array = np.array([rep], dtype=np.int64)
#         loops = np.repeat(loop_array, step_count)
#         stop_condition_list = np.array([])
#
#         # makes the waveforms repeat in loops
#     self._pulser.load_sequence(  # digital_sequences=dosequence,
#         loops_list=loops,
#         segment_map=segment_map,
#         stop_condition_list=stop_condition_list,
#     )
#     self._pulser.start_card(card_idx)
#     self._pulser.arm_trigger(card_idx)


# def delay_sweep_ref_old(
        #         self,
        #         seq_len,
        #         laser_times,
        #         apd_times,
        #         apd_ref,
        #         mw_times,
        #         rep=100,
        #         trigger=bool,
        #         mw_pulse=bool,
        # ):
        #     """
        #     seq_len =   whole length of the sequence in microseconds
        #     laser_times = [laser_in, wait, laser_re] all in microseconds
        #             laser_in:   length of initialisation laser pulse (it starts right at the beginning) in microseconds
        #             wait:       time between the end of the initialisation and the reinitialisation in microseconds
        #             laser_re:   length of the reinitialisation pulse in microseconds
        #     apd_times = [length, min_start, max_start, steps] this will be connected to the reference apd pulse but still doing the sweep
        #             length:      time to wait until the pulse should start
        #             min_start:    minimal length of the apd pulse
        #             max_start:    maximal length of the apd pulse
        #             steps:      step size of the apd pulse sweep
        #     apd_ref = [time to start, length of the pulse] both in microseconds this will be connected to the real apd put it somewhere and dont sweep it
        #     mw_times = [mw_wait_time, mw_pulse_time] all in microseconds
        #             mw_wait_time: time between the end of the first laser pulse and the beginning of the microwave (delay)
        #             mw_pulse_time: length of the mw pulse in between the two laser pulses
        #     rep =       repetition time when no trigger is used
        #     mw_pulse =  With or without MW pulse in the middle
        #     trigger =   either True: software trigger or False: no trigger at all
        #     """
        #     seq_len_s = seq_len * 1e-6
        #     laser_times_s = np.multiply(laser_times, 1e-6)
        #     mw_times_s = np.multiply(mw_times, 1e-6)
        #     apd_times_s = np.multiply(apd_times, 1e-6)
        #     apd_ref_s = np.multiply(apd_ref, 1e-6)
        #
        #     # Sofar we just use card1
        #     card_idx = 1
        #
        #     # Set up the right channels:
        #     mw_outchan = 2  # analog     A0 for card 1   Input at MW generator is +- 0.5 V
        #     laser_outchan = 0  # digital    X0 Input of the laser is 0-1 V
        #     apd_outchan = 1  # digital    X1
        #     apd_ref_outchan = 2  # ADP reference this happens on Card1
        #
        #     clk_rate = self._pulser.get_sample_rate(card_idx)
        #     segment_map = []
        #
        #     # build array with all the desired apd_start values
        #     apd_start = []
        #     step_count = int(((apd_times[2] - apd_times[1]) / apd_times[3]) + 1)
        #     print("step_count: ", step_count)
        #
        #     for i in range(step_count):
        #         apd_start.append((apd_times_s[1] + (i * apd_times_s[3])))
        #     # creating empty arrays
        #     dosequence = []
        #     aosequence = []
        #
        #     # Create laser and apd waveform as they are constant (put it in microseconds)
        #     laser_do_waveform = self.laser_seq(seq_len, laser_times)
        #     apd_ref_do_waveform = self.apd_ref(
        #         seq_len, apd_ref
        #     )
        #
        #     for i in apd_start:
        #         apd_wait_time = np.zeros(int(round(clk_rate * i)))
        #         apd_on = np.ones(int(round(clk_rate * apd_times_s[0])))
        #         apd_off = np.zeros(
        #             int(round(clk_rate * (seq_len_s - (i + apd_times_s[0]))))
        #         )
        #         apd_do_waveform = np.concatenate((apd_wait_time, apd_on, apd_off))
        #         len_apd_do_waveform = self._pulser.waveform_padding(len(apd_do_waveform))
        #         difference = len_apd_do_waveform - (clk_rate * seq_len_s)
        #         apd_do_waveform = np.concatenate(
        #             (apd_do_waveform, np.zeros(int(difference)))
        #         )
        #         if mw_pulse:
        #             mw_wait = np.zeros(
        #                 int(round(clk_rate * (laser_times_s[0] + mw_times_s[0])))
        #             )
        #             mw_pulse_array = np.ones(int(round(clk_rate * mw_times_s[1])))
        #             mw_off = np.zeros(
        #                 int(
        #                     round(
        #                         clk_rate
        #                         * (
        #                                 seq_len_s
        #                                 - (laser_times_s[0] + mw_times_s[0] + mw_times_s[1])
        #                         )
        #                     )
        #                 )
        #             )
        #             mw_ao_waveform = np.concatenate((mw_wait, mw_pulse_array, mw_off))
        #             len_mw_ao_waveform = self._pulser.waveform_padding(len(mw_ao_waveform))
        #             difference = len_mw_ao_waveform - (clk_rate * seq_len_s)
        #             mw_ao_waveform = np.concatenate(
        #                 (mw_ao_waveform, np.zeros(int(difference)))
        #             )
        #             analogs = {"i_chan": mw_ao_waveform, "q_chan": np.zeros(len(mw_ao_waveform))}
        #
        #         else:  # if no MW is needed
        #             d_off = np.zeros(int(round(clk_rate * seq_len_s)))
        #             len_d_off = self._pulser.waveform_padding(len(d_off))
        #             difference = len_d_off - (clk_rate * seq_len_s)
        #             d_off = np.concatenate((d_off, np.zeros(int(difference))))
        #             # Creates a dict for the analog channel if no MW is needed
        #             analogs = {"i_chan": mw_ao_waveform, "q_chan": np.zeros(len(mw_ao_waveform))}
        #
        #         # Now all the waveformes are prepared and have the same length
        #         do_dict = {
        #             laser_outchan: laser_do_waveform,
        #             apd_outchan: apd_ref_do_waveform,
        #             apd_ref_outchan: apd_do_waveform,
        #         }
        #         dosequence.append(do_dict)
        #         aosequence.append(ao_dict)
        #
        #     ###
        #     # Otherwise its not possible to send 3.5 and 3V at the same time
        #     # digital_output_map = {0: [laser_outchan], 1: [apd_outchan]}
        #     digital_output_map = {0: [laser_outchan, apd_outchan], 1: [apd_ref_outchan]}
        #     for aos, dos in zip(aosequence, dosequence):
        #         # Loads waveforms to the awg
        #         self._pulser.load_waveform(
        #             iq_dictionary=aos,
        #             digital_pulses=dos,
        #             digital_output_map=digital_output_map,
        #         )
        #
        #     if trigger == True:
        #         # This sequence waits for a software trigger to start playing and moving to the next step.
        #         self._pulser.configure_ORmask(card_idx, None)
        #         self._pulser.configure_ANDmask(card_idx, None)
        #         array = np.array([0x40000000], dtype=np.int64)
        #         stop_condition_list = np.repeat(array, step_count)
        #         loops = np.ones((len(aosequence)), dtype=np.int64)
        #
        #     else:
        #         # This sequence immediately starts after the sequences are loaded
        #         self._pulser.configure_ORmask(card_idx, "immediate")
        #         self._pulser.configure_ANDmask(card_idx, None)
        #         loop_array = np.array([rep], dtype=np.int64)
        #         loops = np.repeat(loop_array, step_count)
        #         stop_condition_list = np.array([])
        #
        #         # makes the waveforms repeat in loops
        #     self._pulser.load_sequence(  # digital_sequences=dosequence,
        #         loops_list=loops,
        #         segment_map=segment_map,
        #         stop_condition_list=stop_condition_list,
        #     )
        #     self._pulser.start_card(card_idx)
        #     self._pulser.arm_trigger(card_idx)


#     def do_cw_old(self, seq_len):
#         '''
#         Play ones to the APD channel to make sure the counts get through and can be used in a cw measurement.
#         Here I send zeros to all the other channels: Laser, mw, and APD ref
#         '''
#         rep = 100
#         dosequence = []
#         aosequence = []
#         segment_map = []
#         #Channels
#         laser_outchan = 0  # Laser
#         apd_outchan = 1
#         apd_ref_outchan = 2 # ADP reference this happens on Card1
#         mw_outchan = 2  # Laser
#
#         card_idx = 1
#         # convert to seconds
#         seq_len_s = seq_len * 1e-6
#         clk_rate = self._pulser.get_sample_rate(card_idx)
#         #Create laser and apd waveform as they are costant
#         array_zeros = self.just_zeros(seq_len)
#         array_ones = self.just_ones(seq_len)
#         # Now build the  do_dict and ao_dict
#
#         do_dict = {laser_outchan: array_zeros, apd_outchan: array_ones, apd_ref_outchan: array_zeros}
#         ao_dict = {mw_outchan: array_zeros}
#         aosequence.append(ao_dict)
#         dosequence.append(do_dict)
#         aosequence.append(ao_dict)
#         dosequence.append(do_dict)
#         # print('aosequence:', aosequence)
#         # print('dosequence:', dosequence)
#         # digital_output_map = {0: [laser_outchan, apd_outchan]}
#         digital_output_map = {0: [laser_outchan, apd_outchan], 1: [apd_ref_outchan]}
#         # print('digital_output_map: ', digital_output_map)
#         for aos, dos in zip(aosequence, dosequence):
#             # Loads waveforms to the awg
#             #check_out = self.check_len(laser_do_waveform, apd_do_waveform, mw_ao_waveform, mw_pulse=True)
#             self._pulser.load_waveform(ao_waveform_dictionary=aos, do_waveform_dictionary=dos,
#                                    digital_output_map=digital_output_map)
#
#         # This sequence waits for a software trigger to start playing and moving to the next step.
#         self._pulser.configure_ORmask(card_idx, None)
#         self._pulser.configure_ANDmask(card_idx, None)
#         array = np.array([0x40000000], dtype=np.int64)
#         stop_condition_list = np.repeat(array, 2)
#         loops = np.ones((len(aosequence)), dtype=np.int64)
#
#
#         self._pulser.load_sequence(  # digital_sequences=dosequence,
#             loops_list=loops, segment_map=segment_map, stop_condition_list=stop_condition_list)
#         self._pulser.start_card(card_idx)
#         self._pulser.arm_trigger(card_idx)
#         self.trigger()
#
#     def post_measurement_old(self, seq_len):
#         '''
#         Play zeros to both of the APD channels to make sure they don't do weird stuff when the measurement is done.
#         Here I send zeros to all channels: Laser, mw, APD and APD ref
#         '''
#         rep = 100
#         dosequence = []
#         aosequence = []
#         segment_map = []
#         # Channels
#         laser_outchan = 0  # Laser
#         apd_outchan = 1
#         apd_ref_outchan = 2  # ADP reference this happens on Card1
#         mw_outchan = 2  # Laser
#
#         card_idx = 1
#         # convert to seconds
#         seq_len_s = seq_len * 1e-6
#         clk_rate = self._pulser.get_sample_rate(card_idx)
#         # Create laser and apd waveform as they are costant
#         array_zeros = self.just_zeros(seq_len)
#         # print('array_zeros:', array_zeros)
#         # Now build the  do_dict and ao_dict
#
#         do_dict = {laser_outchan: array_zeros, apd_outchan: array_zeros, apd_ref_outchan: array_zeros}
#         ao_dict = {mw_outchan: array_zeros}
#         aosequence.append(ao_dict)
#         dosequence.append(do_dict)
#         aosequence.append(ao_dict)
#         dosequence.append(do_dict)
#         # print('aosequence:', aosequence)
#         # print('dosequence:', dosequence)
#         # digital_output_map = {0: [laser_outchan, apd_outchan]}
#         digital_output_map = {0: [laser_outchan, apd_outchan], 1: [apd_ref_outchan]}
#         # print('digital_output_map: ', digital_output_map)
#         for aos, dos in zip(aosequence, dosequence):
#             # Loads waveforms to the awg
#             # check_out = self.check_len(laser_do_waveform, apd_do_waveform, mw_ao_waveform, mw_pulse=True)
#             self._pulser.load_waveform(ao_waveform_dictionary=aos, do_waveform_dictionary=dos,
#                                        digital_output_map=digital_output_map)
#
#         # This sequence waits for a software trigger to start playing and moving to the next step.
#         self._pulser.configure_ORmask(card_idx, None)
#         self._pulser.configure_ANDmask(card_idx, None)
#         array = np.array([0x40000000], dtype=np.int64)
#         stop_condition_list = np.repeat(array, 2)
#         loops = np.ones((len(aosequence)), dtype=np.int64)
#
#         self._pulser.load_sequence(  # digital_sequences=dosequence,
#             loops_list=loops, segment_map=segment_map, stop_condition_list=stop_condition_list)
#         self._pulser.start_card(card_idx)
#         self._pulser.arm_trigger(card_idx)
#         self.trigger()
#
#
#     def setup_and_play(self, laser_times, mw_times, apd_times, rep=100, rest=100, clk_rate=100, mw_pulse=bool,
#                        trigger=bool):
#         '''
#             laser_times, mw_times, apd_times should be arrays with the right values inside
#             laser_times = [laser_in, wait, laser_re] all in microseconds
#                 laser_in:   length of initialisation laser pulse (it starts right at the beginning) in microseconds
#                 wait:       time between the end of the initialisation and the reinitialisation in microseconds
#                 laser_re:   length of the reinitialisation pulse in microseconds
#             mw_times = [mw_wait_time, mw_pulse_time] all in microseconds
#                 mw_wait_time: time between the end of the first laser pulse and the beginning of the microwave (delay)
#                 mw_pulse_time: length of the mw pulse in between the two laser pulses
#             apd_times = [steps, apd_width]
#                 steps:      step size of the apd pulse sweep
#                 apd_width:   time in microseconds where the apd is on (length of the pulse)
#             rep:        If there is no software trigger, rep gives the number of repetitions the waveform is played
#             rest:       Time between single sequences (just to be able to distinguish between beginning of 1 laser pulse and end of second laser lulse)
#             clk_rate:  set to 100 if None
#             mw_pulse:   True: plays MW pulse and False sets the MW sequence to zero
#             trigger:    True if the software trigger should be used and False if the loop array should define the length
#
#         These parameters work:
#         laser_in=2,wait=1,laser_re=7,
#         steps=1, apd_width=2
#          rest=100,
#          mw_wait_time=0.5, mw_pulse_time=0.25,
#          rep=10,
#          mw_pulse=True, trigger=True)
#         '''
#
#         laser_times = np.multiply(laser_times, 1e-6)
#         mw_times = np.multiply(mw_times, 1e-6)
#         apd_times = np.multiply(apd_times, 1e-6)
#
#         # Sofar we just use card1
#         card_idx = 1
#
#         # Set up the right channels:
#         # Microwave: analog channel 2 and 3 for I and Q parameter. For now we just need I
#         # Laser: digital channel X0
#         # Photon counter: digital channel X1 and X2 (for reference)
#         mw_I = 2  # analog     A0 for card 1   Input at MW generator is +- 0.5 V
#         mw_Q = 3  # analog     A1 for card 1   --> add later
#         laser_channel = 0  # digital    X0 Input of the laser is 0-1 V
#         apd_channel1 = 1  # digital    X1
#         # build an array for that [0,2,1]
#         channels = np.array([laser_channel, mw_I, apd_channel1])
#
#         # if not given differently put the clk_rate to 100ms
#         if clk_rate == 100:
#             clk_rate = MEGA(100)
#             # print('clk_rate not changed')
#             print('clk_rate in 1/s: ', clk_rate)
#         else:
#             clk_rate = MEGA(clk_rate)
#             # print('clk_rate changed')
#             print('clk_rate in 1/s: ', clk_rate)
#
#         # This one sets up the right clk_rate
#         self.set_up_pulser(card_idx, clk_rate)
#         print('AWG is ready')
#
#         if rest == 100:
#             rest = 100 * 1e-6
#             print('break in between single sequences set to 100ms')
#         else:
#             rest *= 1e-6
#             print('break in between single sequences set to the new value')
#         if rep == 100:
#             print('set rep to 100')
#         else:
#             print('set rep to the new value')
#
#         segment_map = []
#         # calculate the whole length of the waveform --> seq_len
#         # The parameters for the laser determine how long the waveform is
#         seq_len = laser_times[0] + laser_times[1] + laser_times[2] + rest  # in s
#
#         # use all these parameters and put it to the general sweep function:
#         self.delay_sweep_general(card_idx, clk_rate, seq_len, laser_times[0], laser_times[1], laser_times[2], rest,
#                                  apd_times[0], apd_times[1], mw_times[0],
#                                  mw_times[1], rep, channels, mw_pulse, trigger, segment_map)
#
#     # This is the main general block...works fine with laser, mw and apd
#     def delay_sweep_general(self, card_idx, clk_rate, seq_len, laser_in, wait, laser_re, rest, steps, apd_width,
#                             mw_wait_time, mw_pulse_time, rep, channels, mw_pulse=bool, trigger=bool, segment_map=[]):
#         '''
#             This can play a laser sequence, microwaves and the apd readout at the same time
#             card_idx:   so far just 1
#             clk_rate:   something like 100
#             seq_len:    length of one whole sequence
#             laser times:
#                 laser_in:   length of initialisation laser pulse (it starts right at the beginning) microseconds
#                 wait:       time between the end of the initialisation and the reinitialisation in microseconds
#                 laser_re:   length of the reinitialisation pulse in microseconds
#             rest:       time in microseconds after the reinitialisation pulse to have a break before the next sequences starts
#             apd times
#                 start:      where in time should the apd pulse start shifting --> this should be the time as the start of the laser reinitialisation = laser_in + wait
#                 steps:      step size of the apd pulse sweep ( this is in microseconds right?)
#                 apd_width:   time in microseconds where the apd is on (length of the pulse)
#             MW times
#                 mw_wait_time: time between the end of the first laser pulse and the beginning of the microwave (delay)
#                 mw_pulse_time: length of the mw pulse in between the two laser pulses
#             rep:        If there is no software trigger, rep gives the number of repetitions the waveform is played
#             channels:   Connectors of the AWG (hardcoded)
#             mw_pulse:   True: sends MW pulse, False: sends zeros instead
#             trigger:    True if the software trigger should be used and False if the loop array should define the length
#             '''
#
#         apd_start_val = laser_in + wait
#         stop = (seq_len - rest) - apd_width
#         laser_outchan = channels[0]  # Laser
#         apd_outchan = channels[2]  # APD
#         mw_outchan = channels[1]
#
#         # build array with all the desired apd_start values
#         apd_start = []
#         step_count = math.floor((stop - apd_start_val) / steps)
#         for i in range(step_count):
#             apd_start.append((apd_start_val + (i * steps)))
#
#         # creating empty arrays
#         dosequence = []
#         aosequence = []
#
#         # Laser waveform does not change the whole time
#         first_pulse = np.ones(self._pulser.waveform_padding(clk_rate * laser_in))
#         laser_wait_time = np.zeros(self._pulser.waveform_padding(clk_rate * wait))
#         sec_pulse = np.ones(self._pulser.waveform_padding(clk_rate * laser_re))
#         rest_array = np.zeros(self._pulser.waveform_padding(clk_rate * rest))
#         laser_do_waveform = np.concatenate((first_pulse, laser_wait_time, sec_pulse, rest_array))
#
#         for i in apd_start:
#             apd_wait_time = np.zeros(self._pulser.waveform_padding(clk_rate * i))
#             apd_on = np.ones(self._pulser.waveform_padding(clk_rate * apd_width))
#             apd_off = np.zeros(self._pulser.waveform_padding(clk_rate * (seq_len - (i + apd_width))))
#             apd_do_waveform = np.concatenate((apd_wait_time, apd_on, apd_off))
#
#             if mw_pulse == True:
#                 mw_wait = np.zeros(self._pulser.waveform_padding(clk_rate * (laser_in + mw_wait_time)))
#                 mw_pulse_array = np.ones(self._pulser.waveform_padding(clk_rate * mw_pulse_time))
#                 mw_off = np.zeros(self._pulser.waveform_padding(clk_rate * seq_len - (laser_in + mw_wait_time + mw_pulse_time)))
#                 mw_pulse_array = np.subtract(mw_pulse_array, 0.5)  # because the MW takes only +-0.5V
#                 mw_ao_waveform = np.concatenate((mw_wait, mw_pulse_array, mw_off))
#
#                 # Make sure the arrays have the same length:
#                 # Check which array is the longest and append the two others:
#                 max_val = np.max(np.array([len(laser_do_waveform), len(apd_do_waveform), len(mw_ao_waveform)]))
#                 if max_val == 0:
#                     # print('laser waveforms is the longest')
#                     difference_mw = len(laser_do_waveform) - len(mw_ao_waveform)
#                     difference_apd = len(laser_do_waveform) - len(apd_do_waveform)
#                     rest_array_mw = np.zeros(int(difference_mw))
#                     rest_array_apd = np.zeros(int(difference_apd))
#                     mw_ao_waveform = np.concatenate((mw_ao_waveform, rest_array_mw))
#                     apd_do_waveform = np.concatenate((apd_do_waveform, rest_array_apd))
#
#                 elif max_val == 1:
#                     # print('apd waveforms is the longest')
#                     difference_mw = len(apd_do_waveform) - len(mw_ao_waveform)
#                     difference_laser = len(apd_do_waveform) - len(laser_do_waveform)
#                     rest_array_mw = np.zeros(int(difference_mw))
#                     rest_array_laser = np.zeros(int(difference_laser))
#                     mw_ao_waveform = np.concatenate((mw_ao_waveform, rest_array_mw))
#                     laser_do_waveform = np.concatenate((laser_do_waveform, rest_array_laser))
#
#                 else:
#                     # print('mw waveforms is the longest')
#                     difference_apd = len(mw_ao_waveform) - len(apd_do_waveform)
#                     difference_laser = len(mw_ao_waveform) - len(laser_do_waveform)
#                     rest_array_apd = np.zeros(int(difference_apd))
#                     rest_array_laser = np.zeros(int(difference_laser))
#                     apd_do_waveform = np.concatenate((apd_do_waveform, rest_array_apd))
#                     laser_do_waveform = np.concatenate((laser_do_waveform, rest_array_laser))
#
#                 ao_dict = {mw_outchan: mw_ao_waveform}
#
#             else:  # if no MW is needed
#                 # check which array is longer laser or apd
#                 if len(apd_do_waveform) < len(laser_do_waveform):
#                     # print('laser waveform is longer than apd waveform')
#                     difference = len(laser_do_waveform) - len(apd_do_waveform)
#                     rest_array = np.zeros(int(difference))
#                     apd_do_waveform = np.concatenate((apd_do_waveform, rest_array))
#                     # Create zero array of same size:
#                     d_off = np.zeros(len(laser_do_waveform))
#                 elif len(laser_do_waveform) < len(apd_do_waveform):
#                     # print('apd waveform is the longer')
#                     difference = len(apd_do_waveform) - len(laser_do_waveform)
#                     rest_array = np.zeros(int(difference))
#                     laser_do_waveform = np.concatenate((laser_do_waveform, rest_array))
#                     # Create zero array of same size:
#                     d_off = np.zeros(len(apd_do_waveform))
#                 else:
#                     print('apd waveform and laser_waveform have the same size')
#                     # d_off = np.zeros(len(apd_do_waveform))
#                 # Creates a dict for the analog channel if no MW is needed
#                 ao_dict = {mw_outchan: d_off}
#
#             # Now all the waveformes are prepared and have the same length
#             do_dict = {laser_outchan: laser_do_waveform, apd_outchan: apd_do_waveform}
#             dosequence.append((do_dict))
#             aosequence.append((ao_dict))
#
#         digital_output_map = {0: [laser_outchan, apd_outchan]}  # directs to the right output channel analog channel0 directs to digital x0 and x1
#
#         for aos, dos in zip(aosequence, dosequence):
#             # Loads waveforms to the awg
#             self._pulser.load_waveform(ao_waveform_dictionary=aos, do_waveform_dictionary=dos,
#                                        digital_output_map=digital_output_map)
#
#         if trigger == True:
#             # This sequence waits for a software trigger to start playing and moving to the next step.
#             self._pulser.configure_ORmask(card_idx, None)
#             self._pulser.configure_ANDmask(card_idx, None)
#             array = np.array([0x40000000], dtype=np.int64)
#             stop_condition_list = np.repeat(array, step_count)
#             loops = np.ones((len(aosequence)), dtype=np.int64)
#
#         else:
#             # This sequence immediately starts after the sequences are loaded
#             self._pulser.configure_ORmask(card_idx, 'immediate')
#             self._pulser.configure_ANDmask(card_idx, None)
#             loop_array = np.array([rep], dtype=np.int64)
#             loops = np.repeat(loop_array, step_count)
#             stop_condition_list = np.array([])
#
#             # makes the waveforms repeat in loops
#         self._pulser.load_sequence(  # digital_sequences=dosequence,
#             loops_list=loops, segment_map=segment_map, stop_condition_list=stop_condition_list)
#         self._pulser.start_card(card_idx)
#         self._pulser.arm_trigger(card_idx)
#
#     def setup_and_play_new(self, laser_times, mw_times, apd_times, rep=100, rest=100, clk_rate=100, mw_pulse=bool,
#                        trigger=bool):
#         '''
#             laser_times, mw_times, apd_times should be arrays with the right values inside
#             laser_times = [laser_in, wait, laser_re] all in microseconds
#                 laser_in:   length of initialisation laser pulse (it starts right at the beginning) in microseconds
#                 wait:       time between the end of the initialisation and the reinitialisation in microseconds
#                 laser_re:   length of the reinitialisation pulse in microseconds
#             mw_times = [mw_wait_time, mw_pulse_time] all in microseconds
#                 mw_wait_time: time between the end of the first laser pulse and the beginning of the microwave (delay)
#                 mw_pulse_time: length of the mw pulse in between the two laser pulses
#             apd_times = [steps, apd_width]
#                 steps:      step size of the apd pulse sweep
#                 apd_width:   time in microseconds where the apd is on (length of the pulse)
#             rep:        If there is no software trigger, rep gives the number of repetitions the waveform is played
#             rest:       Time between single sequences (just to be able to distinguish between beginning of 1 laser pulse and end of second laser lulse)
#             clk_rate:  set to 100 if None
#             mw_pulse:   True: plays MW pulse and False sets the MW sequence to zero
#             trigger:    True if the software trigger should be used and False if the loop array should define the length
#
#         These parameters work:
#         laser_in=2,wait=1,laser_re=7,
#         steps=1, apd_width=2
#          rest=100,
#          mw_wait_time=0.5, mw_pulse_time=0.25,
#          rep=10,
#          mw_pulse=True, trigger=True)
#         '''
#
#         laser_times = np.multiply(laser_times, 1e-6)
#         mw_times = np.multiply(mw_times, 1e-6)
#         apd_times = np.multiply(apd_times, 1e-6)
#
#         # Sofar we just use card1
#         card_idx = 1
#
#         # Set up the right channels:
#         # Microwave: analog channel 2 and 3 for I and Q parameter. For now we just need I
#         # Laser: digital channel X0
#         # Photon counter: digital channel X1 and X2 (for reference)
#         mw_I = 2  # analog     A0 for card 1   Input at MW generator is +- 0.5 V
#         mw_Q = 3  # analog     A1 for card 1   --> add later
#         laser_channel = 0  # digital    X0 Input of the laser is 0-1 V
#         apd_channel1 = 1  # digital    X1
#         # build an array for that [0,2,1]
#         channels = np.array([laser_channel, mw_I, apd_channel1])
#
#         # if not given differently put the clk_rate to 100ms
#         if clk_rate == 100:
#             clk_rate = MEGA(100)
#             # print('clk_rate not changed')
#             print('clk_rate in 1/s: ', clk_rate)
#         else:
#             clk_rate = MEGA(clk_rate)
#             # print('clk_rate changed')
#             print('clk_rate in 1/s: ', clk_rate)
#
#         # This one sets up the right clk_rate
#         self.set_up_pulser(card_idx, clk_rate)
#         print('AWG is ready')
#
#         if rest == 100:
#             rest = 100 * 1e-6
#             print('break in between single sequences set to 100ms')
#         else:
#             rest *= 1e-6
#             print('break in between single sequences set to the new value')
#         if rep == 100:
#             print('set rep to 100')
#         else:
#             print('set rep to the new value')
#
#         segment_map = []
#         # calculate the whole length of the waveform --> seq_len
#         # The parameters for the laser determine how long the waveform is
#         seq_len = laser_times[0] + laser_times[1] + laser_times[2] + rest  # in s
#         apd_start_val = laser_times[0] + laser_times[1]
#         stop = (seq_len - rest) - apd_times[1]
#
#         # build array with all the desired apd_start values
#         apd_start = []
#         step_count = math.floor((stop - apd_start_val) / apd_times[0])
#         for i in range(step_count):
#             apd_start.append((apd_start_val + (i * apd_times[0])))
#
#         # use all these parameters and put it to the general sweep function:
#         self.delay_sweep_general_new(card_idx, clk_rate, seq_len, laser_times[0], laser_times[1], laser_times[2], rest,
#                                  apd_start, apd_times[1], mw_times[0],
#                                  mw_times[1],step_count, rep, channels, mw_pulse, trigger, segment_map)
#
#
#     # This is the main general block...works fine with laser, mw and apd NEW
#     def delay_sweep_general_new(self, card_idx, clk_rate, seq_len, laser_in, wait, laser_re, rest, apd_start, apd_width,
#                             mw_wait_time, mw_pulse_time, step_count, rep, channels, mw_pulse=bool, trigger=bool, segment_map=[]):
#         '''
#             This can play a laser sequence, microwaves and the apd readout at the same time
#             card_idx:   sofar just 1
#             clk_rate:   something like 100
#             seq_len:    length of one whole sequence
#             laser times:
#                 laser_in:   length of initialisation laser pulse (it starts right at the beginning) in ms
#                 wait:       time between the end of the initialisation and the reinitialisation in ms
#                 laser_re:   length of the reinitialisation pulse in ms
#             rest:       time in ms after the reinitialisation pulse to have a break before the next sequences starts
#             apd times
#                 start:      where in time should the apd pulse start shifting --> this should be the time as the start of the laser reinitialisation = laser_in + wait
#                 steps:      step size of the apd pulse sweep ( this is in ms right?)
#                 apd_width:   time in ms where the apd is on (length of the pulse)
#             MW times
#                 mw_wait_time: time between the end of the first laser pulse and the beginning of the microwave (delay)
#                 mw_pulse_time: length of the mw pulse in between the two laser pulses
#             rep:        If there is no software trigger, rep gives the number of repetitions the waveform is played
#             channels:   Connetors of the AWG (hardcoded)
#             mw_pulse:   True: sends MW pulse, False: sends zeros instead
#             trigger:    True if the software trigger should be used and False if the loop array should define the length
#             segment_map:If the order of different segments should be mixed (not there yet)
#             '''
#
#         laser_outchan = channels[0]  # Laser
#         apd_outchan = channels[2]  # APD
#         mw_outchan = channels[1]
#
#         # creating empty arrays
#         dosequence = []
#         aosequence = []
#
#         # Laser waveform does not change the whole time
#         first_pulse = np.ones(self._pulser.waveform_padding(clk_rate * laser_in))
#         laser_wait_time = np.zeros(self._pulser.waveform_padding(clk_rate * wait))
#         sec_pulse = np.ones(self._pulser.waveform_padding(clk_rate * laser_re))
#         rest_array = np.zeros(self._pulser.waveform_padding(clk_rate * rest))
#         laser_do_waveform = np.concatenate((first_pulse, laser_wait_time, sec_pulse, rest_array))
#
#         for i in apd_start:
#             apd_wait_time = np.zeros(self._pulser.waveform_padding(clk_rate * i))
#             apd_on = np.ones(self._pulser.waveform_padding(clk_rate * apd_width))
#             apd_off = np.zeros(self._pulser.waveform_padding(clk_rate * (seq_len - (i + apd_width))))
#             apd_do_waveform = np.concatenate((apd_wait_time, apd_on, apd_off))
#
#             if mw_pulse == True:
#                 mw_wait = np.zeros(self._pulser.waveform_padding(clk_rate * (laser_in + mw_wait_time)))
#                 mw_pulse_array = np.ones(self._pulser.waveform_padding(clk_rate * mw_pulse_time))
#                 mw_off = np.zeros(
#                     self._pulser.waveform_padding(clk_rate * seq_len - (laser_in + mw_wait_time + mw_pulse_time)))
#                 mw_pulse_array = np.subtract(mw_pulse_array, 0.5)  # because the MW takes only +-0.5V
#                 mw_ao_waveform = np.concatenate((mw_wait, mw_pulse_array, mw_off))
#
#                 # Make sure the arrays have the same length:
#                 # Check which array is the longest and append the two others:
#                 check_out = self.check_len(laser_do_waveform,apd_do_waveform,mw_ao_waveform, mw_pulse= True)
#                 # print ('check_out: ', check_out)
#                 # order: mw_ao_waveform, apd_do_waveform, laser_do_waveform, d_off
#                 mw_ao_waveform = check_out[0]
#                 apd_do_waveform = check_out[1]
#                 laser_do_waveform = check_out[2]
#                 d_off = check_out[3]
#                 ao_dict = {mw_outchan: mw_ao_waveform}
#
#             else:  # if no MW is needed
#                 # check which array is longer laser or apd
#                 check_out = self.check_len(laser_do_waveform, apd_do_waveform, mw_ao_waveform=[], mw_pulse=False)
#                 # print('check_out: ', check_out)
#                 # order: mw_ao_waveform, apd_do_waveform, laser_do_waveform, d_off
#                 mw_ao_waveform = check_out[0]
#                 apd_do_waveform = check_out[1]
#                 laser_do_waveform = check_out[2]
#                 d_off = check_out[3]
#                 # Creates a dict for the analog channel if no MW is needed
#                 ao_dict = {mw_outchan: d_off}
#
#             # Now all the waveformes are prepared and have the same length
#             do_dict = {laser_outchan: laser_do_waveform, apd_outchan: apd_do_waveform}
#             dosequence.append((do_dict))
#             aosequence.append((ao_dict))
#
#         digital_output_map = {0: [laser_outchan,
#                                   apd_outchan]}  # directs to the right output channel analog channel0 directs to digital x0 and x1
#
#         for aos, dos in zip(aosequence, dosequence):
#             # Loads waveforms to the awg
#             self._pulser.load_waveform(ao_waveform_dictionary=aos, do_waveform_dictionary=dos,
#                                        digital_output_map=digital_output_map)
#
#         if trigger == True:
#             # This sequence waits for a software trigger to start playing and moving to the next step.
#             self._pulser.configure_ORmask(card_idx, None)
#             self._pulser.configure_ANDmask(card_idx, None)
#             array = np.array([0x40000000], dtype=np.int64)
#             stop_condition_list = np.repeat(array, step_count)
#             loops = np.ones((len(aosequence)), dtype=np.int64)
#
#         else:
#             # This sequence immediately starts after the sequences are loaded
#             self._pulser.configure_ORmask(card_idx, 'immediate')
#             self._pulser.configure_ANDmask(card_idx, None)
#             loop_array = np.array([rep], dtype=np.int64)
#             loops = np.repeat(loop_array, step_count)
#             stop_condition_list = np.array([])
#
#             # makes the waveforms repeat in loops
#         self._pulser.load_sequence(  # digital_sequences=dosequence,
#             loops_list=loops, segment_map=segment_map, stop_condition_list=stop_condition_list)
#         self._pulser.start_card(card_idx)
#         self._pulser.arm_trigger(card_idx)
#
#
# ###############     Test sequenece for testing single channels  ###############
#
#     def test_sequences(self, seq_len, clk_rate=100, rep=100, laser=bool, mw=bool, apd=bool, trigger=bool):
#         '''
#         This is just for testing the different channels. Just sends ones ore zeros to the 3 different channels
#         clk_rate:   something like 100
#         msplay:     time of one sequence in ms
#         trigger:    True if the software trigger should be used and False if the loop array should define the length
#         rep:        If there is no software trigger, rep gives the number of repetitions the sequence is played
#         '''
#
#         # Sofar we just use card1
#         card_idx = 1
#         seq_len *= 1e-3
#
#         if clk_rate == 100:
#             clk_rate = MEGA(100)
#             # print('clk_rate not changed')
#             print('clk_rate in 1/s: ', clk_rate)
#         else:
#             clk_rate = MEGA(clk_rate)
#             # print('clk_rate changed')
#             print('clk_rate in 1/s: ', clk_rate)
#
#         # This one sets up the right clk_rate
#         self.set_up_pulser(card_idx, clk_rate)
#         print('setup of the pulser is done')
#
#         # Set up the right channels:
#         # Microwave: analog channel 2 and 3 for I and Q parameter. For now we just need I
#         # Laser: digital channel X0
#         # Photon counter: digital channel X1 and X2 (for reference)
#         mw_I = 2  # analog     A0 for card 1   Input at MW generator is +- 0.5 V
#         mw_Q = 3  # analog     A1 for card 1   --> add later
#         laser_channel = 0  # digital    X0 Input of the laser is 0-1 V
#         apd_channel1 = 1  # digital    X1
#         channels = np.array([laser_channel, mw_I, apd_channel1])    #[0,2,1]
#
#         if rep == 100:
#             print('set rep to 100')
#         else:
#             print('set rep to something else')
#
#         segment_map = []
#
#         # Create arrays:
#         # For Laser
#         if laser == True:
#             laser_sequence = np.ones(self._pulser.waveform_padding(clk_rate * seq_len))
#         else:
#             laser_sequence = np.zeros(self._pulser.waveform_padding(clk_rate * seq_len))
#
#         # For MW
#         if mw == True:
#             mw_sequence = np.multiply(np.ones(self._pulser.waveform_padding(clk_rate * seq_len)), 0.5)
#         else:
#             mw_sequence = np.zeros(self._pulser.waveform_padding(clk_rate * seq_len))
#
#         # For APD
#         if apd == True:
#             apd_sequence = np.ones(self._pulser.waveform_padding(clk_rate * seq_len))*3
#         else:
#             apd_sequence = np.zeros(self._pulser.waveform_padding(clk_rate * seq_len))
#
#         #create some zeros
#         zero = np.zeros(self._pulser.waveform_padding(clk_rate * seq_len))
#
#         aosequence = [{channels[1]: mw_sequence}, {channels[0]: zero}]
#         dosequence = [{channels[0]: laser_sequence, channels[2]: apd_sequence},
#                       {channels[0]: zero, channels[2]: zero}]
#
#         # Play it 2 times because it doesn't work with only one waveform
#         digital_output_map = {0: [channels[0]], 1: [channels[2]]}  # directs to the right output channel, analog channel0 directs to digital x0 and x1
#         # print('digital output map: ', digital_output_map)
#         for aos, dos in zip(aosequence, dosequence):
#             # Loads waveforms to the awg
#             self._pulser.load_waveform(ao_waveform_dictionary=aos, do_waveform_dictionary=dos,
#                                        digital_output_map=digital_output_map)
#         if trigger == True:
#             # This sequence waits for a software trigger to start playing and moving to the next step.
#             self._pulser.configure_ORmask(card_idx, None)
#             self._pulser.configure_ANDmask(card_idx, None)
#             array = np.array([0x40000000], dtype=np.int64)
#             loops = np.ones((len(aosequence)), dtype=np.int64)
#             stop_condition_list = np.repeat(array,len(aosequence))
#             # print('loops: ', loops)
#             # print('stop con list: ', stop_condition_list)
#
#         else:
#         # This sequence immediately starts after the sequences are loaded
#             self._pulser.configure_ORmask(card_idx, 'immediate')
#             self._pulser.configure_ANDmask(card_idx, None)
#             array = np.array([0x40000000], dtype=np.int64)
#             stop_condition_list = np.repeat(array, len(aosequence))
#             loops= np.array([100,100], dtype=np.int64)
#             # print ('loops: ', loops)
#
#         # makes the waveforms repeat in loops
#         self._pulser.load_sequence(  # digital_sequences=dosequence,
#             loops_list=loops, segment_map=segment_map, stop_condition_list=stop_condition_list)
#         self._pulser.start_card(card_idx)
#         self._pulser.arm_trigger(card_idx)
#
#
# ###############     Old Stuff      ###############
#
#     def check_len(self, laser_do_waveform, apd_do_waveform, mw_ao_waveform, mw_pulse=bool):
#         # checks and adjusts the length of the different waveform because waveform_padding creates a mess
#         if mw_pulse == True:
#             val_list = [len(laser_do_waveform), len(apd_do_waveform), len(mw_ao_waveform)]
#             max_val = max(val_list)
#             # print ('max_val: ', max_value)
#             max_index = val_list.index(max_val)
#             print(max_index)
#             if max_index == 0:
#                 print('laser waveforms is the longest')
#                 difference_mw = len(laser_do_waveform) - len(mw_ao_waveform)
#                 difference_apd = len(laser_do_waveform) - len(apd_do_waveform)
#                 rest_array_mw = np.zeros(int(difference_mw))
#                 rest_array_apd = np.zeros(int(difference_apd))
#                 mw_ao_waveform = np.concatenate((mw_ao_waveform, rest_array_mw))
#                 apd_do_waveform = np.concatenate((apd_do_waveform, rest_array_apd))
#                 d_off = []
#             elif max_index == 1:
#                 print('apd waveforms is the longest')
#                 difference_mw = len(apd_do_waveform) - len(mw_ao_waveform)
#                 difference_laser = len(apd_do_waveform) - len(laser_do_waveform)
#                 rest_array_mw = np.zeros(int(difference_mw))
#                 rest_array_laser = np.zeros(int(difference_laser))
#                 mw_ao_waveform = np.concatenate((mw_ao_waveform, rest_array_mw))
#                 laser_do_waveform = np.concatenate((laser_do_waveform, rest_array_laser))
#                 d_off = []
#             else:
#                 print('mw waveforms is the longest')
#                 difference_apd = len(mw_ao_waveform) - len(apd_do_waveform)
#                 difference_laser = len(mw_ao_waveform) - len(laser_do_waveform)
#                 rest_array_apd = np.zeros(int(difference_apd))
#                 rest_array_laser = np.zeros(int(difference_laser))
#                 apd_do_waveform = np.concatenate((apd_do_waveform, rest_array_apd))
#                 laser_do_waveform = np.concatenate((laser_do_waveform, rest_array_laser))
#                 d_off = []
#         else:
#             if len(apd_do_waveform) < len(laser_do_waveform):
#                 print('laser waveform is longer than apd waveform')
#                 difference = len(laser_do_waveform) - len(apd_do_waveform)
#                 rest_array = np.zeros(int(difference))
#                 apd_do_waveform = np.concatenate((apd_do_waveform, rest_array))
#                 # Create zero array of same size:
#                 d_off = np.zeros(len(laser_do_waveform))
#             elif len(laser_do_waveform) < len(apd_do_waveform):
#                 print('apd waveform is the longer')
#                 difference = len(apd_do_waveform) - len(laser_do_waveform)
#                 rest_array = np.zeros(int(difference))
#                 laser_do_waveform = np.concatenate((laser_do_waveform, rest_array))
#                 # Create zero array of same size:
#                 d_off = np.zeros(len(apd_do_waveform))
#             else:
#                 print('apd_waveform and laser_waveform have the same size')
#                 d_off = np.zeros(len(apd_do_waveform))
#         return mw_ao_waveform, apd_do_waveform, laser_do_waveform, d_off
#
#
#
#     # This is the main block...not that general but works with the laser
#     def delay_sweep_mw(self, clk_rate, laser_in, wait, laser_re, apd_width, rest, start, steps, mw_wait_time,
#                        mw_pulse_time, rep, mw_pulse=bool, trigger=bool, segment_map=[]):
#
#         '''
#         check for the longest trace for sure
#         This can play a laser sequence and the apd readout at the same time
#         clk_rate:   something like 100
#         laser_in:   length of initialisation laser pulse (it starts right at the beginning) in ms
#         wait:       time between the end of the initialisation and the reinitialisation in ms
#         laser_re:   length of the reinitialisation pulse in ms
#
#         apd_width:   time in ms where the apd is on (this needs to be varied to do the delay measurement)
#         rest:       time in ms after the reinitialisation pulse to have a break before the next cycle starts
#         start:      where in time should the apd pulse start shifting
#         steps:      step size of the apd pulse sweep
#         mw_wait_time: when should the mw pulse should happen?
#         mw_pulse_time: length of the mw pulse in between the two laser pulses
#         trigger:    True if the software trigger should be used and False if the loop array should define the length
#         rep:        If there is no software trigger, rep gives the number of repetitions the waveform is played
#         segment_map:If the order of different segments should be mixed (not there yet)
#         '''
#         # building the array to change the length of the apd pulse
#         msplay = laser_in + wait + laser_re + rest
#         clk_rate = MEGA(clk_rate)
#         stop = (msplay - rest) - apd_width
#         # print('stop value: ', stop)
#
#         # build array with all the desired apd_start values]
#         apd_start = []
#         step_count = math.floor((stop - start) / steps)
#         # print('step_count: ', step_count)
#         for i in range(step_count):
#             apd_start.append((start + (i * steps)))
#         # print('adpstart array: ', apd_start)
#         laser_in *= 1e-3
#         mw_wait_time *= 1e-3
#         mw_pulse_time *= 1e-3
#         print('mw specs: ', mw_wait_time)
#         print('mw specs: ', mw_pulse_time)
#         laser_re *= 1e-3
#         wait *= 1e-3
#         apd_width *= 1e-3
#         rest *= 1e-3
#         start *= 1e-3
#         steps *= 1e-3
#         msplay *= 1e-3
#         apd_start = np.multiply(1e-3, apd_start)
#         # print('apd_start: ', apd_start)
#
#         # Laser waveform does not change the whole time
#         first_pulse = np.ones(self._pulser.waveform_padding(clk_rate * laser_in))
#         laser_wait_time = np.zeros(self._pulser.waveform_padding(clk_rate * wait))
#         sec_pulse = np.ones(self._pulser.waveform_padding(clk_rate * laser_re))
#         rest_array = np.zeros(self._pulser.waveform_padding(clk_rate * rest))
#         laser_do_waveform = np.concatenate((first_pulse, laser_wait_time, sec_pulse, rest_array))
#         dosequence = []
#         aosequence = []
#         laser_outchan = 0  # Laser
#         apd_outchan = 1  # APD
#         card_idx = 1
#         for i in apd_start:
#             apd_wait_time = np.zeros(self._pulser.waveform_padding(clk_rate * i))
#             apd_on = np.ones(self._pulser.waveform_padding(clk_rate * apd_width))
#             apd_off = np.zeros(self._pulser.waveform_padding(clk_rate * (msplay - (i + apd_width))))
#             apd_do_waveform = np.concatenate((apd_wait_time, apd_on, apd_off))
#
#             if mw_pulse == True:
#                 # if a mw pulse is desired:
#                 mw_wait = np.zeros(self._pulser.waveform_padding(clk_rate * mw_wait_time))
#                 pi_pulse_array = np.ones(self._pulser.waveform_padding(clk_rate * mw_pulse_time))
#                 mw_off = np.zeros(self._pulser.waveform_padding(clk_rate * msplay - (mw_wait_time + mw_pulse_time)))
#                 mw_ao_waveform = np.concatenate((mw_wait, pi_pulse_array, mw_off))
#                 mw_ao_waveform = np.subtract(mw_ao_waveform, 0.5)  # because the MW takes only +-0.5V
#                 print('mw_ao_waveform max: ', np.max(mw_ao_waveform))
#                 # print('len of mw_ao_waveform: ', len(mw_ao_waveform))
#
#                 # Make sure the arrays have the same length:
#                 # Check which array is the longest and append the two others:
#                 max_val = np.max(np.array([len(laser_do_waveform), len(apd_do_waveform), len(mw_ao_waveform)]))
#                 if max_val == 0:
#                     print('laser waveforms is the longest')
#                     difference_mw = len(laser_do_waveform) - len(mw_ao_waveform)
#                     difference_apd = len(laser_do_waveform) - len(apd_do_waveform)
#                     rest_array_mw = np.zeros(int(difference_mw))
#                     rest_array_apd = np.zeros(int(difference_apd))
#                     mw_ao_waveform = np.concatenate((mw_ao_waveform, rest_array_mw))
#                     apd_do_waveform = np.concatenate((apd_do_waveform, rest_array_apd))
#                     print('len of mw_ao_waveform: ', len(mw_ao_waveform))
#                     print('len of laser_do_waveform: ', len(laser_do_waveform))
#                     print('len of apd_do_waveform: ', len(apd_do_waveform))
#                 elif max_val == 1:
#                     print('apd waveforms is the longest')
#                     difference_mw = len(apd_do_waveform) - len(mw_ao_waveform)
#                     difference_laser = len(apd_do_waveform) - len(laser_do_waveform)
#                     rest_array_mw = np.zeros(int(difference_mw))
#                     rest_array_laser = np.zeros(int(difference_laser))
#                     mw_ao_waveform = np.concatenate((mw_ao_waveform, rest_array_mw))
#                     laser_do_waveform = np.concatenate((laser_do_waveform, rest_array_laser))
#                     print('len of mw_ao_waveform: ', len(mw_ao_waveform))
#                     print('len of laser_do_waveform: ', len(laser_do_waveform))
#                     print('len of apd_do_waveform: ', len(apd_do_waveform))
#                 else:
#                     print('mw waveforms is the longest')
#                     difference_apd = len(mw_ao_waveform) - len(apd_do_waveform)
#                     difference_laser = len(mw_ao_waveform) - len(laser_do_waveform)
#                     rest_array_apd = np.zeros(int(difference_apd))
#                     rest_array_laser = np.zeros(int(difference_laser))
#                     apd_do_waveform = np.concatenate((apd_do_waveform, rest_array_apd))
#                     laser_do_waveform = np.concatenate((laser_do_waveform, rest_array_laser))
#                     print('len of mw_ao_waveform: ', len(mw_ao_waveform))
#                     print('len of laser_do_waveform: ', len(laser_do_waveform))
#                     print('len of apd_do_waveform: ', len(apd_do_waveform))
#
#                 ao_dict = {2: mw_ao_waveform}
#
#             else:  # if no MW is needed
#                 # check which array is longer from laser and apd and use that one
#                 if len(apd_do_waveform) < len(laser_do_waveform):
#                     print('laser waveform is the longer')
#                     difference = len(laser_do_waveform) - len(apd_do_waveform)
#                     rest_array = np.zeros(int(difference))
#                     apd_do_waveform = np.concatenate((apd_do_waveform, rest_array))
#                     # Create zero array of same size:
#                     d_off = np.zeros(len(laser_do_waveform))
#                     print('len of mw_ao_waveform: ', len(d_off))
#                     print('len of laser_do_waveform: ', len(laser_do_waveform))
#                     print('len of apd_do_waveform: ', len(apd_do_waveform))
#                 elif len(laser_do_waveform) < len(apd_do_waveform):
#                     print('apd waveform is the longer')
#                     difference = len(apd_do_waveform) - len(laser_do_waveform)
#                     rest_array = np.zeros(int(difference))
#                     laser_do_waveform = np.concatenate((laser_do_waveform, rest_array))
#                     # Create zero array of same size:
#                     d_off = np.zeros(len(apd_do_waveform))
#                     print('len of mw_ao_waveform: ', len(d_off))
#                     print('len of laser_do_waveform: ', len(laser_do_waveform))
#                     print('len of apd_do_waveform: ', len(apd_do_waveform))
#                 else:
#                     print('apd waveform and laser_waveform have the same size')
#                     d_off = np.zeros(len(apd_do_waveform))
#                     print('len of mw_ao_waveform: ', len(d_off))
#                     print('len of laser_do_waveform: ', len(laser_do_waveform))
#                     print('len of apd_do_waveform: ', len(apd_do_waveform))
#                 # Creates a dict for the analog channel if no MW is needed
#                 ao_dict = {2: d_off}
#
#             # Now all the waveformes are prepared and have the same length
#             do_dict = {laser_outchan: laser_do_waveform, apd_outchan: apd_do_waveform}
#             dosequence.append((do_dict))
#             aosequence.append((ao_dict))
#
#         # print('dosequence: ', dosequence)
#         # print('aosequence: ', aosequence)
#         # load and play
#         # print('len of dosequence: ', len(dosequence))
#         # print('len of aosequence: ', len(aosequence))
#         digital_output_map = {
#             0: [0, 1]}  # directs to the right output channel analog channel0 directs to digital x0 and x1
#         for aos, dos in zip(aosequence, dosequence):
#             # Loads waveforms to the awg
#             self._pulser.load_waveform(ao_waveform_dictionary=aos, do_waveform_dictionary=dos,
#                                        digital_output_map=digital_output_map)
#
#         if trigger == True:
#             # This sequence waits for a software trigger to start playing and moving to the next step.
#             self._pulser.configure_ORmask(card_idx, None)
#             self._pulser.configure_ANDmask(card_idx, None)
#             array = np.array([0x40000000], dtype=np.int64)
#             stop_condition_list = np.repeat(array, step_count)
#             # print('stop_condition_list: ', stop_condition_list)
#             loops = np.ones((len(aosequence)), dtype=np.int64)
#
#         else:
#             # This sequence immediately starts after the sequences are loaded
#             self._pulser.configure_ORmask(card_idx, 'immediate')
#             self._pulser.configure_ANDmask(card_idx, None)
#             loop_array = np.array([rep], dtype=np.int64)
#             loops = np.repeat(loop_array, step_count)
#             print('loops: ', loops)
#             stop_condition_list = np.array([])
#
#         # makes the waveforms repeat in loops
#         self._pulser.load_sequence(  # digital_sequences=dosequence,
#             loops_list=loops, segment_map=segment_map, stop_condition_list=stop_condition_list)
#         self._pulser.start_card(card_idx)
#         self._pulser.arm_trigger(card_idx)
#
#     def mw_test(self, msecondsplay, loops, segment_map=[]):
#         """
#         Function for early debugging of the awg. Remove from the final class.
#
#         The steps to a successful initialization of a sequence are:
#         - Define the clock rate
#         - create the arrays that will represent the sequences
#         - configure the trigger masks
#         - load the sequence
#         - start the card
#         - arm the trigger
#         """
#         msecondsplay *= 1e-3
#
#         clk_rate = MEGA(100)
#
#         card_idx = 1
#         self.set_sample_rate(card_idx, int(clk_rate))
#
#         samples = self.waveform_padding(msecondsplay * clk_rate)
#         time_ax = np.linspace(0, samples / clk_rate, samples)
#         first_seq_ch0 = np.sin(2 * np.pi * time_ax / msecondsplay)
#         second_seq_ch0 = np.sin(3 * 2 * np.pi * time_ax / msecondsplay)
#         first_seq_ch1 = np.linspace(0, 1, len(time_ax))
#         # func_1= 0.5 * (signal.square(2 * np.pi * 20 * time_ax)) + 0.5
#         #
#         sigma = (time_ax[-1]) / 20
#         second_seq_ch1 = np.exp(-(time_ax - time_ax[-1] / 2) ** 2 / (2 * sigma ** 2)) * np.sin(
#             20 * 2 * np.pi * time_ax / msecondsplay)
#
#         aosequence = [{2: first_seq_ch0, 3: first_seq_ch1},
#                       {2: first_seq_ch1, 3: first_seq_ch0},
#                       {2: first_seq_ch0, 3: first_seq_ch1},
#                       {2: first_seq_ch1, 3: first_seq_ch0}
#                       ]
#
#         do_chan = 1
#         do1 = np.zeros(first_seq_ch0.shape)
#         do2 = np.zeros(second_seq_ch0.shape)
#         do3 = np.zeros(first_seq_ch1.shape)
#         do4 = np.copy(do3)
#         do1[first_seq_ch0 > 0] = 1
#         do2[second_seq_ch0 > 0] = 1
#         do3[:] = 1
#         print('len do1', len(do1))
#         print('len do2', len(do2))
#         print('len do3', len(do3))
#         print('len do4', len(do4))
#         outchan = 0
#         dosequence = [{outchan: do1},
#                       {outchan: do2},
#                       {outchan: do3},
#                       {outchan: do4}]
#
#         digital_output_map = {1: [outchan]}
#         for aos, dos in zip(aosequence, dosequence):
#             self.load_waveform(ao_waveform_dictionary=aos, do_waveform_dictionary=dos,
#                                digital_output_map=digital_output_map)
#
#         # This sequence immediately starts after the sequences are loaded
#         # self.configure_ORmask(card_idx, 'immediate')
#         # self.configure_ANDmask(card_idx, None)
#         # stop_condition_list = np.array([])
#
#         # This sequence waits for a software trigger to start playing and moving to the next step.
#         self.configure_ORmask(card_idx, None)
#         self.configure_ANDmask(card_idx, None)
#         loops = np.ones((len(aosequence)), dtype=np.int64)
#         # print('loops: ', loops)
#         stop_condition_list = np.array(
#             [SPCSEQ_ENDLOOPONTRIG, SPCSEQ_ENDLOOPONTRIG, SPCSEQ_ENDLOOPONTRIG, SPCSEQ_ENDLOOPONTRIG], dtype=np.int64)
#         # print('segment map: ', segment_map)
#         # print('stop cond.:', stop_condition_list)
#
#         self.load_sequence(  # digital_sequences=dosequence,
#             loops_list=loops, segment_map=segment_map, stop_condition_list=stop_condition_list)
#
#         self.start_card(card_idx)
#         self.arm_trigger(card_idx)

    # def calc_len_error(self, t_centerwidths, seq_len):
    #     # print(len(t_centerwidths))
    #     # print(t_centerwidths[0][1])
    #     real_lenghts = []
    #     for i in range(len(t_centerwidths)):
    #         if (t_centerwidths[i][1] % 2) == 0: # even
    #             if t_centerwidths[i][1] == 0:
    #                 pulse_len = t_centerwidths[i][1]
    #             else:
    #                 pulse_len = t_centerwidths[i][1] + 1
    #         else: # odd
    #             pulse_len = t_centerwidths[i][1]
    #         real_lenghts.append(pulse_len)
    #     return real_lenghts

