# -*- coding: utf-8 -*-
"""
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

import time
import numpy as np
from scipy.interpolate import interp1d

from logic.generic_logic import GenericLogic
from core.configoption import ConfigOption
from core.connector import Connector
from core.statusvariable import StatusVar
from interface.pulser_interface import PulserInterface
from interface.microwave_interface import MicrowaveInterface


class IQPulserInterfuse(GenericLogic, MicrowaveInterface, PulserInterface):
    """
    An interfuse that merges together the functionality of AWG and IQ-modulated microwave source, in order to
    facilitate single sideband upconversion scheme.
    """

    mw_source = Connector(interface="MicrowaveInterface")
    awg = Connector(interface="PulserInterface")

    if_modulation_freq = StatusVar(default=100e6)

    calibfile_dir = ConfigOption("calibration_file_directory", missing="error")
    calibration_interpolation_method = ConfigOption(
        "calibration_file_interpolation_method", "linear", missing="warn"
    )
    laser_channel = ConfigOption("laser_channel", "x0", missing="warn")
    apd_signal_channel = ConfigOption("apd_signal_channel", "x1", missing="warn")
    apd_read_channel = ConfigOption("apd_read_channel", "x2", missing="warn")
    imod_channel = ConfigOption("imod_channel", "ch2", missing="warn")
    qmod_channel = ConfigOption("qmod_channel", "ch3", missing="warn")

    def on_activate(self):
        self._if_freq = self.if_modulation_freq
        self._laser = self.laser_channel
        self._apdsig = self.apd_signal_channel
        self._apdread = self.apd_read_channel
        self._ichan = self.imod_channel
        self._qchan = self.qmod_channel

        # Make a dictionary to map the channel names to the awg channel names.
        self._channel_dictionary = {
            "laser": int(self.laser_channel[-1]),
            "apd_sig": int(self.apd_signal_channel[-1]),
            "apd_read": int(self.apd_read_channel[-1]),
            "i_chan": int(self.imod_channel[-1]),
            "q_chan": int(self.qmod_channel[-1]),
        }

        # This dictionary builds up the map between user-friendly names and functions used to generate
        # the pulse envelopes
        self._envelope_dictionary = {"square": self.box_envelope}

        # The calibration file should be order in columns of:
        # frequency | i amplitude | i offset | q amplitude | q offset | phase_imbalance
        try:
            self._calibration = np.loadtxt(self.calibfile_dir)
            frequency, calibration_values = (
                self._calibration[:, 0],
                self._calibration[:, 1:],
            )
        except:
            self.log.warning(
                "Calibration file is empty or contains only one value. Using default offsets and amplitudes."
            )
            frequency = np.array([2.87e9, 2.88e9])
            calibration_values = np.zeros((2, 5))
            calibration_values[:, [0, 2]] = 1.0

        # Store the minimum and maximum frequency of the calibration file for warning purposes
        self._freqminmax = frequency[[0, -1]]

        # This function spits out the I amplitude, I offset, Q amplitude, Q offset and phase imbalance in this order
        # for a given frequency. If a requested frequency falls outside the range of frequencies in the calibration file,
        # the function will return the calibration parameters of the closest frequency in the calibration file.
        self._iq_interpfunc = interp1d(
            frequency,
            calibration_values,
            kind=self.calibration_interpolation_method,
            axis=0,
            bounds_error=False,
            fill_value=(calibration_values[0], calibration_values[-1]),
        )

        self._mwsource = self.mw_source()
        self._awg = self.awg()

    def on_deactivate(self):
        pass

    def set_frequency(self, freq=0.0):
        """
        Sets the frequency of the mw source. The actual value is set to be freq - if modulation frequency.
        @param float freq: The desired output frequency.
        @return:
        """
        lo_frequency = freq - self._if_freq
        self._mwsource.set_frequency(lo_frequency)

    def load_waveform(
        self, iq_dictionary=dict(), digital_pulses=dict(), digital_output_map=dict()
    ):
        """
        Method overloading the awg load_waveform.
        @param dict iq_dictionary: the dictionary for i and q channels. Valid keys are i_chan and q_chan.
        @param dict digital_pulses: the dictionary for the digital pulses. Valid keys are laser, apd_sig, apd_ref.
        @param dict digital_output_map: the digital output map, as used in the overloaded method.

        @return int: Error code. 0: ok, -1: error occurred
        """
        try:
            ao_waveform_dictionary = {
                self._channel_dictionary[chan_name]: waveform
                for chan_name, waveform in iq_dictionary.items()
            }
        except:
            self.log.error("The requested analog channel name is not valid.")
            return -1

        try:
            do_waveform_dictionary = {
                self._channel_dictionary[chan_name]: waveform
                for chan_name, waveform in digital_pulses.items()
            }
        except:
            self.log.error("The requested digital channel name is not valid.")

        self._awg.load_waveform(
            ao_waveform_dictionary=ao_waveform_dictionary,
            do_waveform_dictionary=do_waveform_dictionary,
            digital_output_map=digital_output_map,
        )

        return 0

    def iq_pulses(
        self,
        timeaxis,
        t_centrewidth_list,
        output_frequency,
        phases=None,
        pulsenvelope="square",
    ):
        """
        Method used to generate the iq modulation pulses.

        @param np.ndarray timeaxis: the time axis.
        @param list t_centrewidth_list: a list of matrices containing the central times and widths of the pulses, as
                                        given by the envelope functions. Each element of the list corresponds to a
                                        waveform with a given phase.
        @param output_frequency: the desired frequency of the pulses (i.e. local oscillator + if modulation).
        @param np.ndarray phases: The array containing the phases
        @param str pulsenvelope:
        @return tuple: the I and Q waveforms.
        """
        envelope_function = self._envelope_dictionary[pulsenvelope]

        if phases is None:
            phases = np.zeros((len(t_centrewidth_list),))

        phases_number = len(phases)

        lo_frequency = output_frequency - self._if_freq
        if not self._freqminmax[0] < lo_frequency < self._freqminmax[-1]:
            self.log.warning(
                "Requested local oscillator frequency: {:.3f} GHz, falls outside the calibration file ranges. "
                "Using the closest frequency instead.".format(lo_frequency / 1e9)
            )

        # Get the calibrated amplitudes and offset from the interpolation function.
        Iamp, Ioff, Qamp, Qoff, phase_inbalance = self._iq_interpfunc(lo_frequency)

        # Create a matrix with as many rows (in the second index) as the number of requested phases. Essentially,
        # create sinusoidal oscillations (for I and Q) at each phase for the whole time axis duration, then chop
        # them up with the pulse sequence. At the final step the pulses at different phases are summed up,
        # to combine into two waveforms, one for I and one for Q.
        pulse_matrix = np.zeros((2, phases_number, timeaxis.size))

        for phase_idx, phase, t_centrewidth in zip(
            range(phases_number), phases, t_centrewidth_list
        ):
            Imodulation = (
                Iamp * np.cos(2 * np.pi * self.if_modulation_freq * timeaxis + phase)
                + Ioff
            )
            Qmodulation = (
                Qamp
                * np.cos(
                    2 * np.pi * self.if_modulation_freq * timeaxis
                    + phase
                    + phase_inbalance
                    + np.pi / 2
                )
                + Qoff
            )

            pulse_sequence = envelope_function(timeaxis, t_centrewidth)

            Icomp = pulse_sequence * Imodulation
            Qcomp = pulse_sequence * Qmodulation

            pulse_matrix[0, phase_idx] = Icomp
            pulse_matrix[1, phase_idx] = Qcomp

        Ipulses, Qpulses = pulse_matrix.sum(axis=1)

        return Ipulses, Qpulses

    @staticmethod
    def box_envelope(timeaxis, t_centrewidths):
        """
        Method used to generate a square envelope.

        @param np.ndarray timeaxis: the time axis
        @param np.ndarray t_centrewidths: 2D matrix. Each row contains the central time and width of the square pulse.
        @return np.ndarray: the pulse train.
        """

        tstart = t_centrewidths[:, 0] - t_centrewidths[:, 1] / 2
        tend = t_centrewidths[:, 0] + t_centrewidths[:, 1] / 2
        start_edges = timeaxis[None, :] >= tstart[:, None]
        end_edges = timeaxis[None, :] <= tend[:, None]

        pulses = np.logical_and(start_edges, end_edges).astype(np.float64)

        return pulses.sum(axis=0)

    ################################################################################################
    # More custom methods that are not related to the original MicrowaveInterface and PulserInterface
    ################################################################################################
    def set_chan_amplitude(self, channels, amplitudes):
        err_code = self._awg.set_chan_amplitude(channels, amplitudes)
        return err_code

    def set_output_filters(self, channels, filter_active):
        err_code = self._awg.set_output_filters(channels, filter_active)
        return err_code
    def configure_ORmask(self, card_idx, *masks_to_enable):
        err_code = self._awg.configure_ORmask(card_idx, *masks_to_enable)

    def configure_ANDmask(self, card_idx, *masks_to_enable):
        err_code = self._awg.configure_ANDmask(card_idx, *masks_to_enable)
        return err_code

    def send_software_trig(self, card_idx):
        self._awg.send_software_trig(card_idx)

    def start_card(self, card_idx):
        self._awg.start_card(card_idx)

    def arm_trigger(self, card_idx):
        self._awg.arm_trigger(card_idx)

    def disable_trigger(self, card_idx):
        self._awg.disable_trigger(card_idx)

    def configuration_sequence_mode(self, card_idx):
        self._awg.configuration_sequence_mode(card_idx)

    def configuration_single_mode(self, card_idx):
        self._awg.configuration_single_mode(card_idx)

    def stop_replay(self, card_idx):
        self._awg.stop_replay(card_idx)

    def waveform_padding(self, waveform_len):
        return self._awg.waveform_padding(waveform_len)

    def set_IQmod(self, on):
        self._mwsource.set_IQmod(on)

    ################################################################################################
    # Overloading section for the MicrowaveInterface
    ################################################################################################
    def off(self):
        return self._mwsource.off()

    def get_status(self):
        return self._mwsource.get_status()

    def get_power(self):
        return self._mwsource.get_power()

    def get_frequency(self):
        return self._mwsource.get_frequency()

    def cw_on(self):
        return self._mwsource.cw_on()

    def set_cw(self, frequency=None, power=None, useinterleave=None):
        return self._mwsource.set_cw(
            frequency=frequency, power=power, useinterleave=useinterleave
        )

    def list_on(self):
        return self._mwsource.list_on()

    def set_list(self, frequency=None, power=None):
        return self._mwsource.set_list(frequency=frequency, power=power)

    def reset_listpos(self):
        return self._mwsource.reset_listpos()

    def sweep_on(self):
        return self._mwsource.sweep_on()

    def set_sweep(self, start=None, stop=None, step=None, power=None):
        return self._mwsource.set_sweep(start=start, stop=stop, step=step, power=power)

    def reset_sweeppos(self):
        return self._mwsource.reset_sweeppos()

    def set_ext_trigger(self, pol, timing):
        return self._mwsource.set_ext_trigger(pol, timing)

    def trigger(self):
        return self._mwsource.trigger()

    def get_limits(self):
        return self._mwsource.get_limits()

    ################################################################################################
    # Overloading section for the PulserInterface
    ################################################################################################
    def get_constraints(self):
        return self._awg.get_constraints()

    def pulser_on(self):
        self._awg.pulser_on()

    def pulser_off(self):
        self._awg.pulser_off()

    def load_sequence(
        self,
        waveform_list=None,
        segment_map=np.array([]),
        loops_list=np.array([]),
        stop_condition_list=np.array([]),
    ):
        error_out = self._awg.load_sequence(
            waveform_list=None,
            segment_map=np.array([]),
            loops_list=np.array([]),
            stop_condition_list=np.array([]),
        )
        return error_out

    def get_loaded_assets(self):
        return self._awg.get_loaded_assets()

    def clear_all(self):
        return self._awg.clear_all()

    def get_status(self):
        return self._awg.get_status()

    def get_sample_rate(self, card_idx):
        return self._awg.get_sample_rate(card_idx)

    def set_sample_rate(self, card_idx, clk_rate):
        error_code = self._awg.set_sample_rate(card_idx, clk_rate)
        return error_code

    def get_analog_level(self):
        return self._awg.get_analog_level()

    def set_analog_level(self):
        self._awg.set_analog_level()

    def get_digital_level(self):
        return self._awg.get_digital_level()

    def set_digital_level(self):
        self._awg.set_digital_level()

    def get_active_channels(self, ch=None):
        return self._awg.get_active_channels()

    def set_active_channels(self, *channels):
        self._awg.set_active_channels(*channels)

    def write_waveform(self):
        self._awg.write_waveform()

    def write_sequence(self):
        self._awg.write_waveform()

    def get_waveform_names(self):
        return self._awg.get_waveform_names()

    def get_sequence_names(self):
        return self._awg.get_sequence_names()

    def delete_waveform(self):
        return self._awg.delete_waveform()

    def delete_sequence(self):
        return self._awg.delete_sequence

    def get_interleave(self):
        return self._awg.get_interleave()

    def set_interleave(self):
        self._awg.set_interleave()

    def reset(self):
        self._awg.reset()
