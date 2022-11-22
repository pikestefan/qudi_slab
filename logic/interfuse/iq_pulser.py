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

from core.module import Base
from core.configoption import ConfigOption
from core.connector import Connector
from core.statusvariable import StatusVar
from interface.pulser_interface import PulserInterface
from interface.microwave_interface import MicrowaveInterface
from scipy.interpolate import interp1d


class IQPulserInterfuse(Base, PulserInterface, MicrowaveInterface):
    """
    An interfuse that merges together the functionality of AWG and IQ-modulated microwave source, in order to
    facilitate single sideband upconversion scheme.
    """

    mw_source = Connector(interface="PulserInterface")
    awg = Connector(interface="MicrowaveInterface")

    if_modulation_freq = StatusVar(default=100e6)

    calibfile_dir = ConfigOption("calibration_file_directory", missing="error")
    calibration_interpolation_method = ConfigOption(
        "calibration_file_interpolation_method", "linear", missing="warning"
    )

    def on_activate(self):
        self._mwsource = self.mw_source
        self._awg = self.awg

        self._if_freq = self.if_modulation_freq

        # This dictionary builds up the map between user-friendly names and functions used to generate
        # the pulse envelopes
        self._envelope_dictionary = {"square": self.box_envelope}

        # The calibration file should be order in columns of:
        # frequency | i amplitude | i offset | q amplitude | q offset
        self._calibration = np.loadtxt(self.calibfile_dir)
        if self._calibration.size <= 1:
            self.log.warning(
                "Calibration file is empty or contains only one value. Using default offsets and amplitudes."
            )

            frequency = np.array([2.87e9, 2.88e9])
            calibration_values = np.zeros((2, 4))
            calibration_values[:, [0, 2]] = 1.0

        else:
            frequency, calibration_values = (
                self._calibration[:, 0],
                self._calibration[:, 1:],
            )

        # Store the minimum and maximum frequency of the calibration file for warning purposes
        self._freqminmax = frequency[[0, -1]]

        # This function spits out the I amplitude, I offset, Q amplitude, Q offset in this order for a given
        # frequency. If a requested frequency falls outside the range of frequencies in the calibration file, the
        # function will return the calibration parameters of the closest frequency in the calibration file.
        self._iq_interpfunc = interp1d(
            frequency,
            calibration_values,
            kind=self.calibration_interpolation_method,
            axis=0,
            bounds_error=False,
            fill_value=(calibration_values[0], calibration_values[-1]),
        )

    def set_frequency(self, freq=0.0):
        lo_frequency = freq - self._if_freq
        self._mwsource.set_frequency(lo_frequency)

    def iq_pulses(
        self,
        timeaxis,
        t_centrewidth_list,
        lo_frequency,
        phases=None,
        pulsenvelope="square",
    ):
        """
        Method used to generate the iq modulation pulses.

        @param np.ndarray timeaxis: the time axis.
        @param list t_centrewidth_list: a list of matrices containing the central times and widths of the pulses, as
                                        given by the envelope functions. Each element of the list corresponds to a
                                        waveform with a given phase.
        @param np.ndarray phases: The
        @param str pulsenvelope:
        @return: the I and Q waveforms.
        """
        envelope_function = self._envelope_dictionary[pulsenvelope]

        if phases is None:
            phases = np.zeros((len(t_centrewidth_list),))

        phases_number = len(phases)

        if not self._freqminmax[0] < lo_frequency < self._freqminmax[-1]:
            self.log.warning(
                "Requested local oscillator frequency: {:.3f} GHz, falls outside the calibration file ranges. "
                "Using the closest frequency instead."
            )

        # Get the calibrated amplitudes and offset from the interpolation function.
        Iamp, Ioff, Qamp, Qoff = self._iq_interpfunc(lo_frequency)

        # Create a matrix with as many rows (in the second index) as the number of requested phases. Essentially,
        # create sinusoidal oscillations (for I and Q) at each phase for the whole time axis duration, then chop
        # them up with the pulse sequence. At the final step the pulses at different phases are summed up,
        # to combine into two waveforms, one for I and one for Q.
        pulse_matrix = np.zeros((2, phases_number, timeaxis.size))

        for phase_idx, phase, t_centrewidth in zip(
            range(phases_number), phases, t_centrewidth_list
        ):
            Imodulation = (
                Iamp * np.cos(2 * np.pi * self.if_modulation_freq + phase) + Ioff
            )
            Qmodulation = (
                Qamp * np.cos(2 * np.pi * self.if_modulation_freq + phase + np.pi / 2)
                + Qoff
            )
            Icomp = envelope_function(timeaxis, t_centrewidth) * Imodulation
            Qcomp = envelope_function(timeaxis, t_centrewidth) * Qmodulation

            pulse_matrix[0, phase_idx] = Icomp
            pulse_matrix[1, phase_idx] = Qcomp

        Ipulses, Qpulses = pulse_matrix.sum(axis=1)

        return Ipulses, Qpulses

    def box_envelope(self, timeaxis, t_centrewidths):
        """
        Method used to generate a square envelope.

        @param timeaxis: the time axis
        @param t_centrewidths: 2D matrix. Each row contains the central time and width of the square pulse.
        @return: the pulse train.
        """

        tstart = t_centrewidths[:, 0] - t_centrewidths[:, 1] / 2
        tend = t_centrewidths[:, 0] + t_centrewidths[:, 1] / 2
        start_edges = timeaxis[None, :] >= tstart[:, None]
        end_edges = timeaxis[None, :] <= tend[:, None]

        pulses = np.logical_and(start_edges, end_edges).astype(np.float64)

        return pulses.sum(axis=0)
