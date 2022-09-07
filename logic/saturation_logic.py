# -*- coding: utf-8 -*-

"""
This file contains the Qudi Logic module base class.

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
from interface.microwave_interface import MicrowaveMode
from interface.microwave_interface import TriggerEdge
import numpy as np
import time
import datetime
import matplotlib.pyplot as plt

from logic.generic_logic import GenericLogic
from core.util.mutex import Mutex
from core.connector import Connector
from core.configoption import ConfigOption
from core.statusvariable import StatusVar


class SaturationLogic(GenericLogic):
    daq_card = Connector(interface='NationalInstrumentsXSeriesPxScan')
    laser = Connector(interface='NationalInstrumentsUSB6000x')
    plotlogic = Connector(interface='QDPlotLogic')

    # Status Vars
    max_laser_voltage = ConfigOption('max_laser_voltage', 0.95)
    laser_voltage = StatusVar('laser_voltage', 0)
    integration_time = StatusVar('intgration_time', default=30e-3)

    sigMeasuredSaturation = QtCore.Signal(int, np.ndarray, np.ndarray)

    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)
        self.threadlock = Mutex()

    def on_activate(self):
        self._daq_card = self.daq_card()
        self._laser = self.laser()
        self._plotlogic = self.plotlogic()

    def on_deactivate(self):
        pass

    @QtCore.Slot(int, float, float, int, float)
    def measure_saturation(self, plot_index, start, stop, steps, integration_time):
        self.integration_time = integration_time
        self._prepare_devices()

        laser_voltages = np.linspace(start, stop, steps)
        photo_diode_voltages = np.zeros(steps)
        counts = np.zeros(steps)

        for i, v in enumerate(laser_voltages):
            self._set_laser_voltage(v)
            photo_diode_voltages[i] = self._read_photo_diode()
            counts[i] = self._read_counts()

        self._set_laser_voltage(0.62)
        self._close_devices()

        self._plotlogic.set_data(plot_index=plot_index, x=photo_diode_voltages, y=counts)
        self.sigMeasuredSaturation.emit(plot_index, photo_diode_voltages, counts)
        return



    def _prepare_devices(self):
        """
        Initialize the counter and the settings of the MW device prior to starting scanning.
        """

        self._daq_card.module_state.lock()
        self._photon_samples = self._integration_time_to_samples()

        self._daq_card.prepare_counters(samples_to_acquire=self._photon_samples)
        self._daq_card.prepare_photo_diode()
        self._daq_card.module_state.unlock()

    def _close_devices(self):
        self._daq_card.close_photo_diode()
        self._daq_card.close_counters()


    def _read_counts(self):
        counts, _ = self._daq_card.read_pixel(self._photon_samples)
        if len(counts) > 0:
            counts = counts[0] / self.integration_time

        return counts

    def _read_photo_diode(self):
        voltage = self._daq_card.read_voltage()
        if len(voltage) > 0:
            voltage = voltage[0]

        return voltage

    def _set_laser_voltage(self, voltage):
        if voltage > self.max_laser_voltage:
            self.log.warn(
                f'Preventing the setting of a Laser modulation voltage greater than {self.max_laser_voltage}'
            )
            return

        self.laser_voltage = voltage
        self._laser.set_voltage(voltage)

    def _integration_time_to_samples(self):
        return round(self.integration_time * self._daq_card.get_counter_clock_frequency())