# -*- coding: utf-8 -*-

"""
Interface file to control processes in PID control.

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

from core.interface import abstract_interface_method
from core.meta import InterfaceMetaclass


class PowerSupplyInterface(metaclass=InterfaceMetaclass):
    @abstract_interface_method
    def set_voltage(self, setvoltage, channel):
        """
        set the voltage of specified channel to a specific value
        """
    @abstract_interface_method
    def set_current(self, setcurrent, channel):
        """
        set the current of specified channel to a specific value
        """

    @abstract_interface_method
    def get_real_current(self,channel):
        """
        queries the current of target channel
        """

    @abstract_interface_method
    def get_theo_current_all(self, field_x, field_y, field_z):
        """
        estimates currents based on hardware defined coefficients
        """
    @abstract_interface_method
    def set_current_all(self, setcurrent_1, setcurrent_2, setcurrent_3):
        """
        set the current for all channels at the same time
        """

    @abstract_interface_method
    def set_device_state(self, state):
        """
        turn device on or off
        """

    @abstract_interface_method
    def get_coefficients(self):
        """
        returns the hardware specific coeeficients
        """
