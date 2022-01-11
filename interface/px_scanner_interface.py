# -*- coding: utf-8 -*-

"""
This module contains the Qudi interface file for confocal scanner.

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


class PixelScannerInterface(metaclass=InterfaceMetaclass):
    """ This is the Interface class to define the controls for a scanner device

    A scanner device is a hardware that can move on up to 4 multiple axis with repeatability in position and measure
    something at each point of a line in the parameter space trajectory.
    The key idea of this interface is that the hardware handle the timing of the acquisition to remove as much delay as
    possible.

    A typical example of such hardware in the lab is the use with the NI card to move scanner/mirrors and record the
    luminescence at each position to build, line by line, a raster scan.

    ---
    This code, while trying to be compatible with any number of axis, is as of now only used with 3 or 4 axis.
    Using less will raise errors in confocal_logic.

    """

    @abstract_interface_method
    def reset_hardware(self):
        """ Resets the hardware, so the connection is lost and other programs can access it.

        @return int: error code (0:OK, -1:error)
        """
        pass

    @abstract_interface_method
    def get_position_range(self):
        """ Returns the physical range of the scanner.

        @return float [N][2]: array of N ranges with an array containing lower and upper limit, preferably in SI unit.

        """
        pass

    @abstract_interface_method
    def set_position_range(self, myrange=None):
        """ Sets the physical range of the scanner.

        Deprecated : This range should not be accessible by logic. TODO: Discuss and remove from interface ?

        @param float [N][2] myrange: array of N ranges with an array containing lower and upper limit

        @return int: error code (0:OK, -1:error)
        """
        pass

    @abstract_interface_method
    def set_voltage_range(self, myrange=None):
        """ Sets the voltage range of the NI Card.

        Deprecated : This range should not be accessible by logic. TODO: Discuss and remove from interface ?

        @param float [2] myrange: array containing lower and upper limit

        @return int: error code (0:OK, -1:error)
        """
        pass

    @abstract_interface_method
    def get_scanner_count_channels(self):
        """ Returns the list of channels that are recorded while scanning an image.

        @return list(str): channel names

        Most methods calling this might just care about the number of channels.
        """
        pass

    @abstract_interface_method
    def scanner_set_position(self, x=None, y=None, z=None, a=None):
        """ Move stage to x, y, z, a (where a is the fourth channel).

        @param float x: position in x-direction (in axis unit)
        @param float y: position in y-direction (in axis unit)
        @param float z: position in z-direction (in axis unit)
        @param float a: position in a-direction (in axis unit)

        @return int: error code (0:OK, -1:error)

        If a value is not set or set to None, the actual value is implied.
        """
        pass

    @abstract_interface_method
    def get_scanner_position(self):
        """ Get the current position of the scanner hardware.

        @return tuple(float): current position as a tuple. Ex : (x, y, z, a).
        """
        pass

    @abstract_interface_method
    def close_scanner(self):
        """ Closes the scanner and cleans up afterwards.

        @return int: error code (0:OK, -1:error)

        TODO: Give a detail explanation how it is used in practice and why it is necessary.
        """
        pass

    @abstract_interface_method
    def close_scanner_clock(self, power=0):
        """ Closes the clock and cleans up afterwards.

        @return int: error code (0:OK, -1:error)

        TODO: Give a detail explanation how it is used in practice and why it is necessary.
        """
        pass

    @abstract_interface_method
    def read_pixel(self):

        pass





