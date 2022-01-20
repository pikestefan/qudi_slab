# -*- coding: utf-8 -*-

"""
This file contains the Qudi Hardware module NICard class.

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
import re

import PyDAQmx as daq

from core.module import Base
from core.configoption import ConfigOption
from interface.slow_counter_interface import CountingMode
from interface.odmr_counter_interface import ODMRCounterInterface
from interface.px_scanner_interface import PixelScannerInterface


class SnvmDummy(Base, PixelScannerInterface):

    def on_activate(self):
        pass

    def on_deactivate(self):
        pass

    def get_scanner_position(self):
        pass

    def get_position_range(self, stack=None):
        pass

    def get_scanner_count_channels(self):
        pass

    def read_pixel(self, samples=1):
        pass

    def scanner_set_position(self, x=None, y=None, stack=None):
        pass

    def close_scanner(self, stack=None):
        pass

    def close_scanner_clock(self):
        pass

    def reset_hardware(self):
        pass


