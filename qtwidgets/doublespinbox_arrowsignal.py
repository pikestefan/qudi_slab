# -*- coding: utf-8 -*-

"""
QCheckBox with a callback function for acceptance or denial of state change.

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

from qtpy import QtWidgets, QtCore, QtGui


class DoubleSpinBoxArrowSignal(QtWidgets.QDoubleSpinBox):

    stepClicked = QtCore.Signal(float)

    def stepBy(self, step):
        value = self.value()
        super(QtWidgets.QDoubleSpinBox, self).stepBy(step)
        if self.value() != value:
            self.stepClicked.emit(self.value())
