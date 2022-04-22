from qtpy.QtWidgets import QSpinBox, QApplication
"""
This widget is a spin box that displays strings in place of the integer value

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
"""

class StringSpinBox(QSpinBox):
    def __init__(self, *args, **kwargs):
        super(StringSpinBox, self).__init__(*args, **kwargs)
        self._strings = None


    def setStrings(self, strings):
        if not isinstance(strings, list):
            raise(Exception("The input must be a list or a np.ndarray"))

        self._strings = tuple(strings)
        self._values = dict(zip(strings, range(len(strings))))
        self.setRange(0, len(self._values)-1)

    def textFromValue(self, value):
        if self._strings is not None:
            return self._strings[value]
        else:
            return super().textFromValue(value)

if __name__ == '__main__':
    import sys
    app = QApplication(sys.argv)
    gne = StringSpinBox()
    gne.setStrings(['abc','bcd'])
    gne.show()
    app.exec_()