# -*- coding: utf-8 -*-

"""
This file contains the Qudi GUI for general Confocal control.

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
import os
import pyqtgraph as pg
import time
from math import ceil

from core.connector import Connector
from core.configoption import ConfigOption
from core.statusvariable import StatusVar
from qtwidgets.scan_plotwidget import ScanImageItem
from gui.guibase import GUIBase
from gui.guiutils import ColorBar
from gui.colordefs import ColorScaleInferno, BlackAndWhite
from gui.colordefs import QudiPalettePale as palette
from gui.fitsettings import FitParametersWidget
from qtpy import QtCore
from qtpy import QtGui
from qtpy import QtWidgets
from qtpy import uic
from qtpy.QtWidgets import QColorDialog, QFontDialog


class VectorMagnetMainWindow(QtWidgets.QMainWindow):
    """ Create the Main window based on the corresponding *.ui file. """

    sigPressKeyBoard = QtCore.Signal(QtCore.QEvent)
    sigDoubleClick = QtCore.Signal()

    def __init__(self):
        # Get the path to the *.ui file
        this_dir = os.path.dirname(__file__)
        ui_file = os.path.join(this_dir, 'ui_vector_gui.ui')
        self._doubleclicked = False

        # Load it
        super(VectorMagnetMainWindow, self).__init__()
        uic.loadUi(ui_file, self)
        self.show()

class VectorMagnetGui(GUIBase):
    #declare connectors
    power_supply_logic = Connector(interface='PowerSupplyLogic')

    def __init__(self,config, **kwargs):
        super().__init__(config=config, **kwargs)


    def on_activate(self):
        self.initMainUI()
        self._vector_logic = self.power_supply_logic()

        # Connecting user interactions
        self._mainwindow.ApplyButton.clicked.connect(self.ApplyButton)
        self._mainwindow.ResetButton.clicked.connect(self.ResetButton)
        self._mainwindow.Magnitude.valueChanged.connect(self.update_current)
        self._mainwindow.Theta.valueChanged.connect(self.update_current)
        self._mainwindow.Phi.valueChanged.connect(self.update_current)
        self._mainwindow.ChannelsOnOffBox.stateChanged.connect(self.check_channels_off)
        self._mainwindow.LimitSlider.valueChanged.connect(self.LimitSlider_changed)

    def on_deactivate(self):
        """ Deinitialisation performed during deactivation of the module.
        """

    def initMainUI(self):
        self._mainwindow = VectorMagnetMainWindow()
        self.show()

    def show(self):
        # Make main window visible and put it above all other windows.
        # Show the Main Confocal GUI:
        self._mainwindow.show()
        self._mainwindow.activateWindow()
        self._mainwindow.raise_()

    def update_display(self):
        """
        triggers the logic to update the currents
        """
        # update display
        self._vector_logic.applied_current = np.asarray([self._vector_logic.get_real_currents(1), self._vector_logic.get_real_currents(2), self._vector_logic.get_real_currents(3)], dtype = int)
        self._mainwindow.Disp1.display(str(round(self._vector_logic.applied_current[0], 4)))
        self._mainwindow.Disp2.display(str(round(self._vector_logic.applied_current[1], 4)))
        self._mainwindow.Disp3.display(str(round(self._vector_logic.applied_current[2], 4)))

    def ApplyButton(self):
        """
        apply settings to power supply
        """
        print('I was here')
        if self._mainwindow.ChannelsOnOffBox.checkState() == 0:
            # call logic
            self._vector_logic.apply_magnetic_field(1)
            # update GUI
            self.update_display()

    def ResetButton(self):
        """
        reset ALL current outputs to 0A
        """
        self._vector_logic.apply_magnetic_field(0)
        self.update_display()


    def update_current(self):
        """
        calls logic to transform get currents and updates currents spinboxes
        """
        #update field containers
        self._vector_logic.field_mag = self._mainwindow.Magnitude.value()
        self._vector_logic.field_theta = self._mainwindow.Theta.value()
        self._vector_logic.field_phi = self._mainwindow.Phi.value()
        # call logic
        self._vector_logic.Magn_to_Curr()
        # update current
        self._mainwindow.ChannelV1.setPlainText(str(round(self._vector_logic.current_x, 4)))
        self._mainwindow.ChannelV2.setPlainText(str(round(self._vector_logic.current_y, 4)))
        self._mainwindow.ChannelV3.setPlainText(str(round(self._vector_logic.current_z, 4)))

    def check_channels_off(self):
        #shut down channel outputs
        if self._mainwindow.ChannelsOnOffBox.checkState() == 2: #0::unchecked, 1::partially checked, 2::checked
            self._vector_logic.shut_down_channels("Off")
            self._mainwindow.Disp1.setStyleSheet("""QLCDNumber { 
                                                    background-color: darkred; 
                                                    color: white; }""")
            self._mainwindow.Disp2.setStyleSheet("""QLCDNumber { 
                                                    background-color: darkred; 
                                                    color: white; }""")
            self._mainwindow.Disp3.setStyleSheet("""QLCDNumber { 
                                                    background-color: darkred; 
                                                    color: white; }""")
        elif self._mainwindow.ChannelsOnOffBox.checkState() == 0: #0::unchecked, 1::partially checked, 2::checked
            self._vector_logic.shut_down_channels("On")
            self._mainwindow.Disp1.setStyleSheet("""""")
            self._mainwindow.Disp2.setStyleSheet("""""")
            self._mainwindow.Disp3.setStyleSheet("""""")
    def LimitSlider_changed(self):
        #set V-limits according to slider position
        Vlimit = self._mainwindow.LimitSlider.value()*10
        self._vector_logic.set_Vlimit(Vlimit)
        self._mainwindow.LimitValue.setText('{}'.format(Vlimit))


    """ 
        global value
        value = ""
        def __init__(self): #config, **kwargs):
            super().__init__(): #config=config, **kwargs)

        def on_activate(self):
            self.initMainUI()
            # Connecting user interactions
            self._mainwindow.ReadButton.clicked.connect(self.read_clicked)
            self._mainwindow.DeleteButton.clicked.connect(self.delete_clicked)
            self._mainwindow.colourBox.currentIndexChanged.connect(self.change_colour)

        def initMainUI(self):
            self._mainwindow = VectorMagnetMainWindow()
            self.show()

        def show(self):
            #Make main window visible and put it above all other windows.
            # Show the Main Confocal GUI:
            self._mainwindow.show()
            self._mainwindow.activateWindow()
            self._mainwindow.raise_()

        def read_clicked(self):
            global value
            value = self._mainwindow.sendBox.displayText()
            self._mainwindow.readBox.setPlainText(value)
        def delete_clicked(self):
            self._mainwindow.readBox.setPlainText("")

        def change_colour(self):
            global value
            color = self._mainwindow.colourBox.currentText()
            if color == "blue":
                self._mainwindow.readBox.setTextColor(QtGui.QColor(0,255,0,255))
            elif color == "red":
                self._mainwindow.readBox.setTextColor(QtGui.QColor(255, 0, 0, 255))
            elif color == "green":
                self._mainwindow.readBox.setTextColor(QtGui.QColor(0, 0, 255, 255))
            elif color == "colour":
                self._mainwindow.readBox.setTextColor(QtGui.QColor(0, 0, 0, 255))
            #self.read_clicked()
            self._mainwindow.readBox.setPlainText(value)
        


if __name__ == '__main__':
            from qtpy.QtWidgets import QApplication
            import sys

            app = QApplication(sys.argv)
            gne = VectorMagnetGui()
            gne.on_activate()
            sys.exit(app.exec_())
            """













