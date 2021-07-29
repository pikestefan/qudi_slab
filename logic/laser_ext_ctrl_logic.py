import time
import numpy as np
from qtpy import QtCore

from core.connector import Connector
from core.configoption import ConfigOption
from logic.generic_logic import GenericLogic
from core.util.mutex import Mutex
from interface.simple_laser_interface import ControlMode, ShutterState, LaserState

class LaserExtCtrlLogic(GenericLogic):

    laser = Connector(interface='ExtCtrlLaserInterface')
    option_placeholder = ConfigOption('option_placeholder', 0)

    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)
        """
        self.log.debug('The following configuration was found.')

        # checking for the right configuration
        for key in config.keys():
            self.log.info('{0}: {1}'.format(key, config[key]))
        """
        self.threadlock = Mutex()

    def on_activate(self):
        self._laser = self.laser()


    def on_deactivate(self):
        pass

    def set_power(self, power):
        print("Got called")
        self._laser.set_power_extctrl(power)
