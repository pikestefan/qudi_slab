import time
import numpy as np
from qtpy import QtCore

from core.connector import Connector
from core.configoption import ConfigOption
from logic.generic_logic import GenericLogic
from core.util.mutex import Mutex


class SnvmLogic(GenericLogic):

    laser = Connector(interface='PixelScannerInterface')

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
        pass

    def on_deactivate(self):
        pass