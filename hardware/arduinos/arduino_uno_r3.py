import time
import visa
from serial import Serial
from core.module import Base
from core.configoption import ConfigOption
from interface.process_control_interface import ProcessControlInterface
from interface.arduino_interface import ArduinoInterface

class ArduinoUnoR3(Base, ArduinoInterface):
    """
    Example config:
        voltage_generator:
            module.Class: 'power_supply.power_supply_dummy.PowerSupplyDummy'
            voltage_min: 0
            voltage_max_1: 30
            voltage_max_2: 30
            voltage_max_3: 5
            current_max: 3

    """
    _usb_address = ConfigOption('usb_address', missing='error')
    _baud_rate = ConfigOption('baud_rate', missing='error')

    def on_activate(self):
        """
        Initialisation performed during activation of the module.
        """
        self.rm = visa.ResourceManager()
        try:
            self._serialcomm = Serial(self._usb_address, 115200)
            self._serialcomm.timeout = 0.05
        except:
            self.log.error('Could not connect to port "{}". Check '
                           'whether address exists and reload '
                           'module!'.format(self._usb_address))
            raise


    def on_deactivate(self):
        """ Deinitialisation performed during deactivation of the module.
        """
        self._serialcomm.close()

    def write_and_read(self, text):
        """
        pass string to arduino and receive return
        """
        self._serialcomm.write(text.encode())
        time.sleep(0.1)
        output = self._serialcomm.readline().decode('ascii')
        return output