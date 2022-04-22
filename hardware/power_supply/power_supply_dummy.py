import time
import visa

from core.module import Base
from core.configoption import ConfigOption
from interface.process_control_interface import ProcessControlInterface
from interface.power_supply_interface import PowerSupplyInterface

class PowerSupplyDummy(Base, PowerSupplyInterface):
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


    _voltage_min = ConfigOption('voltage_min', missing='error')
    _voltage_max_1 = ConfigOption('voltage_max_1', missing='error')
    _voltage_max_2 = ConfigOption('voltage_max_2', missing='error')
    _voltage_max_3 = ConfigOption('voltage_max_3', missing='error')
    _current_max = ConfigOption('current_max', missing='error')

    def on_activate(self):
        """ Initialisation performed during activation of the module.
        """


    def on_deactivate(self):
        """ Deinitialisation performed during deactivation of the module.
        """
        pass

    def set_current_all(self, setcurrent_x, setcurrent_y, setcurrent_z):
        #[TO BE INCLUDED]
        print('x set to ', setcurrent_x)
        print('y set to ', setcurrent_y)
        print('z set to ', setcurrent_z)


    def get_real_current(self, channel):
        #MEAS:CURR?
        current = 123*channel#FETCh:CURRent?
        return current

    def get_theo_current_all(self, field_x, field_y, field_z):
        K_x = 3
        K_y = 2
        K_z = 1
        setcurrent_x = field_x / K_x
        setcurrent_y = field_y / K_y
        setcurrent_z = field_z / K_z
        return setcurrent_x, setcurrent_y, setcurrent_z

    def get_coefficients(self):
        """
        returns the hardware specific coeeficients
        """
    def set_device_state(self, state):
        """
        turn all channels on or off
        """
        print('Device state: ', state)
    def set_current(self, setcurrent, channel):
        """
        set the current of specified channel to a specific value
        """
    def set_voltage(self, setvoltage, channel):
        """
        set the voltage of specified channel to a specific value
        """
    ##### [TO BE INCLUDED]
    def _ask(self, question):
        """ Ask wrapper.

        @param str question: a question to the device

        @return: the received answer
        """

    def _write(self, command, wait=True):
        """ Write wrapper.

        @param str command: a command to the device
        @param bool wait: optional, is the wait statement should be skipped.

        @return: str: the statuscode of the write command.
        """
