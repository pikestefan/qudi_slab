import time
import visa
import numpy as np
from core.module import Base
from core.configoption import ConfigOption
from interface.process_control_interface import ProcessControlInterface
from interface.power_supply_interface import PowerSupplyInterface

class KeythleyPowerSupply(Base, PowerSupplyInterface):
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
    _usb_timeout = ConfigOption('usb_timeout', 10, missing='warn')

    ###### CLEAN THIS UP?
    _voltage_min = ConfigOption('voltage_min', missing='error')
    _voltage_max_1 = ConfigOption('voltage_max_1', missing='error')
    _voltage_max_2 = ConfigOption('voltage_max_2', missing='error')
    _voltage_max_3 = ConfigOption('voltage_max_3', missing='error')
    _current_max = ConfigOption('current_max', missing='error')

    _negative_polarity = ConfigOption('negative_polarity', missing='error')

    def on_activate(self):
        """
        Initialisation performed during activation of the module.
        """
        self.rm = visa.ResourceManager()
        self._usb_connection = self.rm.open_resource(
                self._usb_address,
                timeout=self._usb_timeout * 1000)
        #print(self._ask('*IDN?'))
        self._usb_connection.write('SYST:REM')
        self._write('OUTP ON')
        # To set the voltage limits. This is saved in the settings
        self.set_Vlimit(5)


        # except:
        #     self.log.error('Could not connect to the usb address "{}". Check '
        #                    'whether address exists and reload '
        #                    'module!'.format(self._usb_address))

    def on_deactivate(self):
        """ Deinitialisation performed during deactivation of the module.
        """
        self._write('APP:CURR {},{},{}'.format(0, 0, 0))
        self._write('OUTP OFF')
        self._usb_connection.close()

    def set_current_all(self, setcurrent):
        # when currents are negative, call arduino to switch
        print('[x,y,z] set to [', round(setcurrent[0], 4), ', ', round(setcurrent[1], 4), ', ', round(setcurrent[2], 4), '] A')
        self._write('APP:CURR {},{},{}'.format(np.absolute(setcurrent[0]),np.absolute(setcurrent[1]),np.absolute(setcurrent[2])))



    def get_real_current(self, channel):
        """
        get the current that is actually applied by the device, internally measured
        """
        self._write('INST:NSEL {}'.format(channel))
        current = self._ask('MEAS:CURR?')
        return current

    def get_theo_current_all(self, field_x, field_y, field_z):
        """
        K_xyz include device dependent parameters, that need to be characterized, as well as the physical parameters
        """
        K_x = 2.939
        K_y = 2.802
        K_z = 3.019
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
        if state == "Off":
            self._write('Outp OFF')
            print('Device state: ', state)
        elif state == "On":
            self._write('Outp ON')
            print('Device state: ', state)
        else:
            print('Cannot change state of device!')
            raise

    def _ask(self, question):
        """ Ask wrapper.

        @param str question: a question to the device

        @return: the received answer
        """
        return self._usb_connection.query(question)

    def _write(self, command, wait=True):
        """ Write wrapper.

        @param str command: a command to the device
        @param bool wait: optional, is the wait statement should be skipped.

        @return: str: the statuscode of the write command.
        """
        statuscode = self._usb_connection.write(command)
        if wait:
            self._usb_connection.write('*WAI')
        return statuscode

    def set_Vlimit(self, Vlimit):
        #sets limits to all channels
        self._write('INST:NSEL 1')
        self._write('Volt:limit {}'.format(Vlimit))
        self._write('Volt:LIMit:STATe 1')
        self._write('Appl ch1,max,min')
        self._write('INST:NSEL 2')
        self._write('Volt:limit {}'.format(Vlimit))
        self._write('Volt:LIMit:STATe 1')
        self._write('Appl ch2,max,min')
        self._write('INST:NSEL 3')
        self._write('Volt:limit {}'.format(Vlimit))
        self._write('Volt:LIMit:STATe 1')
        self._write('Appl ch3,max,min')

    ##### these two functions need to be here!
    def set_current(self, setcurrent, channel):
        """

        """

    def set_voltage(self, setvoltage, channel):
        """

        """
# if __name__ == '__main__':
    # print("Called Power Supply")
    # rm = visa.ResourceManager()
    # print(rm.list_resources())
    # usb = rm.open_resource(
    #     'ASRL6::INSTR',
    #     timeout=10*1000
    # )
    # print(usb)
    # print(usb.write('SYSTem:REMote'))
    # usb.write('*WAI')
    # print(usb.write('Outp ON'))
    # usb.write('*WAI')
