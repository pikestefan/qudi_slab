from logic.generic_logic import GenericLogic
import numpy as np
from core.connector import Connector
import time

class PowerSupplyLogic(GenericLogic):
    """This is the Logic class for ODMR."""

    # declare connectors
    powersupply = Connector(interface='PowerSupplyInterface')
    arduino = Connector(interface='ArduinoInterface')
    def on_activate(self):
        """
        Initialisation performed during activation of the module.
        """
        self._powersupply = self.powersupply() #access interface functions
        self._arduino = self.arduino()  # access interface functions
        # variable containers for the currents and the field components
        self.field_x = 0
        self.field_y = 0
        self.field_z = 0
        self.field_mag = 0
        self.field_theta = 0
        self.field_phi = 0
        self.current_x = 0
        self.current_y = 0
        self.current_z = 0
        self.Vswitch_state = [0, 0, 0] #for all three axes, 0=off, 1=on
        self.applied_current = [0, 0, 0]

    def on_deactivate(self):
        """ Deinitialisation performed during deactivation of the module.
        """
    def Magn_to_Curr(self):
        """
        transforms spherical field components to cartesian currents
        """
        # transform field coordinates
        self.field_x = self.field_mag * np.cos(self.field_phi*2*np.pi/360) * np.sin(self.field_theta*2*np.pi/360)
        self.field_y = self.field_mag * np.sin(self.field_phi*2*np.pi/360) * np.cos(self.field_theta*2*np.pi/360)
        self.field_z = self.field_mag * np.cos(self.field_theta*2*np.pi/360)
        current = self._powersupply.get_theo_current_all(self.field_x, self.field_y, self.field_z)
        self.current_x = current[0]
        self.current_y = current[1]
        self.current_z = current[2]

    def sign_func(self, array):
        j = 0
        out = np.zeros(np.size(array))
        for i in array:
            if np.sign(i) >= 0:
                out[j] = 1
            else:
                out[j] = -1
            j += 1
        return out

    def apply_magnetic_field(self, task):
        """
        call hardware and pass the required fields
        Check, if hardware allows negative outputs, otherwise call V-switch (default: off)
        """
        if self._powersupply._negative_polarity == False:
            sign_init = self.sign_func(self.applied_current)
            sign_final = self.sign_func([self.current_x, self.current_y, self.current_z])
            self._arduino.write_and_read('{:d},{:d},{:d}'.format(*sign_init != sign_final))
            """
            for i in [0, 1, 2]:
                if sign_final[i] != sign_init[i]:
                    self._arduino.write_and_read(str(i)+",on")
                    #print('switched')
                elif sign_final[i] == sign_init[i]:
                    self._arduino.write_and_read(str(i)+",off")
                    #print('not switched')
                else:
                    self.log.error('Could not evaluate the expected V-switch state!')
                    raise
                #time.sleep(1.5)
            """
        #1 = use current containers, 0=reset currents to 0
        if task == 1:
            self._powersupply.set_current_all(self.current_x, self.current_y, self.current_z)
        elif task == 0:
            self._powersupply.set_current_all(0, 0, 0)
        else:
            self.log.error('task for power supplied not given. Must be 0 or 1!')
            raise

def get_real_currents(self, channel):
    # updates the current containers directly from the power supply
    current= self._powersupply.get_real_current(channel)
    return current


def shut_down_channels(self, state):
    # state = 2 :disable all outputs; state = 0 : enable all outputs
    self._powersupply.set_device_state(state)

def set_Vlimit(self, Vlimit):
    #call the hardware and set limit to value
    self._powersupply.set_Vlimit(Vlimit)