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
        self._powersupply = self.powersupply()  # access interface functions
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
        self.Vswitch_state = [0, 0, 0]  # for all three axes, 0=off, 1=on
        self.applied_current = [0, 0, 0]
        self.polarity = [1, 1, 1]

    def on_deactivate(self):
        """ Deinitialisation performed during deactivation of the module.
        """
    def Magn_to_Curr(self):
        """
        transforms spherical field components to cartesian currents
        """
        # transform field coordinates
        self.field_x = self.field_mag * np.cos(self.field_phi*2*np.pi/360) * np.sin(self.field_theta*2*np.pi/360)
        self.field_y = self.field_mag * np.sin(self.field_phi*2*np.pi/360) * np.sin(self.field_theta*2*np.pi/360)
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
            final_current = [self.current_x, self.current_y, self.current_z]
            sign_init = self.sign_func(np.multiply(self.applied_current, self.polarity))
            sign_final = self.sign_func(final_current)
            self.polarity = sign_final

        # 1 = use current containers, 0=reset currents to 0
        if task == 1:
            pol_change = sign_init != sign_final
            if 1 in pol_change:
                flip_current = np.where(pol_change == 0, 1, 0)
                self._powersupply.set_current_all(np.multiply(final_current,flip_current)) #before triggering the switch, set those channels to 0
                self._arduino.write_and_read('{:d},{:d},{:d}'.format(*[task, task, task] != sign_final))    # #False = 0 = pos. pol., True = 1 = neg.pol.
                time.sleep(1.2)
            self._powersupply.set_current_all(final_current) # set to final current
        elif task == 0:
            print('reset to 0!')
            self._powersupply.set_current_all([0, 0, 0])
            self._arduino.write_and_read('{:d},{:d},{:d}'.format(*[task, task, task]))
        else:
            self.log.error('task for power supplied not given. Must be 0 or 1!')
            raise

    def get_real_current(self, channel):
        # updates the current containers directly from the power supply
        current = self._powersupply.get_real_current(channel)
        return current


    def shut_down_channels(self, state):
        # state = 2 :disable all outputs; state = 0 : enable all outputs
        self._powersupply.set_device_state(state)

    def set_Vlimit(self, Vlimit):
        #call the hardware and set limit to value
        self._powersupply.set_Vlimit(Vlimit)