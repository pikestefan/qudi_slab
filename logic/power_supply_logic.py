from logic.generic_logic import GenericLogic
import numpy as np
from core.connector import Connector

class PowerSupplyLogic(GenericLogic):
    """This is the Logic class for ODMR."""

    # declare connectors
    powersupply = Connector(interface='PowerSupplyInterface')

    def on_activate(self):
        """
        Initialisation performed during activation of the module.
        """
        self._powersupply = self.powersupply() #access interface functions

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

    def apply_magnetic_field(self):
        """
        call hardware and pass the required fields
        """
        # call hardware
        self._powersupply.set_current_all(self.current_x, self.current_y, self.current_z)

    def get_real_currents(self, channel):
        # updates the current containers directly from the power supply
        current= self._powersupply.get_real_current(channel)
        return current


    def shut_down_channels(self, state):
        # state = 2 :disable all outputs; state = 0 : enable all outputs
        self._powersupply.set_device_state(state)