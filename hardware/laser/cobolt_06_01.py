from core.module import Base
from core.configoption import ConfigOption
from interface.simple_laser_interface import SimpleLaserInterface
from interface.simple_laser_interface import LaserState
from interface.simple_laser_interface import ShutterState

import visa

class CoboltLaser(Base, SimpleLaserInterface):

    _com_port = ConfigOption('com_port', missing='error')
    #_timeout = ConfigOption('timeout', missing = 10)

    def on_activate(self):
        try:
            self.cobolt = visa.ResourceManager().open_resource(self._com_port)
            self.cobolt.timeout = self._timeout
        except:
            self.log.exception('Could not open given COM port.')

    def on_deactivate(self):
        self.off()
        self.cobolt.close()

    def get_power(self):
        actual_power = self._send_msg('pa?')
        return actual_power

    def set_power(self, power):
        self._send_msg('p')

    def get_power_setpoint(self):
        sp = self._send_msg('p?')
        return sp

    def get_current_setpoint(self):
        curr = self._send_msg('glc?')
        return curr

    def set_current(self, current):
        self._send_msg('slc')

    def on(self):
        self._send_msg('I1')

    def off(self):
        self._send_msg('I0')

    def get_power_range(self):
        pass

    def get_current_unit(self):
        pass

    def get_current_range(self):
        pass

    def get_current(self):
        pass

    def allowed_control_modes(self):
        pass

    def get_control_mode(self):
        pass

    def set_control_mode(self, control_mode):
        pass

    def get_laser_state(self):
        pass

    def set_laser_state(self, state):
        pass

    def get_shutter_state(self):
        pass

    def set_shutter_state(self, state):
        pass

    def get_temperatures(self):
        pass

    def get_temperature_setpoints(self):
        pass

    def set_temperatures(self, temps):
        pass

    def get_extra_info(self):
        pass

    def _send_msg(self, message):
        ret_val = 0
        try:
            if message[-1] == '?':
                ret_val = self.cobolt.query(message)
            else:
                self.cobolt.write(message)
        except:
            self.log.exception('Error sending the message.')
            ret_val = -1

        return ret_val
