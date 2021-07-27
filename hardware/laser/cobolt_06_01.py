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
