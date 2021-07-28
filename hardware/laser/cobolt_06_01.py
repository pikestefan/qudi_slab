from core.module import Base
from core.configoption import ConfigOption
from interface.external_control_laser_interface import ExtCtrlLaserInterface
from interface.external_control_laser_interface import ControlMode
from interface.simple_laser_interface import ShutterState

import re
import visa
import PyDAQmx as daq

class CoboltLaser(Base, ExtCtrlLaserInterface):

    _com_port = ConfigOption('com_port', missing='error')
    _enable_am = ConfigOption('analog_mod_enabled', True, missing='info')
    _enable_dm = ConfigOption('digital_mod_enabled', True, missing='info')
    _laser_com_timeout = ConfigOption('laser_com_timeout', default=10)

    #Requires full address of the output, e.g. /dev(numberhere)/ao0
    _NI_analog_channel = ConfigOption('ao_output_channel', missing='error')
    _NI_voltage_range = ConfigOption('ao_voltage_range', [0, 10], missing='info')
    _RW_timeout = ConfigOption('daq_read_write_timeout', default=10)

    def on_activate(self):
        self._laser_ao_task = None
        try:
            self.cobolt = visa.ResourceManager().open_resource(self._com_port)
            self.cobolt.timeout = self._laser_com_timeout
        except:
            self.log.exception('Could not open given COM port.')

        self.laser_initialize()
        self.set_up_laser_ao_channel()

    def on_deactivate(self):
        self.set_control_mode(ControlMode.POWER)
        self.cobolt.close()

        self._stop_analog_output()
        try:
            daq.DAQmxClearTask(self._laser_ao_task)
            self._scanner_ao_task = None
        except:
            self.log.exception('Could not clear AO output task.')

        self._reset_hardware()

    def get_power(self):
        actual_power = self._send_msg('pa?')
        return float(actual_power)

    def laser_initialize(self, am_on=True, dm_on=True):
        try:
            self.set_control_mode(ControlMode.MODULATION)
            self.set_analog_low_impedance(False)
            self.enable_analog_mod(self._enable_am)
            self.enable_digital_mod(self._enable_dm)
        except:
            self.log.exception('Error occured during laser setup.')

    def set_power(self, power):
        """
        Set power in mW
        """
        self._send_msg(f'p {power*1e-3}')

    def set_power_extctrl(self, voltage):
        voltage = np.array([voltage])
        self._write_ao(self, voltage, start=True)

    def get_power_setpoint(self):
        sp = self._send_msg('p?')
        return float(sp)

    def get_current_setpoint(self):
        curr = self._send_msg('glc?')
        return float(curr)

    def set_current(self, current):
        """
        Set current in mA
        """
        self._send_msg(f'slc {current*1e-3}')

    def on(self):
        self._send_msg('l1')

    def off(self):
        self._send_msg('l0')

    def get_current(self):
        actual_curr = self._send_msg('rlc')
        return float(actual_curr)

    def get_current_unit(self):
        return 'mA'

    def set_analog_low_impedance(self, low_active=False):
        val = int(low_active)
        self._send_msg(f'salis {val}')

    def get_analog_low_impedance(self):
        val = self._send_msg('galis?')
        return bool(val)

    def enable_analog_mod(self, set_enabled = True):
        val = int(set_enabled)
        self._send_msg(f'sames {val}')

    def enable_digital_mod(self, set_enabled = True):
        val = int(set_enabled)
        self._send_msg(f'sdmes {val}')

    def set_control_mode(self, ctrl_mode):
        if ctrl_mode == ControlMode.POWER:
            msg = 'pm'
        elif ctrl_mode == ControlMode.CURRENT:
            msg = 'cm'
        elif ctrl_mode == ControlMode.MODULATION:
            msg = 'em'
        else:
            self.log.exception('Invalid laser mode.')

        self._send_msg(msg)

    def _send_msg(self, message):
        try:
            self.cobolt.write(message)
            ret_val = self.cobolt.read().strip()
        except:
            self.log.exception('Error sending the message.')
            ret_val = -1
        return ret_val

    # Here there are the daq controls
    def _start_analog_output(self):
        try:
            if self._laser_ao_task is not None:
                daq.DAQmxStopTask(self._laser_ao_task)

                daq.DAQmxClearTask(self._laser_ao_task)

                self._laser_ao_task = None

            self._laser_ao_task = daq.TaskHandle()

            daq.DAQmxCreateTask('usbAO', daq.byref(self._laser_ao_task))
            daq.DAQmxCreateAOVoltageChan(
                self._laser_ao_task,
                self._NI_analog_channel,
                'Laser ao channel',
                self._NI_voltage_range[0],
                self._NI_voltage_range[1],
                daq.DAQmx_Val_Volts,
                '')
        except:
            self.log.exception('Error starting analog output task.')
            return -1
        return 0

    def _stop_analog_output(self):
        if self._laser_ao_task is None:
            return -1
        retval = 0
        try:
            daq.DAQmxStopTask(self._laser_ao_task)
        except:
            self.log.exception('Error stopping analog output.')
            retval = -1
        try:
            daq.DAQmxSetSampTimingType(self._laser_ao_task, daq.DAQmx_Val_OnDemand)
        except:
            self.log.exception('Error changing the analog output mode.')

        return retval

    def _write_ao(self, voltages, length=1, start=False):
        self._AONWritten = daq.int32()
        try:
            daq.DAQmxWriteAnalogF64(
                self._laser_ao_task,
                length,
                start,
                self._RWTimeout,
                daq.DAQmx_Val_GroupByChannel,
                voltages,
                daq.byref(self._AONWritten),
                None
            )
        except:
            self.log.exception('Error writing the analog output to the channels.')

        return self._AONWritten.value

    def set_up_laser_ao_channel(self):
        retval = 0
        try:
            daq.DAQmxSetSampTimingType(self._scanner_ao_task, daq.DAQmx_Val_OnDemand)
        except:
            self.log.exception('Error while setting up the channel')
            retval = -1
        return retval

    def _reset_hardware(self):
        retval = 0
        chan_list = [self._laser_ao_channel]
        deviceList = []

        for channel in chan_list:
            if channel is None:
                continue
            match = re.match('^/(?P<dev>[0-9A-Za-z\- ]+[0-9A-Za-z\-_ ]*)/(?P<chan>[0-9A-Za-z]+)',
                             channel)
            if match:
                deviceList.append(match.group('dev'))
            else:
                self.log.error('Did not find device name in {0}'.format(channel)
                               )

        for device in set(deviceList):
            self.log.info('Reset device {0}'.format(device))
            try:
                daq.DAQmxResetDevice(device)
            except:
                self.log.exception('Could not reset NI device {0}'.format(device))
                retval = -1
        return retval

    # Here lies the graveyard of options not enabled by Cobolt.
    def get_power_range(self):
        pass

    def get_current_range(self):
        pass

    def allowed_control_modes(self):
        pass

    def get_control_mode(self):
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
