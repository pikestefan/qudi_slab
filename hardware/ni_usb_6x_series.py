import numpy as np
import re

import PyDAQmx as daq

from core.module import Base
from core.configoption import ConfigOption

class NationalInstrumentsUSB6000x(Base):

    # config options
    _laser_ao_channel = ConfigOption('laser_ao_channel', missing='error')
    _laser_ao_voltage_range = ConfigOption('laser_ao_voltage_range', missing='error')

    _RWTimeout = ConfigOption('read_write_timeout', default=10)
    def on_activate(self):
        self._laser_ao_task = None

        self.active_ao_tasks = []

    def on_deactivate(self):
        self._stop_analog_output()
        try:
            daq.DAQmxClearTask(self._laser_ao_task)
            self._scanner_ao_task = None
        except:
            self.log.exception('Could not clear AO output task.')

        self._reset_hardware()

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
                self._laser_ao_channel,
                'Scanner AO channel 0',
                self._laser_ao_voltage_range[0],
                self._laser_ao_voltage_range[1],
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

    def _write_ao(self, voltage, length=1, start=False):
        self._AONWritten = daq.int32()
        try:
            daq.DAQmxWriteAnalogF64(
                self._laser_ao_task,
                length,
                start,
                self._RWTimeout,
                daq.DAQmx_Val_GroupByChannel,
                voltage,
                daq.byref(self._AONWritten),
                None
            )
        except:
            self.log.exception('Error writing the analog output to the channels.')

        return self._AONWritten.value

    def set_up_laser(self, ao_channel = None):
        retval = 0
        try:
            daq.DAQmxSetSampTimingType(self._scanner_ao_task, daq.DAQmx_Val_OnDemand)
        except:
            self.log.exception('Error while setting up the channel')
            retval = -1
        return retval
