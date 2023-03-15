import numpy as np
import re

import PyDAQmx as daq

from core.module import Base
from core.configoption import ConfigOption
from interface.simple_laser_interface import SimpleLaserInterface


class NationalInstrumentsUSB6000x(Base, SimpleLaserInterface):

    # config options
    _laser_ao_channel = ConfigOption("laser_ao_channel", missing="error")
    _laser_ao_voltage_range = ConfigOption("laser_ao_voltage_range", missing="error")

    _RWTimeout = ConfigOption("read_write_timeout", default=10)

    def on_activate(self):
        self._laser_ao_task = None

        self.active_ao_tasks = []
        self._reset_hardware()
        self._start_analog_output()

    def on_deactivate(self):
        self._stop_analog_output()
        try:
            daq.DAQmxClearTask(self._laser_ao_task)
            self._scanner_ao_task = None
        except:
            self.log.exception("Could not clear AO output task.")

        self._reset_hardware()

    def _reset_hardware(self):
        retval = 0
        chan_list = [self._laser_ao_channel]
        deviceList = []

        for channel in chan_list:
            if channel is None:
                continue
            match = re.match(
                "^/(?P<dev>[0-9A-Za-z\- ]+[0-9A-Za-z\-_ ]*)/(?P<chan>[0-9A-Za-z]+)",
                channel,
            )
            if match:
                deviceList.append(match.group("dev"))
            else:
                self.log.error("Did not find device name in {0}".format(channel))

        for device in set(deviceList):
            self.log.info("Reset device {0}".format(device))
            try:
                daq.DAQmxResetDevice(device)
            except:
                self.log.exception("Could not reset NI device {0}".format(device))
                retval = -1
        return retval

    def _start_analog_output(self):
        try:
            if self._laser_ao_task is not None:
                daq.DAQmxStopTask(self._laser_ao_task)

                daq.DAQmxClearTask(self._laser_ao_task)

                self._laser_ao_task = None

            self._laser_ao_task = daq.TaskHandle()

            daq.DAQmxCreateTask("usbAO", daq.byref(self._laser_ao_task))
            daq.DAQmxCreateAOVoltageChan(
                self._laser_ao_task,
                self._laser_ao_channel,
                "Scanner AO channel 0",
                self._laser_ao_voltage_range[0],
                self._laser_ao_voltage_range[1],
                daq.DAQmx_Val_Volts,
                "",
            )
        except:
            self.log.exception("Error starting analog output task.")
            return -1
        return 0

    def _stop_analog_output(self):
        if self._laser_ao_task is None:
            return -1
        retval = 0
        try:
            daq.DAQmxStopTask(self._laser_ao_task)
        except:
            self.log.exception("Error stopping analog output.")
            retval = -1
        try:
            daq.DAQmxSetSampTimingType(self._laser_ao_task, daq.DAQmx_Val_OnDemand)
        except:
            self.log.exception("Error changing the analog output mode.")

        return retval

    def _write_ao(self, voltage, length=1, start=True):
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
                None,
            )
        except:
            self.log.exception("Error writing the analog output to the channels.")

        return self._AONWritten.value

    def set_voltage(self, voltage):
        arr = np.array([voltage], dtype='float64')
        self._write_ao(arr)

    def set_up_laser(self, ao_channel=None):
        retval = 0
        try:
            daq.DAQmxSetSampTimingType(self._scanner_ao_task, daq.DAQmx_Val_OnDemand)
        except:
            self.log.exception("Error while setting up the channel")
            retval = -1
        return retval


    # Random methods for the interface
    def get_power_range(self):
        """ Return laser power
        @return tuple(p1, p2): Laser power range in watts
        """
        pass


    def get_power(self):
        """ Return laser power
        @return float: Actual laser power in watts
        """
        pass


    def set_power(self, power):
        """ Set laer power ins watts
          @param float power: laser power setpoint in watts

          @return float: laser power setpoint in watts
        """
        pass


    def get_power_setpoint(self):
        """ Return laser power setpoint
        @return float: Laser power setpoint in watts
        """
        pass


    def get_current_unit(self):
        """ Return laser current unit
        @return str: unit
        """
        pass


    def get_current(self):
        """ Return laser current
        @return float: actual laser current as ampere or percentage of maximum current
        """
        pass


    def get_current_range(self):
        """ Return laser current range
        @return tuple(c1, c2): Laser current range in current units
        """
        pass


    def get_current_setpoint(self):
        """ Return laser current
        @return float: Laser current setpoint in amperes
        """
        pass


    def set_current(self, current):
        """ Set laser current
        @param float current: Laser current setpoint in amperes
        @return float: Laser current setpoint in amperes
        """
        pass


    def allowed_control_modes(self):
        """ Get available control mode of laser
          @return list: list with enum control modes
        """
        pass


    def get_control_mode(self):
        """ Get control mode of laser
          @return enum ControlMode: control mode
        """
        pass


    def set_control_mode(self, control_mode):
        """ Set laser control mode.
          @param enum control_mode: desired control mode
          @return enum ControlMode: actual control mode
        """
        pass


    def on(self):
        """ Turn on laser. Does not open shutter if one is present.
          @return enum LaserState: actual laser state
        """
        pass


    def off(self):
        """ Turn off laser. Does not close shutter if one is present.
          @return enum LaserState: actual laser state
        """
        pass


    def get_laser_state(self):
        """ Get laser state.
          @return enum LaserState: laser state
        """
        pass


    def set_laser_state(self, state):
        """ Set laser state.
          @param enum state: desired laser state
          @return enum LaserState: actual laser state
        """
        pass


    def get_shutter_state(self):
        """ Get shutter state. Has a state for no shutter present.
          @return enum ShutterState: actual shutter state
        """
        pass


    def set_shutter_state(self, state):
        """ Set shutter state.
          @param enum state: desired shutter state
          @return enum ShutterState: actual shutter state
        """
        pass


    def get_temperatures(self):
        """ Get all available temperatures from laser.
          @return dict: dict of name, value for temperatures
        """
        pass


    def get_temperature_setpoints(self):
        """ Get all available temperature setpoints from laser.
          @return dict: dict of name, value for temperature setpoints
        """
        pass


    def set_temperatures(self, temps):
        """ Set laser temperatures.
          @param temps: dict of name, value to be set
          @return dict: dict of name, value of temperatures that were set
        """
        pass


    def get_extra_info(self):
        """ Show dianostic information about lasers.
          @return str: diagnostic info as a string
        """
        pass