# -*- coding: utf-8 -*-

"""
This file contains the Qudi Hardware module NICard class.

Qudi is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Qudi is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Qudi. If not, see <http://www.gnu.org/licenses/>.

Copyright (c) the Qudi Developers. See the COPYRIGHT.txt file at the
top-level directory of this distribution and at <https://github.com/Ulm-IQO/qudi/>
"""

import numpy as np
import re

import PyDAQmx as daq

from core.module import Base
from core.configoption import ConfigOption
from interface.slow_counter_interface import CountingMode
from interface.odmr_counter_interface import ODMRCounterInterface
from interface.snvm_scanner_interface import SnvmScannerInterface


class NationalInstrumentsXSeriesPxScan(Base, SnvmScannerInterface, ODMRCounterInterface):
    """ A National Instruments device that can count and control microvave generators.

    !!!!!! NI USB 63XX, NI PCIe 63XX and NI PXIe 63XX DEVICES ONLY !!!!!!

    See [National Instruments X Series Documentation](@ref nidaq-x-series) for details.

    stable: Kay Jahnke, Alexander Stark

    Example config for copy-paste:

    nicard_6343:
        module.Class: 'national_instruments_px_scan.NationalInstrumentsXSeriesPxScan'
        photon_sources:
            - /Dev1/PFI8
        #    - /Dev1/PFI9
        counter_clock: '100kHzTimebase'
        counter_channels:
            - /Dev1/Ctr1
        counter_ai_channels:
            - /Dev1/AI0
        counter_voltage_range: [-10, 10]
        sample_scanner_ao_channels:
            - /Dev1/AO0
            - /Dev1/AO1
        tip_scanner_ao_channels:
            - /Dev1/AO2
            - /Dev1/AO3
        sample_scanner_ai_channels:
            - /Dev1/AI1
        sample_scanner_counter_channels:
            - /Dev1/Ctr1
        tip_scanner_ai_channels:
            - /Dev1/AI2
        tip_scanner_counter_channels:
            - /Dev1/Ctr2
        sample_scanner_voltage_ranges:
            - [0, 10]
            - [0, 10]
        tip_scanner_voltage_ranges:
            - [0, 10]
            - [0, 10]
        sample_scanner_position_ranges:
            - [0e-6, 40e-6]
            - [0e-6, 40e-6]
        tip_scanner_position_ranges:
            - [0e-6, 40e-6]
            - [0e-6, 40e-6]
        default_samples_number: 50
        max_counts: 3e7
        read_write_timeout: 10
        counting_edge_rising: True

    """

    # config options
    _photon_sources = ConfigOption('photon_sources', list(), missing='warn')

    #Clock used when moving the scanners from A to B, without acquiring anything along the way
    _motion_clock_channel = ConfigOption('motion_clock_channel', missing='error')
    _default_motion_clock_frequency = ConfigOption('default_motion_clock_frequency', 100, missing='info')

    # Photon counting settings
    _counter_clock = ConfigOption('counter_clock', '100kHzTimebase', missing='info')
    _counter_clock_frequency = ConfigOption('counter_clock_frequency', 100e3, missing='info')
    _counter_channels = ConfigOption('counter_channels', missing='error')
    _counter_ai_channels = ConfigOption('counter_ai_channels', list(), missing='info')
    _counter_voltage_range = ConfigOption('counter_voltage_range', [-10, 10], missing='info')

    # Sample scanner
    _sample_scanner_ao_channels = ConfigOption('sample_scanner_ao_channels', missing='error')
    _sample_scanner_ai_channels = ConfigOption('sample_scanner_ai_channels', list(), missing='info')
    _sample_scanner_counter_channels = ConfigOption('sample_scanner_counter_channels', list(), missing='warn')
    _sample_scanner_voltage_ranges = ConfigOption('sample_scanner_voltage_ranges', missing='error')
    _sample_scanner_position_ranges = ConfigOption('sample_scanner_position_ranges', missing='error')

    # Tip scanner
    _tip_scanner_ao_channels = ConfigOption('tip_scanner_ao_channels', missing='error')
    _tip_scanner_ai_channels = ConfigOption('tip_scanner_ai_channels', list(), missing='info')
    _tip_scanner_counter_channels = ConfigOption('tip_scanner_counter_channels', list(), missing='warn')
    _tip_scanner_voltage_ranges = ConfigOption('tip_scanner_voltage_ranges', missing='error')
    _tip_scanner_position_ranges = ConfigOption('tip_scanner_position_ranges', missing='error')


    # used as a default for expected maximum counts
    _max_counts = ConfigOption('max_counts', default=3e7)
    # timeout for the Read or/and write process in s
    _RWTimeout = ConfigOption('read_write_timeout', default=10)

    _sample_stack_name = ConfigOption('Sample stack name', default='sample')
    _tip_stack_name = ConfigOption('Tip stack names', default='tip')

    def on_activate(self):
        """ Starts up the NI Card at activation.
        """
        self._stack_names = [self._sample_stack_name, self._tip_stack_name]

        #Regroup the initial parameters for the scanner into dictionaries for ease of access
        self._scanner_ao_channels = dict(zip(self._stack_names,
                                             [self._sample_scanner_ao_channels,
                                              self._tip_scanner_ao_channels]))
        self._scanner_ai_channels = dict(zip(self._stack_names,
                                             [self._sample_scanner_ai_channels,
                                              self._tip_scanner_ai_channels]))
        self._scanner_counter_channels = dict(zip(self._stack_names,
                                              [self._sample_scanner_counter_channels,
                                               self._tip_scanner_counter_channels]))
        self._scanner_voltage_ranges = dict(zip(self._stack_names,
                                                [self._sample_scanner_voltage_ranges,
                                                 self._tip_scanner_voltage_ranges]))

        self._scanner_position_ranges = dict(zip(self._stack_names,
                                                 [self._sample_scanner_position_ranges,
                                                  self._tip_scanner_position_ranges]))

        # the tasks used on that hardware device:
        self._counter_daq_tasks = list()
        self._counter_ai_daq_task = list()

        self._scanner_ao_tasks = dict().fromkeys(self._stack_names, None)
        self._motion_clock_task = None

        self._photon_sources = self._photon_sources if self._photon_sources is not None else list()

        if len(self._sample_scanner_ao_channels) < len(self._sample_scanner_voltage_ranges):
            self.log.error(
                'Specify at least as many sample_scanner_voltage_ranges as sample_scanner_ao_channels!')

        if len(self._tip_scanner_ao_channels) < len(self._tip_scanner_voltage_ranges):
            self.log.error(
                'Specify at least as many tip_scanner_voltage_ranges as tip_scanner_ao_channels!')

        if len(self._sample_scanner_ao_channels) < len(self._sample_scanner_position_ranges):
            self.log.error(
                'Specify at least as many sample_scanner_position_ranges as sample_scanner_ao_channels!')

        if len(self._tip_scanner_ao_channels) < len(self._tip_scanner_position_ranges):
            self.log.error(
                'Specify at least as many tip_scanner_position_ranges as tip_scanner_ao_channels!')
        """
        if len(self._scanner_counter_channels) + len(self._scanner_ai_channels) < 1:
            self.log.error(
                'Specify at least one counter or analog input channel for the scanner!')
        """
        # Analog output is always needed and it does not interfere with the
        # rest, so start it always and leave it running
        if self._start_analog_outputs() < 0:
            self.log.error('Failed to start analog output.')
            raise Exception('Failed to start NI Card module due to analog output failure.')
        self._current_position = [0,0]

    def on_deactivate(self):
        """ Shut down the NI card.
        """
        self._stop_analog_outputs()
        # clear the task
        try:
            for stack_name, scanner_ao_task in self._scanner_ao_tasks.items():
                daq.DAQmxClearTask(scanner_ao_task)
                self._scanner_ao_tasks[stack_name] = None
        except:
            self.log.exception('Could not clear AO Out Task.')

        self.reset_hardware()

    def prepare_counters(self,
                         counter_channels=None,
                         counter_ai_channels=None,
                         sources=None,
                         counter_clock=None,
                         samples_to_acquire=None):
        """ Configures the actual counter with a given clock.

        @param list(str) counter_channels: optional, physical channel of the counter
        @param list(str) sources: optional, physical channel where the photons
                                  are to count from
        @param str clock_channel: optional, specifies the clock channel for the
                                  counter
        @param int counter_buffer: optional, a buffer of specified integer
                                   length, where in each bin the count numbers
                                   are saved.

        @return int: error code (0:OK, -1:error)
        """

        my_counter_channels = counter_channels if counter_channels else self._counter_channels
        my_photon_sources = sources if sources else self._photon_sources
        my_clock_channel = counter_clock if counter_clock else self._counter_clock
        #If no AI is specified, then create an empty array (and do not create AI tasks)
        my_counter_ai_channels = counter_ai_channels# if counter_ai_channels else []
        self._counter_ai_channels = my_counter_ai_channels

        if len(my_photon_sources) < len(my_counter_channels):
            self.log.error('You have given {0} sources but {1} counting channels.'
                           'Please give an equal or greater number of sources.'
                           ''.format(len(my_photon_sources), len(my_counter_channels)))
            return -1
        else:
            counter_num = len(my_counter_channels)

        try:
            for i, counter_chan, ph_source in zip(range(counter_num),
                                                  my_counter_channels,
                                                  my_photon_sources):

                task = daq.TaskHandle()

                daq.DAQmxCreateTask('Counter{0}'.format(i), daq.byref(task))

                daq.DAQmxCreateCICountEdgesChan(
                    task,
                    counter_chan,
                    'Counter Channel {0}'.format(i),
                    daq.DAQmx_Val_Rising, #The trigger is a rising edge
                    0,
                    daq.DAQmx_Val_CountUp, #When edge detected, add one count
                    )

                daq.DAQmxCfgSampClkTiming(task,
                                          my_clock_channel,
                                          self._counter_clock_frequency,
                                          daq.DAQmx_Val_Rising,
                                          daq.DAQmx_Val_FiniteSamps,
                                          samples_to_acquire
                                          )

                #Connect the counter to the physical terminal
                daq.DAQmxSetCICountEdgesTerm(task,
                                             counter_chan,
                                             ph_source
                                             )

                # add task to counter task list
                self._counter_daq_tasks.append(task)

                # Counter analog input task
                if my_counter_ai_channels is not None and len(my_counter_ai_channels) > 0:
                    atask = daq.TaskHandle()

                    daq.DAQmxCreateTask('CounterAnalogIn', daq.byref(atask))

                    daq.DAQmxCreateAIVoltageChan(
                        atask,
                        ', '.join(my_counter_ai_channels),
                        'Counter Analog In',
                        daq.DAQmx_Val_RSE,
                        self._counter_voltage_range[0],
                        self._counter_voltage_range[1],
                        daq.DAQmx_Val_Volts,
                        ''
                    )
                    # Analog in channel timebase
                    daq.DAQmxCfgSampClkTiming(
                        atask,
                        my_clock_channel,
                        self._counter_clock_frequency,
                        daq.DAQmx_Val_Rising,
                        daq.DAQmx_Val_FiniteSamps,
                        samples_to_acquire
                    )
                    self._counter_ai_daq_task = atask
                else:
                    self._counter_ai_daq_task = None

        except:
            self.log.exception('Error while setting up counting task.')
            return -1

        return 0

    def close_counters(self):
        """

        @param bool scanner: specifies if the counter- or scanner- function
                             will be excecuted to close the device.
                                True = scanner
                                False = counter

        @return int: error code (0:OK, -1:error)
        """
        error = 0

        for i, task in enumerate(self._counter_daq_tasks):
            try:
                # stop the counter task
                daq.DAQmxStopTask(task)
                # after stopping delete all the configuration of the counter
                daq.DAQmxClearTask(task)
                # set the task handle to None as a safety
            except:
                self.log.exception('Could not close counter.')
                error = -1
        self._counter_daq_tasks = []

        if self._counter_ai_channels is not None and len(self._counter_ai_channels) > 0:
            try:
                # stop the counter task
                daq.DAQmxStopTask(self._counter_ai_daq_task)
                daq.DAQmxClearTask(self._counter_ai_daq_task)
            except:
                self.log.exception('Could not close counter analog channels.')
                error = -1
            self._counter_ai_daq_task = None
            self._counter_ai_channels = []
        return error

    def prepare_motion_clock(self, clock_frequency = None, clock_channel = None):

        if self._motion_clock_task is not None:
            self.log.error("The motion clock is running, close the other task first.")
            return -1

        clock_task = daq.TaskHandle()
        if clock_frequency is None:
            clock_frequency = self._default_motion_clock_frequency

        if clock_channel is None:
            clock_channel = self._motion_clock_channel

        try:
            task_name = 'MotionClockTask'
            daq.DAQmxCreateTask(task_name, daq.byref(clock_task))

            daq.DAQmxCreateCOPulseChanFreq(
                clock_task,
                clock_channel,
                'Motion clock',
                daq.DAQmx_Val_Hz,
                daq.DAQmx_Val_Low,
                0,
                clock_frequency,
                0.5
            )

            daq.DAQmxCfgImplicitTiming(
                clock_task,
                daq.DAQmx_Val_ContSamps,
                1000)

        except:
            self.log.exception("Error while setting up the motion clock.")
            return -1
        return 0

    def get_counter_channels(self):
        """ Returns the list of counter channel names.

        @return tuple(str): channel names

        Most methods calling this might just care about the number of channels, though.
        """
        ch = self._counter_channels[:]
        return ch

    def get_ai_counter_channels(self, stack_name=None):
        return self._scanner_ai_channels[stack_name]

    # ================ ConfocalScannerInterface Commands =======================
    def _start_analog_outputs(self):
        """ Starts or restarts the analog output.

        @return int: error code (0:OK, -1:error)
        """
        try:
            # If an analog task is already running, kill that one first
            for stack_name, stack_scanner_ao_task in self._scanner_ao_tasks.items():
                if stack_scanner_ao_task is not None:
                    # stop the analog output task
                    daq.DAQmxStopTask(stack_scanner_ao_task)

                    # delete the configuration of the analog output
                    daq.DAQmxClearTask(stack_scanner_ao_task)

                    # set the task handle to None as a safety
                    self._scanner_ao_tasks[stack_name] = None

            # initialize ao channels / task for scanner, should always be active.
            # Define at first the type of the variable as a Task:
            for stack_name in self._scanner_ao_tasks:
                self._scanner_ao_tasks[stack_name] = daq.TaskHandle()

            # create the actual analog output task on the hardware device. Via
            # byref you pass the pointer of the object to the TaskCreation function:
            for stack_name, scanner_ao_task in self._scanner_ao_tasks.items():
                daq.DAQmxCreateTask('{:s}ScannerAO'.format(stack_name.capitalize()),
                                    daq.byref(scanner_ao_task))
                for n, chan in enumerate(self._scanner_ao_channels[stack_name]):
                    # Assign and configure the created task to an analog output voltage channel.
                    daq.DAQmxCreateAOVoltageChan(
                        # The AO voltage operation function is assigned to this task.
                        scanner_ao_task,
                        # use (all) scanner ao_channels for the output
                        chan,
                        # assign a name for that channel
                        'Sample Scanner AO Channel {0}'.format(n),
                        # minimum possible voltage
                        self._scanner_voltage_ranges[stack_name][n][0],
                        # maximum possible voltage
                        self._scanner_voltage_ranges[stack_name][n][1],
                        # units is Volt
                        daq.DAQmx_Val_Volts,
                        # empty for future use
                        '')

                daq.DAQmxSetSampTimingType(scanner_ao_task, daq.DAQmx_Val_OnDemand)

        except:
            self.log.exception('Error starting analog output task.')
            return -1
        return 0

    def _stop_analog_outputs(self):
        """ Stops the analog output.

        @return int: error code (0:OK, -1:error)
        """
        for scanner_ao_task in self._scanner_ao_tasks.values():
            if scanner_ao_task is None:
                return -1
        retval = 0
        try:
            # stop the analog output task
            for scanner_ao_task in self._scanner_ao_tasks.values():
                daq.DAQmxStopTask(scanner_ao_task)
        except:
            self.log.exception('Error stopping analog outputs.')
            retval = -1
        try:
            for scanner_ao_task in self._scanner_ao_tasks.values():
                daq.DAQmxSetSampTimingType(scanner_ao_task, daq.DAQmx_Val_OnDemand)
        except:
            self.log.exception('Error changing analog outputs mode.')
            retval = -1
        return retval

    def get_counter_clock_frequency(self):
        return self._counter_clock_frequency

    def get_scanner_position(self):
        """ Get the current position of the scanner hardware.

        @return float[]: current position in (x, y, z, a).
        """
        return self._current_position.tolist()

    def get_position_range(self, stack=None):
        """ Returns the physical range of the scanner.

        @return float [4][2]: array of 4 ranges with an array containing lower
                              and upper limit. The unit of the scan range is
                              meters.
        """

        return self._scanner_position_ranges[stack]

    def get_voltage_range(self, stack=None):
        """ Returns the physical range of the scanner.

        @return float [4][2]: array of 4 ranges with an array containing lower
                              and upper limit. The unit of the scan range is
                              meters.
        """

        return self._scanner_voltage_ranges[stack]

    def get_stack_names(self):
        """
        Return the names of the stacks, ordered as sample, tip.
        """
        return self._sample_stack_name, self._tip_stack_name

    def get_scanner_count_channels(self):
        """ Return list of counter channels """
        ch = self._scanner_counter_channels[:]
        #ch.extend(self._scanner_ai_channels)
        return ch

    def read_pixel(self, samples=1):
        for i, task in enumerate(self._counter_daq_tasks):
            daq.DAQmxStartTask(task)

        if self._counter_ai_daq_task:
            daq.DAQmxStartTask(self._counter_ai_daq_task)

        count_matrix = np.full((len(self._sample_scanner_counter_channels),
                                samples), 222, dtype=np.uint32)

        final_counts = np.full(len(self._counter_daq_tasks), 222)

        read_samples = daq.int32()
        try:
            for i, task in enumerate(self._counter_daq_tasks):
                daq.DAQmxReadCounterU32(
                    task,
                    samples,
                    self._RWTimeout,
                    count_matrix[i],
                    samples,
                    daq.byref(read_samples),
                    None
                )
                final_counts[i] = count_matrix[i, -1]
                daq.DAQmxStopTask(task)
        except:
            self.log.exception("Failed reading counter samples")

        # FIXME: this part is incorrect: check the ai_channels are the correct ones. Also, the readout needs to be a for
        #  loop, and the data should be stored in each row of the analog_data container.
        if self._counter_ai_daq_task is not None:
            try:
                analog_data = np.full(
                    (len(self._counter_ai_channels), samples),
                    111, dtype=np.float64)

                analog_read_samples = daq.int32()

                daq.DAQmxReadAnalogF64(
                    self._counter_ai_daq_task,
                    samples,
                    self._RWTimeout,
                    daq.DAQmx_Val_GroupByChannel,
                    analog_data,
                    len(self._counter_ai_channels) * samples,
                    daq.byref(analog_read_samples),
                    None
                )

                daq.DAQmxStopTask(self._counter_ai_daq_task)
            except:
                self.log.exception("Failed reading analog samples.")

            analog_data = analog_data.mean(axis=0)
        else:
            analog_data = None

        return final_counts, analog_data


    def write_voltage(self, ao_task, voltages=[], autostart=True):
        AONwritten = daq.int32()
        if ao_task is None:
            self.log.exception('Please specify or create an AO task')

        length = len(voltages)
        try:
            daq.DAQmxWriteAnalogF64(ao_task,
                                    1,
                                    autostart,
                                    self._RWTimeout,
                                    daq.DAQmx_Val_GroupByChannel,
                                    voltages,
                                    daq.byref(AONwritten),
                                    None)
        except:
            self.log.exception('Error writing the analog output to the channels.')
        return AONwritten.value

    def scanner_set_position(self, xypair=None, stack=None):
        """ Move stage to x, y.

        @param float x: position in x-direction (in axis unit)
        @param float y: position in y-direction (in axis unit)

        @return int: error code (0:OK, -1:error)
        """
        x, y = xypair
        if stack is None:
            self.log.error('The scanning stack has not been specified.')
            return -1

        # if self.module_state() == 'locked':
        #     self.log.error('Another scan is already running, close it first.')
        #     return -1

        scanning_task = self._scanner_ao_tasks[stack]

        scanner_position_range = self._scanner_position_ranges[stack]

        if x is not None:
            if not (scanner_position_range[0][0] <= x <= scanner_position_range[0][1]):
                self.log.error('You want to set x out of range: {0:f}.'.format(x))
                return -1
            #self._current_position[0] = np.float(x)

        if y is not None:
            if not (scanner_position_range[1][0] <= y <= scanner_position_range[1][1]):
                self.log.error('You want to set y out of range: {0:f}.'.format(y))
                return -1
            #self._current_position[1] = np.float(y)

        my_position = np.array(xypair)
        # then directly write the position to the hardware
        voltages = self._scanner_position_to_volt(my_position[np.newaxis, :], stack=stack)
        try:
            self.write_voltage(
                scanning_task,
                voltages=voltages,
                autostart=True)
        except:
            return -1
        return 0

    def close_scanner(self, stack=None):
        """ Closes the scanners and cleans up afterwards.

        @return int: error code (0:OK, -1:error)
        """
        a = self._stop_analog_outputs()

        b = 0
        if self._scanner_ai_channels[stack]:
            try:
                # stop the counter task
                daq.DAQmxStopTask(self._counter_ai_daq_task)
                # after stopping delete all the configuration of the counter
                daq.DAQmxClearTask(self._counter_ai_daq_task)
                # set the task handle to None as a safety
                self._scanner_analog_daq_task = None
            except:
                self.log.exception('Could not close analog.')
                b = -1

        c = self.close_counter(scanner=True)
        return -1 if a < 0 or b < 0 or c < 0 else 0

    def close_scanner_clock(self):
        """ Closes the clock and cleans up afterwards.

        @return int: error code (0:OK, -1:error)
        """
        return self.close_clock(scanner=True)

    def reset_hardware(self):
        #FIXME: the resetting is still not working properly.
        """ Resets the NI hardware, so the connection is lost and other
            programs can access it.

        @return int: error code (0:OK, -1:error)
        """
        retval = 0
        chanlist = []
        chanlist.extend(self._scanner_ao_channels)
        chanlist.extend(self._photon_sources)
        chanlist.extend(self._counter_channels)
        chanlist.extend(self._scanner_counter_channels)

        devicelist = []
        for channel in chanlist:
            if channel is None:
                continue
            match = re.match(
                '^/(?P<dev>[0-9A-Za-z\- ]+[0-9A-Za-z\-_ ]*)/(?P<chan>[0-9A-Za-z]+)',
                channel)
            if match:
                devicelist.append(match.group('dev'))
            else:
                self.log.error('Did not find device name in {0}.'.format(channel))
        for device in set(devicelist):
            self.log.info('Reset device {0}.'.format(device))
            try:
                daq.DAQmxResetDevice(device)
            except:
                self.log.exception('Could not reset NI device {0}'.format(device))
                retval = -1
        return retval

    def _scanner_position_to_volt(self, positions=None, stack=None):
        """ Converts a set of position pixels to actual voltages.

        @param np.ndarray positions: array with (n, 2) shape, first column corresponding to x and second one to y.
        @param string stack: the stack name

        @return float[][n]: array of n-part tuples of corresponing voltages
        """
        if not isinstance(positions, np.ndarray):
            self.log.error('Given positions are not and nd array.')
            return np.array([np.NaN])

        if stack is None:
            self.log.error('Please specify a stack to work with.')

        scanner_voltage_ranges = np.array(self._scanner_voltage_ranges[stack])
        scanner_position_ranges = np.array(self._scanner_position_ranges[stack])

        post2volt_coeff = np.diff(scanner_voltage_ranges, axis=1) / np.diff(scanner_position_ranges, axis=1)
        volts = post2volt_coeff.T * positions + scanner_voltage_ranges[:, 0].T

        if (np.any(np.logical_or(volts[:, 0] < scanner_voltage_ranges[0, 0],
                                 volts[:, 0] > scanner_voltage_ranges[0, 1]))
                or
            np.any(np.logical_or(volts[:, 1] < scanner_voltage_ranges[1, 0],
                                 volts[:, 1] > scanner_voltage_ranges[1, 1]))):
            self.log.error('Some of the requested positions are out of the voltage bounds.')
            return np.array([np.nan])
        return volts.astype(np.float64)

    # ================ End ConfocalScannerInterface Commands ===================

    def get_status(self):
        """ Receives the current status of the Fast Counter and outputs it as
            return value.

        0 = unconfigured
        1 = idle
        2 = running
        3 = paused
        -1 = error state
        """
        if self._gated_counter_daq_task is None:
            return 0
        else:
            # return value represents a uint32 value, i.e.
            #   task_done = 0  => False, i.e. device is runnin
            #   task_done !=0  => True, i.e. device has stopped
            task_done = daq.bool32()

            ret_v = daq.DAQmxIsTaskDone(
                # task reference
                self._gated_counter_daq_task,
                # reference to bool value.
                daq.byref(task_done))

            if ret_v != 0:
                return ret_v

            if task_done.value() == 0:
                return 1
            else:
                return 2

    # ======================== Digital channel control ==========================

    def digital_channel_switch(self, channel_name, mode=True):
        """
        Switches on or off the voltage output (5V) of one of the digital channels, that
        can as an example be used to switch on or off the AOM driver or apply a single
        trigger for ODMR.
        @param str channel_name: Name of the channel which should be controlled
                                    for example ('/Dev1/PFI9')
        @param bool mode: specifies if the voltage output of the chosen channel should be turned on or off

        @return int: error code (0:OK, -1:error)
        """
        if channel_name is None:
            self.log.error('No channel for digital output specified')
            return -1
        else:

            self.digital_out_task = daq.TaskHandle()
            if mode:
                self.digital_data = daq.c_uint32(0xffffffff)
            else:
                self.digital_data = daq.c_uint32(0x0)
            self.digital_read = daq.c_int32()
            self.digital_samples_channel = daq.c_int32(1)
            daq.DAQmxCreateTask('DigitalOut', daq.byref(self.digital_out_task))
            daq.DAQmxCreateDOChan(self.digital_out_task, channel_name, "", daq.DAQmx_Val_ChanForAllLines)
            daq.DAQmxStartTask(self.digital_out_task)
            daq.DAQmxWriteDigitalU32(self.digital_out_task, self.digital_samples_channel, True,
                                        self._RWTimeout, daq.DAQmx_Val_GroupByChannel,
                                        np.array(self.digital_data), self.digital_read, None)

            daq.DAQmxStopTask(self.digital_out_task)
            daq.DAQmxClearTask(self.digital_out_task)
            return 0


