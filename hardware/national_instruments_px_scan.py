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
        my_counter_ai_channels = counter_ai_channels if counter_ai_channels else []

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
                if len(my_counter_ai_channels) > 0:
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

        if len(self._sample_scanner_ai_channels) > 0:
            try:
                # stop the counter task
                daq.DAQmxStopTask(self._counter_ai_daq_task)
                # after stopping delete all the configuration of the counter
                daq.DAQmxClearTask(self._counter_ai_daq_task)
                # set the task handle to None as a safety
            except:
                self.log.exception('Could not close counter analog channels.')
                error = -1
            self._counter_ai_daq_task = None
        return error

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

        final_counts = np.full(samples, 222)

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
        if len(self._sample_scanner_ai_channels) > 0:
            try:
                analog_data = np.full(
                    (len(self._sample_scanner_ai_channels), samples),
                    111, dtype=np.float64)

                analog_read_samples = daq.int32()

                daq.DAQmxReadAnalogF64(
                    self._counter_ai_daq_task,
                    samples,
                    self._RWTimeout,
                    daq.DAQmx_Val_GroupByChannel,
                    analog_data,
                    len(self._sample_scanner_ai_channels) * samples,
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

        if self.module_state() == 'locked':
            self.log.error('Another scan is already running, close it first.')
            return -1

        scanning_task = self._scanner_ao_tasks[stack]

        scanner_position_range = self._scanner_position_ranges[stack]

        if x is not None:
            if not (scanner_position_range[0][0] <= x <= scanner_position_range[0][1]):
                self.log.error('You want to set x out of range: {0:f}.'.format(x))
                return -1
            self._current_position[0] = np.float(x)

        if y is not None:
            if not (scanner_position_range[1][0] <= y <= scanner_position_range[1][1]):
                self.log.error('You want to set y out of range: {0:f}.'.format(y))
                return -1
            self._current_position[1] = np.float(y)

        # the position has to be a vstack
        my_position = np.vstack(self._current_position)

        # then directly write the position to the hardware
        try:
            self.write_voltage(
                scanning_task,
                voltages=self._scanner_position_to_volt(my_position),
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
        """ Resets the NI hardware, so the connection is lost and other
            programs can access it.

        @return int: error code (0:OK, -1:error)
        """
        retval = 0
        chanlist = [
            self._clock_channel,
            self._scanner_clock_channel,
            ]
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

    # ==================== ODMRCounterInterface Commands =======================
    def set_up_odmr_clock(self, clock_frequency=None, clock_channel=None):
        """ Configures the hardware clock of the NiDAQ card to give the timing.

        @param float clock_frequency: if defined, this sets the frequency of
                                      the clock
        @param string clock_channel: if defined, this is the physical channel
                                     of the clock

        @return int: error code (0:OK, -1:error)
        """

        return self.set_up_clock(
            clock_frequency=clock_frequency,
            clock_channel=clock_channel,
            scanner=True,
            idle=False)

    def set_up_odmr(self, counter_channel=None, photon_source=None,
                    clock_channel=None, odmr_trigger_channel=None):
        """ Configures the actual counter with a given clock.

        @param string counter_channel: if defined, this is the physical channel
                                       of the counter
        @param string photon_source: if defined, this is the physical channel
                                     where the photons are to count from
        @param string clock_channel: if defined, this specifies the clock for
                                     the counter
        @param string odmr_trigger_channel: if defined, this specifies the
                                            trigger output for the microwave

        @return int: error code (0:OK, -1:error)
        """
        if self._scanner_clock_daq_task is None and clock_channel is None:
            self.log.error('No clock running, call set_up_clock before starting the counter.')
            return -1
        if self._scanner_counter_daq_tasks:
            self.log.error('Another counter is already running, close this one first.')
            return -1
        if self._scanner_ai_channels and self._scanner_analog_daq_task is not None:
            self.log.error('Another analog is already running, close this one first.')
            return -1

        my_clock_channel = clock_channel if clock_channel else self._scanner_clock_channel

        if self._scanner_counter_channels and self._photon_sources:
            my_counter_channel = counter_channel if counter_channel else self._scanner_counter_channels[0]
            my_photon_source = photon_source if photon_source else self._photon_sources[0]

            # this task will count photons with binning defined by the clock_channel
            task = daq.TaskHandle()
            try:
                # create task for the counter
                daq.DAQmxCreateTask('ODMRCounter', daq.byref(task))

                # set up semi period width measurement in photon ticks, i.e. the width
                # of each pulse (high and low) generated by pulse_out_task is measured
                # in photon ticks.
                #   (this task creates a channel to measure the time between state
                #    transitions of a digital signal and adds the channel to the task
                #    you choose)
                daq.DAQmxCreateCISemiPeriodChan(
                    # define to which task to# connect this function
                    task,
                    # use this counter channel
                    my_counter_channel,
                    # name to assign to it
                    'ODMR Counter',
                    # Expected minimum count value
                    0,
                    # Expected maximum count value
                    self._max_counts / self._scanner_clock_frequency,
                    # units of width measurement, here photon ticks
                    daq.DAQmx_Val_Ticks,
                    '')

                # connect the pulses from the clock to the counter
                daq.DAQmxSetCISemiPeriodTerm(
                    task,
                    my_counter_channel,
                    my_clock_channel + 'InternalOutput')

                # define the source of ticks for the counter as self._photon_source
                daq.DAQmxSetCICtrTimebaseSrc(
                    task,
                    my_counter_channel,
                    my_photon_source)

                self._scanner_counter_daq_tasks.append(task)
            except:
                self.log.exception('Error while setting up the digital counter of ODMR scan.')
                return -1

        try:
            # Analog task
            if self._scanner_ai_channels:
                atask = daq.TaskHandle()
                daq.DAQmxCreateTask('ODMRAnalog', daq.byref(atask))

                daq.DAQmxCreateAIVoltageChan(
                    atask,
                    ', '.join(self._scanner_ai_channels),
                    'ODMR Analog',
                    daq.DAQmx_Val_RSE,
                    -10,
                    10,
                    daq.DAQmx_Val_Volts,
                    ''
                )
                self._scanner_analog_daq_task = atask

            # start and stop pulse task to correctly initiate idle state high voltage.
            daq.DAQmxStartTask(self._scanner_clock_daq_task)
            # otherwise, it will be low until task starts, and MW will receive wrong pulses.
            daq.DAQmxStopTask(self._scanner_clock_daq_task)

            if self.lock_in_active:
                ptask = daq.TaskHandle()
                daq.DAQmxCreateTask('ODMRPulser', daq.byref(ptask))
                daq.DAQmxCreateDOChan(
                    ptask,
                    '{0:s}, {1:s}'.format(self._odmr_trigger_line, self._odmr_switch_line),
                    "ODMRPulserChannel",
                    daq.DAQmx_Val_ChanForAllLines)

                self._odmr_pulser_daq_task = ptask

            # connect the clock to the trigger channel to give triggers for the
            # microwave
            daq.DAQmxConnectTerms(
                self._scanner_clock_channel + 'InternalOutput',
                self._odmr_trigger_channel,
                daq.DAQmx_Val_DoNotInvertPolarity)
        except:
            self.log.exception('Error while setting up ODMR scan.')
            return -1
        return 0

    def set_odmr_length(self, length=100):
        """ Sets up the trigger sequence for the ODMR and the triggered microwave.

        @param int length: length of microwave sweep in pixel

        @return int: error code (0:OK, -1:error)
        """
        if self._scanner_counter_channels and len(self._scanner_counter_daq_tasks) < 1:
            self.log.error('No counter is running, cannot do ODMR without one.')
            return -1

        if self._scanner_ai_channels and self._scanner_analog_daq_task is None:
            self.log.error('No analog task is running, cannot do ODMR without one.')
            return -1

        self._odmr_length = length
        try:
            # set timing for odmr clock task to the number of pixel.
            daq.DAQmxCfgImplicitTiming(
                # define task
                self._scanner_clock_daq_task,
                # only a limited number of counts
                daq.DAQmx_Val_FiniteSamps,
                # count twice for each voltage +1 for starting this task.
                # This first pulse will start the count task.
                self._odmr_length + 1)

            # Digital
            if self._scanner_counter_channels:
                # set timing for odmr count task to the number of pixel.
                daq.DAQmxCfgImplicitTiming(
                    # define task
                    self._scanner_counter_daq_tasks[0],
                    # only a limited number of counts
                    daq.DAQmx_Val_ContSamps,
                    # count twice for each voltage +1 for starting this task.
                    # This first pulse will start the count task.
                    2 * (self._odmr_length + 1))

                # read samples from beginning of acquisition, do not overwrite
                daq.DAQmxSetReadRelativeTo(
                    self._scanner_counter_daq_tasks[0],
                    daq.DAQmx_Val_CurrReadPos)

                # do not read first sample
                daq.DAQmxSetReadOffset(
                    self._scanner_counter_daq_tasks[0],
                    0)

                # unread data in buffer will be overwritten
                daq.DAQmxSetReadOverWrite(
                    self._scanner_counter_daq_tasks[0],
                    daq.DAQmx_Val_DoNotOverwriteUnreadSamps)

            # Analog
            if self._scanner_ai_channels:
                # Analog in channel timebase
                daq.DAQmxCfgSampClkTiming(
                    self._scanner_analog_daq_task,
                    self._scanner_clock_channel + 'InternalOutput',
                    self._scanner_clock_frequency,
                    daq.DAQmx_Val_Rising,
                    daq.DAQmx_Val_ContSamps,
                    self._odmr_length + 1
                )

            if self._odmr_pulser_daq_task:
                # pulser channel timebase
                daq.DAQmxCfgSampClkTiming(
                    self._odmr_pulser_daq_task,
                    self._scanner_clock_channel + 'InternalOutput',
                    self._scanner_clock_frequency,
                    daq.DAQmx_Val_Rising,
                    daq.DAQmx_Val_ContSamps,
                    self._odmr_length + 1
                )
        except:
            self.log.exception('Error while setting up ODMR counter.')
            return -1
        return 0

    @property
    def oversampling(self):
        return self._oversampling

    @oversampling.setter
    def oversampling(self, val):
        if not isinstance(val, (int, float)):
            self.log.error('oversampling has to be int of float.')
        else:
            self._oversampling = int(val)

    @property
    def lock_in_active(self):
        return self._lock_in_active

    @lock_in_active.setter
    def lock_in_active(self, val):
        if not isinstance(val, bool):
            self.log.error('lock_in_active has to be boolean.')
        else:
            self._lock_in_active = val
            if self._lock_in_active:
                self.log.warn('You just switched the ODMR counter to Lock-In-mode. \n'
                              'Please make sure you connected all triggers correctly:\n'
                              '  {0:s} is the microwave trigger channel\n'
                              '  {1:s} is the switching channel for the lock in\n'
                              ''.format(self._odmr_trigger_line, self._odmr_switch_line))

    def count_odmr(self, length=100):
        """ Sweeps the microwave and returns the counts on that sweep.

        @param int length: length of microwave sweep in pixel

        @return float[]: the photon counts per second
        """
        if len(self._scanner_counter_daq_tasks) < 1 and self._scanner_counter_channels:
            self.log.error(
                'No counter is running, cannot scan an ODMR line without one.')
            return True, np.array([-1.])

        if self._scanner_ai_channels and self._scanner_analog_daq_task is None:
            self.log.error('No analog task is running, cannot do ODMR without one.')
            return True, np.array([-1.])

        # check if length setup is correct, if not, adjust.
        if self._odmr_pulser_daq_task:
            odmr_length_to_set = length * self.oversampling * 2
        else:
            odmr_length_to_set = length

        if self.set_odmr_length(odmr_length_to_set) < 0:
            self.log.error('An error arose while setting the odmr lenth to {}.'.format(odmr_length_to_set))
            return True, np.array([-1.])

        try:
            # start the scanner counting task that acquires counts synchronously
            if self._scanner_counter_channels:
                daq.DAQmxStartTask(self._scanner_counter_daq_tasks[0])
            if self._scanner_ai_channels:
                daq.DAQmxStartTask(self._scanner_analog_daq_task)
        except:
            self.log.exception('Cannot start ODMR counter.')
            return True, np.array([-1.])

        if self._odmr_pulser_daq_task:
            try:

                # The pulse pattern is an alternating 0 and 1 on the switching channel (line0),
                # while the first half of the whole microwave pulse is 1 and the other half is 0.
                # This way the beginning of the microwave has a rising edge.
                pulse_pattern = np.zeros(self.oversampling * 2, dtype=np.uint32)
                pulse_pattern[:self.oversampling] += 1
                pulse_pattern[::2] += 2

                daq.DAQmxWriteDigitalU32(self._odmr_pulser_daq_task,
                                         len(pulse_pattern),
                                         0,
                                         self._RWTimeout * self._odmr_length,
                                         daq.DAQmx_Val_GroupByChannel,
                                         pulse_pattern,
                                         None,
                                         None)

                daq.DAQmxStartTask(self._odmr_pulser_daq_task)
            except:
                self.log.exception('Cannot start ODMR pulser.')
                return True, np.array([-1.])

        try:
            daq.DAQmxStartTask(self._scanner_clock_daq_task)

            # wait for the scanner clock to finish
            daq.DAQmxWaitUntilTaskDone(
                # define task
                self._scanner_clock_daq_task,
                # maximal timeout for the counter times the positions
                self._RWTimeout * 2 * self._odmr_length)

            # Digital
            if self._scanner_counter_channels:
                # count data will be written here
                odmr_data = np.full(
                    (2 * self._odmr_length + 1, ),
                    222,
                    dtype=np.uint32)

                #number of samples which were read will be stored here
                n_read_samples = daq.int32()

                # actually read the counted photons
                daq.DAQmxReadCounterU32(
                    # read from this task
                    self._scanner_counter_daq_tasks[0],
                    # Read number of double the# number of samples
                    2 * self._odmr_length + 1,
                    # Maximal timeout for the read # process
                    self._RWTimeout,
                    # write into this array
                    odmr_data,
                    # length of array to write into
                    2 * self._odmr_length + 1,
                    # number of samples which were actually read
                    daq.byref(n_read_samples),
                    # Reserved for future use. Pass NULL (here None) to this parameter.
                    None)

            # Analog
            if self._scanner_ai_channels:
                odmr_analog_data = np.full(
                    (len(self._scanner_ai_channels), self._odmr_length + 1),
                    222,
                    dtype=np.float64)

                analog_read_samples = daq.int32()

                daq.DAQmxReadAnalogF64(
                    self._scanner_analog_daq_task,
                    self._odmr_length + 1,
                    self._RWTimeout,
                    daq.DAQmx_Val_GroupByChannel,
                    odmr_analog_data,
                    len(self._scanner_ai_channels) * (self._odmr_length + 1),
                    daq.byref(analog_read_samples),
                    None
                )

            # stop the counter task
            daq.DAQmxStopTask(self._scanner_clock_daq_task)
            if self._scanner_counter_channels:
                daq.DAQmxStopTask(self._scanner_counter_daq_tasks[0])
            if self._scanner_ai_channels:
                daq.DAQmxStopTask(self._scanner_analog_daq_task)
            if self._odmr_pulser_daq_task:
                daq.DAQmxStopTask(self._odmr_pulser_daq_task)

            # prepare array to return data
            all_data = np.full((len(self.get_odmr_channels()), length),
                               222,
                               dtype=np.float64)
            start_index = 0
            if self._scanner_counter_channels:
                # create a new array for the final data (this time of the length
                # number of samples)
                real_data = np.zeros((self._odmr_length, ), dtype=np.uint32)

                # add upp adjoint pixels to also get the counts from the low time of
                # the clock:

                real_data += odmr_data[1:-1:2]
                real_data += odmr_data[:-1:2]

                if self._odmr_pulser_daq_task:
                    differential_data = np.zeros((self.oversampling * length, ), dtype=np.float64)

                    differential_data += real_data[1::2]
                    differential_data -= real_data[::2]
                    differential_data = np.divide(differential_data, real_data[::2],
                                                  np.zeros_like(differential_data),
                                                  where=real_data[::2] != 0)

                    all_data[0] = np.median(np.reshape(differential_data,
                                                       (-1, self.oversampling)),
                                            axis=1
                                            )
                else:
                    all_data[0] = np.array(real_data * self._scanner_clock_frequency, np.float64)
                start_index += 1

            if self._scanner_ai_channels:
                if self._odmr_pulser_daq_task:
                    for i, analog_data in enumerate(odmr_analog_data):
                        differential_data = np.zeros((self.oversampling * length, ), dtype=np.float64)

                        differential_data += analog_data[1:-1:2]
                        differential_data -= analog_data[:-1:2]
                        differential_data = np.divide(differential_data, analog_data[:-1:2],
                                                      np.zeros_like(differential_data),
                                                      where=analog_data[:-1:2] != 0)

                        all_data[i+start_index] = np.median(np.reshape(differential_data,
                                                                       (-1, self.oversampling)),
                                                            axis=1
                                                            )

                else:
                    all_data[start_index:] = odmr_analog_data[:, :-1]

            return False, all_data
        except:
            self.log.exception('Error while counting for ODMR.')
            return True, np.full((len(self.get_odmr_channels()), 1), [-1.])

    def close_odmr(self):
        """ Closes the odmr and cleans up afterwards.

        @return int: error code (0:OK, -1:error)
        """
        retval = 0
        try:
            # disconnect the trigger channel
            daq.DAQmxDisconnectTerms(
                self._scanner_clock_channel + 'InternalOutput',
                self._odmr_trigger_channel)

        except:
            self.log.exception('Error while disconnecting ODMR clock channel.')
            retval = -1

        if self._scanner_ai_channels:
            try:
                # stop the counter task
                daq.DAQmxStopTask(self._scanner_analog_daq_task)
                # after stopping delete all the configuration of the counter
                daq.DAQmxClearTask(self._scanner_analog_daq_task)
                # set the task handle to None as a safety
                self._scanner_analog_daq_task = None
            except:
                self.log.exception('Could not close analog.')
                retval = -1

        if self._odmr_pulser_daq_task:
            try:
                # stop the pulser task
                daq.DAQmxStopTask(self._odmr_pulser_daq_task)
                # after stopping delete all the configuration of the pulser
                daq.DAQmxClearTask(self._odmr_pulser_daq_task)
                # set the task handle to None as a safety
                self._odmr_pulser_daq_task = None
            except:
                self.log.exception('Could not close pulser.')
                retval = -1

        retval = -1 if self.close_counter(scanner=True) < 0 or retval < 0 else 0
        return retval

    def get_odmr_channels(self):
        ch = list()
        if self._scanner_counter_channels:
            ch.append(self._scanner_counter_channels[0])
        ch.extend(self._scanner_ai_channels)
        return ch

    def close_odmr_clock(self):
        """ Closes the odmr and cleans up afterwards.

        @return int: error code (0:OK, -1:error)
        """
        return self.close_clock(scanner=True)

    # ================== End ODMRCounterInterface Commands ====================

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


