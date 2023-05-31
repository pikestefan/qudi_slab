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


class NationalInstrumentsXSeriesPxScan(Base, SnvmScannerInterface):
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
        tip_scanner_ai_channels:
            - /Dev1/AI2
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
    _scanning_photon_sources = ConfigOption('scanning_photon_sources', list(), missing='warn')
    _pulsing_photon_sources = ConfigOption('pulsing_photon_sources', list(), missing='warn')

    #Clock used when moving the scanners from A to B, without acquiring anything along the way
    _motion_clock_channel = ConfigOption('motion_clock_channel', missing='error')
    _default_motion_clock_frequency = ConfigOption('default_motion_clock_frequency', 100, missing='info')
    _motion_speed_conf = ConfigOption('motion_speed_conf', 1e-6, missing='info')  # in m/s
    _motion_speed_snvm = ConfigOption('motion_speed_snvm', 1e-6, missing='info')  # in m/s

    # Photon counting settings
    _clock_counter = ConfigOption('counter_clock', '100kHzTimebase', missing='info')
    _counter_clock_frequency = ConfigOption('counter_clock_frequency', 100e3, missing='info')

    _scanning_counter_channels = ConfigOption('scanning_counter_channels', missing='error')
    _pulsing_counter_channels = ConfigOption('pulsing_counter_channels', missing='error')
    _counter_ai_channels = ConfigOption('counter_ai_channels', list(), missing='info')
    _counter_voltage_range = ConfigOption('counter_voltage_range', [-10, 10], missing='info')

    # Photo diode setting
    _photo_diode_channel = ConfigOption('photo_diode_channel', missing='error')
    _photo_diode_voltage_range = ConfigOption('photo_diode_voltage_range', [0, 10], missing='info')

    # Sample scanner
    _sample_scanner_ao_channels = ConfigOption('sample_scanner_ao_channels', missing='error')
    _sample_scanner_ai_channels = ConfigOption('sample_scanner_ai_channels', list(), missing='info')
    _sample_scanner_voltage_ranges = ConfigOption('sample_scanner_voltage_ranges', missing='error')
    _sample_scanner_position_ranges = ConfigOption('sample_scanner_position_ranges', missing='error')

    # Tip scanner
    _tip_scanner_ao_channels = ConfigOption('tip_scanner_ao_channels', missing='error')
    _tip_scanner_ai_channels = ConfigOption('tip_scanner_ai_channels', list(), missing='info')
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
        self._scanner_voltage_ranges = dict(zip(self._stack_names,
                                                [self._sample_scanner_voltage_ranges,
                                                 self._tip_scanner_voltage_ranges]))

        self._scanner_position_ranges = dict(zip(self._stack_names,
                                                 [self._sample_scanner_position_ranges,
                                                  self._tip_scanner_position_ranges]))

        # the tasks used on that hardware device:
        self._counter_daq_tasks = list()
        self._counter_ai_daq_task = list()
        self._photo_diode_ai_daq_task = None

        self._scanner_ao_tasks = dict().fromkeys(self._stack_names, None)
        self._motion_clock_task = None
        self._motion_clock_frequency = self._default_motion_clock_frequency

        self._clk_task = None

        self._scanning_photon_sources = self._scanning_photon_sources if self._scanning_photon_sources is not None else list()
        self._pulsing_photon_sources = self._pulsing_photon_sources if self._pulsing_photon_sources is not None else list()

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

        self._current_position = dict(zip(self._stack_names, [np.zeros((2,)), np.zeros((2,))]))

    def on_deactivate(self):
        """ Shut down the NI card.
        """
        try:
            for stack_name in self._scanner_ao_tasks:
                self.clear_ao_task(stack_name)
        except:
            self.log.exception('Could not clear AO Out Task.')

        self.reset_hardware()

    def prepare_clock(self,
                      samples_to_acquire=None,
                      clock_channel=None,
                      clock_frequency=None,
                      ):
        """
        @param samples_to_acquire: the number of clock cycles
        @param clock_channel: the counter channel to be used as clock output
        @param clock_frequency: the frequency of the clock
        @return:
        """
        if self._clk_task is not None:
            self.log.error("A clock task is already running, the requested task cannot be started. Operation aborted.")
            return -1
        if samples_to_acquire is None:
            self.log.error("Th")

        if clock_frequency is None:
            clock_frequency = self._counter_clock_frequency

        if clock_channel is None:
            clock_channel = self._clock_counter

        #try:
        taskname = 'PhotonCountingClock'
        clock_task = daq.TaskHandle()
        daq.DAQmxCreateTask(taskname, daq.byref(clock_task))
        daq.DAQmxCreateCOPulseChanFreq(clock_task,  # taskHandle
                                       clock_channel,  # counter
                                       taskname, # nameToAssignToChannel
                                       daq.DAQmx_Val_Hz, # units
                                       daq.DAQmx_Val_Low, # idleState
                                       0., # initialDelay
                                       clock_frequency, # frequency
                                       0.5 # dutyCycle
                                       )

        daq.DAQmxCfgImplicitTiming(clock_task, #taskHandle
                                   daq.DAQmx_Val_FiniteSamps, #sampleMode
                                   samples_to_acquire #sampsPerChanToAcquire
                                   )
        self._clk_task = clock_task
        # except:
        #     self.log.error("The clock initializaton has failed.")
        #     return -1

    def prepare_counters(self,
                         counter_channels=None,
                         counter_ai_channels=None,
                         sources=None,
                         counter_clock=None,
                         samples_to_acquire=None,
                         mode='scanning',
                         auto_clock=True):
        """ Configures the actual counter with a given clock.

        @param list(str) counter_channels: optional, physical channel of the counter
        @param list(str) counter_ai_channels: optional, the analog inputs that should be acquired along
                                              with the counting.
        @param list(str) sources: optional, physical channel where the photons
                                  are to count from
        @param str counter_clock: optional, specifies the clock channel for the
                                  counter
        @param int samples_to_acquire: specifies the number of samples that should be acquired for finite acquisiton.
        @param str mode: specifies the mode of the counter preparation.
        @param bool auto_clock: if True, prepare_counters also prepares the clock counter using the defaults.

        @return int: error code (0:OK, -1:error)
        """

        allowed_modes = ['pulsing', 'scanning']

        if mode not in allowed_modes:
            self.log.error("The requested counter mode is not currently specified.")

        if counter_channels:
            my_counter_channels = counter_channels
        elif mode == 'scanning':
            my_counter_channels = self._scanning_counter_channels
        elif mode == 'pulsing':
            my_counter_channels = self._pulsing_counter_channels

        if sources:
            my_photon_sources = sources
        elif mode == 'scanning':
            my_photon_sources = self._scanning_photon_sources
        elif mode == 'pulsing':
            my_photon_sources = self._pulsing_photon_sources
        my_clock_channel = counter_clock if counter_clock else self._clock_counter

        #Here prepare the clock, if in auto_clock mode
        if auto_clock:
            self.prepare_clock(samples_to_acquire=samples_to_acquire,
                               clock_channel=my_clock_channel,
                               clock_frequency=self._counter_clock_frequency
                               )

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
                                          my_clock_channel + 'InternalOutput',
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
                        my_clock_channel + 'InternalOutput',
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

    def pause_tasks(self):
        try:
            self.close_counters()
        except Exception as e:
            self.log.exception('Could not close the scanning tasks.')

        try:
            self.clear_ao_task('sample')
        except:
            self.log.exception("Could not clear the ao task.")

        try:
            self.module_state.unlock()
        except Exception as e:
            self.log.exception('Could not unlock scanning device.')

        if self._motion_clock_task is not None:
            self.clear_motion_clock()

    def resume_tasks(self):
        for task in self._counter_daq_tasks:
            daq.DAQmxStartTask(task)

        # for task in self._counter_ai_daq_task:
        #     daq.DAQmxStartTask(task)

        for key, task in self._scanner_ao_tasks.items():
            if task is not None:
                daq.DAQmxStartTask(task)

        if self._counter_ai_daq_task is not None:
            daq.DAQmxStartTask(self._counter_ai_daq_task)

        if self._photo_diode_ai_daq_task is not None:
            daq.DAQmxStartTask(self._photo_diode_ai_daq_task)

        if self._motion_clock_task is not None:
            daq.DAQmxStartTask(self._motion_clock_task)



    def prepare_photo_diode(self):
        task = daq.TaskHandle()
        daq.DAQmxCreateTask('PhotoDiode', daq.byref(task))

        daq.DAQmxCreateAIVoltageChan(
            task,
            self._photo_diode_channel,
            'Photo Diode Analog In',
            daq.DAQmx_Val_RSE,
            self._photo_diode_voltage_range[0],
            self._photo_diode_voltage_range[1],
            daq.DAQmx_Val_Volts,
            ''
        )

        self._photo_diode_ai_daq_task = task

    def close_photo_diode(self):

        try:
            daq.DAQmxStopTask(self._photo_diode_ai_daq_task)
            daq.DAQmxClearTask(self._photo_diode_ai_daq_task)
        except:
            self.log.exception('Could not close counter.')

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

        if self._clk_task is not None:
            try:
                daq.DAQmxStopTask(self._clk_task)
                daq.DAQmxClearTask(self._clk_task)
                self._clk_task = None
            except:
                self.log.exception("Could not close the clock task.")
        return error

    def prepare_motion_clock(self, clock_channel=None):

        if self._motion_clock_task is not None:
            self.log.error("The motion clock is running, close the other task first.")
            return -1

        clock_task = daq.TaskHandle()

        if clock_channel is not None:
            self._motion_clock_channel = clock_channel

        try:
            task_name = 'MotionClockTask'
            daq.DAQmxCreateTask(task_name, daq.byref(clock_task))

            daq.DAQmxCreateCOPulseChanFreq(
                clock_task,
                self._motion_clock_channel,
                'Motion clock',
                daq.DAQmx_Val_Hz,
                daq.DAQmx_Val_Low,
                0,
                self._motion_clock_frequency,
                0.5
            )

            self._motion_clock_task = clock_task
        except:
            self.log.exception("Error while setting up the motion clock.")
            return -1
        return 0

    def clear_motion_clock(self):
        try:
            daq.DAQmxStopTask(self._motion_clock_task)

            daq.DAQmxClearTask(self._motion_clock_task)

            self._motion_clock_task = None
        except:
            self.log.exception('Could not close clock.')
            return -1
        return 0

    def set_up_linemotion(self, points=1000, stack = None):
        if self._motion_clock_task is None:
            self.log.error("The motion clock has not been set up yet.")
            return -1

        if stack == None:
            self.log.error("You need to specify a stack to set up the line motion.")

        try:

            daq.DAQmxCfgSampClkTiming(
                self._scanner_ao_tasks[stack],
                self._motion_clock_channel + 'InternalOutput',
                self._motion_clock_frequency,
                daq.DAQmx_Val_Rising,
                daq.DAQmx_Val_FiniteSamps,
                points)

            daq.DAQmxCfgImplicitTiming(
                self._motion_clock_task,
                daq.DAQmx_Val_FiniteSamps,
                points + 1)
        except:
            self.log.exception("Error when setting up the line scanner.")

    def move_along_line(self, position_array=None, stack=None):
        """
        Position array needs to have shape (2, n)
        """

        motion_task = self._scanner_ao_tasks[stack]

        daq.DAQmxSetSampTimingType(motion_task, daq.DAQmx_Val_SampClk)
        self.set_up_linemotion(points=position_array.shape[1], stack=stack)

        voltages = self._scanner_position_to_volt(position_array, stack)
        self.write_voltage(motion_task, voltages=voltages, autostart=False)

        daq.DAQmxStartTask(motion_task)
        daq.DAQmxStartTask(self._motion_clock_task)

        daq.DAQmxWaitUntilTaskDone(self._motion_clock_task,
                                   2 * self._RWTimeout * position_array.shape[1])

        daq.DAQmxWaitUntilTaskDone(motion_task,
                                   2 * self._RWTimeout * position_array.shape[1])

        try:
            daq.DAQmxStopTask(self._motion_clock_task)
        except:
            self.log.exception("stopping clock didnt' work.")
        try:
            daq.DAQmxStopTask(motion_task)
        except:
            self.log.exception("Stopping ao didn't work.")

        try:
            daq.DAQmxSetSampTimingType(motion_task, daq.DAQmx_Val_OnDemand)
            daq.DAQmxStopTask(motion_task)
        except:
            self.log.exception("Changing ao mode didn't work.")

        self._current_position[stack] = position_array[:, -1]

    def get_counter_channels(self):
        """ Returns the list of counter channel names.

        @return tuple(str): channel names

        Most methods calling this might just care about the number of channels, though.
        """
        ch = self._counter_channels[:]
        return ch

    def get_ai_counter_channels(self, stack_name=None):
        if stack_name is None:
            return None
        else:
            return self._scanner_ai_channels[stack_name]

    # ================ Scanning Commands =======================
    def create_ao_task(self, stack=None):
        """
        Implemented this function to start and kill the ao every time it's required. Used in combination with
        clear_ao_task, it helps avoiding resource allocation problems in the DAQ card.
        """
        ao_task = self._scanner_ao_tasks[stack]
        try:
            # If an analog task is already running, kill that one first

            if ao_task is not None:
                # stop the analog output task
                daq.DAQmxStopTask(ao_task)

                # delete the configuration of the analog output
                daq.DAQmxClearTask(ao_task)

                # set the task handle to None as a safety
                self._scanner_ao_tasks[stack] = None

            # initialize ao channels / task for scanner, should always be active.
            # Define at first the type of the variable as a Task:
            ao_task = daq.TaskHandle()

            # create the actual analog output task on the hardware device. Via
            # byref you pass the pointer of the object to the TaskCreation function:

            daq.DAQmxCreateTask('{:s}ScannerAO'.format(stack.capitalize()),
                                daq.byref(ao_task))
            for n, chan in enumerate(self._scanner_ao_channels[stack]):
                # Assign and configure the created task to an analog output voltage channel.
                daq.DAQmxCreateAOVoltageChan(
                    # The AO voltage operation function is assigned to this task.
                    ao_task,
                    # use (all) scanner ao_channels for the output
                    chan,
                    # assign a name for that channel
                    'Sample Scanner AO Channel {0}'.format(n),
                    # minimum possible voltage
                    self._scanner_voltage_ranges[stack][n][0],
                    # maximum possible voltage
                    self._scanner_voltage_ranges[stack][n][1],
                    # units is Volt
                    daq.DAQmx_Val_Volts,
                    # empty for future use
                    '')

            self._scanner_ao_tasks[stack] = ao_task

            daq.DAQmxSetSampTimingType(ao_task, daq.DAQmx_Val_OnDemand)

        except:
            self.log.exception('Error starting analog output task.')
            return -1
        return 0

    def clear_ao_task(self, stack):
        ao_task = self._scanner_ao_tasks[stack]
        if ao_task is None:
            self.log.info(f"The AO task for the {stack} stack is already cleared.")
            return 0

        try:
            daq.DAQmxStopTask(ao_task)
            daq.DAQmxClearTask(ao_task)

            self._scanner_ao_tasks[stack] = None
        except:
            self.log.exception(f"Error clearing the AO {stack} stack task.")
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

    def get_motion_clock_frequency(self):
        return self._motion_clock_frequency

    def get_motion_speed(self, stack):
        if stack == 'sample':
            return self._motion_speed_snvm
        elif stack == 'tip':
            return self._motion_speed_conf
        else:
            raise ValueError(f'Trying to get speed for stack {stack}')

    def get_scanner_position(self, stack):
        """ Get the current position of the scanner hardware.
        """
        #TODO: one day, this is converted into an actual voltage reading
        return self._current_position[stack]

    def get_position_range(self, stack=None):
        """ Returns the physical range of the scanner.

        @return float [4][2]: array of 4 ranges with an array containing lower
                              and upper limit. The unit of the scan range is
                              meters.
        """
        return list(self._scanner_position_ranges[stack])

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
        ch = self._counter_channels[:]
        return ch

    def set_motion_speed(self, speed, stack):
        if stack == 'sample':
            self._motion_speed_snvm = speed
        elif stack == 'tip':
            self._motion_speed_conf = speed
        else:
            raise ValueError(f'Trying to set speed for stack {stack}')

    def set_motion_clock_frequency(self, frequency):
        self._motion_clock_frequency = frequency

    def read_voltage(self, samples=1):
        daq.DAQmxStartTask(self._photo_diode_ai_daq_task)

        if self._photo_diode_ai_daq_task is not None:
            try:
                analog_data = np.full(samples, 111, dtype=np.float64)
                analog_read_samples = daq.int32()

                daq.DAQmxReadAnalogF64(
                    self._photo_diode_ai_daq_task,
                    samples,
                    self._RWTimeout,
                    daq.DAQmx_Val_GroupByChannel,
                    analog_data,
                    samples,
                    daq.byref(analog_read_samples),
                    None
                )

                daq.DAQmxStopTask(self._photo_diode_ai_daq_task)
            except:
                self.log.exception("Failed reading photo diode.")

        else:
            analog_data = None

        return analog_data

    def read_pixel(self, samples=1):
        if self._counter_ai_daq_task:
            daq.DAQmxStartTask(self._counter_ai_daq_task)

        count_matrix = np.full((len(self._counter_daq_tasks),
                                samples), 0, dtype=np.uint32)

        final_counts = np.full(len(self._counter_daq_tasks), 0, dtype=np.uint32)

        read_samples = daq.int32()

        for i, task in enumerate(self._counter_daq_tasks):
            daq.DAQmxStartTask(task)

        daq.DAQmxStartTask(self._clk_task)

        for i, task in enumerate(self._counter_daq_tasks):
            # wait for the scanner counter to finish
            daq.DAQmxWaitUntilTaskDone(
                # define task
                task,
                # Maximum timeout for the counter times the positions. Unit is seconds.
                self._RWTimeout * 2 * samples)

        daq.DAQmxWaitUntilTaskDone(
            self._clk_task,
            # Maximum timeout for the counter times the positions. Unit is seconds.
            self._RWTimeout * 2 * samples)

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

        daq.DAQmxStopTask(self._clk_task)

        return final_counts, analog_data

    def write_voltage(self, ao_task, voltages=[], autostart=True):
        AONwritten = daq.int32()
        if ao_task is None:
            self.log.exception('Please specify or create an AO task')
        try:
            daq.DAQmxWriteAnalogF64(ao_task,
                                    voltages.shape[1],
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
            self._current_position[stack][0] = np.float(x)

        if y is not None:
            if not (scanner_position_range[1][0] <= y <= scanner_position_range[1][1]):
                self.log.error('You want to set y out of range: {0:f}.'.format(y))
                return -1
            self._current_position[stack][1] = np.float(y)

        my_position = np.array(xypair)
        # then directly write the position to the hardware
        voltages = self._scanner_position_to_volt(my_position[:, np.newaxis], stack=stack)
        try:
            self.write_voltage(
                scanning_task,
                voltages=voltages,
                autostart=True)

            self._current_position[stack] = my_position
        except:
            return -1
        return 0

    def set_voltage(self, voltage_x, voltage_y):
        self.create_ao_task(self._tip_stack_name)
        scanning_task = self._scanner_ao_tasks[self._tip_stack_name]
        voltages = np.zeros((2, 1), dtype='float64')
        voltages[0] = voltage_x
        voltages[1] = voltage_y
        self.write_voltage(
            scanning_task,
            voltages=voltages,
            autostart=True)
        self.clear_ao_task(self._tip_stack_name)

    def scanner_slow_motion(self, end_xy, speed=None, stack=None, clear_ao_whenfinished = True):
        """
        Used to move from A to B slowly only once, at the defined speed
        """
        if self._motion_clock_task is None:  # If a motion clock already exists, use that to move to the position.
            self.prepare_motion_clock()
            close_at_end = True
        else:
            close_at_end = False

        #Create the ao task for the scanning
        self.create_ao_task(stack)

        #Move first along y coordinate, then x
        start_xy = self.get_scanner_position(stack=stack)
        x_start, x_end = start_xy[0], end_xy[0]
        y_start, y_end = start_xy[1], end_xy[1]

        if speed is None:
            speed = self.get_motion_speed(stack)

        points_ymotion = int(self._motion_clock_frequency * abs(y_end - y_start) / speed)
        points_xmotion = int(self._motion_clock_frequency * abs(x_end - x_start) / speed)

        if points_ymotion < 2:
            points_ymotion = 2
        if points_xmotion < 2:
            points_xmotion = 2

        yline = np.zeros((2, points_ymotion))
        yline[0], yline[1] = x_start, np.linspace(y_start, y_end, points_ymotion)

        xline = np.zeros((2, points_xmotion))
        xline[0], xline[1] = np.linspace(x_start, x_end, points_xmotion), y_end

        final_line = np.hstack((yline, xline))

        self.move_along_line(final_line, stack=stack)

        if close_at_end:
            self.clear_motion_clock()

        if clear_ao_whenfinished:
            self.clear_ao_task(stack)

        self._current_position[stack] = np.array([x_end, y_end])

    def close_scanner(self, stack=None):
        """ Closes the scanners and cleans up afterwards.

        @return int: error code (0:OK, -1:error)
        """
        a = self.clear_ao_task(stack)

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

        c = self.close_counters()
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
        chanlist = []

        ao_chans = [chan for sublist in self._scanner_ao_channels.values() for chan in sublist]
        ai_chans = [chan for sublist in self._scanner_ai_channels.values() for chan in sublist]

        chanlist.extend(ao_chans)
        chanlist.extend(ai_chans)
        chanlist.extend(self._scanning_photon_sources)
        chanlist.extend(self._pulsing_photon_sources)
        chanlist.extend([self._motion_clock_channel])

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

        @param np.ndarray positions: array with (2, n) shape, first column corresponding to x and second one to y.
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
        offset = scanner_voltage_ranges[:, 0]
        offset = offset[:, np.newaxis]
        volts = post2volt_coeff * positions + offset
        if (np.any(np.logical_or(volts[0, :] < scanner_voltage_ranges[0, 0],
                                 volts[0, :] > scanner_voltage_ranges[0, 1]))
                or
            np.any(np.logical_or(volts[1, :] < scanner_voltage_ranges[1, 0],
                                 volts[1, :] > scanner_voltage_ranges[1, 1]))):
            self.log.error('Some of the requested positions are out of the voltage bounds.')
            return np.array([np.nan])
        return volts.astype(np.float64)

    def _scanner_volt_to_position(self, volts=None, stack=None):
        """ Converts a set of voltage pixels to positionss.

        @param np.ndarray positions: array with (2, n) shape, first column corresponding to x and second one to y.
        @param string stack: the stack name

        @return float[][n]: array of n-part tuples of corresponing positions
        """
        if not isinstance(volts, np.ndarray):
            self.log.error('Given positions are not and nd array.')
            return np.array([np.NaN])

        if stack is None:
            self.log.error('Please specify a stack to work with.')

        scanner_voltage_ranges = np.array(self._scanner_voltage_ranges[stack])
        scanner_position_ranges = np.array(self._scanner_position_ranges[stack])

        post2volt_coeff = np.diff(scanner_voltage_ranges, axis=1) / np.diff(scanner_position_ranges, axis=1)
        offset = scanner_voltage_ranges[:, 0]
        offset = offset[:, np.newaxis]
        positions = (volts - offset) / post2volt_coeff
        if (np.any(np.logical_or(positions[0, :] < scanner_position_ranges[0, 0],
                                 positions[0, :] > scanner_position_ranges[0, 1]))
                or
            np.any(np.logical_or(positions[1, :] < scanner_position_ranges[1, 0],
                                 positions[1, :] > scanner_position_ranges[1, 1]))):
            self.log.error('Some of the positions are out of the positions bounds.')
            return np.array([np.nan])
        return positions.astype(np.float64)
