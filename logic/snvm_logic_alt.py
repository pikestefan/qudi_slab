from qtpy import QtCore
from collections import OrderedDict
from copy import copy
import time
import datetime
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from logic.generic_logic import GenericLogic
from core.util.mutex import Mutex
from core.connector import Connector
from core.statusvariable import StatusVar

"""
This first part is a copy-paste of the original code in qudi.
"""

class OldConfigFileError(Exception):
    """ Exception that is thrown when an old config file is loaded.
    """

    def __init__(self):
        super().__init__('Old configuration file detected. Ignoring confocal history.')


class ConfocalHistoryEntry(QtCore.QObject):
    """ This class contains all relevant parameters of a Confocal scan.
        It provides methods to extract, restore and serialize this data.
    """

    def __init__(self, snvm):
        """ Make a confocal data setting with default values. """
        super().__init__()

        #####
        # Setting up the scanning initial parameters, and get the two stack names.
        #####
        self.sampleStackName, self.tipStackName = snvm._scanning_device.get_stack_names()

        self._active_stack = self.sampleStackName  # Default stack is sample

        # Get the maximum scanning ranges, and the position to voltage conversion factors, and put them in a dictionary
        self.x_maxrange = dict()
        self.y_maxrange = dict()
        self.position_to_voltage_arrays = dict()

        for stack in [self.sampleStackName, self.tipStackName]:
            x_range, y_range = snvm._scanning_device.get_position_range(stack=stack)
            x_volt_range, y_volt_range = snvm._scanning_device.get_voltage_range(stack=stack)
            self.x_maxrange[stack] = x_range
            self.y_maxrange[stack] = y_range

            self.position_to_voltage_arrays[stack] = [(x_volt_range[1] - x_volt_range[0]) / (x_range[1] - x_range[0]),
                                                      (y_volt_range[1] - y_volt_range[0]) / (y_range[1] - y_range[0])]

        # These are the scanning ranges that will be used for the scanning
        self.scanning_x_range = [0,1e-6]
        self.scanning_y_range = [0,1e-6]

        self.scanning_x_resolution = 10
        self.scanning_y_resolution = 10

        # Integration time per pixel
        self.px_time = 30e-3
        self._photon_samples = 0
        self.backward_pixels = 0  # The number is determined by the clock frequency and the bw_speed

        self.store_retrace = False

        self.invalid = np.nan  # Number corresponding to invalid data points.

        ####
        # Set up the ODMR scanner parameters
        ####
        self.start_freq = 2.82e9
        self.stop_freq = 2.92e9
        self.freq_resolution = 1e6
        self.mw_power = -100
        self.odmr_averages = 1

    def restore(self, snvm):
        """ Write data back into confocal logic and pull all the necessary strings """
        snvm.sampleStackName = self.sampleStackName
        snvm.tipStackName = self.tipStackName
        snvm._active_stack = self._active_stack
        snvm.x_maxrange = self.x_maxrange
        snvm.y_maxrange = self.y_maxrange
        snvm.position_to_voltage_arrays = self.position_to_voltage_arrays
        snvm.scanning_x_range = self.scanning_x_range
        snvm.scanning_y_range = self.scanning_y_range
        snvm.scanning_x_resolution = self.scanning_x_resolution
        snvm.scanning_y_resolution = self.scanning_y_resolution
        snvm.px_time = self.px_time
        snvm._photon_samples = self._photon_samples
        snvm.backward_pixels = self.backward_pixels
        snvm.store_retrace = self.store_retrace
        snvm.invalid = self.invalid
        snvm.start_freq = self.start_freq
        snvm.stop_freq = self.stop_freq
        snvm.freq_resolution = self.freq_resolution
        snvm.mw_power = self.mw_power
        snvm.odmr_averages = self.odmr_averages

    def snapshot(self, snvm):
        """ Extract all necessary data from a confocal logic and keep it for later use """
        self.sampleStackName = snvm.sampleStackName
        self.tipStackName = snvm.tipStackName
        self._active_stack = snvm._active_stack
        self.x_maxrange = snvm.x_maxrange
        self.y_maxrange = snvm.y_maxrange
        self.position_to_voltage_arrays = snvm.position_to_voltage_arrays
        self.scanning_x_range = snvm.scanning_x_range
        self.scanning_y_range = snvm.scanning_y_range
        self.scanning_x_resolution = snvm.scanning_x_resolution
        self.scanning_y_resolution = snvm.scanning_y_resolution
        self.px_time = snvm.px_time
        self._photon_samples = snvm._photon_samples
        self.backward_pixels = snvm.backward_pixels
        self.store_retrace = snvm.store_retrace
        self.invalid = snvm.invalid
        self.start_freq = snvm.start_freq
        self.stop_freq = snvm.stop_freq
        self.freq_resolution = snvm.freq_resolution
        self.mw_power = snvm.mw_power
        self.odmr_averages = snvm.odmr_averages

    def serialize(self):
        """ Give out a dictionary that can be saved via the usual means """
        serialized = dict()
        serialized['sampleStackName'] = self.sampleStackName
        serialized['tipStackName'] = self.tipStackName
        serialized['active_stack'] = self._active_stack
        serialized['x_maxrange'] = self.x_maxrange
        serialized['y_maxrange'] = self.y_maxrange
        serialized['position_to_voltage_arrays'] = self.position_to_voltage_arrays
        serialized['scanning_x_range'] = self.scanning_x_range
        serialized['scanning_y_range'] = self.scanning_y_range
        serialized['scanning_x_resolution'] = self.scanning_x_resolution
        serialized['scanning_y_resolution'] = self.scanning_y_resolution
        serialized['px_time'] = self.px_time
        serialized['_photon_samples'] = self._photon_samples
        serialized['backward_pixels'] = self.backward_pixels
        serialized['store_retrace'] = self.store_retrace
        serialized['invalid'] = self.invalid
        serialized['start_freq'] = self.start_freq
        serialized['stop_freq'] = self.stop_freq
        serialized['freq_resolution'] = self.freq_resolution
        serialized['mw_power'] = self.mw_power
        serialized['odmr_averages'] = self.odmr_averages
        return serialized

    def deserialize(self, serialized):
        """ Restore Confocal history object from a dict """
        if 'sampleStackName' in serialized:
            self.sampleStackName = serialized['sampleStackName']
        if 'tipStackName' in serialized:
            self.tipStackName = serialized['tipStackName']
        if 'active_stack' in serialized:
            self._active_stack = serialized['active_stack']
        if 'x_maxrange' in serialized:
            self.x_maxrange = serialized['x_maxrange']
        if 'y_maxrange' in serialized:
            self.y_maxrange = serialized['y_maxrange']
        if 'position_to_voltage_arrays' in serialized:
            self.position_to_voltage_arrays = serialized['position_to_voltage_arrays']
        if 'scanning_x_range' in serialized:
            self.scanning_x_range = serialized['scanning_x_range']
        if 'scanning_y_range' in serialized:
            self.scanning_y_range = serialized['scanning_y_range']
        if 'scanning_x_resolution' in serialized:
            self.scanning_x_resolution = serialized['scanning_x_resolution']
        if 'scanning_y_resolution' in serialized:
            self.scanning_y_resolution = serialized['scanning_y_resolution']
        if 'px_time' in serialized:
            self.px_time = serialized['px_time']
        if '_photon_samples' in serialized:
            self._photon_samples = serialized['_photon_samples']
        if 'backward_pixels' in serialized:
            self.backward_pixels = serialized['backward_pixels']
        if 'store_retrace' in serialized:
            self.store_retrace = serialized['store_retrace']
        if 'invalid' in serialized:
            self.invalid = serialized['invalid']
        if 'start_freq' in serialized:
            self.start_freq = serialized['start_freq']
        if 'stop_freq' in serialized:
            self.stop_freq = serialized['stop_freq']
        if 'freq_resolution' in serialized:
            self.freq_resolution = serialized['freq_resolution']
        if 'mw_power' in serialized:
            self.mw_power = serialized['mw_power']
        if 'odmr_averages' in serialized:
            self.odmr_averages = serialized['odmr_averages']


class SnvmLogic(GenericLogic):

    doublescanner = Connector(interface='SnvmScannerInterface')
    odmrscanner = Connector(interface='MicrowaveInterface')
    savelogic = Connector(interface='HDF5SaveLogic')
    optimizer_logic = Connector(interface='OptimizerLogicPxScan')

    slow_motion_clock_rate = StatusVar('slow_motion_clock_rate', 10)
    backward_speed = StatusVar('slow_motion_speed', 1)
    max_history_length = StatusVar(default=10)

    # Initialize the optimization-while-scanning settings
    optimize_while_scanning = StatusVar('optimize_while_scanning', False)
    every_N_pixels = StatusVar('every_N_pixels', 10)

    # signals
    signal_start_snvm = QtCore.Signal()
    signal_continue_snvm = QtCore.Signal()
    signal_start_confocal = QtCore.Signal()
    signal_continue_confocal = QtCore.Signal()
    signal_stop_scan = QtCore.Signal()
    signal_start_optimizer = QtCore.Signal()

    signal_snvm_image_updated = QtCore.Signal()
    signal_xy_image_updated = QtCore.Signal()
    signal_freq_px_acquired = QtCore.Signal(int) #Emits the row of the temporary data matrix to store the odmr data
    signal_xy_px_acquired = QtCore.Signal()
    signal_scan_finished = QtCore.Signal(bool) #Emits True if the scan was snvm, False otherwise

    signal_goto_start = QtCore.Signal(list, str, bool, str)
    signal_moved_to_point = QtCore.Signal(str)

    signal_snvm_initialized = QtCore.Signal()
    signal_confocal_initialized = QtCore.Signal()

    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)

        self.threadlock = Mutex()
        self.stopRequested = False

    def on_activate(self):
        self._scanning_device = self.doublescanner()
        self._odmrscanner = self.odmrscanner()
        self._savelogic = self.savelogic()
        self._optimizer_logic = self.optimizer_logic()

        self.set_slowmotion_clockrate(self.slow_motion_clock_rate)
        self.set_motion_speed(self.backward_speed)

        self.history = []
        for i in reversed(range(1, self.max_history_length)):
            try:
                new_history_item = ConfocalHistoryEntry(self)
                new_history_item.deserialize(
                    self._statusVariables['history_{0}'.format(i)])
                self.history.append(new_history_item)
            except KeyError:
                pass
            except OldConfigFileError:
                self.log.warning(
                    'Old style config file detected. History {0} ignored.'.format(i))
            except:
                self.log.warning(
                    'Restoring history {0} failed.'.format(i))
        try:
            new_state = ConfocalHistoryEntry(self)
            new_state.deserialize(self._statusVariables['history_0'])
            new_state.restore(self)
        except:
            new_state = ConfocalHistoryEntry(self)
            new_state.restore(self)
        finally:
            self.history.append(new_state)

        self.history_index = len(self.history) - 1

        #These are the indices which will be used to scan through the arrays of the frequencies and position pairs
        self._x_scanning_index = 0
        self._y_scanning_index = 0
        self._x_index_step = 1 #This coefficient is decided to decide the direction of the x scanning
        self._freq_scanning_index = 0
        self._odmr_rep_index = 0 #To keep track of the averages
        self._is_retracing = False
        self.stopRequested = False

        #Initalize the attributes that will be the scan data containers
        self._x_scanning_axis = None
        self._y_scanning_axis = None
        self.xy_scan_matrix = None
        self.snvm_matrix = None
        self._temp_afm_matrix = None #Matrix used to store the AFM values while scanning the ODMR.
        self.temp_freq_matrix = None #This matrix is used to store the ODMR traces to be averaged.
        self.xy_scan_matrix_retrace = None
        self.snvm_matrix_retrace = None
        self.freq_axis = None
        self._xy_matrix_width = None # This variable is used only to keep the code clean and reduce the calls to np.shape
        self._freq_axis_length = None #Same here
        self.average_odmr_trace = None

        self._snvm_active = False

        self._curr_snvm_image = []
        self._curr_afm_image = []
        self._curr_cfc_image = []

        #Now connect all the signals
        self.signal_continue_snvm.connect(self.continue_snvm_scanning, QtCore.Qt.QueuedConnection)
        self.signal_stop_scan.connect(self.stop_scanning, QtCore.Qt.QueuedConnection)
        self.signal_xy_px_acquired.connect(self.move_to_xy_pixel, QtCore.Qt.QueuedConnection)
        self.signal_freq_px_acquired.connect(self.move_to_freq_pixel, QtCore.Qt.QueuedConnection)
        self.signal_continue_confocal.connect(self.continue_confocal_scanning, QtCore.Qt.QueuedConnection)
        self.signal_goto_start.connect(self._go_to_point, QtCore.Qt.QueuedConnection)
        self.signal_start_optimizer.connect(self.start_optimizer, QtCore.Qt.QueuedConnection)
        self._optimizer_logic.sigRefocusFinished.connect(self._optimization_complete)

    def on_deactivate(self):
        """ Reverse steps of activation
        @return int: error code (0:OK, -1:error)
        """
        closing_state = ConfocalHistoryEntry(self)
        closing_state.snapshot(self)
        self.history.append(closing_state)
        histindex = 0
        for state in reversed(self.history):
            self._statusVariables['history_{0}'.format(histindex)] = state.serialize()
            histindex += 1
        return 0

    def _prepare_data_matrices(self):
        #Clip the ranges if they are out of bound
        self.check_xy_ranges()

        #Generate axes and matrices to store the data
        x_axis = np.linspace(self.scanning_x_range[0], self.scanning_x_range[1], self.scanning_x_resolution)
        y_axis = np.linspace(self.scanning_y_range[0], self.scanning_y_range[1], self.scanning_y_resolution)

        #Now generate the matrices to store the data
        xy_scan_matrix = np.zeros((len(y_axis), len(x_axis)), dtype=np.float64)
        # FIXME: for now the stack scanner is the one that's assumed to have the ESR sequence. Maybe consider a flexible
        #  way of doing this
        if self._snvm_active:
            step_number = 1 + round((self.stop_freq - self.start_freq) / self.freq_resolution)
            freq_axis = np.linspace(self.start_freq, self.stop_freq, step_number)
            snvm_matrix = np.zeros((xy_scan_matrix.shape[0], xy_scan_matrix.shape[1],
                                    len(freq_axis), self.odmr_averages), dtype=np.float64)
            temp_freq_matrix = np.full((self.odmr_averages, len(freq_axis)), self.invalid)
            temp_afm_matrix = np.copy(temp_freq_matrix)
            average_odmr_trace = np.zeros((temp_freq_matrix.shape[1],))
        else:
            snvm_matrix = np.copy(xy_scan_matrix[:, :, np.newaxis])
            freq_axis = None
            temp_freq_matrix = None
            average_odmr_trace = None
            temp_afm_matrix = None

        if self.store_retrace:
            xy_scan_matrix_retrace = np.copy(xy_scan_matrix)
            snvm_matrix_retrace = np.copy(snvm_matrix)
        else:
            xy_scan_matrix_retrace = np.zeros(xy_scan_matrix.shape)
            snvm_matrix_retrace = np.zeros(snvm_matrix.shape)

        self._x_scanning_axis = x_axis
        self._y_scanning_axis = y_axis
        self.xy_scan_matrix = xy_scan_matrix
        self.snvm_matrix = snvm_matrix
        self.temp_freq_matrix = temp_freq_matrix
        self._temp_afm_matrix = temp_afm_matrix
        self.average_odmr_trace = average_odmr_trace
        self.xy_scan_matrix_retrace = xy_scan_matrix_retrace
        self.snvm_matrix_retrace = snvm_matrix_retrace
        self.freq_axis = freq_axis
        self._xy_matrix_width = self.xy_scan_matrix.shape[0]
        self._freq_axis_length = len(self.freq_axis) if isinstance(self.freq_axis, np.ndarray) else None

    def check_xy_ranges(self):
        """
        Check that the requested scanning ranges are within the maximum scanning ranges set by the hardware.
        If they are out of bounds, clip the values.
        """
        curr_x_minrange, curr_x_maxrange = self.x_maxrange[self._active_stack]
        curr_y_minrange, curr_y_maxrange = self.y_maxrange[self._active_stack]
        # TODO: emit a signal when the clipping happens, and update the GUI limits accordingly
        if not ((curr_x_minrange <= self.scanning_x_range[0] <= curr_x_maxrange) or
           not (curr_x_minrange <= self.scanning_x_range[1] <= curr_x_maxrange)):
            self.scanning_x_range = np.clip(self.scanning_x_range, curr_x_minrange, curr_x_maxrange)
            self.log.warning("x scanning range limits are out of bounds, clipped back to the maximum values.")

        if not ((curr_y_minrange <= self.scanning_y_range[0] <= curr_y_maxrange) or
           not (curr_y_minrange <= self.scanning_y_range[1] <= curr_y_maxrange)):
            self.scanning_y_range = np.clip(self.scanning_y_range, curr_y_minrange, curr_y_maxrange)
            self.log.warning("y scanning range limits are out of bounds, clipped back to the maximum values.")
        return 0

    def prepare_devices(self):
        """
        Initalize the scanning device, prior to starting the scan.
        """

        self._scanning_device.module_state.lock()
        self._photon_samples= self.pxtime_to_samples()
        if self._snvm_active:
            analog_channels = self._scanning_device.get_ai_counter_channels(stack_name=self._active_stack)
        else:
            analog_channels = None
        self._scanning_device.prepare_counters(samples_to_acquire=self._photon_samples,
                                               counter_ai_channels=analog_channels)

        self._scanning_device.create_ao_task(self._active_stack)
        if self.store_retrace is False:
            self._scanning_device.prepare_motion_clock()
            clk_freq = self._scanning_device.get_motion_clock_frequency()
            speed = self._scanning_device.get_motion_speed()
            self.backward_pixels = int(((self._x_scanning_axis.max() - self._x_scanning_axis.min()) / speed) *
                                       clk_freq)
            if self.backward_pixels < 2:
                self.backward_pixels = 2

        if self._snvm_active:
            self._odmrscanner.module_state.lock()
            try:
                pass
                # FIXME: look into how to operate the srs in list mode
                #self._odmrscanner.set_list(frequency=self.freq_axis, power=self.mw_power)
            except:
                self.log.error("Failed loading the frequency axis into the ODMR scanner. Aborted execution.")
                #self._scanning_device.module_state.unlock()
                #self._odmrscanner.module_state.unlock()

    def start_snvm_scanning(self):

        self.stopRequested = False

        self.module_state.lock()

        self._active_stack = self.sampleStackName
        self._snvm_active = True
        self._prepare_data_matrices()

        self._initialize_scanning_statuses()
        self.prepare_devices()

        self.signal_snvm_image_updated.emit()
        self.signal_snvm_initialized.emit()

        #self._scanning_device.scanner_set_position([self._x_scanning_axis[self._x_scanning_index],
        #                                           self._y_scanning_axis[self._y_scanning_index]],
        #                                           stack=self._active_stack)
        """
        self._scanning_device.scanner_slow_motion([self._x_scanning_axis[self._x_scanning_index],
                                                   self._y_scanning_axis[self._y_scanning_index]],
                                                  stack=self._active_stack, clear_ao_whenfinished=False)
        """
        self.go_to_point([self._x_scanning_axis[self._x_scanning_index],
                          self._y_scanning_axis[self._y_scanning_index]],
                         stack=self._active_stack, clear_ao_whenfinished=False,
                         caller="logic")
        #FIXME: look into how to operate the srs in list mode
        self._odmrscanner.set_frequency(self.freq_axis[self._freq_scanning_index])
        self._odmrscanner.set_power(self.mw_power)
        self._odmrscanner.on()

        self.signal_continue_snvm.emit()

    def start_confocal_scanning(self):
        self.module_state.lock()

        self._active_stack = self.tipStackName
        self._snvm_active = False
        self._prepare_data_matrices()

        self._initialize_scanning_statuses()
        self.prepare_devices()

        self.signal_xy_image_updated.emit()
        self.signal_confocal_initialized.emit()

        #self._scanning_device.scanner_set_position([self._x_scanning_axis[self._x_scanning_index],
        #                                            self._y_scanning_axis[self._y_scanning_index]],
        #                                           stack=self._active_stack)
        """
        self._scanning_device.scanner_slow_motion([self._x_scanning_axis[self._x_scanning_index],
                                                   self._y_scanning_axis[self._y_scanning_index]],
                                                  stack=self._active_stack, clear_ao_whenfinished=False)
        """
        self.go_to_point([self._x_scanning_axis[self._x_scanning_index],
                          self._y_scanning_axis[self._y_scanning_index]],
                         stack=self._active_stack, clear_ao_whenfinished=False,
                         caller="logic")
        self.signal_continue_confocal.emit()

    def continue_snvm_scanning(self):
        acquire_data = False if (self.store_retrace is False) and (self._is_retracing is True) else True
        if acquire_data:
            #If the index of the ODMR is less than the averages, keep acquiring
            if self._odmr_rep_index < self.odmr_averages:
                counts, ainput = self.acquire_pixel()
                counts = counts[0] / self.px_time
                self.temp_freq_matrix[self._odmr_rep_index, self._freq_scanning_index] = counts
                self._temp_afm_matrix[self._odmr_rep_index, self._freq_scanning_index] = ainput.mean()

                if self._odmr_rep_index > 0:
                    self.average_odmr_trace = np.nanmean(self.temp_freq_matrix, axis=0)

                #TODO: in the GUI, when the user changes between trace and retrace, or changes the frequency slice,
                # the refresh plotting takes the average, along the last axis of the snvm_matrix. This means that it will
                # display an extra pixel with lower value, until the move_to_next_pixel is called.
                if self._is_retracing and self.store_retrace:
                    self.snvm_matrix_retrace[self._y_scanning_index, self._x_scanning_index,
                                             self._freq_scanning_index, self._odmr_rep_index] = counts
                else:
                    self.snvm_matrix[self._y_scanning_index, self._x_scanning_index,
                                     self._freq_scanning_index, self._odmr_rep_index] = counts
                self.signal_freq_px_acquired.emit(self._odmr_rep_index)

            #Else, the OMDR acquisition for the pixel has finished. Store the data and ask for the next pixel. If also the
            #Scanning is done, tell that the scanning has finished.
            else:
                if self._is_retracing and self.store_retrace:
                    # FIXME: I am not acquiring and storing properly the analog input, find a way after basic debugging done
                    self.xy_scan_matrix_retrace[self._y_scanning_index, self._x_scanning_index] = self._temp_afm_matrix.mean()
                else:
                    self.xy_scan_matrix[self._y_scanning_index, self._x_scanning_index] = self._temp_afm_matrix.mean()

                #The ODMR sequence has finished. Update the indices accordingly
                self._x_scanning_index += self._x_index_step
                self._odmr_rep_index = 0
                self.temp_freq_matrix[:] = self.invalid
                self.average_odmr_trace[:] = self.invalid

                self.signal_snvm_image_updated.emit()
                self.signal_xy_px_acquired.emit()

    def continue_confocal_scanning(self):
        acquire_data = False if (self.store_retrace is False) and (self._is_retracing is True) else True
        if acquire_data:
            counts, _ = self.acquire_pixel()
            counts = counts[0] / self.px_time
            if self._is_retracing and self.store_retrace:
                self.xy_scan_matrix_retrace[self._y_scanning_index, self._x_scanning_index] = counts
            else:
                self.xy_scan_matrix[self._y_scanning_index, self._x_scanning_index] = counts
        self.signal_xy_image_updated.emit()
        self._x_scanning_index += self._x_index_step
        self.signal_xy_px_acquired.emit()

    def pxtime_to_samples(self):
        return round(self.px_time * self._scanning_device.get_counter_clock_frequency())

    def acquire_pixel(self):
        data = self._scanning_device.read_pixel(self._photon_samples)
        return data

    def move_to_xy_pixel(self):
        if not self._is_retracing and self._x_scanning_index == len(self._x_scanning_axis):
            self._is_retracing = True
            self._x_index_step = -1
            self._x_scanning_index += self._x_index_step
            if (self._y_scanning_index == len(self._y_scanning_axis)-1 ) and not self.store_retrace:
                self.stopRequested = True
        elif self._is_retracing and self._x_scanning_index < 0:
            self._is_retracing = False
            self._x_index_step = 1
            self._x_scanning_index += self._x_index_step
            self._y_scanning_index += 1
            if self._y_scanning_index == len(self._y_scanning_axis):
                self.stopRequested = True

        if not self.stopRequested:
            if self._is_retracing and not self.store_retrace:
                retrace_line = np.linspace(self._x_scanning_axis.max(), self._x_scanning_axis.min(),
                                           self.backward_pixels)
                retrace_line = np.vstack((retrace_line, np.full(retrace_line.shape,
                                                                self._y_scanning_axis[self._y_scanning_index])))
                self._scanning_device.move_along_line(position_array=retrace_line, stack=self._active_stack)
                self._x_scanning_index = 0
                self._y_scanning_index += 1
                self._is_retracing = False
                self._x_index_step = 1

            new_x_pos = self._x_scanning_axis[self._x_scanning_index]
            new_y_pos = self._y_scanning_axis[self._y_scanning_index]
            self._scanning_device.scanner_set_position([new_x_pos, new_y_pos], stack=self._active_stack)

            if self._snvm_active:
                # Check if an optimization needs to be done
                total_scanning_index = self._x_scanning_index + len(self._x_scanning_axis) * self._y_scanning_index
                print(f'total_scanning_index: {total_scanning_index}')

                if self.optimize_while_scanning and total_scanning_index % self.every_N_pixels == 0:
                    self.signal_start_optimizer.emit()
                else:
                    self.signal_continue_snvm.emit()
            else:
                self.signal_continue_confocal.emit()
        else:
            self.signal_stop_scan.emit()

    def start_optimizer(self):
        print('optimizing...')
        print('pausing tasks')
        self._scanning_device.pause_tasks()

        print('starting refocus')
        self._optimizer_logic.start_refocus()

    def _optimization_complete(self, coords):
        if self.stopRequested:
            self.signal_stop_scan.emit()

        elif self._snvm_active:
            self.go_to_point(coords, stack=self._optimizer_logic.optimizer_stack)

            print(f'self._odmrscanner.module_state: {self._odmrscanner.module_state()}')
            if self._odmrscanner.module_state() == 'locked':
                print('unlocking odmrscanner')
                self._odmrscanner.module_state.unlock()

            if self._scanning_device.module_state() == 'locked':
                print('unlocking scanner')
                self._scanning_device.module_state.unlock()

            print('_optimization_complete repreparing devices')
            self.prepare_devices()
            self.signal_continue_snvm.emit()

    def move_to_freq_pixel(self):
        if not self.stopRequested:
            self._freq_scanning_index += 1
            if self._freq_scanning_index == self._freq_axis_length:
                self._odmr_rep_index += 1
                self._freq_scanning_index = 0
            self._odmrscanner.set_frequency(self.freq_axis[self._freq_scanning_index])
            self.signal_continue_snvm.emit()
        else:
            self.signal_stop_scan.emit()

    def _store_data_matrices(self, is_snvm_data=True):
        if is_snvm_data:
            self._curr_snvm_image = [self._x_scanning_axis, self._y_scanning_axis,
                                     self.freq_axis, np.arange(self.odmr_averages),
                                     self.snvm_matrix, self.snvm_matrix_retrace]
            self._curr_afm_image = [self._x_scanning_axis, self._y_scanning_axis,
                                    self.xy_scan_matrix, self.xy_scan_matrix_retrace]
        else:
            self._curr_cfc_image = [self._x_scanning_axis, self._y_scanning_axis,
                                    self.xy_scan_matrix, self.xy_scan_matrix_retrace]

    def stop_scanning(self):
        if self.stopRequested:
            with self.threadlock:
                self._scanning_device.close_counters()
                self.stopRequested = False
                self.stop_xy_scanner()
                if self._snvm_active:
                    self.stop_freq_scanner()
                if not self.store_retrace:
                    self._scanning_device.clear_motion_clock()
                self.module_state.unlock()

            self._store_data_matrices(self._snvm_active)
            self.signal_scan_finished.emit(self._snvm_active)
            self._snvm_active = False

    def stop_xy_scanner(self):
        """Closing the scanner device.

        @return int: error code (0:OK, -1:error)
        """
        try:
            self._scanning_device.close_counters()
        except Exception as e:
            self.log.exception('Could not close the scanning tasks.')
        try:
            self._scanning_device.clear_ao_task(self._active_stack)
        except:
            self.log.exception("Could not clear the ao task.")
        try:
            if self._sca.module_state() == 'locked':
                self._scanning_device.module_state.unlock()
        except Exception as e:
            self.log.exception('Could not unlock scanning device.')

        return 0

    def stop_freq_scanner(self):
        self._odmrscanner.off()
        try:
            self._odmrscanner.module_state.unlock()
        except Exception as e:
            self.log.exception('Could not unlock scanning device.')

        return 0

    def go_to_point(self, xy_coord, stack=None, clear_ao_whenfinished=True,
                    caller="logic"):
        """
        This public method emits only the signal that calls the private
        _go_to_point, which calls the scanner_slow_motion of the scanning device.
        This is done to prevent the GUI from freezing when the scanner_slow_motion
        is running.
        """
        self.signal_goto_start.emit(xy_coord, stack, clear_ao_whenfinished, caller)

    def _go_to_point(self, xy_coord, stack, clear_ao_whenfinished, caller):
        """
        Private method that calls the scanner_slow_motion. To be called only by the public
        go_to_point.
        """
        self._scanning_device.scanner_slow_motion(xy_coord, stack=stack,
                                                  clear_ao_whenfinished=clear_ao_whenfinished)

        self.signal_moved_to_point.emit(caller)

    def get_xy_image_range(self, multiplier=1):
        """
        Multiplier is an optional parameter to convert the range to the desired units
        """
        return [[self._x_scanning_axis[0]*multiplier, self._x_scanning_axis[-1]*multiplier],
                [self._y_scanning_axis[0]*multiplier, self._y_scanning_axis[-1]*multiplier]]

    def get_xy_step_size(self, multiplier = 1):
        return [(self._x_scanning_axis[1]-self._x_scanning_axis[0])*multiplier,
                (self._y_scanning_axis[1]-self._y_scanning_axis[0])*multiplier]

    def get_stack_names(self):
        return self._scanning_device.get_stack_names()

    def get_slowmotion_clockrate(self):
        return self._scanning_device.get_motion_clock_frequency()

    def set_slowmotion_clockrate(self, clockrate):
        # FIXME: dirty trick to keep the clock rate as a status variable which sets the clock rate when reloading Qudi.
        self.slow_motion_clock_rate = clockrate
        self._scanning_device.set_motion_clock_frequency(clockrate)

    def get_motion_speed(self):
        return self._scanning_device.get_motion_speed()

    def set_motion_speed(self, speed):
        #FIXME: dirty trick to keep the motion_speed as a status variable which sets the speed when reloading Qudi.
        self.backward_speed = speed
        self._scanning_device.set_motion_speed(speed * 1e-6)

    def _initialize_scanning_statuses(self):
        self._x_scanning_index = 0
        self._y_scanning_index = 0
        self._x_index_step = 1
        self._freq_scanning_index = 0
        self._odmr_rep_index = 0
        self._is_retracing = False
        self.stopRequested = False

    def save_snvm(self):
        if len(self._curr_snvm_image) > 0:
            xaxsnvm, yaxsnvm, freqax, averages, snvm_trace, snvm_retrace = self._curr_snvm_image
            xaxafm, yaxafm, afm_trace, afm_retrace = self._curr_afm_image
            data = {'snvm/x_axis': xaxsnvm,
                    'snvm/y_axis': yaxsnvm,
                    'snvm/frequency_axis': freqax,
                    'snvm/averages': averages,
                    'snvm/snvm_trace': snvm_trace,
                    'snvm/snvm_retrace': snvm_retrace,
                    'afm/x_axis': xaxafm,
                    'afm/y_axis': yaxafm,
                    'afm/afm_trace': afm_trace,
                    'afm/afm_retrace': afm_retrace}
            self._savelogic.save_hdf5_data(data)
        else:
            self.log.info("The SNVM arrays are empty")
    def save_confocal(self):
        if len(self._curr_cfc_image) > 0:
            xax, yax, cfc_trace, cfc_retrace = self._curr_cfc_image
            data = {'confocal/x_axis': xax,
                    'confocal/y_axis': yax,
                    'confocal/confocal_trace': cfc_trace,
                    'confocal/confocal_retrace': cfc_retrace}
            self._savelogic.save_hdf5_data(data)
        else:
            self.log.info("The confocal arrays are empty")










