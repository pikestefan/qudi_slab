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

class SnvmLogic(GenericLogic):

    doublescanner = Connector(interface='SnvmScannerInterface')
    odmrscanner = Connector(interface='MicrowaveInterface')

    # signals
    signal_start_snvm = QtCore.Signal()
    signal_continue_snvm = QtCore.Signal()
    signal_start_confocal = QtCore.Signal()
    signal_continue_confocal = QtCore.Signal()
    signal_stop_scan = QtCore.Signal()

    signal_snvm_image_updated = QtCore.Signal()
    signal_xy_image_updated = QtCore.Signal()
    signal_freq_px_acquired = QtCore.Signal()
    signal_xy_px_acquired = QtCore.Signal()
    signal_scan_finished = QtCore.Signal()


    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)

        self.threadlock = Mutex()
        self.stopRequested = False

    def on_activate(self):
        self._scanning_device = self.doublescanner()
        self._odmrscanner = self.odmrscanner()

        #TODO: figure out a smart way of storing all these data in another class
        #####
        # Setting up the scanning initial parameters, and get the two stack names.
        #####
        self._sampleStackName, self._tipStackName = self._scanning_device.get_stack_names()

        self._active_stack = self._sampleStackName #Default stack is sample

        #Get the maximum scanning ranges, and the position to voltage conversion factors, and put them in a dictionary
        self.x_maxrange = dict()
        self.y_maxrange = dict()
        self.position_to_voltage_arrays = dict()

        for stack in [self._sampleStackName, self._tipStackName]:
            x_range, y_range = self._scanning_device.get_position_range(stack=stack)
            x_volt_range, y_volt_range = self._scanning_device.get_voltage_range(stack=stack)
            self.x_maxrange[stack] = x_range
            self.y_maxrange[stack] = y_range

            self.position_to_voltage_arrays[stack] = [(x_volt_range[1]-x_volt_range[0])/(x_range[1]-x_range[0]),
                                                      (y_volt_range[1]-y_volt_range[0])/(y_range[1]-y_range[0])]

        #These are the scanning ranges that will be used for the scanning
        self.scanning_x_range = 0 #Initalize to zero
        self.scanning_y_range = 0

        self.scanning_x_resolution = 0
        self.scanning_y_resolution = 0

        #Integration time per pixel
        self.px_time = 0
        self._photon_samples = 0

        #These are the indices which will be used to scan through the arrays of the frequencies and position pairs
        self._x_scanning_index = 0
        self._y_scanning_index = 0
        self._x_index_step = 1 #This coefficient is decided to decide the direction of the x scanning
        self._freq_scanning_index = 0
        self._odmr_rep_index = 0 #To keep track of the averages
        self._is_retracing = False
        self.stopRequested = False

        self.store_retrace = False

        ####
        # Set up the ODMR scanner parameters
        ####
        self.start_freq = 0
        self.stop_freq = 0
        self.freq_resolution = 0
        self.mw_power = -100
        self.odmr_averages = 0

        #Initalize the attributes that will be the scan data containers
        self._x_scanning_axis = None
        self._y_scanning_axis = None
        self.xy_scan_matrix = None
        self.snvm_matrix = None
        self._temp_afm_matrix = None #Matrix used to store the AFM values while scanning the ODMR.
        self._temp_freq_matrix = None #This matrix is used to store the ODMR traces to be averaged.
        self.xy_scan_matrix_retrace = None
        self.snvm_matrix_retrace = None
        self.freq_axis = None
        self._xy_matrix_width = None # This variable is used only to keep the code clean and reduce the calls to np.shape
        self._freq_axis_length = None #Same here

        self._snvm_active = False

        #Now connect all the signals
        self.signal_continue_snvm.connect(self.continue_snvm_scanning, QtCore.Qt.QueuedConnection)
        self.signal_stop_scan.connect(self.stop_scanning, QtCore.Qt.QueuedConnection)
        self.signal_xy_px_acquired.connect(self.move_to_xy_pixel, QtCore.Qt.QueuedConnection)
        self.signal_freq_px_acquired.connect(self.move_to_freq_pixel, QtCore.Qt.QueuedConnection)
        self.signal_continue_confocal.connect(self.continue_confocal_scanning, QtCore.Qt.QueuedConnection)

    def on_deactivate(self):
        pass

    def prepare_data_matrices(self):

        #Clip the ranges if they are out of bound
        self.check_xy_ranges()

        #Generate axes and matrices to store the data
        x_axis = np.linspace(self.scanning_x_range[0], self.scanning_x_range[1], self.scanning_x_resolution)
        y_axis = np.linspace(self.scanning_y_range[0], self.scanning_y_range[1], self.scanning_y_resolution)

        #Now generate the matrices to store the data
        xy_scan_matrix = np.full((len(x_axis), len(y_axis)), 0.)
        # FIXME: for now the stack scanner is the one that's assumed to have the ESR sequence. Maybe consider a flexible
        #  way of doing this
        if self._snvm_active:
            step_number = 1 + round((self.stop_freq - self.start_freq) / self.freq_resolution)
            freq_axis = np.linspace(self.start_freq, self.stop_freq, step_number)

            snvm_matrix = np.repeat(xy_scan_matrix[:, :, np.newaxis], len(freq_axis), axis=-1)
            temp_freq_matrix = np.zeros((self.odmr_averages, len(freq_axis)))
            temp_afm_matrix = np.copy(temp_freq_matrix)
        else:
            snvm_matrix = np.tile(xy_scan_matrix, 1)
            freq_axis = None
            temp_freq_matrix = None

        if self.store_retrace:
            xy_scan_matrix_retrace = np.copy(xy_scan_matrix)
            snvm_matrix_retrace = np.copy(snvm_matrix)
        else:
            xy_scan_matrix_retrace = None
            snvm_matrix_retrace = None

        self._x_scanning_axis = x_axis
        self._y_scanning_axis = y_axis
        self.xy_scan_matrix = xy_scan_matrix
        self.snvm_matrix = snvm_matrix
        self._temp_freq_matrix = temp_freq_matrix
        self._temp_afm_matrix = temp_afm_matrix
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

        if self._snvm_active:
            self._odmrscanner.module_state.lock()
            try:
                pass
                #self._odmrscanner.set_list(frequency=self.freq_axis, power=self.mw_power)
            except:
                self.log.error("Failed loading the frequency axis into the ODMR scanner. Aborted execution.")
                #self._scanning_device.module_state.unlock()
                #self._odmrscanner.module_state.unlock()

    def start_snvm_scanning(self):
        self.module_state.lock()

        self._active_stack = self._sampleStackName
        self._snvm_active = True
        self.prepare_data_matrices()

        self._initialize_scanning_statuses()
        self.prepare_devices()

        self._scanning_device.scanner_set_position([self._x_scanning_axis[self._x_scanning_index],
                                                    self._y_scanning_axis[self._y_scanning_index]],
                                                   stack=self._active_stack)
        #FIXME: look into how to operate the srs in list mode
        self._odmrscanner.set_frequency(self.freq_axis[self._freq_scanning_index])
        self._odmrscanner.on()

        self.signal_continue_snvm.emit()

    def start_confocal_scanning(self):
        self.module_state.lock()

        self._active_stack = self._tipStackName
        self._snvm_active = False
        self.prepare_data_matrices()

        self._initialize_scanning_statuses()
        self.prepare_devices()

        self._scanning_device.scanner_set_position([self._x_scanning_axis[self._x_scanning_index],
                                                    self._y_scanning_axis[self._y_scanning_index]],
                                                   stack=self._active_stack)
        self.signal_continue_confocal.emit()

    def continue_snvm_scanning(self):
        acquire_data = False if (self.store_retrace is False) and (self._is_retracing is True) else True
        if acquire_data:
            #If the index of the ODMR is less than the averages, keep acquiring
            if self._odmr_rep_index < self.odmr_averages:
                counts, ainput = self.acquire_pixel()
                self._temp_freq_matrix[self._odmr_rep_index, self._freq_scanning_index] = counts
                self._temp_afm_matrix[self._odmr_rep_index, self._freq_scanning_index] = ainput.mean()
                self.signal_freq_px_acquired.emit()

            #Else, the OMDR acquisition for the pixel has finished. Store the data and ask for the next pixel. If also the
            #Scanning is done, tell that the scanning has finished.
            else:
                odmr_average_array = self._temp_freq_matrix.mean(axis=0)
                if self._is_retracing and self.store_retrace:
                    self.snvm_matrix_retrace[self._x_scanning_index, self._y_scanning_index, :] = odmr_average_array
                    # FIXME: I am not acquiring and storing properly the analog input, find a way after basic debugging done
                    self.xy_scan_matrix_retrace[self._x_scanning_index, self._y_scanning_index] = np.nan
                else:
                    self.snvm_matrix[self._x_scanning_index, self._y_scanning_index, :] = odmr_average_array
                    self.xy_scan_matrix[self._x_scanning_index, self._y_scanning_index] = self._temp_afm_matrix.mean()

                #The ODMR sequence has finished. Update the indices accordingly
                self._x_scanning_index += self._x_index_step
                self._odmr_rep_index = 0
                self.signal_snvm_image_updated.emit()
                self.signal_xy_px_acquired.emit()
        else:
            self._x_scanning_index += self._x_index_step
            self.signal_xy_px_acquired.emit()

    def continue_confocal_scanning(self):
        acquire_data = False if (self.store_retrace is False) and (self._is_retracing is True) else True
        if acquire_data:
            counts, _ = self.acquire_pixel()
            if self._is_retracing and self.store_retrace:
                self.xy_scan_matrix_retrace[self._x_scanning_index, self._y_scanning_index] = counts
            else:
                self.xy_scan_matrix[self._x_scanning_index, self._y_scanning_index] = counts
        self.self.signal_xy_image_updated.emit()
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
            new_x_pos = self._x_scanning_axis[self._x_scanning_index]
            new_y_pos = self._y_scanning_axis[self._y_scanning_index]

            self._scanning_device.scanner_set_position([new_x_pos, new_y_pos], stack=self._active_stack)
            if self._snvm_active:
                self.signal_continue_snvm.emit()
            else:
                self.signal_continue_confocal.emit()
        else:
            self.signal_stop_scan.emit()

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

    def stop_scanning(self):
        if self.stopRequested:
            with self.threadlock:
                self._scanning_device.close_counters()
                self.stopRequested = False
                self.stop_xy_scanner()
                if self._snvm_active:
                    self.stop_freq_scanner()
                self.module_state.unlock()
            self.signal_scan_finished.emit()

    def stop_xy_scanner(self):
        """Closing the scanner device.

        @return int: error code (0:OK, -1:error)
        """
        try:
            self._scanning_device.close_counters()
        except Exception as e:
            self.log.exception('Could not close the scanning tasks.')
        try:
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

    def _initialize_scanning_statuses(self):
        self._x_scanning_index = 0
        self._y_scanning_index = 0
        self._x_index_step = 1
        self._freq_scanning_index = 0
        self._odmr_rep_index = 0
        self._is_retracing = False
        self.stopRequested = False










