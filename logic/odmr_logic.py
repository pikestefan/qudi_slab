# -*- coding: utf-8 -*-

"""
This file contains the Qudi Logic module base class.

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

from qtpy import QtCore
from collections import OrderedDict
from interface.microwave_interface import MicrowaveMode
from interface.microwave_interface import TriggerEdge
import numpy as np
import time
import datetime
import matplotlib.pyplot as plt

from logic.generic_logic import GenericLogic
from core.util.mutex import Mutex
from core.connector import Connector
from core.configoption import ConfigOption
from core.statusvariable import StatusVar

class ODMRPxLogic(GenericLogic):
    photon_counter = Connector(interface='SnvmScannerInterface')
    fitlogic = Connector(interface='FitLogic')
    mw_source = Connector(interface='MicrowaveInterface')
    savelogic = Connector(interface='HDF5SaveLogic')

    # config option
    mw_scanmode = ConfigOption(
        'scanmode',
        'LIST',
        missing='warn',
        converter=lambda x: MicrowaveMode[x.upper()])

    # ranges = StatusVar('ranges', 1)
    ranges = 1
    range_to_fit = 0
    mw_frequency = StatusVar('cw_mw_frequency', 2870e6)
    mw_power = StatusVar('cw_mw_power', -30)
    # mw_starts = StatusVar('mw_start', [2820e6])
    # mw_stops = StatusVar('mw_stop', [2920e6])
    # mw_steps = StatusVar('mw_steps', [1e6])
    mw_starts = [2820e6]
    mw_stops = [2920e6]
    mw_steps = [1e6]
    averages = StatusVar('averages', 5)  # Averages
    fc = StatusVar('fits', None)
    integration_time = StatusVar('integration_time', 30e-3)  # In ms
    _oversampling = StatusVar('oversampling', default=10)

    # Internal signals
    sigContinueOdmr = QtCore.Signal()
    sigFreqPxAcquired = QtCore.Signal()
    sigOdmrTraceAcquired = QtCore.Signal()
    sigStopOdmr = QtCore.Signal()
    sigOdmrFinished = QtCore.Signal()

    # Update signals, e.g. for GUI module
    sigOdmrFitUpdated = QtCore.Signal(np.ndarray, np.ndarray, dict, str)

    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)
        self.threadlock = Mutex()

    def on_activate(self):
        self._mw_source = self.mw_source()
        self._photon_counter = self.photon_counter()
        self._savelogic = self.savelogic()
        self._fitlogic = self.fitlogic()

        self._mw_source.off()

        ####
        # Set up the ODMR scanner parameters
        ####
        self._freq_scanning_index = 0
        self._average_index = 0

        self.range_to_fit = 0
        self.fits_performed = {}

        # Initalize the attributes that will be the scan data containers
        self.count_matrix = None  # This matrix is used to store the ODMR traces to be averaged.
        self.freq_axis = None
        self._freq_axis_length = None  # Same here
        self.average_odmr_trace = None
        self.curr_odmr_trace = None

        self.invalid = np.nan  # Number corresponding to invalid data points.

        self._photon_samples = None #Number of samples to acquire in each pixel.

        self.stopRequested = False

        self.sigContinueOdmr.connect(self.continue_odmr, QtCore.Qt.QueuedConnection)
        self.sigFreqPxAcquired.connect(self._next_freq_pixel, QtCore.Qt.QueuedConnection)
        self.sigStopOdmr.connect(self.stop_odmr)

    def on_deactivate(self):
        pass

    def get_hw_constraints(self):
        """ Return the names of all ocnfigured fit functions.
        @return object: Hardware constraints object
        """
        constraints = self._mw_source.get_limits()
        return constraints

    def _prepare_count_matrix(self):
        # Here it is assumed that the GUI has already ensured that the frequency step fits an integer amount of times
        # in the desired range.
        freq_axis_list = []
        range_indices_list = []

        for i in range(self.ranges):
            step_number = 1 + round((self.mw_stops[i] - self.mw_starts[i]) / self.mw_steps[i])
            freq_axis_piece = np.linspace(self.mw_starts[i], self.mw_stops[i], step_number)
            freq_axis_list.append(freq_axis_piece)
        freq_axis = np.concatenate(freq_axis_list)
        count_matrix = np.full((self.averages, len(freq_axis)), self.invalid)

        # Populate the range_indices_list. This list will be used to access part of the data if it is split into
        # different mw frequency ranges
        for freq_axis_piece in freq_axis_list:
            range_indices_list.append(np.isin(freq_axis, freq_axis_piece))

        self._freq_scanning_index = 0
        self._average_index = 0
        self.freq_axis = freq_axis
        self.freq_axis_list = freq_axis_list
        self.range_indices_list = range_indices_list
        self.count_matrix = count_matrix
        self.average_odmr_trace = np.zeros(freq_axis.shape)
        self._freq_axis_length = len(freq_axis) if isinstance(self.freq_axis, np.ndarray) else None

    def _prepare_devices(self):
        """
        Initialize the counter and the settings of the MW device prior to starting scanning.
        """

        self._photon_counter.module_state.lock()
        self._photon_samples = self._pxtime_to_samples()

        self._photon_counter.prepare_counters(samples_to_acquire=self._photon_samples)
        self._mw_source.set_mod(False)
        self._mw_source.module_state.lock()
        try:
            pass
            # FIXME: look into how to operate the srs in list mode
            #self._odmrscanner.set_list(frequency=self.freq_axis, power=self.mw_power)
        except:
            self.log.error("Failed loading the frequency axis into the ODMR scanner. Aborted execution.")
            #self._odmrscanner.module_state.unlock()


    def start_odmr(self):
        self.module_state.lock()

        self._prepare_count_matrix()
        self._prepare_devices()

        self._mw_source.set_frequency(self.freq_axis[self._freq_scanning_index])
        self._mw_source.set_power(self.sweep_mw_power)
        self._mw_source.on()

        self.sigContinueOdmr.emit()

    def continue_odmr(self):
        if not self.stopRequested:
            counts = self._acquire_pixel()
            # FIXME: Why is counts sometimes empty?
            if len(counts) == 0:
                return

            self.count_matrix[self._average_index, self._freq_scanning_index] = counts / self.integration_time
            self.curr_odmr_trace = self.count_matrix[self._average_index]

            if self._average_index > 0:
                self.average_odmr_trace = np.nanmean(self.count_matrix, axis=0)

        self.sigFreqPxAcquired.emit()

    def _next_freq_pixel(self):
        self._freq_scanning_index += 1
        if self._freq_scanning_index == self._freq_axis_length:
            self._average_index += 1
            self._freq_scanning_index = 0
            self.sigOdmrTraceAcquired.emit()

        self._mw_source.set_frequency(self.freq_axis[self._freq_scanning_index])
        if (self._average_index == self.averages) or self.stopRequested:
            self.stopRequested = True
            self.sigStopOdmr.emit()
        else:
            self.sigContinueOdmr.emit()

    def stop_odmr(self):
        if self.stopRequested:
            with self.threadlock:
                self._mw_source.off()
                self._photon_counter.close_counters()
                try:
                    self._photon_counter.module_state.unlock()
                except Exception as e:
                    self.log.exception('Could not unlock counter device.')
                try:
                    self._mw_source.module_state.unlock()
                except Exception as e:
                    self.log.exception('Could not unlock mw source.')
                self.module_state.unlock()
            self.stopRequested = False
            self.sigOdmrFinished.emit()

    def save_data(self):
        odmr = {
            'odmr/frequency_axis': self.freq_axis,
            'odmr/averages': self.averages,
            'odmr/odmr_traces': self.count_matrix,
            'odmr/odmr_average_trace': self.average_odmr_trace,
        }

        # Check if fits have been performed. If yes, then store them
        # The fit parameters, like dip position and width will be stored as .h5 attributes to the fit data
        fit = {}
        fit_attributes = {}
        for key, fit_record in self.fits_performed.items():
            odmr_fit_x, odmr_fit_y, result, fit_type = fit_record
            odmr_fit = np.stack((odmr_fit_x, odmr_fit_y))
            name = f'fit/{key}/odmr_fit'
            fit[name] = odmr_fit

            result_dict = result.params.valuesdict()
            # Add the fit type to the saved attributes
            result_dict['fit_type'] = fit_type
            fit_attributes[name] = result_dict

        data = {**odmr, **fit}
        self._savelogic.save_hdf5_data(data, attributes=fit_attributes)

    @fc.constructor
    def sv_set_fits(self, val):
        # Setup fit container
        fc = self.fitlogic().make_fit_container('ODMR sum', '1d')
        fc.set_units(['Hz', 'c/s'])
        if isinstance(val, dict) and len(val) > 0:
            fc.load_from_dict(val)
        else:
            d1 = OrderedDict()
            d1['Lorentzian dip'] = {
                'fit_function': 'lorentzian',
                'estimator': 'dip'
            }
            d1['Two Lorentzian dips'] = {
                'fit_function': 'lorentziandouble',
                'estimator': 'dip'
            }
            d1['N14'] = {
                'fit_function': 'lorentziantriple',
                'estimator': 'N14'
            }
            d1['N15'] = {
                'fit_function': 'lorentziandouble',
                'estimator': 'N15'
            }
            d1['Two Gaussian dips'] = {
                'fit_function': 'gaussiandouble',
                'estimator': 'dip'
            }
            default_fits = OrderedDict()
            default_fits['1d'] = d1
            fc.load_from_dict(default_fits)
        return fc

    @fc.representer
    def sv_get_fits(self, val):
        """ save configured fits """
        if len(val.fit_list) > 0:
            return val.save_to_dict()
        else:
            return None

    def do_fit(self, fit_function=None, x_data=None, y_data=None, fit_range=0):
        """
        Execute the currently configured fit on the measurement data. Optionally on passed data
        """
        if (x_data is None) or (y_data is None):
            x_data = self.freq_axis_list[fit_range]
            # Pick y_data according to the chosen fit range
            range_indices = self.range_indices_list[fit_range]
            y_data = self.average_odmr_trace[range_indices]

        if fit_function is not None and isinstance(fit_function, str):
            if fit_function in self.get_fit_functions():
                self.fc.set_current_fit(fit_function)
            else:
                self.fc.set_current_fit('No Fit')
                if fit_function != 'No Fit':
                    self.log.warning('Fit function "{0}" not available in ODMRLogic fit container.'
                                     ''.format(fit_function))

        self.odmr_fit_x, self.odmr_fit_y, result = self.fc.do_fit(x_data, y_data)
        key = f'range: {fit_range}'
        if fit_function != 'No Fit':
            self.fits_performed[key] = (self.odmr_fit_x, self.odmr_fit_y, result, self.fc.current_fit)
        else:
            if key in self.fits_performed:
                self.fits_performed.pop(key)

        if result is None:
            result_str_dict = {}
        else:
            result_str_dict = result.result_str_dict
        self.sigOdmrFitUpdated.emit(
            self.odmr_fit_x, self.odmr_fit_y, result_str_dict, self.fc.current_fit)
        return

    def get_fit_functions(self):
        """ Return the hardware constraints/limits
        @return list(str): list of fit function names
        """
        return list(self.fc.fit_list)

    def _acquire_pixel(self):
        counts, _ = self._photon_counter.read_pixel(self._photon_samples)
        return counts

    def _pxtime_to_samples(self):
        return round(self.integration_time * self._photon_counter.get_counter_clock_frequency())





