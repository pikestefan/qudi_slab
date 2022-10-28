# -*- coding: utf-8 -*
"""
This file contains the Qudi logic class for optimizing scanner position.

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
import numpy as np
import time

from logic.generic_logic import GenericLogic
from core.connector import Connector
from core.statusvariable import StatusVar
from core.configoption import ConfigOption
from core.util.mutex import Mutex


class OptimizerLogic(GenericLogic):

    """This is the Logic class for optimizing scanner position on bright features. Lazily copied-pasted from the
    original class, then modified to fit the SnvmScannerInterface.
    """

    # declare connectors
    scanner = Connector(interface='SnvmScannerInterface')
    fitlogic = Connector(interface='FitLogic')

    # declare status vars
    optimizer_stack = StatusVar('optimizer_stack', 'tip')
    refocus_XY_size = StatusVar('xy_size', 0.6e-6)
    optimizer_XY_res = StatusVar('xy_resolution', 10)
    integration_time = StatusVar('integration_time', 20)

    # "private" signals to keep track of activities here in the optimizer logic
    _sigScanNextXyLine = QtCore.Signal()
    _sigOptimizationComplete = QtCore.Signal()
    _sigCompletedXyOptimizerScan = QtCore.Signal()

    # public signals
    sigImageUpdated = QtCore.Signal()
    sigRefocusStarted = QtCore.Signal(str)
    sigRefocusXySizeChanged = QtCore.Signal()
    sigRefocusFinished = QtCore.Signal(list)
    sigPositionChanged = QtCore.Signal(float, float, float)
    sigGoToPoint = QtCore.Signal(list)

    def __init__(self, config, **kwargs):
        super().__init__(config=config, **kwargs)

        # locking for thread safety
        self.threadlock = Mutex()

        self.stopRequested = False
        self.is_crosshair = True

        # Keep track of who called the refocus
        self._caller_tag = ''

    def on_activate(self):
        """ Initialisation performed during activation of the module.

        @return int: error code (0:OK, -1:error)
        """
        self._scanning_device = self.scanner()
        self._fit_logic = self.fitlogic()

        #First check the naming is valid
        if self.optimizer_stack not in self._scanning_device.get_stack_names():
            self.log.error("The optimizer stack name is not compatible with the "
                           "scanner stack names.")
        # Reads in the maximal scanning range. The unit of that scan range is micrometer!
        self.x_range = self._scanning_device.get_position_range(stack=self.optimizer_stack)[0]
        self.y_range = self._scanning_device.get_position_range(stack=self.optimizer_stack)[1]
        self._initial_pos_x = 0.
        self._initial_pos_y = 0.
        self.optim_pos_x = self._initial_pos_x
        self.optim_pos_y = self._initial_pos_y
        self.optim_sigma_x = 0.
        self.optim_sigma_y = 0.

        self._max_offset = 3.

        # Sets the current position to the center of the maximal scanning range
        self._current_x = (self.x_range[0] + self.x_range[1]) / 2
        self._current_y = (self.y_range[0] + self.y_range[1]) / 2

        # Samples per pixel (linked to the clock rate)
        self.samps_per_px = 0 #Set it to zero by default (will raise an error if not initialized properly

        ###########################
        # Fit Params and Settings #
        model, params = self._fit_logic.make_gaussianlinearoffset_model()
        self.use_custom_params = {name: False for name, param in params.items()}

        # Sets connections between signals and functions
        self._sigScanNextXyLine.connect(self._refocus_xy_line, QtCore.Qt.QueuedConnection)
        self._sigCompletedXyOptimizerScan.connect(self._set_optimized_xy_from_fit, QtCore.Qt.QueuedConnection)

        self._sigOptimizationComplete.connect(self.finish_refocus)

        self.sigGoToPoint.connect(self._move_to_start_pos, QtCore.Qt.QueuedConnection)
        return 0

    def on_deactivate(self):
        """ Reverse steps of activation

        @return int: error code (0:OK, -1:error)
        """
        return 0

    def get_scanner_count_channels(self):
        """ Get lis of counting channels from scanning device.
          @return list(str): names of counter channels
        """
        return self._scanning_device.get_scanner_count_channels()

    def set_refocus_XY_size(self, size):
        """ Set the number of pixels in the refocus image for X and Y directions

            @param int size: XY image size in pixels
        """
        self.refocus_XY_size = size
        self.sigRefocusXySizeChanged.emit()

    def start_refocus(self, initial_pos=None, caller_tag='unknown', tag='logic'):
        """ Starts the optimization scan around initial_pos

            @param list initial_pos: with the structure [float, float, float]
            @param str caller_tag:
            @param str tag:
        """
        # checking if refocus corresponding to crosshair or corresponding to initial_pos

        if isinstance(initial_pos, (np.ndarray, list, tuple)):
            self._initial_pos_x, self._initial_pos_y = initial_pos
        elif initial_pos is None:
            scpos = self._scanning_device.get_scanner_position(stack=self.optimizer_stack)
            self._initial_pos_x, self._initial_pos_y = scpos
        else:
            pass  # TODO: throw error

        # Keep track of where the start_refocus was initiated
        self._caller_tag = caller_tag
        # Set the optim_pos values to match the initial_pos values.
        # This means we can use optim_pos in subsequent steps and ensure
        # that we benefit from any completed optimization step.
        self.optim_pos_x = self._initial_pos_x
        self.optim_pos_y = self._initial_pos_y
        self.optim_sigma_x = 0.
        self.optim_sigma_y = 0.
        #
        self._xy_scan_line_count = 0
        self._optimization_step = 0

        self.stopRequested = False

        # move to the start of the first line
        x0 = self.optim_pos_x
        y0 = self.optim_pos_y
        xmin = np.clip(x0 - 0.5 * self.refocus_XY_size, self.x_range[0], self.x_range[1])
        ymin = np.clip(y0 - 0.5 * self.refocus_XY_size, self.y_range[0], self.y_range[1])

        self.sigGoToPoint.emit([xmin, ymin])
        # self._move_to_start_pos([xmin, ymin])


    def _launch_optimizer(self):
        scanner_status = self.start_scanner()
        if scanner_status < 0:
            self.sigRefocusFinished.emit(
                self._caller_tag,
                [self.optim_pos_x, self.optim_pos_y])
            return

        self._initialize_xy_refocus_image()

        self.sigRefocusStarted.emit('logic')
        self._sigScanNextXyLine.emit()

    def stop_refocus(self):
        """Stops refocus."""
        with self.threadlock:
            self.stopRequested = True

    def _initialize_xy_refocus_image(self):
        """Initialisation of the xy refocus image."""
        self._xy_scan_line_count = 0
        #
        #         # Take optim pos as center of refocus image, to benefit from any previous
        # optimization steps that have occurred.
        x0 = self.optim_pos_x
        y0 = self.optim_pos_y

        # defining position intervals for refocus
        xmin = np.clip(x0 - 0.5 * self.refocus_XY_size, self.x_range[0], self.x_range[1])
        xmax = np.clip(x0 + 0.5 * self.refocus_XY_size, self.x_range[0], self.x_range[1])
        ymin = np.clip(y0 - 0.5 * self.refocus_XY_size, self.y_range[0], self.y_range[1])
        ymax = np.clip(y0 + 0.5 * self.refocus_XY_size, self.y_range[0], self.y_range[1])

        self._X_values = np.linspace(xmin, xmax, num=self.optimizer_XY_res)
        self._Y_values = np.linspace(ymin, ymax, num=self.optimizer_XY_res)

        self._slowmove_clock = self._scanning_device.get_motion_clock_frequency()
        speed = self._scanning_device.get_motion_speed(stack=self.optimizer_stack)
        slow_pixels = int((xmax - xmin) * self._slowmove_clock / speed)
        if slow_pixels < 2:
            slow_pixels = 2
        self._return_X_values = np.linspace(xmax, xmin, num=slow_pixels)


        self.xy_refocus_image = np.zeros((
            len(self._Y_values),
            len(self._X_values),
            3))
        self.xy_refocus_image[:, :, 0] = np.full((len(self._Y_values), len(self._X_values)), self._X_values)
        y_value_matrix = np.full((len(self._X_values), len(self._Y_values)), self._Y_values)
        self.xy_refocus_image[:, :, 1] = y_value_matrix.transpose()

    def _move_to_start_pos(self, start_pos):
        """Moves the scanner from its current position to the start position of the optimizer scan.

        @param start_pos float[]: 2-point vector giving x, y position to go to.
        """
        try:
            self._scanning_device.scanner_slow_motion(start_pos, stack=self.optimizer_stack,
                                                      clear_ao_whenfinished=False,)
        except Exception as e:
            self.log.error('Error during move to starting point.')
            self.stop_refocus()
            self._sigScanNextXyLine.emit()
            return

        self._launch_optimizer()


    def _refocus_xy_line(self):
        """Scanning a line of the xy optimization image.
        This method repeats itself using the _sigScanNextXyLine
        until the xy optimization image is complete.
        """

        # stop scanning if instructed
        if self.stopRequested:
            with self.threadlock:
                self.stopRequested = False
                self.finish_refocus()
                self.sigImageUpdated.emit()
                self.sigRefocusFinished.emit(
                    [self.optim_pos_x, self.optim_pos_y])
                return

        lsx = self.xy_refocus_image[self._xy_scan_line_count, :, 0]
        lsy = self.xy_refocus_image[self._xy_scan_line_count, :, 1]

        #scan a line of the xy optimization image
        line = np.vstack((lsx, lsy))

        line_counts = self._scan_line(line, motion='fw')
        if np.any(line_counts == -1):
            self.log.error('The scan went wrong, killing the scanner.')
            self.stop_refocus()
            self._sigScanNextXyLine.emit()
            return

        lsx = self._return_X_values
        lsy = self.xy_refocus_image[self._xy_scan_line_count, 0, 1] * np.ones(lsx.shape)

        return_line = np.vstack((lsx, lsy))

        self._scan_line(return_line, motion='bw')

        self.xy_refocus_image[self._xy_scan_line_count, :, 2] = line_counts
        self.sigImageUpdated.emit()

        self._xy_scan_line_count += 1

        if self._xy_scan_line_count < np.size(self._Y_values):
            self._sigScanNextXyLine.emit()
        else:
            self._sigCompletedXyOptimizerScan.emit()

    def _scan_line(self, line, motion='fw'):
        """
        Used to combine the px scanner with the slow scanner in one function. Motion is either fw or bw. If forward,
        the motion happens px by px, and returns a line of counts at the end. Note that the fw scanning is a for loop,
        and thus is a blocking command. This mimics the behaviour of the original line scanner, which does essentially
        the same

        Line is a (2, N) matrix of coordinates.
        """

        if motion == 'fw':
            count_line = np.zeros(line.shape[1]) - 1. #Initialize to array of -1 for error checking
            for ii, px in enumerate(line.T): #For loops of arrays loop on the outermost index (i.e. rows)
                self._scanning_device.scanner_set_position(px, stack=self.optimizer_stack)
                counts, _ = self._scanning_device.read_pixel(self.samps_per_px)
                counts = counts.mean()
                count_line[ii] = counts
        elif motion == 'bw':
            self._scanning_device.move_along_line(line, stack=self.optimizer_stack)
            count_line = None

        return count_line

    def _set_optimized_xy_from_fit(self):
        """Fit the completed xy optimizer scan and set the optimized xy position."""
        fit_x, fit_y = np.meshgrid(self._X_values, self._Y_values)
        xy_fit_data = self.xy_refocus_image[:, :, 2].ravel()
        axes = np.empty((len(self._X_values) * len(self._Y_values), 2))
        axes = (fit_x.flatten(), fit_y.flatten())
        result_2D_gaus = self._fit_logic.make_twoDgaussian_fit(
            xy_axes=axes,
            data=xy_fit_data,
            estimator=self._fit_logic.estimate_twoDgaussian_MLE
        )
        # print(result_2D_gaus.fit_report())

        if result_2D_gaus.success is False:
            self.log.error('Error: 2D Gaussian Fit was not successful!.')
            self.optim_pos_x = self._initial_pos_x
            self.optim_pos_y = self._initial_pos_y
            self.optim_sigma_x = 0.
            self.optim_sigma_y = 0.
        else:
            #                @reviewer: Do we need this. With constraints not one of these cases will be possible....
            if abs(self._initial_pos_x - result_2D_gaus.best_values['center_x']) < self._max_offset and abs(
                    self._initial_pos_x - result_2D_gaus.best_values['center_x']) < self._max_offset:
                if self.x_range[0] <= result_2D_gaus.best_values['center_x'] <= self.x_range[1]:
                    if self.y_range[0] <= result_2D_gaus.best_values['center_y'] <= self.y_range[1]:
                        self.optim_pos_x = result_2D_gaus.best_values['center_x']
                        self.optim_pos_y = result_2D_gaus.best_values['center_y']
                        self.optim_sigma_x = result_2D_gaus.best_values['sigma_x']
                        self.optim_sigma_y = result_2D_gaus.best_values['sigma_y']
            else:
                self.optim_pos_x = self._initial_pos_x
                self.optim_pos_y = self._initial_pos_y
                self.optim_sigma_x = 0.
                self.optim_sigma_y = 0.

        # emit image updated signal so crosshair can be updated from this fit
        self.sigImageUpdated.emit()
        self._sigOptimizationComplete.emit()

    def finish_refocus(self):
        """ Finishes up and releases hardware after the optimizer scans."""
        self.kill_scanner()

        self.log.info(
            'Optimised from ({0:.3e},{1:.3e}) to local '
            'maximum at ({2:.3e},{3:.3e}).'.format(
                self._initial_pos_x,
                self._initial_pos_y,
                self.optim_pos_x,
                self.optim_pos_y,))

        # Signal that the optimization has finished, and "return" the optimal position along with
        # caller_tag
        self.sigRefocusFinished.emit([self.optim_pos_x, self.optim_pos_y])

    def start_scanner(self):
        """Setting up the scanner device.

        @return int: error code (0:OK, -1:error)
        """
        self.module_state.lock()

        motion_clock_status = self._scanning_device.prepare_motion_clock()

        self.samps_per_px = self._samps_per_pixel(self.integration_time*1e-3)

        self._scanning_device.create_ao_task(self.optimizer_stack)
        self._scanning_device.prepare_counters(samples_to_acquire=self.samps_per_px,
                                               counter_ai_channels=None)

        if motion_clock_status < 0:
            self.module_state.unlock()
            return -1

        return 0

    def kill_scanner(self):
        """Closing the scanner device.

        @return int: error code (0:OK, -1:error)
        """
        self._scanning_device.close_counters()
        self._scanning_device.clear_ao_task(self.optimizer_stack)
        self._scanning_device.clear_motion_clock()
        self.module_state.unlock()

    def set_position(self, tag, x=None, y=None):
        """ Set focus position.

            @param str tag: sting indicating who caused position change
            @param float x: x axis position in m
            @param float y: y axis position in m
        """
        if x is not None:
            self._current_x = x
        if y is not None:
            self._current_y = y
        self.sigPositionChanged.emit(self._current_x, self._current_y)

    def _samps_per_pixel(self, integration_time):
        return round(integration_time * self._scanning_device.get_counter_clock_frequency())
