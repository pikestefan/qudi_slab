from qtpy import QtCore
import numpy as np

from logic.generic_logic import GenericLogic
from core.util.mutex import Mutex
from core.connector import Connector
from core.statusvariable import StatusVar

"""
This first part is a copy-paste of the original code in qudi.
"""


class OldConfigFileError(Exception):
    """Exception that is thrown when an old config file is loaded."""

    def __init__(self):
        super().__init__("Old configuration file detected. Ignoring confocal history.")


class ConfocalHistoryEntry(QtCore.QObject):
    """This class contains all relevant parameters of a Confocal scan.
    It provides methods to extract, restore and serialize this data.
    """

    def __init__(self, snvm):
        """Make a confocal data setting with default values."""
        super().__init__()

        #####
        # Setting up the scanning initial parameters, and get the two stack names.
        #####
        (
            self.sampleStackName,
            self.tipStackName,
        ) = snvm._scanning_device.get_stack_names()

        self._active_stack = self.sampleStackName  # Default stack is sample

        # Get the maximum scanning ranges, and the position to voltage conversion factors, and put them in a dictionary
        self.scanning_x_range = dict()
        self.scanning_y_range = dict()
        self.scanning_x_resolution = dict()
        self.scanning_y_resolution = dict()
        self.px_time = dict()
        self.store_retrace = dict()

        for stack in [self.sampleStackName, self.tipStackName]:
            # These are the scanning ranges that will be used for the scanning
            self.scanning_x_range[stack] = [0, 1e-6]
            self.scanning_y_range[stack] = [0, 1e-6]

            self.scanning_x_resolution[stack] = 10
            self.scanning_y_resolution[stack] = 10

            # Integration time per pixel
            self.px_time[stack] = 30e-3

            self.store_retrace[stack] = False

        self._photon_samples = 0
        self.backward_pixels = (
            0  # The number is determined by the clock frequency and the bw_speed
        )
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
        """Write data back into confocal logic and pull all the necessary strings"""
        snvm.sampleStackName = self.sampleStackName
        snvm.tipStackName = self.tipStackName
        snvm._active_stack = self._active_stack
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
        """Extract all necessary data from a confocal logic and keep it for later use"""
        self.sampleStackName = snvm.sampleStackName
        self.tipStackName = snvm.tipStackName
        self._active_stack = snvm._active_stack
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
        """Give out a dictionary that can be saved via the usual means"""
        serialized = dict()
        serialized["sampleStackName"] = self.sampleStackName
        serialized["tipStackName"] = self.tipStackName
        serialized["active_stack"] = self._active_stack
        serialized["scanning_x_range"] = self.scanning_x_range
        serialized["scanning_y_range"] = self.scanning_y_range
        serialized["scanning_x_resolution"] = self.scanning_x_resolution
        serialized["scanning_y_resolution"] = self.scanning_y_resolution
        serialized["px_time"] = self.px_time
        serialized["_photon_samples"] = self._photon_samples
        serialized["backward_pixels"] = self.backward_pixels
        serialized["store_retrace"] = self.store_retrace
        serialized["invalid"] = self.invalid
        serialized["start_freq"] = self.start_freq
        serialized["stop_freq"] = self.stop_freq
        serialized["freq_resolution"] = self.freq_resolution
        serialized["mw_power"] = self.mw_power
        serialized["odmr_averages"] = self.odmr_averages
        return serialized

    def deserialize(self, serialized):
        """Restore Confocal history object from a dict"""
        if "sampleStackName" in serialized:
            self.sampleStackName = serialized["sampleStackName"]
        if "tipStackName" in serialized:
            self.tipStackName = serialized["tipStackName"]
        if "active_stack" in serialized:
            self._active_stack = serialized["active_stack"]
        if "scanning_x_range" in serialized:
            self.scanning_x_range = serialized["scanning_x_range"]
        if "scanning_y_range" in serialized:
            self.scanning_y_range = serialized["scanning_y_range"]
        if "scanning_x_resolution" in serialized:
            self.scanning_x_resolution = serialized["scanning_x_resolution"]
        if "scanning_y_resolution" in serialized:
            self.scanning_y_resolution = serialized["scanning_y_resolution"]
        if "px_time" in serialized:
            self.px_time = serialized["px_time"]
        if "_photon_samples" in serialized:
            self._photon_samples = serialized["_photon_samples"]
        if "backward_pixels" in serialized:
            self.backward_pixels = serialized["backward_pixels"]
        if "store_retrace" in serialized:
            self.store_retrace = serialized["store_retrace"]
        if "invalid" in serialized:
            self.invalid = serialized["invalid"]
        if "start_freq" in serialized:
            self.start_freq = serialized["start_freq"]
        if "stop_freq" in serialized:
            self.stop_freq = serialized["stop_freq"]
        if "freq_resolution" in serialized:
            self.freq_resolution = serialized["freq_resolution"]
        if "mw_power" in serialized:
            self.mw_power = serialized["mw_power"]
        if "odmr_averages" in serialized:
            self.odmr_averages = serialized["odmr_averages"]


class SnvmLogic(GenericLogic):

    doublescanner = Connector(interface="SnvmScannerInterface")
    odmrscanner = Connector(interface="MicrowaveInterface")
    pulser = Connector(interface="PulserInterface")
    savelogic = Connector(interface="HDF5SaveLogic")
    optimizer_logic = Connector(interface="OptimizerLogic")
    master_pulselogic = Connector(interface="MasterPulse")

    slow_motion_clock_rate = StatusVar("slow_motion_clock_rate", 10)
    backward_speed_conf = StatusVar("slow_motion_speed_conf", 1)
    backward_speed_snvm = StatusVar("slow_motion_speed_snvm", 1)
    max_history_length = StatusVar(default=10)

    # Pulsed ODMR settings
    podmr_active = StatusVar("podmr_active", False)
    podmr_start_delay = StatusVar("podmr_start_delay", 0)
    podmr_laser_init = StatusVar("podmr_laser_init", 1500e-9)
    podmr_laser_read = StatusVar("podmr_laser_read", 500e-9)
    podmr_init_delay = StatusVar("podmr_init_delay", 350e-9)
    podmr_read_delay = StatusVar("podmr_read_delay", 350e-9)
    podmr_apd_delay = StatusVar("podmr_apd_delay", 0)
    podmr_apd_read = StatusVar("podmr_apd_read", 500e-9)
    podmr_final_delay = StatusVar("podmr_final_delay", 0)
    podmr_clkrate = StatusVar("podmr_clkrate", 1250e6)
    podmr_pipulse = StatusVar("pdomr_pipulse", 500e-9)

    # Initialize the optimization-while-scanning settings
    optimize_while_scanning = StatusVar("optimize_while_scanning", False)
    every_N_pixels = StatusVar("every_N_pixels", 10)

    # signals
    signal_start_snvm = QtCore.Signal()
    signal_continue_snvm = QtCore.Signal()
    signal_start_confocal = QtCore.Signal()
    signal_continue_confocal = QtCore.Signal()
    signal_stop_scan = QtCore.Signal()
    signal_start_optimizer = QtCore.Signal()

    signal_snvm_image_updated = QtCore.Signal()
    signal_xy_image_updated = QtCore.Signal()
    signal_freq_px_acquired = QtCore.Signal(
        int
    )  # Emits the row of the temporary data matrix to store the odmr data
    signal_odmr_trace_updated = QtCore.Signal(
        int
    )  # Used to refresh the odmr plot in the scanning GUI
    signal_odmr_line_acquired = (
        QtCore.Signal()
    )  # Emitted only if plotting ODMR at the end of scan is requested
    signal_xy_px_acquired = QtCore.Signal()
    signal_scan_finished = QtCore.Signal(
        bool
    )  # Emits True if the scan was snvm, False otherwise

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
        self._master_pulselogic = self.master_pulselogic()
        self._pulser = self.pulser()

        self._master_pulselogic.cw()
        self.set_slowmotion_clockrate(self.slow_motion_clock_rate)
        self.set_motion_speed(self.backward_speed_conf, stack="tip")
        self.set_motion_speed(self.backward_speed_snvm, stack="sample")

        self.history = []
        for i in reversed(range(1, self.max_history_length)):
            try:
                new_history_item = ConfocalHistoryEntry(self)
                new_history_item.deserialize(
                    self._statusVariables["history_{0}".format(i)]
                )
                self.history.append(new_history_item)
            except KeyError:
                pass
            except OldConfigFileError:
                self.log.warning(
                    "Old style config file detected. History {0} ignored.".format(i)
                )
            except:
                self.log.warning("Restoring history {0} failed.".format(i))
        try:
            new_state = ConfocalHistoryEntry(self)
            new_state.deserialize(self._statusVariables["history_0"])
            new_state.restore(self)
        except:
            new_state = ConfocalHistoryEntry(self)
            new_state.restore(self)
        finally:
            self.history.append(new_state)

        self.history_index = len(self.history) - 1

        # Get the maximum scanning ranges and the volt-position mapping from the scanning hardware.
        self.x_maxrange = dict()
        self.y_maxrange = dict()
        self.position_to_voltage_arrays = dict()
        for stack in [self.sampleStackName, self.tipStackName]:
            x_range, y_range = self._scanning_device.get_position_range(stack=stack)
            x_volt_range, y_volt_range = self._scanning_device.get_voltage_range(
                stack=stack
            )
            self.x_maxrange[stack] = x_range
            self.y_maxrange[stack] = y_range

            self.position_to_voltage_arrays[stack] = [
                (x_volt_range[1] - x_volt_range[0]) / (x_range[1] - x_range[0]),
                (y_volt_range[1] - y_volt_range[0]) / (y_range[1] - y_range[0]),
            ]

        # These are the indices which will be used to scan through the arrays of the frequencies and position pairs
        self._x_scanning_index = 0
        self._y_scanning_index = 0
        self._tot_xy_scanning_index = 0
        self._x_index_step = (
            1  # This coefficient is decided to decide the direction of the x scanning
        )
        self._freq_scanning_index = 0
        self._odmr_rep_index = 0  # To keep track of the averages
        self._is_retracing = False
        self.stopRequested = False

        # Settings for the pulsed odmr
        self.use_podmr = False

        # Initalize the attributes that will be the scan data containers
        self._x_scanning_axis = None
        self._y_scanning_axis = None
        self.xy_scan_matrix = None
        self.snvm_matrix = None
        self._temp_afm_matrix = (
            None  # Matrix used to store the AFM values while scanning the ODMR.
        )
        self.temp_freq_matrix = (
            None  # This matrix is used to store the ODMR traces to be averaged.
        )
        self.xy_scan_matrix_retrace = None
        self.snvm_matrix_retrace = None
        self.freq_axis = None
        self._xy_matrix_width = None  # This variable is used only to keep the code clean and reduce the calls to np.shape
        self._freq_axis_length = None  # Same here
        self.average_odmr_trace = None

        self._snvm_active = False

        self._curr_snvm_image = []
        self._curr_afm_image = []
        self._curr_cfc_image = []

        self.completed_pixels_matrix = (
            None  # Matrix to be used as a mask for the plotting
        )
        self.pxbypx_odmr = True

        # Now connect all the signals
        self.signal_continue_snvm.connect(
            self.continue_snvm_scanning, QtCore.Qt.QueuedConnection
        )
        self.signal_stop_scan.connect(self.stop_scanning, QtCore.Qt.QueuedConnection)
        self.signal_xy_px_acquired.connect(
            self.move_to_xy_pixel, QtCore.Qt.QueuedConnection
        )
        self.signal_freq_px_acquired.connect(
            self.move_to_freq_pixel, QtCore.Qt.QueuedConnection
        )
        self.signal_continue_confocal.connect(
            self.continue_confocal_scanning, QtCore.Qt.QueuedConnection
        )
        self.signal_goto_start.connect(self._go_to_point, QtCore.Qt.QueuedConnection)
        self.signal_start_optimizer.connect(
            self.start_optimizer, QtCore.Qt.QueuedConnection
        )
        self._optimizer_logic.sigRefocusFinished.connect(self._optimization_complete)

    def on_deactivate(self):
        """Reverse steps of activation
        @return int: error code (0:OK, -1:error)
        """
        closing_state = ConfocalHistoryEntry(self)
        closing_state.snapshot(self)
        self.history.append(closing_state)
        histindex = 0
        for state in reversed(self.history):
            self._statusVariables["history_{0}".format(histindex)] = state.serialize()
            histindex += 1

        return 0

    def _prepare_data_matrices(self):
        # Clip the ranges if they are out of bound
        self.check_xy_ranges()

        # Generate axes and matrices to store the data
        x_axis = np.linspace(
            self.scanning_x_range[self._active_stack][0],
            self.scanning_x_range[self._active_stack][1],
            self.scanning_x_resolution[self._active_stack],
        )
        y_axis = np.linspace(
            self.scanning_y_range[self._active_stack][0],
            self.scanning_y_range[self._active_stack][1],
            self.scanning_y_resolution[self._active_stack],
        )

        # Now generate the matrices to store the data
        xy_scan_matrix = np.zeros((len(y_axis), len(x_axis)), dtype=np.float64)
        # FIXME: for now the stack scanner is the one that's assumed to have the ESR sequence. Maybe consider a flexible
        #  way of doing this
        if self._snvm_active:
            step_number = 1 + round(
                (self.stop_freq - self.start_freq) / self.freq_resolution
            )
            freq_axis = np.linspace(self.start_freq, self.stop_freq, step_number)
            snvm_matrix = np.zeros(
                (
                    xy_scan_matrix.shape[0],
                    xy_scan_matrix.shape[1],
                    len(freq_axis),
                    self.odmr_averages,
                ),
                dtype=np.float64,
            )
            temp_freq_matrix = np.full(
                (self.odmr_averages, len(freq_axis)), self.invalid
            )
            temp_afm_matrix = np.copy(temp_freq_matrix)
            average_odmr_trace = np.zeros((temp_freq_matrix.shape[1],))
        else:
            snvm_matrix = np.copy(xy_scan_matrix[:, :, np.newaxis])
            freq_axis = None
            temp_freq_matrix = None
            average_odmr_trace = None
            temp_afm_matrix = None

        if self.store_retrace[self._active_stack]:
            xy_scan_matrix_retrace = np.copy(xy_scan_matrix)
            snvm_matrix_retrace = np.copy(snvm_matrix)
        else:
            xy_scan_matrix_retrace = np.zeros(xy_scan_matrix.shape)
            snvm_matrix_retrace = np.zeros(snvm_matrix.shape)

        self._tot_xy_scanning_index = 0
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
        self._freq_axis_length = (
            len(self.freq_axis) if isinstance(self.freq_axis, np.ndarray) else None
        )
        self.completed_pixels_matrix = np.zeros((2,) + xy_scan_matrix.shape, dtype=int)

    def check_xy_ranges(self):
        """
        Check that the requested scanning ranges are within the maximum scanning ranges set by the hardware.
        If they are out of bounds, clip the values.
        """
        stk = self._active_stack
        curr_x_minrange, curr_x_maxrange = self.x_maxrange[stk]
        curr_y_minrange, curr_y_maxrange = self.y_maxrange[stk]
        # TODO: emit a signal when the clipping happens, and update the GUI limits accordingly
        if not (
            (curr_x_minrange <= self.scanning_x_range[stk][0] <= curr_x_maxrange)
            or not (curr_x_minrange <= self.scanning_x_range[stk][1] <= curr_x_maxrange)
        ):
            self.scanning_x_range[stk] = np.clip(
                self.scanning_x_range[stk], curr_x_minrange, curr_x_maxrange
            )
            self.log.warning(
                "x scanning range limits are out of bounds, clipped back to the maximum values."
            )

        if not (
            (curr_y_minrange <= self.scanning_y_range[stk][0] <= curr_y_maxrange)
            or not (curr_y_minrange <= self.scanning_y_range[stk][1] <= curr_y_maxrange)
        ):
            self.scanning_y_range[stk] = np.clip(
                self.scanning_y_range[stk], curr_y_minrange, curr_y_maxrange
            )
            self.log.warning(
                "y scanning range limits are out of bounds, clipped back to the maximum values."
            )
        return 0

    def prepare_devices(self):
        """
        Initalize the scanning device, prior to starting the scan.
        """

        self._scanning_device.module_state.lock()
        self._photon_samples = self.pxtime_to_samples()
        if self._snvm_active:
            analog_channels = self._scanning_device.get_ai_counter_channels(
                stack_name=self._active_stack
            )
        else:
            analog_channels = None
        self._scanning_device.prepare_counters(
            samples_to_acquire=self._photon_samples, counter_ai_channels=analog_channels
        )

        self._scanning_device.create_ao_task(self._active_stack)
        if self.store_retrace[self._active_stack] is False:
            self._scanning_device.prepare_motion_clock()
            clk_freq = self._scanning_device.get_motion_clock_frequency()
            speed = self._scanning_device.get_motion_speed(self._active_stack)
            self.backward_pixels = int(
                (
                    (self._x_scanning_axis.max() - self._x_scanning_axis.min())
                    / speed
                )
                * clk_freq
            )
            if self.backward_pixels < 2:
                self.backward_pixels = 2

            self.x_backward = np.zeros((2, self.backward_pixels))
            self.x_backward[0] = np.linspace(
                self._x_scanning_axis.max(),
                self._x_scanning_axis.min(),
                self.backward_pixels,
            )

        if self._snvm_active:
            self._odmrscanner.module_state.lock()
            self._odmrscanner.set_IQmod(self.podmr_active)
            try:
                pass
                # FIXME: look into how to operate the srs in list mode
                # self._odmrscanner.set_list(frequency=self.freq_axis, power=self.mw_power)
            except:
                self.log.error(
                    "Failed loading the frequency axis into the ODMR scanner. Aborted execution."
                )
                # self._scanning_device.module_state.unlock()
                # self._odmrscanner.module_state.unlock()

    def load_pulse_sequence(self):
        # Convert all the natural units in ns and GHz
        clkrate = round(self.podmr_clkrate * 1e-9)

        start_delay = self.podmr_start_delay * 1e9
        laser_init = self.podmr_laser_init * 1e9
        laser_read = self.podmr_laser_read * 1e9
        init_delay = self.podmr_init_delay * 1e9
        read_delay = self.podmr_read_delay * 1e9
        apd_delay = self.podmr_apd_delay * 1e9
        apd_read = self.podmr_apd_read * 1e9
        final_delay = self.podmr_final_delay * 1e9

        pipulse = self.podmr_pipulse * 1e9
        init_center = start_delay + laser_init / 2
        mw_center = init_center + laser_init / 2 + init_delay + pipulse / 2
        read_center = mw_center + pipulse / 2 + read_delay + laser_read / 2
        apd_center = mw_center + pipulse / 2 + read_delay + apd_delay + apd_read / 2

        tot_time = (
            start_delay
            + laser_init
            + init_delay
            + pipulse
            + read_delay
            + max(laser_read + final_delay, apd_delay + apd_read)
        )

        tot_samples = self._pulser.waveform_padding(tot_time * clkrate)
        timeax = np.linspace(0, tot_samples / clkrate, tot_samples)
        blank_signal = np.zeros_like(timeax)
        high_signal = blank_signal + 1

        mw_waveform = self._pulser.box_envelope(
            timeax, np.array([[mw_center, pipulse]])
        )
        laser_waveform = self._pulser.box_envelope(
            timeax, np.array([[init_center, laser_init], [read_center, laser_read]])
        )
        apd_waveform = self._pulser.box_envelope(
            timeax, np.array([[apd_center, apd_read]])
        )

        card_idx = 1
        self._pulser.stop_all()
        self._pulser.set_sample_rate(card_idx, int(self.podmr_clkrate))

        steps = 2
        do_map = {1: [0, 1, 2]}
        for step in range(steps):
            if step == 0:
                self._pulser.load_waveform(
                    iq_dictionary={"i_chan": mw_waveform, "q_chan": blank_signal},
                    digital_pulses={
                        "laser": laser_waveform,
                        "apd_sig": apd_waveform,
                        "apd_read": blank_signal,
                    },
                    digital_output_map=do_map,
                )
            else:
                # This waveform is going to be used when the optimizer is called
                self._pulser.load_waveform(
                    iq_dictionary={"i_chan": blank_signal, "q_chan": blank_signal},
                    digital_pulses={
                        "laser": high_signal,
                        "apd_sig": high_signal,
                        "apd_read": blank_signal,
                    },
                    digital_output_map=do_map,
                )

        self._pulser.configure_ORmask(card_idx, None)
        self._pulser.configure_ANDmask(card_idx, None)

        loops = np.ones(steps, dtype=np.int64)

        stop_condition_list = np.array(
            [
                0x40000000,
            ]
            * steps,
            dtype=np.int64,
        )

        self._pulser.load_sequence(
            loops_list=loops,
            stop_condition_list=stop_condition_list,
        )

        self._pulser.start_card(card_idx)
        self._pulser.arm_trigger(card_idx)

        self._pulser.send_software_trig(card_idx)

    def start_snvm_scanning(self):

        self.stopRequested = False

        self.module_state.lock()

        self._active_stack = self.sampleStackName
        self._snvm_active = True
        self._prepare_data_matrices()

        self._initialize_scanning_statuses()
        self.prepare_devices()

        # This part needs to stay outside prepare_devices(), otherwise load_pulse_sequence() is going to be called
        # every time the optimizer is finished.
        if self.podmr_active:
            self.load_pulse_sequence()
        else:
            # With this line, make sure that the AWG is playing in CW mode if normal ODMR is requested.
            self._master_pulselogic.cw()

        self.signal_snvm_image_updated.emit()
        self.signal_snvm_initialized.emit()

        self.go_to_point(
            [
                self._x_scanning_axis[self._x_scanning_index],
                self._y_scanning_axis[self._y_scanning_index],
            ],
            stack=self._active_stack,
            clear_ao_whenfinished=False,
            caller="logic",
        )
        # FIXME: look into how to operate the srs in list mode
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

        self.go_to_point(
            [
                self._x_scanning_axis[self._x_scanning_index],
                self._y_scanning_axis[self._y_scanning_index],
            ],
            stack=self._active_stack,
            clear_ao_whenfinished=False,
            caller="logic",
        )
        self.signal_continue_confocal.emit()

    def start_optimizer(self):
        if self.podmr_active:
            self._pulser.send_software_trig(1)
        self._scanning_device.pause_tasks()
        self._optimizer_logic.start_refocus()

    def continue_snvm_scanning(self):
        acquire_data = (
            False
            if (self.store_retrace[self._active_stack] is False)
            and (self._is_retracing is True)
            else True
        )
        if acquire_data:
            # If the index of the ODMR is less than the averages, keep acquiring
            if self._odmr_rep_index < self.odmr_averages:
                counts, ainput = self.acquire_pixel()
                counts = counts[0] / self.px_time[self._active_stack]
                self.temp_freq_matrix[
                    self._odmr_rep_index, self._freq_scanning_index
                ] = counts
                self._temp_afm_matrix[
                    self._odmr_rep_index, self._freq_scanning_index
                ] = ainput.mean()

                if self._odmr_rep_index > 0 and self.pxbypx_odmr:
                    self.average_odmr_trace = np.nanmean(self.temp_freq_matrix, axis=0)

                if self._is_retracing and self.store_retrace[self._active_stack]:
                    self.snvm_matrix_retrace[
                        self._y_scanning_index,
                        self._x_scanning_index,
                        self._freq_scanning_index,
                        self._odmr_rep_index,
                    ] = counts
                else:
                    self.snvm_matrix[
                        self._y_scanning_index,
                        self._x_scanning_index,
                        self._freq_scanning_index,
                        self._odmr_rep_index,
                    ] = counts
                self.signal_freq_px_acquired.emit(self._odmr_rep_index)
                if self.pxbypx_odmr:
                    self.signal_odmr_trace_updated.emit(self._odmr_rep_index)

            # Else, the OMDR acquisition for the pixel has finished. Store the data and ask for the next pixel. If also the
            # Scanning is done, tell that the scanning has finished.
            else:
                if self._is_retracing and self.store_retrace[self._active_stack]:
                    # FIXME: I am not acquiring and storing properly the analog input, find a way after basic debugging done
                    self.xy_scan_matrix_retrace[
                        self._y_scanning_index, self._x_scanning_index
                    ] = self._temp_afm_matrix.mean()
                else:
                    self.xy_scan_matrix[
                        self._y_scanning_index, self._x_scanning_index
                    ] = self._temp_afm_matrix.mean()

                self.completed_pixels_matrix[
                    int(self._is_retracing),
                    self._y_scanning_index,
                    self._x_scanning_index,
                ] = 1
                if not self.pxbypx_odmr:
                    if self._is_retracing and self.store_retrace[self._active_stack]:
                        thetrace = self.snvm_matrix_retrace[
                            self._y_scanning_index, self._x_scanning_index
                        ]
                    else:
                        thetrace = self.snvm_matrix[
                            self._y_scanning_index, self._x_scanning_index
                        ]
                    self.last_odmr_trace = thetrace.mean(axis=-1)
                    self.signal_odmr_line_acquired.emit()

                # The ODMR sequence has finished. Update the indices accordingly
                # self._x_scanning_index += self._x_index_step
                self._odmr_rep_index = 0

                self.temp_freq_matrix[:] = self.invalid
                self.average_odmr_trace[:] = self.invalid

                self.signal_snvm_image_updated.emit()
                self.signal_xy_px_acquired.emit()

    def continue_confocal_scanning(self):
        stk = self._active_stack
        acquire_data = (
            False
            if (self.store_retrace[stk] is False) and (self._is_retracing is True)
            else True
        )
        if acquire_data:
            counts, _ = self.acquire_pixel()
            counts = counts[0] / self.px_time[stk]
            if self._is_retracing and self.store_retrace[stk]:
                self.xy_scan_matrix_retrace[
                    self._y_scanning_index, self._x_scanning_index
                ] = counts
            else:
                self.xy_scan_matrix[
                    self._y_scanning_index, self._x_scanning_index
                ] = counts
            self.completed_pixels_matrix[
                int(self._is_retracing), self._y_scanning_index, self._x_scanning_index
            ] = 1
        self.signal_xy_image_updated.emit()
        # self._x_scanning_index += self._x_index_step
        self.signal_xy_px_acquired.emit()

    def acquire_pixel(self):
        data = self._scanning_device.read_pixel(self._photon_samples)
        return data

    def stop_scanning(self):
        if self.stopRequested:
            with self.threadlock:
                self._scanning_device.close_counters()
                self.stopRequested = False
                self.stop_xy_scanner()
                if self._snvm_active:
                    self.stop_freq_scanner()
                    if self.podmr_active:
                        self._master_pulselogic.cw()
                if not self.store_retrace[self._active_stack]:
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
            self.log.exception("Could not close the scanning tasks.")
        try:
            self._scanning_device.clear_ao_task(self._active_stack)
        except:
            self.log.exception("Could not clear the ao task.")
        try:
            if self._scanning_device.module_state() == "locked":
                self._scanning_device.module_state.unlock()
        except Exception as e:
            self.log.exception("Could not unlock scanning device.")

        return 0

    def stop_freq_scanner(self):
        self._odmrscanner.off()
        self._odmrscanner.set_IQmod(self.podmr_active)
        try:
            self._odmrscanner.module_state.unlock()
        except Exception as e:
            self.log.exception("Could not unlock scanning device.")

        return 0

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

    def move_to_xy_pixel(self):

        self._update_indices()

        if not self.stopRequested:
            if self._is_retracing and not self.store_retrace[self._active_stack]:
                self._clocked_retrace()

                self._x_scanning_index = 0
                self._y_scanning_index += 1
                self._is_retracing = False
                self._x_index_step = 1

            new_x_pos = self._x_scanning_axis[self._x_scanning_index]
            new_y_pos = self._y_scanning_axis[self._y_scanning_index]
            self._scanning_device.scanner_set_position(
                [new_x_pos, new_y_pos], stack=self._active_stack
            )

            if self._snvm_active:
                # Check if an optimization needs to be done
                if (
                    self.optimize_while_scanning
                    and self._tot_xy_scanning_index % self.every_N_pixels == 0
                ):
                    self.signal_start_optimizer.emit()
                else:
                    self.signal_continue_snvm.emit()
            else:
                self.signal_continue_confocal.emit()
        else:
            self.signal_stop_scan.emit()

    def _update_indices(self):
        self._x_scanning_index += self._x_index_step

        # This clause is activated when the last x index is reached.
        if not self._is_retracing and self._x_scanning_index == len(
            self._x_scanning_axis
        ):
            self._is_retracing = True
            self._x_index_step = -1  # If retracing starts, flip the stepping direction
            self._x_scanning_index += self._x_index_step

            # If the maximum y index has been reached and the re-trace is not stored, the scan is done.
            if (
                self._y_scanning_index == len(self._y_scanning_axis) - 1
            ) and not self.store_retrace[self._active_stack]:
                self.stopRequested = True

        # When the index becomes negative, it means that the re-trace has stepped back to the leftmost x position.
        elif self._is_retracing and self._x_scanning_index < 0:
            self._is_retracing = False
            self._x_index_step = 1  # The stepping direction now points forward.
            self._x_scanning_index += self._x_index_step
            self._y_scanning_index += 1

            # If the maximum y index has been reached, the scan is done.
            if self._y_scanning_index == len(self._y_scanning_axis):
                self.stopRequested = True

        self._tot_xy_scanning_index += 1


    def _optimization_complete(self, coords):
        if self._snvm_active:
            self.go_to_point(coords, stack=self._optimizer_logic.optimizer_stack)

            if self._odmrscanner.module_state() == "locked":
                self._odmrscanner.module_state.unlock()

            if self._scanning_device.module_state() == "locked":
                self._scanning_device.module_state.unlock()

            if self.podmr_active:
                self._pulser.send_software_trig(1)
            self.prepare_devices()
            if not self.stopRequested:
                self.signal_continue_snvm.emit()
            else:
                self.signal_stop_scan.emit()

    def pxtime_to_samples(self):
        return round(
            self.px_time[self._active_stack]
            * self._scanning_device.get_counter_clock_frequency()
        )

    def _clocked_retrace(self):

        self.x_backward[1] = self._y_scanning_axis[self._y_scanning_index]

        self._scanning_device.move_along_line(
            position_array=self.x_backward, stack=self._active_stack
        )

    def _store_data_matrices(self, is_snvm_data=True):
        if is_snvm_data:
            self._curr_snvm_image = [
                self._x_scanning_axis,
                self._y_scanning_axis,
                self.freq_axis,
                np.arange(self.odmr_averages),
                self.snvm_matrix,
                self.snvm_matrix_retrace,
            ]
            self._curr_afm_image = [
                self._x_scanning_axis,
                self._y_scanning_axis,
                self.xy_scan_matrix,
                self.xy_scan_matrix_retrace,
            ]
        else:
            self._curr_cfc_image = [
                self._x_scanning_axis,
                self._y_scanning_axis,
                self.xy_scan_matrix,
                self.xy_scan_matrix_retrace,
            ]

    def go_to_point(
        self, xy_coord, stack=None, clear_ao_whenfinished=True, caller="logic"
    ):
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
        self._scanning_device.scanner_slow_motion(
            xy_coord, stack=stack, clear_ao_whenfinished=clear_ao_whenfinished
        )

        self.signal_moved_to_point.emit(caller)

    def get_xy_image_range(self, multiplier=1):
        """
        Multiplier is an optional parameter to convert the range to the desired units
        """
        return [
            [
                self._x_scanning_axis[0] * multiplier,
                self._x_scanning_axis[-1] * multiplier,
            ],
            [
                self._y_scanning_axis[0] * multiplier,
                self._y_scanning_axis[-1] * multiplier,
            ],
        ]

    def get_xy_step_size(self, multiplier=1):
        return [
            (self._x_scanning_axis[1] - self._x_scanning_axis[0]) * multiplier,
            (self._y_scanning_axis[1] - self._y_scanning_axis[0]) * multiplier,
        ]

    def get_stack_names(self):
        return self._scanning_device.get_stack_names()

    def get_slowmotion_clockrate(self):
        return self._scanning_device.get_motion_clock_frequency()

    def get_maxranges(self):
        return self._scanning_device.get_position_range(
            "sample"
        ), self._scanning_device.get_position_range("tip")

    def set_slowmotion_clockrate(self, clockrate):
        self.slow_motion_clock_rate = clockrate
        self._scanning_device.set_motion_clock_frequency(clockrate)

    def get_motion_speed(self, stack):
        return self._scanning_device.get_motion_speed(stack)

    def set_motion_speed(self, speed, stack):
        if stack == "tip":
            self.backward_speed_conf = speed
        elif stack == "sample":
            self.backward_speed_snvm = speed
        else:
            raise ValueError(f"Trying to set speed for stack: {stack}")
        self._scanning_device.set_motion_speed(speed * 1e-6, stack)

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
            (
                xaxsnvm,
                yaxsnvm,
                freqax,
                averages,
                snvm_trace,
                snvm_retrace,
            ) = self._curr_snvm_image
            xaxafm, yaxafm, afm_trace, afm_retrace = self._curr_afm_image
            data = {
                "snvm/x_axis": xaxsnvm,
                "snvm/y_axis": yaxsnvm,
                "snvm/frequency_axis": freqax,
                "snvm/averages": averages,
                "snvm/snvm_trace": snvm_trace,
                "snvm/snvm_retrace": snvm_retrace,
                "afm/x_axis": xaxafm,
                "afm/y_axis": yaxafm,
                "afm/afm_trace": afm_trace,
                "afm/afm_retrace": afm_retrace,
            }
            self._savelogic.save_hdf5_data(data)
        else:
            self.log.info("The SNVM arrays are empty")

    def save_confocal(self):
        if len(self._curr_cfc_image) > 0:
            xax, yax, cfc_trace, cfc_retrace = self._curr_cfc_image
            data = {
                "confocal/x_axis": xax,
                "confocal/y_axis": yax,
                "confocal/confocal_trace": cfc_trace,
                "confocal/confocal_retrace": cfc_retrace,
            }
            self._savelogic.save_hdf5_data(data)
        else:
            self.log.info("The confocal arrays are empty")


# class SnvmLogic_resting(GenericLogic):
#
#     doublescanner = Connector(interface="SnvmScannerInterface")
#     odmrscanner = Connector(interface="MicrowaveInterface")
#     pulser = Connector(interface="PulserInterface")
#     savelogic = Connector(interface="HDF5SaveLogic")
#     optimizer_logic = Connector(interface="OptimizerLogic")
#     master_pulselogic = Connector(interface="MasterPulse")
#
#     slow_motion_clock_rate = StatusVar("slow_motion_clock_rate", 10)
#     backward_speed_conf = StatusVar("slow_motion_speed_conf", 1)
#     backward_speed_snvm = StatusVar("slow_motion_speed_snvm", 1)
#     max_history_length = StatusVar(default=10)
#
#     # Pulsed ODMR settings
#     podmr_active = StatusVar("podmr_active", False)
#     podmr_start_delay = StatusVar("podmr_start_delay", 0)
#     podmr_laser_init = StatusVar("podmr_laser_init", 1500e-9)
#     podmr_laser_read = StatusVar("podmr_laser_read", 500e-9)
#     podmr_init_delay = StatusVar("podmr_init_delay", 350e-9)
#     podmr_read_delay = StatusVar("podmr_read_delay", 350e-9)
#     podmr_apd_delay = StatusVar("podmr_apd_delay", 0)
#     podmr_apd_read = StatusVar("podmr_apd_read", 500e-9)
#     podmr_final_delay = StatusVar("podmr_final_delay", 0)
#     podmr_clkrate = StatusVar("podmr_clkrate", 1250e6)
#     podmr_pipulse = StatusVar("pdomr_pipulse", 500e-9)
#
#     # Initialize the optimization-while-scanning settings
#     optimize_while_scanning = StatusVar("optimize_while_scanning", False)
#     every_N_pixels = StatusVar("every_N_pixels", 10)
#
#     # signals
#     signal_start_snvm = QtCore.Signal()
#     signal_continue_snvm = QtCore.Signal()
#     signal_start_confocal = QtCore.Signal()
#     signal_continue_confocal = QtCore.Signal()
#     signal_stop_scan = QtCore.Signal()
#     signal_start_optimizer = QtCore.Signal()
#
#     signal_snvm_image_updated = QtCore.Signal()
#     signal_xy_image_updated = QtCore.Signal()
#     signal_freq_px_acquired = QtCore.Signal(
#         int
#     )  # Emits the row of the temporary data matrix to store the odmr data
#     signal_odmr_trace_updated = QtCore.Signal(
#         int
#     )  # Used to refresh the odmr plot in the scanning GUI
#     signal_odmr_line_acquired = (
#         QtCore.Signal()
#     )  # Emitted only if plotting ODMR at the end of scan is requested
#     signal_xy_px_acquired = QtCore.Signal()
#     signal_scan_finished = QtCore.Signal(
#         bool
#     )  # Emits True if the scan was snvm, False otherwise
#
#     signal_goto_start = QtCore.Signal(list, str, bool, str)
#     signal_moved_to_point = QtCore.Signal(str)
#
#     signal_snvm_initialized = QtCore.Signal()
#     signal_confocal_initialized = QtCore.Signal()
#
#     def __init__(self, config, **kwargs):
#         super().__init__(config=config, **kwargs)
#
#         self.threadlock = Mutex()
#         self.stopRequested = False
#
#     def on_activate(self):
#         self._scanning_device = self.doublescanner()
#         self._odmrscanner = self.odmrscanner()
#         self._savelogic = self.savelogic()
#         self._optimizer_logic = self.optimizer_logic()
#         self._master_pulselogic = self.master_pulselogic()
#         self._pulser = self.pulser()
#
#         self._master_pulselogic.cw()
#         self.set_slowmotion_clockrate(self.slow_motion_clock_rate)
#         self.set_motion_speed(self.backward_speed_conf, stack="tip")
#         self.set_motion_speed(self.backward_speed_snvm, stack="sample")
#
#         # FIXME: internal variable for testing, remove once finished
#         self.slowstepping = True
#         # FIXME: end
#
#         self.history = []
#         for i in reversed(range(1, self.max_history_length)):
#             try:
#                 new_history_item = ConfocalHistoryEntry(self)
#                 new_history_item.deserialize(
#                     self._statusVariables["history_{0}".format(i)]
#                 )
#                 self.history.append(new_history_item)
#             except KeyError:
#                 pass
#             except OldConfigFileError:
#                 self.log.warning(
#                     "Old style config file detected. History {0} ignored.".format(i)
#                 )
#             except:
#                 self.log.warning("Restoring history {0} failed.".format(i))
#         try:
#             new_state = ConfocalHistoryEntry(self)
#             new_state.deserialize(self._statusVariables["history_0"])
#             new_state.restore(self)
#         except:
#             new_state = ConfocalHistoryEntry(self)
#             new_state.restore(self)
#         finally:
#             self.history.append(new_state)
#
#         self.history_index = len(self.history) - 1
#
#         # Get the maximum scanning ranges and the volt-position mapping from the scanning hardware.
#         self.x_maxrange = dict()
#         self.y_maxrange = dict()
#         self.position_to_voltage_arrays = dict()
#         for stack in [self.sampleStackName, self.tipStackName]:
#             x_range, y_range = self._scanning_device.get_position_range(stack=stack)
#             x_volt_range, y_volt_range = self._scanning_device.get_voltage_range(
#                 stack=stack
#             )
#             self.x_maxrange[stack] = x_range
#             self.y_maxrange[stack] = y_range
#
#             self.position_to_voltage_arrays[stack] = [
#                 (x_volt_range[1] - x_volt_range[0]) / (x_range[1] - x_range[0]),
#                 (y_volt_range[1] - y_volt_range[0]) / (y_range[1] - y_range[0]),
#             ]
#
#         # These are the indices which will be used to scan through the arrays of the frequencies and position pairs
#         self._x_scanning_index = 0
#         self._y_scanning_index = 0
#         self._tot_xy_scanning_index = 0
#         self._x_index_step = (
#             1  # This coefficient is decided to decide the direction of the x scanning
#         )
#         self._freq_scanning_index = 0
#         self._odmr_rep_index = 0  # To keep track of the averages
#         self._is_retracing = False
#         self.stopRequested = False
#
#         # Settings for the pulsed odmr
#         self.use_podmr = False
#
#         # Initalize the attributes that will be the scan data containers
#         self._x_scanning_axis = None
#         self._y_scanning_axis = None
#         self.xy_scan_matrix = None
#         self.snvm_matrix = None
#         self._temp_afm_matrix = (
#             None  # Matrix used to store the AFM values while scanning the ODMR.
#         )
#         self.temp_freq_matrix = (
#             None  # This matrix is used to store the ODMR traces to be averaged.
#         )
#         self.xy_scan_matrix_retrace = None
#         self.snvm_matrix_retrace = None
#         self.freq_axis = None
#         self._xy_matrix_width = None  # This variable is used only to keep the code clean and reduce the calls to np.shape
#         self._freq_axis_length = None  # Same here
#         self.average_odmr_trace = None
#
#         self._snvm_active = False
#
#         self._curr_snvm_image = []
#         self._curr_afm_image = []
#         self._curr_cfc_image = []
#
#         self.completed_pixels_matrix = (
#             None  # Matrix to be used as a mask for the plotting
#         )
#         self.pxbypx_odmr = True
#
#         # Now connect all the signals
#         self.signal_continue_snvm.connect(
#             self.continue_snvm_scanning, QtCore.Qt.QueuedConnection
#         )
#         self.signal_stop_scan.connect(self.stop_scanning, QtCore.Qt.QueuedConnection)
#         self.signal_xy_px_acquired.connect(
#             self.move_to_xy_pixel, QtCore.Qt.QueuedConnection
#         )
#         self.signal_freq_px_acquired.connect(
#             self.move_to_freq_pixel, QtCore.Qt.QueuedConnection
#         )
#         self.signal_continue_confocal.connect(
#             self.continue_confocal_scanning, QtCore.Qt.QueuedConnection
#         )
#         self.signal_goto_start.connect(self._go_to_point, QtCore.Qt.QueuedConnection)
#         self.signal_start_optimizer.connect(
#             self.start_optimizer, QtCore.Qt.QueuedConnection
#         )
#         self._optimizer_logic.sigRefocusFinished.connect(self._optimization_complete)
#
#     def on_deactivate(self):
#         """Reverse steps of activation
#         @return int: error code (0:OK, -1:error)
#         """
#         closing_state = ConfocalHistoryEntry(self)
#         closing_state.snapshot(self)
#         self.history.append(closing_state)
#         histindex = 0
#         for state in reversed(self.history):
#             self._statusVariables["history_{0}".format(histindex)] = state.serialize()
#             histindex += 1
#
#         return 0
#
#     def _prepare_data_matrices(self):
#         # Clip the ranges if they are out of bound
#         self.check_xy_ranges()
#
#         # Generate axes and matrices to store the data
#         x_axis = np.linspace(
#             self.scanning_x_range[self._active_stack][0],
#             self.scanning_x_range[self._active_stack][1],
#             self.scanning_x_resolution[self._active_stack],
#         )
#         y_axis = np.linspace(
#             self.scanning_y_range[self._active_stack][0],
#             self.scanning_y_range[self._active_stack][1],
#             self.scanning_y_resolution[self._active_stack],
#         )
#
#         # Now generate the matrices to store the data
#         xy_scan_matrix = np.zeros((len(y_axis), len(x_axis)), dtype=np.float64)
#         # FIXME: for now the stack scanner is the one that's assumed to have the ESR sequence. Maybe consider a flexible
#         #  way of doing this
#         if self._snvm_active:
#             step_number = 1 + round(
#                 (self.stop_freq - self.start_freq) / self.freq_resolution
#             )
#             freq_axis = np.linspace(self.start_freq, self.stop_freq, step_number)
#             snvm_matrix = np.zeros(
#                 (
#                     xy_scan_matrix.shape[0],
#                     xy_scan_matrix.shape[1],
#                     len(freq_axis),
#                     self.odmr_averages,
#                 ),
#                 dtype=np.float64,
#             )
#             temp_freq_matrix = np.full(
#                 (self.odmr_averages, len(freq_axis)), self.invalid
#             )
#             temp_afm_matrix = np.copy(temp_freq_matrix)
#             average_odmr_trace = np.zeros((temp_freq_matrix.shape[1],))
#         else:
#             snvm_matrix = np.copy(xy_scan_matrix[:, :, np.newaxis])
#             freq_axis = None
#             temp_freq_matrix = None
#             average_odmr_trace = None
#             temp_afm_matrix = None
#
#         if self.store_retrace[self._active_stack]:
#             xy_scan_matrix_retrace = np.copy(xy_scan_matrix)
#             snvm_matrix_retrace = np.copy(snvm_matrix)
#         else:
#             xy_scan_matrix_retrace = np.zeros(xy_scan_matrix.shape)
#             snvm_matrix_retrace = np.zeros(snvm_matrix.shape)
#
#         self._tot_xy_scanning_index = 0
#         self._x_scanning_axis = x_axis
#         self._y_scanning_axis = y_axis
#         self.xy_scan_matrix = xy_scan_matrix
#         self.snvm_matrix = snvm_matrix
#         self.temp_freq_matrix = temp_freq_matrix
#         self._temp_afm_matrix = temp_afm_matrix
#         self.average_odmr_trace = average_odmr_trace
#         self.xy_scan_matrix_retrace = xy_scan_matrix_retrace
#         self.snvm_matrix_retrace = snvm_matrix_retrace
#         self.freq_axis = freq_axis
#         self._xy_matrix_width = self.xy_scan_matrix.shape[0]
#         self._freq_axis_length = (
#             len(self.freq_axis) if isinstance(self.freq_axis, np.ndarray) else None
#         )
#         self.completed_pixels_matrix = np.zeros((2,) + xy_scan_matrix.shape, dtype=int)
#
#     def check_xy_ranges(self):
#         """
#         Check that the requested scanning ranges are within the maximum scanning ranges set by the hardware.
#         If they are out of bounds, clip the values.
#         """
#         stk = self._active_stack
#         curr_x_minrange, curr_x_maxrange = self.x_maxrange[stk]
#         curr_y_minrange, curr_y_maxrange = self.y_maxrange[stk]
#         # TODO: emit a signal when the clipping happens, and update the GUI limits accordingly
#         if not (
#             (curr_x_minrange <= self.scanning_x_range[stk][0] <= curr_x_maxrange)
#             or not (curr_x_minrange <= self.scanning_x_range[stk][1] <= curr_x_maxrange)
#         ):
#             self.scanning_x_range[stk] = np.clip(
#                 self.scanning_x_range[stk], curr_x_minrange, curr_x_maxrange
#             )
#             self.log.warning(
#                 "x scanning range limits are out of bounds, clipped back to the maximum values."
#             )
#
#         if not (
#             (curr_y_minrange <= self.scanning_y_range[stk][0] <= curr_y_maxrange)
#             or not (curr_y_minrange <= self.scanning_y_range[stk][1] <= curr_y_maxrange)
#         ):
#             self.scanning_y_range[stk] = np.clip(
#                 self.scanning_y_range[stk], curr_y_minrange, curr_y_maxrange
#             )
#             self.log.warning(
#                 "y scanning range limits are out of bounds, clipped back to the maximum values."
#             )
#         return 0
#
#     def prepare_devices(self):
#         """
#         Initalize the scanning device, prior to starting the scan.
#         """
#
#         self._scanning_device.module_state.lock()
#         self._photon_samples = self.pxtime_to_samples()
#         if self._snvm_active:
#             analog_channels = self._scanning_device.get_ai_counter_channels(
#                 stack_name=self._active_stack
#             )
#         else:
#             analog_channels = None
#         self._scanning_device.prepare_counters(
#             samples_to_acquire=self._photon_samples, counter_ai_channels=analog_channels
#         )
#
#         # FIXME: Also this if-else is for debugging of the slow stepping, remove once done.
#         if not self.slowstepping:
#             self.signal_xy_px_acquired.disconnect()
#             self.signal_xy_px_acquired.connect(
#                 self.move_to_xy_pixel, QtCore.Qt.QueuedConnection
#             )
#             self._scanning_device.create_ao_task(self._active_stack)
#             if self.store_retrace[self._active_stack] is False:
#                 self._scanning_device.prepare_motion_clock()
#                 clk_freq = self._scanning_device.get_motion_clock_frequency()
#                 speed = self._scanning_device.get_motion_speed(self._active_stack)
#                 self.backward_pixels = int(
#                     (
#                         (self._x_scanning_axis.max() - self._x_scanning_axis.min())
#                         / speed
#                     )
#                     * clk_freq
#                 )
#                 if self.backward_pixels < 2:
#                     self.backward_pixels = 2
#
#                 self.x_backward = np.zeros((2, self.backward_pixels))
#                 self.x_backward[0] = np.linspace(
#                     self._x_scanning_axis.max(),
#                     self._x_scanning_axis.min(),
#                     self.backward_pixels,
#                 )
#         else:
#             self.signal_xy_px_acquired.disconnect()
#             self.signal_xy_px_acquired.connect(
#                 self.move_to_xy_pixel_slow, QtCore.Qt.QueuedConnection
#             )
#             self._scanning_device.create_ao_task(self._active_stack)
#             self._scanning_device.prepare_motion_clock()
#             clk_freq = self._scanning_device.get_motion_clock_frequency()
#             speed = self._scanning_device.get_motion_speed(self._active_stack)
#             if self.store_retrace[self._active_stack] is False:
#                 self.backward_pixels = int(
#                     (
#                         (self._x_scanning_axis.max() - self._x_scanning_axis.min())
#                         / speed
#                     )
#                     * clk_freq
#                 )
#                 if self.backward_pixels < 2:
#                     self.backward_pixels = 2
#
#                 self.x_backward = np.zeros((2, self.backward_pixels))
#                 self.x_backward[0] = np.linspace(
#                     self._x_scanning_axis.max(),
#                     self._x_scanning_axis.min(),
#                     self.backward_pixels,
#                 )
#
#             fw_x_step = self._x_scanning_axis[1] - self._x_scanning_axis[0]
#             fw_y_step = self._y_scanning_axis[1] - self._y_scanning_axis[0]
#             self.fw_x_pixels = max(int(clk_freq * fw_x_step / speed), 2)
#             self.fw_y_pixels = max(int(clk_freq * fw_y_step / speed), 2)
#
#             self.x_forward = np.zeros((2, self.fw_x_pixels))
#             self.y_forward = np.zeros((2, self.fw_y_pixels))
#         # FIXME: end
#
#         if self._snvm_active:
#             self._odmrscanner.module_state.lock()
#             self._odmrscanner.set_IQmod(self.podmr_active)
#             try:
#                 pass
#                 # FIXME: look into how to operate the srs in list mode
#                 # self._odmrscanner.set_list(frequency=self.freq_axis, power=self.mw_power)
#             except:
#                 self.log.error(
#                     "Failed loading the frequency axis into the ODMR scanner. Aborted execution."
#                 )
#                 # self._scanning_device.module_state.unlock()
#                 # self._odmrscanner.module_state.unlock()
#
#     def load_pulse_sequence(self):
#         # Convert all the natural units in ns and GHz
#         clkrate = round(self.podmr_clkrate * 1e-9)
#
#         start_delay = self.podmr_start_delay * 1e9
#         laser_init = self.podmr_laser_init * 1e9
#         laser_read = self.podmr_laser_read * 1e9
#         init_delay = self.podmr_init_delay * 1e9
#         read_delay = self.podmr_read_delay * 1e9
#         apd_delay = self.podmr_apd_delay * 1e9
#         apd_read = self.podmr_apd_read * 1e9
#         final_delay = self.podmr_final_delay * 1e9
#
#         pipulse = self.podmr_pipulse * 1e9
#         init_center = start_delay + laser_init / 2
#         mw_center = init_center + laser_init / 2 + init_delay + pipulse / 2
#         read_center = mw_center + pipulse / 2 + read_delay + laser_read / 2
#         apd_center = mw_center + pipulse / 2 + read_delay + apd_delay + apd_read / 2
#
#         tot_time = (
#             start_delay
#             + laser_init
#             + init_delay
#             + pipulse
#             + read_delay
#             + max(laser_read + final_delay, apd_delay + apd_read)
#         )
#
#         tot_samples = self._pulser.waveform_padding(tot_time * clkrate)
#         timeax = np.linspace(0, tot_samples / clkrate, tot_samples)
#         blank_signal = np.zeros_like(timeax)
#         high_signal = blank_signal + 1
#
#         mw_waveform = self._pulser.box_envelope(
#             timeax, np.array([[mw_center, pipulse]])
#         )
#         laser_waveform = self._pulser.box_envelope(
#             timeax, np.array([[init_center, laser_init], [read_center, laser_read]])
#         )
#         apd_waveform = self._pulser.box_envelope(
#             timeax, np.array([[apd_center, apd_read]])
#         )
#
#         card_idx = 1
#         self._pulser.stop_all()
#         self._pulser.set_sample_rate(card_idx, int(self.podmr_clkrate))
#
#         steps = 2
#         do_map = {1: [0, 1, 2]}
#         for step in range(steps):
#             if step == 0:
#                 self._pulser.load_waveform(
#                     iq_dictionary={"i_chan": mw_waveform, "q_chan": blank_signal},
#                     digital_pulses={
#                         "laser": laser_waveform,
#                         "apd_sig": apd_waveform,
#                         "apd_read": blank_signal,
#                     },
#                     digital_output_map=do_map,
#                 )
#             else:
#                 # This waveform is going to be used when the optimizer is called
#                 self._pulser.load_waveform(
#                     iq_dictionary={"i_chan": blank_signal, "q_chan": blank_signal},
#                     digital_pulses={
#                         "laser": high_signal,
#                         "apd_sig": high_signal,
#                         "apd_read": blank_signal,
#                     },
#                     digital_output_map=do_map,
#                 )
#
#         self._pulser.configure_ORmask(card_idx, None)
#         self._pulser.configure_ANDmask(card_idx, None)
#
#         loops = np.ones(steps, dtype=np.int64)
#
#         stop_condition_list = np.array(
#             [
#                 0x40000000,
#             ]
#             * steps,
#             dtype=np.int64,
#         )
#
#         self._pulser.load_sequence(
#             loops_list=loops,
#             stop_condition_list=stop_condition_list,
#         )
#
#         self._pulser.start_card(card_idx)
#         self._pulser.arm_trigger(card_idx)
#
#         self._pulser.send_software_trig(card_idx)
#
#     def start_snvm_scanning(self):
#
#         self.stopRequested = False
#
#         self.module_state.lock()
#
#         self._active_stack = self.sampleStackName
#         self._snvm_active = True
#         self._prepare_data_matrices()
#
#         self._initialize_scanning_statuses()
#         self.prepare_devices()
#
#         # This part needs to stay outside prepare_devices(), otherwise load_pulse_sequence() is going to be called
#         # every time the optimizer is finished.
#         if self.podmr_active:
#             self.load_pulse_sequence()
#         else:
#             # With this line, make sure that the AWG is playing in CW mode if normal ODMR is requested.
#             self._master_pulselogic.cw()
#
#         self.signal_snvm_image_updated.emit()
#         self.signal_snvm_initialized.emit()
#
#         self.go_to_point(
#             [
#                 self._x_scanning_axis[self._x_scanning_index],
#                 self._y_scanning_axis[self._y_scanning_index],
#             ],
#             stack=self._active_stack,
#             clear_ao_whenfinished=False,
#             caller="logic",
#         )
#         # FIXME: look into how to operate the srs in list mode
#         self._odmrscanner.set_frequency(self.freq_axis[self._freq_scanning_index])
#         self._odmrscanner.set_power(self.mw_power)
#         self._odmrscanner.on()
#
#         self.signal_continue_snvm.emit()
#
#     def start_confocal_scanning(self):
#         self.module_state.lock()
#
#         self._active_stack = self.tipStackName
#         self._snvm_active = False
#         self._prepare_data_matrices()
#
#         self._initialize_scanning_statuses()
#         self.prepare_devices()
#
#         self.signal_xy_image_updated.emit()
#         self.signal_confocal_initialized.emit()
#
#         self.go_to_point(
#             [
#                 self._x_scanning_axis[self._x_scanning_index],
#                 self._y_scanning_axis[self._y_scanning_index],
#             ],
#             stack=self._active_stack,
#             clear_ao_whenfinished=False,
#             caller="logic",
#         )
#         self.signal_continue_confocal.emit()
#
#     def start_optimizer(self):
#         if self.podmr_active:
#             self._pulser.send_software_trig(1)
#         self._scanning_device.pause_tasks()
#         self._optimizer_logic.start_refocus()
#
#     def continue_snvm_scanning(self):
#         acquire_data = (
#             False
#             if (self.store_retrace[self._active_stack] is False)
#             and (self._is_retracing is True)
#             else True
#         )
#         if acquire_data:
#             # If the index of the ODMR is less than the averages, keep acquiring
#             if self._odmr_rep_index < self.odmr_averages:
#                 counts, ainput = self.acquire_pixel()
#                 counts = counts[0] / self.px_time[self._active_stack]
#                 self.temp_freq_matrix[
#                     self._odmr_rep_index, self._freq_scanning_index
#                 ] = counts
#                 self._temp_afm_matrix[
#                     self._odmr_rep_index, self._freq_scanning_index
#                 ] = ainput.mean()
#
#                 if self._odmr_rep_index > 0 and self.pxbypx_odmr:
#                     self.average_odmr_trace = np.nanmean(self.temp_freq_matrix, axis=0)
#
#                 if self._is_retracing and self.store_retrace[self._active_stack]:
#                     self.snvm_matrix_retrace[
#                         self._y_scanning_index,
#                         self._x_scanning_index,
#                         self._freq_scanning_index,
#                         self._odmr_rep_index,
#                     ] = counts
#                 else:
#                     self.snvm_matrix[
#                         self._y_scanning_index,
#                         self._x_scanning_index,
#                         self._freq_scanning_index,
#                         self._odmr_rep_index,
#                     ] = counts
#                 self.signal_freq_px_acquired.emit(self._odmr_rep_index)
#                 if self.pxbypx_odmr:
#                     self.signal_odmr_trace_updated.emit(self._odmr_rep_index)
#
#             # Else, the OMDR acquisition for the pixel has finished. Store the data and ask for the next pixel. If also the
#             # Scanning is done, tell that the scanning has finished.
#             else:
#                 if self._is_retracing and self.store_retrace[self._active_stack]:
#                     # FIXME: I am not acquiring and storing properly the analog input, find a way after basic debugging done
#                     self.xy_scan_matrix_retrace[
#                         self._y_scanning_index, self._x_scanning_index
#                     ] = self._temp_afm_matrix.mean()
#                 else:
#                     self.xy_scan_matrix[
#                         self._y_scanning_index, self._x_scanning_index
#                     ] = self._temp_afm_matrix.mean()
#
#                 self.completed_pixels_matrix[
#                     int(self._is_retracing),
#                     self._y_scanning_index,
#                     self._x_scanning_index,
#                 ] = 1
#                 if not self.pxbypx_odmr:
#                     if self._is_retracing and self.store_retrace[self._active_stack]:
#                         thetrace = self.snvm_matrix_retrace[
#                             self._y_scanning_index, self._x_scanning_index
#                         ]
#                     else:
#                         thetrace = self.snvm_matrix[
#                             self._y_scanning_index, self._x_scanning_index
#                         ]
#                     self.last_odmr_trace = thetrace.mean(axis=-1)
#                     self.signal_odmr_line_acquired.emit()
#
#                 # The ODMR sequence has finished. Update the indices accordingly
#                 # self._x_scanning_index += self._x_index_step
#                 self._odmr_rep_index = 0
#
#                 self.temp_freq_matrix[:] = self.invalid
#                 self.average_odmr_trace[:] = self.invalid
#
#                 self.signal_snvm_image_updated.emit()
#                 self.signal_xy_px_acquired.emit()
#
#     def continue_confocal_scanning(self):
#         stk = self._active_stack
#         acquire_data = (
#             False
#             if (self.store_retrace[stk] is False) and (self._is_retracing is True)
#             else True
#         )
#         if acquire_data:
#             counts, _ = self.acquire_pixel()
#             counts = counts[0] / self.px_time[stk]
#             if self._is_retracing and self.store_retrace[stk]:
#                 self.xy_scan_matrix_retrace[
#                     self._y_scanning_index, self._x_scanning_index
#                 ] = counts
#             else:
#                 self.xy_scan_matrix[
#                     self._y_scanning_index, self._x_scanning_index
#                 ] = counts
#             self.completed_pixels_matrix[
#                 int(self._is_retracing), self._y_scanning_index, self._x_scanning_index
#             ] = 1
#         self.signal_xy_image_updated.emit()
#         # self._x_scanning_index += self._x_index_step
#         self.signal_xy_px_acquired.emit()
#
#     def acquire_pixel(self):
#         data = self._scanning_device.read_pixel(self._photon_samples)
#         return data
#
#     def stop_scanning(self):
#         if self.stopRequested:
#             with self.threadlock:
#                 self._scanning_device.close_counters()
#                 self.stopRequested = False
#                 self.stop_xy_scanner()
#                 if self._snvm_active:
#                     self.stop_freq_scanner()
#                     if self.podmr_active:
#                         self._master_pulselogic.cw()
#
#                 # FIXME: this needs to be removed as well when slowstepping is debugged.
#                 if not self.store_retrace[self._active_stack] or self.slowstepping:
#                     self._scanning_device.clear_motion_clock()
#                 self.module_state.unlock()
#
#             self._store_data_matrices(self._snvm_active)
#             self.signal_scan_finished.emit(self._snvm_active)
#             self._snvm_active = False
#
#     def stop_xy_scanner(self):
#         """Closing the scanner device.
#
#         @return int: error code (0:OK, -1:error)
#         """
#         try:
#             self._scanning_device.close_counters()
#         except Exception as e:
#             self.log.exception("Could not close the scanning tasks.")
#         try:
#             self._scanning_device.clear_ao_task(self._active_stack)
#         except:
#             self.log.exception("Could not clear the ao task.")
#         try:
#             if self._scanning_device.module_state() == "locked":
#                 self._scanning_device.module_state.unlock()
#         except Exception as e:
#             self.log.exception("Could not unlock scanning device.")
#
#         return 0
#
#     def stop_freq_scanner(self):
#         self._odmrscanner.off()
#         self._odmrscanner.set_IQmod(self.podmr_active)
#         try:
#             self._odmrscanner.module_state.unlock()
#         except Exception as e:
#             self.log.exception("Could not unlock scanning device.")
#
#         return 0
#
#     def move_to_freq_pixel(self):
#         if not self.stopRequested:
#             self._freq_scanning_index += 1
#             if self._freq_scanning_index == self._freq_axis_length:
#                 self._odmr_rep_index += 1
#                 self._freq_scanning_index = 0
#             self._odmrscanner.set_frequency(self.freq_axis[self._freq_scanning_index])
#             self.signal_continue_snvm.emit()
#         else:
#             self.signal_stop_scan.emit()
#
#     def move_to_xy_pixel(self):
#
#         self._update_indices()
#
#         if not self.stopRequested:
#             if self._is_retracing and not self.store_retrace[self._active_stack]:
#                 self._clocked_retrace()
#
#                 self._x_scanning_index = 0
#                 self._y_scanning_index += 1
#                 self._is_retracing = False
#                 self._x_index_step = 1
#
#             new_x_pos = self._x_scanning_axis[self._x_scanning_index]
#             new_y_pos = self._y_scanning_axis[self._y_scanning_index]
#             self._scanning_device.scanner_set_position(
#                 [new_x_pos, new_y_pos], stack=self._active_stack
#             )
#
#             if self._snvm_active:
#                 # Check if an optimization needs to be done
#                 if (
#                     self.optimize_while_scanning
#                     and self._tot_xy_scanning_index % self.every_N_pixels == 0
#                 ):
#                     self.signal_start_optimizer.emit()
#                 else:
#                     self.signal_continue_snvm.emit()
#             else:
#                 self.signal_continue_confocal.emit()
#         else:
#             self.signal_stop_scan.emit()
#
#     def move_to_xy_pixel_slow(self):
#
#         curr_x = self._x_scanning_axis[self._x_scanning_index]
#         curr_y = self._y_scanning_axis[self._y_scanning_index]
#
#         self._update_indices()
#
#         if not self.stopRequested:
#             if self._is_retracing and not self.store_retrace[self._active_stack]:
#                 self._clocked_retrace()
#
#                 self._x_scanning_index = 0
#                 curr_x = self._x_scanning_axis[0]
#                 self._y_scanning_index += 1
#                 self._is_retracing = False
#                 self._x_index_step = 1
#
#         if not self.stopRequested:
#             new_x_pos = self._x_scanning_axis[self._x_scanning_index]
#             new_y_pos = self._y_scanning_axis[self._y_scanning_index]
#
#             if new_y_pos == curr_y and new_x_pos != curr_x:
#                 # The following is actually faster than (new_x_pos-curr_x)*array_stored_in_memory + curr_x
#                 self.x_forward[0] = np.linspace(curr_x, new_x_pos, self.fw_x_pixels)
#                 self.x_forward[1] = curr_y
#                 self._scanning_device.move_along_line(
#                     position_array=self.x_forward, stack=self._active_stack
#                 )
#             elif new_y_pos != curr_y:
#                 self.y_forward[0] = curr_x
#                 self.y_forward[1] = np.linspace(curr_y, new_y_pos, self.fw_y_pixels)
#                 self._scanning_device.move_along_line(
#                     position_array=self.y_forward, stack=self._active_stack
#                 )
#
#             if self._snvm_active:
#                 # Check if an optimization needs to be done
#                 if (
#                     self.optimize_while_scanning
#                     and self._tot_xy_scanning_index % self.every_N_pixels == 0
#                 ):
#                     self.signal_start_optimizer.emit()
#                 else:
#                     self.signal_continue_snvm.emit()
#             else:
#                 self.signal_continue_confocal.emit()
#         else:
#             self.signal_stop_scan.emit()
#
#     def _update_indices(self):
#         self._x_scanning_index += self._x_index_step
#
#         # This clause is activated when the last x index is reached.
#         if not self._is_retracing and self._x_scanning_index == len(
#             self._x_scanning_axis
#         ):
#             self._is_retracing = True
#             self._x_index_step = -1  # If retracing starts, flip the stepping direction
#             self._x_scanning_index += self._x_index_step
#
#             # If the maximum y index has been reached and the re-trace is not stored, the scan is done.
#             if (
#                 self._y_scanning_index == len(self._y_scanning_axis) - 1
#             ) and not self.store_retrace[self._active_stack]:
#                 self.stopRequested = True
#
#         # When the index becomes negative, it means that the re-trace has stepped back to the leftmost x position.
#         elif self._is_retracing and self._x_scanning_index < 0:
#             self._is_retracing = False
#             self._x_index_step = 1  # The stepping direction now points forward.
#             self._x_scanning_index += self._x_index_step
#             self._y_scanning_index += 1
#
#             # If the maximum y index has been reached, the scan is done.
#             if self._y_scanning_index == len(self._y_scanning_axis):
#                 self.stopRequested = True
#
#         self._tot_xy_scanning_index += 1
#
#
#     def _optimization_complete(self, coords):
#         if self._snvm_active:
#             self.go_to_point(coords, stack=self._optimizer_logic.optimizer_stack)
#
#             if self._odmrscanner.module_state() == "locked":
#                 self._odmrscanner.module_state.unlock()
#
#             if self._scanning_device.module_state() == "locked":
#                 self._scanning_device.module_state.unlock()
#
#             if self.podmr_active:
#                 self._pulser.send_software_trig(1)
#             self.prepare_devices()
#             if not self.stopRequested:
#                 self.signal_continue_snvm.emit()
#             else:
#                 self.signal_stop_scan.emit()
#
#     def pxtime_to_samples(self):
#         return round(
#             self.px_time[self._active_stack]
#             * self._scanning_device.get_counter_clock_frequency()
#         )
#
#     def _clocked_retrace(self):
#
#         self.x_backward[1] = self._y_scanning_axis[self._y_scanning_index]
#
#         self._scanning_device.move_along_line(
#             position_array=self.x_backward, stack=self._active_stack
#         )
#
#     def _store_data_matrices(self, is_snvm_data=True):
#         if is_snvm_data:
#             self._curr_snvm_image = [
#                 self._x_scanning_axis,
#                 self._y_scanning_axis,
#                 self.freq_axis,
#                 np.arange(self.odmr_averages),
#                 self.snvm_matrix,
#                 self.snvm_matrix_retrace,
#             ]
#             self._curr_afm_image = [
#                 self._x_scanning_axis,
#                 self._y_scanning_axis,
#                 self.xy_scan_matrix,
#                 self.xy_scan_matrix_retrace,
#             ]
#         else:
#             self._curr_cfc_image = [
#                 self._x_scanning_axis,
#                 self._y_scanning_axis,
#                 self.xy_scan_matrix,
#                 self.xy_scan_matrix_retrace,
#             ]
#
#     def go_to_point(
#         self, xy_coord, stack=None, clear_ao_whenfinished=True, caller="logic"
#     ):
#         """
#         This public method emits only the signal that calls the private
#         _go_to_point, which calls the scanner_slow_motion of the scanning device.
#         This is done to prevent the GUI from freezing when the scanner_slow_motion
#         is running.
#         """
#         self.signal_goto_start.emit(xy_coord, stack, clear_ao_whenfinished, caller)
#
#     def _go_to_point(self, xy_coord, stack, clear_ao_whenfinished, caller):
#         """
#         Private method that calls the scanner_slow_motion. To be called only by the public
#         go_to_point.
#         """
#         self._scanning_device.scanner_slow_motion(
#             xy_coord, stack=stack, clear_ao_whenfinished=clear_ao_whenfinished
#         )
#
#         self.signal_moved_to_point.emit(caller)
#
#     def get_xy_image_range(self, multiplier=1):
#         """
#         Multiplier is an optional parameter to convert the range to the desired units
#         """
#         return [
#             [
#                 self._x_scanning_axis[0] * multiplier,
#                 self._x_scanning_axis[-1] * multiplier,
#             ],
#             [
#                 self._y_scanning_axis[0] * multiplier,
#                 self._y_scanning_axis[-1] * multiplier,
#             ],
#         ]
#
#     def get_xy_step_size(self, multiplier=1):
#         return [
#             (self._x_scanning_axis[1] - self._x_scanning_axis[0]) * multiplier,
#             (self._y_scanning_axis[1] - self._y_scanning_axis[0]) * multiplier,
#         ]
#
#     def get_stack_names(self):
#         return self._scanning_device.get_stack_names()
#
#     def get_slowmotion_clockrate(self):
#         return self._scanning_device.get_motion_clock_frequency()
#
#     def get_maxranges(self):
#         return self._scanning_device.get_position_range(
#             "sample"
#         ), self._scanning_device.get_position_range("tip")
#
#     def set_slowmotion_clockrate(self, clockrate):
#         self.slow_motion_clock_rate = clockrate
#         self._scanning_device.set_motion_clock_frequency(clockrate)
#
#     def get_motion_speed(self, stack):
#         return self._scanning_device.get_motion_speed(stack)
#
#     def set_motion_speed(self, speed, stack):
#         if stack == "tip":
#             self.backward_speed_conf = speed
#         elif stack == "sample":
#             self.backward_speed_snvm = speed
#         else:
#             raise ValueError(f"Trying to set speed for stack: {stack}")
#         self._scanning_device.set_motion_speed(speed * 1e-6, stack)
#
#     def _initialize_scanning_statuses(self):
#         self._x_scanning_index = 0
#         self._y_scanning_index = 0
#         self._x_index_step = 1
#         self._freq_scanning_index = 0
#         self._odmr_rep_index = 0
#         self._is_retracing = False
#         self.stopRequested = False
#
#     def save_snvm(self):
#         if len(self._curr_snvm_image) > 0:
#             (
#                 xaxsnvm,
#                 yaxsnvm,
#                 freqax,
#                 averages,
#                 snvm_trace,
#                 snvm_retrace,
#             ) = self._curr_snvm_image
#             xaxafm, yaxafm, afm_trace, afm_retrace = self._curr_afm_image
#             data = {
#                 "snvm/x_axis": xaxsnvm,
#                 "snvm/y_axis": yaxsnvm,
#                 "snvm/frequency_axis": freqax,
#                 "snvm/averages": averages,
#                 "snvm/snvm_trace": snvm_trace,
#                 "snvm/snvm_retrace": snvm_retrace,
#                 "afm/x_axis": xaxafm,
#                 "afm/y_axis": yaxafm,
#                 "afm/afm_trace": afm_trace,
#                 "afm/afm_retrace": afm_retrace,
#             }
#             self._savelogic.save_hdf5_data(data)
#         else:
#             self.log.info("The SNVM arrays are empty")
#
#     def save_confocal(self):
#         if len(self._curr_cfc_image) > 0:
#             xax, yax, cfc_trace, cfc_retrace = self._curr_cfc_image
#             data = {
#                 "confocal/x_axis": xax,
#                 "confocal/y_axis": yax,
#                 "confocal/confocal_trace": cfc_trace,
#                 "confocal/confocal_retrace": cfc_retrace,
#             }
#             self._savelogic.save_hdf5_data(data)
#         else:
#             self.log.info("The confocal arrays are empty")