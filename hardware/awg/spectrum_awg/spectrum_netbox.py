# -*- coding: utf-8 -*-

"""
This file contains the Qudi hardware module for the Keysight M3202A PXIe AWG device.
(previously Signadyne SD1).

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

from core.module import Base
from core.configoption import ConfigOption
from collections import defaultdict
from .spcm_tools import *
from .pyspcm import *
from math import log2
from interface.pulser_interface import PulserInterface, PulserConstraints
from qtpy import QtCore


class SpectrumNetbox(Base, PulserInterface):
    _card_ip = ConfigOption("card_ip", missing="error")
    _netbox_type = ConfigOption("netbox_type", "DN2.663-04", missing="info")
    _max_frequency = ConfigOption("max_frequency", 100e6, missing="info")
    _max_ao_voltages = ConfigOption(
        "max_ao_voltage_mV", [500, 500, 500, 500], missing="info"
    )

    _ip_addon = "::inst{}::INSTR"

    # Initalize some lists of parameters that are going to be used in the methods
    __amplitude_chans = [SPC_AMP0, SPC_AMP1, SPC_AMP2, SPC_AMP3]
    __filter_chans = [SPC_FILTER0, SPC_FILTER1, SPC_FILTER2, SPC_FILTER3]
    __stoplevel_chans = [
        SPC_CH0_STOPLEVEL,
        SPC_CH1_STOPLEVEL,
        SPC_CH2_STOPLEVEL,
        SPC_CH3_STOPLEVEL,
    ]
    __output_chan_enable = [
        SPC_ENABLEOUT0,
        SPC_ENABLEOUT1,
        SPC_ENABLEOUT2,
        SPC_ENABLEOUT3,
    ]
    __card_ormask_trigger_dict = {
        None: SPC_TMASK_NONE,
        "immediate": SPC_TMASK_SOFTWARE,  # The cards default, plays immediately after arming
        "ext0": SPC_TMASK_EXT0,
        "ext1": SPC_TMASK_EXT1,
    }

    __card_andmask_trigger_dict = {
        None: SPC_TMASK_NONE,
        "ext0": SPC_TMASK_EXT0,
        "ext1": SPC_TMASK_EXT1,
    }
    __chan_digout_sources = [
        SPCM_XMODE_DIGOUTSRC_CH0,
        SPCM_XMODE_DIGOUTSRC_CH1,
        SPCM_XMODE_DIGOUTSRC_CH2,
    ]
    __bit_digout_sources = [
        SPCM_XMODE_DIGOUTSRC_BIT15,
        SPCM_XMODE_DIGOUTSRC_BIT14,
        SPCM_XMODE_DIGOUTSRC_BIT13,
    ]
    __digout_channels = [SPCM_X0_MODE, SPCM_X1_MODE, SPCM_X2_MODE]

    sigStepLoaded = QtCore.Signal(int, int)

    def on_activate(self):
        # TODO: maybe implement a network discovery. If so, it should be fast.
        if self._netbox_type == "DN2.663-04":
            self._card_num = 2
            self._ao_chans_percard = 2
            self._min_clock_rate = MEGA(50)
        else:
            self.log.error("Unknown digitizerNETBOX/generatorNETBOX.")
            return -1

        self._error_msg = create_string_buffer(ERRORTEXTLEN)

        self._netbox = CardCollection()

        for ii in range(self._card_num):
            card_address = bytes(self._card_ip + self._ip_addon.format(ii), "utf-8")
            card = spcm_hOpen(create_string_buffer(card_address))
            if card is None:
                self.log.error("Could not load card.")
                return -1
            # Make sure that the card is an AO card.
            card_function = self._command_get(card, SPC_FNCTYPE)
            cardfeatures = self._command_get(card, SPC_PCIFEATURES)
            if (card_function != SPCM_TYPE_AO) or (
                cardfeatures & SPCM_FEAT_SEQUENCE
            ) == 0:
                self.log.error(
                    "Something is wrong, you are supposed to have an AO card with sequence mode."
                )
                return -1

            # The card with the starhub is also the one that controls the synchronization, do and triggers.
            has_starhub = self._command_get(card, SPCM_FW_MODEXTRA)
            if has_starhub:
                is_master = True
            else:
                is_master = False
            chan_num = self._command_get(card, SPC_MIINST_CHPERMODULE)
            self._netbox.add_card(card, is_master=is_master)
            self._netbox.set_chan_num(chan_num)
        self.log.info("Loaded all the cards successfully.")
        if self._netbox.master() is None:
            self.log.warning("No card with starhub found.")
        self._netbox.get_info()

        # Now open the starhub communication
        self._netbox.starHub = spcm_hOpen(create_string_buffer(b"sync0"))

        if self._netbox.netbox_type != self._netbox_type:
            self.log.error(
                f"The expected netbox type is {self._netbox_type}, "
                f"but the device returned {self._netbox.netbox_type}. Aborting."
            )
            self.on_deactivate()

        self._waveform_container = []
        self._sequence_container = dict()

    def on_deactivate(self):
        for ii in range(self._netbox.card_amount()):
            spcm_vClose(self._netbox.card(ii))
            self.log.info(f"Successfully close card {ii} of netbox {self._netbox_type}")
        if self._netbox.starHub is not None:
            spcm_vClose(self._netbox.starHub)
            self.log.info("Successfully closed starhub.")

    def set_chan_amplitude(self, channels, amplitudes):
        """
        Set the channel max output amplitude (i.e., set the gain of the output amplifier).
        @param list channels: the list of the channels, numbered from 0 to 3
        @param list amplitudes: the list of the amplitudes in mV required for each channel.
        @return int: Error code. 0: ok, -1: error occurred
        """
        # Make sure everything is sorted
        channels, amplitudes = (
            list(srlist) for srlist in zip(*sorted(zip(channels, amplitudes)))
        )
        for channel, amplitude, max_ao_voltage in zip(
            channels, amplitudes, self._max_ao_voltages
        ):
            card_idx = channel // self._card_num
            if amplitude > max_ao_voltage:
                amplitude = max_ao_voltage
                self.log.warning(
                    f"Requested an ao voltage exceeding the {max_ao_voltage} bound."
                    "Clipping applied."
                )
            command = self.__amplitude_chans[channel % self._card_num]
            spcm_dwSetParam_i64(self._netbox.card(card_idx), command, int64(amplitude))
            errorout = self._get_error_msg(self._netbox.card(card_idx))
            if errorout:
                return -1
        return 0

    def set_output_filters(self, channels, filter_active):
        """
        Set the output filters.

        @param list channels: the list of channels, numbered from 0 to 3
        @param list filter_active: 0 or 1, specifies if the filter needs to be activated or not.
        @return int: Error code. 0: ok, -1: error occurred
        """
        # Make sure everything is sorted
        channels, filter_active = (
            list(srlist) for srlist in zip(*sorted(zip(channels, filter_active)))
        )
        for channel, value in zip(channels, filter_active):
            card_idx = channel // self._card_num
            spcm_dwSetParam_i64(
                self._netbox.card(card_idx),
                self.__filter_chans[channel % self._card_num],
                int64(value),
            )
            errorout = self._get_error_msg(self._netbox.card(card_idx))
            if errorout:
                return -1
        return 0

    def set_sample_rate(self, card_idx, clk_rate=MEGA(50)):
        """
        Set the sampling rate for one of the internal cards of the awg.

        @param int card_idx: The card index.
        @param int clk_rate: The clock rate.
        @return int: Error code. 0: ok, -1: error occurred
        """
        card = self._netbox.card(card_idx)
        if not isinstance(clk_rate, int):
            clk_rate = int(clk_rate)
            self.log.warning(
                f"Given clock rate {clk_rate} is not an integer. Rounding will occurr."
            )
        if clk_rate < self._min_clock_rate:
            self.log.error(
                f"Minimum clock rate is {self._min_clock_rate}, you provided {clk_rate}."
            )
            return -1

        spcm_dwSetParam_i64(card, SPC_SAMPLERATE, clk_rate)
        errorout = self._get_error_msg(card)
        if errorout:
            return -1
        else:
            return 0

    def get_sample_rate(self, card_idx):
        """
        Get the sampling rate for one of the internal cards of the awg.
        @param card_idx: The card index.
        """
        card = self._netbox.card(card_idx)
        return self._command_get(card, SPC_SAMPLERATE)

    def set_active_channels(self, *channels):
        """
        Activates or deactivates disables the physical output channels. *channels represents the list of requested channels,
        indexed from 0 to n (to stay consistent with python indexing).
        If the physical channels of each card in the netbox are less than the required channels, but if there are enough
        outputs in the rest of the netbox, the analog output is adjusted accordingly.

        @params int channels: unpacked list of channel integers, from 0 to n.
        @return int: Error code. 0: ok, -1: error occurred
        """

        # Assuming here that all cards have the same number of channels.
        all_channels = range(self._netbox.get_chan_num(0) * self._netbox.card_amount())
        for channel in all_channels:
            card_idx = channel // self._card_num
            ch_idx = channel % self._card_num
            command = self.__output_chan_enable[ch_idx]
            value = 1 if channel in channels else 0
            spcm_dwSetParam_i64(self._netbox.card(card_idx), command, value)
        return 0

    def get_active_channels(self, ch=None):
        """
        Check if a list of channels is enabled.
        @param list ch: The list of channels to be checked. Default is None (all the channels will be checked).
        @return list: The list indicating whether the channels are active (element=1) or not (element=0).
        """
        if ch is None:
            ch = range(self._netbox.get_chan_num(0) * self._netbox.card_amount())

        channels_active = [
            0,
        ] * len(ch)
        for ii, channel in enumerate(ch):
            card_idx = channel // self._card_num
            ch_idx = channel % self._card_num
            command = self.__output_chan_enable[ch_idx]
            channels_active[ii] = self._command_get(
                self._netbox.card(card_idx), command
            )
        return channels_active

    def load_sequence(
        self,
        waveform_list=None,
        segment_map=np.array([]),
        loops_list=np.array([]),
        stop_condition_list=np.array([]),
    ):
        """
        Method that loads the waveforms into the AWG. Assumes that all the different sequences will be played
        synchronously, meaning that all the requested sequences need to have the same length.
        @param list waveform_list: A list of AbstractWaveform. If None, the sequences will
                                   be loaded by the array of abstract waveforms created using the load_waveform method.
        @param np.ndarray segment_map: A list containing the order of the memory segments associated with the sequence
                                       steps. E.g. for three sequences, segment_map = [0,2,1] will associate step 0 to
                                       segment 0, step 1 to segment 2, step 2 to segment 1. Default is sequential.
        @param np.ndarray loops_list: A list containing the number of times each sequence step will be repeated.
        @param np.ndarray stop_condition_list: A list setting the stop condition for each sequence step. The default
                                               waits for each step to finish the loops, until the last one which stops
                                               the whole sequence.
        @param dict digital_output_map: optional dictionary that maps the source of the waveforms (i.e. the ao channels)
                                        to the physical do channels. I.e., digital_output_map = {0: [2,0]}, assigns
                                        the two waveforms that will be contained in ch0 waveform to the physical do x2
                                        and x0 channels, in this order.
        """

        if waveform_list is None:
            waveform_list = self._waveform_container

        # Here it is assumed that all the sequence steps use the same cards
        first_wform = waveform_list[0]
        cards = first_wform.required_cards()
        required_channels = first_wform.required_channels()

        self.set_active_channels(*required_channels)
        self._chan_enable(*required_channels)

        steps_num = len(waveform_list)

        if not len(segment_map) > 0:
            segment_map = np.arange(0, steps_num, dtype=np.int64)
        if not len(loops_list) > 0:
            loops_list = np.ones((steps_num,), dtype=np.int64)
        if not len(stop_condition_list) > 0:
            stop_condition_list = np.zeros((steps_num,), dtype=np.int64)
            stop_condition_list[:-1] = SPCSEQ_ENDLOOPALWAYS
            stop_condition_list[-1] = SPCSEQ_END

        # The segments in the memory need to be a power of 2.
        max_segments = 1 << (int(log2(steps_num - 1)) + 1)
        if max_segments > self._netbox.maxsegments:
            self.log.error(
                "The requested steps exceed the maximum allowed segmentation."
            )
            return -1

        # Get the memory requirements for each segment, then set up the sequence settings
        # TODO: for now, I'm not using the segment_memory_requirements to do the error checking on the length of each
        #  segment. Implement it in the near future.
        segment_memory_requirements = dict()
        for card_idx in cards:
            current_card = self._netbox.card(card_idx)
            active_chans = self._count_active_channels(current_card)
            min_segment_memory = 384 // active_chans
            max_segment_memory = int(self._netbox.memsize / active_chans / max_segments)

            segment_memory_requirements[card_idx] = [
                min_segment_memory,
                max_segment_memory,
            ]

            # Set up the sequence settings
            spcm_dwSetParam_i64(current_card, SPC_CARDMODE, SPC_REP_STD_SEQUENCE)
            spcm_dwSetParam_i64(
                current_card, SPC_SEQMODE_MAXSEGMENTS, int64(max_segments)
            )
            spcm_dwSetParam_i64(current_card, SPC_SEQMODE_STARTSTEP, 0)

        self.sigStepLoaded.emit(0, steps_num)
        for step in range(steps_num):
            current_waveform = waveform_list[step]
            (
                current_ao_waveform,
                current_do_waveform,
                digital_output_map,
            ) = current_waveform.generate_awg_waveform()

            self._assign_digital_output_channels(
                digital_waveforms=current_do_waveform,
                digital_output_map=digital_output_map,
            )

            self._create_waveform_buffers(
                analog_waveforms=current_ao_waveform,
                digital_waveforms=current_do_waveform,
                sequence_step=step,
            )
            self.sigStepLoaded.emit(step, steps_num)
        # Now combine all the segment settings into a value, which is then written to the card
        step_array = np.arange(
            0, steps_num, dtype=np.int64
        )  # By default, the steps in memory are loaded from 0 to step_num

        # The steps map should always start with the first step to be played. With this, it is assumed that if no stop
        # event occurs, the sequence will start replaying from the first step after going through the rest of the steps
        next_steps = np.roll(step_array, -1)
        # Now combine all the step, next step, loops and stop conditions into a single value
        values_for_stepmem = np.stack(
            (segment_map, next_steps << 16, loops_list << 32, stop_condition_list << 32)
        )
        values_for_stepmem = np.bitwise_or.reduce(
            values_for_stepmem, axis=0, dtype=np.int64
        )

        # Write the value to the card
        for step, value in zip(step_array, values_for_stepmem):
            for card_idx in cards:
                current_card = self._netbox.card(card_idx)
                spcm_dwSetParam_i64(
                    current_card, SPC_SEQMODE_STEPMEM0 + int(step), int64(value)
                )
                errorout = self._get_error_msg(current_card)
                if errorout:
                    return -1

        for card_idx in cards:
            current_card = self._netbox.card(card_idx)

            # This step is fundamental. Without it, the card won't be able to apply the settings properly.
            spcm_dwSetParam_i64(current_card, SPC_M2CMD, M2CMD_CARD_WRITESETUP)

    def play_waveform(self, waveform=None, loops=0):
        """
        Writes a waveform to the awg memory, without breaking the latter into chunks as it's done for sequences.

        @param AbstractWaveform waveform: An AbstractWaveform.
        @param int loops: Number of loops for the waveform replay. Default is 0 (i.e. infinite replay).
        @return int: Error code. 0: ok, -1: error occurred
        """
        if waveform is None:
            self.log.error("Please provide a waveform as input to this method.")
            return -1
        cards = waveform.required_cards()
        required_channels = waveform.required_channels()

        self.set_active_channels(*required_channels)
        self._chan_enable(*required_channels)

        for card_idx in cards:
            current_card = self._netbox.card(card_idx)
            spcm_dwSetParam_i64(current_card, SPC_LOOPS, loops)
            spcm_dwSetParam_i64(current_card, SPC_CARDMODE, SPC_REP_STD_SINGLE)

            (
                ao_waveform,
                do_waveform,
                digital_output_map,
            ) = waveform.generate_awg_waveform()
            self._assign_digital_output_channels(
                digital_waveforms=do_waveform, digital_output_map=digital_output_map
            )

            self._create_waveform_buffers(
                analog_waveforms=ao_waveform, digital_waveforms=do_waveform
            )

        for card_idx in cards:
            current_card = self._netbox.card(card_idx)

            # This step is fundamental. Without it, the card won't be able to apply the settings properly.
            spcm_dwSetParam_i64(current_card, SPC_M2CMD, M2CMD_CARD_WRITESETUP)

        return 0

    def configure_ORmask(self, card_idx, *masks_to_enable):
        """
        Configure the OR mask of an internal awg card.
        @param card_idx: The index of the card
        @param masks_to_enable: Unpacked list of non-keyworded arguments. The values should correspond to the values
                                that can be passed to the SPC_TRIG_ORMASK register (check the Spectrum Netbox docs).
        @return int: Error code. 0: ok, -1: error occurred
        """
        final_mask = 0
        for mask in masks_to_enable:
            if mask not in self.__card_ormask_trigger_dict:
                self.log.error("Requested register is not valid. Operation aborted.")
                return -1
            final_mask |= self.__card_ormask_trigger_dict[mask]

        spcm_dwSetParam_i64(self._netbox.card(card_idx), SPC_TRIG_ORMASK, final_mask)

        return 0

    def configure_ANDmask(self, card_idx, *masks_to_enable):
        """
        Configure the AND mask of an internal awg card.
        @param card_idx: The index of the card
        @param masks_to_enable: Unpacked list of non-keyworded arguments. The values should correspond to the values
                                that can be passed to the SPC_TRIG_ANDMASK register (check the Spectrum Netbox docs).
        @return int: Error code. 0: ok, -1: error occurred
        """
        final_mask = 0
        for mask in masks_to_enable:
            if mask not in self.__card_andmask_trigger_dict:
                self.log.error("Requested register is not valid. Operation aborted.")
                return -1
            final_mask |= self.__card_ormask_trigger_dict[mask]

        spcm_dwSetParam_i64(self._netbox.card(card_idx), SPC_TRIG_ANDMASK, final_mask)

    def send_software_trig(self, card_idx):
        """
        Send a software trigger to one of the internal cards of the awg.

        @param int card_idx: The card index.
        """
        spcm_dwSetParam_i64(
            self._netbox.card(card_idx), SPC_M2CMD, M2CMD_CARD_FORCETRIGGER
        )

    def start_card(self, card_idx):
        """
        Starts one of the internal cards of the awg. NB: this method starts the card, but
        the arm_trigger method needs to be called to actually start playing the loaded waveforms.

        @param int card_idx:
        """
        spcm_dwSetParam_i64(self._netbox.card(card_idx), SPC_M2CMD, M2CMD_CARD_START)

    def arm_trigger(self, card_idx):
        """
        Arms one of the internal cards of the awg to start listening for digital or software triggers.

        @param int card_idx: The index of the card.
        """
        spcm_dwSetParam_i64(
            self._netbox.card(card_idx), SPC_M2CMD, M2CMD_CARD_ENABLETRIGGER
        )

    def disable_trigger(self, card_idx):
        """
        Disables the trigger listening for one of the internal cards of the awg.

        @param int card_idx: The index of the card.
        """
        spcm_dwSetParam_i64(
            self._netbox.card(card_idx), SPC_M2CMD, M2CMD_CARD_DISABLETRIGGER
        )

    def configuration_sequence_mode(self, card_idx):
        """
        Configure an internal awg card for sequence mode.

        @param int card_idx: The index of the card.
        """
        spcm_dwSetParam_i64(
            self._netbox.card(card_idx), SPC_CARDMODE, SPC_REP_STD_SEQUENCE
        )

    def configuration_single_mode(self, card_idx):
        """
        Configure an internal awg card for single waveform replay mode.

        @param int card_idx: The index of the card.
        """
        spcm_dwSetParam_i64(
            self._netbox.card(card_idx), SPC_CARDMODE, SPC_REP_STD_SINGLE
        )

    def stop_replay(self, card_idx):
        """
        Stops the replay for an internal awg card.

        @param int card_idx:
        """
        spcm_dwSetParam_i64(self._netbox.card(card_idx), SPC_M2CMD, M2CMD_CARD_STOP)

    def reset(self):
        """
        Resets all the awg internal cards.
        """
        for card_idx in range(self._netbox.card_amount()):
            self._card_reset(self._netbox.card(card_idx))
        self._waveform_container = []

    def get_interleave(self):
        """
        Not used
        """
        return False

    def set_interleave(self, state=False):
        """
        Not used
        """
        pass

    @staticmethod
    def waveform_padding(waveform_len):
        """
        Takes a desired waveform length, and calculates the correct value to be used to generate the waveform arrays.
        This is done to ensure that the the waveform length is a multiple of 32, as required by the awg.

        @param int waveform_len: The length
        """
        if waveform_len % 32 > 0:
            outlen = 32 * (waveform_len // 32) + 32
        else:
            outlen = waveform_len
        return int(outlen)

    def get_constraints(self):
        """
        Gets the awg constraints.
        @return:
        """
        return self._netbox.constraints

    def pulser_on(self):
        """
        Activates all the AWG channels.
        """
        all_chans = self._netbox.get_chan_num(0) * self._netbox.card_amount()
        self.set_active_channels([ii for ii in range(all_chans)])

    def pulser_off(self):
        """
        Disasbles all the awg channels.
        """
        self.set_active_channels([])

    def load_waveform(
        self,
        ao_waveform_dictionary=dict(),
        do_waveform_dictionary=dict(),
        digital_output_map=dict(),
    ):
        """
        Creates an AbstractWaveform and appends it to a private array. The private array is going to be used by the
        load_waveform method (check the method docstring for details).

        @param dict ao_waveform_dictionary: A dictionary containing the analog waveforms. Each key is an integer which
                                            indicates the  awg channel (i.e. 0=ch0, 1=ch1, ...).
        @param dict do_waveform_dictionary: A dictionary containing the digital waveforms. Each key is an integer which
                                            indicates the digital output channel (i.e. 0=x0, 1=x1, ...).
        @param dict digital_output_map: Dictionary that maps the source associated to the digital waveforms (i.e. the
                                        ao channels) to the physical do channels. I.e., digital_output_map = {0: [2,0]},
                                         assigns the two waveforms that will be contained in ch0 waveform to the physical
                                         do x2 and x0 channels, in this order.
        """
        waveform = AbstractWaveform(self._netbox.card_amount(), self._netbox.masteridx)
        waveform.assign_ao_waveforms(ao_waveform_dictionary)
        waveform.assign_do_waveforms(do_waveform_dictionary, digital_output_map)
        self._waveform_container.append(waveform)

    def get_waveform(self, waveform_idx=-1):
        """
        Return a waveform loaded in the waveform buffer, for a given index. Default is -1 (last waveform).
        """

        return self._waveform_container[waveform_idx]

    def get_loaded_assets(self):
        """
        Not used
        """
        return {}, ""

    def clear_all(self):
        """
        Not used
        """
        self._waveform_container = []
        return 0

    def get_status(self):
        """
        Not used
        """
        return 0, {}

    def get_analog_level(self):
        """
        Not used
        """
        return {}, {}

    def set_analog_level(self):
        """
        Not used
        """
        pass

    def get_digital_level(self):
        """
        Not used
        """
        return {}, {}

    def set_digital_level(self):
        """
        Not used
        """
        pass

    def write_waveform(self):
        """
        Not used
        """
        self.log.error("Write waveform not implemented.")

    def write_sequence(self):
        """
        Not used
        """
        self.log.error("Write sequence not implemented.")

    def get_waveform_names(self):
        """
        Not used
        """
        return []

    def get_sequence_names(self):
        """
        Not used
        """
        return []

    def delete_waveform(self):
        """
        Not used
        """
        self.log.error("delete_waveform not implemented.")
        return []

    def delete_sequence(self):
        """
        Not used
        """
        self.log.error("delete_sequence not implemented.")
        return []

    def get_interleave(self):
        """
        Not used
        """
        self.log.error("Interleave not available on this device.")
        return False

    def set_interleave(self):
        """
        Not used
        """
        self.log.error("Interleave not available on this device.")

    def _assign_digital_output_channels(
        self, digital_waveforms=dict(), digital_output_map=dict()
    ):
        """
        Function that assigns the digital waveforms to their correct AO source in the master card, and also assigns
        them to the physical output channels.
        @param dict digital_waveforms: the dictionary containing the do waveforms, ordered as in _create_waveform_buffers
        @param dict digital_output_map: optional dictionary that maps the source of the waveforms (i.e. the ao channels)
                                        to the physical do channels. I.e., digital_output_map = {0: [2,0]}, assigns
                                        the two waveforms that will be contained in ch0 waveform to the physical do x2
                                        and x0 channels, in this order.
        """
        if digital_waveforms:
            mastercard = self._netbox.master()
            if not digital_output_map:
                digital_output_map = dict().fromkeys(digital_waveforms.keys())
                last_used_chan = 0
                for chan in digital_output_map:
                    digital_output_map[chan] = [
                        last_used_chan + ii
                        for ii in range(len(digital_waveforms[chan]))
                    ]
                    last_used_chan += digital_output_map[chan][-1] + 1

            for channel, xoutchannels in digital_output_map.items():
                for ii, xoutput in enumerate(xoutchannels):
                    x_chan_mode = (
                        SPCM_XMODE_DIGOUT
                        | self.__chan_digout_sources[channel]
                        | self.__bit_digout_sources[ii]
                    )
                    spcm_dwSetParam_i64(
                        mastercard, self.__digout_channels[xoutput], x_chan_mode
                    )

    def _create_waveform_buffers(
        self, analog_waveforms=dict(), digital_waveforms=dict(), sequence_step=None
    ):
        """
        Given a bunch of waveforms, it prepares the buffers for data transfer and loads them into the memory.
        Channels have to be activated somewhere else.

        @param dict analog_waveforms: A dictionary of subdictionaries. The keys of the dictionary need to be
                                      integer numbers, and they describe which is the required card. Each
                                      subdictionary is indexed with the channels indices (i.e. 0,1,2...)
        @param dict digital_waveforms:  A dictionary containing the digital waveforms. The keys are 0,1,2,...
                                        corresponding to which output channel the digital waveforms will be assigned
                                        to. Each element of the dictionary is a numpy array with shape MxN, where M
                                        corresponds to number of digital waveforms assigned to the channel. If you have
                                        a device with 16 bits, they will first be assigned to the 15th bit, then the
                                        14th, and so on.
        @param int sequence_step: Optional parameter to be fed when loading sequences. Defines the segment nummber.
        @param bool rescale_maxbit_all: Optional. If True, check which of the ao channels has the least available bits,
                                        and scales the amplitude of the remaining channels accordingly. I.e., if
                                        ch0 has 14-bit resolution (i.e. two digital channels are assigned to ch0),
                                        rescale_max_bit_all ensures that the waveforms assigned to the other channels
                                        are rescaled to a maximum value of 14 bits.
        """

        masteridx = self._netbox.masteridx
        maxADC = self._netbox.adc_resolution

        # First calculate how many cards are needed.
        cards_required = self._get_required_cards(analog_waveforms, digital_waveforms)

        # Now start preparing the buffers
        for card_idx in cards_required:
            card_ao_waveforms = None
            req_ao_channels = []
            req_do_channels = []
            waveform_ao_size = None
            waveform_do_size = None

            # First some check on the ao and do waveforms
            if analog_waveforms and card_idx in analog_waveforms:
                card_ao_waveforms = analog_waveforms[card_idx]
                req_ao_channels = list(card_ao_waveforms.keys())

                # Check the lengths of the analog waveforms to be loaded, they need to be the same.
                ao_sizes = list(map(len, card_ao_waveforms.values()))
                if len(set(ao_sizes)) > 1:
                    self.log.error(
                        "The analog waveforms loaded on the same card need to have equal length."
                    )
                    return -1
                else:
                    waveform_ao_size = ao_sizes[0]
            if card_idx == masteridx and digital_waveforms:
                req_do_channels = list(digital_waveforms.keys())

                # Check the lengths of the digital waveforms to be loaded, they need to be the same.
                do_sizes = [wform.shape[1] for wform in digital_waveforms.values()]
                if len(set(do_sizes)) > 1:
                    self.log.error(
                        "The digital waveforms loaded on the same card need to have equal length."
                    )
                    return -1
                else:
                    waveform_do_size = do_sizes[0]
            # Check that both waveforms sizes match for the same card.
            if (
                waveform_ao_size is not None
                and waveform_do_size is not None
                and (waveform_ao_size != waveform_do_size)
            ):
                self.log.error(
                    "The length of the analog waveforms must match the length of the digital waveforms on "
                    "the same card."
                )
                return -1
            waveform_size = (
                waveform_ao_size if waveform_ao_size is not None else waveform_do_size
            )

            tot_channels = sorted(set(req_ao_channels + req_do_channels))

            # Now combine all the waveforms and channels for the same card into one array.
            channel_matrix = np.zeros((waveform_size, len(tot_channels)), dtype=c_int16)

            for chan_idx, channel in enumerate(tot_channels):
                # By default, use the full resolution of the instrument
                maxADC_chan = maxADC - 1

                if (
                    card_idx == masteridx
                    and digital_waveforms
                    and channel in digital_waveforms
                ):

                    # If there are digital channels, you need to sacrifice one bit for each digital channel
                    # (unless rescale_maxbit_all was requested)
                    maxADC_chan = (maxADC // 2 ** len(digital_waveforms[channel])) - 1

                    channel_do_matrix = (
                        digital_waveforms[channel].T
                        * 2 ** (15 - np.arange(len(digital_waveforms[channel])))
                    ).astype(c_int16)
                    # The line above shifts puts each do signal in its right bit
                    channel_matrix[:, chan_idx] = np.sum(channel_do_matrix, axis=1)
                if card_ao_waveforms and channel in card_ao_waveforms:
                    # Apply rescaling and adjust depending on whether the waveforms have 16, 15, or 14 usable bits
                    rescaled_ao = self._waveform_16_to_nbits(
                        card_ao_waveforms[channel], maxADC_chan
                    )
                    channel_matrix[:, chan_idx] = (
                        channel_matrix[:, chan_idx] + rescaled_ao
                    )

            # Reshape the matrix to interleave the sample arrays
            channel_matrix = channel_matrix.reshape(np.prod(channel_matrix.shape))

            buffer_size = self._netbox.bytes_persample * len(channel_matrix)
            pvData = pvAllocMemPageAligned(
                buffer_size
            )  # The page-aligned buffer that will be transferred
            # FIXME: here I'm casting to a 16 bit pointer. This is not general for any netbox.
            pnData = cast(
                pvData, ptr16
            )  # The pointer to the buffer, used to fill it up

            # Here a nice ctypes/C trick to transfer the whole array to the pnData pointer, without using for loops
            memmove(
                pnData,
                channel_matrix.ctypes.data,
                sizeof(c_int16) * len(channel_matrix),
            )

            if sequence_step is not None:
                spcm_dwSetParam_i64(
                    self._netbox.card(card_idx),
                    SPC_SEQMODE_WRITESEGMENT,
                    int64(sequence_step),
                )
                spcm_dwSetParam_i64(
                    self._netbox.card(card_idx),
                    SPC_SEQMODE_SEGMENTSIZE,
                    int64(waveform_size),
                )

            spcm_dwDefTransfer_i64(
                self._netbox.card(card_idx),
                SPCM_BUF_DATA,
                SPCM_DIR_PCTOCARD,
                0,
                pvData,
                0,
                int64(buffer_size),
            )
            spcm_dwSetParam_i64(
                self._netbox.card(card_idx),
                SPC_M2CMD,
                M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA,
            )

    def _get_required_cards(self, analog_waveforms=dict(), digital_waveforms=dict()):
        """
        Get which cards are needed to play the required analog_waveforms and digital_waveforms.

        @param analog_waveforms: A dictionary of analog_waveforms, as passed to _create_waveforms_buffers.
        @param digital_waveforms: A dictionary of digital_waveforms, as passed to _create_waveforms_buffers.
        @return list: The list of required cards. Returns -1 if an error occurred.
        """
        masteridx = self._netbox.masteridx

        cards_required = []
        if analog_waveforms:
            cards_required += list(analog_waveforms.keys())
        if digital_waveforms:
            cards_required.append(
                masteridx
            )  # The card used for digital output is always the master card.
            for chan, do_waveform in digital_waveforms.items():
                if len(do_waveform.shape) != 2:
                    self.log.error(
                        f"Failed at channel {chan}. The digital output waveform matrices need to be "
                        f"two-dimensional. "
                    )
                    return -1

        return sorted(set(cards_required))

    def _get_required_channels(self, analog_waveforms=dict(), digital_waveforms=dict()):
        """
        Get which channels need to be activated to play the desired waveforms.

        @param analog_waveforms: A dictionary of analog_waveforms, as passed to _create_waveforms_buffers.
        @param digital_waveforms: A dictionary of digital_waveforms, as passed to _create_waveforms_buffers.
        @return list: The list of required cards.
        """
        masteridx = self._netbox.masteridx

        # Assuming same number of channels per card
        chns_per_card = self._netbox.get_chan_num(0)

        channels_required = []
        if analog_waveforms:
            for card, card_waveforms in analog_waveforms.items():
                channels_required += [
                    chn + card * chns_per_card for chn in card_waveforms.keys()
                ]
        if digital_waveforms:
            channels_required += [
                chn + masteridx * chns_per_card for chn in digital_waveforms.keys()
            ]

        return sorted(set(channels_required))

    def _waveform_16_to_nbits(self, waveform, max_value):
        out_waveform = (max_value * waveform).astype(c_int16)
        if max_value < (2**15 - 1):
            sign_waveform = np.zeros(out_waveform.shape, dtype=c_int16)
            sign_waveform[out_waveform < 0] = 1
            out_waveform = (out_waveform & max_value) + (
                sign_waveform * (max_value + 1)
            )
        return out_waveform

    def _chan_enable(self, *channels):
        """
        Enables or disables channels. This affects the maximum available memory that each channel has for the waveforms.
        *channels represents the list of requested channels, indexed from 0 to n (to stay consistent with python indexing).
        If the physical channels of each card in the netbox are less than the required channels, but if there are enough
        outputs in the rest of the netbox, the analog output is adjusted accordingly.

        @params int channels: unpacked list of channel integers, from 0 to n.
        @return int: -1 if errors occured, 0 otherwise
        """
        channels = sorted(channels)
        chans_to_activate = defaultdict(int)
        for channel in channels:
            card_idx = channel // self._card_num
            chans_to_activate[card_idx] += 2 ** (channel % self._card_num)

        for card_idx, chan_value in chans_to_activate.items():
            if chan_value != 0:
                spcm_dwSetParam_i64(
                    self._netbox.card(card_idx), SPC_CHENABLE, int64(chan_value)
                )
                errorout = self._get_error_msg(self._netbox.card(card_idx))
                if errorout:
                    return -1
        return 0

    def _sequence_max_memory(self, card):
        """
        Gets the maximum memory for the sequence, which depends on the number of active channels.

        @param card: The card onto which the sequence will be loaded.
        """
        active_chans = self._count_active_channels(card)
        return self._netbox.memsize // active_chans

    def _set_stop_level(self, channels, stoplevels):
        """
        Set the stop voltage level at the end of replay.
        @param list channels: the list of channels, numbered from 0 to 3
        @param list stoplevels: the list of available stoplevels, as listed in the documentation.
        @return int: -1 if errors occured, 0 otherwise
        """
        # Make sure everything is sorted
        channels, stoplevels = (
            list(srlist) for srlist in zip(*sorted(zip(channels, stoplevels)))
        )

        for channel, stoplevel in zip(channels, stoplevels):
            card_idx = channel // self._card_num
            spcm_dwSetParam_i64(
                self._netbox.card(card_idx),
                self.__stoplevel_chans[channel % self._card_num],
                int64(stoplevel),
            )

            errorout = self._get_error_msg(self._netbox.card(card_idx))
            if errorout:
                return -1
        return 0

    def _count_active_channels(self, card):
        """
        Return the active channel count for a card.
        @param int card: The card index.
        @return int: the amount of active channels.
        """
        return self._command_get(card, SPC_CHCOUNT)

    @staticmethod
    def _command_get(card, command):
        """
        Function to stop keeping writing numbers and passing their pointers to the dwGetParam.
        """
        output = int64(0)
        spcm_dwGetParam_i64(card, command, byref(output))
        return output.value

    @staticmethod
    def _card_reset(card):
        """
        Resets one of the awg internal cards.
        @param card: The card index.
        """
        spcm_dwSetParam_i64(card, SPC_M2CMD, M2CMD_CARD_RESET)

    def _get_error_msg(self, card):
        """
        Gets the error message from an internal awg card, translates it to a readable string and prints it with the
        logger.

        @param card: The card index.
        @return int: Error code. 0: no error occurred, 1: error occured.
        """
        error_reg, error_val = uint32(0), int32(0)
        error_code = spcm_dwGetErrorInfo_i32(
            card, byref(error_reg), byref(error_val), byref(self._error_msg)
        )
        if error_code != ERR_OK:
            message_values = (
                card,
                error_reg.value,
                error_val.value,
                self._error_msg.value.decode("utf-8"),
            )
            self.log.error(
                "Error occured on card {}, at register {} with value {}: {}.".format(
                    *message_values
                )
            )
            return -1
        else:
            return 0


class AbstractWaveform(object):
    """
    Helpers class used to bridge between user-friendly dictionaries and the dictionaries that need to be passed
    to the _create_waveform_buffers method of the SpectrumNetbox class. The number of cards contained in the awg and
    the index of the awg master card need to be passed when initializing the class.
    """

    def __init__(self, awg_card_num, mastercard_idx):
        self._awg_card_num = awg_card_num
        self._mastercard_idx = mastercard_idx
        self._analog_waveforms = dict()
        self._digital_waveforms = dict()
        self._digital_map = dict()
        self._required_channels = []

    def assign_ao_waveforms(self, waveform_dict):
        """
        A dictionary of waveforms. The keys 0,1,2,... correspond to the analog output channels.

        @param dict waveform_dict: The dictionary of analog waveforms.
        """
        self._analog_waveforms.update(waveform_dict)
        self._required_channels += list(self._analog_waveforms.keys())

    def assign_do_waveforms(self, waveform_dict, digital_out_map):
        """
        A dictionary of waveforms. The keys 0,1,2,... correspond to the digital output channels.

        @param dict waveform_dict: The dictionary of analog waveforms.
        @param dict digital_output_map: The digital_output_map, as specified in the _create_waveform_buffers method.
        """
        self._digital_waveforms.update(waveform_dict)
        self._digital_map = digital_out_map
        self._required_channels += [
            self._mastercard_idx + ao_chan + 1 for ao_chan in self._digital_map.keys()
        ]

    def generate_awg_waveform(self):
        """
        Generates the waveform dictionaries as required by the SpectrumNetbox._create_waveform_buffers method

        @return tuple: Returns the analog waveform dictionary and the digital waveform dictionary in awg-friendly format,
                       and the digital_output_map.
        """
        if self._analog_waveforms:
            analog_out = defaultdict(dict)
            for key, waveform in self._analog_waveforms.items():
                card_index, channel_index = (
                    key // self._awg_card_num,
                    key % self._awg_card_num,
                )
                analog_out[card_index].update({channel_index: waveform})
        else:
            analog_out = dict()

        digital_out = dict()
        if self._digital_waveforms and self._digital_map:
            for channel, do_out_list in self._digital_map.items():
                # All the do waveforms assigned to the same channel need to have the same length
                do_matrix = np.zeros(
                    (len(do_out_list), len(self._digital_waveforms[do_out_list[0]])),
                    dtype=np.int64,
                )
                for ii, do_out in enumerate(do_out_list):
                    do_matrix[ii, :] = self._digital_waveforms[do_out]
                digital_out[channel] = do_matrix

        return analog_out, digital_out, self._digital_map

    def required_cards(self):
        """
        Gets the cards used by the loaded waveforms.

        @return list: The list of required cards.
        """
        cardifier = lambda chan: chan // self._awg_card_num
        return sorted(set(map(cardifier, self._required_channels)))

    def required_channels(self):
        """
        Gets the channels used by the loaded waveforms.

        @return list: The list of required channels.
        """
        return sorted(set(self._required_channels))

    def analog_waveforms(self):
        """
        Returns the currently loaded analog waveforms
        @return dict: The dictionary of loaded analog waveforms.
        """
        return self._analog_waveforms

    def digital_waveforms(self):
        """
        Returns the currently loaded digital waveforms
        @return dict: The dictionary of loaded digital waveforms.
        """
        return self._digital_waveforms


class CardCollection(object):
    """
    Small class that is used for code readability of the loaded spectrum cards.
    """

    def __init__(self):
        self._masteridx = None
        self._cards = []
        self._chans_per_card = []
        self.max_samprate = None
        self.memsize = None
        self.bytes_persample = None
        self.adc_resolution = None
        self.maxsegments = None
        self.netbox_type = None
        self.starHub = None
        self._constraints = None

    def add_card(self, card, is_master):
        self._cards.append(card)
        if is_master:
            self._masteridx = len(self._cards) - 1

    def master(self):
        """
        Return the master card.
        """
        return self._cards[self._masteridx] if self._masteridx else None

    def card(self, index):
        """
        Return any card
        """
        return self._cards[index] if index < len(self._cards) else None

    def slavecard(self, index):
        """
        Return the slave card
        """
        if index == self._masteridx:
            self.log.info("Requested the index of the master card.")
            return -1
        else:
            return self._cards[index]

    def card_amount(self):
        """
        Return number of loaded cards
        """
        return len(self._cards)

    def set_chan_num(self, num):
        self._chans_per_card.append(num)

    def get_chan_num(self, card_idx):
        return self._chans_per_card[card_idx]

    def get_info(self):
        if self.master():
            card = self.master()
        elif self.card(0):
            card = self.card(0)
        else:
            self.log.error("No cards are loaded.")
            return -1
        self.max_samprate = self._command_get(card, SPC_PCISAMPLERATE)
        self.memsize = self._command_get(card, SPC_PCIMEMSIZE)
        self.bytes_persample = self._command_get(card, SPC_MIINST_BYTESPERSAMPLE)
        self.adc_resolution = self._command_get(card, SPC_MIINST_MAXADCVALUE)
        self.maxsegments = self._command_get(card, SPC_SEQMODE_AVAILMAXSEGMENT)

        card_type = self._command_get(card, SPC_NETBOX_TYPE)
        if card_type != 0:
            series, family, speed_grade, channels = (
                (card_type >> 24) & 0xFF,
                (card_type >> 16) & 0xFF,
                (card_type >> 8) & 0xFF,
                card_type & 0xFF,
            )
            series = hex(series)[2:]
            family = hex(family)[2:]
            self.netbox_type = (
                f"DN{series}.{family}{speed_grade}-{str(channels).zfill(2)}"
            )

    @staticmethod
    def _command_get(card, command):
        """
        Function to stop keeping to write numbers and passing their pointers to the dwGetParam.
        """
        output = int64(0)
        spcm_dwGetParam_i64(card, command, byref(output))
        return output.value

    def _set_constraints(self):
        self._constraints = PulserConstraints()
        self._constraints.sample_rate.min = 50e6
        self._constraints.sample_rate.max = self.max_samprate
        self._constraints.a_ch_amplitude.min = 80e-3
        self._constraints.a_ch_amplitude.max = 2.5

    @property
    def constraints(self):
        return self._constraints

    @property
    def masteridx(self):
        return self._masteridx
