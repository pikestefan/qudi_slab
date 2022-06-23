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
import sys
from pathlib import Path
from .spcm_tools import *
from .pyspcm import *

class CardCollection(object):
    """
    Small class that is used for code readability of the loaded spectrum cards.
    """
    def __init__(self):
        self.masteridx = None
        self._cards = []
        self._chans_per_card = []
        self.max_samprate = None
        self.memsize = None
        self.bytes_persample = None
        self.adc_resolution = None
        self.maxsegments = None
        self.netbox_type = None
        self.starHub = None

    def add_card(self, card, is_master):
        self._cards.append(card)
        if is_master:
            self.masteridx = len(self._cards) - 1

    def master(self):
        """
        Return the master card.
        """
        return self._cards[self.masteridx] if self.masteridx else None

    def card(self, index):
        """
        Return any card
        """
        return self._cards[index] if index < len(self._cards) else None

    def slavecard(self, index):
        """
        Return the slave card
        """
        if index==self.masteridx:
            self.log.info("Requested the index of the master card. Aborted.")
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


# FIXME: for now, the inheritance gets only Base. Once everything is settled, check the pulser interface class and see
#  if the commands are -or can be made- compatible with this module.
class SpectrumNetbox(Base):
    # TODO: one day, maybe implement an autodetection if the ip is not provided.
    _card_ip = ConfigOption('card_ip', missing='error')
    _netbox_type = ConfigOption('netbox_type', 'DN2.663-04', missing='info')
    _max_frequency = ConfigOption('max_frequency', 100e6, missing='info')
    _max_ao_voltage = ConfigOption('max_ao_voltage_mV', 500, missing='info')

    _ip_addon = '::inst{}::INSTR'


    #Initalize some lists of parameters that are going to be used in the methods
    __amplitude_chans = [SPC_AMP0, SPC_AMP1, SPC_AMP2, SPC_AMP3]
    __filter_chans = [SPC_FILTER0, SPC_FILTER1, SPC_FILTER2, SPC_FILTER3]
    __stoplevel_chans = [SPC_CH0_STOPLEVEL, SPC_CH1_STOPLEVEL, SPC_CH2_STOPLEVEL, SPC_CH3_STOPLEVEL]
    __output_chan_enable = [SPC_ENABLEOUT0, SPC_ENABLEOUT1, SPC_ENABLEOUT2, SPC_ENABLEOUT3]

    def on_activate(self):
        #FIXME: maybe implement a network discovery. If so, it should be fast.
        if self._netbox_type == 'DN2.663-04':
            self._card_num = 2
            self._ao_chans_percard = 2
            self._min_clock_rate = MEGA(50)
        else:
            self.log.error("Unknown digitizerNETBOX/generatorNETBOX.")
            return -1

        self._netbox = CardCollection()

        for ii in range(self._card_num):
            card_address = bytes(self._card_ip + self._ip_addon.format(ii), 'utf-8')
            card = spcm_hOpen(create_string_buffer(card_address))
            if card is None:
                self.log.error("Could not load card.")
                return -1
            # Make sure that the card is an AO card.
            card_function = self._command_get(card, SPC_FNCTYPE)
            cardfeatures = self._command_get(card, SPC_PCIFEATURES)
            if (card_function != SPCM_TYPE_AO) or (cardfeatures & SPCM_FEAT_SEQUENCE) == 0:
                self.log.error("Something is wrong, you are supposed to have an AO card with sequence mode.")
                return -1
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
        self._get_device_info()

        self._netbox.starHub = spcm_hOpen(create_string_buffer(b"sync0"))

        if self._netbox.netbox_type != self._netbox_type:
            self.log.error(f"The expected netbox type is {self._netbox_type}, "
                           f"but the device returned {self._netbox.netbox_type}. Aborting.")
            self.on_deactivate()

    def on_deactivate(self):
        for ii in range(self._netbox.card_amount()):
            spcm_vClose(self._netbox.card(ii))
            self.log.info(f"Successfully close card {ii} of netbox {self._netbox_type}")
        if self.starHub is not None:
            spcm_vClose(self.starHub)
            self.log.info("Successfully closed starhub.")

    def load_sequence(self, waveform_sequences=None, digital_sequences=None, step_map=None):
        """
        Method that loads the waveforms into the AWG. Assumes that all the different sequences will be played
        synchronously, meaning that all the requested sequences need to have the same length.
        @param dict waveform_sequences: A dictionary containing the list of ao waveforms. The keys must be called:
                                        ch0, ch1, etc... and refer to the ao channels of the AWG.
        @param dict digital_sequences: A dictionary containing the list of digital waveforms. The keys must be called:
                                       x0, x1, etc... and refer to the digital I/O channels of the AWG.
        @param list step_map: A list containing the order of sequence steps to be replayed. For example, if there are
                              three sequences, step_map = [0,2,1] will play the first, the last and the middle one.
                              Default is sequential.
        """
        # Check that all the requested channels have the same number of steps.
        seq_lengths = []
        for dictionary in [waveform_sequences, digital_sequences]:
            if dictionary is not None:
                seq_lengths += list(map(len, dictionary.values()))
        seq_lengths = set(seq_lengths)
        if len(seq_lengths) > 1:
            self.log.error("All the requested sequences need to have the same number of segments.")
            return -1

        if step_map is None:
            step_map = [ii for ii in range(seq_lengths[0])]

        # Get which channels need to be activated for the ao waveforms.
        get_chan_num = lambda key: int(key[-1])
        chans_to_activate = list(map(get_chan_num, waveform_sequences.keys()))

        # Now check also the digital inputs. They will require the activation of the master card. They will first be
        # distributed on each analog channel to minimize the sacrificed bits.
        if digital_sequences is not None:
            masteridx = self._netbox.masteridx
            if len(digital_sequences.keys()) > 1:
                do_chans = [masteridx+1 + chnum for chnum in range(self._netbox.get_chan_num(masteridx))]
            else:
                do_chans = [masteridx+1]
        else:
            do_chans = []

        #Merge and get the unique numbers
        chans_to_activate = list(set(chans_to_activate + do_chans))
        self._chan_enable(*chans_to_activate)
        pass

    def set_chan_amplitude(self, channels, amplitudes):
        """
        Set the channel max output amplitude (i.e., set the gain of the amplifier).
        @param list channels: the list of the channels, numbered from 0 to 3
        @param list amplitudes: the list of the amplitudes in mV required for each channel.

        @return int: -1 if errors occured, 0 otherwise
        """
        # Make sure everything is sorted
        channels, amplitudes = (list(srlist) for srlist in zip(*sorted(zip(channels, amplitudes))))
        for channel, amplitude in zip(channels, amplitudes):
            card_idx = channel // self._card_num
            # TODO: this is not general. Apply the clipping only to the channels that need it (e.g., the channels
            #  used for IQ modulation on the SRS.
            if amplitude > self._max_ao_voltage:
                amplitude = self._max_ao_voltage
                self.log.warning(f"Requested an ao voltage exceeding the {self._max_ao_voltage} bound."
                                 "Clipping applied.")
            command = self.__amplitude_chans[channel % self._card_num]
            errorout = spcm_dwSetParam_i64(self._netbox.card(card_idx),
                                           command,
                                           int64(amplitude))
            errorout = self._get_error_msg(self._netbox.card(card_idx), errorout)
            if errorout:
                return -1
        return 0

    def set_output_filters(self, channels, filter_active):
        """
       Set the output filters.
       @param list channels: the list of channels, numbered from 0 to 3
       @param list filter_active: 0 or 1, specifies if the filter needs to be activated or not.

       @return int: -1 if errors occured, 0 otherwise
       """
        # Make sure everything is sorted
        channels, filter_active = (list(srlist) for srlist in zip(*sorted(zip(channels, filter_active))))
        for channel, value in zip(channels, filter_active):
            card_idx = channel // self._card_num
            errorout = spcm_dwSetParam_i64(self._netbox.card(card_idx), self.__filter_chans[channel % self._card_num],
                                           int64(value))
            errorout = self._get_error_msg(self._netbox.card(card_idx), errorout)
            if errorout:
                return -1
        return 0

    def set_clock_rate(self, card, clk_rate = MEGA(50)):
        if not isinstance(clk_rate, int):
            clk_rate = int(clk_rate)
            self.log.warning(f"Given clock rate {clk_rate} is not an integer. Rounding will occurr.")
        if clk_rate < self._min_clock_rate:
            self.log.error(f"Minimum clock rate is {self._min_clock_rate}, you provided {clk_rate}.")
            return -1

        errorout = spcm_dwSetParam_i64(card, SPC_SAMPLERATE, clk_rate)
        errorout = self._get_error_msg(card, errorout)
        if errorout:
            return -1

    def get_clock_rate(self, card):
        return self._command_get(card, SPC_SAMPLERATE)

    def activate_outputs(self, *channels):
        """
        Enables or disables channels. This affects the maximum available memory that each channel has for the waveforms.
        *channels represents the list of requested channels, indexed from 0 to n (to stay consistent with python indexing).
        If the physical channels of each card in the netbox are less than the required channels, but if there are enough
        outputs in the rest of the netbox, the analog output is adjusted accordingly.

        @params int channels: unpacked list of channel integers, from 0 to n.

        @return int: -1 if errors occured, 0 otherwise
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

    def play_single_waveform(self, analog_waveforms = dict(), digital_waveforms = dict(), clk_rate = MEGA(50)):

        cards = self._get_required_cards(analog_waveforms, digital_waveforms)
        if len(cards) > 1:
            synchmaster = self._netbox.starHub
            spcm_dwSetParam_i64(synchmaster, SPC_SYNC_ENABLEMASK, 0x3)
        else:
            synchmaster = cards[0]

        for card in cards:
            spcm_dwSetParam_i64(card, SPC_CARDMODE, SPC_REP_STD_SINGLE)

    def simple_test(self, mseconds_play, loops, clk_rate=MEGA(50) ,amps = [100,100]):
        card = self._netbox.card(1)

        spcm_dwSetParam_i64(card, SPC_M2CMD, M2CMD_CARD_RESET)

        #self.set_clock_rate(card, int(clk_rate))

        samples = self._waveform_padding(mseconds_play * 1e-3 * clk_rate)
        tax = np.linspace(0, mseconds_play, samples)
        out = np.zeros((samples,))
        out[tax < mseconds_play/2] = 1

        out2 = np.zeros((samples,))
        out2[tax < mseconds_play / 2] = -1

        analog_waveforms = {1: {0: out}}#, 1: out2}}
        do_out = np.zeros((1, samples), dtype=c_int16) + 1
        #do_out[0] = np.where(out >= 0.5, 1, 0)
        digital_waveforms = {0: do_out}

        spcm_dwSetParam_i64(card, SPC_SAMPLERATE, clk_rate)
        spcm_dwSetParam_i32(card, SPC_CARDMODE, SPC_REP_STD_SINGLE)
        self._chan_enable(2)
        spcm_dwSetParam_i32(card, SPC_MEMSIZE, samples)
        spcm_dwSetParam_i32(card, SPC_LOOPS, loops)
        self.activate_outputs(2)
        self.set_chan_amplitude([2], amps)

        #dwXMode = int32(SPCM_XMODE_DIGOUT | SPCM_XMODE_DIGOUTSRC_CH0 | SPCM_XMODE_DIGOUTSRC_BIT15)
        #spcm_dwSetParam_i32(card, SPCM_X0_MODE, dwXMode)

        spcm_dwSetParam_i64(card, SPC_TRIG_ORMASK, SPC_TMASK_SOFTWARE)  # SPC_TMASK_NONE)
        error = spcm_dwSetParam_i64(card, SPC_TRIG_ANDMASK, 0)
        self.log.info("Transferring...")
        self._create_waveform_buffers(analog_waveforms=analog_waveforms)
        self.log.info("Finished.")

        self.log.info(self._get_error_msg(card, error))

        spcm_dwSetParam_i32(card, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA)
        spcm_dwSetParam_i32(card, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER)

    def _create_waveform_buffers(self, analog_waveforms=dict(), digital_waveforms=dict()):
        """
        Given a bunch of waveforms, it prepares the buffers for data transfer. Channels have to be activated somewhere
        else, this function just prepares and returns the buffers for data transfer.
        @param dict waveform_sequences: A dictionary containing subdictionaries. The keys of the dictionary need to be
                                        integer numbers, and they describe which is the required card. Each
                                        subdictionary is indexed with the channels indices (i.e. 0,1,2...)
        @param dict digital_sequences:  A dictionary containing the digital waveforms. The keys are 0,1,2,...
                                        corresponding to which digital channels the digital waveforms will be assigned
                                        to. Each element of the dictionary is a numpy array with shape MxN, where M
                                        corresponds to number of digital waveforms assigned to the channel. If you have
                                        a device with 16 bits, they will first be assigned to the 15th bit, then the
                                        14th, and so on.
        """

        masteridx = self._netbox.masteridx
        maxADC = self._netbox.adc_resolution

        #First calculate how many cards are needed.
        cards_required = self._get_required_cards(analog_waveforms, digital_waveforms)

        # Now start preparing the buffers
        card_buffers = dict().fromkeys(cards_required)
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
                    self.log.error("The analog waveforms loaded on the same card need to have equal length.")
                    return -1
                else:
                    waveform_ao_size = ao_sizes[0]
            if card_idx == masteridx and digital_waveforms:
                req_do_channels = list(digital_waveforms.keys())

                # Check the lengths of the digital waveforms to be loaded, they need to be the same.
                do_sizes = [wform.shape[1] for wform in digital_waveforms.values()]
                if len(set(do_sizes)) > 1:
                    self.log.error("The digital waveforms loaded on the same card need to have equal length.")
                    return -1
                else:
                    waveform_do_size = do_sizes[0]
            # Check that both waveforms sizes match for the same card.
            if waveform_ao_size is not None and waveform_do_size is not None and (waveform_ao_size != waveform_do_size):
                self.log.error("The lengh of the analog waveforms must match the length of the digital waveforms on "
                               "the same card.")
                return -1
            waveform_size = waveform_ao_size if waveform_ao_size is not None else waveform_do_size

            tot_channels = sorted(set(req_ao_channels + req_do_channels))

            # Now combine all the waveforms and channels for the same card into one array.
            channel_matrix = np.zeros((waveform_size, len(tot_channels)), dtype=c_int16)
            for channel in tot_channels:
                # By default, use the full resolution of the instrument
                maxADC_chan = maxADC - 1
                if card_idx == masteridx and digital_waveforms and channel in digital_waveforms:
                    # If there are digital channels, you need to sacrifice one bit for each digital channel
                    maxADC_chan = (maxADC // 2 ** (len(digital_waveforms[channel]))) - 1
                    channel_do_matrix = (digital_waveforms[channel].T *
                                         2 ** (15 - np.arange(len(digital_waveforms[channel])))).astype(c_int16)
                    #                    The line above shifts puts each do signal in its right bit
                    channel_matrix[:, channel] = np.sum(channel_do_matrix, axis=1)
                if card_ao_waveforms and channel in card_ao_waveforms:
                    # Apply mask to avoid overriding the digital output samples
                    rescaled_ao = self._waveform_16_to_nbits(card_ao_waveforms[channel], maxADC_chan)
                    channel_matrix[:, channel] = channel_matrix[:, channel] + rescaled_ao
            # Reshape the matrix to interleave the sample arrays
            channel_matrix = channel_matrix.reshape(np.prod(channel_matrix.shape))

            buffer_size = (self._netbox.bytes_persample * waveform_size * len(tot_channels))
            pvData = pvAllocMemPageAligned(buffer_size)  # The page-aligned buffer that will be transferred
            # FIXME: here I'm casting to a 16 bit pointer. This is not general for any netbox.
            pnData = cast(pvData, ptr16)  # The pointer to the buffer, used to fill it up

            # Here a nice ctypes/C trick to transfer the whole array to the pnData pointer, without using for loops
            memmove(pnData, channel_matrix.ctypes.data, sizeof(c_int16) * len(channel_matrix))

            spcm_dwDefTransfer_i64(self._netbox.card(card_idx), SPCM_BUF_DATA, SPCM_DIR_PCTOCARD,
                                   0, pvData, 0, buffer_size)

    def _get_required_cards(self, analog_waveforms = dict(), digital_waveforms = dict()):
        # TODO: get also the required channels e.g. 0,1,2,3 from the this action
        masteridx = self._netbox.masteridx
        cards_required = []
        if analog_waveforms:
            cards_required += list(analog_waveforms.keys())
        if digital_waveforms:
            cards_required.append(masteridx)  # The card used for digital output is always the master card.
            for chan, do_waveform in digital_waveforms.items():
                if len(do_waveform.shape) != 2:
                    self.log.error(f"Failed at channel {chan}. The digital output waveform matrices need to be "
                                   f"two-dimensional. ")
                    return -1

        return sorted(set(cards_required))

    def _waveform_16_to_nbits(self, waveform, max_value):
        out_waveform = (max_value * waveform).astype(c_int16)
        if max_value < (2**15-1):
            sign_waveform = np.zeros(out_waveform.shape, dtype=c_int16)
            sign_waveform[out_waveform < 0] = 1
            out_waveform = (out_waveform & max_value) + (sign_waveform * (max_value + 1))
        return out_waveform

    def _waveform_padding(self, waveform_len):
        """
        Calculate the padding required to take the waveform to the correct length. The waveform length needs to be a
        multiple of 32.
        """
        outlen = 0
        if waveform_len % 32 > 0:
            outlen = 32 * (waveform_len // 32) + 32
        else:
            outlen = waveform_len
        return int(outlen)

    def _get_device_info(self):
        #FIXME: consider moving this to the CardCollection class, rather than doing it here.
        if self._netbox.master():
            card = self._netbox.master()
        elif self._netbox.card(0):
            card = self._netbox.card(0)
        else:
            self.log.error("No cards are loaded.")
            return -1
        self._netbox.max_samprate = self._command_get(card, SPC_PCISAMPLERATE)
        self._netbox.memsize = self._command_get(card, SPC_PCIMEMSIZE)
        self._netbox.bytes_persample = self._command_get(card, SPC_MIINST_BYTESPERSAMPLE)
        self._netbox.adc_resolution = self._command_get(card, SPC_MIINST_MAXADCVALUE)
        self._netbox.maxsegments = self._command_get(card, SPC_SEQMODE_MAXSEGMENTS)

        card_type = self._command_get(card, SPC_NETBOX_TYPE)
        if card_type != 0:
            series, family, speed_grade, channels = (card_type >> 24) & 0xff, (card_type >> 16) & 0xff, \
                                                    (card_type >> 8) & 0xff, card_type & 0xff
            series = hex(series)[2:]
            family = hex(family)[2:]
            self._netbox.netbox_type = f"DN{series}.{family}{speed_grade}-{str(channels).zfill(2)}"

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
            chans_to_activate[card_idx] += 2**(channel % self._card_num)

        for card_idx, chan_value in chans_to_activate.items():
            if chan_value != 0:
                errorout = spcm_dwSetParam_i64(self._netbox.card(card_idx), SPC_CHENABLE, chan_value)
                errorout = self._get_error_msg(self._netbox.card(card_idx), errorout)
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
        channels, stoplevels = (list(srlist) for srlist in zip(*sorted(zip(channels, stoplevels))))

        for channel, stoplevel in zip(channels, stoplevels):
            card_idx = channel // self._card_num
            errorout = spcm_dwSetParam_i64(self._netbox.card(card_idx),
                                           self.__stoplevel_chans[channel % self._card_num],
                                           int64(stoplevel))

            errorout = self._get_error_msg(self._netbox.card(card_idx), errorout)
            if errorout:
                return -1
        return 0

    def _count_active_channels(self, card):
        """
        Return the active channel count for a card.
        """
        return self._command_get(card, SPC_CHCOUNT)

    @staticmethod
    def _command_get(card, command):
        """
        Function to stop keeping to write numbers and passing their pointers to the dwGetParam.
        """
        output = int64(0)
        spcm_dwGetParam_i64(card, command, byref(output))
        return output.value

    @staticmethod
    def _card_reset(card):
        spcm_dwSetParam_i64(card, SPC_M2CMD, M2CMD_CARD_RESET)

    def _get_error_msg(self, card, err):
        if err == ERR_OK:
            return 0
        else:
            error_text_buffer = create_string_buffer(ERRORTEXTLEN)
            spcm_dwGetErrorInfo_i32(card, 0, 0, byref(error_text_buffer))
            message = error_text_buffer.value
            self.log.error("Spectrum Netbox error." + message.decode("utf-8"))
            return -1



