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

    def on_activate(self):
        #FIXME: maybe implement a network discovery. If so, it should be fast.
        if self._netbox_type == 'DN2.663-04':
            self._card_num = 2
            self._ao_chans_percard = 2
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

        self.starHub = spcm_hOpen(create_string_buffer(b"sync0"))

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
        cards_required = []
        if analog_waveforms:
            cards_required.append(list(analog_waveforms.keys()))
        if digital_waveforms:
            cards_required.append(masteridx)  # The card used for digital output is always the master card.
        cards_required = sorted(set(cards_required))

        card_buffers = dict().fromkeys(cards_required)
        for card_idx in cards_required:
            waveform_size = 0
            card_ao_waveforms = None
            req_ao_channels = []
            req_do_channels = []
            if analog_waveforms[card_idx]:
                card_ao_waveforms = analog_waveforms[card_idx]
                req_ao_channels = list(card_ao_waveforms.keys())
                waveform_size = len(card_ao_waveforms[req_ao_channels[0]])  # It's the same for all the waveforms in the card
            if card_idx == masteridx and digital_waveforms:
                req_do_channels = list(digital_waveforms.keys())
                waveform_size = len(req_do_channels[req_do_channels[0]])  # It's the same for all the waveforms in the card

            tot_channels = max(len(req_ao_channels), len(req_do_channels))
            buffer_size = (self._netbox.bytes_persample * waveform_size * tot_channels)

            pvData = pvAllocMemPageAligned(buffer_size)  # The page-aligned buffer that will be transferred
            #FIXME: here I'm casting to a 16 bit pointer. This is not general for any netbox.
            pnData = cast(pvData, ptr16)  # The pointer to the buffer, used to fill it up
            for buff_idx in range(0, waveform_size):
                for ii, channel in enumerate(tot_channels):
                    sample = 0
                    maxADC_chan = maxADC-1 # By default, use the full resolution of the instrument
                    if card_idx == masteridx and digital_waveforms and channel in digital_waveforms:
                        # If there are digital channels, you need to sacrifice one bit for each digital channel
                        maxADC_chan = (maxADC // 2**(len(digital_waveforms[channel]))) - 1
                        for jj, do_line in digital_waveforms[channel]:
                            sample |= (do_line[buff_idx] << (15-jj))  # Add the digital samples to the correct bit
                    if card_ao_waveforms and channel in card_ao_waveforms:
                        ao_sample = maxADC_chan * int(card_ao_waveforms[channel][buff_idx])
                        ao_sample &= 2**maxADC_chan  # For safety, to avoid overriding the already written do samples
                        sample |= ao_sample  # Combine the sample with the ao_sample

                    pnData[buff_idx + ii] = int16(sample)
            card_buffers[card_idx] = [pvData, pnData]
        return card_buffers

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
        chans_to_activate = [0, ] * self._card_num
        for channel in channels:
            card_idx = channel // self._card_num
            chans_to_activate[card_idx] += (channel % self._card_num+1)
            # The +1 is because the value for the FPGA register goes from 1 to 4, not from 0 to 3

        for ii, chan_value in enumerate(chans_to_activate):
            if chan_value != 0:
                errorout = spcm_dwSetParam_i64(self._netbox.card(ii), SPC_CHENABLE, chan_value)
                errorout = self._get_error_msg(self._netbox.card(ii), errorout)
                if errorout:
                    return -1
        return 0

    def set_chan_amplitude(self, channels, amplitudes):
        """
        Set the channel max output amplitude (i.e., set the gain of the amplifier).
        @param list channels: the list of the channels, numbered from 0 to 3
        @param list amplitudes: the list of the amplitudes required for each channel.

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
            errorout = spcm_dwSetParam_i64(self._netbox.card(card_idx),
                                           self.__amplitude_chans[channel % self._card_num],
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

    def set_clock_rate(self, card, clk_rate=50):
        if not isinstance(clk_rate, int):
            clk_rate = int(clk_rate)
            self.log.warning(f"Given clock rate {clk_rate} is not an integer. Rounding will occurr.")

        errorout = spcm_dwSetParam_i64(card, SPC_SAMPLERATE, clk_rate)
        errorout = self._get_error_msg(card, errorout)
        if errorout:
            return -1

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



