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
Test Push
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
    __card_ormask_trigger_dict = {None: SPC_TMASK_NONE,
                                  'immediate': SPC_TMASK_SOFTWARE,  # The cards default, plays immediately after arming
                                  'ext0': SPC_TMASK_EXT0,
                                  'ext1': SPC_TMASK_EXT1
                                  }

    __card_andmask_trigger_dict = {None: SPC_TMASK_NONE,
                                   'ext0': SPC_TMASK_EXT0,
                                   'ext1': SPC_TMASK_EXT1,
                                   }
    __chan_digout_sources = [SPCM_XMODE_DIGOUTSRC_CH0, SPCM_XMODE_DIGOUTSRC_CH1, SPCM_XMODE_DIGOUTSRC_CH2]
    __bit_digout_sources = [SPCM_XMODE_DIGOUTSRC_BIT15, SPCM_XMODE_DIGOUTSRC_BIT14, SPCM_XMODE_DIGOUTSRC_BIT13]
    __digout_channels = [SPCM_X0_MODE, SPCM_X1_MODE, SPCM_X2_MODE]

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
        self._get_device_info()

        # Now open the starhub communication
        self._netbox.starHub = spcm_hOpen(create_string_buffer(b"sync0"))

        if self._netbox.netbox_type != self._netbox_type:
            self.log.error(f"The expected netbox type is {self._netbox_type}, "
                           f"but the device returned {self._netbox.netbox_type}. Aborting.")
            self.on_deactivate()

    def on_deactivate(self):
        for ii in range(self._netbox.card_amount()):
            spcm_vClose(self._netbox.card(ii))
            self.log.info(f"Successfully close card {ii} of netbox {self._netbox_type}")
        if self._netbox.starHub is not None:
            spcm_vClose(self._netbox.starHub)
            self.log.info("Successfully closed starhub.")

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

    def set_clock_rate(self, card, clk_rate=MEGA(50)):
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

    def load_sequence(self, analog_sequences=None, digital_sequences=None, step_map=list(), loops_list=list(),
                      digital_output_map=dict()):
        """
        Method that loads the waveforms into the AWG. Assumes that all the different sequences will be played
        synchronously, meaning that all the requested sequences need to have the same length.
        @param list analog_sequences: A list of dictionaries. Each dictionary containing the subdictionaries of ao waveforms.
                                        The keys must be called: 0, 1, etc..., and they refer to the AWG cards.
                                        The subdictionary keys, numbered 0, 1, etc..., refer to the ao channels of the
                                        single card.
        @param list digital_sequences: A dictionary containing the list of digital waveforms. The keys must be called:
                                       0, 1, etc... and refer to the digital I/O channels of the AWG.
        @param list step_map: A list containing the order of sequence steps to be replayed. For example, if there are
                              three sequences, step_map = [0,2,1] will play the first, the last and the middle one.
                              Default is sequential.
        @param list loops_list: A list containing the number of times each sequence step will be repeated.
        @param dict digital_output_map: optional dictionary that maps the source of the waveforms (i.e. the ao channels)
                                        to the physical do channels. I.e., digital_output_map = {0: [2,0]}, assigns
                                        the two waveforms that will be contained in ch0 waveform to the physical do x2
                                        and x0 channels, in this order.
        """
        # Here I assume that all the sequence steps use the same cards
        cards = self._get_required_cards(analog_sequences[0], digital_sequences[0])

        # Check the requested sequences
        len_ao, len_do = len(analog_sequences), len(digital_sequences)
        if (len_ao > 0 and len_do > 0) and (len_ao != len_do):
            self.log.error("The number of AO steps and DO steps must be the same.")
            return -1

        if len_ao > 0:
            steps_num = len_ao
        elif len_do > 0:
            steps_num = len_do

        if steps_num > self._netbox.maxsegments:
            self.log.error("The requested steps exceed the maximum allowed segmentation.")
            return -1

        # Get the memory requirements for each segment, then set up the sequence settings
        segment_memory_requirements = dict()
        for card_idx in cards:
            current_card = self._netbox.card(card_idx)
            active_chans = self._count_active_channels(current_card)
            min_segment_memory = 384 // active_chans
            max_segment_memory = int(self._netbox.memsize / active_chans / steps_num)

            segment_memory_requirements[card_idx] = [min_segment_memory, max_segment_memory]

            # Set up the sequence settings
            spcm_dwSetParam_i64(current_card, SPC_CARDMODE, SPC_REP_STD_SEQUENCE)
            spcm_dwSetParam_i64(current_card, SPC_SEQMODE_MAXSEGMENTS, steps_num)

        for step in range(steps_num):
            if len_ao:
                current_ao_waveform = analog_sequences[step]
            else:
                current_ao_waveform = dict()
            if len_do:
                current_do_waveform = digital_sequences[step]
            else:
                current_ao_waveform = dict()

            self._create_waveform_buffers(analog_waveforms=current_ao_waveform,
                                          digital_waveforms=current_do_waveform,
                                          sequence_step=step)





    def play_single_waveform(self, analog_waveforms=dict(), digital_waveforms=dict(), digital_output_map=dict(),
                             repeats=0, clk_rate=MEGA(50), channel_amplitudes=[], trigger_or_mask=list(),
                             trigger_and_mask=list()):

        card_indices = self._get_required_cards(analog_waveforms, digital_waveforms)
        channels_required = self._get_required_channels(analog_waveforms, digital_waveforms)

        if not channel_amplitudes:
            channel_amplitudes = [self._max_ao_voltage, ] * len(channels_required)

        # Check that all the waveforms have the same length
        for card_idx in card_indices:
            ao_wforms = []
            do_wforms = []
            if analog_waveforms and card_idx in analog_waveforms:
                ao_wforms += list(map(len, analog_waveforms[card_idx].values()))
        if digital_waveforms:
            do_wforms += [doform.shape[1] for doform in digital_waveforms.values()]
        wform_lengths = list(set(ao_wforms + do_wforms))
        if len(wform_lengths) > 1:
            self.log.error("All the analog and digital waveforms need to have the same length in this mode.")
            return -1
        else:
            memsize = wform_lengths[0]

        if (memsize % 32) != 0:
            self.log.error("The waveform length needs to be a multiple of 32. Use the method"
                           "SpectrumNetbox.waveform_padding to help you get the correct waveform length.")
            return -1

        if len(card_indices) > 1:
            synchmaster = self._netbox.starHub
            spcm_dwSetParam_i64(synchmaster, SPC_SYNC_ENABLEMASK, 0x3) # Harcoded, connects by default the only two cards we have
        else:
            synchmaster = self._netbox.card(card_indices[0])

        if self._netbox.masteridx in card_indices:
            mask_card = self._netbox.master()
        else:
            mask_card = self._netbox.card(card_indices[0])
            trigger_or_mask = [trig for trig in trigger_or_mask if trig in [None, 'immediate']]
            trigger_and_mask = [None]

        self.activate_outputs(*channels_required)

        self._chan_enable(*channels_required)
        self.set_chan_amplitude(channels_required, channel_amplitudes)

        for card_idx in card_indices:
            card = self._netbox.card(card_idx)
            self.set_clock_rate(card, clk_rate=clk_rate)
            spcm_dwSetParam_i64(card, SPC_CARDMODE, SPC_REP_STD_SINGLE)
            spcm_dwSetParam_i64(card, SPC_MEMSIZE, memsize)
            spcm_dwSetParam_i64(card, SPC_LOOPS, repeats)

            if card == mask_card:
                self.configure_ORmask(card, *trigger_or_mask)
                self.configure_ANDmask(card, *trigger_and_mask)
            else:
                self.configure_ORmask(card, *[])
                self.configure_ANDmask(card, *[])

        self._assign_digital_output_channels(digital_waveforms=digital_waveforms, digital_output_map=digital_output_map)

        self._create_waveform_buffers(analog_waveforms=analog_waveforms, digital_waveforms=digital_waveforms)

        self.start_card(synchmaster)
        self.arm_trigger(synchmaster)

    def simple_test(self, mseconds_play, loops, clk_rate=MEGA(50), amps=[500, 500]):

        spcm_dwSetParam_i64(self._netbox.card(0), SPC_M2CMD, M2CMD_CARD_RESET)
        spcm_dwSetParam_i64(self._netbox.card(1), SPC_M2CMD, M2CMD_CARD_RESET)

        #self.set_clock_rate(self._netbox.card(1), int(clk_rate))

        samples = self.waveform_padding(mseconds_play * 1e-3 * clk_rate)
        tax = np.linspace(0, mseconds_play, samples)
        out = np.zeros((samples,))

        halftime = len(tax[tax < mseconds_play/2])
        tinterval = tax[-1] - tax[halftime]
        out = np.cos(tax * 2*np.pi/tinterval)
        out2 = np.sin(tax * 2 * np.pi / tinterval)

        analog_waveforms = {1:{0:out2}}
        do_out = np.zeros((1, samples), dtype=c_int16)
        #do_out[0:len(tax)//2] = 1
        do_out[:, out2<=0] = 1

        do_out2 = ~do_out
        digital_waveforms = {1: do_out}

        self.play_single_waveform(analog_waveforms=analog_waveforms, digital_waveforms=digital_waveforms,
                                  repeats=loops, channel_amplitudes=amps, trigger_or_mask=['immediate'],
                                  digital_output_map={1:[2]},
                                  trigger_and_mask=[None], clk_rate=clk_rate)


    def configure_ORmask(self, card, *masks_to_enable):
        final_mask = 0
        for mask in masks_to_enable:
            if mask not in self.__card_ormask_trigger_dict:
                self.log.error("Requested register is not valid. Operation aborted.")
                return -1
            final_mask |= self.__card_ormask_trigger_dict[mask]

        spcm_dwSetParam_i64(card, SPC_TRIG_ORMASK, final_mask)

    def configure_ANDmask(self, card, *masks_to_enable):
        final_mask = 0
        for mask in masks_to_enable:
            if mask not in self.__card_andmask_trigger_dict:
                self.log.error("Requested register is not valid. Operation aborted.")
                return -1
            final_mask |= self.__card_ormask_trigger_dict[mask]

        spcm_dwSetParam_i64(card, SPC_TRIG_ANDMASK, final_mask)

    def send_software_trig(self, card):
        spcm_dwSetParam_i64(card, SPC_M2CMD, M2CMD_CARD_FORCETRIGGER)

    def start_card(self, card):
        spcm_dwSetParam_i64(card, SPC_M2CMD, M2CMD_CARD_START)

    def arm_trigger(self, card):
        spcm_dwSetParam_i64(card, SPC_M2CMD, M2CMD_CARD_ENABLETRIGGER)

    def _assign_digital_output_channels(self, digital_waveforms=dict(), digital_output_map=dict()):
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
                    digital_output_map[chan] = [last_used_chan + ii for ii in range(len(digital_waveforms[chan]))]
                    last_used_chan += digital_output_map[chan][-1] + 1

            for channel, xoutchannels in digital_output_map.items():
                for ii, xoutput in enumerate(xoutchannels):
                    x_chan_mode = SPCM_XMODE_DIGOUT | self.__chan_digout_sources[channel] | self.__bit_digout_sources[ii]
                    spcm_dwSetParam_i64(mastercard, self.__digout_channels[xoutput], x_chan_mode)

    def _create_waveform_buffers(self, analog_waveforms=dict(), digital_waveforms=dict(),
                                 sequence_step=None):
        """
        Given a bunch of waveforms, it prepares the buffers for data transfer. Channels have to be activated somewhere
        else, this function just prepares and returns the buffers for data transfer.
        @param dict analog_waveforms: A dictionary containing subdictionaries. The keys of the dictionary need to be
                                        integer numbers, and they describe which is the required card. Each
                                        subdictionary is indexed with the channels indices (i.e. 0,1,2...)
        @param dict digital_waveforms:  A dictionary containing the digital waveforms. The keys are 0,1,2,...
                                        corresponding to which digital channels the digital waveforms will be assigned
                                        to. Each element of the dictionary is a numpy array with shape MxN, where M
                                        corresponds to number of digital waveforms assigned to the channel. If you have
                                        a device with 16 bits, they will first be assigned to the 15th bit, then the
                                        14th, and so on.
        @param int sequence_step: Optional parameter to be fed when loading sequences. Defines the segment nummber.
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
            for chan_idx, channel in enumerate(tot_channels):
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
                    # Apply rescaling and adjust depending if the waveforms has 16, 15, or 14 usable bits
                    rescaled_ao = self._waveform_16_to_nbits(card_ao_waveforms[channel], maxADC_chan)
                    channel_matrix[:, chan_idx] = channel_matrix[:, chan_idx] + rescaled_ao
            # Reshape the matrix to interleave the sample arrays
            channel_matrix = channel_matrix.reshape(np.prod(channel_matrix.shape))

            buffer_size = (self._netbox.bytes_persample * waveform_size * len(tot_channels))
            pvData = pvAllocMemPageAligned(buffer_size)  # The page-aligned buffer that will be transferred
            # FIXME: here I'm casting to a 16 bit pointer. This is not general for any netbox.
            pnData = cast(pvData, ptr16)  # The pointer to the buffer, used to fill it up

            # Here a nice ctypes/C trick to transfer the whole array to the pnData pointer, without using for loops
            memmove(pnData, channel_matrix.ctypes.data, sizeof(c_int16) * len(channel_matrix))

            if sequence_step is not None:
                spcm_dwSetParam_i64(self._netbox.card(card_idx), SPC_SEQMODE_WRITESEGMENT, sequence_step)
                spcm_dwSetParam_i64(self._netbox.card(card_idx), SPC_SEQMODE_SEGMENTSIZE,
                                    waveform_size * len(tot_channels))

            error = spcm_dwDefTransfer_i64(self._netbox.card(card_idx), SPCM_BUF_DATA, SPCM_DIR_PCTOCARD,
                                           0, pvData, 0, buffer_size)
            error = spcm_dwSetParam_i64(self._netbox.card(card_idx), SPC_M2CMD,
                                        M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA)

    def _get_required_cards(self, analog_waveforms=dict(), digital_waveforms=dict()):
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

    def _get_required_channels(self, analog_waveforms=dict(), digital_waveforms=dict()):
        masteridx = self._netbox.masteridx

        # Assuming same number of channels per card
        chns_per_card = self._netbox.get_chan_num(0)

        channels_required = []
        if analog_waveforms:
            for card, card_waveforms in analog_waveforms.items():
                channels_required += [chn + card*chns_per_card for chn in card_waveforms.keys()]
        if digital_waveforms:
            channels_required += [chn + masteridx*chns_per_card for chn in digital_waveforms.keys()]

        return sorted(set(channels_required))

    def _waveform_16_to_nbits(self, waveform, max_value):
        out_waveform = (max_value * waveform).astype(c_int16)
        if max_value < (2**15-1):
            sign_waveform = np.zeros(out_waveform.shape, dtype=c_int16)
            sign_waveform[out_waveform < 0] = 1
            out_waveform = (out_waveform & max_value) + (sign_waveform * (max_value + 1))
        return out_waveform

    def waveform_padding(self, waveform_len):
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



