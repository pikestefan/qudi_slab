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

from core.module import Base
from core.configoption import ConfigOption
import sys
from pathlib import Path
from .spcm_tools import *
from .pyspcm import *

class CardCollection():
    """
    Small class that is used for readable access of the loaded spectrum cards.
    """
    def __init__(self):
        self.master = None
        self.slaves = []

    def add_slave(self, card):
        self.slaves.append(card)

#FIXME: for now, the inheritance gets only Base. Once everything is settled, check the pulser interface class and see
# if the commands are -or can be made- compatible with this module.
class SpectrumNetbox(Base):
    # TODO: one day, maybe implement an autodetection if the ip is not provided
    _card_ip = ConfigOption('card_ip', missing='error')
    _netbox_type = ConfigOption('netbox_type', 'DN2.663-04', missing='info')

    _ip_addon = '::inst{}::INSTR'
    def on_activate(self):

        if self._netbox_type == 'DN2.663-04':
            self._card_num = 2
        else:
            self.log.error("Unknown digitizerNETBOX/generatorNETBOX.")
            return -1

        self._card_num = 2
        self._cardlist = CardCollection()
        #Now load the cards for real, and get their info
        for ii in range(self._card_num):
            card_address = bytes(self._card_ip + self._ip_addon.format(ii), 'utf-8')
            card = spcm_hOpen(create_string_buffer(card_address))
            if card is None:
                self.log.error("Could not load card.")
                return -1
            #Make sure that the card is an AO card.
            card_function = self._command_get(card, SPC_FNCTYPE)
            cardfeatures = self._command_get(card, SPC_PCIFEATURES)
            if (card_function != SPCM_TYPE_AO) or (cardfeatures & SPCM_FEAT_SEQUENCE) == 0:
                self.log.error("Something is wrong, you are supposed to have an AO card with sequence mode.")
                return -1
            has_starhub = self._command_get(card, SPCM_FW_MODEXTRA)
            if has_starhub:
                self._cardlist.master = card
            else:
                self._cardlist.add_slave(card)
        self.log.info("Loaded all the cards successfully.")
        if self._cardlist.master is None:
            self.log.warning("No card with starhub found.")
        self._get_device_info()

    def on_deactivate(self):
        spcm_vClose(self._cardlist.master)
        for card in self._cardlist.slaves:
            spcm_vClose(card)

    def _get_device_info(self):
        if self._cardlist.master:
            card = self._cardlist.master
        elif len(self._cardlist.slaves) > 0:
            card = self._cardlist.slaves[0]
        else:
            self.log.error("No cards are not loaded.")
            return -1
        max_samprate = self._command_get(card, SPC_PCISAMPLERATE)
        memsize = self._command_get(card, SPC_PCIMEMSIZE)
        bytes_persample = self._command_get(card, SPC_MIINST_BYTESPERSAMPLE)
        adc_resolution = self._command_get(card, SPC_MIINST_MAXADCVALUE)

        self._cardlist.max_samprate = max_samprate
        self._cardlist.memsize = memsize
        self._cardlist.bytes_persample = bytes_persample
        self._cardlist.adc_resolution = adc_resolution

    def _get_error_msg(self, card):
        error_text_buffer = create_string_buffer(ERRORTEXTLEN)
        spcm_dwGetErrorInfo_i32(card, 0, 0, byref(error_text_buffer))
        message = error_text_buffer.value
        self.log.error(message)

    def _command_get(self, card, command):
        output = int64(0)
        spcm_dwGetParam_i64(card, command, byref(output))
        return output.value

    def _card_reset(self, card):
        spcm_dwSetParam_i64(card, SPC_M2CMD, M2CMD_CARD_RESET)
