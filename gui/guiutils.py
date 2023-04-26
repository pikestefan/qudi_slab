# -*- coding: utf-8 -*-

"""
This file contains the Qudi GUI module utility classes.

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

import pyqtgraph as pg
from pyqtgraph import QtCore
from pyqtgraph import functions as fn
import numpy as np
import weakref


class ColorBar(pg.GraphicsObject):
    """Create a ColorBar according to a previously defined color map.

    @param object pyqtgraph.ColorMap cmap: a defined colormap
    @param float width: width of the colorbar in x direction, starting from
                        the origin.
    @param numpy.array ticks: optional, definition of the relative ticks marks
    """

    def __init__(self, cmap, width, cb_min, cb_max):

        pg.GraphicsObject.__init__(self)

        # handle the passed arguments:
        self.stops, self.colors = cmap.getStops("float")
        self.stops = (self.stops - self.stops.min()) / self.stops.ptp()
        self.width = width

        # Constructs an empty picture which can be altered by QPainter
        # commands. The picture is a serialization of painter commands to an IO
        # device in a platform-independent format.
        self.pic = pg.QtGui.QPicture()

        self.refresh_colorbar(cb_min, cb_max)

    def refresh_colorbar(
        self, cb_min, cb_max, width=None, height=None, xMin=None, yMin=None
    ):
        """Refresh the appearance of the colorbar for a changed count range.

        @param float cb_min: The minimal count value should be passed here.
        @param float cb_max: The maximal count value should be passed here.
        @param float width: optional, with that you can change the width of the
                            colorbar in the display.
        """

        if width is None:
            width = self.width
        else:
            self.width = width

        #       FIXME: Until now, if you want to refresh the colorbar, a new QPainter
        #              object has been created, but I think that it is not necassary.
        #              I have to figure out how to use the created object properly.
        p = pg.QtGui.QPainter(self.pic)
        p.drawRect(self.boundingRect())
        p.setPen(pg.mkPen("k"))
        grad = pg.QtGui.QLinearGradient(
            width / 2.0, cb_min * 1.0, width / 2.0, cb_max * 1.0
        )
        for stop, color in zip(self.stops, self.colors):
            grad.setColorAt(1.0 - stop, pg.QtGui.QColor(*[255 * c for c in color]))
        p.setBrush(pg.QtGui.QBrush(grad))
        if xMin is None:
            p.drawRect(pg.QtCore.QRectF(0, cb_min, width, cb_max - cb_min))
        else:
            # If this picture whants to be set in a plot, which is going to be
            # saved:
            p.drawRect(pg.QtCore.QRectF(xMin, yMin, width, height))
        p.end()

        vb = self.getViewBox()
        # check whether a viewbox is already created for this object. If yes,
        # then it should be adjusted according to the full screen.
        if vb is not None:
            vb.updateAutoRange()
            vb.enableAutoRange()

    def paint(self, p, *args):
        """Overwrite the paint method from GraphicsObject.

        @param object p: a pyqtgraph.QtGui.QPainter object, which is used to
                         set the color of the pen.

        Since this colorbar object is in the end a GraphicsObject, it will
        drop an implementation error, since you have to write your own paint
        function for the created GraphicsObject.
        """
        # paint colorbar
        p.drawPicture(0, 0, self.pic)

    def boundingRect(self):
        """Overwrite the paint method from GraphicsObject.

        Get the position, width and hight of the displayed object.
        """
        return pg.QtCore.QRectF(self.pic.boundingRect())


class InteractiveColBar(pg.ImageItem):
    """
    Author: Lucio Stefan

    A partially hacky way to get an interactive colorbar with pyqtgraph=0.10.
    """

    def __init__(
        self, image, parent, lut=None, levels=None, multiplier=1, label="", width=None
    ):
        super().__init__()

        self.multiplier = multiplier

        self.parentref = weakref.ref(parent)
        self.parent = self.parentref()
        self.parent.disableAutoRange()
        self.parent.hideButtons()
        self.parent.setMouseEnabled(x=False, y=False)
        self.parent.setMenuEnabled(False)
        self.parent.hideAxis("bottom")
        self.parent.setLabel("left", label)
        if width is not None:
            self.parent.plotItem.setFixedWidth(width)
            self.parent.setFixedWidth(width)

        self.image_weakref = weakref.ref(image)
        self.reference_image = self.image_weakref()

        gradient = np.linspace(0, 1, 255)
        gradient = np.hstack((gradient[:, None], gradient[:, None]))
        self.setImage(image=gradient, axisOrder="row-major")

        if lut is not None:
            self.setLookupTable(lut)
        else:
            lut = self.reference_image.lut
            self.setLookupTable(lut)

        self.minlvl, self.maxlvl = None, None
        if levels is None:
            levels = [0, 1 / self.multiplier]

        self.regionAutoUpdate = True

        self._prepare_region(*levels)
        self.set_levels(*levels)
        self.parent.addItem(self)

        self.region.signalLineDragged.connect(self.update_image_levels)
        self.region.sigRegionChanged.connect(self.update_image_levels)

    def set_levels(self, minLevel, maxLevel):
        if minLevel != self.minlvl or maxLevel != self.maxlvl:
            self.minlvl = minLevel * self.multiplier
            self.maxlvl = maxLevel * self.multiplier
            self.parent.setYRange(self.minlvl, self.maxlvl, padding=0.01)
            self.setRect(QtCore.QRectF(0, self.minlvl, 1, self.maxlvl - self.minlvl))
            self.region.setBounds([self.minlvl, self.maxlvl])
            if self.regionAutoUpdate:
                self.region.setRegion([self.minlvl, self.maxlvl])

    def get_levels(self):
        minlvl, maxlvl = self.region.getRegion()
        return minlvl / self.multiplier, maxlvl / self.multiplier

    def update_image_levels(self):
        self.regionAutoUpdate = False
        minval, maxval = self.region.getRegion()
        self.reference_image.setLevels(
            [minval / self.multiplier, maxval / self.multiplier]
        )

    def _prepare_region(self, minlvl, maxlvl):
        self.region = DoubleClickableRegion(
            [minlvl, maxlvl],
            "horizontal",
            swapMode="block",
            movable=True,
            brush=pg.mkBrush(None),
            pen=pg.mkPen("r", width=2),
            hoverPen=pg.mkPen("g", width=2),
            hoverBrush=pg.mkBrush(37, 150, 190, 30),
            bounds=[minlvl, maxlvl],
        )
        self.region.setZValue(1000)
        self.region.lines[0].addMarker("|>", size=6)
        self.region.lines[1].addMarker("<|", size=6)

        self.region.signalDoubleClicked.connect(self._reset_region)

        self.parent.addItem(self.region)

    def _reset_region(self):
        self.regionAutoUpdate = True
        self.region.setBounds([self.minlvl, self.maxlvl])
        self.region.setRegion([self.minlvl, self.maxlvl])
        self.reference_image.setLevels(
            [self.minlvl / self.multiplier, self.maxlvl / self.multiplier]
        )


class DoubleClickableRegion(pg.LinearRegionItem):
    """
    Author: Lucio Stefan
    Custom LinearRegionItem class that can detect double clicks and use them to reset the bounds and position of the region.
    """

    signalDoubleClicked = QtCore.Signal()
    signalLineDragged = QtCore.Signal()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.lines[0].sigDragged.connect(self.signalLineDragged)
        self.lines[1].sigDragged.connect(self.signalLineDragged)

    def mouseDoubleClickEvent(self, ev):
        self.signalDoubleClicked.emit()

    def lineMoved(self, i):
        """
        Overriden copy-paste of the class original method, but doesn't emit the signalRegionChanged at the end.
        """
        if self.blockLineSignal:
            return

        # lines swapped
        if self.lines[0].value() > self.lines[1].value():
            if self.swapMode == "block":
                self.lines[i].setValue(self.lines[1 - i].value())
            elif self.swapMode == "push":
                self.lines[1 - i].setValue(self.lines[i].value())

        self.prepareGeometryChange()
