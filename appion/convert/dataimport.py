# **************************************************************************
# *
# * Authors:  Laura del Cano (ldelcano@cnb.csic.es)
# *           Yunior C. Fonseca Reyna (cfonseca@cnb.csic.es)
# *
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from os.path import exists

from pwem.objects import Coordinate
import pwem.emlib.metadata as md


class DogpickerImport:

    def __init__(self, protocol):
        self.protocol = protocol

    def importCoordinates(self, fileName, addCoordinate):
        print("In importCoordinates Appion with filename=%s" % fileName)
        if exists(fileName):
            mdata = md.MetaData()
            mdata.readPlain(fileName, 'xcoor ycoor')
            for objId in mdata:
                x = mdata.getValue(md.MDL_XCOOR, objId)
                y = mdata.getValue(md.MDL_YCOOR, objId)
                coord = Coordinate()
                coord.setPosition(x, y)
                addCoordinate(coord)
