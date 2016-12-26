#!/usr/bin/env ipython
# coding: utf-8
'''
Copyright (C) 2014 Ömer Özak
This program defines functions that are useful for working with GIS data
Usage:

import georasters as gr

======================================================
Author:  Ömer Özak, 2013--2014 (ozak at smu.edu)
Website: http://omerozak.com
GitHub:  https://github.com/ozak/
======================================================

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

from __future__ import division
import georasters as gr
import hmi
import geopandas as gp
import geostats
import hmi
import os

# Path to data
pathtestdata = gr.__path__[0]+'/../tests/data/'
# Load raster data
data = gr.from_file(pathtestdata+'pre1500.tif')
# Load country geometries
col = gp.read_file(pathtestdata+'COL.shp')
df = data.clip(col, keep=True)
print(df)

# Select clipped raster
colraster = df.GeoRaster[0]
colraster.plot(cmap='Reds')

# Compute Global autocorrelation stats
colraster.pysal_G()
colraster.pysal_Gamma()
colraster.pysal_Geary()
colraster.pysal_Join_Counts()
colraster.pysal_Moran()
colraster.pysal_Geary()

# Print some Global autocorrelation stats
print(colraster.G.EG)
print(colraster.G.EG2)
print(colraster.G.EG_sim)
print(colraster.G.p_norm)
print(colraster.G.p_sim)
print(colraster.G.p_z_sim)

# Similar for other measures

# Compute Local autocorrelation stats
colraster.pysal_Moran_Local()
colraster.pysal_G_Local()

# Plot Local Moran
colraster.Moran_Local.Is.plot()
colraster.Moran_Local.p_z_sim.plot()

# Similar for other measures

