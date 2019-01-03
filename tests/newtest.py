# Test GeoRasters
from __future__ import division
import os, sys
import pytest
import georasters
import numpy as np
from osgeo import osr


DATA = os.path.join(os.path.dirname(georasters.__path__[0]), "tests/data")
raster = os.path.join(DATA, 'pre1500.tif')

# Create rasters
projection = osr.SpatialReference()
projection.ImportFromEPSG(4326)

# Example that works correctly
# Raster 1
A = np.array([[1]])
A1 = georasters.GeoRaster(A, (0, 1, 0, 0, 0,-1), nodata_value=-1, projection=projection)
# Raster 2
B = np.array([[3]])
B1 = georasters.GeoRaster(B, (2, 1, 0, -1, 0,-1), nodata_value=-1, projection=projection)

# Union
C1 = A1.union(B1)

# Expected
C = np.array([[1, -1, -1], [-1, -1, 3]])
CE = georasters.GeoRaster(C, (0, 1, 0, 0, 0,-1), nodata_value=-1, projection=projection)

assert((C1.raster.data==CE.raster.data).all())
assert((C1.raster.mask==CE.raster.mask).all())

# Example that works may not correctly
# Raster 1
A = np.array([[1]])
A1 = georasters.GeoRaster(A, (0, 1, 0, 0, 0,-1), nodata_value=-1, projection=projection)
# Raster 2
B = np.array([[3]])
B1 = georasters.GeoRaster(B, (2.5, 1, 0, -1, 0,-1), nodata_value=-1, projection=projection)

# Union
C1 = A1.union(B1)

# Expected
C = np.array([[1, -1, -1], [-1, -1, 3]])
CE = georasters.GeoRaster(C, (0, 1, 0, 0, 0,-1), nodata_value=-1, projection=projection)

assert((C1.raster.data==CE.raster.data).all())
assert((C1.raster.mask==CE.raster.mask).all())
