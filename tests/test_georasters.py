# test zonal stats
from __future__ import division
import os, sys
import pytest
import numpy as np
import pandas as pd
from osgeo import gdal, gdalnumeric, ogr, osr
from gdalconst import *
from skimage.measure import block_reduce
import matplotlib.pyplot as plt
import georasters

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
raster = os.path.join(DATA, 'slope.tif')

def test_main():
    import georasters as gr
    A = gr.from_file(raster)
    assert A.mean()==56.813046574133502
    assert A.raster.count() == 6232
    assert A.min() == 0
    assert A.max() == 737.9317+1.6601562720097718e-06
    assert A.projection.ExportToProj4() == '+proj=aea +lat_1=43 +lat_2=48 +lat_0=34 +lon_0=-120 +x_0=600000 +y_0=0 +ellps=GRS80 +units=m +no_defs '

