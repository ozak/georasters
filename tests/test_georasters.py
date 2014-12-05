# test zonal stats
from __future__ import division
import os, sys
import pytest
import georasters
from rasterstats import zonal_stats
'''
import numpy as np
import pandas as pd
from osgeo import gdal, gdalnumeric, ogr, osr
from gdalconst import *
from skimage.measure import block_reduce
import matplotlib.pyplot as plt
'''

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
raster = os.path.join(DATA, 'slope.tif')

def test_main():
    import georasters as gr
    A = gr.from_file(raster)
    assert A.raster.count() == 6232
    assert A.min() == 0
    assert A.projection.ExportToProj4() == '+proj=aea +lat_1=43 +lat_2=48 +lat_0=34 +lon_0=-120 +x_0=600000 +y_0=0 +ellps=GRS80 +units=m +no_defs '

def test_extract():
    import georasters as gr
    raster = os.path.join(DATA, 'pre1500.tif')
    data = gr.from_file(raster)
    (xmin,xsize,x,ymax,y,ysize)=data.geot
    (x,y)=(xmin+data.shape[1]/2*xsize, ymax+data.shape[0]/2*ysize)
    assert data.raster[gr.map_pixel(x,y,data.x_cell_size,data.y_cell_size,data.xmin,data.ymax)]==data.extract(x,y).max()
    assert data.raster[gr.map_pixel(x,y,data.x_cell_size,data.y_cell_size,data.xmin,data.ymax)]==data.map_pixel(x,y)

def test_union():
    import georasters as gr
    raster = os.path.join(DATA, 'pre1500.tif')
    data = gr.from_file(raster)
    (xmin,xsize,x,ymax,y,ysize)=data.geot
    data1 = gr.GeoRaster(data.raster[:data.shape[0]/2,:], data.geot, 
                          nodata_value=data.nodata_value, projection=data.projection, datatype=data.datatype)
    data2 = gr.GeoRaster(data.raster[data.shape[0]/2:,:], (xmin,xsize,x,ymax+ysize*data.shape[0]/2,y,ysize), 
                          nodata_value=data.nodata_value, projection=data.projection, datatype=data.datatype)
    assert (data1.union(data2).raster==data.raster).sum()==data.count()

