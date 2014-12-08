# test zonal stats
from __future__ import division
import os, sys
import pytest
import georasters
import numpy as np
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
    assert A.raster.count() == A.count()
    assert A.min() == 0
    assert A.projection.ExportToProj4() == '+proj=longlat +datum=WGS84 +no_defs '

def test_extract():
    import georasters as gr
    raster = os.path.join(DATA, 'pre1500.tif')
    data = gr.from_file(raster)
    (xmin,xsize,x,ymax,y,ysize)=data.geot
    (x,y)=(xmin+2507*xsize, ymax+1425*ysize)
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
    '''
    import matplotlib.pyplot as plt
    plt.figure()
    data1.plot()
    plt.savefig(os.path.join(DATA,'data1.png'))
    
    plt.figure()
    data2.plot()
    plt.savefig(os.path.join(DATA,'data2.png'))
    
    from rasterstats import zonal_stats
    import geopandas as gp
    import pandas as pd
    
    # Import shapefiles
    pathshp = os.path.join(DATA, 'COL.shp')
    dfcol=gp.GeoDataFrame.from_file(pathshp)
    pathshp = os.path.join(DATA, 'TUR.shp')
    dftur=gp.GeoDataFrame.from_file(pathshp)
    
    # Joint geopandas df
    df=dfcol.append(dftur)
    df.reset_index(drop=True,inplace=True)
    
    stats = zonal_stats(df, raster, copy_properties=True, all_touched=True, raster_out=True, opt_georaster=True)
    dfcol=pd.merge(dfcol,pd.DataFrame(data=stats), 
    
    
    
    
    
    '''
    assert (data1.union(data2).raster==data.raster).sum()==data.count()

def test_stats():
    import georasters as gr
    raster = os.path.join(DATA, 'pre1500.tif')
    data = gr.from_file(raster)
    assert data.argmax() == data.raster.argmax()
    assert data.argmin() == data.raster.argmin()
    assert data.sum() == data.raster.sum()
    assert data.max() == data.raster.max()
    assert data.min() == data.raster.min()
    assert data.median() == np.ma.median(data.raster)
    assert data.std() == data.raster.std()
    assert data.var() == data.raster.var()
    assert data.prod() == data.raster.prod()
    assert data.count() == data.raster.count()
