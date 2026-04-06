# Test GeoRasters
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
raster = os.path.join(DATA, 'pre1500.tif')

def test_main():
    import georasters as gr
    raster = os.path.join(DATA, 'pre1500.tif')
    A = gr.from_file(raster)
    assert A.count() == 2277587
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
    data1 = gr.GeoRaster(data.raster[:int(data.shape[0]/2),:], data.geot,
                          nodata_value=data.nodata_value, projection=data.projection, datatype=data.datatype)
    data2 = gr.GeoRaster(data.raster[int(data.shape[0]/2):,:], (xmin,xsize,x,ymax+ysize*data.shape[0]/2,y,ysize),
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

def test_stats2():
    import georasters as gr
    raster = os.path.join(DATA, 'pre1500.tif')
    data = gr.from_file(raster)
    assert data.argmin() == data.raster.argmin()

def test_stats3():
    import georasters as gr
    raster = os.path.join(DATA, 'pre1500.tif')
    data = gr.from_file(raster)
    assert data.sum() == data.raster.sum()

def test_stats4():
    import georasters as gr
    raster = os.path.join(DATA, 'pre1500.tif')
    data = gr.from_file(raster)
    assert data.max() == data.raster.max()

def test_stats5():
    import georasters as gr
    raster = os.path.join(DATA, 'pre1500.tif')
    data = gr.from_file(raster)
    assert data.min() == data.raster.min()

def test_stats6():
    import georasters as gr
    raster = os.path.join(DATA, 'pre1500.tif')
    data = gr.from_file(raster)
    assert data.median() == np.ma.median(data.raster)

def test_stats7():
    import georasters as gr
    raster = os.path.join(DATA, 'pre1500.tif')
    data = gr.from_file(raster)
    assert data.std() == data.raster.std()

def test_stats8():
    import georasters as gr
    raster = os.path.join(DATA, 'pre1500.tif')
    data = gr.from_file(raster)
    assert data.var() == data.raster.var()

def test_stats9():
    import georasters as gr
    raster = os.path.join(DATA, 'pre1500.tif')
    data = gr.from_file(raster)
    assert data.prod() == data.raster.prod()

def test_stats10():
    import georasters as gr
    raster = os.path.join(DATA, 'pre1500.tif')
    data = gr.from_file(raster)
    assert data.count() == data.raster.count()

def test_aggregate_preserves_mask():
    """Issues #72, #51: aggregate() must propagate the mask to the reduced raster."""
    import georasters as gr
    import numpy as np
    # 4x4 raster with a 2x2 block of masked values in the top-left
    data = np.ma.array(
        [[1, 2, 3, 4],
         [5, 6, 7, 8],
         [9, 10, 11, 12],
         [13, 14, 15, 16]],
        mask=[[True, True, False, False],
              [True, True, False, False],
              [False, False, False, False],
              [False, False, False, False]],
        dtype=float
    )
    geo = gr.GeoRaster(data, (0.0, 1.0, 0.0, 4.0, 0.0, -1.0), nodata_value=-9999.0)
    agg = geo.aggregate((2, 2))
    # Result should be 2x2; top-left block was all masked → still masked
    assert agg.raster.shape == (2, 2)
    assert agg.raster.mask[0, 0], "Top-left block (all masked) should remain masked"
    assert not agg.raster.mask[0, 1], "Top-right block (no masked pixels) should be unmasked"
    assert not agg.raster.mask[1, 0], "Bottom-left block should be unmasked"
    assert not agg.raster.mask[1, 1], "Bottom-right block should be unmasked"

def test_block_reduce_preserves_mask():
    """Issues #72, #51: block_reduce() must propagate the mask to the reduced raster."""
    import georasters as gr
    import numpy as np
    data = np.ma.array(
        [[1, 2, 3, 4],
         [5, 6, 7, 8],
         [9, 10, 11, 12],
         [13, 14, 15, 16]],
        mask=[[True, True, False, False],
              [True, True, False, False],
              [False, False, False, False],
              [False, False, False, False]],
        dtype=float
    )
    geo = gr.GeoRaster(data, (0.0, 1.0, 0.0, 4.0, 0.0, -1.0), nodata_value=-9999.0)
    reduced = geo.block_reduce((2, 2), how=np.mean)
    assert reduced.raster.shape == (2, 2)
    assert reduced.raster.mask[0, 0], "Top-left block (all masked) should remain masked"
    assert not reduced.raster.mask[0, 1], "Top-right block should be unmasked"

def test_nodata_none_roundtrip(tmp_path):
    """Issues #47, #34: opening and immediately saving a raster with nodata=None should not crash."""
    import georasters as gr
    import warnings
    from osgeo import gdal, gdal_array
    # Build a small raster with no nodata value set (simulate files without NDV)
    driver = gdal.GetDriverByName('GTiff')
    path = str(tmp_path / 'no_ndv.tif')
    ds = driver.Create(path, 4, 4, 1, gdal.GDT_Int16)
    ds.SetGeoTransform((0.0, 1.0, 0.0, 4.0, 0.0, -1.0))
    from osgeo import osr
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    ds.SetProjection(srs.ExportToWkt())
    import numpy as np
    ds.GetRasterBand(1).WriteArray(np.arange(16, dtype=np.int16).reshape(4, 4))
    # Intentionally do NOT call SetNoDataValue
    ds.FlushCache()
    ds = None
    # Load — should warn but not crash
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        data = gr.from_file(path)
        assert any(issubclass(warning.category, UserWarning) for warning in w), \
            "Expected a UserWarning about missing nodata value"
    # Save — should not crash
    out = str(tmp_path / 'saved.tif')
    data.to_tiff(out)
    assert (tmp_path / 'saved.tif').exists()

def test_to_tiff_no_duplicate_extension(tmp_path):
    """Issue #46: to_tiff should not double-append .tif extension."""
    import georasters as gr
    data = gr.from_file(os.path.join(DATA, 'pre1500.tif'))
    # Pass filename WITH .tif extension
    out = str(tmp_path / 'output.tif')
    data.to_tiff(out)
    assert (tmp_path / 'output.tif').exists(), "output.tif should exist"
    assert not (tmp_path / 'output.tif.tif').exists(), "output.tif.tif should NOT exist"
    # Pass filename WITHOUT extension — should still produce .tif
    out2 = str(tmp_path / 'output2')
    data.to_tiff(out2)
    assert (tmp_path / 'output2.tif').exists(), "output2.tif should exist"
