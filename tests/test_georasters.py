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
    assert A.projection.ExportToProj4().strip() == '+proj=longlat +datum=WGS84 +no_defs'

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

# ---------------------------------------------------------------------------
# Issue #4 — reproject
# ---------------------------------------------------------------------------

def test_reproject_changes_crs():
    """Issue #4: reproject() should return a GeoRaster in the requested CRS."""
    import georasters as gr
    data = gr.from_file(os.path.join(DATA, 'pre1500.tif'))
    # Source is WGS84 (EPSG:4326); reproject to Web Mercator (EPSG:3857)
    reprojected = data.reproject(3857)
    src_wkt  = data.projection.ExportToProj4()
    dst_wkt  = reprojected.projection.ExportToProj4()
    assert src_wkt != dst_wkt, "Projection should have changed"
    assert 'merc' in dst_wkt.lower() or '3857' in reprojected.projection.ExportToWkt(), \
        "Output should be Web Mercator"

def test_reproject_epsg_string():
    """reproject() should accept 'EPSG:NNNN' strings."""
    import georasters as gr
    data = gr.from_file(os.path.join(DATA, 'pre1500.tif'))
    reprojected = data.reproject('EPSG:3857')
    assert 'merc' in reprojected.projection.ExportToProj4().lower() or \
        '3857' in reprojected.projection.ExportToWkt()

def test_reproject_preserves_nodata():
    """reproject() must carry the nodata value through to the output."""
    import georasters as gr
    data = gr.from_file(os.path.join(DATA, 'pre1500.tif'))
    reprojected = data.reproject(3857)
    assert reprojected.nodata_value == data.nodata_value

def test_reproject_output_is_masked():
    """reproject() output raster should be a masked array with some masked cells."""
    import georasters as gr
    data = gr.from_file(os.path.join(DATA, 'pre1500.tif'))
    reprojected = data.reproject(3857)
    assert isinstance(reprojected.raster, np.ma.MaskedArray)
    assert reprojected.raster.mask.any(), "Reprojected raster should have some masked (nodata) cells"

def test_reproject_bounds_are_valid():
    """reproject() output should have a valid, non-zero spatial extent."""
    import georasters as gr
    data = gr.from_file(os.path.join(DATA, 'pre1500.tif'))
    reprojected = data.reproject(3857)
    assert reprojected.xmin < reprojected.xmax
    assert reprojected.ymin < reprojected.ymax

def test_reproject_invalid_resampling_raises():
    """reproject() should raise ValueError for unknown resampling method."""
    import georasters as gr
    import pytest
    data = gr.from_file(os.path.join(DATA, 'pre1500.tif'))
    with pytest.raises(ValueError):
        data.reproject(3857, resampling='bogus')

def test_reproject_osr_input():
    """reproject() should accept an osr.SpatialReference as dst_crs."""
    import georasters as gr
    from osgeo import osr
    data = gr.from_file(os.path.join(DATA, 'pre1500.tif'))
    dst = osr.SpatialReference()
    dst.ImportFromEPSG(3857)
    reprojected = data.reproject(dst)
    assert 'merc' in reprojected.projection.ExportToProj4().lower() or \
        '3857' in reprojected.projection.ExportToWkt()

def test_reproject_roundtrip_shape():
    """reproject() result should have non-trivial shape matching expected output dims."""
    import georasters as gr
    data = gr.from_file(os.path.join(DATA, 'pre1500.tif'))
    reprojected = data.reproject(3857)
    assert reprojected.shape[0] > 0 and reprojected.shape[1] > 0

def test_reproject_nan_nodata_masked():
    """reproject() must correctly mask border pixels when nodata_value is nan (float raster)."""
    import georasters as gr
    import numpy as np
    data = gr.from_file(os.path.join(DATA, 'pre1500.tif'))
    # Cast to float so nan is a valid fill value, then set nodata to nan
    float_raster = np.ma.array(data.raster.data.astype(np.float64),
                               mask=data.raster.mask)
    geo = gr.GeoRaster(float_raster, data.geot, nodata_value=np.nan,
                       projection=data.projection, datatype=data.datatype)
    reprojected = geo.reproject(3857)
    assert isinstance(reprojected.raster, np.ma.MaskedArray)
    assert reprojected.raster.mask.any(), \
        "Border fill pixels (nan) should be masked when nodata_value=nan"

def test_reproject_against_gdal_cea():
    """
    Ground-truth test for reproject(): compare georasters output against an
    independent GDAL warp to CEA (ESRI:54034).

    Both pipelines must produce:
      - the same output shape
      - the same geotransform (within floating-point tolerance)
      - pixel values that agree at all valid (non-nodata) locations
    """
    import georasters as gr
    import numpy as np
    from osgeo import gdal, osr
    from rasterio.crs import CRS as RioCRS
    from rasterio.warp import reproject as rio_reproject, Resampling, calculate_default_transform
    from affine import Affine

    src_path = os.path.join(DATA, 'pre1500.tif')
    CEA = 'ESRI:54034'

    # ------------------------------------------------------------------
    # Reference: reproject with raw rasterio (no georasters code path)
    # ------------------------------------------------------------------
    src_ds  = gdal.Open(src_path)
    ndv     = src_ds.GetRasterBand(1).GetNoDataValue()
    geot    = src_ds.GetGeoTransform()
    src_srs = osr.SpatialReference()
    src_srs.ImportFromWkt(src_ds.GetProjectionRef())
    src_crs = RioCRS.from_wkt(src_srs.ExportToWkt())
    dst_crs = RioCRS.from_string(CEA)

    xmin, xsize, _, ymax, _, ysize = geot
    rows, cols = src_ds.RasterYSize, src_ds.RasterXSize
    ymin  = ymax + ysize * rows
    xmax  = xmin + xsize * cols

    ref_transform, ref_width, ref_height = calculate_default_transform(
        src_crs, dst_crs, cols, rows,
        left=xmin, bottom=ymin, right=xmax, top=ymax)

    import numpy as np
    src_arr  = src_ds.GetRasterBand(1).ReadAsArray().astype(np.float64)
    ref_arr  = np.full((ref_height, ref_width), ndv, dtype=np.float64)
    rio_reproject(
        src_arr, ref_arr,
        src_transform=Affine.from_gdal(*geot),
        src_crs=src_crs,
        dst_transform=ref_transform,
        dst_crs=dst_crs,
        src_nodata=ndv, dst_nodata=ndv,
        resampling=Resampling.bilinear,
    )
    ref_mask = (ref_arr == ndv)

    # ------------------------------------------------------------------
    # Under test: georasters reproject()
    # ------------------------------------------------------------------
    geo          = gr.from_file(src_path)
    reprojected  = geo.reproject(CEA, resampling='bilinear')
    gr_arr       = reprojected.raster.filled(ndv)
    gr_mask      = np.ma.getmaskarray(reprojected.raster)

    # --- shape must match exactly ---
    assert reprojected.shape == (ref_height, ref_width), (
        "Shape mismatch: georasters {} vs GDAL reference {}".format(
            reprojected.shape, (ref_height, ref_width)))

    # --- geotransform must match to sub-metre precision ---
    ref_geot = ref_transform.to_gdal()
    for i, (g, r) in enumerate(zip(reprojected.geot, ref_geot)):
        assert abs(g - r) < 1.0, (
            "geot[{}] mismatch: georasters {:.4f} vs reference {:.4f}".format(i, g, r))

    # --- masks must agree at every cell ---
    np.testing.assert_array_equal(
        gr_mask, ref_mask,
        err_msg="Nodata mask differs between georasters and GDAL reference")

    # --- valid pixel values must agree within bilinear resampling tolerance ---
    valid = ~ref_mask
    np.testing.assert_allclose(
        gr_arr[valid], ref_arr[valid], rtol=0, atol=1.0,
        err_msg="Pixel values differ by more than 1 unit at valid locations")

# ---------------------------------------------------------------------------
# Issue #13 — Spatial autocorrelation methods
#
# Strategy: build a tiny synthetic 5x5 GeoRaster (no masked cells), compute
# stats via georasters, then verify against:
#   (a) direct esda calls on the same array
#   (b) esda calls on the geopandas representation of the same raster
# ---------------------------------------------------------------------------

@pytest.fixture(scope='module')
def tiny_geo():
    """Synthetic 5x5 GeoRaster, no masked cells, WGS84, for fast SA tests."""
    import georasters as gr
    from osgeo import osr
    # Spatially autocorrelated pattern: gradient + noise
    data = np.array([
        [10, 12, 11, 13, 10],
        [20, 22, 21, 23, 20],
        [30, 32, 31, 33, 30],
        [40, 42, 41, 43, 40],
        [50, 52, 51, 53, 50],
    ], dtype=np.float32)
    raster = np.ma.array(data, mask=False)
    proj = osr.SpatialReference()
    proj.ImportFromEPSG(4326)
    return gr.GeoRaster(raster,
                        (0.0, 1.0, 0.0, 5.0, 0.0, -1.0),
                        nodata_value=-9999.0,
                        projection=proj,
                        datatype='Float32')

# --- helpers ----------------------------------------------------------------

def _direct_weights(nrows, ncols, rook=False):
    """Build row-standardised lat2W weights matching georasters' raster_weights."""
    from libpysal.weights import lat2W
    w = lat2W(nrows, ncols, rook=rook)
    w.transform = 'r'
    return w

def _geopandas_queen_weights(gdf):
    """Build row-standardised Queen weights from the geopandas representation."""
    from libpysal.weights import Queen
    w = Queen.from_dataframe(gdf, use_index=False)
    w.transform = 'r'
    return w

# --- _compressed ------------------------------------------------------------

def test_compressed_excludes_masked(tiny_geo):
    """_compressed() returns a plain ndarray of all non-masked values."""
    valid = tiny_geo._compressed()
    n_valid = int((~np.ma.getmaskarray(tiny_geo.raster)).sum())
    assert len(valid) == n_valid
    assert isinstance(valid, np.ndarray)
    assert not isinstance(valid, np.ma.MaskedArray)

# --- Moran's I --------------------------------------------------------------

def test_pysal_Moran_matches_direct_esda(tiny_geo):
    """Moran.I via georasters must equal esda.Moran on the same array/weights."""
    from esda import Moran
    tiny_geo.weights = None
    tiny_geo.pysal_Moran(permutations=0)
    y = tiny_geo._compressed()
    w = _direct_weights(5, 5)
    ref = Moran(y, w, permutations=0)
    np.testing.assert_almost_equal(tiny_geo.Moran.I, ref.I, decimal=10)

def test_pysal_Moran_matches_geopandas(tiny_geo):
    """Moran.I via georasters must equal Moran.I computed from geopandas representation."""
    from esda import Moran
    tiny_geo.weights = None
    tiny_geo.pysal_Moran(permutations=0)
    gdf = tiny_geo.to_geopandas(name='value')
    w_gp = _geopandas_queen_weights(gdf)
    ref = Moran(gdf['value'].values.astype(float), w_gp, permutations=0)
    np.testing.assert_almost_equal(tiny_geo.Moran.I, ref.I, decimal=5)

def test_pysal_Moran_seed_reproducible(tiny_geo):
    """Same seed must produce identical p_sim across two calls."""
    tiny_geo.weights = None
    tiny_geo.pysal_Moran(permutations=19, seed=42)
    p1 = tiny_geo.Moran.p_sim
    tiny_geo.weights = None
    tiny_geo.pysal_Moran(permutations=19, seed=42)
    p2 = tiny_geo.Moran.p_sim
    assert p1 == p2

def test_pysal_Moran_seed_does_not_change_I(tiny_geo):
    """Seeded permutations must not alter the point estimate I."""
    from esda import Moran
    tiny_geo.weights = None
    tiny_geo.pysal_Moran(permutations=0)
    I_base = tiny_geo.Moran.I
    tiny_geo.weights = None
    tiny_geo.pysal_Moran(permutations=19, seed=7)
    np.testing.assert_almost_equal(tiny_geo.Moran.I, I_base, decimal=10)

# --- Geary's C --------------------------------------------------------------

def test_pysal_Geary_matches_direct_esda(tiny_geo):
    """Geary.C via georasters must equal esda.Geary on the same array/weights."""
    from esda import Geary
    tiny_geo.weights = None
    tiny_geo.pysal_Geary(permutations=0)
    y = tiny_geo._compressed()
    w = _direct_weights(5, 5)
    ref = Geary(y, w, permutations=0)
    np.testing.assert_almost_equal(tiny_geo.Geary.C, ref.C, decimal=10)

def test_pysal_Geary_matches_geopandas(tiny_geo):
    """Geary.C via georasters must equal Geary.C from geopandas representation."""
    from esda import Geary
    tiny_geo.weights = None
    tiny_geo.pysal_Geary(permutations=0)
    gdf = tiny_geo.to_geopandas(name='value')
    w_gp = _geopandas_queen_weights(gdf)
    ref = Geary(gdf['value'].values.astype(float), w_gp, permutations=0)
    np.testing.assert_almost_equal(tiny_geo.Geary.C, ref.C, decimal=5)

def test_pysal_Geary_seed_reproducible(tiny_geo):
    """Same seed must produce identical p_sim across two calls."""
    tiny_geo.weights = None
    tiny_geo.pysal_Geary(permutations=19, seed=7)
    p1 = tiny_geo.Geary.p_sim
    tiny_geo.weights = None
    tiny_geo.pysal_Geary(permutations=19, seed=7)
    p2 = tiny_geo.Geary.p_sim
    assert p1 == p2

# --- G (no seed — esda.G does not support it) ------------------------------

def test_pysal_G_matches_direct_esda(tiny_geo):
    """G statistic via georasters must equal esda.G on the same array/weights."""
    from esda import G
    tiny_geo.weights = None
    tiny_geo.pysal_G(permutations=0)
    y = tiny_geo._compressed()
    w = _direct_weights(5, 5)
    ref = G(y, w, permutations=0)
    np.testing.assert_almost_equal(tiny_geo.G.G, ref.G, decimal=10)

# --- Local Moran (seed passed to esda natively) ----------------------------

def test_pysal_Moran_Local_mapped_to_grid(tiny_geo):
    """Moran_Local.Is must be mapped back as a GeoRaster of the same shape."""
    tiny_geo.weights = None
    tiny_geo.pysal_Moran_Local(permutations=0, seed=42)
    assert hasattr(tiny_geo.Moran_Local, 'Is')
    assert isinstance(tiny_geo.Moran_Local.Is, type(tiny_geo))
    assert tiny_geo.Moran_Local.Is.shape == tiny_geo.shape

def test_pysal_Moran_Local_values_match_direct_esda(tiny_geo):
    """Local Moran Is values must match esda.Moran_Local on same data."""
    from esda import Moran_Local
    tiny_geo.weights = None
    tiny_geo.pysal_Moran_Local(permutations=0, seed=None)
    y = tiny_geo._compressed()
    w = _direct_weights(5, 5)
    ref = Moran_Local(y, w, permutations=0)
    gr_Is = tiny_geo.Moran_Local.Is._compressed()
    # float32 raster → ~7 sig figs; decimal=6 matches that precision
    np.testing.assert_array_almost_equal(gr_Is, ref.Is, decimal=6)

def test_pysal_Moran_Local_seed_accepted(tiny_geo):
    """Moran_Local seed kwarg must be accepted without error (passed to esda)."""
    tiny_geo.weights = None
    tiny_geo.pysal_Moran_Local(permutations=0, seed=42)
    assert tiny_geo.Moran_Local is not None

# --- Local G (seed passed to esda natively) --------------------------------

def test_pysal_G_Local_mapped_to_grid(tiny_geo):
    """G_Local.Zs must be mapped back as a GeoRaster of the same shape."""
    tiny_geo.weights = None
    tiny_geo.pysal_G_Local(permutations=0, seed=42)
    assert hasattr(tiny_geo.G_Local, 'Zs')
    assert isinstance(tiny_geo.G_Local.Zs, type(tiny_geo))
    assert tiny_geo.G_Local.Zs.shape == tiny_geo.shape

def test_pysal_G_Local_seed_accepted(tiny_geo):
    """G_Local seed kwarg must be accepted without error (passed to esda)."""
    tiny_geo.weights = None
    tiny_geo.pysal_G_Local(permutations=0, seed=3)
    assert tiny_geo.G_Local is not None

# --- kwargs isolation -------------------------------------------------------

def test_rook_kwarg_does_not_bleed_to_esda(tiny_geo):
    """rook= is consumed by raster_weights() and must not reach esda constructors."""
    tiny_geo.weights = None
    tiny_geo.pysal_Moran(permutations=0, rook=True)
    assert tiny_geo.Moran is not None
