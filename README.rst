GeoRaters
===========

The ``Georasters`` is a python module that provides a fast and flexible
tool to work with GIS raster files. It includes tools to 

- Given a point (lat,lon) find its location in a raster
- Aggregate rasters to lower resolutions
- Align two rasters of different sizes to common area and size
- Get all the geographical information of raster
- Create GeoTiff files easily
- Load GeoTiff files as masked numpy rasters

Install
-------

.. code-block:: python
    
    pip install gisrastertools
   
Example Usage
-------------

.. code-block:: python
    
    from gisrastertools import *
    # Get info on raster
    NDV, xsize, ysize, GeoT, Projection, DataType = get_geo_info(raster)
    
    # Load raster
    data = load_tiff(raster)
       
    # Find location of point (x,y) on raster, e.g. to extract info at that location
    col, row = map_pixel(x,y,GeoT[1],GeoT[-1], GeoT[0],GeoT[3])
    value = data[row,col]
    
    # Agregate raster by summing over cells in order to increase pixel size by e.g. 10
    aggregate(data,NDV,(10,10))
    
    # Align two rasters
    data2 = load_tiff(raster2)
    (alignedraster_o, alignedraster_a, GeoT_a) = align_rasters(raster, raster2, how=np.mean)
   