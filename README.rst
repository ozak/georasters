GeoRasters
===========

|BuildStatus|_ 
|CoverageStatus|_
|PyPiVersion|_
|PyPiDownloads|_

The ``GeoRasters`` package is a python module that provides a fast and flexible
tool to work with GIS raster files. It provides the GeoRaster class, which makes working with rasters quite transparent and easy.
In a way it tries to do for rasters what GeoPandas does for geometries.

It includes tools to 

- Merge rasters
- Plot rasters
- Extract information from rasters
- Given a point (lat,lon) find its location in a raster
- Aggregate rasters to lower resolutions
- Align two rasters of different sizes to common area and size
- Get all the geographical information of raster
- Create GeoTiff files easily
- Load GeoTiff files as masked numpy rasters

Install
-------

.. code-block:: python
    
    pip install georasters
   
Example Usage
-------------

.. code-block:: python
    
    import georasters as gr
    # Get info on raster
    NDV, xsize, ysize, GeoT, Projection, DataType = gr.get_geo_info(raster)
    
    # Load raster
    data = load_tiff(raster)
       
    # Find location of point (x,y) on raster, e.g. to extract info at that location
    col, row = gr.map_pixel(x,y,GeoT[1],GeoT[-1], GeoT[0],GeoT[3])
    value = data[row,col]
    
    # Agregate raster by summing over cells in order to increase pixel size by e.g. 10
    gr.aggregate(data,NDV,(10,10))
    
    # Align two rasters
    data2 = load_tiff(raster2)
    (alignedraster_o, alignedraster_a, GeoT_a) = gr.align_rasters(raster, raster2, how=np.mean)
    
    # Create GeoRaster
    A=gr.GeoRaster(data, GeoT, nodata_value=NDV)

    # Load another raster
    NDV, xsize, ysize, GeoT, Projection, DataType = gr.get_geo_info(raster2)
    data = load_tiff(raster2)
    B=gr.GeoRaster(data2, GeoT, nodata_value=NDV)
    
    # Plot Raster
    A.plot()
    
    # Merge both rasters and plot
    C=B.merge(A)
    C.plot()
    
Issues
------

Find a bug? Report it via github issues by providing

- a link to download the smallest possible raster and vector dataset necessary to reproduce the error
- python code or command to reproduce the error
- information on your environment: versions of python, gdal and numpy and system memory

.. |BuildStatus| image:: https://api.travis-ci.org/ozak/georasters.png
.. _BuildStatus: https://travis-ci.org/ozak/georasters

.. |CoverageStatus| image:: https://coveralls.io/repos/ozak/georasters/badge.png
.. _CoverageStatus: https://coveralls.io/r/ozak/georasters

.. |PyPiVersion| image:: https://pypip.in/v/georasters/badge.png
.. _PyPiVersion: http://pypi.python.org/pypi/georasters

.. |PyPiDownloads| image:: https://pypip.in/d/georasters/badge.png
.. _PyPiDownloads: http://pypi.python.org/pypi/georasters