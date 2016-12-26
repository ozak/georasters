.. PyPI download statistics are broken, but the new PyPI warehouse makes PyPI
.. download statistics available through Google BigQuery
.. (https://bigquery.cloud.google.com).
.. Query to list package downloads by version:
..
   SELECT
     file.project,
     file.version,
     COUNT(*) as total_downloads,
     SUM(CASE WHEN REGEXP_EXTRACT(details.python, r"^([^\.]+\.[^\.]+)") = "2.6" THEN 1 ELSE 0 END) as py26_downloads,
     SUM(CASE WHEN REGEXP_EXTRACT(details.python, r"^([^\.]+\.[^\.]+)") = "2.7" THEN 1 ELSE 0 END) as py27_downloads,
     SUM(CASE WHEN REGEXP_EXTRACT(details.python, r"^([^\.]+)\.[^\.]+") = "3" THEN 1 ELSE 0 END) as py3_downloads,
   FROM
     TABLE_DATE_RANGE(
       [the-psf:pypi.downloads],
       TIMESTAMP("19700101"),
       CURRENT_TIMESTAMP()
     )
   WHERE
     file.project = 'georasters'
   GROUP BY
     file.project, file.version
   ORDER BY
     file.version DESC

GeoRasters
===========

|PyPiVersion|_
|PyPiDownloads|_
|BuildStatus|_ 
|CoverageStatus|_

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
- Clip raster using geometries
- Get zonal statistics using geometries
- Spatial analysis tools

Install
-------

.. code-block:: python
    
    pip install git+git://github.com/ozak/georasters.git
    pip install georasters
   
Example Usage: GeoRasters
-------------------------

.. code-block:: python
    
    import georasters as gr
    
    # Load data
    raster = './data/slope.tif'
    data = gr.from_file(raster)
    
    # Plot data
    data.plot()
    
    # Get some stats
    data.mean()
    data.sum()
    data.std()
    
    # Convert to Pandas DataFrame
    df = data.to_pandas()
    
    # Save transformed data to GeoTiff
    data2 = data**2
    data2.to_tiff('./data2')
    
    # Algebra with rasters
    data3 = np.sin(data.raster) / data2
    data3.plot()
    
    # Notice that by using the data.raster object, 
    # you can do any mathematical operation that handles 
    # Numpy Masked Arrays
    
    # Find value at point (x,y) or at vectors (X,Y)
    value = data.map_pixel(x,y)
    Value = data.map_pixel(X,Y)
    
Example Merge GeoRasters:
-------------------------

.. code-block:: python

    import os
    import georasters as gr
    import matplotlib.pyplot as plt

    DATA = "/path/to/tiff/files"

    # Import raster
    raster = os.path.join(DATA, 'pre1500.tif')
    data = gr.from_file(raster)
    (xmin, xsize, x, ymax, y, ysize) = data.geot

    # Split raster in two
    data1 = gr.GeoRaster(
        data.raster[:data.shape[0] / 2, :],
        data.geot,
        nodata_value=data.nodata_value,
        projection=data.projection,
        datatype=data.datatype,
    )
    data2 = gr.GeoRaster(
        data.raster[data.shape[0] / 2:, :],
        (xmin, xsize, x, ymax + ysize * data.shape[0] / 2, y, ysize),
        nodata_value=data.nodata_value,
        projection=data.projection,
        datatype=data.datatype,
    )

    # Plot both parts and save them
    plt.figure(figsize=(12, 8))
    data1.plot()
    plt.savefig(os.path.join(DATA, 'data1.png'), bbox_inches='tight')

.. image :: ./tests/data/data1.png
    
.. code-block:: python

    plt.figure(figsize=(12,8))
    data2.plot()
    plt.savefig(os.path.join(DATA,'data2.png'), bbox_inches='tight')
    
.. image :: ./tests/data/data2.png
    
.. code-block:: python

    # Generate merged raster
    
    data3 = data1.union(data2)
    
    # Plot it and save the figure
    
    plt.figure(figsize=(12,8))
    data3.plot()
    plt.savefig(os.path.join(DATA,'data3.png'), bbox_inches='tight')
    
.. image :: ./tests/data/data3.png
    

Another Merge:
--------------


Example Usage: Other functions
------------------------------

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

.. |CoverageStatus| image:: https://img.shields.io/coveralls/ozak/georasters.svg
.. _CoverageStatus: https://coveralls.io/r/ozak/georasters

.. |PyPiVersion| image:: https://img.shields.io/pypi/v/georasters.svg
.. _PyPiVersion: :target: https://pypi.python.org/pypi/georasters/
    :alt: Version on Pypi

.. |PyPiDownloads| image:: https://img.shields.io/pypi/dm/georasters.svg
.. _PyPiDownloads:     :target: https://pypi.python.org/pypi/georasters/
     :alt: Pypi downloads
