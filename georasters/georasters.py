#!/usr/bin/env python
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
import numpy as np
from osgeo import gdal, gdalnumeric, ogr, osr, gdal_array
from gdalconst import *
from skimage.measure import block_reduce
from skimage.transform import resize
import skimage.graph as graph
import matplotlib.pyplot as plt
import pandas as pd
from fiona.crs import from_string
import geopandas as gp

# Function to read the original file's projection:
def get_geo_info(FileName):
    ''' Gets information from a Raster data set
    '''
    SourceDS = gdal.Open(FileName, GA_ReadOnly)
    NDV = SourceDS.GetRasterBand(1).GetNoDataValue()
    xsize = SourceDS.RasterXSize
    ysize = SourceDS.RasterYSize
    GeoT = SourceDS.GetGeoTransform()
    Projection = osr.SpatialReference()
    Projection.ImportFromWkt(SourceDS.GetProjectionRef())
    DataType = SourceDS.GetRasterBand(1).DataType
    DataType = gdal.GetDataTypeName(DataType)
    return NDV, xsize, ysize, GeoT, Projection, DataType

# Function to map location in pixel of raster array
def map_pixel(point_x, point_y, cellx, celly, xmin, ymax):
    '''
    Usage: map_pixel(xcoord, ycoord, x_cell_size, y_cell_size, xmin, ymax)
    where: 
            xmin is leftmost X coordinate in system
            ymax is topmost Y coordinate in system
    Example:
            raster = HMISea.tif'
            NDV, xsize, ysize, GeoT, Projection, DataType = GetGeoInfo(raster)
            row, col = map_pixel(x,y,GeoT[1],GeoT[-1], GeoT[0],GeoT[3])
    '''
    point_x=np.asarray(point_x)
    point_y=np.asarray(point_y)
    col = np.floor((point_x - xmin) / cellx).astype(int)
    row = np.floor((point_y - ymax) / celly).astype(int)
    return row,col

def map_pixel_inv(row, col, cellx, celly, xmin, ymax):
    '''
    Usage: map_pixel(xcoord, ycoord, x_cell_size, y_cell_size, xmin, ymax)
    where: 
            xmin is leftmost X coordinate in system
            ymax is topmost Y coordinate in system
    Example:
            raster = HMISea.tif'
            NDV, xsize, ysize, GeoT, Projection, DataType = GetGeoInfo(raster)
            row, col = map_pixel(x,y,GeoT[1],GeoT[-1], GeoT[0],GeoT[3])
    '''
    col=np.asarray(col)
    row=np.asarray(row)
    point_x = xmin+col*cellx 
    point_y = ymax+row*celly
    return point_x, point_y

# Aggregate raster to higher resolution using sums
def aggregate(raster,NDV,block_size):
    '''
    Aggregate raster to smaller resolution, by adding cells.
    Usage:
            aggregate(raster,NDV,block_size)
    where
            raster is a Numpy array created by importing the raster (e.g. GeoTiff)
            NDV is the NoData Value for the raster (can be read using the GetGeoInfo function)
            block_size is a duple of factors by which the raster will be shrinked
    Example:
            raster = HMISea.tif'
            NDV, xsize, ysize, GeoT, Projection, DataType = GetGeoInfo(raster)
            costs = load_tiff(raster)
            costs2=aggregate(costs,NDV,(10,10))
    '''
    raster2=block_reduce(raster,block_size,func=np.ma.sum)
    return raster2

# Function to write a new file.
def create_geotiff(Name, Array, driver, NDV, xsize, ysize, GeoT, Projection, DataType):
    '''
    Creates new GeoTiff from array
    '''
    if type(DataType)!=np.int:
        if DataType.startswith('gdal.GDT_')==False:
            DataType=eval('gdal.GDT_'+DataType)
    NewFileName = Name+'.tif'
    # Set nans to the original No Data Value
    Array[np.isnan(Array)] = NDV
    # Set up the dataset
    DataSet = driver.Create( NewFileName, xsize, ysize, 1, DataType)
    # the '1' is for band 1.
    DataSet.SetGeoTransform(GeoT)
    DataSet.SetProjection( Projection.ExportToWkt() )
    # Write the array
    DataSet.GetRasterBand(1).WriteArray( Array )
    DataSet.GetRasterBand(1).SetNoDataValue(NDV)
    return NewFileName

# Function to aggregate and align rasters
def align_rasters(raster,alignraster,how=np.ma.mean,cxsize=None,cysize=None,masked=False):
    '''
    Align two rasters so that data overlaps by geographical location
    Usage: (alignedraster_o, alignedraster_a, GeoT_a) = AlignRasters(raster, alignraster, how=np.mean)
    where 
        raster: string with location of raster to be aligned
        alignraster: string with location of raster to which raster will be aligned
        how: function used to aggregate cells (if the rasters have different sizes)
    It is assumed that both rasters have the same size
    '''
    NDV1, xsize1, ysize1, GeoT1, Projection1, DataType1=GetGeoInfo(raster)
    NDV2, xsize2, ysize2, GeoT2, Projection2, DataType2=GetGeoInfo(alignraster)
    if Projection1.ExportToMICoordSys()==Projection2.ExportToMICoordSys():
        blocksize=(np.round(GeoT2[1]/GeoT1[1]),np.round(GeoT2[-1]/GeoT1[-1]))
        mraster=gdalnumeric.LoadFile(raster)
        mraster=np.ma.masked_array(mraster, mask=mraster==NDV1, fill_value=NDV1)
        mmin=mraster.min()
        mraster=block_reduce(mraster,blocksize,func=how)
        araster=gdalnumeric.LoadFile(alignraster)
        araster=np.ma.masked_array(araster, mask=araster==NDV2, fill_value=NDV2)
        amin=araster.min()
        if GeoT1[0]<=GeoT2[0]:
            row3,mcol=map_pixel(GeoT2[0], GeoT2[3], GeoT1[1] *blocksize[0],GeoT1[-1]*blocksize[1], GeoT1[0], GeoT1[3])
            acol=0
        else:
            row3,acol=map_pixel(GeoT1[0], GeoT1[3], GeoT2[1],GeoT2[-1], GeoT2[0], GeoT2[3])
            mcol=0
        if GeoT1[3]<=GeoT2[3]:
            arow,col3=map_pixel(GeoT1[0], GeoT1[3], GeoT2[1],GeoT2[-1], GeoT2[0], GeoT2[3])
            mrow=0
        else:
            mrow,col3=map_pixel(GeoT2[0], GeoT2[3], GeoT1[1] *blocksize[0],GeoT1[-1]*blocksize[1], GeoT1[0], GeoT1[3])
            arow=0
        '''
        col3,row3=map_pixel(GeoT1[0], GeoT1[3], GeoT2[1],GeoT2[-1], GeoT2[0], GeoT2[3])
        col3=max(0,col3)
        row3=max(0,row3)
        araster=araster[row3:,col3:]
        col3,row3=map_pixel(GeoT2[0], GeoT2[3], GeoT1[1] *blocksize[0],GeoT1[-1]*blocksize[1], GeoT1[0], GeoT1[3])
        col3=max(0,abs(col3))
        row3=max(0,np.abs(row3))
        mraster=mraster[row3:,col3:]
        '''
        mraster=mraster[mrow:,mcol:]
        araster=araster[arow:,acol:]
        if cxsize and cysize:
            araster=araster[:cysize,:cxsize]
            mraster=mraster[:cysize,:cxsize]
        else:
            rows = min(araster.shape[0],mraster.shape[0])
            cols = min(araster.shape[1],mraster.shape[1])
            araster=araster[:rows,:cols]
            mraster=mraster[:rows,:cols]
        #mraster=mraster[row3:rows+row3,col3:cols+col3]
        if masked:
            mraster=np.ma.masked_array(mraster,mask=mraster<mmin, fill_value=NDV1)
            araster=np.ma.masked_array(araster,mask=araster<amin, fill_value=NDV2)
        GeoT=(max(GeoT1[0],GeoT2[0]), GeoT1[1]*blocksize[0], GeoT1[2], min(GeoT1[3],GeoT2[3]), GeoT1[4] ,GeoT1[-1]*blocksize[1])
        return (mraster,araster,GeoT)
    else:
        print("Rasters need to be in same projection")
        return (-1,-1,-1)

# Load GeoTif raster data
def load_tiff(file):
    """
    Load a GeoTiff raster keeping NDV values using a masked array
    Usage:
            data=LoadTiffRaster(file)
    """
    NDV, xsize, ysize, GeoT, Projection, DataType=get_geo_info(file)
    data=gdalnumeric.LoadFile(file)
    data=np.ma.masked_array(data, mask=data==NDV,fill_value=NDV)
    return data

class RasterGeoTError(Exception):
    pass

class RasterGeoError(Exception):
    pass

class RasterGeoTWarning(Exception):
    pass


# GeoRaster Class
class GeoRaster(object):
    '''
    GeoRaster class to create and handle GIS rasters
    Eash GeoRaster object is a numpy masked array + geotransfrom + nodata_value
    Usage:
        geo=GeoRaster(raster, geot, nodata_value = ndv)
    where
        raster: Numpy masked array with the raster data, which could be loaded with the load_tiff(file)
        geot: GDAL Geotransformation
        nodata_value: No data value in raster, optional
    '''
    def __init__(self, raster, geot, nodata_value = np.nan, projection = None, datatype=None):
        '''
        Initialize Georaster
        Usage:
            geo=GeoRaster(raster, geot, nodata_value = ndv)
        where
            raster: Numpy masked array with the raster data, which could be loaded with from_file(file) or load_tiff(file)
            geot: GDAL Geotransformation
            nodata_value: No data value in raster, optional
        '''
        super(GeoRaster, self).__init__()
        if isinstance(raster, np.ma.core.MaskedArray):
            self.raster = raster
        else:
            self.raster = np.ma.masked_array(raster, mask= raster==nodata_value, fill_value=nodata_value)
        self.geot = geot
        self.nodata_value = nodata_value
        self.shape = raster.shape
        self.x_cell_size = geot[1]
        self.y_cell_size = geot[-1]
        self.xmin = geot[0]
        self.ymax = geot[3]
        self.xmax = self.xmin + self.x_cell_size * self.shape[1]
        self.ymin = self.ymax + self.y_cell_size * self.shape[0]
        self.bounds = (self.xmin, self.ymin, self.xmax, self.ymax)
        self.projection = projection
        self.datatype=datatype
        self.mcp_cost=None


    def __getitem__(self,indx):
        rast = self.raster.__getitem__(indx)
        proj = self.projection
        nodata = self.nodata_value
        datatype = self.datatype
        geot = list(self.geot)
        geot[0] += indx[0].start*geot[1]
        geot[3] += indx[1].start*geot[-1]
        geot = tuple(geot)
        return GeoRaster(rast,geot,nodata,proj,datatype)
    
    def __getslice__(self,i,j):
        return self.raster.__getslice__(i,j)

#    def __getattribute__(self, attr):
#        return eval('self.'+attr)

    def __lt__(self,other):
        if isinstance(other, GeoRaster):
            return self.raster<other.raster
        elif isinstance(other, np.ndarray):
            return self.raster<other
        else:
            return self.raster<other

    def __le__(self,other):
        if isinstance(other, GeoRaster):
            return self.raster<=other.raster
        elif isinstance(other, np.ndarray):
            return self.raster<=other
        else:
            return self.raster<=other

    def __gt__(self,other):
        if isinstance(other, GeoRaster):
            return self.raster>other.raster
        elif isinstance(other, np.ndarray):
            return self.raster>other
        else:
            return self.raster>other

    def __ge__(self,other):
        if isinstance(other, GeoRaster):
            return self.raster>=other.raster
        elif isinstance(other, np.ndarray):
            return self.raster>=other
        else:
            return self.raster>=other

    def __eq__(self,other):
        if isinstance(other, GeoRaster):
            return self.raster==other.raster
        elif isinstance(other, np.ndarray):
            return self.raster==other
        else:
            return self.raster==other

    def __ne__(self,other):
        if isinstance(other, GeoRaster):
            return self.raster!=other.raster
        elif isinstance(other, np.ndarray):
            return self.raster!=other
        else:
            return self.raster!=other

    def __pos__(self):
        return self

    def __neg__(self):
        return GeoRaster(-self.raster, self.geot, nodata_value=self.nodata_value, projection = self.projection, datatype = self.datatype)

    def __add__(self, other):
        if isinstance(other,GeoRaster):
            if self.geot!=other.geot:
                raise RasterGeoTWarning("Rasters do not have same geotransform. If needed first create union or allign them.")
            if self.nodata_value==other.nodata_value:
                ndv=self.nodata_value
            else:
                ndv=np.nan
            return GeoRaster(self.raster+other.raster, self.geot, nodata_value=ndv, projection = self.projection, datatype = self.datatype)
        else:
            return GeoRaster(self.raster+other, self.geot, nodata_value=self.nodata_value, projection = self.projection, datatype = self.datatype)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self+other.__neg__()

    def __rsub__(self, other):
        return self.__sub__(other)

    def __mul__(self, other):
        if isinstance(other,GeoRaster):
            if self.geot!=other.geot:
                raise RasterGeoTWarning("Rasters do not have same geotransform. If needed first create union or allign them.")
            if self.nodata_value==other.nodata_value:
                ndv=self.nodata_value
            else:
                ndv=np.nan
            return GeoRaster(self.raster*other.raster, self.geot, nodata_value=ndv, projection = self.projection, datatype = self.datatype)
        else:
            return GeoRaster(self.raster*other, self.geot, nodata_value=self.nodata_value, projection = self.projection, datatype = self.datatype)

    def __rmul__(self, other):
        return self.__mul__(other) 

    def __truediv__(self, other):
        if isinstance(other,GeoRaster):
            if self.geot!=other.geot:
                raise RasterGeoTWarning("Rasters do not have same geotransform. If needed first create union or allign them.")
            if self.nodata_value==other.nodata_value:
                ndv=self.nodata_value
            else:
                ndv=np.nan
            return GeoRaster(self.raster/other.raster, self.geot, nodata_value=ndv, projection = self.projection, datatype = self.datatype)
        else:
            return GeoRaster(self.raster/other, self.geot, nodata_value=self.nodata_value, projection = self.projection, datatype = self.datatype)

    def __rtruediv__(self, other):
        if isinstance(other,GeoRaster):
            return other.__truediv__(self)
        else:
            return GeoRaster(other/self.raster, self.geot, nodata_value=self.nodata_value, projection = self.projection, datatype = self.datatype)

    def __floordiv__(self, other):
        A=self/other
        A.raster=A.raster.astype(int)
        return A

    def __rfloordiv__(self, other):
        if isinstance(other,GeoRaster):
            if self.geot!=other.geot:
                raise RasterGeoTWarning("Rasters do not have same geotransform. If needed first create union or allign them.")
            if self.nodata_value==other.nodata_value:
                ndv=self.nodata_value
            else:
                ndv=np.nan
            return GeoRaster((other.raster/self.raster).astype(int), self.geot, nodata_value=ndv, projection = self.projection, datatype = self.datatype)
        else:
            return GeoRaster((other/self.raster).astype(int), self.geot, nodata_value=self.nodata_value, projection = self.projection, datatype = self.datatype)

    def __pow__(self,other):
        if isinstance(other,GeoRaster):
            if self.geot!=other.geot:
                raise RasterGeoTWarning("Rasters do not have same geotransform. If needed first create union or allign them.")
            if self.nodata_value==other.nodata_value:
                ndv=self.nodata_value
            else:
                ndv=np.nan
            return GeoRaster(self.raster**other.raster, self.geot, nodata_value=ndv, projection = self.projection, datatype = self.datatype)
        else:
            return GeoRaster(self.raster**other, self.geot, nodata_value=self.nodata_value, projection = self.projection, datatype = self.datatype)

    def copy(self):
        """Returns copy of itself"""
        return GeoRaster(self.raster.copy(), self.geot, nodata_value = self.nodata_value, projection = self.projection, datatype=self.datatype)

    def to_tiff(self, filename):
        '''
        geo.to_tiff(filename)
        
        Saves GeoRaster as GeoTiff filename.tif with type datatype
        
        If GeoRaster does not have datatype, then it tries to assign a type.
        You can assign the type yourself by setting
         geo.datatype = 'gdal.GDT_'+type
        '''
        if self.datatype is None:
            self.datatype = gdal_array.NumericTypeCodeToGDALTypeCode(self.raster.data.dtype)
            if self.datatype is None:
                if self.raster.data.dtype.name.find('int')!=-1:
                    self.raster = self.raster.astype(np.int32)
                    self.datatype = gdal_array.NumericTypeCodeToGDALTypeCode(self.raster.data.dtype)
                else:
                    self.raster = self.raster.astype(np.float64)
                    self.datatype = gdal_array.NumericTypeCodeToGDALTypeCode(self.raster.data.dtype)
        self.raster.data[self.raster.mask] = self.nodata_value
        create_geotiff(filename, self.raster, gdal.GetDriverByName('GTiff'), self.nodata_value, self.shape[1], self.shape[0], self.geot, self.projection, self.datatype)

    def to_pandas(self):
        """
        Convert GeoRaster to Pandas DataFrame, which can be easily exported to other types of files
        The DataFrame has the row, col, value, x, and y values for each cell
        """
        df = to_pandas(self)
        return df

    def to_geopandas(self):
        """
        Convert GeoRaster to GeoPandas DataFrame, which can be easily exported to other types of files and used to 
        do other types of operations.
        The DataFrame has the geometry (Polygon), row, col, value, x, and y values for each cell
        """
        df = to_geopandas(self)
        return df

    def plot(self):
        '''
        geo.plot()

        Returns plot of raster data
        '''
        plt.matshow(self.raster)

    def union(self, other):
        '''
        geo.union(Georaster)

        Returns union of GeoRaster with another one
        '''
        return union([self,other])

    def merge(self, other):
        '''
        geo.merge(Georaster)

        Returns merge of GeoRaster with another one
        '''
        return merge([self,other])

    def mean(self, *args, **kwargs):
        '''
        geo.mean(axis=None, dtype=None, out=None)

        Returns the average of the array elements along given axis.

        Refer to `numpy.mean` for full documentation.

        See Also
        --------
        numpy.mean : equivalent function
        '''
        return self.raster.mean(*args, **kwargs)

    def max(self, *args, **kwargs):
        '''
        geo.max(axis=None, out=None)

        Return the maximum along a given axis.

        Refer to `numpy.amax` for full documentation.

        See Also
        --------
        numpy.amax : equivalent function
        '''
        return self.raster.max(*args, **kwargs)

    def min(self, *args, **kwargs):
        '''
        geo.min(axis=None, out=None)

        Return the minimum along a given axis.

        Refer to `numpy.amin` for full documentation.

        See Also
        --------
        numpy.amin : equivalent function
        '''
        return self.raster.min(*args, **kwargs)

    def median(self, *args, **kwargs):
        '''
        geo.median(axis=None, out=None, overwrite_input=False)
        
        axis : int, optional
            Axis along which the medians are computed. The default (axis=None)
            is to compute the median along a flattened version of the array.
        out : ndarray, optional
            Alternative output array in which to place the result. It must have
            the same shape and buffer length as the expected output, but the
            type (of the output) will be cast if necessary.
        overwrite_input : bool, optional
           If True, then allow use of memory of input array (a) for
           calculations. The input array will be modified by the call to
           median. This will save memory when you do not need to preserve the
           contents of the input array. Treat the input as undefined, but it
           will probably be fully or partially sorted. Default is False. Note
           that, if `overwrite_input` is True and the input is not already an
           ndarray, an error will be raised.
        '''
        return np.ma.median(self.raster, *args, **kwargs)

    def std(self, *args, **kwargs):
        '''
        geo.std(axis=None, dtype=None, out=None, ddof=0)

        Returns the standard deviation of the array elements along given axis.

        Refer to `numpy.std` for full documentation.

        See Also
        --------
        numpy.std : equivalent function
        '''
        return self.raster.std(*args, **kwargs)

    def argmax(self, *args, **kwargs):
        '''
        geo.argmax(axis=None, out=None)

        Return indices of the maximum values along the given axis.

        Refer to `numpy.argmax` for full documentation.

        See Also
        --------
        numpy.argmax : equivalent function
        '''
        return self.raster.argmax(*args, **kwargs)

    def argmin(self, *args, **kwargs):
        '''
        geo.argmin(axis=None, out=None)

        Return indices of the minimum values along the given axis of `a`.

        Refer to `numpy.argmin` for detailed documentation.

        See Also
        --------
        numpy.argmin : equivalent function
        '''
        return self.raster.argmin(*args, **kwargs)

    def sum(self, *args, **kwargs):
        '''
        geo.sum(axis=None, dtype=None, out=None)

        Return the sum of the array elements over the given axis.

        Refer to `numpy.sum` for full documentation.

        See Also
        --------
        numpy.sum : equivalent function
        '''
        return self.raster.sum(*args, **kwargs)

    def prod(self, *args, **kwargs):
        '''
        geo.prod(axis=None, dtype=None, out=None)

        Return the product of the array elements over the given axis

        Refer to `numpy.prod` for full documentation.

        See Also
        --------
        numpy.prod : equivalent function
        '''
        return self.raster.prod(*args, **kwargs)

    def var(self, *args, **kwargs):
        '''
        geo.var(axis=None, dtype=None, out=None, ddof=0)

        Returns the variance of the array elements, along given axis.

        Refer to `numpy.var` for full documentation.

        See Also
        --------
        numpy.var : equivalent function
        '''
        return self.raster.var(*args, **kwargs)

    def count(self, *args, **kwargs):
        '''
        geo.count(axis=None)
        Count the non-masked elements of the array along the given axis.
        '''
        return self.raster.count(*args, **kwargs)

    def gini(self):
        """
        geo.gini()
        
        Return computed Gini coefficient.
        """
        if self.count()>1:
            xsort = sorted(self.raster.data[self.raster.mask==False].flatten()) # increasing order
            y = np.cumsum(xsort)
            B = sum(y) / (y[-1] * len(xsort))
            return 1 + 1./len(xsort) - 2*B
        else:
            return 1

    def flatten(self, *args, **kwargs):
        '''
        geo.flatten(order='C')

        Return a copy of the array collapsed into one dimension.

        Parameters
        ----------
        order : {'C', 'F', 'A'}, optional
            Whether to flatten in C (row-major), Fortran (column-major) order,
            or preserve the C/Fortran ordering from `a`.
            The default is 'C'.
        '''
        return self.raster.flatten(*args, **kwargs)

    def apply(self, func, *args, **kwargs):
        '''
        geo.apply(func, *args, **kwargs)
        
        Returns the value of applying function func on the raster data
        
        func: Python function
        *args: Arguments of function
        **kwargs: Additional arguments of function
        '''
        return func(self.raster, *args, **kwargs)

    def map_pixel(self, point_x, point_y):
        '''
        geo.map_pixel(point_x, point_y)
        
        Return value of raster in location 
        '''
        row, col =map_pixel(point_x, point_y, self.x_cell_size, self.y_cell_size, self.xmin, self.ymax)
        return self.raster[row, col]

    def map_pixel_location(self, point_x, point_y):
        '''
        geo.map_pixel(point_x, point_y)
        
        Return value of raster in location 
        '''
        row, col =map_pixel(point_x, point_y, self.x_cell_size, self.y_cell_size, self.xmin, self.ymax)
        return np.array([row, col])

    def extract(self, point_x, point_y, radius=0):
        '''
        geo.extract(x, y, radius=r)
        
        Return subraster of raster geo around location (x,y) with radius r
        where (x,y) and r are in the same coordinate system as geo
        '''
        row, col =map_pixel(point_x, point_y, self.x_cell_size, self.y_cell_size, self.xmin, self.ymax)
        col2 = np.abs(radius/self.x_cell_size)
        row2 = np.abs(radius/self.y_cell_size)
        return GeoRaster(self.raster[max(row-row2, 0):min(row+row2+1, self.shape[0]), \
                        max(col-col2, 0):min(col+col2+1, self.shape[1])], self.geot, nodata_value = self.nodata_value,\
                        projection=self.projection, datatype = self.datatype)

    # Align GeoRasters
    def align(self,alignraster,how=np.mean,cxsize=None,cysize=None):
        '''
        geo.align(geo2, how=np.mean)

        Returns both georasters aligned and with the same pixelsize
        '''
        return align_georasters(self,alignraster,how=how,cxsize=cxsize,cysize=cysize)

    def aggregate(self,block_size):
        '''
        geo.aggregate(block_size)

        Returns copy of raster aggregated to smaller resolution, by adding cells.
        '''
        raster2=block_reduce(self.raster,block_size,func=np.ma.sum)
        geot = self.geot
        geot = (geot[0], block_size[0] * geot[1], geot[2], geot[3],geot[4], block_size[1] * geot[-1])
        return GeoRaster(raster2, geot, nodata_value=self.nodata_value,\
                        projection=self.projection, datatype = self.datatype)

    def block_reduce(self,block_size, how=np.ma.mean):
        '''
        geo.block_reduce(block_size, how=func)

        Returns copy of raster aggregated to smaller resolution, by adding cells.
        Default: func=np.ma.mean
        '''
        raster2=block_reduce(self.raster,block_size,func=how)
        return GeoRaster(raster2, self.geot, nodata_value=self.nodata_value,\
                        projection=self.projection, datatype = self.datatype)

    def resize(self, block_size, order=0, mode='constant', cval=False, preserve_range=True):
        '''
        geo.resize(new_shape, order=0, mode='constant', cval=np.nan, preserve_range=True)
        
        Returns resized georaster
        '''
        if not cval:
            cval=np.nan
        raster2=resize(self.raster.data, block_size, order=order, mode=mode, cval=cval, preserve_range=preserve_range)
        mask=resize(self.raster.mask, block_size, order=order, mode=mode, cval=cval, preserve_range=preserve_range)
        raster2=np.ma.masked_array(raster2, mask=mask, fill_value=self.raster.fill_value)
        raster2[raster2.mask]=self.nodata_value
        raster2.mask=np.logical_or(np.isnan(raster2.data),raster2.data==self.nodata_value)
        geot=list(self.geot)
        [geot[-1],geot[1]]=np.array([geot[-1],geot[1]])*self.shape/block_size
        return GeoRaster(raster2, tuple(geot), nodata_value=self.nodata_value,\
                        projection=self.projection, datatype = self.datatype)

    def resize_old(self, block_size, order=0, mode='constant', cval=False):
        '''
        geo.resize(new_shape, order=0, mode='constant', cval=np.nan, preserve_range=True)
        
        Returns resized georaster
        '''
        if not cval:
            cval=np.nan
        if self.raster.dtype.name.find('float')!=-1 and np.max(np.abs([self.max(),self.min()]))>1:
            raster2=(self.raster-self.min())/(self.max()-self.min())
        else:
            raster2=self.raster.copy()
        raster2=raster2.astype(float)
        raster2[self.raster.mask]=np.nan
        raster2=resize(raster2,block_size, order=order, mode=mode, cval=cval)
        raster2=np.ma.masked_array(raster2, mask=np.isnan(raster2), fill_value=self.raster.fill_value)
        raster2=raster2*(self.max()-self.min())+self.min()
        raster2[raster2.mask]=self.nodata_value
        raster2.mask=np.logical_or(np.isnan(raster2.data),raster2.data==self.nodata_value)
        geot=list(self.geot)
        [geot[-1],geot[1]]=np.array([geot[-1],geot[1]])*self.shape/block_size
        return GeoRaster(raster2, tuple(geot), nodata_value=self.nodata_value,\
                        projection=self.projection, datatype = self.datatype)

    # Setup Graph for distance computations and provide distance functions
    def mcp(self, *args, **kwargs):
        """
        Setup MCP_Geometric object from skimage for optimal travel time computations
        """
        # Create Cost surface to work on
        self.mcp_cost=graph.MCP_Geometric(self.raster, *args, **kwargs)
    pass

    # Determine minimum travel cost to each location 
    def distance(self, sources, destinations, x='x', y='y', isolation=True, export_raster=False, export_shape=False, routes=False, path='./'):
        """
        Compute cost distance measured from each start point to all end points.
        The function returns the distances between the start point and the end points as a Pandas dataframe.
        Additionally, for each start point it computes the level of isolation, i.e. its average travel distance to 
        all other locations
        """
        start_points=sources.copy()
        end_points=destinations.copy()
        if not isinstance(start_points,pd.core.frame.DataFrame) and not isinstance(start_points,gp.geodataframe.GeoDataFrame):
            raise TypeError('Sources has to be a (Geo)Pandas Data Frame Object.')
        if not isinstance(end_points,pd.core.frame.DataFrame) and not isinstance(end_points,gp.geodataframe.GeoDataFrame):
            raise TypeError('Destinations has to be a (Geo)Pandas Data Frame Object.')
        if not self.mcp_cost:
            self.mcp()
        count=0
        start_points['row'], start_points['col'] = self.map_pixel_location(start_points[x], start_points[y])
        end_points['row'], end_points['col'] = self.map_pixel_location(end_points[x], end_points[y])
        start_points['ID']=start_points.index.values
        end_points['ID']=end_points.index.values+start_points['ID'].max()+1
        
        for i in start_points.iterrows():
            cumulative_costs, traceback = self.mcp_cost.find_costs([[i[1].row,i[1].col]])
            dist=cumulative_costs[end_points.row.values,end_points.col.values].transpose()/(7*24)
            df2=pd.DataFrame(np.array([(i[1]['ID']*np.ones_like(dist)).flatten(),
                                    end_points['ID'],dist.flatten()]).transpose(), 
                                    columns=['ID1','ID2','dist'])
            # Keep only locations that are accessible
            df2 = df2.loc[df2['dist']<np.inf]
            if isolation:
                grisolation = np.ma.masked_array(cumulative_costs, mask = np.logical_or(self.raster.mask, cumulative_costs==np.inf)
                                    , fill_value=np.nan).mean()/(7*24)
                start_points.loc[i[0],'Iso'] = grisolation
            if export_raster:
                cumulative_costs = gr.GeoRaster(np.ma.masked_array(cumulative_costs, 
                                                    mask = np.logical_or(self.raster.mask, cumulative_costs==np.inf), 
                                                    fill_value = np.nan), self.geot, self.nodata_value, 
                                                    projection=self.projection, datatype = self.datatype)
                cumulative_costs.raster.data[cumulative_costs.raster.mask] = cumulative_costs.nodata_value
                cumulative_costs.to_tiff(path+str(i[1]['ID']))
            if df2.size>0:
                if export_shape:
                    routes=True
                if routes:
                    df2['geometry'] = df2['ID2'].apply(lambda x: 
                                                self.mcp_cost.traceback(end_points.loc[end_points['ID']==x][['row','col']].values[0]))
                    df2['geometry'] = df2.geometry.apply(lambda x: [gr.map_pixel_inv(y[0],y[1], self.geot[1], 
                                                    self.geot[-1], self.geot[0], self.geot[-3]) for y in x])
                    df2['geometry'] = df2.geometry.apply(lambda x: LineString(x) if int(len(x)>1) else LineString([x[0],x[0]]))
                    df2 = gp.GeoDataFrame(df2, crs = cea)
                if isolation:
                    df2['Iso'] = grisolation
                if count==0:
                    self.grdist=df2.copy()
                else:
                    self.grdist=self.grdist.append(df2)
                count+=1
        if routes:
            self.grdist = gp.GeoDataFrame(self.grdist, crs = cea)
        if export_shape:
            start_pointscols = sources.columns.values
            end_pointscols = destinations.columns.values
            if 'geometry' in end_pointscols:
                self.grdist = pd.merge(self.grdist, end_points[['ID']+end_pointscols.tolist()].drop('geometry', axis=1), left_on='ID2', right_on='ID', how='left')
            else:
                self.grdist = pd.merge(self.grdist, end_points[['ID']+end_pointscols.tolist()], left_on='ID2', right_on='ID', how='left')
            if 'geometry' in self.start_pointscols:
                self.grdist = pd.merge(self.grdist, start_points[['ID']+start_pointscols.tolist()].drop('geometry', axis=1), left_on='ID1', right_on='ID', how='left',
                             suffixes = ['_2','_1'])
            else:
                self.grdist = pd.merge(self.grdist, start_points[['ID']+start_pointscols.tolist()], left_on='ID1', right_on='ID', how='left',
                             suffixes = ['_2','_1'])
            self.grdist = gp.GeoDataFrame(self.grdist, crs = cea)
            self.grdist.to_file(path+'routes.shp')
    pass

##########################################
# Useful functions defined on GeoRasters
##########################################

# Union of rasters
def union(rasters):
    """
    Union of rasters
    Usage:
        union(rasters)
    where
        rasters is a list of GeoRaster objects
    """
    if sum([rasters[0].x_cell_size==i.x_cell_size for i in rasters])==len(rasters) \
       and sum([rasters[0].y_cell_size==i.y_cell_size for i in rasters])==len(rasters)\
       and sum([rasters[0].projection.ExportToProj4()==i.projection.ExportToProj4() for i in rasters])==len(rasters):
        if sum([rasters[0].nodata_value==i.nodata_value for i in rasters])==len(rasters):
            ndv=rasters[0].nodata_value
        else:
            ndv = np.nan
        if ndv==None:
            ndv = np.nan
        if sum([rasters[0].datatype==i.datatype for i in rasters])==len(rasters):
            datatype=rasters[0].datatype
        else:
            datatype = None
        projection = rasters[0].projection
        lonmin = min([i.xmin for i in rasters])
        lonmax = max([i.xmax for i in rasters])
        latmin = min([i.ymin for i in rasters])
        latmax = max([i.ymax for i in rasters])
        shape = (np.abs(np.round((latmax-latmin)/rasters[0].y_cell_size)).astype(int),np.round((lonmax-lonmin)/rasters[0].x_cell_size).astype(int))
        out = ndv*np.ones(shape)
        outmask = np.ones(shape).astype(bool)
        for i in rasters:
            (row,col) = map_pixel(i.xmin, i.ymax, rasters[0].x_cell_size, rasters[0].y_cell_size, lonmin, latmax)
            out[row:row+i.shape[0],col:col+i.shape[1]] = np.where(i.raster.data!=i.nodata_value, i.raster.data,\
                                                         out[row:row+i.shape[0],col:col+i.shape[1]])#i.raster
            outmask[row:row+i.shape[0],col:col+i.shape[1]] = np.where(i.raster.mask==False, False,\
                                                         outmask[row:row+i.shape[0],col:col+i.shape[1]])#i.raster
        out=np.ma.masked_array(out, mask=outmask, fill_value=ndv)
        return GeoRaster(out, (lonmin, rasters[0].x_cell_size, 0.0, latmax, 0.0, rasters[0].y_cell_size), nodata_value=ndv, projection=projection, datatype=datatype)
    else:
        raise RasterGeoError('Rasters need to have same pixel sizes. Use the aggregate or dissolve functions to generate correct GeoRasters')

def merge(rasters):
    """
    Merges GeoRasters, same as union
    Usage:
        merge(rasters)
    where
        rasters is a list of GeoRaster objects
    """
    return union(rasters)

# Load data into a GeoRaster from file
def from_file(filename):
    """
    Create a GeoRaster object from a file
    """
    NDV, xsize, ysize, GeoT, Projection, DataType = get_geo_info(filename)
    data = gdalnumeric.LoadFile(filename)
    data = np.ma.masked_array(data, mask=data==NDV,fill_value=NDV)
    return GeoRaster(data,GeoT, nodata_value=NDV, projection=Projection, datatype=DataType)

# Convert Pandas DataFrame to raster
def from_pandas(df, value='value', x='x', y='y', cellx= None, celly=None, xmin=None, ymax=None,
                GeoT=None, nodata_value=None, projection=None, datatype=None):
    """
    Creates a GeoRaster from a Pandas DataFrame. Useful to plot or export data to rasters.
    Usage:
        raster = from_pandas(df, value='value', x='x', y='y', cellx= cellx, celly=celly, xmin=xmin, ymax=ymax,
                        GeoT=GeoT, nodata_value=NDV, projection=Projection, datatype=DataType)

    Although it does not require all the inputs, it is highly recommended to include the geographical information, 
    so that the GeoRaster is properly constructed. As usual, the information can be added afterwards directly to the 
    GeoRaster.
    """
    if not cellx:
        cellx = (df.sort(x)[x]-df.sort(x).shift(1)[x]).max()
    if not celly:
        celly = (df.sort(y, ascending=False)[y]-df.sort(y, ascending=False).shift(1)[y]).min()
    if not xmin:
        xmin = df[x].min()
    if not ymax:
        ymax = df[y].max()
    row,col = map_pixel(df[x], df[y], cellx, celly, xmin, ymax)
    dfout = pd.DataFrame(np.array([row, col,df.value]).T, columns=['row','col','value'])
    dfout = dfout.set_index(['row','col']).unstack().values
    if nodata_value:
        dfout[np.isnan(dfout)]=nodata_value
    if not nodata_value:
        nodata_value=np.nan
    if not GeoT:
        GeoT = (xmin, cellx, 0, ymax, 0, celly)
    return GeoRaster(dfout, GeoT, nodata_value=nodata_value, projection=projection, datatype=datatype)

# GDAL conversion types
NP2GDAL_CONVERSION = {
  "uint8": 1,
  "int8": 1,
  "uint16": 2,
  "int16": 3,
  "uint32": 4,
  "int32": 5,
  "float32": 6,
  "float64": 7,
  "complex64": 10,
  "complex128": 11,
}

# Align GeoRasters
def align_georasters(raster,alignraster,how=np.mean,cxsize=None,cysize=None):
    '''
    Align two rasters so that data overlaps by geographical location
    Usage: (alignedraster_o, alignedraster_a) = AlignRasters(raster, alignraster, how=np.mean)
    where 
        raster: string with location of raster to be aligned
        alignraster: string with location of raster to which raster will be aligned
        how: function used to aggregate cells (if the rasters have different sizes)
    It is assumed that both rasters have the same size
    '''
    (NDV1, xsize1, ysize1, GeoT1, Projection1, DataType1)=(raster.nodata_value, raster.shape[1], raster.shape[0], raster.geot, raster.projection, raster.datatype)
    (NDV2, xsize2, ysize2, GeoT2, Projection2, DataType2)=(alignraster.nodata_value, alignraster.shape[1], alignraster.shape[0], alignraster.geot, alignraster.projection, alignraster.datatype)
    if Projection1.ExportToMICoordSys()==Projection2.ExportToMICoordSys():
        blocksize=(np.round(max(GeoT2[1]/GeoT1[1],1)),np.round(max(GeoT2[-1]/GeoT1[-1],1)))
        mraster=raster.raster
        mmin=mraster.min()
        if block_reduce!=(1,1):
            mraster=block_reduce(mraster,blocksize,func=how)
        blocksize=(np.round(max(GeoT1[1]/GeoT2[1],1)),np.round(max(GeoT1[-1]/GeoT2[-1],1)))
        araster=alignraster.raster
        amin=araster.min()
        if block_reduce!=(1,1):
            araster=block_reduce(araster,blocksize,func=how)
        if GeoT1[0]<=GeoT2[0]:
            row3,mcol=map_pixel(GeoT2[0], GeoT2[3], GeoT1[1] *blocksize[0],GeoT1[-1]*blocksize[1], GeoT1[0], GeoT1[3])
            acol=0
        else:
            row3,acol=map_pixel(GeoT1[0], GeoT1[3], GeoT2[1],GeoT2[-1], GeoT2[0], GeoT2[3])
            mcol=0
        if GeoT1[3]<=GeoT2[3]:
            arow,col3=map_pixel(GeoT1[0], GeoT1[3], GeoT2[1],GeoT2[-1], GeoT2[0], GeoT2[3])
            mrow=0
        else:
            mrow,col3=map_pixel(GeoT2[0], GeoT2[3], GeoT1[1] *blocksize[0],GeoT1[-1]*blocksize[1], GeoT1[0], GeoT1[3])
            arow=0
        mraster=mraster[mrow:,mcol:]
        araster=araster[arow:,acol:]
        if cxsize and cysize:
            araster=araster[:cysize,:cxsize]
            mraster=mraster[:cysize,:cxsize]
        else:
            rows = min(araster.shape[0],mraster.shape[0])
            cols = min(araster.shape[1],mraster.shape[1])
            araster=araster[:rows,:cols]
            mraster=mraster[:rows,:cols]
        mraster=np.ma.masked_array(mraster,mask=mraster<mmin, fill_value=NDV1)
        araster=np.ma.masked_array(araster,mask=araster<amin, fill_value=NDV2)
        GeoT=(max(GeoT1[0],GeoT2[0]), GeoT1[1]*blocksize[0], GeoT1[2], min(GeoT1[3],GeoT2[3]), GeoT1[4] ,GeoT1[-1]*blocksize[1])
        mraster=GeoRaster(mraster, GeoT, projection=Projection1, nodata_value=NDV1, datatype=DataType1)
        araster=GeoRaster(araster, GeoT, projection=Projection2, nodata_value=NDV2, datatype=DataType2)
        return (mraster,araster)
    else:
        print("Rasters need to be in same projection")
        return (-1,-1)

# Test if two GeoRasters are in same GeoT and have same projection
def is_geovalid(grasterlist):
    if np.sum(map(isinstance,grasterlist,[GeoRaster for i in range(len(grasterlist))]))==len(grasterlist):
        graster0 = grasterlist[-1]
        while grasterlist !=[]:
            grasterlist = grasterlist[:-1]
            graster1 = grasterlist[-1]
            if graster1.geot == graster0.geot and graster1.projection.ExportToPrettyWkt() == graster0.ExportToPrettyWkt():
                return 0
            else:
                raise RasterGeoTError("Rasters must have same geotransform and projection.")
                return 1
    else:
        raise RasterGeoTError("List must contain only GeoRasters.")

# Convert GeoRaster to Pandas DataFrame, which can be easily exported to other types of files
# Function to
def to_pandas(raster):
    """
    Convert GeoRaster to Pandas DataFrame, which can be easily exported to other types of files
    The DataFrame has the row, col, value, x, and y values for each cell
    Usage:
        df = gr.to_pandas(raster)
    """
    df = pd.DataFrame(raster.raster)
    df = df.stack()
    df = df.reset_index()
    df.columns = ['row','col','value']
    df['x'] = df.col.apply(lambda col: raster.geot[0]+(col)*raster.geot[1])
    df['y'] = df.row.apply(lambda row: raster.geot[3]+(row)*raster.geot[-1])
    return df

# Convert GeoRaster to GeoPandas
def squares(row, georaster=None):
    geometry = Polygon([(row.x,row.y),(row.x+georaster.x_cell_size,row.y),
                       (row.x+georaster.x_cell_size,row.y+georaster.y_cell_size),
                       (row.x,row.y+georaster.y_cell_size)])
    return geometry

def to_geopandas(raster):
    """
    Convert GeoRaster to GeoPandas DataFrame, which can be easily exported to other types of files and used to 
    do other types of operations.
    The DataFrame has the geometry (Polygon), row, col, value, x, and y values for each cell
    Usage:
        df = gr.to_geopandas(raster)
    """
    df = to_pandas(raster)
    df['geometry'] = df.apply(squares, georaster=raster, axis=1)
    df = gp.GeoDataFrame(df, crs=from_string(raster.projection.ExportToProj4()))
    return df
