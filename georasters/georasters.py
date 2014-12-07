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
import matplotlib.pyplot as plt
import pandas as pd

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
    point_x=np.array(point_x)
    point_y=np.array(point_y)
    col = np.around((point_x - xmin) / cellx).astype(int)
    row = np.around((point_y - ymax) / celly).astype(int)
    return row,col

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
    raster2=np.where(raster==NDV,0,raster)
    raster3=block_reduce(raster2,block_size,func=np.sum)
    raster2=np.where(raster==NDV,NDV,0)
    raster4=block_reduce(raster2,block_size,func=np.sum)
    raster2=np.where(raster4<0,NDV,raster3)
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
def align_rasters(raster,alignraster,how=np.mean,cxsize=None,cysize=None,masked=False):
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
    NDV, xsize, ysize, GeoT, Projection, DataType=GetGeoInfo(file)
    data=gdalnumeric.LoadFile(file)
    data=np.ma.masked_array(data, mask=data==NDV,fill_value=-np.inf)
    return data

class RasterGeoTError(Exception):
    pass

class RasterGeoError(Exception):
    pass

# GeoRaster Class
class GeoRaster():
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
            raster: Numpy masked array with the raster data, which could be loaded with the load_tiff(file)
            geot: GDAL Geotransformation
            nodata_value: No data value in raster, optional
        '''
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

    def __getitem__(self,indx):
        return self.raster.__getitem__(indx)

    def __getslice__(self,i,j):
        return self.raster.__getslice__(i,j)

    def __getattribute__(self, attr):
        return eval('self.'+attr)

    def __lt__(self,other):
        if isinstance(other, GeoRaster):
            return self.raster<other.raster
        elif isinstance(other, np.ndarray):
            return self.raster<other

    def __le__(self,other):
        if isinstance(other, GeoRaster):
            return self.raster<=other.raster
        elif isinstance(other, np.ndarray):
            return self.raster<=other

    def __gt__(self,other):
        if isinstance(other, GeoRaster):
            return self.raster>other.raster
        elif isinstance(other, np.ndarray):
            return self.raster>other

    def __ge__(self,other):
        if isinstance(other, GeoRaster):
            return self.raster>=other.raster
        elif isinstance(other, np.ndarray):
            return self.raster>=other

    def __eq__(self,other):
        if isinstance(other, GeoRaster):
            return self.raster==other.raster
        elif isinstance(other, np.ndarray):
            return self.raster==other

    def __ne__(self,other):
        if isinstance(other, GeoRaster):
            return self.raster!=other.raster
        elif isinstance(other, np.ndarray):
            return self.raster!=other

    def __pos__(self):
        return self

    def __neg__(self):
        return GeoRaster(-self.raster, self.geot, nodata_value=self.nodata_value, projection = self.projection, datatype = self.datatype)

    def __add__(self, other):
        if isinstance(other,GeoRaster):
            if self.geot==other.geot:
                if self.nodata_value==other.nodata_value:
                    ndv=self.nodata_value
                else:
                    ndv=np.nan
                return GeoRaster(self.raster+other.raster, self.geot, nodata_value=ndv, projection = self.projection, datatype = self.datatype)
            else:
                raise RasterGeoTError("Rasters must have same geotransform. If needed first create union or allign them.")
        else:
            return GeoRaster(self.raster+other, self.geot, self.geot, nodata_value=self.nodata_value, projection = self.projection, datatype = self.datatype)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self+other.__neg__()

    def __rsub__(self, other):
        return self.__sub__(other)

    def __mul__(self, other):
        if isinstance(other,GeoRaster):
            if self.geot==other.geot:
                if self.nodata_value==other.nodata_value:
                    ndv=self.nodata_value
                else:
                    ndv=np.nan
                return GeoRaster(self.raster*other.raster, self.geot, nodata_value=ndv, projection = self.projection, datatype = self.datatype)
            else:
                raise RasterGeoTError("Rasters must have same geotransform. If needed first create union or allign them.")
        else:
            return GeoRaster(self.raster*other, self.geot, self.geot, nodata_value=self.nodata_value, projection = self.projection, datatype = self.datatype)

    def __rmul__(self, other):
        return self.__mul__(other) 

    def __truediv__(self, other):
        if isinstance(other,GeoRaster):
            if self.geot==other.geot:
                if self.nodata_value==other.nodata_value:
                    ndv=self.nodata_value
                else:
                    ndv=np.nan
                return GeoRaster(self.raster/other.raster, self.geot, nodata_value=ndv, projection = self.projection, datatype = self.datatype)
            else:
                raise RasterGeoTError("Rasters must have same geotransform. If needed first create union or allign them.")
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
            if self.geot==other.geot:
                if self.nodata_value==other.nodata_value:
                    ndv=self.nodata_value
                else:
                    ndv=np.nan
                return GeoRaster((other.raster/self.raster).astype(int), self.geot, nodata_value=ndv, projection = self.projection, datatype = self.datatype)
            else:
                raise RasterGeoTError("Rasters must have same geotransform. If needed first create union or allign them.")
        else:
            return GeoRaster((other/self.raster).astype(int), self.geot, nodata_value=self.nodata_value, projection = self.projection, datatype = self.datatype)

    def __pow__(self,other):
        if isinstance(other,GeoRaster):
            if self.geot==other.geot:
                if self.nodata_value==other.nodata_value:
                    ndv=self.nodata_value
                else:
                    ndv=np.nan
                return GeoRaster(self.raster**other.raster, self.geot, nodata_value=ndv, projection = self.projection, datatype = self.datatype)
            else:
                raise RasterGeoTError("Rasters must have same geotransform. If needed first create union or allign them.")
        else:
            return GeoRaster(self.raster**other, self.geot, nodata_value=self.nodata_value, projection = self.projection, datatype = self.datatype)

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
        create_geotiff(filename, self.raster, gdal.GetDriverByName('GTiff'), self.nodata_value, self.shape[1], self.shape[0], self.geot, self.projection, self.datatype)

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
            xsort = sorted(self.raster.data[self.raster.data!=self.nodata_value].flatten()) # increasing order
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
        Align two rasters so that data overlaps by geographical location
        Usage: (alignedraster_o, alignedraster_a) = AlignRasters(raster, alignraster, how=np.mean)
        where 
            raster: string with location of raster to be aligned
            alignraster: string with location of raster to which raster will be aligned
            how: function used to aggregate cells (if the rasters have different sizes)
        It is assumed that both rasters have the same size
        '''
        return align_georasters(self,alignraster,how=how,cxsize=cxsize,cysize=cysize)

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
        for i in rasters:
            (row,col) = map_pixel(i.xmin, i.ymax, rasters[0].x_cell_size, rasters[0].y_cell_size, lonmin, latmax)
            out[row:row+i.shape[0],col:col+i.shape[1]] = np.where(i.raster.data!=i.nodata_value, i.raster.data,\
                                                         out[row:row+i.shape[0],col:col+i.shape[1]])#i.raster
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

