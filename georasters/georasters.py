#!/usr/bin/env python
# coding: utf-8
'''
Copyright (C) 2014 Ömer Özak
This program defines functions that are useful for working with GIS data
Usage:

import gisrastertools

or

from gisrastertools import *

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
import pandas as pd
from osgeo import gdal, gdalnumeric, ogr, osr
from gdalconst import *
from skimage.measure import block_reduce
import matplotlib.pyplot as plt

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
			col, row = map_pixel(x,y,GeoT[1],GeoT[-1], GeoT[0],GeoT[3])
	'''
	point_x=np.array(point_x)
	point_y=np.array(point_y)
	col = np.around((point_x - xmin) / cellx).astype(int)
	row = np.around((point_y - ymax) / celly).astype(int)
	return col,row

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
			mcol,row3=map_pixel(GeoT2[0], GeoT2[3], GeoT1[1] *blocksize[0],GeoT1[-1]*blocksize[1], GeoT1[0], GeoT1[3])
			acol=0
		else:
			acol,row3=map_pixel(GeoT1[0], GeoT1[3], GeoT2[1],GeoT2[-1], GeoT2[0], GeoT2[3])
			mcol=0
		if GeoT1[3]<=GeoT2[3]:
			col3,arow=map_pixel(GeoT1[0], GeoT1[3], GeoT2[1],GeoT2[-1], GeoT2[0], GeoT2[3])
			mrow=0
		else:
			col3,mrow=map_pixel(GeoT2[0], GeoT2[3], GeoT1[1] *blocksize[0],GeoT1[-1]*blocksize[1], GeoT1[0], GeoT1[3])
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
		print "Rasters need to be in same projection"
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
    '''
    def __init__(self, raster, geot, nodata_value = np.nan):
        self.raster = raster
        self.geot = geot
        self.nodata_value = nodata_value
        self.shape = raster.shape
        self.xmin = geot[0]
        self.ymax = geot[3]
        self.xmax = self.xmin + geot[1] * self.shape[1]
        self.ymin = self.ymax + geot[-1] * self.shape[0]
        self.bounds = (self.xmin, self.ymin, self.xmax, self.ymax)

    def __pos__(self):
        return self

    def __neg__(self):
        return GeoRaster(-self.raster, self.geot, nodata_value=self.nodata_value)

    def __add__(self, other):
        if isinstance(other,GeoRaster):
            if self.geot==other.geot:
                if self.nodata_value==other.nodata_value:
                    ndv=self.nodata_value
                else:
                    ndv=np.nan
                return GeoRaster(self.raster+other.raster, self.geot, nodata_value=ndv)
            else:
                raise RasterGeoTError("Rasters must have same geotransform. If needed first create union or allign them.")
        else:
            return GeoRaster(self.raster+other, self.geot, nodata_value=self.nodata_value)

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
                return GeoRaster(self.raster*other.raster, self.geot, nodata_value=ndv)
            else:
                raise RasterGeoTError("Rasters must have same geotransform. If needed first create union or allign them.")
        else:
            return GeoRaster(self.raster*other, self.geot, nodata_value=self.nodata_value)

    def __rmul__(self, other):
        return self.__mul__(other) 

    def __truediv__(self, other):
        if isinstance(other,GeoRaster):
            if self.geot==other.geot:
                if self.nodata_value==other.nodata_value:
                    ndv=self.nodata_value
                else:
                    ndv=np.nan
                return GeoRaster(self.raster/other.raster, self.geot, nodata_value=ndv)
            else:
                raise RasterGeoTError("Rasters must have same geotransform. If needed first create union or allign them.")
        else:
            return GeoRaster(self.raster/other, self.geot, nodata_value=self.nodata_value)

    def __rtruediv__(self, other):
        if isinstance(other,GeoRaster):
            return other.__truediv__(self)
        else:
            return GeoRaster(other/self.raster, self.geot, nodata_value=self.nodata_value)

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
                return GeoRaster((other.raster/self.raster).astype(int), self.geot, nodata_value=ndv)
            else:
                raise RasterGeoTError("Rasters must have same geotransform. If needed first create union or allign them.")
        else:
            return GeoRaster((other/self.raster).astype(int), self.geot, nodata_value=self.nodata_value)

    def __pow__(self,other):
        if isinstance(other,GeoRaster):
            if self.geot==other.geot:
                if self.nodata_value==other.nodata_value:
                    ndv=self.nodata_value
                else:
                    ndv=np.nan
                return GeoRaster(self.raster**other.raster, self.geot, nodata_value=ndv)
            else:
                raise RasterGeoTError("Rasters must have same geotransform. If needed first create union or allign them.")
        else:
            return GeoRaster(self.raster**other, self.geot, nodata_value=self.nodata_value)

    def plot(self):
        plt.matshow(self.raster)

    def union(self, other):
        return union([self,other])

# Union of rasters
def union(rasters):
    """
    Union of rasters
    Usage:
        union(rasters)
    where
        rasters is a list of GeoRaster objects
    """
    if sum([rasters[0].geot[1]==i.geot[1] for i in rasters])==len(rasters) and sum([rasters[0].geot[-1]==i.geot[-1] for i in rasters])==len(rasters):
        if sum([rasters[0].nodata_value==i.nodata_value for i in rasters])==len(rasters):
            ndv=rasters[0].nodata_value
        else:
            ndv = np.nan
        lonmin = min([i.xmin for i in rasters])
        lonmax = max([i.xmax for i in rasters])
        latmin = min([i.ymin for i in rasters])
        latmax = max([i.ymax for i in rasters])
        shape = (np.round((latmax-latmin)/rasters[0].geot[1]).astype(int),np.round((lonmax-lonmin)/rasters[0].geot[1]).astype(int))
        out = ndv*np.ones(shape)
        for i in rasters:
            (col,row) = map_pixel(i.xmin, i.ymax, rasters[0].geot[1], rasters[0].geot[-1], lonmin, latmax)
            out[row:row+i.shape[0],col:col+i.shape[1]] = np.where(i.raster.data!=ndv, i.raster.data, out[row:row+i.shape[0],col:col+i.shape[1]])#i.raster
        return GeoRaster(out, (lonmin, rasters[0].geot[1], 0.0, latmax, 0.0, rasters[0].geot[-1]), nodata_value=ndv)
    else:
        raise RasterGeoError('Rasters need to have same pixel sizes. Use the aggregate or dissolve functions to generate correct GeoRasters')
