# Copyright (c) 2015, Roger Lew (rogerlew.gmail.com)
# Date: 12/27/2014
# License: BSD (3-clause license)
# 
# The project described was supported by NSF award number IIA-1301792
# from the NSF Idaho EPSCoR Program and by the National Science Foundation.
 
from __future__ import print_function

from collections import OrderedDict

import warnings

import numpy as np
from scipy import ndimage
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource, Normalize

import colormath

from osgeo import gdal, osr, ogr
from osgeo.gdalconst import *


class DEM32:
    """
    make shift class to represent a 32-bit DEM
    """
    
    # I monkey wrench this to make it do what I need as needed
    def __init__(self, fname):

            # open the image
            ds = gdal.Open(fname, GA_ReadOnly)

            if ds is None:
                raise Exception('Could not open "%s"' % fname)

            self.fname = fname
            self.ds = ds
            proj = ds.GetProjectionRef()
            srs = None
            if proj is not None:
                srs = osr.SpatialReference()
                srs.ImportFromWkt(proj)
            self.proj, self.srs = proj, srs

            self.epsg = self._identifyEPSG()
            
            self.band = ds.GetRasterBand(1)
            self.transform = self.ds.GetGeoTransform(can_return_null = True)
            self.xres = self.transform[1]
            self.yres = self.transform[5]

            nx, ny = ds.RasterXSize, ds.RasterYSize
            
            ul = self.ul= self.getLngLat(0, 0)
            ll = self.ll = self.getLngLat(0, ny)
            ur = self.ur= self.getLngLat(nx, 0)
            lr = self.lr = self.getLngLat(nx, ny)
            
            corners = self.corners = [ul, ll, ur, lr]
            self.extents = [min(zip(*corners)[0]), max(zip(*corners)[1]), 
                            max(zip(*corners)[0]), min(zip(*corners)[1])]
          
    def _identifyEPSG(self):
        for L in str(self.srs).split('\n')[::-1]:
            if 'AUTHORITY' in L and 'EPSG' in L:
                return 'epsg:' + \
                       ''.join([c for c in L if c in '0123456789'])
        return None

    def getLngLat(self, x, y):
        assert self.transform != None

        xOrigin, xPixSize, xZero, \
        yOrigin, yZero, yPixSize = self.transform

        assert xZero == yZero == 0.0

        lng = xOrigin + xPixSize*x 
        lat = yOrigin + yPixSize*y  

        return lng, lat

    def getPixelCoords(self, lng, lat):

        if not self.__contains__((lng, lat)):
            warnings.warn('coordinates not in extents')
        
        assert self.transform != None
        
        xOrigin, xPixSize, xZero, \
        yOrigin, yZero, yPixSize = self.transform

        assert xZero == yZero == 0.0
        
        x = (lng - xOrigin) / xPixSize
        y = (lat - yOrigin) / yPixSize

        return x, y
    
    def __contains__(self, (lng, lat)):
        left, upper, right, lower = self.extents
        return left < lng and lng < right and \
               lower < lat and lat < upper

    def getElevation(self, lng, lat, method='cubic'):
        x, y = self.getPixelCoords(lng, lat)
        w, h = self.ds.RasterXSize, self.ds.RasterYSize
        
        if x < 0 or x > w or y < 0 or y > h:
            return float('nan')
        
        if method == 'bilinear':
            _x = [int(floor(x))-1, int(ceil(x))-1]
            if _x[0] < 0 : _x = [0, 1]
            
            _y = [int(floor(y))-1, int(ceil(y))-1]
            if _y[0] < 0 : _y = [0, 1]
            
            data = self.band.ReadAsArray(_x[0], _y[0], 2, 2)
            func = interpolate.interp2d(_x, _y, data, kind='linear')
            z = func(x, y)[0]
                
        elif method == 'cubic':
            xr, yr = int(round(x)), int(round(y))

            if xr - 2 < 0:
                _x = [0, 1, 2, 3, 4]
            elif xr + 2 > w:
                _x = [w-5, w-4, w-3, w-2, w-1]
            else:
                _x = [xr-3, xr-2, xr-1, xr, xr+1]
                
            if yr - 2 < 0:
                _y = [0, 1, 2, 3, 4]
            elif yr + 2 >= h:
                _y = [h-5, h-4, h-3, h-2, h-1]
            else:
                _y = [yr-3, yr-2, yr-1, yr, yr+1]
              
            data = self.band.ReadAsArray(_x[0], _y[0], 5, 5)
            func = interpolate.interp2d(_x, _y, data, kind='cubic')
            
            z = func(x, y)[0]
            
        else:
            x, y = int(round(x)), int(round(y))

            if x == w : x = w-1
            if y == h : y = h-1
            
            z = self.band.ReadAsArray(x, y, 1, 1)

            if z is not None:
                return z[0,0]
            else:
                return float('nan')
            
        return z
    
    def getData(self):
        """
        data isn't loaded by default to conserve memory
        """
        return self.band.ReadAsArray()

    def getElevationLimits(self):
        # read one row at a time to handle huge files
        # this is slower, but saves memory
        
        nx, ny = self.ds.RasterXSize, self.ds.RasterYSize
        band = self.band

        zmin = 1e38
        zmax = -1e38
        
        for i in xrange(ny):
            data = band.ReadAsArray(0, i, nx, 1, ny, 1)

            dmin = np.min(data)
            dmax = np.max(data)

            if dmin < zmin:
                zmin = dmin

            if dmax > zmax:
                zmax = dmax

        return zmin, zmax
    
    def getMeshGrid(self):
        ds = self.ds
        
        nx = ds.RasterXSize
        ny = ds.RasterYSize
        
        xmin, ymax, xmax, ymin = self.extents
        
        xi = np.linspace(xmin, xmax, nx)
        yi = np.linspace(ymin, ymax, ny)
        return np.meshgrid(xi, yi) 
   
    def _calculate_surface_normals(self, scale=1):
        xres, yres = self.xres, self.yres
        
        elevation = self.getData() * scale
        
        # gradient in x and y directions
        dy, dx = np.gradient(elevation)  # find gradient
        dx /=  -xres
        dy /=  -yres
        
        slope = np.arctan(np.hypot(dx, dy))
        dz = np.cos(slope)
        
        d = np.sqrt(dx**2 + dy**2 + dz**2)
        
        return dx/d, dy/d, dz/d
         
    def _write_ds(self, bands, dst_fname, drivername, dtype=GDT_Byte):
        nx, ny = self.ds.RasterXSize, self.ds.RasterYSize
        nb = len(bands)
        
        # initialize raster
        driver = gdal.GetDriverByName(drivername)
        ds = driver.Create(dst_fname, nx, ny, nb, dtype)

        # set projection

        wkt = self.srs.ExportToWkt()
        proj = osr.SpatialReference()
        proj.CopyGeogCSFrom(self.srs)
        status = proj.ImportFromWkt(wkt)

        ds.SetProjection(proj.ExportToWkt())

        # set geotransform
        ds.SetGeoTransform(self.transform)

        # write data
        for i, b in enumerate(bands):
            ds.GetRasterBand(i+1).WriteArray(b)

        ds = None  # Writes and closes file

    def set_nodata_to_zero(self, dst_fname, nodata=0.0, drivername='GTiff'):
        data = self.getData()
        nodata_indx = np.where(data < 0)
        data[nodata_indx] = 0.0
        self._write_ds([data], dst_fname, drivername, dtype=GDT_Float32)
        
    def feet_to_meters(self, dst_fname, drivername='GTiff'):
        data = self.getData()
        data *= 0.3048
        self._write_ds([data], dst_fname, drivername, dtype=GDT_Float32)

    def gaussian_blur(self, dst_fname, sigma=1, drivername='GTiff'):
        data = self.getData()
        data = ndimage.filters.gaussian_filter(data, sigma)
        
        self._write_ds([data], dst_fname, drivername, dtype=GDT_Float32)
                
    def _calculate_slope(self):
        xres, yres = self.xres, self.yres
        elevation = self.getData()
        
        dy, dx = np.gradient(elevation) 
        
        slope = np.arctan(np.hypot(dx/xres, dy/yres)) # slope in radians
        return slope
        
    def _calculate_aspect(self):
        xres, yres = self.xres, self.yres
        elevation = self.getData()
        
        dy, dx = np.gradient(elevation) 
        
        aspect = np.arctan2(dx, dy)

        return aspect

    def gradient_weighted_maps(self):
        slope = self._calculate_slope()
        aspect = self._calculate_aspect()

        slopeW = np.cos(slope)
        slopeWinv = 1.0 - slopeW

        northW = np.cos(aspect)
        northW[np.where(northW < 0)] = 0
        northW *= slopeWinv

        southW = -1 * np.cos(aspect)
        southW[np.where(southW < 0)] = 0
        southW *= slopeWinv

        eastW = -1 * np.sin(aspect)
        eastW[np.where(eastW < 0)] = 0
        eastW *= slopeWinv

        westW = np.sin(aspect)
        westW[np.where(westW < 0)] = 0
        westW *= slopeWinv

        return OrderedDict([('slope', slopeW),
                            ('north', northW),
                            ('south', southW),
                            ('east',  eastW),
                            ('west',  westW)])
        
    def normal_map(self, dst_fname, drivername='GTiff', scale=1.0):
        dx, dy, dz = self._calculate_surface_normals(scale)
        
        r = np.array(127.5 * (dx + 1.0), dtype=np.uint8)
        g = np.array(127.5 * (dy + 1.0), dtype=np.uint8)
        b = np.array(127.5 * (dz + 1.0), dtype=np.uint8)
        
        self._write_ds([r, g, b], dst_fname, drivername)

    def getMaskFromPolyCoords(self, poly_coords):
        """
        build a blurred mask for merging bathymetry into
        DEM

        poly_coords - a list of (lng, lat) coordinates
        src - file name of the DEM used for merging

        returns a 2D np.array with values between [0-1]
        Should be 1 around the convex hull of the zero_coords
        """
        left, top, right, bottom = self.extents
        xs, ys = zip(*poly_coords)

        assert left < np.min(xs)
        assert right > np.max(xs)
        
        assert top > np.max(ys)
        assert bottom < np.min(ys)
        
        ref_data = self.band.ReadAsArray()
        srs = self.srs
        proj = self.proj
        transform = self.transform
        xsize, ysize = self.ds.RasterXSize, self.ds.RasterYSize

        # Create a new raster dataset in memory 
        driver = gdal.GetDriverByName('MEM')
        dst_ds = driver.Create('', xsize, ysize, 
                               1, gdal.GDT_Float32)
        dst_ds.SetGeoTransform(transform)
        dst_ds.SetProjection(proj)

        data = np.zeros((ysize, xsize))
        dst_band = dst_ds.GetRasterBand(1)
        dst_band.WriteArray(data)

        # Create a memory layer to rasterize from.
        rast_ogr_ds = \
                ogr.GetDriverByName('Memory').CreateDataSource( 'wrk' )
        rast_mem_lyr = rast_ogr_ds.CreateLayer( 'poly', srs=srs )

        # Add a polygon.
        coord_str = ','.join(['%f %f'%(lng,lat) for lng,lat in poly_coords])
        wkt_geom = 'POLYGON((' + coord_str + '))'

        feat = ogr.Feature( rast_mem_lyr.GetLayerDefn() )
        feat.SetGeometryDirectly( ogr.Geometry(wkt = wkt_geom) )

        rast_mem_lyr.CreateFeature( feat )

        # Run the algorithm.
        err = gdal.RasterizeLayer( dst_ds, [1], rast_mem_lyr,
                                   burn_values = [1.0] )

        # Pull data back out of the dataset
        data = dst_band.ReadAsArray()
        
        dst_ds = None # close dataset to make sure memory is released

        return data
    
def make_normal(src, dst, drivername='GTiff', scale=1.0):
    gt = DEM32(src)
    gt.normal_map(dst, drivername, scale)

if __name__ == '__main__':

    dem_fn = r"D:\ownCloud\documents\geo_data\DryCreek\Kormos_etal_2013_data\GIS_DATA\DEM.tif"

    dem = DEM32(dem_fn)
    
        
##    import argparse
##
##    parser = argparse.ArgumentParser()
##    parser.add_argument('mode', type=str,
##                        help='Utility has 2 modes: normal, splat')
##    parser.add_argument('src_dataset', type=str,
##                        help='The input raster to be processed')
##    parser.add_argument('dst_dataset', type=str,
##                        help='The output raster produced')
##    parser.add_argument('--drivername',
##                        help='path to elemlist file ("Gtiff")')
##    
##    args = parser.parse_args()
##    mode = args.mode
##    src_dataset = args.src_dataset
##    dst_dataset = args.dst_dataset
##    
##    drivername = (args.drivername, 'GTiff')[args.drivername is None]
##
##    if mode.lower().strip() == 'normal':
##        make_normal(src_dataset, dst_dataset, drivername)
##    
    
