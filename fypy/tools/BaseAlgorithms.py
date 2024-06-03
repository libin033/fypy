# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : BaseAlgorithms.py

@Modify Time :  2024/2/1 15:08   

@Author      : fypy Team    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime
import tempfile
import glob
import time

from osgeo import gdal, gdalconst
from fypy.tools.tifpro import writetiff, GetGDALType
from fypy.tools.ncpro import writenc, writenc_fileinfo
from fypy.tools.hdfpro import writehdf

class BaseAlgorithms :

    def clip(self, srcfile, extent=None, shpname=None, resolution=None,
             dstSRS='EPSG:4326', epsg=None, dstNodata=None,
             resampleAlg=gdalconst.GRA_NearestNeighbour, **kwargs):
        ''' 影像裁剪 '''

        if epsg is not None and dstSRS is None:
            dstSRS = "EPSG:%d" %(int(epsg))

        if shpname is None :
            ds = gdal.Warp("", srcfile, format='MEM',
                           dstSRS=dstSRS, dstNodata=dstNodata,
                           xRes=resolution, yRes=resolution,
                           resampleAlg=resampleAlg,
                           outputBounds=extent, **kwargs)
        else:
            ds = gdal.Warp(
                "", srcfile, format='MEM', dstSRS=dstSRS,
                resampleAlg=resampleAlg, dstNodata=dstNodata,
                cutlineDSName=shpname, cropToCutline=True,
                xRes=resolution, yRes=resolution, **kwargs)

        return ds

    def ds2data(self, srcDS) :
        ''' 读取数据 '''
        if srcDS is None :
            return None
        return srcDS.ReadAsArray()

    def ds2tiff(self, outname, srcDS):
        ''' 保存为GeoTiff文件 '''

        data  = srcDS.ReadAsArray()
        trans = srcDS.GetGeoTransform()
        prj   = srcDS.GetProjection()
        fillvalue = srcDS.GetRasterBand(1).GetNoDataValue()

        writetiff(outname, data, trans, prj, fillvalue=fillvalue)

    def ds2netcdf(self, outname, sdsname, srcDS):
        ''' 保存为NetCDF文件 '''

        if srcDS is None :
            raise Exception('srcDS为None')

        data  = srcDS.ReadAsArray()
        trans = srcDS.GetGeoTransform()
        prj   = srcDS.GetProjection()
        fillvalue = srcDS.GetRasterBand(1).GetNoDataValue()

        if len(data.shape) == 2 :
            level, (height, width) = 1, data.shape
        elif len(data.shape) == 3 :
            level, height, width = data.shape
        else:
            raise Exception('仅暂支持2D或3D数据的输出')

        lon = trans[0] + trans[1] * np.arange(width)
        lat = trans[3] + trans[5] * np.arange(height)


        dictfileinfo = {
            'trans' : trans,
            'prj'   : prj
        }
        writenc_fileinfo(outname, dictfileinfo=dictfileinfo, overwrite=1)
        writenc(outname, 'latitude',  lat, overwrite=0)
        writenc(outname, 'longitude', lon, overwrite=0)
        if level == 1 :
            writenc(outname, sdsname, data, dimension=('latitude', 'longitude'), overwrite=0)
        else:
            writenc(outname, 'level', np.arange(level), overwrite=0)
            writenc(outname, sdsname, data, dimension=('level', 'latitude', 'longitude'),
                    overwrite=0, fill_value=fillvalue)

    def ds2hdf(self, outname, sdsname, srcDS):
        ''' 保存为HDF5文件 '''
        if srcDS is None :
            raise Exception('srcDS为None')

        data  = srcDS.ReadAsArray()
        trans = srcDS.GetGeoTransform()
        prj   = srcDS.GetProjection()
        fillvalue = srcDS.GetRasterBand(1).GetNoDataValue()
        if data is None :
            raise Exception('读取图层数据失败')

        if len(data.shape) == 2 :
            level, (height, width) = 1, data.shape
        elif len(data.shape) == 3 :
            level, height, width = data.shape
        else:
            raise Exception('仅暂支持2D或3D数据的输出')

        lon = trans[0] + trans[1] * np.arange(width)
        lat = trans[3] + trans[5] * np.arange(height)
        dictfileinfo = {
            'trans' : trans,
            'prj'   : prj
        }
        dictsdsinfo = {
            'name' : sdsname,
            'fillvalue' : fillvalue
        }
        writenc_fileinfo(outname, dictfileinfo=dictfileinfo, overwrite=1)
        writehdf(outname, 'latitude',  lat, overwrite=0)
        writehdf(outname, 'longitude',  lon, overwrite=0)
        writehdf(outname, sdsname, data, overwrite=0, dictsdsinfo=dictsdsinfo)

    def __del__(self):
        # 删除系统中由fypy创建的临时文件
        tmp_file = tempfile.NamedTemporaryFile(prefix="tmp_fypy_", delete=True)
        temphdf = tmp_file.name + '.hdf'
        tempdir = os.path.dirname(temphdf)
        filelist = glob.glob(os.path.join(tempdir, 'tmp_fypy_*'))
        for filename in filelist :
            if (time.time() - os.stat(filename).st_mtime ) > 5*24*60*60 :
                if os.path.isfile(filename) :
                    try:
                        os.remove(filename)
                        # print(item)
                    except BaseException as e:
                        print('删除%s失败' %(filename))