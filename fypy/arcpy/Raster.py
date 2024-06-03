# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : Raster.py

@Modify Time :  2024/3/15 9:36   

@Author      : fypy Team    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime

from osgeo import gdal, gdalconst, ogr, osr

gdal.SetConfigOption("GDAL_FILENAME_IS_UTF8", "YES")
gdal.SetConfigOption("SHAPE_ENCODING", 'GBK')
gdal.UseExceptions()
ogr.RegisterAll()

def set_encoding(encoding='GBK') :
    gdal.SetConfigOption("SHAPE_ENCODING", encoding)

def clip(filename, extent=None, shpname=None, resolution=None,
         dstSRS='EPSG:4326', dstNodata=None, resampleAlg=gdalconst.GRA_NearestNeighbour):

    if shpname is None :
        ds = gdal.Warp("", filename, format='MEM',
                       dstSRS=dstSRS, dstNodata=dstNodata,
                       xRes=resolution, yRes=resolution,
                       resampleAlg=resampleAlg,
                       outputBounds=extent)
    else:
        ds = gdal.Warp(
            "", filename, format='MEM', dstSRS=dstSRS,
            resampleAlg=resampleAlg, dstNodata=dstNodata,
            cutlineDSName=shpname, cropToCutline=True,
            xRes=resolution, yRes=resolution)

    if ds is None :
        return None, None, None

    data = ds.ReadAsArray()
    trans = ds.GetGeoTransform()      # 获取仿射矩阵信息
    prj = ds.GetProjection()          # 获取投影信息

    return data, trans, prj


def mask(InRaster=None, MaskShp=None, burn_name=None, burn_values=None, ref_data=None, trans=None, prj=None):
    '''
    根据参考的tiff文件创建同等投影类型和数据维度的掩膜文件
    :param InRaster: 参考模板的geotiff文件
    :param MaskShp: 待转换的掩膜矢量文件
    :param burn_name: 栅格化的属性字段名称
    :param burn_values:
    :return:
    '''

    if ref_data is None or trans is None or prj is None :
        ref_info = gdal.Open(InRaster)
        ref_data = ref_info.GetRasterBand(1).ReadAsArray()

        trans = ref_info.GetGeoTransform()
        prj = ref_info.GetProjection()
        ref_info = None

    y_size, x_size = ref_data.shape

    vector=ogr.Open(MaskShp)
    vl=vector.GetLayer(0)
    dst_ds = gdal.GetDriverByName('MEM').Create('', x_size, y_size,1,gdal.GDT_Int32)
    dst_ds.SetGeoTransform(trans)
    dst_ds.SetProjection(prj)
    if burn_name is not None :
        gdal.RasterizeLayer(dst_ds,[1],vl,options=["ATTRIBUTE=%s" %burn_name])
    if burn_values is not None :
        gdal.RasterizeLayer(dst_ds,[1],vl,burn_values=burn_values)
    vector=None
    bd = dst_ds.GetRasterBand(1)
    data = bd.ReadAsArray()
    # dst_ds.FlushCache()
    dst_ds=None

    return data, trans, prj

def merge(self, filelist, outname=None, resolution=None, fillvalue=None,
          epsg=4326, extent=None, resampleAlg=0):

    if resampleAlg == 0:
        resampleType = gdalconst.GRA_NearestNeighbour
    elif resampleAlg == 1:
        resampleType = gdalconst.GRA_Bilinear
    else:
        resampleType = gdalconst.GRA_Cubic

    if outname is None :
        mode = 'MEM'
        outname=''
    else:
        mode = 'GTiff'
    options = gdal.WarpOptions(
        dstSRS='EPSG:%d' %(epsg),
        format=mode,
        outputBounds=extent,    # (minX, minY, maxX, maxY)
        resampleAlg=resampleType,
        xRes=resolution, yRes=resolution,
        dstNodata=fillvalue,
        creationOptions=["COMPRESS=LZW"])
    ds = gdal.Warp(outname, filelist, options=options)

    if outname is None :
        im_data = ds.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)#获取数据
        ds = None

        return im_data
    else:
        return outname

def to_vector(shpname, rastername, layername=None, field=None, geom_type=ogr.wkbMultiPolygon):
    '''
    栅格转矢量
    :param shpname: 输出的矢量文件名
    :param rastername: 输入的栅格数据名
    :param layername: option, 待转图层名称
    :param field: option，图层属性字段名
    :return:
    '''
    rasterid = gdal.Open(rastername)
    inband = rasterid.GetRasterBand(1)
    maskband = inband.GetMaskBand()

    prj = osr.SpatialReference()
    prj.ImportFromWkt(rasterid.GetProjection())

    drv = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.isfile(shpname) :
        drv.DeleteDataSource(shpname)

    polygon = drv.CreateDataSource(shpname)
    if layername is not None :
        poly_layer = polygon.CreateLayer(layername, srs=prj,
                                         geom_type=geom_type)
    else:
        poly_layer = polygon.CreateLayer(rastername[:-4], srs=prj,
                                         geom_type=geom_type)

    if field is not None :
        newfield = ogr.FieldDefn(field, ogr.OFTReal)
    else:
        newfield = ogr.FieldDefn('value', ogr.OFTReal)
    poly_layer.CreateField(newfield)

    gdal.FPolygonize(inband, maskband, poly_layer, 0)
    polygon.SyncToDisk()
    polygon = None

def __xuanran(self, ds, colorlist):
    ''' GeoTiff 栅格渲染 https://blog.csdn.net/amyniez/article/details/114765006 '''
    ct = gdal.ColorTable()
    for i in range(len(colorlist)):
        value1, color1, label1 = colorlist[i]
        # value2, color2, label2 = colorlist[i+1]
        # ct.CreateColorRamp(value1, color1,
        #                    value2, color2)
        ct.SetColorEntry(value1, color1)
    band = ds.GetRasterBand(1)
    band.SetRasterColorTable(ct)
    band.SetRasterColorInterpretation(gdal.GCI_PaletteIndex)
    # ct.CreateColorRamp(value1, color1,
    #                    value2, color2)
    # band = ds.GetRasterBand(i + 1)
    # band.SetRasterColorTable(ct)

    return ds

def __jingjiaozhun(self, ds, PosInfo):
    ''' https://blog.csdn.net/amyniez/article/details/114765006 '''
    srs = osr.SpatialReference()
    srs.SetWellKnownGeogCS('WGS84')

    gcps = []
    for x, y, z, line, pixel in PosInfo :
        gcps.append(gdal.GCP(x, y, z, line, pixel))

    ds.SetGCPs(gcps, srs.ExportToWkt())
    ds.SetProjection(srs.ExportToWkt())

    ds.SetGeoTransform(gdal.GCPsToGeoTransform(gcps))









