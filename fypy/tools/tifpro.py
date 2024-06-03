# -*- coding:utf-8 -*-
'''
@Project     : code

@File        : tifpro.py

@Modify Time :  2022/11/28 10:19

@Author      : fypy Team

@Version     : 1.0

@Description :

'''



import numpy as np

from osgeo import gdal, gdalconst, osr, ogr


def readtiff(filename):
    '''
    读取TIFF文件
    :param filename: 输入TIFF文件名
    :return: 数据、仿射函数、投影参数
    '''

    dataset = gdal.Open(filename, gdal.GA_ReadOnly)
    if dataset == None:
        print(filename+"文件无法打开")
        return None, None, None
    width = dataset.RasterXSize
    height = dataset.RasterYSize
    bands = dataset.RasterCount

    data = dataset.ReadAsArray()
    trans = dataset.GetGeoTransform()
    prj = dataset.GetProjection()

    # 对填充值进行NAN填充
    for iband in range(bands) :
        fillvalue = dataset.GetRasterBand(iband+1).GetNoDataValue()
        if not fillvalue is None :
            if data.dtype == np.float16 or data.dtype == np.float32 or data.dtype == np.float64 :
                data[data==fillvalue] = np.nan

            # if im_data.dtype == np.int16 or im_data.dtype == np.int32 or im_data.dtype == np.int64 :
            #     im_data = np.array(im_data, dtype=np.float32)
            #     im_data[im_data==fillvalue] = np.nan

        # 获取各个波段的斜率和截距
        scale = dataset.GetRasterBand(iband+1).GetScale()
        offset = dataset.GetRasterBand(iband+1).GetOffset()
        if not scale is None :
            data[iband] *= scale

        if not offset is None :
            data[iband] += offset

    del dataset

    return data, trans, prj

def writetiff(outname, data, trans=None, prj=None, fillvalue=None, scale=None, offset=None):
    '''
    输出geotiff文件
    :param outname: 输出文件名
    :param data: 输出数据
    :param trans: 仿射函数
    :param prj: 投影信息
    :param fillvalue: 填充值，作为Nodata
    :return: None
    '''
    im_data = data.copy()
    im_data = np.array(im_data)
    if fillvalue :
        im_data[np.isnan(im_data)] = fillvalue

    datatype = GetGDALType(im_data.dtype)

    if len(im_data.shape) == 3:
        im_bands, im_height, im_width = im_data.shape
    elif len(im_data.shape) == 2:
        im_bands, (im_height, im_width) = 1,im_data.shape
    else:
        im_bands, (im_height, im_width) = 1,im_data.shape
        #创建文件
    driver = gdal.GetDriverByName("GTiff")
    dataset = driver.Create(outname, im_width, im_height, im_bands, datatype,  options=["COMPRESS=LZW",
                                                                                        'BIGTIFF=IF_SAFER'])
    if dataset is None :
        return None

    #写入仿射变换参数
    if not trans is None :
        dataset.SetGeoTransform(trans)

    #写入投影
    if prj is None :
        # 获取地理坐标系统信息，用于选取需要的地理坐标系统
        sr = osr.SpatialReference()
        sr.SetWellKnownGeogCS("EPSG:4326")
        # 给新建图层赋予投影信息
        dataset.SetProjection(sr.ExportToWkt())
    else:
        dataset.SetProjection(prj)

    if im_bands == 1:
        dataset.GetRasterBand(1).WriteArray(im_data)
        if fillvalue is not None:
            dataset.GetRasterBand(1).SetNoDataValue(float(fillvalue))
    else:
        for i in range(im_bands):
            dataset.GetRasterBand(i+1).WriteArray(im_data[i])
            if fillvalue is not None:
                dataset.GetRasterBand(i+1).SetNoDataValue(float(fillvalue))

    del dataset


def getDateSet(data, trans=None, prj=None, fillvalue=None, epsg=4326, scale=None, offset=None):
    '''

    :param data: 输出数据
    :param trans: 仿射函数
    :param prj: 投影信息
    :param fillvalue: 填充值，作为Nodata
    :return: 返回dataset
    '''
    im_data = data.copy()
    im_data = np.array(im_data)
    if fillvalue :
        im_data[np.isnan(im_data)] = fillvalue

    datatype = GetGDALType(im_data.dtype)

    if len(im_data.shape) == 3:
        im_bands, im_height, im_width = im_data.shape
    elif len(im_data.shape) == 2:
        im_bands, (im_height, im_width) = 1,im_data.shape
    else:
        im_bands, (im_height, im_width) = 1,im_data.shape
    #创建文件
    driver = gdal.GetDriverByName("MEM")
    dataset = driver.Create('', im_width, im_height, im_bands, datatype)
    if dataset is  None:
        return None

    if not trans is None :
        dataset.SetGeoTransform(trans) #写入仿射变换参数

    if prj is None :
        # 获取地理坐标系统信息，用于选取需要的地理坐标系统
        sr = osr.SpatialReference()
        dstSRS = 'EPSG:%d' %(int(epsg))
        sr.SetWellKnownGeogCS(dstSRS)
        # 给新建图层赋予投影信息
        dataset.SetProjection(sr.ExportToWkt())
    else:
        dataset.SetProjection(prj) #写入投影

    if im_bands == 1:
        dataset.GetRasterBand(1).WriteArray(im_data)
        if fillvalue is not None:
            dataset.GetRasterBand(1).SetNoDataValue(float(fillvalue))
    else:
        for i in range(im_bands):
            dataset.GetRasterBand(i+1).WriteArray(im_data[i])
            if fillvalue is not None:
                dataset.GetRasterBand(i+1).SetNoDataValue(float(fillvalue))

    return dataset


def get_fillvalue(filename):
    dataset = gdal.Open(filename, gdal.GA_ReadOnly)
    if dataset == None:
        print(filename+"文件无法打开")
        return None, None, None

    fillvalue = dataset.GetRasterBand(1).GetNoDataValue()

    return fillvalue


def GetGDALType(dtype):
    '''
    根据numpy的数据类型，匹配GDAL中的数据类型
    :param dtype:
    :return: GDAL数据类型
    '''

    if dtype == np.byte or dtype == np.uint8:
        return gdal.GDT_Byte
    elif dtype == np.uint16 :
        return gdal.GDT_UInt16
    elif dtype == np.int16 :
        return gdal.GDT_Int16
    elif dtype == np.uint32 :
        return gdal.GDT_UInt32
    elif dtype == np.int32 :
        return gdal.GDT_Int32
    elif dtype == np.float32 or dtype.str in ['>f4', '<f4']:
        return gdal.GDT_Float32
    elif dtype == np.float64 or dtype.str in ['>f8', '<f8']:
        return gdal.GDT_Float64
    else:
        return gdal.GDT_Unknown

def GetNumpyType(dtype):
    if dtype == gdal.GDT_Byte :
        return np.byte
    elif dtype == gdal.GDT_UInt16 :
        return np.uint16
    elif dtype == gdal.GDT_Int16 :
        return np.int16
    elif dtype == gdal.GDT_UInt32 :
        return np.uint32
    elif dtype == gdal.GDT_Int32 :
        return np.int32
    elif dtype == gdal.GDT_Float32 :
        return np.float32
    elif dtype == gdal.GDT_Float64 :
        return np.float64
    else:
        return None
