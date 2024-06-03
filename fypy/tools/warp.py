# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : warp.py

@Modify Time :  2023/10/10 17:54   

@Author      : fypy Team    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime
import tempfile
from osgeo import gdal, ogr, gdalconst
from fypy import parm
ParmDir = os.path.abspath(list(parm.__path__)[0])



def vectorSplit(tempshp, srcshp, regionID, field='NAME'):
    """
    拆分矢量文件为多个单要素矢量文件,注意拆分后的fid需要重置。
    :param srcshp: 要拆分的矢量文件
    :param tempshp: 生成的矢量文件保存目录
    :param field: 根据提供的字段进行矢量拆分
    :return:
    """
    # if not os.path.isdir(outpath) :
    #     os.makedirs(outpath)
    #

    matchID = getField(srcshp, field, regionID)
    if len(matchID) == 0 :
        return None
    elif len(matchID) > 1 :
        print(matchID)
        raise Exception('匹配到多个行政区划，请指定具体行政区划')


    data = ogr.Open(srcshp)
    layer = data.GetLayer()
    spatial = layer.GetSpatialRef()
    geomType = layer.GetGeomType()
    layerDefn = layer.GetLayerDefn()
    fieldCount = layerDefn.GetFieldCount()
    feature = layer.GetNextFeature()
    while feature:
        fid = feature.GetFieldAsString(field)
        fileName, layerName = str(fid), str(fid)
        if regionID in fileName :
            # outShapeFileName = fileName + ".shp"
            # outShapeFilePath = os.path.join(outpath, outShapeFileName)
            # if os.path.isfile(outShapeFilePath) :
            #     return outShapeFilePath

            driverName = "ESRI Shapefile"
            driver = ogr.GetDriverByName(driverName)
            outData = driver.CreateDataSource(tempshp)
            gdal.SetConfigOption("SHAPE_ENCODING", "CP936")
            outLayer = outData.CreateLayer(layerName, spatial, geomType)
            for fieldIndex in range(fieldCount):
                fieldDefn = layerDefn.GetFieldDefn(fieldIndex)
                outLayer.CreateField(fieldDefn)
            outLayer = outData.GetLayer()
            feature.SetFID(0)
            outLayer.CreateFeature(feature)
            outData.Destroy()

            return tempshp

        feature = layer.GetNextFeature()

    data.Destroy()

    return None

def getField(srcfile, fieldname, regionID):

    fieldvalue = []
    driver = ogr.GetDriverByName('ESRI Shapefile')
    srcds = driver.Open(srcfile, gdalconst.GA_ReadOnly)
    src_layer = srcds.GetLayer()
    for feat in src_layer:
        value = feat.GetField(fieldname)
        if regionID in value :
            fieldvalue.append(value)

    return fieldvalue

def clipByAdmin(regionID, adminlevel='sheng'):
    tmp_file = tempfile.NamedTemporaryFile(prefix="tmp_fypy_adminshp_", delete=True)
    tempshp = tmp_file.name + '.shp'

    if adminlevel in ['province', 'sheng', '省', 1] :
        admin = 'sheng'
    elif adminlevel in ['city', 'shi', '市', 2] :
        admin = 'shi'
    elif adminlevel in ['county', 'xian', '县', 3] :
        admin = 'xian'
    else:
        raise Exception('请输入正确的行政等级【sheng, shi, xian】')

    srcshp = os.path.join(ParmDir, 'shapefile', 'china_%s_polygon.shp' %(admin))

    maskshp = vectorSplit(tempshp, srcshp=srcshp, regionID=regionID)
    if maskshp is None:
        raise Exception('未找到指定区县【%s】，请核查' %(regionID))

    return maskshp

