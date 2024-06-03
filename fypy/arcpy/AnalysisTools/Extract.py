# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : Extract.py

@Modify Time :  2024/4/15   

@Author      : Lee    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime
from osgeo import gdal, gdalconst, osr, ogr
import geopandas as gpd

def clip(outfile, srcfile, maskfile):
    gdf_in = gpd.read_file(srcfile)
    gdf_mask = gpd.read_file(maskfile)

    gdf_diff = gpd.overlay(gdf_in, gdf_mask, how='intersection')
    gdf_diff.to_file(outfile, driver='ESRI Shapefile')

def select():
    pass

def split(outpath, srcshp, field):
    """
    拆分矢量文件为多个单要素矢量文件,注意拆分后的fid需要重置。
    :param srcshp: 要拆分的矢量文件
    :param outpath: 生成的矢量文件保存目录
    :param field: 根据提供的字段进行矢量拆分
    :return:
    """
    if not os.path.isdir(outpath) :
        os.makedirs(outpath)

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
        driverName = "ESRI Shapefile"
        driver = ogr.GetDriverByName(driverName)
        outShapeFileName = fileName + ".shp"
        outShapeFilePath = os.path.join(outpath, outShapeFileName)
        outData = driver.CreateDataSource(outShapeFilePath)
        gdal.SetConfigOption("SHAPE_ENCODING", "CP936")
        outLayer = outData.CreateLayer(layerName, spatial, geomType)
        for fieldIndex in range(fieldCount):
            fieldDefn = layerDefn.GetFieldDefn(fieldIndex)
            outLayer.CreateField(fieldDefn)
        outLayer = outData.GetLayer()
        feature.SetFID(0)
        outLayer.CreateFeature(feature)
        outData.Destroy()
        feature = layer.GetNextFeature()
    data.Destroy()

    return True

def splitByAttr():
    pass
