# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : Overlay.py

@Modify Time :  2024/4/15   

@Author      : Lee    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime
import glob
import geopandas as gpd
from osgeo import gdal, gdalconst, osr, ogr

def erase(outfile, srcfile, erasefile):
    gdf_in = gpd.read_file(srcfile)
    gdf_mask = gpd.read_file(erasefile)

    gdf_diff = gpd.overlay(gdf_in, gdf_mask, how='difference')
    gdf_diff.to_file(outfile, driver='ESRI Shapefile')

def identity():
    pass

def intersect(outname, srcfile, clipfile, layername=None):
    '''
    矢量求交，裁剪

    Parameters
    ----------
    outname : str
        输出裁剪结果文件
    srcfile ： str
        被裁剪对象
    clipfile ：str
        掩膜对象
    layername ： str, optional
        输出文件名图层名称

    Returns
    -------

    '''
    if layername is None :
        layername = os.path.basename(outname).replace('.shp','')

    tempfiles = None
    outdir = os.path.dirname(outname)
    if not os.path.isdir(outdir) :
        os.makedirs(outdir)

    driver = ogr.GetDriverByName('ESRI Shapefile')
    srcds1 = driver.Open(srcfile, gdalconst.GA_ReadOnly)
    srcds2 = driver.Open(clipfile, gdalconst.GA_ReadOnly)
    src_layer1 = srcds1.GetLayer()
    src_layer2 = srcds2.GetLayer()
    srs1 = src_layer1.GetSpatialRef()
    srs2 = src_layer2.GetSpatialRef()
    if srs1.GetAttrValue('AUTHORITY',1) != srs2.GetAttrValue('AUTHORITY',1):
        print("空间参考不一致!将进行投影转换投影转换")

        # 对clip文件进行投影转换
        dstpatialRef = get_spatialref(srcfile)
        tempfile = os.path.join(outdir, os.path.basename(clipfile))
        tempfiles = glob.glob(tempfile.replace('.shp', '.*'))
        for item in tempfiles :
            try:
                os.remove(item)
            except BaseException as e :
                continue
        transform(tempfile, clipfile, dstpatialRef=dstpatialRef)

        srcds2.Destroy()
        # 获取转换后的结果
        srcds2 = driver.Open(tempfile, gdalconst.GA_ReadOnly)
        src_layer2 = srcds2.GetLayer()
        srs2 = src_layer2.GetSpatialRef()
        if srs1.GetAttrValue('AUTHORITY',1) != srs2.GetAttrValue('AUTHORITY',1):
            raise Exception("空间参考不一致!将进行投影转换投影转换")

        tempfiles = glob.glob(tempfile.replace('.shp', '.*'))

    # 创建输出文件
    target_ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource(outname)
    target_layer = target_ds.CreateLayer(layername, srs1,
                                         geom_type=src_layer1.GetGeomType(),
                                         options=["ENCODING=UTF-8"]) # 设置编码为UTF-8，防止中文出现乱码

    # 定义输出属性表信息
    feature = src_layer1.GetFeature(0)
    field_count = feature.GetFieldCount()
    field_names = []
    for attr in range(field_count):
        field_defn = feature.GetFieldDefnRef(attr)
        field_names.append(field_defn.GetName())
        target_layer.CreateField(field_defn)

    for feat1 in src_layer1:
        geom1 = feat1.GetGeometryRef()
        for feat2 in src_layer2:
            geom2 = feat2.GetGeometryRef()
            # 判断两个feature是否相交
            if not geom1.Intersects(geom2):
                continue
            # intersect = geom1.Intersection(geom1)
            feature = ogr.Feature(target_layer.GetLayerDefn())
            feature.SetGeometry(geom1)
            # 将源文件中的字段信息写入匹配feature中
            for field_name in field_names:
                feature.SetField(field_name, feat1.GetField(field_name))
            target_layer.CreateFeature(feature)

            feat1.Destroy()
            feat2.Destroy()
            break

    # 清理引用
    target_ds.Destroy()
    srcds1.Destroy()
    srcds2.Destroy()
    if tempfiles is not None :
        for item in tempfiles :
            try:
                os.remove(item)
            except BaseException as e :
                continue

def union(outfile, srcfile, maskfile):
    gdf_in = gpd.read_file(srcfile)
    gdf_mask = gpd.read_file(maskfile)

    gdf_diff = gpd.overlay(gdf_in, gdf_mask, how='union')
    gdf_diff.to_file(outfile, driver='ESRI Shapefile')

def update(self):
    pass


def transform(outname, srcshp, dstpatialRef):
    # 当前地理参考
    srcpatialRef = get_spatialref(srcshp)
    # dstpatialRef = self.get_spatialref(dstshp)

    ds = ogr.Open(srcshp)
    srclayer = ds.GetLayer(0)

    basename = os.path.basename(outname)
    driver = ogr.GetDriverByName('ESRI Shapefile')
    out_ds = driver.CreateDataSource(outname)
    outlayer = out_ds.CreateLayer(basename[:-4], geom_type=srclayer.GetGeomType())

    # 投影转换
    coordinate_transfor = osr.CoordinateTransformation(srcpatialRef, dstpatialRef)

    # 定义输出属性表信息
    feature = srclayer.GetFeature(0)
    field_count = feature.GetFieldCount()
    field_names = []
    for attr in range(field_count):
        field_defn = feature.GetFieldDefnRef(attr)
        field_names.append(field_defn.GetName())
        outlayer.CreateField(field_defn)

    out_fielddefn = outlayer.GetLayerDefn()

    for feature in srclayer:
        geometry = feature.GetGeometryRef()
        geometry.Transform(coordinate_transfor)

        out_feature = ogr.Feature(out_fielddefn)
        out_feature.SetGeometry(geometry)
        for field_name in field_names:
            out_feature.SetField(field_name,feature.GetField(field_name))
        outlayer.CreateFeature(out_feature)

        feature.Destroy()
        out_feature.Destroy()

    # 清除缓存
    ds.Destroy()
    out_ds.Destroy()
    # 保存投影文件
    dstpatialRef.MorphFromESRI()
    prjname = outname.replace(".shp",".prj")
    fn = open(prjname,'w')
    fn.write(dstpatialRef.ExportToWkt())
    fn.close()

def get_spatialref(shpname):
    ''' 获取矢量图层的空间参考信息 '''
    ds = ogr.Open(shpname)
    layer = ds.GetLayer(0)

    # 当前地理参考
    spatialRef = layer.GetFeature(0).GetGeometryRef().GetSpatialReference()
    # 清除缓存
    ds.Destroy()

    return spatialRef

def get_geomType(shpname):
    ds = ogr.Open(shpname)
    layer = ds.GetLayer(0)

    # 图层类型：ogr.wkbPoint、ogr.wkbLineString、ogr.wkbPolygon
    geomType = layer.GetGeomType()

    # 清除缓存
    ds.Destroy()

    return geomType


