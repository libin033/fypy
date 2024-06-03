# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : Feature.py

@Modify Time :  2024/3/15 9:53   

@Author      : fypy Team    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime
import glob
import shapefile
from tqdm import tqdm
from fypy.tools.jsonpro import writejson
from osgeo import gdal, gdalconst, ogr, osr
import geopandas as gpd


gdal.SetConfigOption("GDAL_FILENAME_IS_UTF8", "YES")
gdal.SetConfigOption("SHAPE_ENCODING", 'GBK')
gdal.UseExceptions()
ogr.RegisterAll()

def set_encoding(encoding='GBK') :
    gdal.SetConfigOption("SHAPE_ENCODING", encoding)

def erase(outfile, srcfile, erasefile, ) :

    gdf_in = gpd.read_file(srcfile)
    gdf_mask = gpd.read_file(erasefile)

    gdf_diff = gpd.overlay(gdf_in, gdf_mask, how='difference')
    gdf_diff.to_file(outfile, driver='ESRI Shapefile')

def clip(outfile, srcfile, maskfile):
    gdf_in = gpd.read_file(srcfile)
    gdf_mask = gpd.read_file(maskfile)

    gdf_diff = gpd.overlay(gdf_in, gdf_mask, how='intersection')
    gdf_diff.to_file(outfile, driver='ESRI Shapefile')

def union(outfile, srcfile, maskfile):
    gdf_in = gpd.read_file(srcfile)
    gdf_mask = gpd.read_file(maskfile)

    gdf_diff = gpd.overlay(gdf_in, gdf_mask, how='union')
    gdf_diff.to_file(outfile, driver='ESRI Shapefile')

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

def to_raster(outname, shpname, field, resolution,
              fillvalue=None, all_touch=False):
    '''
    矢量转栅格
    :param outname: 输出栅格文件的名称
    :param shpname: 输入矢量文件的名称
    :param resolution: 输出栅格数据的分辨率
    :param field: str, 矢量文件中的字段信息转换成栅格数据的像元值
    :param fillvalue: 输出tif的填充值（NoData）
    :param all_touch: bool，是否进行八方向处理
    :return:
    '''
    vec = ogr.Open(shpname)
    layer = vec.GetLayer()

    extent = layer.GetExtent()
    srs = layer.GetSpatialRef()

    feature = layer.GetFeature(0)
    # for feat in layer.schema :
    #     name = feat.GetName()#获取字段名称
    #     if name == field :
    #         type = feat.GetTypeName()#获取字段类型
    #         type = feat.GetField(field)#获取字段类型
    #         print('字段名称：%s ,字段类型：%s'%(name, type))

    fid = feature.GetField(field)
    dtype = str(type(fid))

    dtype = getNPDType(dtype)

    lonmin, lonmax, latmin, latmax = extent

    line = int(np.ceil((latmax - latmin)/resolution))
    pixel = int(np.ceil((lonmax - lonmin)/resolution))
    band = 1

    trans = [lonmin, resolution, 0, latmax, 0, -resolution]

    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(outname, pixel, line, band, dtype)

    ds.SetGeoTransform(trans)
    ds.SetProjection(str(srs))

    bands = ds.GetRasterBand(1)
    if fillvalue is not None :
        bands.SetNoDataValue(fillvalue)

    # bands.FlushCache()

    if all_touch :
        all_touch = 'True'
    else:
        all_touch = 'False'

    if field :
        gdal.RasterizeLayer(ds, [1], layer,
                            options=["ALL_TOUCHED="+all_touch, "ATTRIBUTE="+field])
    else:
        gdal.RasterizeLayer(ds, [1], layer,
                            options=["ALL_TOUCHED="+all_touch])

def to_geojson(outname, shapename, shp_encoding='utf-8', json_encoding='utf-8'):

    reader = shapefile.Reader(shapename, encoding=shp_encoding)
    fields = reader.fields[1:]
    fieldname = [field[0] for field in fields]
    buffer = []
    for sr in tqdm(reader.shapeRecords()) :
        record = sr.record

        record = [r.decode('gb2312', 'ignore') if isinstance(r, bytes) else r for r in record]

        attr = dict(zip(fieldname, record))

        geom = sr.shape.__geo_interface__

        buffer.append(dict(type='Feature',
                           geometry=geom,
                           properties=attr))

    writejson(outname, dict(type= 'FeatureCollection', features=buffer))

def vectorTranslate(shpname, outdir,
                    format="GeoJSON", accessMode=None, dstSrsESPG=4326,
                    selectFields=None, geometryType="POLYGON", dim="XY"):
    '''
    转换矢量文件，包括坐标系，名称，格式，字段，类型，纬度等。

    Parameters
    ----------
    shpname : str
        要转换的矢量文件
    outdir : str
        生成矢量文件保存目录
    format : str
        矢量文件格式，强烈建议不要使用ESRI Shapefile格式。
    accessMode : str
        None代表creation,'update','append','overwrite'
    dstSrsESPG : str
        目标坐标系EPSG代码，4326是wgs84地理坐标系
    selectFields : str
        需要保留的字段列表如果都保留，则为None
    geometryType : str
        几何类型,"POLYGON","POINT"。。。
    dim : str
        新矢量文件坐标纬度,建议查阅官方API。
    Returns
    -------

    '''

    data = ogr.Open(shpname)
    layer = data.GetLayer()
    spatial = layer.GetSpatialRef()
    layerName = layer.GetName()
    data.Destroy()
    dstSRS = osr.SpatialReference()
    dstSRS.ImportFromEPSG(int(dstSrsESPG))
    if format == "GeoJSON":
        destDataName = layerName + ".geojson"
        destDataPath = os.path.join(outdir, destDataName)
    elif format == "ESRI Shapefile":
        destDataName = os.path.join(outdir, layerName)
        flag = os.path.exists(destDataName)
        os.makedirs(destDataName) if not flag else None
        destDataPath = os.path.join(destDataName, layerName + ".shp")
    else:
        print("不支持该格式！")
        return
    options = gdal.VectorTranslateOptions(
        format=format,
        accessMode=accessMode,
        srcSRS=spatial,
        dstSRS=dstSRS,
        reproject=True,
        selectFields=selectFields,
        layerName=layerName,
        geometryType=geometryType,
        dim=dim
    )
    gdal.VectorTranslate(
        destDataPath,
        srcDS=shpname,
        options=options
    )
    return destDataPath

def createPoint(shpname, coords, dictfield={}, layername=None, epsg=4326):
    '''
    利用GDAL 创建点矢量，将point转换成矢量数据

    Parameters
    ----------
    shpname : str
        输出文件名
    coords: list
        eg. [[lon1, lat1], [lon2, lat2],... ]
    epsg: int

    Returns
    -------
        None
    '''

    driver = ogr.GetDriverByName('ESRI Shapefile')
    # 如果文件存在，则删除后再创建
    if os.access(shpname, os.F_OK) :
        driver.DeleteDataSource(shpname)

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(int(epsg))

    if layername is None :
        layername = os.path.basename(shpname).replace('.shp','')

    ds = driver.CreateDataSource(shpname)
    lyr = ds.CreateLayer(layername, srs, ogr.wkbPoint)
    for key in dictfield :
        # fieldName = ogr.FieldDefn(key, ogr.OFTString)
        # fieldName.SetWidth(50)
        dtype = getOGRType(np.array(dictfield[key]).dtype)
        fieldName = ogr.FieldDefn(key, dtype)
        if dtype == ogr.OFTString :
            fieldName.SetWidth(200)

        lyr.CreateField(fieldName)

    index = 0
    if isinstance(coords, list) :
        for item in coords:
            val_lat = float(item[1])
            val_lon = float(item[0])

            wkt = 'POINT(%f %f)' %(val_lon, val_lat)

            geom = ogr.CreateGeometryFromWkt(wkt)
            feature = ogr.Feature(lyr.GetLayerDefn())
            feature.SetGeometry(geom)
            for key in dictfield :
                feature.SetField(key, dictfield[key][index])
            lyr.CreateFeature(feature)
            index += 1
    else:
        raise Exception('coords类型必须是list')

    ds = None
    driver =None

def createPolyline(shpname, coords, dictfield={}, layername=None, epsg=4326):
    '''
    利用GDAL 创建线矢量，将polyline转换成矢量数据

    Parameters
    ----------
    shpname : str
        输出文件名
    coords: list
        eg. [
                [[lon11, lat11], [lon12, lat12],... ],  # line1
                [[lon21, lat21], [lon22, lat22],... ],  # line2
            ]
    epsg: int

    Returns
    -------
        None
    '''

    driver = ogr.GetDriverByName('ESRI Shapefile')
    if os.access(shpname, os.F_OK) :
        driver.DeleteDataSource(shpname)

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(int(epsg))

    if layername is None :
        layername = os.path.basename(shpname).replace('.shp','')

    ds = driver.CreateDataSource(shpname)
    lyr = ds.CreateLayer(layername, srs, ogr.wkbLineString)

    for key in dictfield :
        dtype = getOGRType(np.array(dictfield[key]).dtype)
        fieldName = ogr.FieldDefn(key, dtype)
        if dtype == ogr.OFTString :
            fieldName.SetWidth(200)
        lyr.CreateField(fieldName)

    index = 0
    if isinstance(coords, list) :
        # 逐线创建
        for coord in coords :
            count = len(coord)
            wkt = 'LINESTRING('
            for i in range(count) :
                if i == count-1:
                    wkt += '%f %f' %(coord[i][0], coord[i][1])
                else:
                    wkt += '%f %f,' %(coord[i][0], coord[i][1])
            wkt += ')'
            geom = ogr.CreateGeometryFromWkt(wkt)
            feature = ogr.Feature(lyr.GetLayerDefn())
            feature.SetGeometry(geom)
            for key in dictfield :
                feature.SetField(key, dictfield[key][index])
            lyr.CreateFeature(feature)
            index += 1
    else:
        raise Exception('coords类型必须是list')

    ds = None
    driver =None

def createPolygon(shpname, coords, dictfield={}, layername=None, epsg=4326):
    '''
    利用GDAL 创建面矢量，将polygon转换成矢量数据

    Parameters
    ----------
    shpname : str
        输出文件名
    coords: list
        eg. [
                [[lon11, lat11], [lon12, lat12],... ],  # polygon1
                [[lon21, lat21], [lon22, lat22],... ],  # polygon2
            ]
    epsg: int

    Returns
    -------
        None
    '''

    driver = ogr.GetDriverByName('ESRI Shapefile')
    if os.access(shpname, os.F_OK) :
        driver.DeleteDataSource(shpname)

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(int(epsg))

    if layername is None :
        layername = os.path.basename(shpname).replace('.shp','')

    ds = driver.CreateDataSource(shpname)
    lyr = ds.CreateLayer(layername, srs, ogr.wkbPolygon)

    for key in dictfield :
        dtype = getOGRType(np.array(dictfield[key]).dtype)
        fieldName = ogr.FieldDefn(key, dtype)
        if dtype == ogr.OFTString :
            fieldName.SetWidth(200)
        lyr.CreateField(fieldName)

    index = 0
    if isinstance(coords, list) :
        # 逐个面对象创建
        for coord in coords :
            count = len(coord)
            wkt = 'POLYGON(('
            for i in range(count) :
                if i == count-1:
                    wkt += '%f %f' %(coord[i][0], coord[i][1])
                else:
                    wkt += '%f %f,' %(coord[i][0], coord[i][1])

            # 为了形成闭合圆环，最后一个点必须是第一个点位置
            if coord[count-1][0] != coord[0][0] and coord[count-1][1] != coord[0][1] :
                wkt += ',%f %f' %(coord[0][0], coord[0][1])

            wkt = wkt + '))'

            geom = ogr.CreateGeometryFromWkt(wkt)
            feature = ogr.Feature(lyr.GetLayerDefn())
            feature.SetGeometry(geom)
            for key in dictfield :
                feature.SetField(key, dictfield[key][index])
            lyr.CreateFeature(feature)
            index += 1
    else:
        raise Exception('coords类型必须是list')

    ds = None
    driver =None


def merge(outname, shplist):
    ''' 多个矢量文件合并 '''

    driver = ogr.GetDriverByName("ESRI Shapefile")
    ds_dst = driver.CreateDataSource(outname)
    spatialref = get_spatialref(shplist[0])
    geomtype = ogr.wkbPolygon

    basename = os.path.basename(outname)

    # 创建矢量面图层
    layer=ds_dst.CreateLayer(basename[:-4], srs=spatialref,
                             geom_type=geomtype, options=["ENCODING=UTF-8"])

    Field_Names = {}

    for shpname in shplist:
        ds_src = ogr.Open(shpname)
        Layer_src = ds_src.GetLayer()

        for field in Layer_src.schema :
            fieldname = field.name
            fieldtype = field.GetType()
            if not fieldname in Field_Names :
                Field_Names[fieldname] = fieldtype
        ds_src.Destroy()

    # 创建字段
    for fieldname in Field_Names :
        fieldName = ogr.FieldDefn(fieldname, Field_Names[fieldname])
        layer.CreateField(fieldName)

    for shpname in shplist:
        ds_src = ogr.Open(shpname)
        Layer_src = ds_src.GetLayer()

        for feature in Layer_src:

            geom = feature.GetGeometryRef()

            dst_feature = ogr.Feature(layer.GetLayerDefn())
            dst_feature.SetGeometry(geom)

            # # 定义输出属性表信息
            field_count = feature.GetFieldCount()
            for attr in range(field_count):
                field_defn = feature.GetFieldDefnRef(attr)
                fieldname = field_defn.GetName()
                # fieldtype = field_defn.GetType()

                # 将源文件中的字段信息写入匹配feature中
                dst_feature.SetField(fieldname, feature.GetField(fieldname))

            layer.CreateFeature(dst_feature)

        ds_src.Destroy()

    ds_dst.Destroy()


    # def merge1(self, outname, shplist):
    #     ''' 多个矢量文件合并 '''
    #
    #     driver = ogr.GetDriverByName("ESRI Shapefile")
    #     ds_dst = driver.CreateDataSource(outname)
    #     spatialref = self.get_spatialref(shplist[0])
    #     geomtype = ogr.wkbPolygon
    #
    #     basename = os.path.basename(outname)
    #
    #     # 创建矢量面图层
    #     layer=ds_dst.CreateLayer(basename[:-4], srs=spatialref,
    #                              geom_type=geomtype, options=["ENCODING=UTF-8"])
    #
    #     Field_Names = {}
    #
    #     errorInfo = {
    #         'Id' : 'classId',
    #         'classname' : 'class_name',
    #         'NAME' : 'class_name',
    #     }
    #
    #     for shpname in shplist:
    #         ds_src = ogr.Open(shpname)
    #         Layer_src = ds_src.GetLayer()
    #
    #         for field in Layer_src.schema :
    #             fieldname = field.name
    #             fieldtype = field.GetType()
    #             if not fieldname in ['class_name', 'classId'] :
    #                 continue
    #
    #             if not fieldname in Field_Names :
    #                 Field_Names[fieldname] = fieldtype
    #         ds_src.Destroy()
    #
    #     # 创建字段
    #     for fieldname in Field_Names :
    #         fieldName = ogr.FieldDefn(fieldname, Field_Names[fieldname])
    #         layer.CreateField(fieldName)
    #
    #     for shpname in shplist:
    #         ds_src = ogr.Open(shpname)
    #         Layer_src = ds_src.GetLayer()
    #
    #         for feature in Layer_src:
    #
    #             geom = feature.GetGeometryRef()
    #
    #             dst_feature = ogr.Feature(layer.GetLayerDefn())
    #             dst_feature.SetGeometry(geom)
    #
    #             # # 定义输出属性表信息
    #             field_count = feature.GetFieldCount()
    #             for attr in range(field_count):
    #                 field_defn = feature.GetFieldDefnRef(attr)
    #                 fieldname = field_defn.GetName()
    #                 # fieldtype = field_defn.GetType()
    #
    #                 if fieldname in errorInfo :
    #                     dst_feature.SetField(errorInfo[fieldname], feature.GetField(fieldname))
    #
    #                 if fieldname in Field_Names :
    #                     # 将源文件中的字段信息写入匹配feature中
    #                     dst_feature.SetField(fieldname, feature.GetField(fieldname))
    #
    #             layer.CreateFeature(dst_feature)
    #
    #         ds_src.Destroy()
    #
    #     ds_dst.Destroy()


def getField(srcfile, fieldname):

    fieldvalue = []
    driver = ogr.GetDriverByName('ESRI Shapefile')
    srcds = driver.Open(srcfile, gdalconst.GA_ReadOnly)
    src_layer = srcds.GetLayer()
    for feat in src_layer:
        value = feat.GetField(fieldname)
        fieldvalue.append(value)

    return fieldvalue

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

def fishgrid(outname, extent, xRes, yRes, epsg=4326, xindex=None, yindex=None):
    ''' https://blog.csdn.net/weixin_42924891/article/details/85865260
    https://zhuanlan.zhihu.com/p/403957133
    '''

    #初始化起始格网四角范围
    minX = extent[0]
    minY = extent[1]
    maxX = extent[2]
    maxY = extent[3]

    #计算行数和列数
    Height = int(np.ceil((maxY - minY) / yRes))
    Width = int(np.ceil((maxX - minX) / xRes))

    if xindex is not None and isinstance(xindex, str) :
        xindex = np.full(shape=(Width), fill_value=xindex, dtype='U256')

    if yindex is not None and isinstance(yindex, str) :
        yindex = np.full(shape=(Height), fill_value=yindex, dtype='U256')

    if xindex is not None and len(xindex) != Width :
        raise Exception('请检查输入xindex的维度，该维度长度须为%d' %(Width))

    if yindex is not None and len(yindex) != Height :
        raise Exception('请检查输入yindex的维度，该维度长度须为%d' %(Height))

    # 初始化投影
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(int(epsg))

    #创建输出文件
    outdriver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(outname):
        outdriver.DeleteDataSource(outname)
    outds = outdriver.CreateDataSource(outname)
    outlayer = outds.CreateLayer(outname, srs, geom_type = ogr.wkbPolygon)
    # 不添加属性信息，获取图层属性
    outfielddefn  = outlayer.GetLayerDefn()

    # 根据传入的xlabel和ylabel写入字段信息
    if xindex is None or yindex is None :
        for key in ['xindex', 'yindex'] :
            fieldName = ogr.FieldDefn(key, ogr.OFTInteger64)
            outlayer.CreateField(fieldName)
    else:
        for key in ['xindex', 'yindex'] :
            dtype = getOGRType(np.array(eval(key)).dtype)
            fieldName = ogr.FieldDefn(key, dtype)
            if dtype == ogr.OFTString :
                fieldName.SetWidth(200)
            outlayer.CreateField(fieldName)

    #遍历列，每一列写入格网
    for col in range(Width) :
        RectXLeft  = minX + xRes * col
        RectXRight = minX + xRes * (col+1)
        if RectXRight > maxX :
            RectXRight = maxX

        #遍历行，对这一列每一行格子创建和写入
        for row in range(Height) :
            RectYTop    = maxY - yRes * row
            RectYBottom = maxY - yRes * (row+1)
            if RectYBottom < minY :
                RectYBottom = minY

            #创建左上角第一个格子
            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(RectXLeft,  RectYTop)
            ring.AddPoint(RectXRight, RectYTop)
            ring.AddPoint(RectXRight, RectYBottom)
            ring.AddPoint(RectXLeft,  RectYBottom)
            ring.CloseRings()

            #写入几何多边形
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)
            #创建要素，写入多边形
            outfeat = ogr.Feature(outfielddefn)
            outfeat.SetGeometry(poly)
            if xindex is None or yindex is None :
                outfeat.SetField('xindex', col+1)
                outfeat.SetField('yindex', row+1)
            else:
                outfeat.SetField('xindex', xindex[col])
                outfeat.SetField('yindex', yindex[row])

            # 写入图层
            outlayer.CreateFeature(outfeat)

            outfeat = None
    #写入后清除缓存
    outds = None

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

def grid(srcshp, field, outname=None, format="GTiff", kind='idw',
         outputBounds=None,  # [minX, minY, maxX, maxY]
         xRes=None, yRes=None,
         power=2.0, smoothing=0.0, radius1=0.0, radius2=0.0,
         radius=0.0, angle=0.0, max_points=12, min_points=0,
         max_points_per_quadrant=0, min_points_per_quadrant=0,
         nodata=0.0, **kwargs):
    ''' https://gdal.org/programs/gdal_grid.html#interpolation-algorithms '''

    if kind.lower() in ['invdist', 'idw'] :
        algorithm = f'invdist:power={power}:smoothing={smoothing}' \
                    f':radius1={radius1}:radius2={radius2}:radius={radius}' \
                    f':angle={angle}:max_points={max_points}:min_points={min_points}' \
                    f':max_points_per_quadrant={max_points_per_quadrant}' \
                    f':min_points_per_quadrant={min_points_per_quadrant}' \
                    f':nodata={nodata}'
    elif kind.lower() in ['invdistnn'] :
        algorithm = f'invdistnn:power={power}:smoothing={smoothing}' \
                    f':radius={radius}' \
                    f':max_points={max_points}:min_points={min_points}' \
                    f':max_points_per_quadrant={max_points_per_quadrant}' \
                    f':min_points_per_quadrant={min_points_per_quadrant}' \
                    f':nodata={nodata}'
    elif kind.lower() in ['average', 'avg'] :
        algorithm = f'average:radius1={radius1}:radius2={radius2}:radius={radius}' \
                    f':angle={angle}:max_points={max_points}:min_points={min_points}' \
                    f':max_points_per_quadrant={max_points_per_quadrant}' \
                    f':min_points_per_quadrant={min_points_per_quadrant}' \
                    f':nodata={nodata}'
    elif kind.lower() in ['nearest'] :
        algorithm = f'nearest:radius1={radius1}:radius2={radius2}:radius={radius}' \
                    f':angle={angle}:nodata={nodata}'
    elif kind.lower() in ['linear'] :
        algorithm = f'linear:radius={radius}:nodata={nodata}'
    else:
        raise Exception('暂不支持该插值方法【%s】,请选择插值方案：“%s”' %(
            kind, 'invdist, invdistnn, average, nearest, linear'))

    vec = ogr.Open(srcshp)
    layer = vec.GetLayer()
    shpextent = layer.GetExtent()
    minX, maxX, minY, maxY = shpextent
    if outputBounds is None :
        outputBounds = [minX, minY, maxX, maxY]

    minX, minY, maxX, maxY = outputBounds
    if xRes is not None :
        width = int(np.ceil((maxX - minX) / xRes))
    else:
        width = 0

    if yRes is not None :
        height = int(np.ceil((maxY - minY) / yRes))
    else:
        height = 0

    if outname is None or format in ['MEM']:
        creationOptions=[]
        format = 'MEM'
    else:
        creationOptions=["COMPRESS=LZW"]


    opts = gdal.GridOptions(format=format, algorithm=algorithm, zfield=field,
                            outputBounds=outputBounds, width=width, height=height,
                            creationOptions=creationOptions, noData=nodata, **kwargs)
    ds = gdal.Grid(destName=outname, srcDS=srcshp, options=opts)

    if ds is not None :
        data = ds.ReadAsArray()
        trans = ds.GetGeoTransform()
        prj = ds.GetProjection()
        ds = None
        return data, trans, prj
    else:
        print('插值失败！')
        return None, None, None

def createBuffer(outname, srcshp, bdistance=None):
    """
    :param srcshp: 输入的矢量
    :param outname: 输出的矢量
    :param bdistance: 缓冲区距离
    :return:
    """

    in_ds = ogr.Open(srcshp)
    in_lyr = in_ds.GetLayer()

    if bdistance is None :
        extent = in_lyr.GetExtent()
        minX, maxX, minY, maxY = extent
        xdis = (maxX - minX) * 0.08
        ydis = (maxY - minY) * 0.08
        bdistance = np.nanmin([xdis, ydis])

    # 创建输出Buffer文件
    driver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(outname):
        driver.DeleteDataSource(outname)

    # 新建DataSource，Layer
    out_ds = driver.CreateDataSource(outname)
    out_lyr = out_ds.CreateLayer(outname, in_lyr.GetSpatialRef(), ogr.wkbPolygon, options=["ENCODING=UTF-8"])
    def_feature = out_lyr.GetLayerDefn()

    # 定义输出属性表信息
    feature = in_lyr.GetFeature(0)
    field_count = feature.GetFieldCount()
    field_names = []
    for attr in range(field_count):
        field_defn = feature.GetFieldDefnRef(attr)
        field_names.append(field_defn.GetName())
        out_lyr.CreateField(field_defn)

    # 遍历原始的Shapefile文件给每个Geometry做Buffer操作
    for feature in in_lyr:
        geometry = feature.GetGeometryRef()
        buffer = geometry.Buffer(bdistance)
        out_feature = ogr.Feature(def_feature)
        out_feature.SetGeometry(buffer)
        # 将源文件中的字段信息写入匹配feature中
        for field_name in field_names:
            out_feature.SetField(field_name, feature.GetField(field_name))
        out_lyr.CreateFeature(out_feature)
        out_feature = None

    out_ds.FlushCache()
    del in_ds, out_ds

def __createMultiPolygon_(strVectorFile):
    strDriverName = "ESRI Shapefile"  # 创建数据，这里创建ESRI的shp文件
    oDriver = ogr.GetDriverByName(strDriverName)
    if oDriver == None:
        print("%s 驱动不可用！\n", strDriverName)

    oDS = oDriver.CreateDataSource(strVectorFile)  # 创建数据源
    if oDS == None:
        print("创建文件【%s】失败！", strVectorFile)

    srs = osr.SpatialReference()  # 创建空间参考
    srs.ImportFromEPSG(4326)  # 定义地理坐标系WGS1984
    papszLCO = []
    # 创建图层，创建一个多边形图层,"TestPolygon"->属性表名
    oLayer = oDS.CreateLayer("TestPolygon", srs, ogr.wkbPolygon, papszLCO)
    if oLayer == None:
        print("图层创建失败！\n")

    '''下面添加矢量数据，属性表数据、矢量数据坐标'''
    oFieldID = ogr.FieldDefn("FieldID", ogr.OFTInteger)  # 创建一个叫FieldID的整型属性
    oLayer.CreateField(oFieldID, 1)

    oFieldName = ogr.FieldDefn("FieldName", ogr.OFTString)  # 创建一个叫FieldName的字符型属性
    oFieldName.SetWidth(100)  # 定义字符长度为100
    oLayer.CreateField(oFieldName, 1)

    oDefn = oLayer.GetLayerDefn()  # 定义要素

    # 创建单个面
    oFeatureTriangle = ogr.Feature(oDefn)
    oFeatureTriangle.SetField(0, 0)  # 第一个参数表示第几个字段，第二个参数表示字段的值
    oFeatureTriangle.SetField(1, "单个面")
    ring = ogr.Geometry(ogr.wkbLinearRing)  #  构建几何类型:线
    ring.AddPoint(0, 0)  #  添加点01
    ring.AddPoint(10, 0)  #  添加点02
    ring.AddPoint(10, 10)  #  添加点03
    ring.AddPoint(0, 10)  #  添加点04
    yard = ogr.Geometry(ogr.wkbPolygon)  #  构建几何类型:多边形
    yard.AddGeometry(ring)
    yard.CloseRings()

    geomTriangle = ogr.CreateGeometryFromWkt(str(yard))  # 将封闭后的多边形集添加到属性表
    oFeatureTriangle.SetGeometry(geomTriangle)
    oLayer.CreateFeature(oFeatureTriangle)

    # 第一种方法创建多个面，注意如果两个面有重叠，出现镂空显示
    oFeatureTriangle = ogr.Feature(oDefn)
    oFeatureTriangle.SetField(0, 1)  # 第一个参数表示第几个字段，第二个参数表示字段的值
    oFeatureTriangle.SetField(1, "多个面")
    box1 = ogr.Geometry(ogr.wkbLinearRing)
    box1.AddPoint(15, 0)
    box1.AddPoint(25, 0)
    box1.AddPoint(20, 10)
    garden1 = ogr.Geometry(ogr.wkbPolygon)
    garden1.AddGeometry(box1)

    box2 = ogr.Geometry(ogr.wkbLinearRing)
    box2.AddPoint(30, 0)
    box2.AddPoint(40, 0)
    box2.AddPoint(35, 10)
    garden2 = ogr.Geometry(ogr.wkbPolygon)
    garden2.AddGeometry(box2)

    box3 = ogr.Geometry(ogr.wkbLinearRing)
    box3.AddPoint(32, 1)
    box3.AddPoint(38, 1)
    box3.AddPoint(35, 7)
    garden3 = ogr.Geometry(ogr.wkbPolygon)
    garden3.AddGeometry(box3)

    gardens = ogr.Geometry(ogr.wkbMultiPolygon)
    gardens.AddGeometry(garden1) #分别将三个多边形面添加到总多变形面中
    gardens.AddGeometry(garden2)
    gardens.AddGeometry(garden3)
    gardens.CloseRings()

    geomTriangle = ogr.CreateGeometryFromWkt(str(gardens))  # 将封闭后的多边形集添加到属性表
    oFeatureTriangle.SetGeometry(geomTriangle)
    oLayer.CreateFeature(oFeatureTriangle)

    # 第二种方法创建多个面，注意如果两个面有重叠，出现镂空显示
    oFeatureTriangle = ogr.Feature(oDefn)
    oFeatureTriangle.SetField(0, 2)  # 第一个参数表示第几个字段，第二个参数表示字段的值
    oFeatureTriangle.SetField(1, "多个面")

    lot = ogr.Geometry(ogr.wkbLinearRing)
    lot.AddPoint(58, 38.5)
    lot.AddPoint(53, 6)
    lot.AddPoint(99.5, 19)
    lot.AddPoint(73, 42)

    house = ogr.Geometry(ogr.wkbLinearRing)
    house.AddPoint(67.5, 29)
    house.AddPoint(69, 25.5)
    house.AddPoint(64, 23)
    house.AddPoint(69, 15)
    house.AddPoint(82.5, 22)
    house.AddPoint(76, 31.5)

    yard = ogr.Geometry(ogr.wkbPolygon)
    yard.AddGeometry(lot) #分别将两个多边形线添加到总多边形面中
    yard.AddGeometry(house)
    yard.CloseRings()

    geomTriangle = ogr.CreateGeometryFromWkt(str(yard))  # 将封闭后的多边形集添加到属性表
    oFeatureTriangle.SetGeometry(geomTriangle)
    oLayer.CreateFeature(oFeatureTriangle)

    oDS.Destroy()
    print("数据集创建完成！\n")


def __polygon2polyline(polyfn, linefn):
    """
        This function is used to make polygon convert to line
    :param polyfn: the path of input, the shapefile of polygon
    :param linefn: the path of output, the shapefile of line
    :return:
    """
    driver = ogr.GetDriverByName('ESRI Shapefile')
    polyds = ogr.Open(polyfn, 0)
    polyLayer = polyds.GetLayer()
    spatialref = polyLayer.GetSpatialRef()
    #创建输出文件
    if os.path.exists(linefn):
        driver.DeleteDataSource(linefn)
    lineds =driver.CreateDataSource(linefn)
    linelayer = lineds.CreateLayer(linefn, srs=spatialref, geom_type=ogr.wkbLineString)
    featuredefn = linelayer.GetLayerDefn()
    #获取ring到几何体
    #geomline = ogr.Geometry(ogr.wkbGeometryCollection)
    for feat in polyLayer:
        geom = feat.GetGeometryRef()
        ring = geom.GetGeometryRef(0)
        #geomcoll.AddGeometry(ring)
        outfeature = ogr.Feature(featuredefn)
        outfeature.SetGeometry(ring)
        linelayer.CreateFeature(outfeature)
        outfeature = None


def __smoothPolygon(_input_shapefile, _smoothed_shapefile_output_path, _buffer_distance=0.0005):
    '''平滑矢量文件边界
    _input_shapefile：输入的矢量文件路径
    _smoothed_shapefile_output_path：平滑后输出的矢量文件路径
    _buffer_distance：平滑缓冲区距离，单位是米，可以根据影像分辨率设定，不能低于影像分辨率'''
    in_ds = ogr.Open(_input_shapefile)
    in_lyr = in_ds.GetLayer()
    # feature_number = in_lyr.GetFeatureCount()
    driver = ogr.GetDriverByName('ESRI Shapefile')
    if Path(_smoothed_shapefile_output_path).exists():
        # 如果输出的文件存在，则删除这个文件
        driver.DeleteDataSource(_smoothed_shapefile_output_path)
    out_ds = driver.CreateDataSource(_smoothed_shapefile_output_path)
    out_lyr = out_ds.CreateLayer(_smoothed_shapefile_output_path, in_lyr.GetSpatialRef(), ogr.wkbPolygon)
    def_feature = out_lyr.GetLayerDefn()

    for feature in in_lyr:
        geometry = feature.GetGeometryRef()
        _buffer = geometry.Buffer(_buffer_distance).Buffer(-1 * _buffer_distance)
        out_feature = ogr.Feature(def_feature)
        out_feature.SetGeometry(_buffer)
        out_lyr.CreateFeature(out_feature)
        out_feature = None
    out_ds.FlushCache()
    in_ds.Destroy()
    out_ds.Destroy()


# def contourGenerate(dstLayer, srcBand, contourInterval, contourBase,
#                     fixedLevelCount, useNoData, noDataValue,  idField, elevField):
#     gdal.ContourGenerate()

def contour(srcfile, dstfile, attrname, interval, offset=0, band=1,
            srcNoData=None, polygons=False, format='ESRI Shapefile'):
    ''' 创建等值线 '''

    #  The gdal_contour generates a vector contour file
    #      from the input raster elevation model (DEM).
    #     The contour line-strings are oriented consistently
    #     and the high side will be on the right,
    #     i.e. a line string goes clockwise around a top.
    #

    # gdal_contour [-b <band>] [-a <attribute_name>] [-amin <attribute_name>] [-amax <attribute_name>]
    # [-3d] [-inodata]
    # [-snodata n] [-i <interval>]
    # [-f <formatname>] [[-dsco NAME=VALUE] ...] [[-lco NAME=VALUE] ...]
    # [-off <offset>] [-fl <level> <level>...] [-e <exp_base>]
    # [-nln <outlayername>] [-q] [-p]
    # <src_filename> <dst_filename>

    cmd = 'gdal_contour -f "{format}" {srcfile} {dstfile} ' \
          '-b {band} -a {attrname} -i {interval} -off {offset}'.format(
        format=format,
        srcfile = srcfile,
        dstfile = dstfile,
        band=band,
        attrname=attrname,
        interval=interval,
        offset=offset
    )
    if srcNoData is not None :
        cmd += f' -snodata {srcNoData}'

    if polygons :
        cmd += ' -p'
    print('执行下载命令: 【%s】' %(cmd))
    status = os.system(cmd)
    print(status)

def getGDALType(datatype):
    '''
    根据numpy的数据类型，匹配GDAL中的数据类型

    Parameters
    ----------
    datatype

    Returns
    -------
        GDAL数据类型
    '''

    if datatype == np.byte or datatype == np.uint8:
        return gdal.GDT_Byte
    elif datatype == np.uint16 :
        return gdal.GDT_UInt16
    elif datatype == np.int16 :
        return gdal.GDT_Int16
    elif datatype == np.uint32 :
        return gdal.GDT_UInt32
    elif datatype == np.int32 :
        return gdal.GDT_Int32
    elif datatype == np.float32 :
        return gdal.GDT_Float32
    elif datatype == np.float64 :
        return gdal.GDT_Float64
    else:
        return gdal.GDT_Unknown


def getOGRType(dtype) :
    '''
    根据numpy的数据类型，匹配OGR中的数据类型

    Parameters
    ----------
    dtype

    Returns
    -------
        OGR数据类型
    '''

    if dtype == np.byte or dtype == np.uint8 or 'U' in str(dtype):
        return ogr.OFTString
    elif dtype == np.uint16 or dtype == np.int16 :
        return ogr.OFSTInt16
    elif dtype == np.uint32 or dtype == np.int32  :
        return ogr.OFTInteger
    elif dtype == np.int64 :
        return ogr.OFTInteger64
    elif dtype == np.float32 :
        return ogr.OFSTFloat32
    elif dtype == np.float64 :
        return ogr.OFTReal
    else:
        return ogr.OJUndefined

def getNPDType(datatype) :
    '''将numpy.dtype转换成 gdal Type'''
    if 'int' in datatype:
        return gdal.GDT_Int32
    elif 'float' in datatype:
        return gdal.GDT_Float32
    else:
        return gdal.GDT_Unknown
