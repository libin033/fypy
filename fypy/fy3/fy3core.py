# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : fy3core.py

@Modify Time :  2024/1/26 16:10   

@Author      : fypy Team    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime


from osgeo import gdal
from fypy.tools.tifpro import getDateSet
from fypy.tools.hdfpro import readhdf

from .fy3config import FY3Block10CoefX, FY3Block10CoefY, Prj_Info, FY3ProdInfo


def calref(dn, coef):
    ''' 可见光通道定标 '''

    ref = np.full_like(dn, fill_value=np.nan, dtype=np.float32)
    for i in range(dn.shape[0]) :
        ref[i] = (coef[i, 0] + coef[i, 1] * dn[i] + coef[i, 2] * dn[i] * dn[i]) * 0.01 # 将反射率转为0~1

        ref[i, dn[i] == 65535] = np.nan
        ref[i, ref[i]<=0] = 0

    return ref

def calemiss(dn, wavenum):
    ''' 红外波段定标 '''

    temp = dn * 0.01
    temp[(temp<=0) | (temp >= 600.0)] = np.nan

    bt = np.full_like(temp, fill_value=np.nan, dtype=np.float32)

    for i in range(temp.shape[0]) :
        bt[i] = planck_r2t(temp[i, :, :], wavenum[i])

    return bt

def planck_r2t(rad, wn):
    '''
    普朗克函数：将辐射值转成亮温（K）

    Parameters
    ----------
    rad : numpy.narray
        mW/(m2.cm-1.sr)
    wn : float or numpy.narray
        wave number(cm^-1)

    Returns
    -------
        numpy.narray
        bright temperature
        units: K
    '''

    # 普朗克系数
    Radiation_C1 = 1.191042953E-5
    Radiation_C2 = 1.4387774

    bt = (Radiation_C2 * wn / np.log(Radiation_C1 * wn * wn * wn / (rad)+1.0))

    return bt

def GetSourceInfo(filename, sdsname):
    src_ds = gdal.Open(filename)
    layers = src_ds.GetSubDatasets()

    # 获取sdsname所在的图层栅格索引
    src_raster = GetLayer(layers, sdsname)
    if src_raster is None :
        raise Exception('数据集【%s】不在文件中【%s】' %(sdsname, filename))

    return src_raster

def GetLayer(layers, sdsname):
    '''
    获取指定的图层的索引名
    :param layers: tuple
    :return: str
    '''

    if sdsname:
        for layer in layers :
            l_name = layer[0].split(':')[-1].replace('"','')

            if sdsname in l_name:
                return layer[0]

    return None

def CreateVrt(vrtDir, srcfile, layer, srclon=None, srclat=None):

    if layer is None :
        gdal.Translate(vrtDir,
                       srcfile,
                       format='vrt')
    else:
        gdal.Translate(vrtDir,
                       layer,
                       format='vrt')

    lines = []
    with open(vrtDir, 'r') as f:
        for line in f:
            lines.append(line)
    lines.insert(1,'<Metadata domain="GEOLOCATION">\n')
    lines.insert(2,' <MDI key="LINE_OFFSET">1</MDI>\n')
    lines.insert(3, ' <MDI key="LINE_STEP">1</MDI>\n')
    lines.insert(4, ' <MDI key="PIXEL_OFFSET">1</MDI>\n')
    lines.insert(5, ' <MDI key="PIXEL_STEP">1</MDI>\n')
    lines.insert(6, ' <MDI key="SRS">GEOGCS["WGS84",'
                    'DATUM["WGS_1984",'
                    'SPHEROID["WGS84",6378137,298.257223563,'
                    'AUTHORITY["EPSG","7030"]],'
                    'TOWGS84[0,0,0,0,0,0,0],'
                    'AUTHORITY["EPSG","6326"]],'
                    'PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],'
                    'UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9108"]],'
                    'AUTHORITY["EPSG","4326"]]</MDI>\n')
    lines.insert(7, ' <MDI key="X_BAND">1</MDI>')
    lines.insert(8, ' <MDI key="X_DATASET">HDF5:"{geofile}":/{srclon}</MDI>\n'.format(
        geofile=srcfile, srclon=srclon))
    lines.insert(9, ' <MDI key="Y_BAND">1</MDI>\n')
    lines.insert(10, ' <MDI key="Y_DATASET">HDF5:"{geofile}":/{srclat}</MDI>\n'.format(
        geofile=srcfile, srclat=srclat))
    lines.insert(11, '</Metadata>\n')
    with open(vrtDir, 'w') as f:
        for line in lines:
            f.writelines(line)

def GetNameInfo(filename):
    ''' 根据输入数据的文件名获取信息 '''
    # FY3D_MERSI_GBAL_L1_YYYYMMDD_HHmm_GEO1K_MS.HDF
    # FY3D_MERSI_GBAL_L2_PWV_MLT_GLL_YYYYMMDD_POAD_5000M_MS.HDF
    # FY3D_MERSI_GBAL_L3_NVI_MLT_GLL_YYYYMMDD_AOAM_5000M_MS.HDF

    nameinfo = {}
    basename = os.path.basename(filename)
    if len(basename) != 45 and len(basename) != 57 :
        print('输入文件名为非标准文件名，将不作信息提取【%s】' %(basename))
        return nameinfo

    namelist = basename.split('_')
    nameinfo['SatID']  = namelist[0]
    nameinfo['InstID'] = namelist[1]
    nameinfo['RegionID'] = namelist[2]
    nameinfo['LevelID'] = namelist[3]
    if 'L1' in namelist[3] :
        nameinfo['StartTime'] = datetime.datetime.strptime('%s %s' %(namelist[4], namelist[5]),
                                                           '%Y%m%d %H%M')
        if '1000M' in namelist[6] or '1K' in namelist[6] :
            nameinfo['Resolution'] = 0.01
        elif '250M' in namelist[6] or 'QK' in namelist[6] :
            nameinfo['Resolution'] = 0.0025
    else:
        nameinfo['ProdID'] = namelist[4]
        nameinfo['Proj'] = namelist[6]
        nameinfo['StartTime'] = datetime.datetime.strptime(namelist[7], '%Y%m%d')
        nameinfo['CType'] = namelist[8]
        if 'KM' in namelist[9] :
            nameinfo['Resolution'] = float(namelist[9].replace('KM',''))/100.0
        elif 'M' in namelist[9] :
            nameinfo['Resolution'] = float(namelist[9].replace('M',''))/100.0/1000.0

    return nameinfo

def get_data(filename, ProdID, sdsname, SatID=None, InstID=None):

    if SatID is None or InstID is None :
        nameinfo = GetNameInfo(filename)
        if len(nameinfo) == 0 :
            raise Exception('请提供参数【SatID】和【InstID】')
        SatID = nameinfo['SatID']
        InstID = nameinfo['InstID']

    try:
        if ProdID in FY3ProdInfo[SatID][InstID]:
            val, sdsinfo = readhdf(filename, sdsname, dictsdsinfo={})
            if 'FillValue' in sdsinfo :

                if isinstance(sdsinfo['FillValue'], np.ndarray) :
                    fillvalue = sdsinfo['FillValue'][0]
                else:
                    fillvalue = sdsinfo['FillValue']
            else:
                fillvalue = 65535

            if 'valid_range' in sdsinfo :
                vmin = sdsinfo['valid_range'][0]
                vmax = sdsinfo['valid_range'][1]
                fillflag = (val<vmin) | (val>vmax)

            if 'Slope' in sdsinfo and 'Intercept' in sdsinfo :
                val = val * sdsinfo['Slope'] + sdsinfo['Intercept']

            val[fillflag] = fillvalue
            # if 'valid_range' in sdsinfo :
            #     vmin = sdsinfo['valid_range'][0]
            #     vmax = sdsinfo['valid_range'][1]
            # fillflag = (val<vmin) | (val>vmax)
            # val[fillflag] = self.fillvalue

            if ProdID == 'CLM':
                data_u8 = np.uint8(val)
                data_bit = np.unpackbits(data_u8).reshape(data_u8.size, 8)
                bit21 = data_bit[:, 5] * 4 + data_bit[:, 6] * 2 + data_bit[:, 7]
                val = bit21.reshape(val.shape)
                # 001 = Cloudy （1）
                # 011 = Uncertain （3）
                # 101 = Probably  Clear （5）
                # 111 = Confident  Clear （7）
                val[(val != 1) & (val != 3) & (val != 5) & (val != 7)] = 255
                val[val==1] = 0
                val[val==3] = 1
                val[val==5] = 2
                val[val==7] = 3
                fillvalue = 255

            key = sdsname.replace(' ', '_')
            if '/' in key :
                key = key.split('/')[-1]

            return val
        else:
            print('产品【%s】不在绘图要求列表之内' %(ProdID))
            return None
    except BaseException as  e :
        print('读取失败【%s】：%s' %(ProdID, filename))
        return None

def single_files(filename, productName, sdsname, fillvalue) :
    '''单个数据产品投影拼接'''

    NameInfo = GetNameInfo(filename)
    if len(NameInfo) == 0 :
        raise Exception('请传入官方标准文件名格式')

    resolution = NameInfo['Resolution']

    prjType = NameInfo['Proj']
    if prjType in Prj_Info :
        prjwkt = Prj_Info[prjType]
    else:
        prjwkt = Prj_Info['GLL']

    blockID = NameInfo['RegionID']

    lat = FY3Block10CoefY[blockID[0:2]]
    lon = FY3Block10CoefX[blockID[2:4]]

    mtrans = (float(lon) * 100000, float(resolution) * 100000, 0,
              float(lat) * 100000, 0, -1 * float(resolution) * 100000)

    data = get_data(filename, productName, sdsname)
    ds   = getDateSet(data, mtrans, prj=prjwkt, fillvalue=fillvalue)

    return ds


def xy2latlon(x_s, y_s, x_res, y_res, rows, cols):
    ''' hammer转为WGS84坐标 '''
    # z^2 = 1 - x^2/2 - y^2/2
    # longitude = 2 * atan(sqrt(2) * x * z / (2 * z^2 - 1))
    # latitude = asin(sqrt(2) * y * z)

    x, y = np.meshgrid(np.linspace(x_s, x_s+x_res*(cols-1), num=cols),
                       np.linspace(y_s, y_s+y_res*(rows-1), num=rows))

    # 将degree转为meter: 100.0 * 1000.0
    x = np.where(x > (180.0 * 100.0 * 1000.0), (180.0 * 100.0 * 1000.0) - x, x)
    x = x / (180.0 * 100.0 * 1000.0)

    y = np.where(y > (90.0 * 100.0 * 1000.0), (90.0 * 100.0 * 1000.0) - y, y)
    y = y / (9000.0 * 1000.0)

    z = np.sqrt(1 - np.square(x) / 2.0 - np.square(y) / 2.0)
    lon = 2 * np.arctan(np.sqrt(2) * x * z / (2.0 * (np.square(z)) - 1))
    lat = np.arcsin(np.sqrt(2) * y * z)

    # 弧度转为度
    lon = lon / np.pi * 180.0
    lat = lat / np.pi * 180.0

    return lon, lat

def latlon2xy(x, y, Pixel, Line, resolution):

    lon, lat = np.meshgrid(np.linspace(x, x+resolution*(Pixel-1), num=Pixel),
                           np.linspace(y,y-1 * resolution*(Line-1), num=Line))

    # 将经纬度坐标转化为Hammer坐标
    lon = np.where(lon > 180.0, lon - 360.0, lon)
    lon = lon / 180.0 * np.pi
    lat = lat / 180.0 * np.pi

    newz = np.sqrt(1 + np.cos(lat) * np.cos(lon / 2.0))
    x = np.cos(lat) * np.sin(lon / 2.0) / newz
    y = np.sin(lat) / newz

    del lat, lon

    x = x * (180.0 * 100.0 * 1000.0)
    y = y * (90.0  * 100.0 * 1000.0)

    return x, y



def get_bounding_box(filelist):

    tileindex = []

    for filename in filelist :
        basename = os.path.basename(filename)
        names = basename.split('_')
        blockID = names[2]

        lat = FY3Block10CoefY[blockID[0:2]]
        lon = FY3Block10CoefX[blockID[2:4]]

        tileindex.append([lat, lon])

    tileindex = np.array(tileindex)
    maxY = np.nanmax(tileindex[:,0])
    minY = np.nanmin(tileindex[:,0])

    maxX = np.nanmax(tileindex[:,1])
    minX = np.nanmin(tileindex[:,1])

    extent = [minX, maxX+10, minY, maxY+10]

    return extent




