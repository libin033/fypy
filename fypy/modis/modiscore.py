# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : modiscore.py

@Modify Time :  2024/2/1 16:05   

@Author      : fypy Team    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime

from osgeo import gdal

def GetSourceInfo(filename, sdsname):

    src_ds = gdal.Open(filename)
    layers = src_ds.GetSubDatasets()

    # 获取sdsname所在的图层栅格索引
    if sdsname:
        for layer in layers :
            l_name = layer[0].split(':')[-1].replace('"','')

            if sdsname in l_name:
                return layer[0]

    raise Exception('数据集【%s】不在文件中【%s】' %(sdsname, filename))



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
    lines.insert(6, ' <MDI key="SRS">GEOGCS["WGS 84",'
                    'DATUM["WGS_1984",'
                    'SPHEROID["WGS84",6378137,298.257223563,'
                    'AUTHORITY["EPSG","7030"]],'
                    'AUTHORITY["EPSG","6326"]],'
                    'PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],'
                    'UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],'
                    'AXIS["Latitude",NORTH],AXIS["Longitude",EAST],'
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



