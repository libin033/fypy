# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : sentinelscene.py

@Modify Time :  2024/4/15   

@Author      : Lee    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime

from tqdm import tqdm
import shutil
import json
from osgeo import gdal, osr, ogr
from dateutil.relativedelta import relativedelta
from fypy.tools.BaseAlgorithms import BaseAlgorithms
from fypy import parm
import xml.dom.minidom

EXEPATH = os.path.abspath(list(parm.__path__)[0])


class SentinelScene(BaseAlgorithms) :

    def __init__(self):
        pass

    def calibrate(self, srcfile, BandId, fillvalue=65535, blocksize=1000):
        ''' 辐射定标 '''

        srcdset = gdal.Open(srcfile, gdal.GA_ReadOnly)
        if srcdset == None:
            print("文件无法打开【%s】" % (srcfile))
            return None

        # 获取行列号
        cols = srcdset.RasterXSize                  # 栅格矩阵的列数
        rows = srcdset.RasterYSize                  # 栅格矩阵的行数
        Bands = srcdset.RasterCount                 # 波段数

        srcTrans = srcdset.GetGeoTransform()   # 获取仿射矩阵信息
        srcProj = srcdset.GetProjection()             # 获取投影信息

        # 创建输出结果文件
        driver = gdal.GetDriverByName("MEM")
        dstdataset = driver.Create('', cols, rows, Bands, gdal.GDT_Float32)
        dstdataset.SetGeoTransform(srcTrans)
        dstdataset.SetProjection(srcProj)

        ReadBand = srcdset.GetRasterBand(1)
        outband = dstdataset.GetRasterBand(1)
        outband.SetNoDataValue(fillvalue)

        #进度条参数
        XBlockcount = np.ceil(cols / blocksize)
        YBlockcount = np.ceil(rows / blocksize)
        i = 0
        j = 0
        try:
            with tqdm(total=XBlockcount*YBlockcount, iterable='iterable',
                      desc = '正在进行第%i波段定标' %(BandId), mininterval=1) as pbar:
                while i < rows:
                    while j < cols:
                        # 保存分块大小
                        nXBK = blocksize
                        nYBK = blocksize

                        # 最后一块
                        if i+blocksize > rows:
                            nYBK = rows - i
                        if j+blocksize > cols:
                            nXBK=cols - j

                        # 分块读取影像
                        srcdata = ReadBand.ReadAsArray(j, i, nXBK,nYBK)
                        TOARef =np.where(srcdata>0, (srcdata-1000) / 10000.0, fillvalue)
                        outband.WriteArray(TOARef, j, i)
                        j=j+nXBK
                        pbar.update(1)
                    j=0
                    i=i+nYBK
        except BaseException :
            pbar.close()

        pbar.close()

    def getMetaInfo(self, metafile):
        dom = xml.dom.minidom.parse(metafile)
        dict_meta = {}

        #太阳天顶角、方位角
        SunAngle = dom.getElementsByTagName('Mean_Sun_Angle')
        dict_meta['sunz'] = float(SunAngle[0].getElementsByTagName('ZENITH_ANGLE')[0].firstChild.data)
        dict_meta['suna'] = float(SunAngle[0].getElementsByTagName('AZIMUTH_ANGLE')[0].firstChild.data)

        #卫星天顶角、方位角
        ViewAngles = dom.getElementsByTagName('Mean_Viewing_Incidence_Angle')
        dict_satz = {}
        dict_sata = {}

        for angle in ViewAngles:
            ViewAngle = int(angle.getAttribute('bandId'))
            dict_satz[ViewAngle+1] = float(angle.getElementsByTagName('ZENITH_ANGLE')[0].firstChild.data)
            dict_sata[ViewAngle+1]= float(angle.getElementsByTagName('AZIMUTH_ANGLE')[0].firstChild.data)

        dict_meta["ViewZenithAngle"] = dict_satz
        dict_meta["ViewAzimuthAngle"] = dict_sata
        # 日期:月、日  2023-04-11T03:50:21.042436Z
        dict_meta["Date"] = datetime.datetime.strptime(dom.getElementsByTagName('SENSING_TIME')[0].firstChild.data,
                                                       '%Y-%m-%dT%H:%M:%S.%fZ')

        #求影像中心经纬度
        PointULX = int(dom.getElementsByTagName('ULX')[0].firstChild.data)
        PointULY = int(dom.getElementsByTagName('ULY')[0].firstChild.data)

        Imgsizes = dom.getElementsByTagName('Size')

        for Imgsize in Imgsizes:
            Resolution = Imgsize.getAttribute('resolution')
            if Resolution == '10':
                dict_meta["Nrows"] = int(Imgsize.getElementsByTagName('NROWS')[0].firstChild.data)
                dict_meta["Ncols"] = int(Imgsize.getElementsByTagName('NCOLS')[0].firstChild.data)

        PointBRX = PointULX + 10*dict_meta["Ncols"]
        PointBRY = PointULY - 10*dict_meta["Nrows"]

        # 将投影坐标转为经纬度坐标（具体的投影坐标系由给定数据确定）
        Proj = dom.getElementsByTagName('HORIZONTAL_CS_CODE')[0].firstChild.data
        ProjCode = int(Proj.split(':')[1])

        source = osr.SpatialReference()
        source.ImportFromEPSG(ProjCode)
        target = osr.SpatialReference()
        target.ImportFromEPSG(4326)
        ct = osr.CoordinateTransformation(source,target)
        CoordsUL,CoordsBR = ct.TransformPoints([(PointULX,PointULY),(PointBRX,PointBRY)])

        ULLat = CoordsUL[0]
        ULLon = CoordsUL[1]
        BRLat = CoordsBR[0]
        BRLon = CoordsBR[1]

        sLongitude = (ULLon+BRLon) / 2
        sLatitude = (ULLat+BRLat) / 2

        dict_meta['extent'] = [ULLon, BRLat, BRLon, ULLat]  # minX, minY, maxX, maxY
        dict_meta['centre_pos'] = [sLongitude, sLatitude]

        return dict_meta


def getmw(wl, resp) :
    resp = np.array(resp)
    newwl = np.arange(wl[0], wl[1]+0.0002, 0.0025)

    index = np.argmin(1-resp)

    return newwl[index]

if __name__ == '__main__':

    from fypy.tools.jsonpro import readjson, writejson
    ccfile  = r'D:\pypi\fypy\fypy\RSDP\resp\bak\Sentinel\CalibrationCoefficient.json'
    rsffile = r'D:\pypi\fypy\fypy\RSDP\resp\bak\Sentinel\Sentinel.json'

    data2 = readjson(ccfile)
    data1 = readjson(rsffile)

    for satid in data1 :
        for instid in data1[satid] :
            dict_info = {}
            info1 = data1[satid][instid]

            for item in info1 :
                if item in ['8A'] :
                    continue

                bandid = int(item.replace('B', ''))

                print(satid, instid, bandid)
                dict_info['B%02d' %(bandid)] = {}

                mw = getmw(info1[item]['wl'], info1[item]['SRF'])

                dict_info['B%02d' %(bandid)]['MW'] = float('%.3f' %(mw))

                wl = list([float('%.4f' %(info1[item]['wl'][0])), float('%.4f' %(info1[item]['wl'][1]))])
                dict_info['B%02d' %(bandid)]['wl'] = wl
                dict_info['B%02d' %(bandid)]['SRF'] = info1[item]['SRF']

                if not satid in data2 :
                    continue

                info2 = data2[satid][instid]
                if 'ESUN' in info2 :
                    if '8A' in data1[satid][instid] and bandid > 8 :
                        sun = info2['ESUN'][bandid]
                        dict_info['B%02d' %(bandid)]['ESUN'] = float('%.4f' %(sun))
                    else:
                        sun = info2['ESUN'][bandid-1]

                        dict_info['B%02d' %(bandid)]['ESUN'] = float('%.4f' %(sun))


            outjson = '%s_%s.json' %(satid.upper(), instid.upper())
            writejson(outjson, dict_info)




