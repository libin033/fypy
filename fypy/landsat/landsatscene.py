# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : landsatscene.py

@Modify Time :  2024/4/15   

@Author      : Lee    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime
from osgeo import gdal
from tqdm import tqdm
import tarfile
from dateutil.relativedelta import relativedelta
from fypy.tools.BaseAlgorithms import BaseAlgorithms
from fypy.tools.tifpro import readtiff
from fypy import parm

EXEPATH = os.path.abspath(list(parm.__path__)[0])


class LandsatScene(BaseAlgorithms) :

    def __init__(self):
        pass

    def calibrate(self, srcfile, metafile, BandId=None, fillvalue=65535, blocksize=1000):

        self.metadata = self.getMetaInfo(metafile)

        if BandId is None :
            BandId = int((os.path.basename(srcfile).split('.')[0])[-1])

        srcdset = gdal.Open(srcfile, gdal.GA_ReadOnly)
        if srcdset == None:
            print("文件无法打开【%s】" % (srcfile))
            return None

        cols = srcdset.RasterXSize                  # 栅格矩阵的列数
        rows = srcdset.RasterYSize                  # 栅格矩阵的行数
        Bands = srcdset.RasterCount                 # 波段数

        trans = srcdset.GetGeoTransform()   # 获取仿射矩阵信息
        prj = srcdset.GetProjection()             # 获取投影信息

        # 创建输出结果文件
        driver = gdal.GetDriverByName("MEM")
        dstdset = driver.Create('', cols, rows, Bands, gdal.GDT_Float32)
        # options=["COMPRESS=LZW", "BigTIFF=YES"])
        dstdset.SetGeoTransform(trans)
        dstdset.SetProjection(prj)

        ReadBand = srcdset.GetRasterBand(1)
        outband = dstdset.GetRasterBand(1)
        outband.SetNoDataValue(fillvalue)

        i = 0
        j = 0
        #进度条参数
        XBlockcount = np.ceil(cols / blocksize)
        YBlockcount = np.ceil(rows / blocksize)
        try:
            with tqdm(total=XBlockcount*YBlockcount, iterable='iterable',
                      desc = '正在进行第%i波段定标' %(BandId), mininterval=1) as pbar:
                while i < rows:
                    while j < cols:
                        # 保存分块大小
                        nXBK = blocksize
                        nYBK = blocksize

                        # 最后一块
                        if i+blocksize>rows:
                            nYBK = rows - i
                        if j+blocksize>cols:
                            nXBK=cols - j

                        # 分块读取影像
                        DN = ReadBand.ReadAsArray(j, i, nXBK,nYBK)

                        #辐射校正
                        # radiance = self.GetRadiance(BandId, ImgRasterData, self.metadata)
                        data = self.GetReflectance(BandId, DN, self.metadata, fillvalue=fillvalue)

                        outband.WriteArray(data,j,i)
                        j=j+nXBK
                        pbar.update(1)
                    j=0
                    i=i+nYBK
        except BaseException :
            pbar.close()
        pbar.close()

        srcdset = None

        return dstdset

    def combine(self, outdir, filename) :
        ''' 多通道合成 '''

        basename = os.path.basename(filename)
        basename = basename.replace('.tar', '')
        untardir = os.path.join(outdir, basename)
        if not os.path.isdir(untardir) :
            os.makedirs(untardir)

        namelist = basename.split('_')
        tifname = os.path.join(outdir, basename+'.TIF')
        if os.path.isfile(tifname) :
            print('文件已存在【%s】' %(tifname))
            return tifname

        metafile =  basename + '_MTL.txt'
        metafile = self.untar(untardir, filename, metafile)
        meta_data = self.readMTL(metafile)

        allband = []
        for iband in np.arange(1, 8) :

            name = basename + '_SR_B%d.TIF' %(iband)
            outname = self.untar(untardir, filename, name)

            data, trans, prj = readtiff(outname)
            ref = self.GetReflectance(iband, data, meta_data)

            allband.append(ref)

        combdir = os.path.dirname(tifname)
        if not os.path.isdir(combdir) :
            os.makedirs(combdir)

        print('成功处理【%s】' %(tifname))
        print('='*100)

    def getMetaInfo(self, metafile):
        metadata = {}

        with open(metafile, 'r') as fp :
            lines = fp.readlines()

        for line in lines :
            if 'END_GROUP' in line :
                continue

            if 'GROUP' in line :
                line = line.replace('\n','')
                lineinfo = line.split('=')
                key1 = lineinfo[1].replace(' ', '')
                # key1 = lineinfo[1]
                metadata[key1] = {}
                continue

            if '=' in line :
                line = line.replace('\n','')
                lineinfo = line.split('=')
                key = lineinfo[0].replace(' ', '')
                # value = lineinfo[1].replace(' ', '')
                value = lineinfo[1]
                metadata[key1][key] = value

        return metadata

    def untar(self, outdir, filename, item):

        tar = tarfile.open(filename)

        outname = os.path.join(outdir, item)
        if not os.path.isfile(outname) :
            print('正在解压【%s】' %(outname))
            tar.extract(item, outdir)

        tar.close()

        return outname

    def getExtent(self, metadata):
        '''

        :param nowdate:
        :param BandId:
        :param radiance:
        :param metadata:
        :param dem:
        :return:
        '''

        # 中心经纬度
        point1lat = float(metadata['CORNER_UL_LAT_PRODUCT'])
        point1lon = float(metadata['CORNER_UL_LON_PRODUCT'])
        point2lat = float(metadata['CORNER_UR_LAT_PRODUCT'])
        point2lon = float(metadata['CORNER_UR_LON_PRODUCT'])
        point3lat = float(metadata['CORNER_LL_LAT_PRODUCT'])
        point3lon = float(metadata['CORNER_LL_LON_PRODUCT'])
        point4lat = float(metadata['CORNER_LR_LAT_PRODUCT'])
        point4lon = float(metadata['CORNER_LR_LON_PRODUCT'])

        lon = [point1lon, point2lon, point3lon, point4lon]
        lat = [point1lat, point2lat, point3lat, point4lat]

        sunz = 90-float(metadata['SUN_ELEVATION'])
        suna =    float(metadata['SUN_AZIMUTH'])

        return [np.nanmin(lon), np.nanmin(lat),
                np.nanmax(lon), np.nanmax(lat)]


    def GetRadiance(self, BandId, data, metadata, fillvalue=65535):
        ''' 计算辐射亮度参数：DN -> radiance '''

        if 'LEVEL2_SURFACE_REFLECTANCE_PARAMETERS' in metadata :
            Gain = float(metadata['LEVEL2_SURFACE_REFLECTANCE_PARAMETERS']['RADIANCE_MULT_BAND_%d' %(BandId)])
            Bias = float(metadata['LEVEL2_SURFACE_REFLECTANCE_PARAMETERS']['RADIANCE_ADD_BAND_%d' %(BandId)])
        else:
            Gain = float(metadata['LEVEL1_RADIOMETRIC_RESCALING']['RADIANCE_MULT_BAND_%d' %(BandId)])
            Bias = float(metadata['LEVEL1_RADIOMETRIC_RESCALING']['RADIANCE_ADD_BAND_%d' %(BandId)])

        radiance = np.where(data>0 ,Gain * data + Bias, fillvalue)
        radiance = np.array(radiance, dtype=np.float32)

        return radiance

    def GetReflectance(self, BandId, data, metadata, fillvalue=65535):
        ''' 计算反射率：DN -> Reflectance '''

        if 'LEVEL2_SURFACE_REFLECTANCE_PARAMETERS' in metadata :
            Gain = float(metadata['LEVEL2_SURFACE_REFLECTANCE_PARAMETERS']['REFLECTANCE_MULT_BAND_%d' %(BandId)])
            Bias = float(metadata['LEVEL2_SURFACE_REFLECTANCE_PARAMETERS']['REFLECTANCE_ADD_BAND_%d' %(BandId)])
        else:
            Gain = float(metadata['LEVEL1_RADIOMETRIC_RESCALING']['REFLECTANCE_MULT_BAND_%d' %(BandId)])
            Bias = float(metadata['LEVEL1_RADIOMETRIC_RESCALING']['REFLECTANCE_ADD_BAND_%d' %(BandId)])

        reflectance = np.where(data>0 ,Gain * data + Bias, fillvalue)
        reflectance = np.array(reflectance, dtype=np.float32)

        return reflectance