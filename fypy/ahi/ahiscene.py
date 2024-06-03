# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : ahiscene.py

@Modify Time :  2024/1/29 18:26   

@Author      : fypy Team    

@Version     : 1.0   

@Description :

'''
import os
import sys
import time

import numpy as np
import datetime
from PIL import Image
from tqdm import tqdm

from osgeo import gdal, osr
from fypy.tools.tifpro import writetiff, GetGDALType
from fypy.tools.BaseAlgorithms import BaseAlgorithms
from fypy.ahi.ahisearchtable import ahisearchtable
from fypy.ahi.ahi_read_hsd import ahi_read_hsd
from fypy.ahi.ahiconfig import AreaInfo


class AHIScene(ahi_read_hsd, ahisearchtable, BaseAlgorithms) :

    def __init__(self, subpoint, resolution):
        super().__init__(subpoint, resolution)
        self.Tempfile = []

    def hsdBlock(self, srcHSDfiles, tmppath, fillvalue=65535) :
        ''' 对H8、H9的HSD文件进行解析、拼接成NOM '''

        # HS_H09_20230115_0400_B01_FLDK_R10_S0110.DAT.bz2
        BandID, BlockIDMin, BlockIDMax, SegmentTotal = self.setHSDInfo(srcHSDfiles)
        outdata = None
        BlockIDs = []
        with tqdm(total=len(srcHSDfiles), iterable='iterable',
                  desc = '正在进行第%i波段块合成' %(BandID), mininterval=1) as pbar:
            for hsdname in srcHSDfiles :
                if not os.path.isfile(hsdname):
                    print('文件不存在【%s】' %(hsdname))
                    pbar.update(1)
                    continue

                # 获取文件名信息
                nameinfo = self.getHSDNameInfo(hsdname)
                if nameinfo is None :
                    pbar.update(1)
                    continue
                SegmentNum = nameinfo['SegmemtID']

                # print('正在解压bz2文件【%s】' %(hsdname))
                self._unzipped = self.unzip_file(hsdname, tmppath)
                if self._unzipped:
                    self.is_zipped = True
                    filename = self._unzipped

                    self.Tempfile.append(filename)
                else:
                    filename = hsdname

                if filename.endswith('.bz2') :
                    print('解压bz2文件失败【%s】' %(filename))
                    pbar.update(1)
                    continue

                # 根据块号对数据进行拼接
                data = self.readhsd(filename, SegmentNum)
                if data is None :
                    pbar.update(1)
                    continue

                if outdata is None :
                    line, pixel = data.shape
                    outdata = np.full(shape=(line*SegmentTotal, pixel),
                                      fill_value=fillvalue, dtype=np.uint16)

                data[np.isnan(data)] = fillvalue/100.0
                outdata[(SegmentNum-BlockIDMin)*line:(SegmentNum-BlockIDMin+1)*line, :] \
                    = np.array(data*100.0, dtype=np.uint16)
                BlockIDs.append(SegmentNum)
                pbar.update(1)
        pbar.close()

        return self.project(outdata, BlockIDs, line, srcNodata=fillvalue)

    def project(self, data, blockIDs, blockLine, srcNodata=65535):

        im_data = np.array(data)
        dtype = GetGDALType(im_data.dtype)

        # 只支持2、3维数据处理
        if len(im_data.shape) == 2:
            im_bands, (im_height, im_width) = 1,im_data.shape
        elif len(im_data.shape) == 3:
            im_bands, im_height, im_width = im_data.shape

        Driver = gdal.GetDriverByName('MEM')
        memDs = Driver.Create('', im_width, im_height, im_bands, dtype)

        # 写入数据
        if im_bands == 1:
            memDs.GetRasterBand(1).WriteArray(im_data)
            if srcNodata is not None:
                memDs.GetRasterBand(1).SetNoDataValue(float(srcNodata))
        else:
            for i in range(im_bands):
                memDs.GetRasterBand(i+1).WriteArray(im_data[i])
                if srcNodata is not None:
                    memDs.GetRasterBand(i+1).SetNoDataValue(float(srcNodata))

        # 设置参考投影
        srs = osr.SpatialReference()
        srs.ImportFromProj4('+proj=geos +h=35785863 +a=6378137.0 +b=6356752.3 '
                            '+lon_0={sublon} +no_defs'.format(sublon=self.subpoint))
        memDs.SetProjection(srs.ExportToWkt())

        # 设置仿射变换
        memDs.SetGeoTransform(self.setTrans(BlockIDs=blockIDs, blockLine=blockLine))

        return memDs

    def setTrans(self, BlockIDs, blockLine) :
        ''' 设置仿射变换 '''

        if str(self.resolution) not in AreaInfo['ahi']['ahi'] :
            raise Exception('该空间分辨率【%s】不在处理范围内【0.005，0.01，0.02】' )

        dictinfo = AreaInfo['ahi']['ahi'][str(self.resolution)]

        BlockIDMin = np.nanmin(BlockIDs)
        BlockIDMax = np.nanmax(BlockIDs)

        maxY = dictinfo['extent'][3] - ((BlockIDMin-1)*blockLine) * self.resolution * 100 * 1000
        minX = dictinfo['extent'][0]

        trans = [minX, self.resolution*100*1000, 0,
                 maxY, 0, -self.resolution*100*1000]

        return trans

    def load(self, vis650, vis550, vis450):
        '''
        绘制葵花数据真彩图

        Parameters
        ----------
        outname : string
            输出文件名，PNG或者JPG
        vis650 : numpy.array-2D
            650um通道的反射率，范围在0-1
        vis550 : numpy.array-2D
            550um通道的反射率，范围在0-1
        vis450 : numpy.array-2D
            450um通道的反射率，范围在0-1

        Returns
        -------

        '''

        r = vis650 * 255
        g = vis550 * 255
        b = vis450 * 255

        r[r<0] = 0
        g[g<0] = 0
        b[b<0] = 0

        r[r>255] = 255
        g[g>255] = 255
        b[b>255] = 255

        rgbArray = np.zeros((r.shape[0], r.shape[1], 4), dtype=np.float64)
        rgbArray[..., 0] = r
        rgbArray[..., 1] = g
        rgbArray[..., 2] = b
        rgbArray[..., 3] = 255

        img = Image.fromarray(rgbArray.astype(np.uint8))

    def setHSDInfo(self, filelist):

        BandID = None
        BlockIDs = []
        for filename in filelist :
            nameinfo = self.getHSDNameInfo(filename)
            if nameinfo is None :
                continue

            if BandID is None :
                BandID = nameinfo['BandID']
            elif BandID != nameinfo['BandID'] :
                raise Exception('输入的文件列表中有多个波段的块数据文件【%s】' %(filename))
            BlockIDs.append(nameinfo['SegmemtID'])

        BlockIDMin = np.nanmin(BlockIDs)
        BlockIDMax = np.nanmax(BlockIDs)

        SegmentTotal = int(BlockIDMax-BlockIDMin+1)

        return BandID, BlockIDMin, BlockIDMax, SegmentTotal

    def getHSDNameInfo(self, filename):

        basename = os.path.basename(filename)
        basename = basename.split('.')[0]
        if len(basename) != 39 :
            print('非标准文件名，需要输入文件名【HS_H09_YYYYMMDD_HHMM_BXX_FLDK_R20_S0810】')
            return None

        nameinfo = {}
        namelist = basename.split('_')

        nameinfo['SatID'] = namelist[1]
        nameinfo['StartTime'] = datetime.datetime.strptime('%s %s' %(namelist[2], namelist[3]), '%Y%m%d %H%M')
        nameinfo['BandID'] = int(namelist[4][1:])   # 2-digit band number (varies from "01" to "16");
        nameinfo['ObsType'] = namelist[5]
        nameinfo['Resolution'] = float(namelist[6][1:])/10.0/100    # spatial resolution ("05": 0.5km, "10": 1.0km, "20": 2.0km);
        nameinfo['SegmemtID'] = int(namelist[7][1:3])
        nameinfo['SegmemtTotal'] = int(namelist[7][3:5])    # total number of segments (fixed to "10")

        return nameinfo

    def __del__(self):
        pass
        # for filename in self.Tempfile :
        #     if os.path.isfile(filename) :
        #         try:
        #             os.remove(filename)
        #         except BaseException as e :
        #             time.sleep(1)
        #             try:
        #                 fp = open(filename, 'r')
        #                 fp.close()
        #                 os.remove(filename)
        #             except BaseException as e :
        #                 pass