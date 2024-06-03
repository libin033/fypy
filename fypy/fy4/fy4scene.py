# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : fy4pro.py

@Modify Time :  2022/11/10 14:45

@Author      : fypy Team

@Version     : 1.0

@Description :

'''
import os
import sys
import numpy as np
import datetime
from osgeo import gdal, osr

from fypy.tools.tifpro import GetGDALType
from fypy.tools.ncpro import readnc, readnc_sdsinfo
from fypy.tools.BaseAlgorithms import BaseAlgorithms
from fypy.fy4.fy4core import fy4searchtable

from PIL import Image
RATE = 4


class FY4Scene(fy4searchtable, BaseAlgorithms) :

    def __init__(self, filename=None, satID=None, instID=None, prodID=None,
                 sublon=None, resolution=None, regionID=None, levelID=None,
                 startTime=None, endTime=None, proj=None):
        ''' 通过文件名获取文件相关信息 '''
        self.img = None

        self.Parse(filename, satID, instID, prodID, sublon, regionID,
                   levelID, startTime, endTime, resolution, proj)

        if (not hasattr(self, 'SubLon')) or (not hasattr(self, 'Resolution')) :
            raise Exception('请设置【sublon】和【resolution】参数')

        super().__init__(self.SubLon, self.Resolution)

    def calibrate(self, filename, bandID=1, fillvalue=65535, projFlag=False):
        '''
        读取FY4 L1数据，并完成辐射定标
        Parameters
        ----------
        filename: str
            L1数据文件名
        bandID : int
            波段索引，从1开始

        Returns
        -------
            numpy.array
            辐射定标转换后的ref或bt
            如果projFlag为True, 则对数据进行投影，返回gdal.Dateset
        '''
        if not os.path.isfile(filename) :
            print('文件不存在【%s】' %(filename))
            return None

        import h5py
        fp = h5py.File(filename, 'r')
        # 转换到区域的行列号（考虑去除图像偏移）
        Begin_Line_Number  = fp.attrs['Begin Line Number'][0]
        End_Line_Number    = fp.attrs['End Line Number'][0]
        Begin_Pixel_Number = fp.attrs['Begin Pixel Number'][0]
        End_Pixel_Number   = fp.attrs['End Pixel Number'][0]

        if 'FY4B' in filename :
            cal = fp['Calibration']['CALChannel%02d' %(bandID)][:]
            dn = fp['Data']['NOMChannel%02d' %(bandID)][:]
        else:
            cal = fp['CALChannel%02d' %(bandID)][:]
            dn = fp['NOMChannel%02d' %(bandID)][:]
        fp.close()

        flag = dn>=len(cal)
        dn[flag] = 0

        data = cal[dn]
        data[flag] = fillvalue

        if projFlag :
            return self.project(data, Begin_Line_Number, Begin_Pixel_Number, srcNodata=fillvalue)

        return data

    def getGEOData(self, filename, sdsname, fillvalue=65535, projFlag=False):
        import h5py
        fp = h5py.File(filename, 'r')
        data1 = fp[sdsname][:]

        # 转换到区域的行列号（考虑去除图像偏移）
        Begin_Line_Number = fp.attrs['Begin Line Number'][0]
        End_Line_Number = fp.attrs['End Line Number'][0]
        Begin_Pixel_Number = fp.attrs['Begin Pixel Number'][0]
        End_Pixel_Number = fp.attrs['End Pixel Number'][0]
        fp.close()

        if projFlag :
            return self.project(data1, Begin_Line_Number, Begin_Pixel_Number, srcNodata=fillvalue)
        else:
            return data1

    def project(self, data, sLine, sPixel, srcNodata=65535):

        im_data = np.array(data, dtype=np.float32)
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
                            '+lon_0={sublon} +no_defs'.format(sublon=self.SubLon))
        memDs.SetProjection(srs.ExportToWkt())

        # 设置仿射变换
        memDs.SetGeoTransform(self.SetTrans(sLine, sPixel))

        return memDs

    def load(self, filename, ProdID=None):

        data, dictsdsinfo  = readnc(filename, ProdID, dictsdsinfo={})
        if data is None :
            print('读取文件失败【%s】【%s】' %(ProdID, filename))
            return None

        if 'FillValue' in dictsdsinfo :
            fillvalue = dictsdsinfo['FillValue']
        elif '_FillValue' in dictsdsinfo :
            fillvalue = dictsdsinfo['_FillValue']
        else:
            fillvalue = None

        extentinfo = readnc_sdsinfo(filename, 'geospatial_lat_lon_extent')

        return self.project(data,
                            extentinfo['begin_line_number'],
                            extentinfo['begin_pixel_number'], srcNodata=fillvalue)

    def show(self, filename, ProdID='true_color'):
        vis045 = self.calibrate(filename, bandID=1, projFlag=False)
        vis065 = self.calibrate(filename, bandID=2, projFlag=False)
        vis085 = self.calibrate(filename, bandID=3, projFlag=False)

        rr = vis085.copy()
        rr[vis085<vis065] = vis065[vis085<vis065]

        vis055 = self.set_green(vis045, vis065)

        rr = np.array(rr*255, dtype=np.uint8)
        gg = np.array(vis055*255, dtype=np.uint8)
        bb = np.array(vis045*255, dtype=np.uint8)

        self.img = Image.merge('RGB', [self.Arr2Img(i) for i in (rr, gg, bb)])

        if self.img is None :
            raise Exception('请先加载【load】一个对象后再【show】')

        self.img.show()

    def saveThematic(self, outname):
        if self.img is None :
            raise Exception('请先加载【load】一个对象后再【SaveThematic】')

        self.img.save(outname)

    def SetTrans(self, sLine, sPixel):

        maxY = (self.rowmax/2.0-sLine) * self.Resolution * 100 * 1000
        minX = -(self.colmax/2.0-sPixel) * self.Resolution * 100 * 1000

        trans = [minX, self.Resolution*100*1000, 0,
                 maxY, 0, -self.Resolution*100*1000]

        return trans

    def Arr2Img(self, arr):
        if arr.dtype == np.uint8:
            return Image.fromarray(arr, "L")
        if arr.dtype == np.uint16:
            return Image.fromarray(self.T0_255(arr), "L")
        return Image.fromarray(self.T0_255(arr.astype('u2')), "L")

    def T0_255(self, raw):
        return (raw >> RATE).astype(np.uint8)

    def set_green(self, vis047, vis065, fractions=(1.0, 0.13, 0.87)):
        ''' 用红光和蓝光通道模拟绿光通道 '''

        res = (vis047 * fractions[0] - vis065 * fractions[1]) / fractions[2]

        return res

    def Parse(self, filename=None, SatID=None, InstID=None, ProdID=None, SubLon=None,
              RegionID=None, LevelID=None, StartTime=None, EndTime=None, Resolution=None,
              Proj=None):
        nameinfo = {}
        if filename is not None :
            nameinfo = self.GetNameInfo(filename)

        if SatID is not None :
            self.SatID = SatID
        elif 'SatID' in nameinfo :
            self.SatID = nameinfo['SatID']

        if InstID is not None :
            self.InstID = InstID
        elif 'InstID' in nameinfo :
            self.InstID = nameinfo['InstID']

        if ProdID is not None :
            self.ProdID = ProdID
        elif 'ProdID' in nameinfo :
            self.ProdID = nameinfo['ProdID']

        if Proj is not None :
            self.Proj = Proj
        elif 'Proj' in nameinfo :
            self.Proj = nameinfo['Proj']

        if RegionID is not None :
            self.RegionID = RegionID
        elif 'RegionID' in nameinfo :
            self.RegionID = nameinfo['RegionID']

        if LevelID is not None :
            self.LevelID = LevelID
        elif 'LevelID' in nameinfo :
            self.LevelID = nameinfo['LevelID']

        if StartTime is not None :
            self.StartTime = StartTime
        elif 'StartTime' in nameinfo :
            self.StartTime = nameinfo['StartTime']

        if EndTime is not None :
            self.EndTime = EndTime
        elif 'EndTime' in nameinfo :
            self.EndTime = nameinfo['EndTime']

        if Resolution is not None :
            self.Resolution = Resolution
        elif 'Resolution' in nameinfo :
            self.Resolution = nameinfo['Resolution']

        if SubLon is not None :
            self.SubLon = SubLon
        elif 'SubLon' in nameinfo :
            self.SubLon = float(nameinfo['SubLon'].replace('E', ''))/10.0


    def GetNameInfo(self, filename):
        ''' 根据输入数据的文件名获取信息 '''
        # FY4A-_AGRI--_N_DISK_1047E_L1-_FDI-_MULT_NOM_20211211060000_20211211061459_4000M_V0001.HDF
        # FY4A-_AGRI--_N_DISK_1047E_L1-_GEO-_MULT_NOM_20231231040000_20231231041459_4000M_V0001.HDF
        # FY4A-_AGRI--_N_DISK_1047E_L2-_CLM-_MULT_NOM_20240106030000_20240106031459_4000M_V0001.NC
        # FY4A-_AGRI--_N_REGC_1047E_L2-_CLM-_MULT_NOM_20231230011500_20231230011917_4000M_V0001.NC

        nameinfo = {}
        basename = os.path.basename(filename)
        if len(basename) != 88 and len(basename) != 89 :
            print('输入文件名为非标准文件名，将不作信息提取【%s】' %(basename))
            return nameinfo

        basename = basename.replace('-', '')
        namelist = basename.split('_')
        nameinfo['SatID']  = namelist[0]
        nameinfo['InstID'] = namelist[1]
        nameinfo['ObsType'] = namelist[2]
        nameinfo['RegionID'] = namelist[3]
        nameinfo['SubLon'] = namelist[4]
        nameinfo['LevelID'] = namelist[5]
        nameinfo['ProdID'] = namelist[6]
        nameinfo['Proj'] = namelist[8]
        nameinfo['StartTime'] = datetime.datetime.strptime(namelist[9], '%Y%m%d%H%M%S')
        nameinfo['EndTime'] = datetime.datetime.strptime(namelist[10], '%Y%m%d%H%M%S')
        if 'KM' in namelist[11] :
            nameinfo['Resolution'] = float(namelist[11].replace('KM',''))/100.0
        elif 'M' in namelist[11] :
            nameinfo['Resolution'] = float(namelist[11].replace('M',''))/100.0/1000.0
        nameinfo['Version'] = namelist[12]

        return nameinfo