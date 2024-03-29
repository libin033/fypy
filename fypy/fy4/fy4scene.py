# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : fy4pro.py

@Modify Time :  2022/11/10 14:45

@Author      : Lee

@Version     : 1.0

@Description :

'''
import os
import sys
import numpy as np
from osgeo import gdal, osr

from fypy.tools.tifpro import GetGDALType
from fypy.tools.ncpro import readnc, readnc_sdsinfo
from fypy.tools.BaseAlgorithms import BaseAlgorithms
from fypy.fy4.fy4core import GetNameInfo, fy4searchtable

from PIL import Image
RATE = 4


class fy4scene(fy4searchtable, BaseAlgorithms) :

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

    def Calibration(self, filename, bandID=1, fillvalue=65535):
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

        cal = fp['CALChannel%02d' %(bandID)][:]
        dn = fp['NOMChannel%02d' %(bandID)][:]
        fp.close()

        flag = dn>=len(cal)
        dn[flag] = 0

        data1 = cal[dn]
        data1[flag] = fillvalue

        return self.Reprojection(data1, Begin_Line_Number, Begin_Pixel_Number, srcNodata=fillvalue)

    def GetGEOData(self, filename, sdsname, fillvalue=65535):
        import h5py
        fp = h5py.File(filename, 'r')
        data1 = fp[sdsname][:]

        # 转换到区域的行列号（考虑去除图像偏移）
        Begin_Line_Number = fp.attrs['Begin Line Number'][0]
        End_Line_Number = fp.attrs['End Line Number'][0]
        Begin_Pixel_Number = fp.attrs['Begin Pixel Number'][0]
        End_Pixel_Number = fp.attrs['End Pixel Number'][0]
        fp.close()

        return self.Reprojection(data1, Begin_Line_Number, Begin_Pixel_Number, srcNodata=fillvalue)

    def Clip(self, srcDS, shpname=None, extent=None, srcNodata=65535, dstNodata=None,
             dstSRS='EPSG:4326', resampleAlg='near'):
        '''
        将标称投影转换成等经纬投影（NOM->GLL）

        Parameters
        ----------
        data : numpy.array
            输入数据
        shpname: str, optional
            掩膜的面矢量（polygon）
        extent : list, optional
            output bounds as (minX, minY, maxX, maxY) in target SRS
        dstNodata : float
            数据填充值
        dstSRS : str
            输出投影坐标系

        Returns
        -------

        '''

        # 影像裁剪
        warpDs = gdal.Warp('', srcDS, format='MEM', dstSRS=dstSRS,
                           cutlineDSName=shpname, cropToCutline=True,
                           outputBounds=extent, resampleAlg=resampleAlg,
                           xRes=self.Resolution, yRes=self.Resolution,
                           srcNodata=srcNodata, dstNodata=dstNodata)

        if warpDs is None :
            return None

        return warpDs

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

        return self.Reprojection(data,
                                 extentinfo['begin_line_number'],
                                 extentinfo['begin_pixel_number'], srcNodata=fillvalue)


    def show(self, filename, ProdID='true_color'):
        vis045 = self.Calibration(filename, bandID=1)
        vis065 = self.Calibration(filename, bandID=2)
        vis085 = self.Calibration(filename, bandID=3)

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

    def SaveThematic(self, outname):
        if self.img is None :
            raise Exception('请先加载【load】一个对象后再【SaveThematic】')

        self.img.save(outname)

    def Reprojection(self, data, sLine, sPixel, srcNodata=65535):

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

    def SetTrans(self, sLine, sPixel):

        maxY = (self.rowmax/2.0-sLine) * self.resolution * 100 * 1000
        minX = -(self.colmax/2.0-sPixel) * self.resolution * 100 * 1000

        trans = [minX, self.resolution*100*1000, 0,
                 maxY, 0, -self.resolution*100*1000]

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
            nameinfo = GetNameInfo(filename)

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
