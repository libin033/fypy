# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : modisscene.py

@Modify Time :  2024/2/1 16:01   

@Author      : fypy Team    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime
import tempfile
from osgeo import gdal, gdalconst, osr, ogr

from fypy.tools.BaseAlgorithms import BaseAlgorithms
from fypy.tools.tifpro import getDateSet
from fypy.tools.hdfpro import writehdf
from fypy.modis.modiscore import GetSourceInfo, CreateVrt

class MODISScene(BaseAlgorithms) :

    def __init__(self):
        self.TempFile = []

    def calibrate(self, filename, sdsname):

        srcDS = self._GetSourceInfo(filename, sdsname)
        meta = srcDS.GetMetadata()

        band_names = meta['band_names']
        bands = band_names.split(',')
        fillvalue = meta['_FillValue']
        valid_range = meta['valid_range']
        DN = srcDS.ReadAsArray()
        DstData = np.full_like(DN, fill_value=fillvalue, dtype=np.float32)
        for bandid in bands :
            index = bands.index(bandid)
            if '13' in bandid :
                ID = 13
            elif '14' in bandid :
                ID = 14
            else:
                ID = int(bandid)

            if ID >= 21 and ID <= 36 and ID != 26 :
                scales = meta['radiance_scales']
                offsets = meta['radiance_offsets']
            else:
                scales = meta['reflectance_scales']
                offsets = meta['reflectance_offsets']

            scales = np.array(scales.split(','), dtype=np.float64)
            offsets = np.array(offsets.split(','), dtype=np.float64)

            scale = scales[index]
            offset = offsets[index]

            if ID >= 20 and ID <= 36 and ID != 26 :
                data = scale * (DN[index] - offset)
                data = self.calibrate_bt(data, index)
                data[DN[index]==fillvalue] = fillvalue
            else:
                data = scale * (DN[index] - offset)
                data[DN[index]==fillvalue] = fillvalue

            DstData[index] = data

        return DstData

    def project(self, srcdata, srclat, srclon, resolution=None,
                vmin=None, vmax=None, extent=None, resampleAlg='near',
                srcNodata=-999.0, dstNodata=None, dstSRS="EPSG:4326"):
        '''  '''

        if resolution is None :
            if hasattr(self, 'Resolution') :
                resolution = self.__getattribute__('Resolution')
            else:
                raise Exception('需要指定输出空间分辨率【resolution】')

        if dstNodata is None :
            dstNodata = srcNodata

        flag = (srclon > 180) | (srclon < -180) | (srclat > 90)  | (srclat < -90)
        srclon[flag] = np.nan
        srclat[flag] = np.nan
        # srcdata[flag] = np.nan

        # 获取投影数据的范围
        if extent is None :
            extent = [np.nanmin(srclon), np.nanmin(srclat),
                      np.nanmax(srclon), np.nanmax(srclat)]

        # 获取输入数据的有效值范围
        if vmin is None : vmin = np.nanmin(srcdata)
        if vmax is None : vmax = np.nanmax(srcdata)
        if np.isnan(vmin) or np.isnan(vmax) :
            print('数据有效范围为NAN，将不做投影')
            return None

        data = np.array(srcdata).copy()

        if vmax is not None and vmin is not None :
            data[(srcdata < vmin) | (srcdata > vmax)] = srcNodata

        data[np.isnan(data)] = srcNodata

        tmp_file = tempfile.NamedTemporaryFile(prefix="tmp_fypy_orbit_", delete=True)
        temphdf = tmp_file.name + '.hdf'
        self.TempFile.append(temphdf)
        # 创建临时的数据文件
        writehdf(temphdf, 'data', data, overwrite=1)
        writehdf(temphdf, 'lon', srclon, overwrite=0)
        writehdf(temphdf, 'lat', srclat, overwrite=0)

        layer = GetSourceInfo(temphdf, 'data')

        vrtFile = tmp_file.name + '.vrt'
        self.TempFile.append(vrtFile)

        CreateVrt(vrtFile, temphdf, layer, '/lon', '/lat')

        ds = gdal.Warp('', vrtFile,
                       format='MEM',  geoloc=True,
                       dstSRS=dstSRS,  resampleAlg=resampleAlg,
                       srcNodata= srcNodata, dstNodata=dstNodata,
                       outputBounds=extent,  # (minX, minY, maxX, maxY)
                       xRes=resolution, yRes=resolution)

        if ds is None:
            print('处理失败')
            return None
        print('投影转换成功')

        return ds

    def converModisByGDAL(self, files, sdsname, resolution=None,
                          dstSRS='EPSG:4326',
                          resampleAlg=gdalconst.GRA_NearestNeighbour, maxerror=0.125):

        if isinstance(files, list):
            return self._MultiFiles(files, sdsname, resolution=resolution,
                                    dstSRS=dstSRS,
                                    resampleAlg=resampleAlg, maxerror=maxerror)
        elif isinstance(files, str):
            srcDS = self._GetSourceInfo(files, sdsname)
            return self.Reproject(srcDS, resolution, dstSRS=dstSRS,
                                  resampleAlg=resampleAlg,
                                  maxerror=maxerror)
        else:
            raise Exception('输出参数hdfname错误，仅支持str或者list')

    def quickProject(self, filename, sdsname, resolution=None,
                     dstSRS='EPSG:4326', resampleAlg=gdalconst.GRA_NearestNeighbour):
        ''' 基于MODIS的L1B文件进行快速投影，该投影精度存在一定偏差 '''

        target_dataset = GetSourceInfo(filename, sdsname)
        ds = gdal.Warp('', target_dataset,
                       format='MEM',
                        dstSRS=dstSRS, resampleAlg=resampleAlg,
                       xRes=resolution, yRes=resolution,
                       multithread=True, geoloc=True)
        return ds

    def _GetSourceInfo(self, filename, sdsname):

        # 获取sdsname所在的图层栅格索引
        srcDS = gdal.Open(filename)
        layers = srcDS.GetSubDatasets()

        # 获取sdsname所在的图层栅格索引
        if sdsname:
            for layer in layers :
                l_name = layer[0].split(':')[-1].replace('"','')

                if sdsname in l_name:
                    dstDS = gdal.Open(layer[0])
                    del srcDS

                    return dstDS

        print('文件中不含数据集【%s】<--【%s】' %(sdsname, filename))
        return None


        # src_proj = srcDS.GetProjection()
        # src_trans = srcDS.GetGeoTransform()
        #
        # src_meta = srcDS.GetMetadata()
        # #  列数
        # # self.tileColumns = int(self.src_meta["DATACOLUMNS"])
        # #  行数
        # # self.tileRows = int(self.src_meta["DATAROWS"])
        #
        # tiledata = srcDS.ReadAsArray()
        # # if 'scale_factor' in self.src_meta and 'add_offset' in self.src_meta :
        # #     print('scale_factor: ', np.float(self.src_meta['scale_factor']) , ' add_offset: ', np.float(self.src_meta['add_offset']))
        # #     self.tiledata = self.tiledata * np.float32(self.src_meta['scale_factor']) + np.float32(self.src_meta['add_offset'])
        #
        # datasize = tiledata.shape
        # if len(datasize) == 2 :
        #     #  列数
        #     tileColumns = int(datasize[1])
        #     #  行数
        #     tileRows = int(datasize[0])
        # elif len(datasize) == 3 :
        #     #  列数
        #     tileColumns = int(datasize[2])
        #     #  行数
        #     tileRows = int(datasize[1])
        #
        # # self.src_meta = src_driver.GetMetadata()
        # band = srcDS.GetRasterBand(1)
        #
        # if '_FillValue' in list(src_meta.keys()):
        #     data_fill_value = src_meta['_FillValue']
        # elif band.GetNoDataValue():
        #     data_fill_value = band.GetNoDataValue()
        # else:
        #     data_fill_value = None
        # datatype = band.DataType

    def _MultiFiles(self, files, sdsname, resolution,
                    dstSRS, resampleAlg, maxerror=0.125):

        countfile = len(files)
        if countfile == 0 :
            return None
        elif countfile == 1 :
            srcDS = self._GetSourceInfo(files, sdsname)
            return self.Reproject(srcDS, resolution, dstSRS=dstSRS,
                                  resampleAlg=resampleAlg,
                                  maxerror=maxerror)
        else:
            mosaicds = []
            ulx = None
            uly = None
            lrx = None
            lry = None
            psize_x = None
            psize_y = None
            for filename in files :
                srcds = self._GetSourceInfo(filename, sdsname)
                bandtype = srcds.GetRasterBand(1).DataType
                bandcnt = srcds.RasterCount
                tempds = self.Reproject(srcds, resolution, dstSRS=dstSRS,
                               resampleAlg=resampleAlg,
                               maxerror=maxerror)
                if tempds is not None:
                    mosaicds.append(tempds)

                width  = tempds.RasterXSize
                height = tempds.RasterYSize
                trans  = tempds.GetGeoTransform()
                fillvalue = tempds.GetRasterBand(1).GetNoDataValue()
                if ulx is None or uly is None or lrx is None or lry is None :
                    ulx = trans[0]
                    uly = trans[3]
                    lrx = ulx + trans[1] * width
                    lry = uly + trans[5] * height
                    psize_x = trans[1]
                    psize_y = trans[5]
                else:
                    ulx = min(ulx, trans[0])
                    uly = max(uly, trans[3])
                    lrx = max(lrx, trans[0] + trans[1] * width)
                    lry = min(lry, trans[3] + trans[5] * height)

                if resolution is None :
                    resolution = psize_x

            if len(mosaicds) == 0 :
                return None

            geotransform = [ulx, psize_x, 0, uly, 0, psize_y]
            xsize = int(np.ceil((lrx - ulx) / geotransform[1]))
            ysize = int(np.ceil((lry - uly) / geotransform[5]))
            dr = gdal.GetDriverByName("MEM")
            dstDS = dr.Create('', xsize, ysize, bandcnt, bandtype)
            if fillvalue :
                for i in range(bandcnt) :
                    dstDS.GetRasterBand(i+1).SetNoDataValue(fillvalue)
                    dstDS.GetRasterBand(i+1).Fill(fillvalue)
            dstDS.SetProjection(dstSRS)
            dstDS.SetGeoTransform(geotransform)

            for ds in mosaicds :
                trans = ds.GetGeoTransform()
                width  = ds.RasterXSize
                height = ds.RasterYSize

                xoffset = int(np.abs((trans[0] - geotransform[0])/geotransform[1]))
                yoffset = int(np.fabs((trans[3] - geotransform[3])/geotransform[5]))

                srcdata = ds.ReadAsArray()
                if fillvalue :
                    dstdata = dstDS.ReadAsArray(xoff=xoffset, yoff=yoffset, xsize=width, ysize=height)
                    srcdata[srcdata==fillvalue] = dstdata[srcdata==fillvalue]
                dstDS.WriteArray(srcdata, xoff=xoffset, yoff=yoffset)

            return dstDS

    def Reproject(self, srcDS, resolution=None, dstSRS=None,
                  resampleAlg=gdalconst.GRA_NearestNeighbour,
                  maxerror = 0.125):

        # 投影转换
        tmp_ds = gdal.AutoCreateWarpedVRT(srcDS,
                                          srcDS.GetProjection(),
                                          dstSRS,
                                          resampleAlg,
                                          maxerror)
        if tmp_ds is None :
            return None

        if not resolution:
            width = tmp_ds.RasterXSize
            height = tmp_ds.RasterYSize
            dstTrans = tmp_ds.GetGeoTransform()
        else:
            bbox = self._boundingBox(tmp_ds)
            width = self._calculateRes(bbox[0][0], bbox[1][0],
                                       resolution)
            height = self._calculateRes(bbox[0][1], bbox[1][1],
                                        resolution)

            dstTrans = [bbox[0][0], resolution, 0.0,
                        bbox[1][1], 0.0, -resolution]
        del tmp_ds

        #创建文件
        driver = gdal.GetDriverByName("MEM")
        dstDS = driver.Create('', width, height,
                              srcDS.RasterCount,
                              srcDS.GetRasterBand(1).DataType)
        if dstDS is  None:
            return None

        dstDS.SetProjection(dstSRS)
        dstDS.SetGeoTransform(dstTrans)

        meta = srcDS.GetMetadata()
        if '_FillValue' in list(meta.keys()):
            fillvalue = meta['_FillValue']
        elif srcDS.GetRasterBand(1).GetNoDataValue():
            fillvalue = srcDS.GetRasterBand(1).GetNoDataValue()
        else:
            fillvalue = None

        if fillvalue :
            for i in range(srcDS.RasterCount) :
                dstDS.GetRasterBand(i+1).SetNoDataValue(float(fillvalue))
                dstDS.GetRasterBand(i+1).Fill(float(fillvalue))

        cbk = self._progressCallback
        cbk_user_data = None
        try:
            gdal.ReprojectImage(srcDS, dstDS, srcDS.GetProjection(),
                                dstDS.GetProjection(), resampleAlg,
                                0, maxerror, cbk, cbk_user_data)
        except:
            raise Exception('Not possible to reproject dataset')
        dstDS.SetMetadata(meta)

        return dstDS

    def _boundingBox(self, src):

        src_gtrn = src.GetGeoTransform(can_return_null=True)

        src_bbox_cells = ((0., 0.),
                          (0, src.RasterYSize),
                          (src.RasterXSize, 0),
                          (src.RasterXSize, src.RasterYSize))

        geo_pts_x = []
        geo_pts_y = []
        for x, y in src_bbox_cells:
            x2 = src_gtrn[0] + src_gtrn[1] * x + src_gtrn[2] * y
            y2 = src_gtrn[3] + src_gtrn[4] * x + src_gtrn[5] * y
            geo_pts_x.append(x2)
            geo_pts_y.append(y2)
        return ((min(geo_pts_x), min(geo_pts_y)), (max(geo_pts_x),
                                                   max(geo_pts_y)))

    def _calculateRes(self, minn, maxx, res):
        """Calculate the number of pixel from extent and resolution

           :param float minn: minimum value of extent
           :param float maxx: maximum value of extent
           :param int res: resolution of output raster

           :return: integer number with the number of pixels
        """
        return int(round((maxx - minn) / res))

    def _progressCallback(self, pct, message, user_data):
        """For the progress status"""
        return 1  # 1 to continue, 0 to stop

    def calibrate_bt(self, array,  index):
        """Calibration for the emissive channels."""
        # Planck constant (Joule second)
        h__ = np.float32(6.6260755e-34)

        # Speed of light in vacuum (meters per second)
        c__ = np.float32(2.9979246e+8)

        # Boltzmann constant (Joules per Kelvin)
        k__ = np.float32(1.380658e-23)

        # Derived constants
        c_1 = 2 * h__ * c__ * c__
        c_2 = (h__ * c__) / k__

        # Effective central wavenumber (inverse centimeters)
        cwn = np.array([
            2.641775E+3, 2.505277E+3, 2.518028E+3, 2.465428E+3,
            2.235815E+3, 2.200346E+3, 1.477967E+3, 1.362737E+3,
            1.173190E+3, 1.027715E+3, 9.080884E+2, 8.315399E+2,
            7.483394E+2, 7.308963E+2, 7.188681E+2, 7.045367E+2],
            dtype=np.float32)

        # Temperature correction slope (no units)
        tcs = np.array([
            9.993411E-1, 9.998646E-1, 9.998584E-1, 9.998682E-1,
            9.998819E-1, 9.998845E-1, 9.994877E-1, 9.994918E-1,
            9.995495E-1, 9.997398E-1, 9.995608E-1, 9.997256E-1,
            9.999160E-1, 9.999167E-1, 9.999191E-1, 9.999281E-1],
            dtype=np.float32)

        # Temperature correction intercept (Kelvin)
        tci = np.array([
            4.770532E-1, 9.262664E-2, 9.757996E-2, 8.929242E-2,
            7.310901E-2, 7.060415E-2, 2.204921E-1, 2.046087E-1,
            1.599191E-1, 8.253401E-2, 1.302699E-1, 7.181833E-2,
            1.972608E-2, 1.913568E-2, 1.817817E-2, 1.583042E-2],
            dtype=np.float32)

        # Transfer wavenumber [cm^(-1)] to wavelength [m]
        cwn = 1. / (cwn * 100)

        # Some versions of the modis files do not contain all the bands.
        # emmissive_channels = ["20", "21", "22", "23", "24", "25", "27", "28", "29",
        #                       "30", "31", "32", "33", "34", "35", "36"]
        # global_index = emmissive_channels.index(band_name)

        cwn = cwn[index]
        tcs = tcs[index]
        tci = tci[index]
        array = c_2 / (cwn * np.log(c_1 / (1000000 * array * cwn ** 5) + 1))
        array = (array - tci) / tcs
        return array


    def __del__(self):

        for tempfile in  self.TempFile :
            if os.path.isfile(tempfile) :
                try:
                    os.remove(tempfile)
                except BaseException :
                    print('删除临时文件失败【%s】' %(tempfile))
