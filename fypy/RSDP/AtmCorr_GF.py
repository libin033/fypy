# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : AtmCorr_GF.py

@Modify Time :  2022/12/23 16:38   

@Author      : Lee    

@Version     : 1.0   

@Description :

'''
import os
import shutil
import sys
import numpy as np
import datetime
import json
from tqdm import tqdm
from .AtmCorr import AtmCorr
from osgeo import gdal, gdalconst
from dateutil.relativedelta import relativedelta
import xml.dom.minidom
EXEPATH = os.path.split(os.path.realpath(__file__))[0]

class AtmCorr_GF(AtmCorr):

    def __init__(self):
        super(AtmCorr_GF, self).__init__()

    def FLAASH(self, outname, tiffFile, metedata, nowdate, SatID, InstID,
               InstType=None, fillvalue=65535, blocksize=1000):

        outdir = os.path.dirname(outname)
        if not os.path.isdir(outdir) :
            os.makedirs(outdir)
            print('成功创建路径【%s】' %(outdir))

        if not os.path.isfile(tiffFile) :
            print('文件不存在【%s】' %(tiffFile))
            return None

        srcdataset = gdal.Open(tiffFile, gdal.GA_ReadOnly)
        if srcdataset == None:
            print("文件无法打开【%s】" %(tiffFile))
            return None

        cols = srcdataset.RasterXSize                  # 栅格矩阵的列数
        rows = srcdataset.RasterYSize                  # 栅格矩阵的行数
        Bands = srcdataset.RasterCount                 # 波段数

        geoTransform1 = srcdataset.GetGeoTransform()   # 获取仿射矩阵信息
        proj1 = srcdataset.GetProjection()             # 获取投影信息

        # 创建输出结果文件
        driver = gdal.GetDriverByName("GTiff")
        dstdataset = driver.Create(outname, cols, rows, Bands, gdal.GDT_UInt16,
                                   options=["COMPRESS=LZW", "BigTIFF=YES"])
        dstdataset.SetGeoTransform(geoTransform1)
        dstdataset.SetProjection(proj1)

        #分别读取4个波段
        for BandID in range(0, Bands):
            ReadBand = srcdataset.GetRasterBand(BandID+1)
            outband = dstdataset.GetRasterBand(BandID+1)
            outband.SetNoDataValue(fillvalue)
            #获取对应波段的增益gain和偏移bias
            Gain, Bias = self.RadiometricCalibration(nowdate, SatID, InstID)

            if 'PMS' in InstID :
                if InstType is None :
                    raise Exception('请传入InstType参数【全色：PAN，多光谱：MSS】')

                if 'MSS' in InstType :
                    BandID += 1
                elif 'PAN' in InstType :
                    BandID = BandID
                else:
                    raise Exception('请传入InstType参数【全色：PAN，多光谱：MSS】')

            #获取大气校正系数
            # if 'PMS' in InstID :
            #     if 'PAN' not in InstType :
            #         self.setParam(nowdate, metedata, SatID, InstID, BandID+1)
            #         AtcCofa, AtcCofb, AtcCofc = self.corrCoeff()

            self.setParam(nowdate, metedata, SatID, InstID, BandID+1)
            AtcCofa, AtcCofb, AtcCofc = self.corrCoeff()


            i = 0
            j = 0
            #进度条参数
            XBlockcount = np.ceil(cols / blocksize)
            YBlockcount = np.ceil(rows / blocksize)
            # print("第%d波段校正："%BandID)
            try:
                with tqdm(total=XBlockcount*YBlockcount, iterable='iterable',
                          desc = '正在进行第%i波段校正' %(BandID+1), mininterval=1) as pbar:
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
                            Image = ReadBand.ReadAsArray(j, i, nXBK,nYBK)

                            # if 'PMS' in InstID :
                            #     if 'PAN' in InstType :  # 全色波段将不做大气校正，直接返回
                            #         outImage =np.where(Image>0, (Image*Gain[BandID] + Bias[BandID])*10000, fillvalue)
                            #         outband.WriteArray(outImage,j,i)
                            #         j=j+nXBK
                            #         pbar.update(1)
                            #         continue

                            outImage =np.where(Image>0, Image*Gain[BandID] + Bias[BandID], fillvalue)

                            y = np.where(outImage!=fillvalue, AtcCofa * outImage - AtcCofb, fillvalue)
                            atcImage = np.where(y!=fillvalue, (y / (1 + y * AtcCofc))*10000, fillvalue)

                            outband.WriteArray(atcImage,j,i)
                            j=j+nXBK
                            pbar.update(1)
                        j=0
                        i=i+nYBK
            except BaseException :
                pbar.close()
                continue
            pbar.close()

        dstdataset = None
        srcdataset = None

    def setParam(self, nowdate, metedata, SatID, InstID, BandId, dem=0.010):

        # 读取头文件
        dom = xml.dom.minidom.parse(metedata)

        # 太阳和卫星角度信息
        sunz = 90 - float(dom.getElementsByTagName('SolarZenith')[0].firstChild.data)
        suna = float(dom.getElementsByTagName('SolarAzimuth')[0].firstChild.data)
        # s.geometry.view_z = float(dom.getElementsByTagName('SatelliteZenith')[0].firstChild.data)
        # s.geometry.view_a = float(dom.getElementsByTagName('SatelliteAzimuth')[0].firstChild.data)
        satz = 0
        sata = 0
        # 日期
        DateTimeparm = dom.getElementsByTagName('CenterTime')[0].firstChild.data
        DateTime = datetime.datetime.strptime(DateTimeparm, '%Y-%m-%d %H:%M:%S')

        # 中心经纬度
        TopLeftLat     = float(dom.getElementsByTagName('TopLeftLatitude')[0].firstChild.data)
        TopLeftLon     = float(dom.getElementsByTagName('TopLeftLongitude')[0].firstChild.data)
        TopRightLat    = float(dom.getElementsByTagName('TopRightLatitude')[0].firstChild.data)
        TopRightLon    = float(dom.getElementsByTagName('TopRightLongitude')[0].firstChild.data)
        BottomRightLat = float(dom.getElementsByTagName('BottomRightLatitude')[0].firstChild.data)
        BottomRightLon = float(dom.getElementsByTagName('BottomRightLongitude')[0].firstChild.data)
        BottomLeftLat  = float(dom.getElementsByTagName('BottomLeftLatitude')[0].firstChild.data)
        BottomLeftLon  = float(dom.getElementsByTagName('BottomLeftLongitude')[0].firstChild.data)

        ImageCenterLat = (TopLeftLat + TopRightLat + BottomRightLat + BottomLeftLat) / 4

        self.set_geom(nowdate, sunz=sunz, suna=suna, satz=satz, sata=sata)
        self.set_atm(sLatitude=ImageCenterLat, nowdate=nowdate)
        self.set_aer()
        self.set_vis()
        self.set_altitude(dem=dem)

        #读取辐射校正和大气校正所需参数:增益、偏移和光谱响应函数
        SRFFile = os.path.join(EXEPATH, 'resp', 'GF', "GF.json")
        if not os.path.isfile(SRFFile) :
            raise Exception('高分卫星光谱响应文件不存在【%s】' %(SRFFile))
        SRF = json.load(open(SRFFile))

        minwl = SRF[SatID][InstID]['B%d' %(BandId)]['wl'][0]
        maxwl = SRF[SatID][InstID]['B%d' %(BandId)]['wl'][1]
        response = SRF[SatID][InstID]['B%d' %(BandId)]['SRF']
        self.set_resp(minwl, maxwl, response)

    def RadiometricCalibration(self, nowdate, SatID, InstID):
        ''' 获取高分卫星定标系数 gain和 bias '''

        #读取辐射校正和大气校正所需参数:增益、偏移和光谱响应函数
        CalCoeffFile = os.path.join(EXEPATH, "resp", "GF", "CalibrationCoefficient.json")
        if not os.path.isfile(CalCoeffFile) :
            raise Exception('高分卫星定标系数文件不存在【%s】' %(CalCoeffFile))
        CalCoeff = json.load(open(CalCoeffFile))

        satid = SatID.replace('-', '')
        instid = InstID.replace('-', '')

        if not satid in CalCoeff :
            raise Exception('请确认【%s】是否在定标系数列表，当前仅支持' %(SatID), CalCoeff.keys())

        if not instid in CalCoeff[satid] :
            raise Exception('请确认【%s】是否在定标系数列表，当前仅支持' %(InstID), CalCoeff[satid].keys())

        while nowdate >= datetime.datetime.strptime('2000', '%Y') :
            stryear = nowdate.strftime('%Y')
            if not stryear in CalCoeff[satid][instid] :
                nowdate -= relativedelta(years=1)
            else:
                gain = CalCoeff[satid][instid][stryear]['gain']
                bias = CalCoeff[satid][instid][stryear]['bias']

                return gain, bias

        # raise Exception('未匹配到该年【%s】系数' %(nowdate.strftime('%Y')))
        return None, None

    def ortho_rectification(self, dstfile, srcfile, rpcfile, demfile=None):
        """
        正射校正
        :param input:输入原始影像
        :param output:输出正射影像
        """

        if not os.path.isfile(rpcfile) :
            raise Exception('RPC文件不存在【%s】' %(rpcfile))

        srcdir = os.path.dirname(srcfile)
        shutil.copy(rpcfile, srcdir)

        dataset = gdal.Open(srcfile, gdal.GA_Update)#读入影像
        rpc = dataset.GetMetadata("RPC")#读入影像，rpc

        if demfile is None :
            dst_ds = gdal.Warp(dstfile, dataset, dstSRS='EPSG:4326',
                               # xRes=resolution,
                               # yRes=resolution,
                               resampleAlg=gdal.GRIORA_Bilinear,
                               rpc=True, #使用RPC模型进行校正
                               warpOptions=['INIT_DEST=NO_DATA'],
                               creationOptions=["COMPRESS=LZW", "BigTIFF=YES"])
        else:
            dst_ds = gdal.Warp(dstfile, dataset, dstSRS='EPSG:4326',
                               # xRes=resolution,
                               # yRes=resolution,
                               # resampleAlg=gdal.GRIORA_Bilinear,
                               rpc=True, #使用RPC模型进行校正
                               transformerOptions=[r'RPC_DEM=%s' %(demfile),
                                                   "RPC_DEMINTERPOLATION=bilinear"],
                               creationOptions=["COMPRESS=LZW", "BigTIFF=YES"])

        dst_ds = None

    def pansharpen(self,
                   dst_filename,
                   pan_name,
                   dict_spectrac=None,
                   format = None,
                   resampleAlg= None,
                   nodata_value = None,
                   weights = None,
                   creation_options = None,
                   spat_adjust = None,
                   num_threads = None,
                   bitdepth = None,
                   verbose_vrt = False,
                   progress_callback = gdal.TermProgress_nocb):

        band_nums = []
        weights = weights or []
        creation_options = creation_options or []
        spectral_ds = None
        spectral_bands = None

        if dict_spectrac :
            spectral_ds, spectral_bands = self.parse_spectral_names(dict_spectrac)

        if pan_name is None or spectral_bands is None:
            return 1

        # 打开PAN图像
        pan_ds = gdal.Open(pan_name)
        if pan_ds is None:
            return 1

        if format is None:
            format = self.GetOutputDriverFor(dst_filename)

        if not band_nums:
            band_nums = [j + 1 for j in range(len(spectral_bands))]
        else:
            for band in band_nums:
                if band < 0 or band > len(spectral_bands):
                    print('Invalid band number in -b: %d' % band)
                    return 1

        if weights and len(weights) != len(spectral_bands):
            print('There must be as many -w values specified as input spectral bands')
            return 1

        # 开始构建VRT文件
        vrt_xml = """<VRTDataset subClass="VRTPansharpenedDataset">\n"""
        if band_nums != [j + 1 for j in range(len(spectral_bands))]:
            for i, band in enumerate(band_nums):
                sband = spectral_bands[band - 1]
                datatype = gdal.GetDataTypeName(sband.DataType)
                colorname = gdal.GetColorInterpretationName(sband.GetColorInterpretation())
                vrt_xml += """  <VRTRasterBand dataType="%s" band="%d" subClass="VRTPansharpenedRasterBand">
          <ColorInterp>%s</ColorInterp>
      </VRTRasterBand>\n""" % (datatype, i + 1, colorname)

        vrt_xml += """  <PansharpeningOptions>\n"""

        # 对各个波段增加权重
        if weights:
            vrt_xml += """      <AlgorithmOptions>\n"""
            vrt_xml += """        <Weights>"""
            for i, weight in enumerate(weights):
                if i > 0:
                    vrt_xml += ","
                vrt_xml += "%.16g" % weight
            vrt_xml += "</Weights>\n"
            vrt_xml += """      </AlgorithmOptions>\n"""

        if resampleAlg is not None:
            vrt_xml += f'      <Resampling>{resampleAlg}</Resampling>\n'

        if num_threads is not None:
            vrt_xml += f'      <NumThreads>{num_threads}</NumThreads>\n'

        if bitdepth is not None:
            vrt_xml += f'      <BitDepth>{bitdepth}</BitDepth>\n'

        if nodata_value is not None:
            vrt_xml += f'      <NoData>{nodata_value}</NoData>\n'

        if spat_adjust is not None:
            vrt_xml += f'      <SpatialExtentAdjustment>{spat_adjust}</SpatialExtentAdjustment>\n'

        pan_relative = '0'
        if format.upper() == 'VRT':
            if not os.path.isabs(pan_name):
                pan_relative = '1'
                pan_name = os.path.relpath(pan_name, os.path.dirname(dst_filename))

        vrt_xml += """    <PanchroBand>
          <SourceFilename relativeToVRT="%s">%s</SourceFilename>
          <SourceBand>1</SourceBand>
        </PanchroBand>\n""" % (pan_relative, pan_name)

        for i, sband in enumerate(spectral_bands):
            dstband = ''
            for j, band in enumerate(band_nums):
                if i + 1 == band:
                    dstband = ' dstBand="%d"' % (j + 1)
                    break

            ms_relative = '0'
            ms_name = spectral_ds[i].GetDescription()
            if format.upper() == 'VRT':
                if not os.path.isabs(ms_name):
                    ms_relative = '1'
                    ms_name = os.path.relpath(ms_name, os.path.dirname(dst_filename))

            vrt_xml += """    <SpectralBand%s>
          <SourceFilename relativeToVRT="%s">%s</SourceFilename>
          <SourceBand>%d</SourceBand>
        </SpectralBand>\n""" % (dstband, ms_relative, ms_name, sband.GetBand())

        vrt_xml += """  </PansharpeningOptions>\n"""
        vrt_xml += """</VRTDataset>\n"""

        if format.upper() == 'VRT':
            f = gdal.VSIFOpenL(dst_filename, 'wb')
            if f is None:
                print('Cannot create %s' % dst_filename)
                return 1
            gdal.VSIFWriteL(vrt_xml, 1, len(vrt_xml), f)
            gdal.VSIFCloseL(f)
            if verbose_vrt:
                vrt_ds = gdal.Open(dst_filename, gdal.GA_Update)
                vrt_ds.SetMetadata(vrt_ds.GetMetadata())
            else:
                vrt_ds = gdal.Open(dst_filename)
            if vrt_ds is None:
                return 1

            return 0

        vrt_ds = gdal.Open(vrt_xml)
        out_ds = gdal.GetDriverByName(format).CreateCopy(dst_filename, vrt_ds, 0, creation_options,
                                                         callback=progress_callback)
        if out_ds is None:
            return 1
        return 0

    def DoesDriverHandleExtension(self, drv, ext):
        exts = drv.GetMetadataItem(gdal.DMD_EXTENSIONS)
        return exts is not None and exts.lower().find(ext.lower()) >= 0

    def GetExtension(self, filename):
        ext = os.path.splitext(filename)[1]
        if ext.startswith('.'):
            ext = ext[1:]
        return ext


    def GetOutputDriversFor(self, filename):
        drv_list = []
        ext = self.GetExtension(filename)
        for i in range(gdal.GetDriverCount()):
            drv = gdal.GetDriver(i)
            if (drv.GetMetadataItem(gdal.DCAP_CREATE) is not None or
                drv.GetMetadataItem(gdal.DCAP_CREATECOPY) is not None) and \
                    drv.GetMetadataItem(gdal.DCAP_RASTER) is not None:
                if ext and self.DoesDriverHandleExtension(drv, ext):
                    drv_list.append(drv.ShortName)
                else:
                    prefix = drv.GetMetadataItem(gdal.DMD_CONNECTION_PREFIX)
                    if prefix is not None and filename.lower().startswith(prefix.lower()):
                        drv_list.append(drv.ShortName)

        # GMT is registered before netCDF for opening reasons, but we want
        # netCDF to be used by default for output.
        if ext.lower() == 'nc' and not drv_list and \
                drv_list[0].upper() == 'GMT' and drv_list[1].upper() == 'NETCDF':
            drv_list = ['NETCDF', 'GMT']

        return drv_list


    def GetOutputDriverFor(self,
                           filename=None,
                           is_raster=True,
                           default_raster_format='GTiff',
                           default_vector_format='ESRI Shapefile'):
        if not filename:
            return 'MEM'
        drv_list = self.GetOutputDriversFor(filename)
        ext = self.GetExtension(filename)
        if not drv_list:
            if not ext:
                return default_raster_format if is_raster else default_vector_format
            else:
                raise Exception("Cannot guess driver for %s" % filename)
        elif len(drv_list) > 1:
            print("Several drivers matching %s extension. Using %s" % (ext if ext else '', drv_list[0]))
        return drv_list[0]

    def parse_spectral_names(self, dict_spectrac):
        spectral_ds = []
        spectral_bands = []

        for spectral_name in dict_spectrac :
            ds = gdal.Open(spectral_name)
            if ds is None:
                print('读取文件失败【%s】' %(spectral_name))
                return None,None

            # 根据指定波段获取dataset ID
            for band_num in dict_spectrac[spectral_name] :
                band = ds.GetRasterBand(band_num)
                spectral_ds.append(ds)
                spectral_bands.append(band)

        return spectral_ds, spectral_bands
