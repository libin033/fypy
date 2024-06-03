# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : gfscene.py

@Modify Time :  2024/4/15   

@Author      : Lee    

@Version     : 1.0   

@Description :

'''
import glob
import os
import sys
import numpy as np
import datetime
from osgeo import gdal
import xml.dom.minidom
from tqdm import tqdm
import shutil
import json

from dateutil.relativedelta import relativedelta
from fypy.tools.BaseAlgorithms import BaseAlgorithms
from fypy import parm

EXEPATH = os.path.abspath(list(parm.__path__)[0])


class GFScene(BaseAlgorithms) :

    def __init__(self):
        pass

    def calibrate(self, srcfile, nowdate,
                  SatID, InstID, InstType=None,
                  rpcfile=None, demfile=None,
                  fillvalue=65535, blocksize=1000):
        ''' 对高分卫星数据进行辐射定标 '''

        if rpcfile is not None :
            srcdset = self.ortho_rectification(srcfile, rpcfile=rpcfile, demfile=demfile)
        else:
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
        dstdset = driver.Create('', cols, rows, Bands, gdal.GDT_Float32,)
                                   # options=["COMPRESS=LZW", "BigTIFF=YES"])
        dstdset.SetGeoTransform(trans)
        dstdset.SetProjection(prj)

        #分别读取4个波段
        for BandID in range(0, Bands):
            ReadBand = srcdset.GetRasterBand(BandID+1)
            outband = dstdset.GetRasterBand(BandID+1)
            outband.SetNoDataValue(fillvalue)
            #获取对应波段的增益gain和偏移bias
            Gain, Bias = self.RadiometricCalibration(nowdate, SatID, InstID, BandId=BandID+1)

            if 'PMS' in InstID :
                if InstType is None :
                    raise Exception('请传入InstType参数【全色：PAN，多光谱：MSS】')

                if 'MSS' in InstType :
                    BandID += 1
                elif 'PAN' in InstType :
                    BandID = BandID
                else:
                    raise Exception('请传入InstType参数【全色：PAN，多光谱：MSS】')

            i = 0
            j = 0
            #进度条参数
            XBlockcount = np.ceil(cols / blocksize)
            YBlockcount = np.ceil(rows / blocksize)
            # print("第%d波段校正："%BandID)
            try:
                with tqdm(total=XBlockcount*YBlockcount, iterable='iterable',
                          desc = '正在进行第%i波段定标' %(BandID+1), mininterval=1) as pbar:
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

                            atcImage =np.where(Image>0, (Image*Gain + Bias)*0.01, fillvalue)

                            outband.WriteArray(atcImage,j,i)
                            j=j+nXBK
                            pbar.update(1)
                        j=0
                        i=i+nYBK
            except BaseException as e :
                print(e)
                pbar.close()
                continue
            pbar.close()

        srcdset = None

        return dstdset

    def RadiometricCalibration(self, nowdate, SatID, InstID, BandId):
        ''' 获取高分卫星定标系数 gain和 bias '''

        satid = SatID.replace('-', '')
        instid = InstID.replace('-', '')

        #读取辐射校正所需参数:增益、偏移和光谱响应函数
        CalCoeffFile = os.path.join(EXEPATH, "coef", '%s_%s.json' %(satid.upper(), instid.upper()))
        if not os.path.isfile(CalCoeffFile) :
            coeffiles = glob.glob(os.path.join(EXEPATH, "coef", 'GF*.json'))
            print('定标系数文件不存在【%s】，仅支持-->' %(CalCoeffFile), coeffiles)
            return None

        CalCoeff = json.load(open(CalCoeffFile))

        while nowdate >= datetime.datetime.strptime('2000', '%Y') :
            stryear = nowdate.strftime('%Y')
            if not stryear in CalCoeff['B%02d' %(BandId)] :
                nowdate -= relativedelta(years=1)
            else:
                gain = CalCoeff['B%02d' %(BandId)][stryear]['gain']
                bias = CalCoeff['B%02d' %(BandId)][stryear]['bias']

                return gain, bias

        # raise Exception('未匹配到该年【%s】系数' %(nowdate.strftime('%Y')))
        return None, None

    def ortho_rectification(self, srcfile, rpcfile, demfile=None, dstfile=None):
        """
        正射校正
        :param input:输入原始影像
        :param output:输出正射影像
        """

        if not os.path.isfile(rpcfile) :
            raise Exception('RPC文件不存在【%s】' %(rpcfile))

        # 保证rpc文件与源文件在同一个文件夹下
        srcdir = os.path.dirname(srcfile)
        dstrpcfile = os.path.join(srcdir, os.path.basename(rpcfile))
        if not os.path.isfile(dstrpcfile) :
            shutil.copy(rpcfile, srcdir)

        dataset = gdal.Open(srcfile, gdal.GA_Update)    # 读入影像
        rpc = dataset.GetMetadata("RPC")                # 读入影像，rpc

        if dstfile is None :
            format = 'MEM'
            dstfile = ''
            creationOptions=[]
        else:
            format = 'GTiff'
            dstfile = dstfile
            creationOptions=["COMPRESS=LZW", "BigTIFF=YES"]

        if demfile is None :
            dst_ds = gdal.Warp(dstfile, dataset,
                               dstSRS='EPSG:4326',
                               format=format,
                               resampleAlg=gdal.GRA_NearestNeighbour,
                               rpc=True, #使用RPC模型进行校正
                               warpOptions=['INIT_DEST=NO_DATA'],
                               creationOptions=creationOptions)
        else :
            dst_ds = gdal.Warp(dstfile, dataset,
                               dstSRS='EPSG:4326',
                               format=format,
                               resampleAlg=gdal.GRA_NearestNeighbour,
                               rpc=True, #使用RPC模型进行校正
                               transformerOptions=[r'RPC_DEM=%s' %(demfile),
                                                   "RPC_DEMINTERPOLATION=bilinear"],
                               creationOptions=creationOptions)

        return dst_ds

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

def getmw(wl, resp) :
    resp = np.array(resp)
    newwl = np.arange(wl[0], wl[1]+0.0002, 0.0025)

    index = np.argmin(1-resp)

    return newwl[index]

if __name__ == '__main__':

    from fypy.tools.jsonpro import readjson, writejson
    ccfile = r'D:\pypi\fypy\fypy\RSDP\resp\bak\GF\CalibrationCoefficient.json'
    rsffile = r'D:\pypi\fypy\fypy\RSDP\resp\bak\GF\GF.json'

    data1 = readjson(ccfile)
    data2 = readjson(rsffile)

    for satid in data1 :
        for instid in data1[satid] :
            dict_info = {}
            info1 = data1[satid][instid]
            if 'ESUN' in info1 :
                bandcount = len(info1['ESUN'])
            elif instid in ['IRS'] :
                bandcount = 1
            else:
                bandcount = 4
            for i in range(bandcount) :
                bandid = i + 1
                print(satid, instid, bandid)
                dict_info['B%02d' %(bandid)] = {}

                if 'ESUN' in info1 :
                    sun = info1['ESUN'][i]
                    dict_info['B%02d' %(bandid)]['ESUN'] = float('%.4f' %(sun))

                for item in info1 :
                    if item in ['ESUN', 'BMSI'] :
                        continue
                    if satid == 'GF4' and instid == 'PMS' :
                        dict_info['B%02d' %(bandid)][item] = {}
                        for item1 in info1[item] :
                            dict_info['B%02d' %(bandid)][item][item1] = {}
                            dict_info['B%02d' %(bandid)][item][item1]['gain'] = info1[item][item1]['gain'][i]
                            dict_info['B%02d' %(bandid)][item][item1]['bias'] = info1[item][item1]['bias'][i]
                    else:
                        dict_info['B%02d' %(bandid)][item] = {}
                        dict_info['B%02d' %(bandid)][item]['gain'] = info1[item]['gain'][i]
                        dict_info['B%02d' %(bandid)][item]['bias'] = info1[item]['bias'][i]

                if not satid in data2 :
                    print('缺失卫星【%s】光谱信息' %(satid))
                    continue

                if not instid in data2[satid] :
                    print('【%s】【%s】缺失光谱信息' %(satid, instid))
                else:
                    info2 = data2[satid][instid]['B%d' %(bandid)]

                    mw = getmw(info2['wl'], info2['SRF'])

                    dict_info['B%02d' %(bandid)]['MW'] = float('%.3f' %(mw))

                    wl = list([float('%.4f' %(info2['wl'][0])), float('%.4f' %(info2['wl'][1]))])
                    dict_info['B%02d' %(bandid)]['wl'] = wl
                    dict_info['B%02d' %(bandid)]['SRF'] = info2['SRF']
            outjson = '%s_%s.json' %(satid.upper(), instid.upper())
            writejson(outjson, dict_info)



