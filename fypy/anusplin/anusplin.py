# coding:utf-8
'''
@Project: anusplin
-------------------------------------
@File   : anusplin.py
-------------------------------------
@Modify Time      @Author    @Version    
--------------    -------    --------
2021/10/13 9:24     Lee        1.0         
-------------------------------------
@Desciption
-------------------------------------

'''

import os
import sys
import time
import numpy as np
import tempfile
import platform

from osgeo import gdal, osr
from fypy.tools.tifpro import GetGDALType, readtiff
from fypy.tools.BaseAlgorithms import BaseAlgorithms

EXEPath = os.path.dirname(os.path.realpath(__file__))


class anusplin(BaseAlgorithms) :
    # https://www.cnblogs.com/wkfvawl/p/16647422.html

    def __init__(self, lat, lon, data, resolution,
                 vmax=None, vmin=None,
                 extent=None, dict_cov={}):
        self.tempfile = []

        self.srclat  = np.array(lat).reshape(-1)
        self.srclon  = np.array(lon).reshape(-1)
        self.srcdata = np.array(data).reshape(-1)

        self.checkdictcov(dict_cov)

        self.dict_cov = dict_cov

        # 获取自变量范围
        self.vmax = vmax
        self.vmin = vmin

        self.checkExtent(extent)
        self.resolution = resolution

        # 根据插值区域范围计算行列数
        self.rows = int(np.ceil((self.maxY - self.minY) / self.resolution))
        self.cols = int(np.ceil((self.maxX - self.minX) / self.resolution))
        print('='*100)
        print('输出格网--行数: %d , 列数:%d' %(self.rows, self.cols))

    def interp2D(self, prodname ='temp') :

        tmpfile = tempfile.NamedTemporaryFile(prefix="tmp_fypy_anusplin_", delete=True)

        stime = time.time()

        # 参数设置
        self.SetOptions(tmpfile.name, prodname)

        # 对协变量进行标准格式处理
        covvarlist = []
        for key in self.dict_cov :
            covvar = self.covdataformat(self.covname + self.dict_cov[key]['covname'],
                                        self.dict_cov[key], self.srclat, self.srclon)
            self.tempfile.append(self.covname + self.dict_cov[key]['covname'])
            covvarlist.append(covvar)

        self.initdata(self.srclat, self.srclon, self.srcdata, covvarlist)

        # 添加初始节点
        if not self.selnot(tmpfile.name, os.path.join(EXEPath, 'bin')) :
            return None

        # 任意个独立变量或多个协变量的薄盘样条函数。数据平滑度由 GCV 或 GML 决定。
        # 最多 10000 个站点，2000 个节点。
        if  not self.splinb(tmpfile.name, os.path.join(EXEPath, 'bin')) :
            return None

        # 生成拟合曲面或贝叶斯标准误差曲面
        if not self.lapgrd(tmpfile.name, os.path.join(EXEPath, 'bin')) :
            return None

        etime = time.time()
        print('ANSPLIN共耗时【%.2f秒】' %(etime - stime))

        # 将生成的grib文件转换成TIFF
        return self.GetGrid(self.gridname)

    def SetOptions(self, outdir, prodname):

        covcount= len(self.dict_cov)

        self.prodname = prodname
        self.units = 0           # units code for surface value and fillvalue
        # 0:undefined    1:meters  2:feet
        # 3:kilometers   4:miles   5:degrees
        # 6:radians      7:millimetres   8:megajoules
        self.independentvar = 2  # 自变量（x, y ）, 最多10个
        self.covcount   = covcount  # 协变量 dem
        self.sfc = 0             # 表面样条变量个数， 每个表面采用不同的样条变量
        self.sfcc = 0            # 表面协变量个数， 每个表面采用不同的样条变量
        # 对自变量范围设置
        # 每一个自变量的上下限、有无转换、单位，该项为拟合区域的X的最大值最小值坐标，即最左、最右边界的X值；
        # 0：不转换，1：转换；
        self.exchange = 0

        # units--> 1：m, 5:degree
        self.posunits = 0

        # 对协变量范围设置
        # 高程的上下限，
        # self.covvar1_min = -500
        # self.covvar1_max = 9000
        # self.covvar1_units = 0

        # 自变量转换参数，0代表不转换
        self.independent_var_exchange = 1
        self.covariate_var_exchange = 0

        # 样条次数
        self.splinenum = 3

        # 输出表面个数
        self.sfcnum = 1

        # 所有的表面用相同的权重
        # 0：不需要
        # 1：采用不同权重，需要对各个表面所用权重进行设置
        self.weightforsfc = 0

        # 优化参数指标
        self.optimalparm = 1

        # 平滑参数选择方法
        # 1 ：GCV，2：MSE，3：固定值，4：GML
        self.smoothparm = 1

        # 输入数据
        self.datadir = outdir + '.in'

        # 样本总数
        self.samplenum = 0
        self.knotsnum = int(self.samplenum * 0.5)
        # 站点的label
        self.stationlabel = 0

        # 输入数据的变量存储格式
        self.format = '(%df10.3)' %(covcount+3)  # lat, lon, val, cov1,....

        # 协变量文件
        self.covname = outdir + '.cov'

        # 输出文件
        # 残差文件
        self.resname = outdir + '.res'

        # 表面文件，用于下一步LAPGRD的输入数据
        self.surname = outdir + '.sur'
        self.notname = outdir + '.not'

        self.ourdir  = outdir + '.out'
        self.rejname = outdir + '.rej'
        self.optname = outdir + '.opt'
        ############################################################
        ############################################################

        # 表面个数
        self.sfcnum = 1

        # 1：转换恢复
        self.re_exchange = 1

        # 误差类型
        # 1：模型标准误差，2：预测标准误差
        self.biastype = 1

        # 掩摸方式
        # 0：不提供掩摸文件
        self.masktype = 0

        # 独立的协变量数据格式，2代表输入输出格式大小一致
        self.covformat = 2

        # 空值数据
        self.fillvalue = -9999.0

        # 输出文件名
        self.gridname = outdir + '.grd'

        # 输出数据格式
        self.gridformat = '(%df15.2)' %(self.cols)

        self.tempfile.append(self.datadir)
        self.tempfile.append(self.covname)
        self.tempfile.append(self.resname)
        self.tempfile.append(self.surname)
        self.tempfile.append(self.notname)
        self.tempfile.append(self.ourdir)
        self.tempfile.append(self.rejname)
        self.tempfile.append(self.gridname)
        self.tempfile.append(self.optname)

        # 清理已存在的临时文件
        self.clean()

    def selnot(self, outdir, exedir):

        # step1
        selnot = os.path.join(exedir, 'selnot.exe')
        confname = outdir + '_selnot.conf'
        logname  = outdir + '_selnot.log'

        self.tempfile.append(confname)
        self.tempfile.append(logname)

        self.createConfForSelnot(confname)

        cmd1 = '{exename} <{confname}> {logname}'.format(
            exename=selnot, confname=confname, logname=logname)

        # 在Linux系统运行需要wine
        if platform.system().lower() != 'windows' :
            cmd1 = 'wine ' + cmd1

        print(cmd1)
        status=os.system(cmd1)
        if os.path.isfile(self.notname) :
            print('【selnot】处理成功')
            return True
        else:
            print('【selnot】处理失败')
            return False

    def splinb(self, outdir, exedir):

        # step2
        splinb = os.path.join(exedir, 'splinb.exe')
        confname = outdir + '_splinb.conf'
        logname  = outdir + '_splinb.log'
        self.tempfile.append(confname)
        self.tempfile.append(logname)
        self.createConfForSplineB(confname)

        cmd2 = '{exename} <{confname}> {logname}'.format(
            exename=splinb, confname=confname, logname=logname)
        # 在Linux系统运行需要wine
        if platform.system().lower() != 'windows' :
            cmd2 = 'wine ' + cmd2

        print(cmd2)
        status=os.system(cmd2)
        if os.path.isfile(self.surname) :
            print('【splinb】处理成功')
            return True
        else:
            print('【splinb】处理失败')
            return False

    def lapgrd(self, outdir, exedir):

        # step3
        lapgrd = os.path.join(exedir, 'lapgrd.exe')
        confname = outdir + '_lapgrd.conf'
        logname  = outdir + '_lapgrd.log'
        self.tempfile.append(confname)
        self.tempfile.append(logname)

        self.createConfForLapgrd(confname)
        cmd3 = '{exename} <{confname}> {logname}'.format(
            exename=lapgrd, confname=confname, logname=logname)
        # 在Linux系统运行需要wine
        if platform.system().lower() != 'windows' :
            cmd3 = 'wine ' + cmd3

        print(cmd3)
        status=os.system(cmd3)
        if os.path.isfile(self.gridname) :
            print('【lapgrd】处理成功')
            return True
        else:
            print('【lapgrd】处理失败')
            return False

    def createConfForSelnot(self, confname):
        '''
        为SPLINB设置knots

        Parameters
        ----------
        confname : str
            配置文件名

        Returns
        -------

        '''
        with open(confname, 'w') as fp :
            fp.write('%d\n' %(self.independentvar))   # 独立样条变量的个数（lat, lon）                         # 1
            fp.write('%d\n' % (self.covcount))     # 独立协变量的个数                                      # 2
            fp.write('%d\n' %(self.sfc))              # 表面独立样条变量的个数                                 # 3
            fp.write('%d\n' %(self.sfcc))             # 表面独立协变量数                                    # 4

            # 自变量（lat,lon）下限和上限，转换代码，单位代码，可选边际
            fp.write('%.3f %.3f %d %d\n' % (self.minX, self.maxX, self.exchange, self.posunits))      # 5
            fp.write('%.3f %.3f %d %d\n' % (self.minY, self.maxY, self.exchange, self.posunits))      # 6
            for key in self.dict_cov :
                fp.write('%.3f %.3f %d %d\n' %(self.dict_cov[key]['vmin'], self.dict_cov[key]['vmax'], 1, 0))
                fp.write('%d\n' %(self.independent_var_exchange))                                               # 8
            fp.write('%d\n' %(self.covariate_var_exchange))                                                 # 9
            fp.write('%d\n' %(self.sfcnum))            # 表面的数量                                          # 10
            fp.write('%d\n' %(self.weightforsfc))      # 相对误差方差                                        # 11
            fp.write('%s\n' %(self.datadir))           # 数据文件名                                          # 12
            fp.write('%d\n' %(self.samplenum))         # 数据最多个数                                         # 13
            fp.write('%d\n' %(self.stationlabel))      #  站点的名称                                         # 14
            fp.write('%s\n' %(self.format))            # 数据存储格式                                        # 15
            fp.write('%s\n' %(self.notname))           # 输出knot的文件名                                   # 16
            fp.write('%s\n' %(self.rejname))           # 输出被拒点的文件名                                  # 17
            fp.write('%d\n' %(self.knotsnum))          # knots点数                                            # 18

    def createConfForSplineB(self, confname):
        '''任意个独立变量或多个协变量的薄盘样条函数。数据平滑度由 GCV 或 GML 决定。最多 10000 个站点'''
        with open(confname, 'w') as fp :
            fp.write('%s\n' %(self.prodname))       # 表面变量的标题                                        # 1
            # 表面变量的单位:
            # -----------------------------------------------
            # 0 - UNDEFINED        5 - DEGREES
            # 1 - METRES           6 - RADIANS
            # 2 - FEET             7 - MILLIMETRES
            # 3 - KILOMETRES       8 - MEGAJOULES5
            # 4 - MILES
            fp.write('%d\n' %(self.units))                                                           # 2
            # INDEPENDENT VARIABLES
            # ---------------------
            # 独立变量的个数（lat, lon） (0 TO 10)
            fp.write('%d\n' %(self.independentvar))                                                  # 3
            # 协变量的个数 (0 TO  8)
            fp.write('%d\n' % (self.covcount))                                                    # 4
            # 表面样条变量数 (0 TO  7)
            fp.write('%d\n' %(self.sfc))                                                             # 5
            # 表面协变量数 (0 TO  7)
            fp.write('%d\n' %(self.sfcc))                                                            # 6
            #     转换代码                     REFERENCE UNIT CODES
            # --------------------            --------------------
            # 0 - NO TRANSFORMATION           0 - UNDEFINED
            # 1 - X/A                         1 - METRES
            # 2 - X*A                         2 - FEET
            # 3 - A*LOG(X + B)                3 - KILOMETRES
            # 4 - (X/B)**A                    4 - MILES
            # 5 - A*EXP(X/B)                  5 - DEGREES
            # 6 - A*TANH(X/B)                 6 - RADIANS
            # 7 - ANISOTROPY ANGLE            7 - MILLIMETRES
            # 8 - ANISOTROPY FACTOR           8 - MEGAJOULES

            # LOWER & UPPER LIMITS, TRANSF CODE, REF UNIT, MARGIN(S) FOR VARIABLE
            fp.write('%.3f %.3f %d %d\n' % (self.minX, self.maxX, self.exchange, self.posunits))         # 7
            fp.write('%.3f %.3f %d %d\n' % (self.minY, self.maxY, self.exchange, self.posunits))         # 8
            for key in self.dict_cov :
                fp.write('%.3f %.3f %d %d\n' %(self.dict_cov[key]['vmin'], self.dict_cov[key]['vmax'], 1, 0))         # 9
                # ENTER  1 TRANSFORMATION COEFFICIENT(S)
                # SURFACE DIRECTIVES
                # ------------------
                # DEPENDENT VARIABLE TRANSFORMATION:
                # 0 - NO TRANSFORMATION
                # 1 - NATURAL LOGARITHM
                # 2 - SQUARE ROOT
                # 5 - OCCURRENCE0
                fp.write('%d\n' %(self.independent_var_exchange))                                        # 10
            fp.write('%d\n' %(self.covariate_var_exchange))                                          # 11
            # 样条数 (AT LEAST 2)
            fp.write('%d\n' %(self.splinenum))                                                       # 12
            # 表面变量数 (AT LEAST 1)
            fp.write('%d\n' %(self.sfcnum))                                                          # 13
            # 相对误差变量数 (0 OR 1)
            fp.write('%d\n' %(self.weightforsfc))                                                    # 14
            # OPTIMIZATION DIRECTIVE (NORMALLY 1):
            #   0 - COMMON SMOOTHING PARAMETER FOR ALL SURFACES
            #   1 - COMMON SMOOTHING DIRECTIVE FOR ALL SURFACES
            #   2 - DIFFERENT SMOOTHING DIRECTIVE FOR EACH SURFACE1
            fp.write('%d\n' %(self.optimalparm))                                                     # 15
            # SMOOTHING DIRECTIVE (NORMALLY 1):
            #   0 - FIXED SMOOTHING PARAMETER FOR EACH SURFACE
            #   1 - MINIMIZE GCV FOR EACH SURFACE
            #   2 - MINIMIZE TRUE MEAN SQUARE ERROR FOR EACH SURFACE
            #   3 - FIXED SIGNAL FOR EACH SURFACE
            #   4 - MINIMIZE GML FOR EACH SURFACE1
            fp.write('%d\n' %(self.smoothparm))                                                      # 16
            # DATA FILE NAME
            fp.write('%s\n' %(self.datadir))                                                         # 17
            # MAXIMUM NUMBER OF DATA POINTS (AT LEAST  7)
            fp.write('%d\n' %(self.samplenum))                                                       # 18
            # NO. OF CHARACTERS IN SITE LABEL (0 TO 20)
            fp.write('%d\n' %(self.stationlabel))                                                    # 19
            # DATA FORMAT (SITE LABEL,  3 INDEP VARS,  1 SURFACES,  0 REL VARIANCES)
            fp.write('%s\n' %(self.format))                                                          # 20
            # KNOT INDEX FILE NAME
            fp.write('%s\n' %(self.notname))                                                         # 21
            # MAXIMUM NUMBER OF DATA KNOTS (AT LEAST  1)
            fp.write('%d\n' %(self.knotsnum))                                                            # 22
            # INPUT BAD DATA FLAG FILE NAME (BLANK IF NOT REQUIRED)
            fp.write('\n')                                                                      # 23
            # OUTPUT BAD DATA FLAG FILE NAME (BLANK IF NOT REQUIRED)
            fp.write('\n')                                                                      # 24
            # OUTPUT LARGE RESIDUAL FILE NAME (BLANK IF NOT REQUIRED)
            fp.write('%s\n' %(self.resname))                                                                      # 25
            # OUTPUT OPTIMIZATION PARAMETERS FILE NAME (BLANK IF NOT REQUIRED):
            fp.write('%s\n' %(self.optname))                                                                      # 26
            # OUTPUT SURFACE COEFFICIENTS FILE NAME (BLANK IF NOT REQUIRED):
            fp.write('%s\n' %(self.surname))                                                         # 27
            # OUTPUT DATA LIST FILE NAME (BLANK IF NOT REQUIRED):
            fp.write('\n')                                                                      # 28
            # OUTPUT ERROR COVARIANCE FILE NAME (BLANK IF NOT REQUIRED):
            fp.write('\n')                                                                      # 29
            # OUTPUT ERROR COVARIANCE FILE NAME (BLANK IF NOT REQUIRED):
            fp.write('\n')                                                                      # 30
            # fp.write('%d\n' %(self.samplenum))                                                       # 31
            # fp.write('%d\n' %(self.stationlabel))                                                    # 32
            # fp.write('%s\n' %(self.format))                                                       # 33
            # fp.write('%s\n' %(self.ourdir))                                                          # 34

    def createConfForLapgrd(self, confname):
        with open(confname, 'w') as fp :
            fp.write('%s\n' %(self.surname))
            fp.write('%d\n' %(self.sfcnum))
            fp.write('%d\n' %(self.re_exchange))
            fp.write('\n')
            fp.write('%d\n' %(self.biastype))
            fp.write('%d\n' %(1))
            fp.write('%.3f %.3f %.4f\n' % (self.minX, self.maxX, self.resolution))
            fp.write('%d\n' %(2))
            fp.write('%.3f %.3f %.4f\n' % (self.minY, self.maxY, self.resolution))
            fp.write('%d\n' %(self.masktype))
            fp.write('%d\n' %(self.covformat))
            for key in self.dict_cov :
                fp.write('%s\n' % (self.covname + self.dict_cov[key]['covname']))
            # fp.write('%s\n' % (self.covname))
                fp.write('%d\n' %(2))                   # arcgis grid
            fp.write('%f\n' %(self.fillvalue))
            fp.write('%s\n' %(self.gridname))
            fp.write('%s\n' %(self.gridformat))

    def initdata(self, srclat, srclon, srcdata, covvarlist):
        '''
        按格式生成输入数据ASCII文件

        Args:
            srclat: array, 纬度
            srclon: array, 经度
            srcdata: array, 自变量

        Returns:

        '''

        # 获取指定区域范围内的输入数据
        flag = (srclon <= self.maxX) & (srclon >= self.minX) \
               & (srclat <= self.maxY) & (srclat >= self.minY)
        if self.vmax is not None and self.vmin is not None :
            flag = (srcdata >= self.vmin) & (srcdata <= self.vmax) \
                   & flag

        mdata = srcdata[flag]
        mlon = srclon[flag]
        mlat = srclat[flag]

        # 如果未指定数据的有效值范围，将自动根据输入获取最值
        if self.vmax is  None :
            self.vmax = np.nanmax(mdata)

        if  self.vmin is None :
            self.vmin = np.nanmin(mdata)
        print('数据有效值范围vmin: %.2f  vmax:  %.2f' %(self.vmin, self.vmax))
        print('区域范围【lon: %.2f ~ %.2f】，【lat: %.2f ~ %.2f】' %(self.minX, self.maxX,
                                                                   self.minY, self.maxY))
        # 对输入数据进行标准格式化处理
        Line = int(np.ceil((self.maxY - self.minY) / self.resolution))
        Pixel = int(np.ceil((self.maxX - self.minX) / self.resolution))

        row = np.array((self.maxY - mlat) / self.resolution, dtype=np.int32)
        col = np.array((mlon - self.minX) / self.resolution, dtype=np.int32)
        row[row>=Line] = Line-1
        col[col>=Pixel] = Pixel-1

        row[row<0] = 0
        col[col<0] = 0

        # ####################################################################################
        # # 对输入数据进行输出临时文件，以便后续分析
        # ####################################################################################
        # GlobleData = np.full(shape=(Line, Pixel), dtype=mdata.dtype, fill_value=-1)
        # GlobleData[row, col] = mdata
        # trans = [self.minX, self.resolution, 0,
        #          self.maxY, 0, -self.resolution]
        # # writetiff(self.outname+'_srcdata.tif', GlobleData, trans, fillvalue=-1)
        # del GlobleData

        # 将来优化该模块
        Index = np.arange(len(mdata))
        GlobleData = np.full(shape=(Line, Pixel), dtype=np.int32, fill_value=-1)
        GlobleData[row, col] = Index

        flag1  = GlobleData[GlobleData>=0]
        mdata  = mdata[flag1]
        mlon   = mlon[flag1]
        mlat   = mlat[flag1]

        result = []
        result.append(mlon)
        result.append(mlat)

        for i in range(self.covcount) :
            srccovvar = covvarlist[i]
            covvar = srccovvar[flag]
            covvar = covvar[flag1]
            covvar[covvar<0] = 0
            result.append(covvar)

        result.append(mdata)
        self.samplenum = len(mdata)
        print('samplenum:', self.samplenum)

        # Note:
        # 如果输入控制点数据过多，将导致计算量过大，将控制在3000个点以内
        self.knotsnum = int(self.samplenum * 0.75)
        if self.knotsnum > 3000:
            self.knotsnum = 3000

        result = np.array(result, dtype=np.float64).T

        np.savetxt(self.datadir, result, fmt='%10.3f',  comments='')

        print('成功格式化输入数据')

    def covdataformat(self, covariatename, dict_covdata=None, srclat=None, srclon=None):
        '''
        对高程数据进行格式化处理

        Args:
            covariatename: 协变量文件名
            demfile: 高程静态数据

        Returns:

        '''
        # covdata = readnc(demfile, 'dem')
        # srcdem,_,_ = readtiff(demfile)
        covdata = np.array(dict_covdata['data'])
        minlon = dict_covdata['extent'][0]
        minlat = dict_covdata['extent'][1]
        maxlon = dict_covdata['extent'][2]
        maxlat = dict_covdata['extent'][3]

        # res = (maxY - minY) / line
        res = (maxlat-minlat) / covdata.shape[0]
        Line = covdata.shape[0]
        Pixel = covdata.shape[1]

        lon, lat = np.meshgrid(np.arange(self.minX, self.maxX, self.resolution),
                               np.arange(self.maxY, self.minY, -self.resolution))

        row = np.array((maxlat - lat) / res, dtype=np.int32)
        col = np.array((lon - minlon) / res, dtype=np.int32)
        row[row>=Line] = Line-1
        row[row<0] = 0

        col[col>=Pixel] = Pixel-1
        col[col<0] = 0

        dem = covdata[row, col]
        if dem is not None :
            rows, cols = dem.shape

            header = 'ncols         %d\n' \
                     'nrows         %d\n' \
                     'xllcorner     %f\n' \
                     'yllcorner     %f\n' \
                     'cellsize      %f\n' \
                     'NODATA_value  %f' %(cols, rows, self.minX, self.minY,
                                          self.resolution, self.fillvalue)

            np.savetxt(covariatename, dem, fmt='%4d', header=header, comments='')

        row = np.array((maxlat - srclat) / res, dtype=np.int32)
        col = np.array((srclon - minlon) / res, dtype=np.int32)
        row[row>=Line] = Line-1
        row[row<0] = 0

        col[col>=Pixel] = Pixel-1
        col[col<0] = 0
        covvar = covdata[row, col]

        del covdata

        return covvar

    def GetGrid(self, gridname):
        '''
        将模型输出的grid文件转换成tif文件

        Args:
            outname: 输出的TIFF文件名
            gridname: 输入的grid文件名

        Returns:

        '''
        if not os.path.isfile(gridname) :
            print('格点文件不存在，转换失败【%s】' %(gridname))
            return None

        try:
            dict_data = self.readgrid(gridname)
            toplat = dict_data['bottomlat'] + dict_data['rows'] * dict_data['resolution']
            trans = [dict_data['leftlon'], dict_data['resolution'], 0,
                     toplat, 0, -dict_data['resolution']]
            data = dict_data['data']
            # data[data < self.vmin] = self.vmin
            # data[data > self.vmax] = self.vmax
            # flag = (data < self.vmin) | (data > self.vmax)
            # data[flag] = dict_data['fillvalue']
            fillvalue = dict_data['fillvalue']
        except BaseException as e :
            print(e)
            return None

        data = np.array(data)
        datatype = GetGDALType(data.dtype)

        if len(data.shape) == 3:
            im_bands, im_height, im_width = data.shape
        elif len(data.shape) == 2:
            im_bands, (im_height, im_width) = 1,data.shape
        else:
            im_bands, (im_height, im_width) = 1,data.shape
            #创建文件
        driver = gdal.GetDriverByName("MEM")
        dstDS = driver.Create('', im_width, im_height, im_bands, datatype)

        dstDS.SetGeoTransform(trans) #写入仿射变换参数
        # 获取地理坐标系统信息，用于选取需要的地理坐标系统
        sr = osr.SpatialReference()
        sr.SetWellKnownGeogCS("EPSG:4326")
        # 给新建图层赋予投影信息
        dstDS.SetProjection(sr.ExportToWkt())

        if im_bands == 1:
            dstDS.GetRasterBand(1).WriteArray(data)
        else:
            for i in range(im_bands):
                dstDS.GetRasterBand(i+1).WriteArray(data[i])

        if fillvalue is not None :
            dstDS.GetRasterBand(1).SetNoDataValue(float(fillvalue))

        return dstDS

    def readgrid(self, gridname):

        dict_data = {}

        with open(gridname, 'r') as fp :
            strNCOLS = fp.readline()
            strNROWS = fp.readline()
            strXLLCORNER = fp.readline()
            strYLLCORNER = fp.readline()
            strCELLSIZE = fp.readline()
            strNODATA_VALUE = fp.readline()

        dict_data['cols'] = np.int32(strNCOLS.split()[1])
        dict_data['rows'] = np.int32(strNROWS.split()[1])
        dict_data['leftlon'] = np.float32(strXLLCORNER.split()[1])
        dict_data['bottomlat'] = np.float32(strYLLCORNER.split()[1])
        dict_data['resolution'] = np.float32(strCELLSIZE.split()[1])
        dict_data['fillvalue'] = np.float32(strNODATA_VALUE.split()[1])


        data = np.loadtxt(gridname, dtype=np.float32, skiprows=6)

        # data[(data<390) | (data > 430)] = np.nan

        if (data.shape[0] != dict_data['rows']) | (data.shape[1] != dict_data['cols']) :
            raise Exception('The Row: %d, The Col: %d, but the data region size is [%d, %d]' %(
                dict_data['rows'], dict_data['cols'], data.shape[0], data.shape[1]))

        dict_data['data'] = data

        return dict_data

    def clean(self):
        for filename in self.tempfile :
            if os.path.isfile(filename) :
                try:
                    os.remove(filename)
                except BaseException as e :
                    print('删除失败：%s' %(filename))

    def __del__(self):

        self.clean()



def setcovdem():
    demfile = os.path.join(EXEPath, 'conf', 'Global_DEM_50KM.TIF')
    srcdem, trans, prj = readtiff(demfile)

    dict_cov = {
        'cov1' : {
            'covname' : 'dem',
            'data' : np.zeros_like(srcdem, dtype=np.float32),
            'extent' : [-180.0, -90.0, 180.0, 90.0], #(minX, minY, maxX, maxY)
            'vmin' : 0.0,
            'vmax' : 10000.0
        },
    }

    return dict_cov





