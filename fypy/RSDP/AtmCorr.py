# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : AtmCorr.py

@Modify Time :  2022/12/23 15:43   

@Author      : Lee    

@Version     : 1.0   

@Description :
https://www.l3harrisgeospatial.com/docs/FLAASH.html
'''
import os
import sys
import numpy as np
import datetime
from osgeo import gdal

from fypy.py6s import SixS, Geometry, AtmosProfile, AeroProfile, \
    Altitudes, Wavelength, PredefinedWavelengths, AtmosCorr, GroundReflectance


class AtmCorr() :

    def __init__(self):
        # 6S模型
        self.s = SixS()

    def set_geom(self, nowdate, sunz, suna, satz, sata):
        # 第一步：设置几何参数
        self.s.geometry = Geometry.User()
        self.s.geometry.solar_z = sunz
        self.s.geometry.solar_a = suna
        self.s.geometry.view_z = satz
        self.s.geometry.view_a = sata

        self.s.geometry.month = nowdate.month
        self.s.geometry.day = nowdate.day

    def set_atm(self, sLatitude, nowdate):

        self.s.atmos_profile = AtmosProfile.PredefinedType(
            int(AtmosProfile.FromLatitudeAndDate(sLatitude, nowdate.strftime('%Y-%m-%d %H:%M:%S'))))
        return
        # 大气模式类型
        if sLatitude > -15 and sLatitude <= 15:
            s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.Tropical)

        if sLatitude > 15 and sLatitude <= 45:
            if s.geometry.month > 4 and s.geometry.month <= 9:
                s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeSummer)
            else:
                s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeWinter)

        if sLatitude > 45 and sLatitude <= 60:
            if s.geometry.month > 4 and s.geometry.month <= 9:
                s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.SubarcticSummer)
            else:
                s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.SubarcticWinter)

    def set_aer(self, type=None):
        '''
        6s模型中气溶胶类型

        主要根据影像中地物类型分为：
            NoAerosols、Continental、Maritime、Urban、
            Desert、BiomassBurning、Strat-ospheric。
        在自动化程序中统一选择continental气溶胶类型。

        Parameters
        ----------
        s
        type

        Returns
        -------

        '''
        # 气溶胶类型大陆
        self.s.aero_profile = AtmosProfile.PredefinedType(AeroProfile.Continental)

        # 目标地物
        self.s.ground_reflectance = GroundReflectance.HomogeneousLambertian(0.36)

    def set_vis(self, vis=40, opt=0.14497):
        '''
        能见度指的是影像获取当天大气条件，
        也可以使用550nm气溶胶光学厚度代表，
        实际的气溶胶厚度或能见度难以获得，
        在FLAASH中设置为默认的40km，
        在使用6s模型时与FLAASH保持一致都设置为40km，
        对应的550nm气溶胶厚度为0.14497,
        就是只能用气溶胶厚度，而不能用能见度。

        Parameters
        ----------
        s
        vis
        opt

        Returns
        -------

        '''
        # 550nm气溶胶光学厚度,根据日期从MODIS处获取。
        self.s.visibility=vis
        # self.s.aot550 = opt

    def set_altitude(self, dem):
        '''
        影像的平均高程可以借助DEM和影像头文件中的位置信息计算得到。
        6s模型传感器高度根据遥感平台分为三类：
            * -1000代表卫星传感器；
            * 0代表地面传感器；
            * 0～100km代表航空遥感平台的传感器如无人机、自由气球等。
        在自动化工程中选择卫星传感器。

        Parameters
        ----------
        s
        dem

        Returns
        -------

        '''
        # 研究区海拔、卫星传感器轨道高度
        self.s.altitudes = Altitudes()
        self.s.altitudes.set_target_custom_altitude(dem)
        self.s.altitudes.set_sensor_satellite_level()

    def set_resp(self, minwl=None, maxwl=None, response=None, wavelength=None):
        # 校正波段（根据波段名称）
        if wavelength is None :
            self.s.wavelength = Wavelength(minwl, maxwl, response)
        else:
            self.s.wavelength = Wavelength(wavelength)

        # 下垫面非均一、朗伯体
        self.s.atmos_corr = AtmosCorr.AtmosCorrLambertianFromReflectance(-0.1)

    def corrCoeff(self):
        # 运行6s大气模型
        self.s.run()

        a = self.s.outputs.coef_xa
        b = self.s.outputs.coef_xb
        c = self.s.outputs.coef_xc
        # x = s.outputs.values
        # measured radiance [w/m2/sr/mic]
        # y=xa*(measured radiance)-xb;  acr=y/(1.+xc*y)

        return a, b, c

    def corrImage(self, data, a, b, c, slope=10000):

        y = a * data - b
        result = (y / (1 + y * c)) * slope

        return result

    def MeanDEM(self, pointUL, pointDR):
        '''
        计算影像所在区域的平均高程.
        '''
        script_path = os.path.split(os.path.realpath(__file__))[0]
        dem_path = os.path.join(script_path,"GMTED2km.tif")

        try:
            DEMIDataSet = gdal.Open(dem_path)
        except Exception as e:
            pass

        DEMBand = DEMIDataSet.GetRasterBand(1)
        geotransform = DEMIDataSet.GetGeoTransform()
        # DEM分辨率
        pixelWidth = geotransform[1]
        pixelHight = geotransform[5]

        # DEM起始点：左上角，X：经度，Y：纬度
        originX = geotransform[0]
        originY = geotransform[3]

        # 研究区左上角在DEM矩阵中的位置
        yoffset1 = int((originY - pointUL['lat']) / pixelWidth)
        xoffset1 = int((pointUL['lon'] - originX) / (-pixelHight))

        # 研究区右下角在DEM矩阵中的位置
        yoffset2 = int((originY - pointDR['lat']) / pixelWidth)
        xoffset2 = int((pointDR['lon'] - originX) / (-pixelHight))

        # 研究区矩阵行列数
        xx = xoffset2 - xoffset1
        yy = yoffset2 - yoffset1
        # 读取研究区内的数据，并计算高程
        DEMRasterData = DEMBand.ReadAsArray(xoffset1, yoffset1, xx, yy)

        MeanAltitude = np.mean(DEMRasterData)
        return MeanAltitude


def CalDES(nowdate) :
    ''' 日地距离的计算 '''
    doy = nowdate.timetuple().tm_yday
    diy = 365

    X = 2 * np.pi * (doy-1) / diy
    temp = 1.000109 + 0.033494 * np.cos(X) + 0.001472 * np.sin(X)\
         + 0.000768 * np.cos(2*X) + 0.000079 * np.sin(2*X)

    return 1.0/temp

def radiance2reflectance(rad, Esun, sunz, nowdate):
    '''
    将radiance转换为reflectance
    :param rad: 卫星载荷通道入瞳处等效辐射亮度, W/m2/um
    :param Esun: 太阳辐照度, W/m2/um
    :param sunz: 太阳天顶角, degree
    :param nowdate: datetime, 当前日期
    :return:
    '''

    des = CalDES(nowdate)      # 日地距离，天体单位

    return (np.pi * rad * des**2) / (Esun * np.cos(sunz/180*np.pi))

def reflectance2radiance(ref, Esun, sunz, nowdate):
    '''
    将radiance转换为reflectance
    :param ref: 卫星载荷通道入瞳处等效辐射亮度
    :param Esun: 太阳辐照度
    :param sunz: 太阳天顶角, degree
    :param nowdate: datetime, 当前日期
    :return:
    '''

    des = CalDES(nowdate)      # 日地距离，天体单位

    return  (ref * Esun * np.cos((90-sunz)/180*np.pi)) / (np.pi * des**2)


