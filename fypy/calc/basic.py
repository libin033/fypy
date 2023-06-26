# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : basic.py

@Modify Time :  2022/10/31 15:52   

@Author      : Lee    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime

Deg2Rad = np.pi / 180.0
Rad2Deg = 180.0 / np.pi

def calc_earth_distance(lat1, lon1, lat2, lon2) :
    '''  计算球面两点之间距离（KM） '''
    lon11 = lon1 * Deg2Rad
    lat11 = lat1 * Deg2Rad
    lon21 = lon2 * Deg2Rad
    lat21 = lat2 * Deg2Rad

    dlon = lon21 - lon11
    dlat = lat21 - lat11
    h = np.sin(dlat/2)**2 + np.cos(lat11) * np.cos(lat21) * np.sin(dlon/2)**2
    distance = 2 * 6371.009 * np.arcsin(np.sqrt( h )) # 地球平均半径，6371km

    return distance

def planck_t2r(bt, wn):
    '''
    普朗克函数：将亮温（K）转成辐射值

    Parameters
    ----------
    bt : numpy.narray
        bright temperature, units: K
    wn : float or numpy.narray
        wave number(cm^-1)

    Returns
    -------
        numpy.narray
        卫星观测辐射值，mW/(m2.cm-1.sr)

    '''
    Radiation_C1 = 1.191042953E-5
    Radiation_C2 = 1.4387774
    a = Radiation_C1 * wn * wn * wn
    b = (Radiation_C2 * wn) / bt
    c = np.exp(b) - 1.0
    radi = a / c

    return radi


def planck_r2t(rad, wn):
    '''
    普朗克函数：将辐射值转成亮温（K）

    Parameters
    ----------
    rad : numpy.narray
        mW/(m2.cm-1.sr)
    wn : float or numpy.narray
        wave number(cm^-1)

    Returns
    -------
        numpy.narray
        bright temperature
        units: K
    '''

    Radiation_C1 = 1.191042953E-5
    Radiation_C2 = 1.4387774
    vs = wn
    bt = (Radiation_C2*vs/np.log(Radiation_C1*vs*vs*vs/(rad)+1.0))

    return bt


def k_index(pressure, temperature, dewpoint):
    '''
    Calculate K Index from the pressure temperature and dewpoint.

    K Index formula derived from [George1960]:
        K = (T850 - T500) + Td850 - (T700 - Td700)

    where:
    T850 is the temperature at 850 hPa
    T700 is the temperature at 700 hPa
    T500 is the temperature at 500 hPa
    Td850 is the dewpoint at 850 hPa
    Td700 is the dewpoint at 700 hPa

    Calculation of the K Index is defined as the temperature
    difference between the static instability between 850 hPa and 500 hPa,
     add with the moisture at 850hPa, then subtract from the dryness of the airmass at 700 hPa.

    Parameters
    ----------
    pressure
    temperature
    dewpoint

    Returns
    -------

    '''

def total_totals_index(pressure, temperature, dewpoint):
    '''
    Calculate Total Totals Index from the pressure temperature and dewpoint.

    Total Totals Index formula derived from [Miller1972]:
        TT = (T850 + Td850) - (2 * T500)

    where:

    T850 is the temperature at 850 hPa
    T500 is the temperature at 500 hPa
    Td850 is the dewpoint at 850 hPa

    Calculation of the Total Totals Index is defined as
    the temperature at 850 hPa plus the dewpoint at 850 hPa,
    minus twice the temperature at 500 hPa.
    This index consists of two components,
    the Vertical Totals (VT) and the Cross Totals (CT).

    Parameters
    ----------
    pressure
    temperature
    dewpoint

    Returns
    -------

    '''

def lifted_index(pressure, temperature, parcel_profile):
    '''
    Calculate Lifted Index from the pressure temperature and parcel profile.

    Lifted index formula derived from [Galway1956] and referenced by [DoswellSchultz2006]:
        LI = T500 - Tp500

    where:
    T500 is the measured temperature at 500 hPa
     Tp500 is the temperature of the lifted parcel at 500 hPa

    Calculation of the lifted index is defined as the temperature
     difference between the observed 500 hPa temperature and the
     temperature of a parcel lifted from the surface to 500 hPa.

    Parameters
    ----------
    pressure
    temperature
    parcel_profile

    Returns
    -------

    '''


def wind_direction(u, v, convention='from'):
    '''
    Compute the wind direction from u and v-components.

    Parameters
    ----------
    u
    v
    convention

    Returns
    -------

    '''


def calc_ndvi(vis006, nir008):

    return (nir008 - vis006) / (nir008 + vis006)

def calc_evi(vis004, vis006, nir008):

    G = 2.5     # 增益系数
    c1 = 6.0    # 气溶胶系数
    c2 = 7.5    # 气溶胶系数
    L = 1.0     # 土壤背景的调节系数

    return G * (nir008 - vis006) / (nir008 + c1 * vis006 - c2 * vis004 + L)

def calc_tvi(vis005, vis006, nir008):

    return 0.5 * (120 * (nir008 - vis005) - 200 * (vis006 - vis005))

def calc_rvi(vis006, nir008):

    return nir008 / vis006



