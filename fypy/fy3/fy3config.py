# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : fy3config.py

@Modify Time :  2022/10/27 16:13   

@Author      : fypy Team    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime

# FY3 10度块 编码对应关系
FY3Block10CoefX = {
    "00":0.0,
    "10":10.0,
    "20":20.0,
    "30":30.0,
    "40":40.0,
    "50":50.0,
    "60":60.0,
    "70":70.0,
    "80":80.0,
    "90":90.0,
    "A0":100.0,
    "B0":110.0,
    "C0":120.0,
    "D0":130.0,
    "E0":140.0,
    "F0":150.0,
    "G0":160.0,
    "H0":170.0,
    "I0":-10.0,
    "J0":-20.0,
    "K0":-30.0,
    "L0":-40.0,
    "M0":-50.0,
    "N0":-60.0,
    "O0":-70.0,
    "P0":-80.0,
    "Q0":-90.0,
    "R0":-100.0,
    "S0":-110.0,
    "T0":-120.0,
    "U0":-130.0,
    "V0":-140.0,
    "W0":-150.0,
    "X0":-160.0,
    "Y0":-170.0,
    "Z0":-180.0,
}

FY3Block10CoefY = {
    "80":  90.0,
    "70":  80.0,
    "60":  70.0,
    "50":  60.0,
    "40":  50.0,
    "30":  40.0,
    "20":  30.0,
    "10":  20.0,
    "00":  10.0,
    "90":   0.0,
    "A0": -10.0,
    "B0": -20.0,
    "C0": -30.0,
    "D0": -40.0,
    "E0": -50.0,
    "F0": -60.0,
    "G0": -70.0,
    "H0": -80.0,
}


FY3ProdInfo = {
    'FY3D' :{
        'MERSI' : {
            'CLA': {
                'sdsname': ['Global High Cloud Amount Day', 'Global High Cloud Amount Night'],
            },
            'CLM': {
                'sdsname': ['CLM_DAILY_D', 'CLM_DAILY_N',
                            'MERSI_NDVI_D', 'MERSI_NDVI_N'],
            },
            'GFR': {
                'sdsname': ['FIRES'],
            },
            'LST': {
                'sdsname': ['MERSI_25Km_LST_D', 'MERSI_25Km_LST_N',
                            'MERSI_1Km_LST_D', 'MERSI_1Km_LST_N',
                            'MERSI_obt_LST_D', 'MERSI_obt_LST_N'],
            },
            'NDVI': {
                'sdsname': ['5KM_NDVI', '250m NDVI'],
            },
            'NVI': {
                'sdsname': ['5KM_NDVI', '250m NDVI'],
            },
            'PWV': {
                'sdsname': ['MERSI_PWV'],
            },
            'SST': {
                'sdsname': ['sea_surface_temperature'],
            },
            'SIC': {
                'sdsname': ['Daily_Sea_Ice_Both', 'Both'],
            },
            'TPW': {
                'sdsname': ['MERSI_TPW', 'MERSI_DAY_TPWSDS', 'MERSI_NIGHT_TPWSDS'],
            },
        },
        'MWRI'  : {
            'CLW': {
                'sdsname': [''],
            },
            'MRR': {
                'sdsname': ['RainRate'],
            },
            'SWE': {
                'sdsname': ['SWE_Northern_Daily', 'SWE_Southern_Daily'],
            },
            'SWS': {
                'sdsname': ['SWS_ORBIT', 'SWS_Ascending', 'SWS_Descending'],
            },
        },
        'MWHS'  : {
            'RDT': {
                'sdsname': ['Rain Detection'],
            },
        },
        'TSHS' : {
            'AVP': {
                'sdsname': ['/DATA/TSHS_AT_Prof', '/DATA/TSHS_AH_Prof',
                            '/GEO/Latitude','/GEO/Longitude'],
            },
        },
        'GNOS' : {
            'ATP' : {
                'sdsname': ['Temp'],
            }
        },
    },
}

Prj_Info = {
    'GLL' : 'GEOGCS["WGS84",'
                'DATUM["WGS_1984",'
                    'SPHEROID["WGS84",6378137,298.257223563,'
                        'AUTHORITY["EPSG","7030"]],'
                    'TOWGS84[0,0,0,0,0,0,0],'
                'AUTHORITY["EPSG","6326"]],'
                'PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],'
                'UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9108"]],'
                'AUTHORITY["EPSG","4326"]]',

    'HAM' : 'PROJCRS["World_Hammer_Aitoff",' \
                'BASEGEOGCRS["WGS 84", ' \
                    'DATUM["World Geodetic System 1984",' \
                        'ELLIPSOID["WGS 84",6378137,298.257223563,LENGTHUNIT["metre",1]]],' \
                'PRIMEM["Greenwich",0,ANGLEUNIT["Degree",0.0174532925199433]]],' \
                'CONVERSION["World_Hammer_Aitoff",METHOD["Hammer_Aitoff"],' \
                'PARAMETER["False_Easting",0,LENGTHUNIT["metre",1]],' \
            'PARAMETER["False_Northing",0,LENGTHUNIT["metre",1]],' \
            'PARAMETER["Central_Meridian",0,ANGLEUNIT["Degree",0.0174532925199433]]],' \
            'CS[Cartesian,2],AXIS["(E)",east,ORDER[1],LENGTHUNIT["metre",1]],' \
            'AXIS["(N)",north,ORDER[2],LENGTHUNIT["metre",1]],' \
            'USAGE[SCOPE["Not known."],AREA["World."],BBOX[-90,-180,90,180]],' \
            'ID["ESRI",54044]]',

}