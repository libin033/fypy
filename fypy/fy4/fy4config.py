# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : fy3config.py

@Modify Time :  2022/10/27 16:37   

@Author      : fypy Team    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime

from fypy import parm

exedir = os.path.abspath(list(parm.__path__)[0])

FontTTF = os.path.join(exedir, 'font', 'simsun.ttf')
if not os.path.isfile(FontTTF) :
    raise Exception('字体文件不存在【%s】' %(FontTTF))

FY4ProdInfo = {
    'FY4A':{
        'AGRI' : {
            'CIX': {
                'sdsname': ['Convective_Initiation'],
            },
            'CLM': {
                'sdsname': ['CLM'],
            },
            'CLP': {
                'sdsname': ['CLP'],
            },
            'CLT': {
                'sdsname': ['CLT'],
            },
            'CTH': {
                'sdsname': ['CTH'],
            },
            'CTP': {
                'sdsname': ['CTP'],
            },
            'CTT': {
                'sdsname': ['CTT'],
            },
            'DSD': {
                'sdsname': ['DST'],
            },
            'FHS': {
                'sdsname': ['FPA'],
            },
            'FOG': {
                'sdsname': ['FOG'],
            },
            'LST': {
                'sdsname': ['LST'],
            },
            'SST': {
                'sdsname': ['SST'],
            },
            'QPE': {
                'sdsname': ['Precipitation'],
            },
            'SNC': {
                'sdsname': ['SNC'],
            },
        },
        'LMI' : {
            'LMIE': {
                'sdsname': [''],
            },
            'LMIG': {
                'sdsname': [''],
            },
        },
        'GIIRS' : {
            'AVP': {
                'sdsname': ['AT_Prof'],
            },
        }
    }
}

AreaInfo = {
    'fy4': {
        'agri' : {
            '0.005' : {
                'description' : 'FY-4 full disk area definition at 500m resolution',
                'shape': [21984, 21984],
            },
            '0.01' : {
                'description' : 'FY-4 full disk area definition at 1km resolution',
                'shape': [10992, 10992],
            },
            '0.02' : {
                'description' : 'FY-4 full disk area definition at 2km resolution',
                'shape': [5496, 5496],
            },
            '0.04' : {
                'description' : 'FY-4 full disk area definition at 4km resolution',
                'shape': [2748, 2748],
            }
        }
    },
    'ahi': {
        'ahi' : {
            '0.005' : {
                'description' : 'Himawari-8/9 full disk area definition at 500m resolution',
                'shape': [22000, 22000],
                'extent':
                    [-5499999.9684, -5499999.9684, 5499999.9684, 5499999.9684],        # minX, minY, maxX, maxY
            },
            '0.01' : {
                'description' : 'Himawari-8/9 full disk area definition at 1km resolution',
                'shape': [11000, 11000],
                'extent':
                    [-5500000.0355, -5500000.0355, 5500000.0355, 5500000.0355],        # minX, minY, maxX, maxY
            },
            '0.02' : {
                'description' : 'Himawari-8/9 full disk area definition at 1km resolution',
                'shape': [5500, 5500],
                'extent':
                    [-5499999.9012, -5499999.9012, 5499999.9012, 5499999.9012],        # minX, minY, maxX, maxY
            }
        }
    }
}



