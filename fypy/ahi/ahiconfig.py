# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : ahiconfig.py

@Modify Time :  2024/1/30 9:30   

@Author      : fypy Team    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime


AreaInfo = {
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
