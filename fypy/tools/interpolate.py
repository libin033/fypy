# -*- coding:utf-8 -*-
'''
@Project     : fypy

@File        : interpolate.py

@Modify Time :  2022/12/21 16:44   

@Author      : fypy Team    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime

# import cv2

def griddata_simple(srclat, srclon, srcdata, gridy, gridx):
    lats = np.array(srclat).flatten()
    lons = np.array(srclon).flatten()
    data = np.array(srcdata).flatten()

    xlim = [gridx[0], gridx[-1]]
    ylim = [gridy[0], gridy[-1]]
    # id = np.where((lons >= xlim[0]) & (lons <= xlim[1]) & (lats >= ylim[0]) & (lats <= ylim[1]))
    # data = data[id]
    # lons = lons[id]
    # lats = lats[id]

    dx = gridx[1] - gridx[0]
    dy = gridy[1] - gridy[0]

    lons_around = np.floor((lons - xlim[0]) / dx)
    lats_around = np.floor((lats - ylim[0]) / dy)
    nn = (lats_around * gridx.size + lons_around).reshape(-1)

    nn_sort_index = np.argsort(nn)
    nn = nn[nn_sort_index].astype(np.int)
    nn_index, nn_unique_index = np.unique(nn, True)
    nn_unique_index = nn_unique_index[1:] - 1
    nn_unique_index = np.append(nn_unique_index, nn.shape[0] - 1)

    re = np.zeros(gridx.size*gridy.size)
    means_t = data.reshape(-1)
    means_t = means_t[nn_sort_index]
    means_t = np.cumsum(means_t)
    means_t = means_t[nn_unique_index]

    re[nn[nn_unique_index[1:]]] = (means_t[1:] - means_t[0:-1]) / (nn_unique_index[1:] - nn_unique_index[0:-1])
    re[nn[nn_unique_index[0]]] = means_t[0] / (nn_unique_index[1])
    re = re.reshape([gridy.size, gridx.size])
    # if sn > 0:
    #
    #     for i in range(0, sn):
    #         re[np.isnan(re)] = 0
    #         re_ = cv2.filter2D(re, -1, np.ones((3, 3)))
    #         re01 = np.zeros(re_.shape)
    #         re01[re > 0] = 1
    #         re01_ = cv2.filter2D(re01, -1, np.ones((3, 3)))
    #         ren = re_ / re01_
    #
    #         re[re01 == 0] = ren[re01 == 0]

    return re

