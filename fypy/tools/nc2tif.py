# -*- coding:utf-8 -*-
'''
@Project  : fypy

@File     : nc2tif.py

@Modify Time : 2022/8/11 15:34

@Author      : fypy Team

@Version : 1.0

@Description :
针对netCDF4数据，将数据集进行投影转换成GeoTiff文件

'''
import os
import sys
import numpy as np

import netCDF4
from fypy.tools import writetiff


def nc2tif(outname, filename, sdsname, fillvalue=-999):
    '''
    将NC文件中的数据集根据lat、lon维度进行投影输出geotiff文件

    Parameters
    ----------
    outname: str
        输出tif文件名
    filename: str
        输入的nc文件名
    sdsname: str
        需要转换的数据集名
    fillvalue: float or int
        数据填充值
    Returns
    -------
        bool
        成功返回TRUE，失败返回FALSE
    '''

    try:
        fp = netCDF4.Dataset(filename, 'r')

        dsetid = fp[sdsname]
        dims = dsetid.dimensions

        grid = False
        xdim = None
        ydim = None
        for item in dims :
            if 'lon' in item.lower() :
                grid = True
                xdim = item

            if 'lat' in item.lower() :
                grid = True
                ydim = item

        sdsinfo = dsetid.ncattrs()
        if xdim is not None and ydim is not None :
            x = fp[xdim][:]
            y = fp[ydim][:]
            xres = x[1] - x[0]
            yres = y[1] - y[0]

            trans = [x[0], xres, 0, y[0], 0, yres]
            data = fp[sdsname][:]
            if '_FillValue' in sdsinfo :
                fillvalue = sdsinfo['_FillValue']
            writetiff(outname, data, trans=trans, fillvalue=fillvalue)
        fp.close()

        return True
    except BaseException as e :
        print(e)
        return False

    # 后续更新根据coordinate属性来完成投影设置

    for item in sdsinfo :
        if 'coordinate' in item.lower() :
            grid = True



