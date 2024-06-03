# -*- coding:utf-8 -*-
'''
@Project  : fypy

@File     : ncpro.py

@Modify Time : 2022/8/11 15:34

@Author      : fypy Team

@Version : 1.0

@Description :
Netcdf4文件的读写操作

'''
import os
import netCDF4
import numpy as np
from netCDF4 import num2date, date2num

def readnc(filename, sdsname, dictsdsinfo=None):
    ''' 读取netcdf4 文件数据集和数据集属性 '''
    if not os.path.isfile(filename) :
        print('%s is not exist, please check it...' %(filename))
        if dictsdsinfo is None :
            return None
        else:
            return None, None
    else:
        try:
            with netCDF4.Dataset(filename, 'r') as fp :
                if not sdsname in fp.variables:
                    print('%s is not a dataset in %s' % (sdsname, filename))

                    if dictsdsinfo is None :
                        return None
                    else:
                        return None, None
                else:
                    dsetid = fp.variables[sdsname]
                    data = fp.variables[sdsname][:]

                    if not dictsdsinfo is None :
                        for key in dsetid.ncattrs() :
                            dictsdsinfo.update({key: dsetid.getncattr(key)})


                        return data, dictsdsinfo

        except BaseException as e :
            print(e)
            if dictsdsinfo is None :
                return None
            else:
                return None, None

        return data

def readnc_fileinfo(filename) :
    ''' 读取netcdf4 文件全局属性 '''
    dictfileinfo = {}
    if not os.path.isfile(filename) :
        print('%s is not exist, please check it...' %(filename))
        return dictfileinfo
    else:
        try:
            with netCDF4.Dataset(filename, 'r') as fp :

                for key in fp.ncattrs():
                    dictfileinfo.update({key: fp.getncattr(key)})

                return  dictfileinfo
        except BaseException as e :
            print(e)
            return dictfileinfo

def readnc_sdsinfo(filename, sdsname) :
    ''' 读取数据集属性信息 '''
    dictsdsinfo = {}
    if not os.path.isfile(filename) :
        print('%s is not exist, please check it...' %(filename))
        return dictsdsinfo
    else:
        try:
            with netCDF4.Dataset(filename, 'r') as fp :

                if not sdsname in fp.variables:
                    print('%s is not a dataset in %s' %(sdsname, filename))
                    return dictsdsinfo
                else:
                    dsetid = fp.variables[sdsname]
                    for key in dsetid.ncattrs():
                        dictsdsinfo.update({key: dsetid.getncattr(key)})

                    return dictsdsinfo
        except BaseException as e :
            print(e)
            return dictsdsinfo


def readncfortimes(filename, sdsname,
                   units = 'hours since 1900-01-01 00:00:00.0',
                   calendar = "gregorian"):
    '''
    读取时间维度数据集

    Parameters
    ----------
    filename : str
        文件名
    sdsname : str
        数据集名
    units : str
        时间单位，从什么时间开始计数
    calendar : str

    Returns
    -------
        numpy.array
        时间维度数据集
    '''

    if not os.path.isfile(filename) :
        print('%s is not exist, please check it...' %(filename))
        return None
    else:
        try:
            with netCDF4.Dataset(filename, 'r') as fp :
                if not sdsname in fp.variables:
                    print('%s is not a dataset in %s' % (sdsname, filename))
                    return None
                else:
                    dsetid = fp.variables[sdsname]
                    data = fp.variables[sdsname][:]

                    if 'units' in dsetid.ncattrs() :
                        units = dsetid.getncattr('units')

                    if 'calendar' in dsetid.ncattrs() :
                        calendar = dsetid.getncattr('calendar')

                    data = num2date(data,units=units,calendar=calendar)

                    return data

        except BaseException as e :
            print(e)
            return None

        return data


def writenc(filename, sdsname, srcdata, dimension=None, overwrite=1,
            complevel=9, dictsdsinfo=None, fill_value=None,
            standard_name=None, long_name=None, description=None, units=None,
            valid_range=None,
            scale_factor=None, add_offset=None, **kwargs):
    '''
    写NC文件

    Parameters
    ----------
    filename : str
    sdsname : str
    srcdata : str
    dimension : tuple
        数据集关联的维度名
    overwrite : int
    complevel : int
        压缩等级，0-9，值越大，压力等级越大
    dictsdsinfo ：dict
        数据集属性
    fill_value  : int or float
        填充值 或 无效值
    standard_name
    long_name
    description
    units
    valid_range
    scale_factor
    add_offset
    kwargs

    Returns
    -------

    '''

    data = np.array(srcdata)

    if (dictsdsinfo is None) or not isinstance(dictsdsinfo, dict) :
        dictsdsinfo={}

    dictsdsinfo = checkinfo(dictsdsinfo, 'standard_name', standard_name)
    dictsdsinfo = checkinfo(dictsdsinfo, 'standard_name', sdsname)
    dictsdsinfo = checkinfo(dictsdsinfo, 'long_name', long_name)
    dictsdsinfo = checkinfo(dictsdsinfo, 'long_name', sdsname)
    dictsdsinfo = checkinfo(dictsdsinfo, 'valid_range', valid_range)
    dictsdsinfo = checkinfo(dictsdsinfo, 'units', units)
    dictsdsinfo = checkinfo(dictsdsinfo, 'scale_factor', scale_factor)
    dictsdsinfo = checkinfo(dictsdsinfo, 'add_offset', add_offset)
    dictsdsinfo = checkinfo(dictsdsinfo, 'description', description)

    for key in kwargs :
        if not key in dictsdsinfo :
            dictsdsinfo[key] = kwargs[key]

    if overwrite :
        fp = netCDF4.Dataset(filename, 'w')
    else:
        fp = netCDF4.Dataset(filename, 'r+')

    dsetid = None
    # 如果没有对应的关联维度，则创建该维度
    if dimension is None :
        fp.createDimension(sdsname, data.shape[0])
        dsetid = fp.createVariable(sdsname, data.dtype, dimensions=sdsname,
                                  zlib=True, complevel=9)
        fp.variables[sdsname][:] = data
    else:
        # 判断数据集维度是否存在，如果不存在，则创建该维度，并赋值为0
        k = 0
        for item in dimension :
            if not item in fp.dimensions :
                fp.createDimension(item, data.shape[k])
                dsetid = fp.createVariable(item, np.int32, dimensions=item, zlib=True, complevel=9)
                fp.variables[item][:] = range(data.shape[k])
            k = k + 1
        dsetid = fp.createVariable(sdsname, data.dtype, dimensions=dimension, zlib=True,
                                   complevel=complevel, fill_value=fill_value)
        fp.variables[sdsname][:] = data

    if (not dictsdsinfo is None) and (not dsetid is None) :
        for key in dictsdsinfo :
            dsetid.setncattr(key, dictsdsinfo[key])

    fp.close()
    # print('create %s %s success...' %(filename, sdsname))


def writenc_fileinfo(filename, dictfileinfo, overwrite=1):

    try:
        if overwrite :
            fp = netCDF4.Dataset(filename, 'w')
        else:
            if os.path.isfile(filename) :
                fp = netCDF4.Dataset(filename, 'r+')
            else:
                print('%s is not exist, will create it...' %(filename))
                fp = netCDF4.Dataset(filename, 'w')

        if not dictfileinfo is None :
            for key in dictfileinfo:
                fp.setncattr(key, dictfileinfo[key])
        fp.close()
        return True

    except BaseException as e :
        print('write %s file information error!!!' %(filename))
        return False


def writenc_sdsinfo(filename, sdsname, dictsdsinfo, overwrite=1):
    try:
        if overwrite:
            fp = netCDF4.Dataset(filename, 'w')
        else:
            if os.path.isfile(filename):
                fp = netCDF4.Dataset(filename, 'r+')
            else:
                print('%s is not exist, will create it...' % (filename))
                fp = netCDF4.Dataset(filename, 'w')
        if not sdsname in fp.variables:
            print('%s is not a dataset in %s' %(sdsname, filename))
            fp.close()
            return False
        else:
            dsetid = fp.variables[sdsname]
            for key in dictsdsinfo:
                dsetid.setncattr(key, dictsdsinfo[key])
        fp.close()
        return True
    except BaseException as e:
        print('write %s dataset %s information error!!!' % (filename, sdsname))
        return False

def writencfortimes(filename, sdsname, srcdata, overwrite=1,
                    units = 'hours since 1900-01-01 00:00:00.0',
                    calendar = "gregorian",
                    complevel=9, dictsdsinfo=None):
    '''
    写时间维度数据集

    Parameters
    ----------
    filename : str
    sdsname: str
        datetime array
    srcdata : datetime array

    overwrite: int
    units : str
        时间单位，从什么时间开始计数
    calendar : str
    complevel : int
    dictsdsinfo : dict
        数据集属性

    Returns
    -------

    '''

    data = np.array(srcdata)

    if overwrite :
        fp = netCDF4.Dataset(filename, 'w')
    else:
        fp = netCDF4.Dataset(filename, 'r+')

    # 如果没有对应的关联维度，则创建该维度

    fp.createDimension(sdsname, data.shape[0])
    times = fp.createVariable(sdsname, 'f8', dimensions=sdsname,
                              zlib=True, complevel=complevel)
    times.units = units
    times.calendar = calendar

    fp.variables[sdsname][:] = date2num(data,units=times.units,calendar=times.calendar)


    if (not dictsdsinfo is None) and (not times is None) :
        for key in dictsdsinfo :
            times.setncattr(key, dictsdsinfo[key])

    fp.close()
    print('create %s %s success...' %(filename, sdsname))


def checkinfo(dictinfo, key, value):

    if not key in dictinfo :
        if value is not None :
            dictinfo[key] = value

    return dictinfo

