# -*- coding:utf-8 -*-
'''
@Project  : fypy

@File     : hdf4pro.py

@Modify Time : 2022/8/11 15:34

@Author      : fypy Team

@Version : 1.0

@Description :
对hdf4文件进行读写

'''
import os

def readhdf4(h4file, sdsname, dictsdsattrs=None, dictfileattrs=None):
    ''' 读取hdf4数据集 '''
    data = None

    from pyhdf import SD
    if not os.path.isfile(h4file):
        print('%s is not exist, will continue...' %(h4file))
        return data

    fp4 = SD.SD(h4file, SD.SDC.READ)

    if dictfileattrs is not None :
        dictfileattrs = {}
        attrs = fp4.attributes(full=1)

        for item in fp4.attributes().keys():
            dictfileattrs[item] = attrs[item][0]


    if sdsname in fp4.datasets().keys() :
        try:
            sds4id = fp4.select(sdsname)
            data = sds4id[:]

            if dictsdsattrs is not None :
                dictsdsattrs = {}
                attrs = sds4id.attributes(full=1)
                for key in sds4id.attributes() :
                    dictsdsattrs[key] =  attrs[key][0]

        except BaseException as e:
            print(e)
    fp4.end()

    if dictfileattrs is not None and dictsdsattrs is not None :
        return data, dictfileattrs, dictsdsattrs
    elif dictfileattrs is not None and dictsdsattrs is None :
        return data, dictfileattrs
    elif dictfileattrs is None and dictsdsattrs is not None :
        return data, dictsdsattrs
    else:
        return data


def readhdf4sdsattrs(h4file, sdsname):
    ''' 读取hdf4文件数据集属性信息'''
    from pyhdf import SD

    sdsattrs = {}
    if not os.path.isfile(h4file):
        print('%s is not exist, will continue...' %(h4file))
        return sdsattrs
    try:
        fp4 = SD.SD(h4file, SD.SDC.READ)
        sds4id = fp4.select(sdsname)
        attrs = sds4id.attributes(full=1)
        for key in sds4id.attributes() :
            sdsattrs[key] =  attrs[key][0]
    except BaseException as e:
        print(e)

    return sdsattrs


def readhdf4fileattrs(h4file):
    ''' 读取hdf4文件属性'''
    from pyhdf import SD

    fileattrs = {}

    if not os.path.isfile(h4file):
        print('%s is not exist, will continue...' %(h4file))
        return fileattrs

    fp4 = SD.SD(h4file, SD.SDC.READ)

    attrs = fp4.attributes(full=1)

    for item in fp4.attributes().keys():
        try:
            fileattrs[item] = attrs[item][0]
        except BaseException as e:
            print(e)

    return fileattrs