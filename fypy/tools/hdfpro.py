# -*- coding:utf-8 -*-
'''
@Project  : fypy

@File     : hdfpro.py

@Modify Time : 2022/8/11 15:34

@Author      : fypy Team

@Version : 1.0

@Description :
对hdf5文件进行读写操作

'''
import os
import h5py

def readhdf(filename, sdsname, dictsdsinfo=None) :
    r'''
    读取hdf数据集和属性
    :param filename:
    :param sdsname:
    :param dictsdsinfo:
    :return:
    '''
    if not os.path.isfile(filename) :
        print('%s is not exist, will return None' %(filename))
        return None
    else:
        try:
            fp = h5py.File(filename, 'r')
            dsetid = fp[sdsname]
    
            data = dsetid[:]
            if not dictsdsinfo is None :
                for key in dsetid.attrs :
                    # print(key)
                    dictsdsinfo.update({key : dsetid.attrs[key]})
                fp.close()
                return data, dictsdsinfo
    
            fp.close()
        except BaseException as e :
            print(e)
            return None

        return data

def readhdf_fileinfo(filename) :
    '''
    读取hdf数据集和属性
    :param filename:
    :param sdsname:
    :param dictsdsinfo:
    :param dictfileinfo:
    :return:
    '''

    dictfileinfo = {}
    if not os.path.isfile(filename) :
        print('%s is not exist, will return None' %(filename))
        return dictfileinfo
    else:
        fp = h5py.File(filename, 'r')
        for key in fp.attrs :
            # print(key)
            dictfileinfo.update({key : fp.attrs[key]})
        fp.close()
        return dictfileinfo

def writehdf(filename, sdsname, data, overwrite=True,
             dictsdsinfo=None, dictfileinfo=None,
             compression=9, info=False):
    '''
    创建hdf5文件
    :param filename:
    :param sdsname:
    :param data:
    :param overwrite:
    :param dictsdsinfo:
    :param dictfileinfo:
    :param compression:
    :return:
    '''

    try:

        if overwrite :
            fp = h5py.File(filename, 'w')
        else:
            fp = h5py.File(filename, 'r+')
        if not dictfileinfo is None :
            for key in dictfileinfo :
                fp.attrs[key] = dictfileinfo[key]

        dsetid = fp.create_dataset(sdsname, data=data, compression = compression)

        if not dictsdsinfo is None :
            for key in dictsdsinfo :
                dsetid.attrs[key] = dictsdsinfo[key]

        fp.close()
        if info:
            print('create %s %s success...' %(filename, sdsname))
    except BaseException as e :
        print(e)
        return False

    return True


def writehdf_fileinfo(filename, dictfileinfo, overwrite=True, info=False):
    '''
    创建文件属性
    :param filename:
    :param dictfileinfo:
    :param overwrite:
    :param info: 是否打印信息
    :return:
    '''

    try:

        if overwrite :
            fp = h5py.File(filename, 'w')
        else:
            fp = h5py.File(filename, 'r+')

        if isinstance(dictfileinfo, dict):
            for key in dictfileinfo :
                fp.attrs[key] = dictfileinfo[key]
                if info :
                    print('create %s attr [%s] success...' %(
                        filename, key
                    ))


    except BaseException as e :
        print(e)
        return False

    return True

