# -*- coding:utf-8 -*-
'''
@Project  : fypy

@File     : h4toh5.py

@Modify Time : 2022/8/11 15:34

@Author      : fypy Team

@Version : 1.0

@Description :
将hdf4格式文件转换成hdf5文件

'''
import os
import sys
import h5py



def h4toh5(h4file, h5file):
    '''
    将hdf4格式文件转换成hdf5文件

    Parameters
    ----------
    h4file: str
        输入hdf4文件名
    h5file: str
        输出hdf5文件名
    Returns
    -------
        bool
        成功返回True, 失败返回False
    '''
    from pyhdf import SD

    if not os.path.isfile(h4file):
        print('%s is not exist, will continue...' %(h4file))
        return False

    fp4 = SD.SD(h4file, SD.SDC.READ)
    fp5 = h5py.File(h5file, 'w')
    attrs = fp4.attributes(full=1)

    # 将hdf4中的文件属性写入hdf5
    for item in fp4.attributes().keys():
        # print(item)
        try:
            fp5.attrs[item] = attrs[item][0]
        except BaseException as e:
            print(e)

    # 将hdf4中的数据集和属性写入hdf5
    for sdsname in fp4.datasets().keys() :
        # print(sdsname)
        try:
            sds4id = fp4.select(sdsname)
            data = sds4id[:]
            dsetid = fp5.create_dataset(name= sdsname, data=data, compression=9)
            attrs = sds4id.attributes(full=1)
            for key in sds4id.attributes() :
                dsetid.attrs[key] =  attrs[key][0]
        except BaseException as e:
            print(e)
    fp4.end()
    fp5.close()

    return True











