# -*- coding:utf-8 -*-
'''
@Project     : ANUSPLIN

@File        : readStationFile.py

@Modify Time :  2022/9/27 13:50   

@Author      : fypy Team    

@Version     : 1.0   

@Description :

'''
import os
import sys
import numpy as np
import datetime



def readStationsFile(filename):
    dict_data = {}
    if not os.path.isfile(filename):
        print("%s is not exist, will be return!!!" %filename)
        return dict_data

    temp = []
    datasize = 0

    fin = open(filename, 'r')#, encoding='utf-8')
    data = fin.readlines()
    for item in data[1:]:
        item = item.strip('\n')
        item_split = item.split(' ')
        item_split = ['NULL' if i == '' else i for i in item_split]

        # 记录总共有多少个字段
        if datasize == 0 :
            datasize = len(item_split)

        if len(item_split) != datasize  or datasize == 0:
            continue
        temp.append(item_split)
        #print(item.split())

    temp = np.array(temp)
    print(temp.shape)
    if temp.shape[0] == 0:
        print('Not DownLoad Valid Data, will return...')
        return dict_data


    for i in range(temp.shape[1]):
        key = temp[0, i]
        temp_data = temp[1:, i]
        # if key in ['Station_Name', 'Province', 'City', 'Cnty'] :
        #     temp_data = temp[1:, i]
        # elif "Datetime" in key:
        #     temp_data = np.array(temp_data, dtype=np.int64) # dtype= '<S20')
        # elif "Station_Id_C" in key:
        #     continue
        #     # temp_data = np.array(temp_data, dtype=np.int64) # dtype= '<S20')
        # else:
        #     temp_data = np.array(temp_data, dtype=np.float32)

        if not key in dict_data.keys():
            dict_data.update({ key : temp_data})
        else:
            dict_data[key].extend(temp_data)

    return dict_data

def AnalysisData(txtname, hdfname = None):
    '''
    解析气象观测数据，并输出HDF文件，每个气象观测要素为一个数据集
    :param txtname: 气象数据下载文件
    :param hdfname: 输出HDF文件
    :return:
    '''
    dictdata = readStationsFile(txtname)

    flag = 1
    for key in dictdata.keys():

        # if 'Station_Name' in key:
        #     continue
        # elif "Datetime" in key:
        #     WriteHDF(hdfname, key, np.array(dictdata[key], dtype= '<S20'), overwrite=flag)
        # else:           ### 今后在此做容错方案，对写入HDF数据进行检查，如果有字符串则转换成string型写入HDF文件
        #     WriteHDF(hdfname, key, np.array(dictdata[key], dtype=np.float32), overwrite=flag)
        flag = 0

    return dictdata

    # if 'Station_levl' in dictdata.keys():
    #     # Station_level = np.array(dictdata['Station_levl'], dtype=np.int32)
    #     Station_level = dictdata['Station_levl']
    #     ind = np.where(Station_level <= 12)
    #

    # 开始按站点ID提取气象观测数据
    # for item in ind[0]:
    #     stationID = dictdata['Station_Id_d'][item]
    #     stationname = r'D:\ProductSystem\MeteDownLoad\data\station\%05d.HDF' %(int(stationID))      # 拼接站点文件名
    #     print(stationname)
    #     # if os.path.isfile(stationname):
    #     with h5py.File(stationname, 'a') as fp:
    #         for key in SDS:
    #             if key in dictdata.keys():
    #                 print(stationname, key)
    #                 data = dictdata[key][item]
    #                 if key in fp:
    #                     dset = fp[key]
    #                     print(dset.shape)
    #                     dset.resize(dset.shape[0] + 1)
    #                     dset[dset.shape[0]:dset.shape[0] + data.shape[0]] = data
    #                 else:
    #                     fp.create_dataset(key, data=data)



