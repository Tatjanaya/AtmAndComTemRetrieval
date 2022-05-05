""" 通过查找表获取结果

@functionA: 查找出最符合的LAI值

@lutCalLAI: 输出tif影像
"""
import numpy as np
import pandas as pd
import datetime
from numba import jit
from readtiff import read_tiff
import writeoneRaster

@jit(nopython=True)
def functionA(row, col, data1, array1, result, x, y, z):
    for i in range(row):
        for j in range(col):
            if data1[i, j, x] < 0 or data1[i, j, x] == 65536:
                data1[i, j, :] = np.nan
            elif (data1[i, j, z] + data1[i, j, y]) != 0 and (data1[i, j, z] - data1[i, j, y]) / (data1[i, j, z] + data1[i, j, y]) < 0.35:
                data1[i, j, :] = np.nan
            else:
                temp1 = 10000000
                t = 0
                for k in range(array1.shape[0]):  # 这里修正了以前的一个小bug，即可能出现array[-1,3]的情况，不过np默认array[-1,3]为倒数第一行，不影响结果
                    temp = np.power((data1[i, j, x] / 10000 - array1[k, 1 + 0]), 2) \
                        + np.power((data1[i, j, y] / 10000 - array1[k, 1 + 1]), 2) \
                            + np.power((data1[i, j, z] / 10000 - array1[k, 1 + 2]), 2)
                    if temp < temp1:
                        temp1 = temp
                        t = k
                if array1[t, 0] < 5.0:
                    result[i, j] = round(array1[t, 0] / 0.7, 2)
                else:
                    result[i, j] = array1[t, 0]


def lutCalLAI(picTif, LAIExcelFile, resTif, x, y, z):
    starttime = datetime.datetime.now()
    _, geoTransform1, proj1, _, row, col, data1 = read_tiff(picTif)
    result = np.zeros([row, col])
    array1 = np.array(pd.read_excel(LAIExcelFile, header=None))

    functionA(row, col, data1, array1, result, x, y, z)    

    writeoneRaster.writeOneRaster(resTif, result, geoTransform1, proj1, row, col)
    endtime = datetime.datetime.now()
    print((endtime - starttime).seconds)