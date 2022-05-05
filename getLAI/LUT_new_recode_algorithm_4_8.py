import numpy as np
import pandas as pd
import datetime
from numba import jit
#from searchclose import find_close
from readtiff import read_tiff
import writeoneRaster

starttime = datetime.datetime.now()

@jit(nopython=True)
def find_close(arr , e):
    low = 0
    high = len(arr) - 1
    idx = -1

    while low <= high:
        mid = int((low + high) / 2)
        if e == arr[mid] or mid == low:
            idx = mid
            break
        elif e > arr[mid]:
            low = mid
        elif e < arr[mid]:
            high = mid

    if idx + 1 < len(arr) and abs(e - arr[idx]) > abs(e - arr[idx+1]):
        idx += 1

    return idx

@jit(nopython=True)
def functionA(row, col, data1, array1, m_range, Result):
    for i in range(row):
        for j in range(col):
            if data1[i, j, 0] < 0 or data1[i, j, 0] == 65536:
                data1[i, j, :] = np.nan
            elif (data1[i, j, 2] + data1[i, j, 1]) != 0 and (data1[i, j, 2] - data1[i, j, 1]) / (data1[i, j, 2] + data1[i, j, 1]) < 0.35:
                data1[i, j, :] = np.nan
            else:
                # p1 = find_close(array1[:, 3], data1[i, j, 2] / 10000 - m_range)
                # p2 = find_close(array1[:, 3], data1[i, j, 2] / 10000 + m_range)
                # p = np.arange(p1, p2 + 1, 1)
                temp1 = 10000000
                t = 0
                for k in range(array1.shape[0]):  # 这里修正了以前的一个小bug，即可能出现array[-1,3]的情况，不过np默认array[-1,3]为倒数第一行，不影响结果
                    temp = np.power((data1[i, j, 0] / 10000 - array1[k, 1 + 0]), 2) \
                        + np.power((data1[i, j, 1] / 10000 - array1[k, 1 + 1]), 2) \
                            + np.power((data1[i, j, 2] / 10000 - array1[k, 1 + 2]), 2)
                    if temp < temp1:
                        temp1 = temp
                        t = k
                if array1[t, 0] < 5.0:
                    Result[i, j] = round(array1[t, 0] / 0.7, 2)
                else:
                    Result[i, j] = array1[t, 0]

#@jit(nogil=True, parallel=True)
data, geoTransform1, proj1, band, row, col, data1 = read_tiff(r'res.tif')
Result = np.zeros([row, col])
array1 = np.array(pd.read_excel(r'./res/LAI.xlsx', header=None))
Ban = 3  # 用户要反演的波段数
g_range = array1[:, 0].max() - array1[:, 0].min()
r_range = array1[:, 1].max() - array1[:, 1].min()
m_range = np.sqrt(np.power(g_range, 2)+np.power(r_range, 2))

functionA(row, col, data1, array1, m_range, Result)    

writeoneRaster.writeOneRaster("tcLAI.tif", Result, geoTransform1, proj1, row, col)
endtime = datetime.datetime.now()

print((endtime - starttime).seconds)