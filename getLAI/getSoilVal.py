import gdal
import numpy as np
from scipy import optimize
""" 获取土壤线
@getSoilVal: 获取土壤线
originTif: 目标影像
"""
def getSoilVal(originTif):
    # 读数据
    ds = gdal.Open(originTif)
    col = ds.RasterXSize
    row = ds.RasterYSize
    data = ds.ReadAsArray(0, 0, col, row)
    green = data[0, 0: row, 0: col]
    red = data[1, 0: row, 0: col]
    nir = data[2, 0: row, 0: col]
    
    # 分为256个组
    redFlatten = red.flatten()
    redFlatten.sort()
    # 去掉为0的点
    temp = 0
    for i in range(len(redFlatten)):
        if redFlatten[i] != 0:
            temp = i
            break
    redFlatten = redFlatten[temp: len(redFlatten)]
    # 取这些组最小近红外的1%的点
    calLst = []
    du = int(len(redFlatten) / 256)
    count = 0
    # 由于求取相关系数发现0-50%是理想集合 因此只取前一半的点
    for i in range(128):
        left = i * du
        right = left + du
        if right >= len(redFlatten):
            right = len(redFlatten) - 1
        redSlice = redFlatten[left: right + 1]
        # 去重
        redLst = list(set(redSlice))
        # nirLst (nir, [][])
        nirLst = []
        # redSlice遍历 将对应位置的nir值找出 然后选择1%的点加入calLst calLst为位置
        for red_val in redLst:
            tempLst = np.where(red == red_val)
            for j in range(len(tempLst[0])):
                nirLst.append((nir[tempLst[0][j]][tempLst[1][j]], (tempLst[0][j], tempLst[1][j])))
        # nirLst 排序
        nirLst.sort()
        # 选取前1% 的列入 calLst
        one_per = int(len(nirLst) * 0.01)
        for j in range(one_per):
            calLst.append(nirLst[j])
        
        # 只选第一个
        # calLst.append(nirLst[0])
        count += 1
        print(count)
    
    # 做拟合直线 循环28次 剔除掉离直线最远的点 剩下100个作为最终结果取平均
    x_r = []
    y_nir = []
    z_green = []
    for val in calLst:
        x = val[1][0]
        y = val[1][1]
        x_r.append(red[x][y])
        y_nir.append(nir[x][y])
        z_green.append(green[x][y])
    
    targetLen = int(len(x_r) / 0.9)
    print(targetLen)
    
    while len(x_r) > targetLen:
        print(len(x_r))
        # 直线拟合与绘制
        A, B = optimize.curve_fit(f_1, x_r, y_nir)[0]
        # 计算列表中的点与直线距离 距离最大删除
        indexLoc = -1
        max_val = -1
        for i in range(len(x_r)):
            if calDis(A, B, x_r, y_nir, i) > max_val:
                max_val = calDis(A, B, x_r, y_nir, i)
                indexLoc = i
        
        del x_r[indexLoc]
        del y_nir[indexLoc]
        del z_green[indexLoc]
        
    # 平均
    nir_res = sum(y_nir) / len(y_nir)
    red_res = sum(x_r) / len(x_r)
    green_res = sum(z_green) / len(z_green)
    
    print('green: ' + str(green_res))
    print('res: ' + str(red_res))
    print('nir: ' + str(nir_res))
    
# 直线方程函数
def f_1(x, A, B):
    return A * x + B

# A B 是直线方程系数
# xLst yLst 是x y值列表
# i是列表中的索引
def calDis(A, B, xLst, yLst, i):
    return np.abs(A * xLst[i] - yLst[i] + B) / np.sqrt(A * A + 1)
    
getSoilVal('res.tif')