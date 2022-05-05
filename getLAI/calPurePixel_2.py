import os
import sys
import gdal
import numpy as np
""" 从特征空间中计算纯净像元-纯土壤像元和纯植被像元
@calPurePixel: 计算植被和土壤反射率并保存
resTif: 原始影像
resTifRedBand: 红波段在原始影像顺序序号
resTifNirBand: 近红波段在原始影像顺序序号
ndviTif: 同一区域NDVI影像
bandLst: 波段列表（预处理后的原始影像可能不存在波段信息，需要额外对应）
txtDir: 生成的植被和土壤txt保存文件夹位置
"""
def calPurePixel(resTif, redBandCount, nirBandCount, bandLst):
    # 读取预处理后原始影像数据
    ds = gdal.Open(resTif)
    
    row = ds.RasterXSize
    col = ds.RasterYSize
    band = ds.RasterCount
    
    # 如果波段列表和原始影像波段个数不一致 不再运行
    if len(bandLst) != band:
        print('Something wrong happens')
        sys.exit()
    
    # 读取ndvi影像数据
    # ds_ndvi = gdal.Open(ndviTif)
    
    # row_ndvi = ds_ndvi.RasterXSize
    # col_ndvi = ds_ndvi.RasterYSize
    
    # 如果ndvi和res不等大 则退出程序
    # if row != row_ndvi or col != col_ndvi:
    #     print("Something wrong happens")
    #     sys.exit()
    
    # 读取ndvi
    # ndvi = ds_ndvi.ReadAsArray(0, 0, row, col)
    
    # 读取红波段和近红波段数据
    redBand = ds.GetRasterBand(redBandCount)
    red = redBand.ReadAsArray(0, 0, row, col)
    
    nirBand = ds.GetRasterBand(nirBandCount)
    nir = nirBand.ReadAsArray(0, 0, row, col)
    
    # 将red值和nir值放入一维列表中
    redFlatten = red.flatten()
    nirFlatten = nir.flatten()
    
    # 排序
    redFlatten.sort()
    nirFlatten.sort()
    
    # 获取初始土壤点集
    # 划分256个区间 取每个区间近红外最小值作为初始土壤点
    # 初始土壤点存储列表
    soilPointLst = []
    
    interval = len(redFlatten) / 256
    for i in range(256):
        print('xxxx')
        # 找到每个区间的值
        firLoc = i * interval
        secLoc = firLoc + interval
        nirVal = 10000
        loc_x = -1
        loc_y = -1
        for j in range(int(firLoc), int(secLoc)):
            # 找到位置
            redVal = redFlatten[j]
            redValLoc = np.where(red == redVal)
            nir_x = redValLoc[0][0]
            nir_y = redValLoc[1][0]
            if nir[nir_x][nir_y] < nirVal:
                nirVal = nir[nir_x][nir_y]
                loc_x = nir_x
                loc_y = nir_y
        soilPointLst.append([loc_x, loc_y])
    
    # 计算0-50% 0-75% 0-100% 25-75% 25-100%子集的最小二乘相关系数
    # 0-50%
    lst1 = soilPointLst[0: int(len(soilPointLst) / 2)]
    lst1_r = lstSqCorrCoe(lst1, red, nir)
    print('ok')
    # 0-75%
    lst2 = soilPointLst[0: int(len(soilPointLst) * 0.75)]
    lst2_r = lstSqCorrCoe(lst2, red, nir)
    print('ok')
    # 0-100%
    lst3 = soilPointLst[:]
    lst3_r = lstSqCorrCoe(lst3, red, nir)
    print('ok')
    # 25-75%
    lst4 = soilPointLst[int(len(soilPointLst) * 0.25) : int(len(soilPointLst) * 0.75)]
    lst4_r = lstSqCorrCoe(lst4, red, nir)
    print('ok')
    # 25-100%
    lst5 = soilPointLst[int(len(soilPointLst) * 0.25):]
    lst5_r = lstSqCorrCoe(lst5, red, nir)
    print('ok')
    # 相关系数最大的子集确定为土壤点
    lst_r = [lst1_r, lst2_r, lst3_r, lst4_r, lst5_r]
    lst_t = [lst1, lst2, lst3, lst4, lst5]
    lst_res = lst_t[lst_r.index(max(lst_r))]
    
    # 取这些点平均值作为土壤结果
    soilLst = []
    for i in range(1, band + 1):
        iBand = ds.GetRasterBand(i)
        ival = iBand.ReadAsArray(0, 0, row, col)
        nowBandValTotal = 0
        for val in lst_res:
            nowBandValTotal += ival[val[0]][val[1]]
        nowBandValTotal /= len(lst_res)
        soilLst.append(nowBandValTotal)
        print(nowBandValTotal)
    
    # 植被 ndvi 前10个 nir 前10个的重合点平均值 如果没有 则预先设置内置参数
    vegLst = []
    # ndviFlatten = ndvi.flatten()
    # ndviFlatten.sort()
    # ndviMaxLst = []
    
    # 10个
    ndviNumPer2 = 10
    nirNumPer2 = 10
    # for i in range(10):
        
    #     ndviMaxLst.append(ndviFlatten[-(i + 1)])
    # # 找到ndvi最大值所在位置
    # maxLocLst = []
    # for val in ndviMaxLst:
    #     maxLocLst.append(np.where(ndvi == val))
    # for i in range(1, band + 1):
    #     iBand = ds.GetRasterBand(i)
    #     ival = iBand.ReadAsArray(0, 0, row, col)
    #     nowBandValTotal = 0
    #     for val in maxLocLst:
    #         nowBandValTotal += ival[val[0][0]][val[1][0]]
    #     nowBandValTotal /= len(maxLocLst)
    #     vegLst.append(nowBandValTotal)
    
    # 写入txt 如果没有该文件夹则创建
    if not os.path.exists('./res/'):
        os.makedirs('./res/')
        
    # f = open('./res/leaf.txt', 'w')
    # for i in range(len(bandLst)):
    #     f.write(str(bandLst[i]) + ' ' + str(format(vegLst[i] / 10000, '.4f')) + ' ' + str(format(vegLst[i] / 10000, '.4f')))
    #     if i < len(bandLst) - 1:
    #         f.write('\n')
    # f.close()
    
    f = open('./res/soil.txt', 'w')
    for i in range(len(bandLst)):
        f.write(str(bandLst[i]) + ' ' + str(format(soilLst[i] / 10000, '.4f')))
        if i < len(bandLst) - 1:
            f.write('\n')
    f.close()

# 计算最小二乘相关系数
def lstSqCorrCoe(lst, redTif, nirTif):
    redValLst = []
    nirValLst = []
    # print(lst)
    for temp in lst:
        redValLst.append(redTif[temp[0]][temp[1]])
        nirValLst.append(nirTif[temp[0]][temp[1]])
    # 平均值
    redAvg = sum(redValLst) / len(redValLst)
    nirAvg = sum(nirValLst) / len(nirValLst)
    upVal = 0
    downLeftVal = 0
    downRightVal = 0
    for i in range(len(redValLst)):
        upVal += (redValLst[i] - redAvg) * (nirValLst[i] - nirAvg)
        downLeftVal += (redValLst[i] - redAvg) ** 2
        downRightVal += (nirValLst[i] - nirAvg) ** 2
    
    r = upVal / np.sqrt(downLeftVal) / np.sqrt(downRightVal)
    return r
    
calPurePixel('E:\integrate\getLAI\data\mp_saihanba.tif', 2, 3, [560, 665, 865])