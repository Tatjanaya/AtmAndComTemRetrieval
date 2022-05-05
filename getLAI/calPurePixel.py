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
def calPurePixel(resTif, ndviTif, bandLst):
    # 读取预处理后原始影像数据
    ds = gdal.Open(resTif)
    
    row = ds.RasterXSize
    col = ds.RasterYSize
    band = ds.RasterCount
    
    # 读取ndvi影像数据
    ds_ndvi = gdal.Open(ndviTif)
    
    row_ndvi = ds_ndvi.RasterXSize
    col_ndvi = ds_ndvi.RasterYSize
    
    # 如果ndvi和res不等大 则退出程序
    if row != row_ndvi or col != col_ndvi:
        print("Something wrong happens")
        sys.exit()
    
    # 读取ndvi
    ndvi = ds_ndvi.ReadAsArray(0, 0, row, col)
    
    # 找到区域内ndvi的最大值和最小值
    ndvi_min = ndvi.min()
    ndvi_max = ndvi.max()
    
    # 像元二分法模型基础上计算植被覆盖度
    vfc = (ndvi - ndvi_min) / (ndvi_max - ndvi_min)
    
    # 认为植被覆盖度最高的为该区域纯植被像元 植被覆盖度最低的为该区域纯土壤像元
    # 防止异常像元代入 使用累计百分比 取1% 和 99% 分别为纯土壤像元和纯植被像元
    flaVfc = vfc.flatten()
    valLoc = int(len(flaVfc) * 0.01)
    flaVfc.sort()
    maxVal = flaVfc[-valLoc]
    minVal = flaVfc[valLoc]
    
    locMax = np.where(vfc == maxVal)
    locMin = np.where(vfc == minVal)

    max_x = locMax[0][0]
    max_y = locMax[1][0]
    
    min_x = locMin[0][0]
    min_y = locMin[1][0]
    
    # 如果波段列表和原始影像波段个数不一致 不再运行
    if len(bandLst) != band:
        print('Something wrong happens')
        sys.exit()
    
    # 获得纯植被像元反射率和纯土壤像元反射率
    vegLst = []
    soilLst = []
    for i in range(1, band + 1):
        iBand = ds.GetRasterBand(i)
        ival = iBand.ReadAsArray(0, 0, row, col)
        vegLst.append(ival[max_x][max_y])
        soilLst.append(ival[min_x][min_y])
        
    # 写入txt 如果没有该文件夹则创建
    if not os.path.exists('./res/'):
        os.makedirs('./res/')
        
    f = open('./res/leaf.txt', 'w')
    for i in range(len(bandLst)):
        f.write(str(bandLst[i]) + ' ' + str(vegLst[i] / 10000) + ' ' + str(vegLst[i] / 10000))
        if i < len(bandLst) - 1:
            f.write('\n')
    f.close()
    
    f = open('./res/soil.txt', 'w')
    for i in range(len(bandLst)):
        f.write(str(bandLst[i]) + ' ' + str(soilLst[i] / 10000))
        if i < len(bandLst) - 1:
            f.write('\n')
    f.close()

# calPurePixel('res.tif', 'ndvi.tif', [554, 659, 868])