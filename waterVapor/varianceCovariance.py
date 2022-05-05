import numpy as np
import os
import gdal

'''
求最后结果
paraLst: 参数 形式为a + b * x 的 a 和 b
s8nBtToa s9nBtToa: 两波段Nadir角度TOA亮温
'''
def getViaNihe(paraLst, s8nBtToa, s9nBtToa):
    # 如果不存在文件夹 则创建
    dirs = './res/'
    if not os.path.exists(dirs):
        os.makedirs(dirs)
    # 读取文件
    s8_data, row, col, ras_datatype, geoTransform, proj = tifPara(s8nBtToa)
    s9_data, _, _, _, _, _ = tifPara(s9nBtToa)
    # 计算Rji
    res_rji = np.zeros((row, col))
    # n为5
    for x in range(row):
        for y in range(col):
            if x < 2 or y < 2 or x >= row - 2 or y >= col - 2:
                continue
            res_rji[x][y] = calRji(s8_data, s9_data, x, y)
    # 周边两条用隔壁补
    for x in range(row):
        for y in range(col):
            if x < 2 and y < 2:
                res_rji[x][y] = res_rji[2][2]
            elif x < 2 and y >= col - 2:
                res_rji[x][y] = res_rji[2][col - 3]
            elif x >= row - 2 and y < 2:
                res_rji[x][y] = res_rji[row - 3][2]
            elif x >= row - 2 and y >= col - 2:
                res_rji[x][y] = res_rji[row - 3][col - 3]
            elif x < 2 and y >= 2 and y < col - 2:
                res_rji[x][y] = res_rji[2][y]
            elif x >= row - 2 and y >= 2 and y < col - 2:
                res_rji[x][y] = res_rji[row - 3][y]
            elif y < 2 and x >= 2 and x < row - 2:
                res_rji[x][y] = res_rji[x][2]
            elif y >= col - 2 and x >= 2 and x < row - 2:
                res_rji[x][y] = res_rji[x][col - 3]
            else:
                continue
    # 得到水汽
    wv_cal = paraLst[0] + paraLst[1] * res_rji
    # 如果有负数的 划为0
    wv_cal = np.where(wv_cal > 0, wv_cal, 0)
    
    # 再生成透过率和上行辐射
    # -0.005x*x+-0.083x+0.981	0.060x*x+0.640x+-0.064
    s8_trans = -0.005 * wv_cal * wv_cal - 0.083 * wv_cal + 0.981
    s8_lup = 0.060 * wv_cal * wv_cal + 0.640 * wv_cal - 0.064
    
    # 保存
    saveTif(wv_cal, ras_datatype, './res/wv_cal.tif', geoTransform, proj)
    saveTif(s8_trans, ras_datatype, './res/s8_trans.tif', geoTransform, proj)
    saveTif(s8_lup, ras_datatype, './res/s8_lup.tif', geoTransform, proj)
    
            
# 计算每一个像元位置处的Rji
def calRji(s8_data, s9_data, x, y):
    # 计算平均值
    s8_avg = 0
    s9_avg = 0
    for i in range(x - 2, x + 3):
        for j in range(y - 2, y + 3):
            s8_avg += float(s8_data[i][j])
            s9_avg += float(s9_data[i][j])
    s8_avg /= 25
    s9_avg /= 25
    # 计算分子和分母
    # 分子 numerator
    # 分母 denominator
    numerator = 0
    denominator = 0
    for i in range(x - 2, x + 3):
        for j in range(y - 2, y + 3):
            numerator += (s8_data[i][j] - s8_avg) * (s9_data[i][j] - s9_avg)
            denominator += (s8_data[i][j] - s8_avg) ** 2
    return numerator / denominator

# 打开影像文件
def tifPara(tifFile):
    ds = gdal.Open(tifFile)
    col = ds.RasterXSize
    row = ds.RasterYSize
    geoTransform = ds.GetGeoTransform()
    proj = ds.GetProjection()
    rasBand = ds.GetRasterBand(1)
    img = ds.ReadAsArray(0, 0, col, row)
    return img, row, col, rasBand.DataType, geoTransform, proj

def saveTif(nArray, nInfo, targetLoc, geoTransform, proj):
    driver = gdal.GetDriverByName('GTiff')
    op = driver.Create(targetLoc, nArray.shape[1], nArray.shape[0], 1, nInfo)
    op.SetGeoTransform(geoTransform)
    op.SetProjection(proj)
    op.GetRasterBand(1).WriteArray(nArray)

paraLst = [14.7806, -14.6851]
s8nBtToa = './data/s8nbt.tif'
s9nBtToa = './data/s9nbt.tif'
getViaNihe(paraLst, s8nBtToa, s9nBtToa)