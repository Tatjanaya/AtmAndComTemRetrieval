import numpy as np
import gdal

def upscale(originTif, cankaoTif, saveLoc):
    img, row, col, _, _, _ = tifPara(originTif)
    _, row_cankao, col_cankao, rasDataType, geoTransform, proj = tifPara(cankaoTif)
    res = np.zeros((row_cankao, col_cankao))
    row_interval = int(row / row_cankao)
    col_interval = int(col / col_cankao)
    for i in range(row_cankao):
        for j in range(col_cankao):
            if i == row_cankao - 1 and j == col_cankao - 1:
                res[i][j] = calAvgVal(img, i * row_interval, row - 1, j * col_interval, col - 1)
            elif i == row_cankao - 1:
                res[i][j] = calAvgVal(img, i * row_interval, row - 1, j * col_interval, j * col_interval + col_interval)
            elif j == col_cankao - 1:
                res[i][j] = calAvgVal(img, i * row_interval, i * row_interval + row_interval, j * col_interval, col - 1)
            else: 
                res[i][j] = calAvgVal(img, i * row_interval, i * row_interval + row_interval, j * col_interval, j * col_interval + col_interval)
    saveTif(res, rasDataType, saveLoc, geoTransform, proj)
                
def calAvgVal(data, xminRan, xmaxRan, yminRan, ymaxRan):
    res = 0
    count = (xmaxRan - xminRan) * (ymaxRan - yminRan)
    for x in range(xminRan, xmaxRan):
        for y in range(yminRan, ymaxRan):
            # 如果有异常点
            if data[x][y] < 0.5:
                count -= 1
                continue
            res += data[x][y]
    res /= count
    return res

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
    
upscale('./res/s8_trans.tif', '../res/s8_0_tran_init.tif', './res/upscale_tran.tif')