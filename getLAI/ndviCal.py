import sys
import gdal
from readtiff import read_tiff
import writeoneRaster
""" 计算NDVI

@ndviCal: 计算NDVI值并保存
tifLoc: 原始影像数据位置
redBandLoc: 红波段位置
nirBandLoc: 近红波段位置
ndviLoc: 生成的NDVI影像保存位置
"""
# redBandLoc 2 / nirBandLoc 3
def ndviCal(tifLoc, redBandLoc, nirBandLoc, ndviLoc):
    
    _, geoTransform1, proj1, bandCount, col, row, _ = read_tiff(tifLoc)
    ds = gdal.Open(tifLoc)
    
    # 如果 redBandLoc 和 nirBandLoc 超过了bandCount 那么退出程序
    # GetRasterBand 是从1开始计数的
    if redBandLoc < 1 or nirBandLoc < 1 or redBandLoc > bandCount or nirBandLoc > bandCount:
        print('Something wrong happens')
        sys.exit()
    
    # 读取红波段
    redBand = ds.GetRasterBand(redBandLoc)
    red = redBand.ReadAsArray(0, 0, row, col)
    
    # 读取近红波段
    nirBand = ds.GetRasterBand(nirBandLoc)
    nir = nirBand.ReadAsArray(0, 0, row, col)
    
    # 计算NDVI
    ndvi = (nir - red) / (nir + red)
    
    # 保存
    writeoneRaster.writeOneRaster(ndviLoc, ndvi, geoTransform1, proj1, col, row)