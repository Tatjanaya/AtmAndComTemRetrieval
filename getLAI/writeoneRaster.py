import gdal
import numpy as np
import pandas as pd

def writeOneRaster(outpath, data, geoTransform, proj, row, col):
    # 判断类型
    if 'int8' in data.dtype.name:
        datatype = gdal.GDT_Byte
    elif 'int16' in data.dtype.name:
        datatype = gdal.GDT_UInt16
    else:
        datatype = gdal.GDT_Float32
    # 判断维数
    if len(data.shape) == 3:
        bands, im_height, im_width = data.shape
    elif len(data.shape) == 2:
        im_data = np.array([data])
    else:
        im_bands, (im_height, im_width) = 1, data.shape
    driver = gdal.GetDriverByName("GTiff")
    outRaster = driver.Create(outpath, col, row, 1, datatype)
    outRaster.SetGeoTransform(geoTransform)
    outRaster.SetProjection(proj)  # 直接用现有的proj
    outRaster.GetRasterBand(1).WriteArray(data) #写入数组数据
    del outRaster