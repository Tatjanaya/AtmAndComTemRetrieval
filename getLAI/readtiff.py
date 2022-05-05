import gdal
import numpy as np

def read_tiff(inpath):
    ds = gdal.Open(inpath)  # 一般的遥感数据打开方式
    row = ds.RasterYSize  # 获取行数
    col = ds.RasterXSize  # 获取列数
    band = ds.RasterCount  # 获取波段数
    geoTransform = ds.GetGeoTransform()  # 仿射矩阵，共有六个参数，分表代表左上角x坐标；东西方向上图像的分辨率；如果北边朝上，地图的旋转角度，0表示图像的行与x轴平行；左上角y坐标；如果北边朝上，地图的旋转角度，0表示图像的列与y轴平行；南北方向上地图的分辨率
    proj = ds.GetProjection()  # 地图投影信息
    data = ds.ReadAsArray(0, 0, col, row)  # 用ReadAsArray(<xoff>, <yoff>, <xsize>, <ysize>)，读出从(xoff,yoff)开始，大小为(xsize,ysize)的矩阵
    data1 = np.zeros([row, col, band])
    for i in range(band):
        data1[:, :, i] = ds.GetRasterBand(i+1).ReadAsArray(0, 0, col, row)
    del ds
    return data, geoTransform, proj, band, row, col, data1