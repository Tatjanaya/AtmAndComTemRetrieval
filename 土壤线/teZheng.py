from cProfile import label
import gdal
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 
import getSoilValNew

# 设置中文画图
plt.rcParams['font.sans-serif'] = ['SimHei']  # 中文字体设置-黑体
plt.rcParams['axes.unicode_minus'] = False  # 解决保存图像是负号'-'显示为方块的问题
sns.set(style="whitegrid", font='SimHei')

def tifPara(tifFile):
    ds = gdal.Open(tifFile)
    col = ds.RasterXSize
    row = ds.RasterYSize
    rasBand = ds.GetRasterBand(1)
    img = ds.ReadAsArray(0, 0, col, row)
    geoTransform = ds.GetGeoTransform()
    proj = ds.GetProjection()
    return img, row, col, rasBand.DataType, geoTransform, proj

red, row, col, _, _, _ = tifPara('r.tif')
nir, _, _, _, _, _ = tifPara('n.tif')

red = red.astype(np.float)
nir = nir.astype(np.float)

# red /= 10000
# nir /= 10000
num = 10000

zhibei_lst_r = []
zhibei_lst_n = []
red_lst = []
nir_lst = []
for x in range(row):
    for y in range(col):
        if (nir[x][y] - red[x][y]) / (nir[x][y] + red[x][y]) > 0.68:
            zhibei_lst_r.append(red[x][y])
            zhibei_lst_n.append(nir[x][y])
        red_lst.append(red[x][y])
        nir_lst.append(nir[x][y])
soil_r_lst, soil_n_lst = getSoilValNew.getSoilVal(nir, red, row, col, 2000000, 88)

soil_r_lst = np.array(soil_r_lst)
soil_n_lst = np.array(soil_n_lst)
zhibei_lst_r = np.array(zhibei_lst_r)
zhibei_lst_n = np.array(zhibei_lst_n)
print(np.mean(soil_r_lst))
print(np.mean(soil_n_lst))
print(np.mean(zhibei_lst_r))
print(np.mean(zhibei_lst_n))

plt.xlabel('红波段反射率')
plt.ylabel('近红波段反射率')

plt.xlim(xmin=0, xmax=0.4)
plt.ylim(ymin=0, ymax=0.6)

plt.scatter(np.array(red_lst) / num, np.array(nir_lst) / num, s=1, c='dodgerblue', label='像元点')
plt.scatter(np.array(zhibei_lst_r) / num, np.array(zhibei_lst_n) / num, s=1, c='forestgreen', label='纯植被像元')
plt.scatter(np.array(soil_r_lst) / num, np.array(soil_n_lst) / num, s=1, c='darkred', label='纯土壤像元')
plt.legend(loc='best')
plt.show()