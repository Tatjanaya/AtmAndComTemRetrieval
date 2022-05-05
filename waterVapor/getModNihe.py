import os
import numpy as np
import gdal
import xlrd
import matplotlib.pyplot as plt
import seaborn as sns 
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score

'''
获取 Sentinel-3 s8 Nadir t / s9 Nadir t 透过率比值和水汽含量之间的关系 线性
tranDataFile: chn透过率计算结果文件
waterVaporFile: 大气廓线 水汽对照表
'''
def getModNihe(tranDataFile, waterVaporFile):
    # 如果不存在文件夹 则创建
    dirs = './res/'
    if not os.path.exists(dirs):
        os.makedirs(dirs)
        
    content = xlrd.open_workbook(waterVaporFile)
    sh = content.sheet_by_name('page_1')
    # 水汽廓线列表
    waterVapor = []
    for i in range(1, sh.nrows):
        waterVapor.append(float(sh.cell(i, 1).value))
    # 打开chn文件 读取两个波段透过率并且求比值 s9 / s8
    transRatio = []
    f = open(tranDataFile, 'r',encoding='utf-8')
    chnData = f.readlines()
    f.close()
    for row in range(946):
        transRatio.append(float(1.0 - float(chnData[row * 7 + 6].split()[2])) / float(1.0 - float(chnData[row * 7 + 5].split()[2])))
    
    # 设置中文画图
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 黑体
    plt.rcParams['axes.unicode_minus'] = False    # 解决无法显示符号的问题
    # 坐标轴的刻度设置向内(in)或向外(out)
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    sns.set(style="white", font='SimHei')
    # 进行拟合
    para, _ = curve_fit(Fun, transRatio, waterVapor)
    
    plt.title('水汽含量与S9 S8波段透过率比值拟合图')
    plt.xlabel('τ9/τ8')
    plt.ylabel('Water vapor content(g cm-2)')
    plt.grid()
    plt.scatter(transRatio, waterVapor, s=1, c='#df03fc', marker='o')
    xx = np.linspace(0, max(transRatio), 100)
    yy = Fun(xx, para[0], para[1])
    plt.xlim(xmax=1.0, xmin=0.4)
    plt.ylim(ymax=8.0, ymin=0.0)
    plt.plot(xx, yy, '-g')
    plt.legend(["线性拟合结果", "MODTRAN结果散点"])
    plt.savefig('./res/水汽和透过率比值拟合结果.png')
    # 把拟合参数和RMSE结果打印到txt中
    predictRes = para[0] + para[1] * np.array(transRatio)
    mse = mean_squared_error(waterVapor, predictRes)
    rmse = np.sqrt(mse)
    r2 = r2_score(waterVapor, predictRes)
    
    f = open("./res/方差协方差比值法.txt", "a")
    f.write('参数 -> ' + '形式 y = ' + str(format(para[0], '.4f')) + ' + ' + str(format(para[1], '.4f')) + ' * x ')
    f.write('\n')
    f.write('r2 ' + str(format(r2, '.4f')))
    f.write('\n')
    f.write('rmse ' + str(format(rmse, '.4f')))
    
def Fun(x, a, b):
    return a + b * x

# 得到近似的角度 整10的 ->
def getApproxVza(vzaTif):
    img_data, row, col, _ = tifPara(vzaTif)
    res = np.zeros((row, col))
    for x in range(row):
        for y in range(col):
            res[x][y] = getApproxVal(float(img_data[x][y]))
    return res
            
def getApproxVal(num):
    tenVal = int(num) / 10
    sinVal = num % 10
    if sinVal >= 5:
        return (tenVal + 1) * 10
    else:
        return tenVal * 10

# 打开影像文件
def tifPara(tifFile):
    ds = gdal.Open(tifFile)
    col = ds.RasterXSize
    row = ds.RasterYSize
    rasBand = ds.GetRasterBand(1)
    img = ds.ReadAsArray(0, 0, col, row)
    return img, row, col, rasBand.DataType

getModNihe('../../modtran/tape5/trans_10.chn', './data/廓线-水汽对照表.xlsx')