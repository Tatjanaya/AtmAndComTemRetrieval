from math import ceil, floor
import sys
import xlrd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 
import matplotlib.font_manager as fm
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
""" 从chn文件中提取有用信息 并储存在一个excel文件中

@drawLineChart: 读取存储数据的excel文件，并拟合任意两列数据
excelFile: 存储数据的excel文件位置
baseCol: 拟合的自变量 第几列
targetCol: 拟合的因变量 第几列 

@colName2Chinese: 将列名转换为中文名
axisName: 列名

@transNumber: 转换字符串为列号
col: 列号或列名
sh: xlrd读取的excel文件组织格式

@Fun: 拟合结构 设置为二次方程

@error: 拟合偏差

@getIndexes: 获取RMSE和R方
"""
def drawLineChart(excelFile, baseCol, targetCol, weidu, wv, duibiexcel):
    # 设置中文画图
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 中文字体设置-黑体
    plt.rcParams['axes.unicode_minus'] = False  # 解决保存图像是负号'-'显示为方块的问题
    sns.set(style="white", font='SimHei')
    
    # 如果不存在文件夹 则创建
    dirs = './res/'
    if not os.path.exists(dirs):
        os.makedirs(dirs)
    
    # 打开excel表
    xlsxFile = xlrd.open_workbook(excelFile)
    # 按Sheet定位工作表
    sh = xlsxFile.sheet_by_name('page_1')
    
    # 提取对比excel
    duibiFile = xlrd.open_workbook(duibiexcel)
    shduibi = duibiFile.sheet_by_name('page_1')
    bianhao_lst = []
    wv_lst = []
    for i in range(1, shduibi.nrows):
        bianhao_lst.append(int(shduibi.cell(i, 0).value))
        wv_lst.append(float(shduibi.cell(i, 1).value))
    
    # 转换列表名
    baseCol = transNumber(baseCol, sh)
    targetCol = transNumber(targetCol, sh)
    
    # 如果目标两列列号大于文件的列号 终止程序 列号是最后一列也终止
    if baseCol >= sh.ncols or targetCol >= sh.ncols:
        print("Something wrong happens")
        sys.exit()
        
    # 首先获取两列的名称 例如 trans_0_S8 ——> S8波段0°透过率
    xAxis = sh.cell(0, baseCol).value
    yAxis = sh.cell(0, targetCol).value
    
    # 两列中文名
    xAxisChinese = colName2Chinese(xAxis)
    yAxisChinese = colName2Chinese(yAxis)
    
    # 根据weidu判断取出范围
    if weidu == 'trop':
        minHao = 1
        maxHao = 872
    elif weidu == 'midsum':
        minHao = 873
        maxHao = 1260
    elif weidu == 'midwin':
        minHao = 1261
        maxHao = 1614
    elif weidu == 'polar':
        minHao = 1615
        maxHao = 2311
    else:
        minHao = 873
        maxHao = 1614
    
    # 两列数据取出
    xData = []
    yData = []
    for i in range(1, sh.nrows):
        if bianhao_lst[i - 1] >= minHao and bianhao_lst[i - 1] <= maxHao and wv_lst[i - 1] >= wv[0] and wv_lst[i - 1] <= wv[1]:
            xData.append(float(sh.cell(i, baseCol).value))
            yData.append(float(sh.cell(i, targetCol).value))
        
    # 拟合
    x1 = np.array(xData)
    y1 = np.array(yData)
    para, _ = curve_fit(Fun, x1, y1)
    
    plt.title("Sentinel3_SLSTR_" + xAxisChinese + "和" + yAxisChinese)
    # x y范围需要根据两列数据情况调整
    xMax = max(xData)
    yMax = max(yData)
    xMin = min(xData)
    yMin = min(yData)
    
    plt.xlim(xmax=int(ceil(xMax)), xmin=int(floor(xMin)))
    plt.ylim(ymax=int(ceil(yMax)), ymin=int(floor(yMin)))
    
    plt.xlabel(xAxisChinese)
    plt.ylabel(yAxisChinese)
    plt.grid()
    plt.scatter(xData, yData, s=1, c='#FF0000', marker='o')
    
    # 确定x 边界
    xxSec = 1
    if max(xMax, yMax) > 1:
        xxSec = max(xMax, yMax) + 1
        
    xx = np.linspace(0, xxSec, 100)
    yy = Fun(xx, para[0], para[1], para[2])
    plt.plot(xx, yy, '-g')
    plt.legend(["二次拟合曲线", "MODTRAN模拟散点"])
    
    predictRes = para[0] * x1 * x1 + para[1] * x1 + para[2] # predict
    realResult = y1 # real
    mse = mean_squared_error(realResult, predictRes)
    rmse = np.sqrt(mse)
    r2 = r2_score(realResult, predictRes)
    # 将结果写入文本文件
    f = open("./res/MODTRANResult-2.txt", "a")
    f.write("\n" + str(xAxis) + " " + str(yAxis) + "\n" + str(format(para[0], '.3f')) + 'x*x+' + str(format(para[1], '.3f')) + 'x+' + str(format(para[2], '.3f')) + "\n" \
        + " RMSE: " + str(format(rmse, '.3f')) + "\n")
    # 保存
    plt.savefig("./res/" + xAxis + '_' + yAxis + ".png")
    plt.cla()
    # plt.show()
    
def waterVaporDraw(excelFile, atmCol, lowWv, highWv):
    # 设置中文画图
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 中文字体设置-黑体
    plt.rcParams['axes.unicode_minus'] = False  # 解决保存图像是负号'-'显示为方块的问题
    sns.set(style="white", font='SimHei')
    
    # 如果不存在文件夹 则创建
    dirs = './res/'
    if not os.path.exists(dirs):
        os.makedirs(dirs)
    
    # 打开excel表
    xlsxFile = xlrd.open_workbook(excelFile)
    # 按Sheet定位工作表
    sh = xlsxFile.sheet_by_name('page_1')
    # 转换列表名
    atmCol = transNumber(atmCol, sh)
    
    # 找到大气参数一列
    yAxis = sh.cell(0, atmCol).value
    yAxisChinese = colName2Chinese(yAxis)
    yData = []
    # 找到水汽含量一列
    xAxis = sh.cell(0, sh.ncols - 1).value
    xAxisChinese = '水汽含量(g/cm2)'
    xData = []
    for i in range(1, sh.nrows):
        if float(sh.cell(i, sh.ncols - 1).value) >= lowWv and float(sh.cell(i, sh.ncols - 1).value) <= highWv:
            xData.append(float(sh.cell(i, sh.ncols - 1).value))
            yData.append(float(sh.cell(i, atmCol).value))
    
    # 拟合
    x1 = np.array(xData)
    y1 = np.array(yData)
    para, _ = curve_fit(Fun, x1, y1)
    
    plt.title("Sentinel3_SLSTR_" + xAxisChinese + "和" + yAxisChinese)
    # x y范围需要根据两列数据情况调整
    xMax = max(xData)
    yMax = max(yData)
    xMin = min(xData)
    yMin = min(yData)
    
    plt.xlim(xmax=int(ceil(xMax)), xmin=int(floor(xMin)))
    plt.ylim(ymax=int(ceil(yMax)), ymin=int(floor(yMin)))
    
    plt.xlabel(xAxisChinese)
    plt.ylabel(yAxisChinese)
    plt.grid()
    plt.scatter(xData, yData, s=1, c='#FF0000', marker='o')
    
    # 确定x 边界
    xxSec = 1
    if max(xMax, yMax) > 1:
        xxSec = max(xMax, yMax) + 1
        
    xx = np.linspace(0, xxSec, 100)
    yy = Fun(xx, para[0], para[1], para[2])
    plt.plot(xx, yy, '-g')
    plt.legend(["二次拟合曲线", "MODTRAN模拟散点"])
    
    predictRes = para[0] * x1 * x1 + para[1] * x1 + para[2] # predict
    realResult = y1 # real
    mse = mean_squared_error(realResult, predictRes)
    rmse = np.sqrt(mse)
    r2 = r2_score(realResult, predictRes)
    # 将结果写入文本文件
    f = open("./res/WaterVaporResult.txt", "a")
    f.write("\n" + str(xAxis) + " " + str(yAxis) + "\n" + str(format(para[0], '.3f')) + 'x*x+' + str(format(para[1], '.3f')) + 'x+' + str(format(para[2], '.3f')) + "\n" \
        + " RMSE: " + str(format(rmse, '.3f')) + "\n")
    # 保存
    plt.savefig("./res/" + xAxis + '_' + yAxis + ".png")
    plt.clf()
    
    
def colName2Chinese(axisName):
    # 下划线分割 ——> trans 0 S8
    nameLst = axisName.split('_')
    
    tempStr = ''
    if nameLst[0] == 'trans':
        tempStr = '透过率'
    elif nameLst[0] == 'up':
        tempStr = '上行辐射'
    elif nameLst[0] == 'down':
        tempStr = '下行辐射'
    else:
        tempStr = ''
    
    return nameLst[2] + '波段' + nameLst[1] + '°' + tempStr

# 可能输入列号 也可能输入列名 如果是列名 则进行转换
def transNumber(col, sh):
    if isinstance(col, int):
        return col
    elif isinstance(col, float):
        return int(np.floor(col))
    else:
        for i in range(sh.ncols):
            if col == str(sh.cell(0, i).value):
                return i
    print('Something wrong happens for transNumber')
    sys.exit()

def Fun(x, a1, a2, a3):
    return a1 * x**2 + a2 * x + a3

def error(p, x, y):
    return Fun(p, x)-y

def getIndexes(y_pre, y_data):
    n = len(y_data)
    SSE = ((y_data - y_pre)**2).sum()
    MSE = SSE / n
    RMSE = np.sqrt(MSE)

    u = y_data.mean()
    SST = ((y_data - u)**2).sum()
    SSR = SST - SSE
    R_square = SSR / SST
    return RMSE, R_square

# drawLineChart.drawLineChart(excelFile, 'trans_10_S8', 'up_10_S8')
excelFile = r'E:\integrate\modtran\result\res.xlsx'
baseCol = 'up_10_S8'
targetCol = 'down_55_S9'
weidu = 'midsum'
wv = [0, 2.5]
duibiexcel = r'E:/integrate/data/廓线-水汽对照表.xlsx'
# drawLineChart(excelFile, baseCol, targetCol, weidu, wv, duibiexcel)