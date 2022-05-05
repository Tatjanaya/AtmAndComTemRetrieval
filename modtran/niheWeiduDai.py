import numpy as np
import xlrd
import sys
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score

# 拟合不同纬度带水汽关系

def calwvRele(excelFile, baseCol, targetCol, weidu, duibiexcel):
    if weidu == 'trop':
        minHao = 1
        maxHao = 872
    elif weidu == 'midsum':
        minHao = 873
        maxHao = 1260
    elif weidu == 'midwin':
        minHao = 1261
        maxHao = 1614
    elif weidu == 'polarsum':
        minHao = 1615
        maxHao = 1718
    elif weidu == 'polarwin':
        minHao = 1719
        maxHao = 2311
    else:
        minHao = 873
        maxHao = 1614
    
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
    
    # 两列数据取出
    xData = []
    yData = []
    for i in range(1, sh.nrows):
        if bianhao_lst[i - 1] >= minHao and bianhao_lst[i - 1] <= maxHao:
            xData.append(float(sh.cell(i, baseCol).value))
            yData.append(float(sh.cell(i, targetCol).value))
        
    # 拟合
    x1 = np.array(xData)
    y1 = np.array(yData)
    para, _ = curve_fit(Fun, x1, y1)
    
    predictRes = para[0] * x1 * x1 + para[1] * x1 + para[2] # predict
    realResult = y1 # real
    mse = mean_squared_error(realResult, predictRes)
    rmse = np.sqrt(mse)
    
    f = open("./res/通道相关性计算所需.txt", "a")
    f.write("\n" + str(xAxis) + " " + str(yAxis) + " " + str(weidu) + "\n" + str(format(para[0], '.3f')) + ' ' + str(format(para[1], '.3f')) + ' ' + str(format(para[2], '.3f')) + "\n" \
        + " RMSE: " + str(format(rmse, '.3f')) + "\n")
    
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

excelFile = r'E:/integrate/modtran/result/res.xlsx'
duibiexcel = r'E:/integrate/data/廓线-水汽对照表.xlsx'
weidu_lst = ['trop', 'midsum', 'midwin', 'polarsum', 'polarwin']
vza_lst = range(0, 60, 5)

for weidu in weidu_lst:
    for vza in vza_lst:
        calwvRele(excelFile, 'trans_' + str(vza) + '_S8', 'up_' + str(vza) + '_S8', weidu, duibiexcel)
        calwvRele(excelFile, 'trans_' + str(vza) + '_S9', 'up_' + str(vza) + '_S9', weidu, duibiexcel)
#         calwvRele(excelFile, 'waterVapor', 'trans_' + str(vza) + '_S8', weidu, duibiexcel)
#         calwvRele(excelFile, 'waterVapor', 'trans_' + str(vza) + '_S9', weidu, duibiexcel)
        
#         calwvRele(excelFile, 'waterVapor', 'up_' + str(vza) + '_S8', weidu, duibiexcel)
#         calwvRele(excelFile, 'waterVapor', 'up_' + str(vza) + '_S9', weidu, duibiexcel)
        
#         calwvRele(excelFile, 'waterVapor', 'down_' + str(vza) + '_S8', weidu, duibiexcel)
#         calwvRele(excelFile, 'waterVapor', 'down_' + str(vza) + '_S9', weidu, duibiexcel)

# 提取对比excel
# duibiFile = xlrd.open_workbook(duibiexcel)
# shduibi = duibiFile.sheet_by_name('page_1')
# bianhao_lst = []
# wv_lst = []
# trop_sum = 0
# midsum_sum = 0
# midwin_sum = 0
# polarsum_sum = 0
# polarwin_sum = 0
# q1 = 0
# q2 = 0
# q3 = 0
# q4 = 0
# q5 = 0
# for i in range(1, shduibi.nrows):
#     bianhao_lst.append(int(shduibi.cell(i, 0).value))
#     wv_lst.append(float(shduibi.cell(i, 1).value))
    
# 算平均
# for i in range(len(bianhao_lst)):
#     if bianhao_lst[i] >= 1 and bianhao_lst[i] <= 872:
#         trop_sum += wv_lst[i]
#         q1 += 1
#     elif bianhao_lst[i] >= 873 and bianhao_lst[i] <= 1260:
#         midsum_sum += wv_lst[i]
#         q2 += 1
#     elif bianhao_lst[i] >= 1261 and bianhao_lst[i] <= 1614:
#         midwin_sum += wv_lst[i]
#         q3 += 1
#     elif bianhao_lst[i] >= 1615 and bianhao_lst[i] <= 1718:
#         polarsum_sum += wv_lst[i]
#         q4 += 1
#     else:
#         polarwin_sum += wv_lst[i]
#         q5 += 1
        
# print(trop_sum / q1)
# print(midsum_sum / q2)
# print(midwin_sum / q3)
# print(polarsum_sum / q4)
# print(polarwin_sum / q5)        