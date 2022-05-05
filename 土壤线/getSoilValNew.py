from cProfile import label
from tkinter import N
from turtle import color, right
import gdal
import numpy as np
from scipy import optimize
import sys
import seaborn as sns
import matplotlib.pyplot as plt

'''
获取土壤线结果
'''
# nirband redband greenband 以0为开头计算波段位置
def getSoilVal(nir, red, row, col, N, targetLen):
    # 读取三波段图像
    # nir, red, green, row, col = getNirRedGreenData(originTif, nirband, redband, greenband)
    # 计算步长 这里没使用步长 直接取1
    step = int(getStep(row, col, N))
    print('step: ' + str(step))
    step = 1
    # 绘制特征空间 -> 组成一个列表 [(r, n, g), (), (), ...]
    # 元组顺序解释 -> sorted 按照元组第一个元素进行排序 而 red是横坐标
    data_lst = getTupleLst(red, nir, row, col, step)
    print('getTupleLst finish')
    # 获取初始土壤点集 根据横坐标 r 最大值最小值分为256组 此处不这么干的原因是会导致过多的空子集出现
    # group_256 = groupDatas(data_lst)
    # print('groupDatas finish')
    # 获取相关系数 确定选择哪个分组
    rangeStr = getRelePara(data_lst)
    # rangeStr = '25-75'
    print('rangeStr: ' + str(rangeStr))
    print('getRelePara finish')
    # 根据返回的字符串 截取group_256子列表
    if rangeStr == '0-100':
        cal_group_lst = data_lst[:]
        groupNum = 256
    elif rangeStr == '0-50':
        cal_group_lst = data_lst[: int(len(data_lst) * 0.5) + 1]
        groupNum = 128
    elif rangeStr == '0-75':
        cal_group_lst = data_lst[: int(len(data_lst) * 0.75) + 1]
        groupNum = 192
    elif rangeStr == '25-75':
        cal_group_lst = data_lst[int(len(data_lst) * 0.25): int(len(data_lst) * 0.75) + 1]
        groupNum = 128
    elif rangeStr == '25-100':
        cal_group_lst = data_lst[int(len(data_lst) * 0.25):]
        groupNum = 192
    else:
        print('Rele part has something wrong')
        sys.exit()
    # 接下来的计算以cal_group_lst为准
    # nir red green lst     r n g
    if len(cal_group_lst) == 0:
        print('cal_group_lst is []')
        sys.exit()
    
    # 将cal_group_lst 分为256组
    group_256 = groupDatas(cal_group_lst, groupNum)
    # 去掉空列表
    group_256 = list(filter(None, group_256))
    # 输出每组nir最小值
    min_nir_lst = findMinNirAsResLst(group_256)
    print('findMinNirAsResLst finish')
    
    red_lst = [x[0] for x in min_nir_lst]
    nir_lst = [x[1] for x in min_nir_lst]
    
    while len(red_lst) > targetLen:
        A, B = optimize.curve_fit(funFit, red_lst, nir_lst)[0]
        indexLoc = -1
        max_val = -1
        for i in range(len(red_lst)):
            if calDis(A, B, red_lst, nir_lst, i) > max_val:
                max_val = calDis(A, B, red_lst, nir_lst, i)
                indexLoc = i
        
        del red_lst[indexLoc]
        del nir_lst[indexLoc]
        
    print('targetLen finish')
    
    # 平均
    nir_res = sum(nir_lst) / len(nir_lst)
    red_res = sum(red_lst) / len(red_lst)
    
    # print('red: ' + str(red_res))
    # print('nir: ' + str(nir_res))
    
    # 绘图
    # plt.rcParams['font.sans-serif'] = ['SimHei']  # 黑体
    # plt.rcParams['axes.unicode_minus'] = False    # 解决无法显示符号的问题
    # # 坐标轴的刻度设置向内(in)或向外(out)
    # plt.rcParams['xtick.direction'] = 'in'
    # plt.rcParams['ytick.direction'] = 'in'
    # sns.set(style="white", font='SimHei')
    
    # drawPic(red_lst, nir_lst, red.flatten(), nir.flatten())
    return red_lst, nir_lst

# 找到每个分组里nir最小的值 参数 group_lst = [[(r, n, g), (r, n, g)], []]
def findMinNirAsResLst(group_lst):
    res_lst = []
    for group_zi_lst in group_lst:
        minNir = 1000000
        loc = -1
        for i in range(len(group_zi_lst)):
            if group_zi_lst[i][1] < minNir:
                minNir = group_zi_lst[i][1]
                loc = i
        res_lst.append(group_zi_lst[loc])
    return res_lst
        
# A B 是直线方程系数
# xLst yLst 是x y值列表
# i是列表中的索引
def calDis(A, B, xLst, yLst, i):
    return np.abs(A * xLst[i] - yLst[i] + B) / np.sqrt(A * A + 1)

# 相关系数获取
def getRelePara(data_lst):
    xLst_0_50 = []
    xLst_0_75 = []
    xLst_0_100 = []
    xLst_25_75 = []
    xLst_25_100 = []
    
    yLst_0_50 = []
    yLst_0_75 = []
    yLst_0_100 = []
    yLst_25_75 = []
    yLst_25_100 = []
    # 位置
    loc_25 = int(len(data_lst) * 0.25)
    loc_50 = int(len(data_lst) * 0.50)
    loc_75 = int(len(data_lst) * 0.75)
    
    for i in range(len(data_lst)):
        xLst_0_100.append(data_lst[i][0])
        yLst_0_100.append(data_lst[i][1])
        if i < loc_50:
            xLst_0_50.append(data_lst[i][0])
            yLst_0_50.append(data_lst[i][1])
        if i < loc_75:
            xLst_0_75.append(data_lst[i][0])
            yLst_0_75.append(data_lst[i][1])
        if i >= loc_25 and i < loc_75:
            xLst_25_75.append(data_lst[i][0])
            yLst_25_75.append(data_lst[i][1])
        if i >= loc_25:
            xLst_25_100.append(data_lst[i][0])
            yLst_25_100.append(data_lst[i][1])
            
    lst = []
    lst.append(str(rela(xLst_0_100, yLst_0_100)))
    lst.append(str(rela(xLst_0_50, yLst_0_50)))
    lst.append(str(rela(xLst_0_75, yLst_0_75)))
    lst.append(str(rela(xLst_25_75, yLst_25_75)))
    lst.append(str(rela(xLst_25_100, yLst_25_100)))
    
    # rela_0_100 rela_0_50 rela_0_75 rela_25_75 rela_25_100
    # 找最大值
    res_loc = getRangeStr(lst)
    if res_loc == 0:
        return '0-100'
    elif res_loc == 1:
        return '0-50'
    elif res_loc == 2:
        return '0-75'
    elif res_loc == 3:
        return '25-75'
    elif res_loc == 4:
        return '25-100'
    else:
        return ''

# 找到斜率小于1中最小的releLst元素的位置
def getRangeStr(releLst):
    res_loc = -1
    minVal = -1
    for i in range(len(releLst)):
        if float(releLst[i]) > minVal:
            res_loc = i
            minVal = float(releLst[i])
    return res_loc

def funFit(x, A, B):
    return A * x + B

# 拟合结果
def fitRes(xLst, yLst):
    k, _ = optimize.curve_fit(funFit, xLst, yLst)[0]
    return k

def rela(xLst, yLst):
    avgX = sum(xLst) / len(xLst)
    avgY = sum(yLst) / len(yLst)
    upTotal = 0
    downXTotal = 0
    downYTotal = 0
    for i in range(len(xLst)):
        upTotal += (xLst[i] - avgX) * (yLst[i] - avgY)
        downXTotal += (xLst[i] - avgX) ** 2
        downYTotal += (yLst[i] - avgY) ** 2
    return upTotal / np.sqrt(downXTotal) / np.sqrt(downYTotal)    

# 分组
def groupDatas(data_lst, groupNum):
    # 取最大值和最小值
    minRedVal = data_lst[0][0]
    maxRedVal = data_lst[-1][0]
    # 得到分组步长
    groupStep = (maxRedVal - minRedVal) / (groupNum - 1)
    # 填充
    res_lst = []
    mark = 0
    left = minRedVal
    for i in range(groupNum):
        # 处理
        temp_lst = []
        if i == groupNum - 1:
            right = maxRedVal + 1
        else:
            right = left + groupStep
        while True:
            if mark >= len(data_lst) or data_lst[mark][0] >= right:
                break
            temp_lst.append(data_lst[mark])
            mark += 1
        res_lst.append(temp_lst)
        left = right
    # res_lst 结构 [[(r, n, g), ()], [], ...]
    return res_lst

# 获取元组列表
def getTupleLst(redData, nirData, row, col, step):
    lst = []
    if step == 1:
        for i in range(row):
            for j in range(col):
                # 含0的为无影像区域 需要全部去除
                # 现在改成都小于500 没参考价值应该去除
                if redData[i][j] > nirData[i][j]:
                    continue
                if redData[i][j] <= 500 and nirData[i][j] <= 500:
                    continue
                lst.append((redData[i][j], nirData[i][j]))
    else:
        # 按照step步长求取
        ix = 0
        iy = 0
        while ix < row and iy < col:
            if redData[ix][iy] != 0 and nirData[ix][iy] != 0:
                lst.append((redData[ix][iy], nirData[ix][iy]))
            ix += int(step)
            iy += int(step)
    return sorted(lst)        

def getNirRedGreenData(originTif, nirband, redband, greenband):
    # 读数据
    ds = gdal.Open(originTif)
    col = ds.RasterXSize
    row = ds.RasterYSize
    
    data = ds.ReadAsArray(0, 0, col, row)
    nir = data[nirband, 0: row, 0: col]
    red = data[redband, 0: row, 0: col]
    green = data[greenband, 0: row, 0: col]
    
    return nir, red, green, row, col

# 计算步长 计算方法如下
# 若 M ≥ N 根号 M / N
# 若 M ＜ N 1
# 通常情况下 N 取 2 000 000
def getStep(row, col, bestN):
    M = row * col
    if M >= bestN:
        return np.sqrt(M / bestN)
    else:
        return 1
    
def drawPic(selectLstR, selectLstNir, totalLstR, totalLstNir):
    plt.xlabel("R")
    plt.ylabel("NIR")
    plt.grid()
    plt.scatter(totalLstR, totalLstNir, color='red', s=5)
    plt.scatter(selectLstR, selectLstNir, color='green', s=5)
    # plt.legend(loc='upper right')
    plt.show()
    plt.clf()
