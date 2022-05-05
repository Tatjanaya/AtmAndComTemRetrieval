from fnmatch import translate
import os
import sys
import linecache
import pandas as pd
import numpy as np
""" 从chn文件中提取有用信息 并储存在一个excel文件中

@extractChnData: 从所有chn文件中提取数据 并存储到一个excel文件中
excelDir: excel存放位置 所有数据都只存到一个excel文件中
chnDir: 透过率chn文件储存位置
chn文件默认设定为Sentinel-3 SLSTR S8 S9
originalFile: 原始数据，用于提取水汽含量
"""
def extractChnData(excelDir, chnDir, originalFile):
    # 提取chnDir文件夹下所有后缀是chn且前缀是trans up down的文件
    transLst = []
    upLst = []
    downLst = []
    for _, _, k in os.walk(chnDir):
        for fileName in k:
            if fileName.endswith('.chn') and fileName.startswith('trans'):
                transLst.append(fileName)
            if fileName.endswith('.chn') and fileName.startswith('up'):
                upLst.append(fileName)
            if fileName.endswith('.chn') and fileName.startswith('down'):
                downLst.append(fileName)
    
    # 根据角度大小排序
    transLst = sortByAngle(transLst)
    upLst = sortByAngle(upLst)
    downLst = sortByAngle(downLst)
    
    # 合并列表 按透过率 上行辐射 下行辐射的顺序排列
    totalLst = transLst + upLst + downLst
    
    # 总的数量 两个波段 * 2
    total = len(totalLst) * 2 + 1
    
    # 分配空间 待会要存储在excel表里的矩阵大小
    excelMatrix = np.zeros((946, total))
    
    # 水汽含量
    watervapor = extractWaterVapor(originalFile)
    
    # 表头 表头以 trans/up/down_angle_S8/S9 来命名
    excelHeaders = []
    for ele in totalLst:
        excelHeaders.append(ele.split('.')[0] + '_S8')
        excelHeaders.append(ele.split('.')[0] + '_S9')
        
    excelHeaders.append('waterVapor')
    
    # location 目前写入数据的列数
    col = 0
    
    # 将这些文件的数据写入excel表中
    writer = pd.ExcelWriter(excelDir + "\\" + "res.xlsx")
    for chnFile in totalLst:
        # 打开文件 读取全部行
        f = open(chnDir + "\\" + chnFile, 'r',encoding='utf-8')
        chnData = f.readlines()
        f.close()
        
        # 判断读取类型 一类为透过率 一类为上下行辐射
        if chnFile.startswith('trans'):
            for row in range(946):
                excelMatrix[row][col] = 1.0 - float(chnData[row * 7 + 5].split()[2])
                excelMatrix[row][col + 1] = 1.0 - float(chnData[row * 7 + 6].split()[2])
        elif chnFile.startswith('down') or chnFile.startswith('up'):
            for row in range(946):
                excelMatrix[row][col] = float(chnData[row * 7 + 5].split()[3]) * 10000
                excelMatrix[row][col + 1] = float(chnData[row * 7 + 6].split()[3]) * 10000
        else:
            print('Something wrong happens')
            sys.exit()
        
        col = col + 2
    
    for row in range(946):
        excelMatrix[row][-1] = float(watervapor[row])
    
    # 保存
    data_df = pd.DataFrame(data=excelMatrix, columns=excelHeaders)
    data_df.to_excel(writer, 'page_1', float_format='%.5f', index=False, header=True)
    writer.save()
    
def extractWaterVapor(originalFile):
    # 获得文件行数
    linecount = lineCount(originalFile)
    # 水汽存储列表
    waterLst = []
    index = 0
    while index * 43 + 1 < linecount:
        line = linecache.getline(originalFile, index * 43 + 1).strip()
        # 找到TCWV位置
        tcwvLoc = line.find('TCWV')
        # 添加
        waterLst.append(line[tcwvLoc + 5: tcwvLoc + 9])
        index += 1
    return waterLst
    
# 根据角度大小排序
def sortByAngle(lst):
    return sorted(lst, key=lambda x: int(x.split('.')[0].split('_')[1]))
    
# 统计文件行数 不用readlines一次性读取是为了防止文件过大导致运行速度过慢
def lineCount(fileName):
    count = 0
    for _, _ in enumerate(open(fileName, 'r')):
        count += 1
    return count

# extractWaterVapor('TIGR_atmProfilesSelection.txt')