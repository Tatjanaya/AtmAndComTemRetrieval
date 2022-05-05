import os
import linecache
import pandas as pd
""" 提取各条廓线水汽含量

@extractWaterVapor: 提取廓线编号和水汽含量 保存在一个excel文件中
tigrFile: TIGR大气廓线文件位置

@lineCount: 统计文件行数
fileName: 文件名字
"""
def extractWaterVapor(tigrFile):
    # 如果不存在文件夹 则创建
    dirs = './result/'
    if not os.path.exists(dirs):
        os.makedirs(dirs)
    
    start = 0
    tigrNum = []
    waterVpor = []
    
    excelHeaders = ["TIGR编号", "水汽含量(g cm-2)"]
    
    rowNum = lineCount(tigrFile)
    while start * 43 + 1 < rowNum:
        line = linecache.getline(tigrFile, start * 43 + 1).strip()
        # 找到TIGR和TCWV位置
        tigrLoc = line.find('TIGR', 1)
        tcwvLoc = line.find('TCWV')
        
        # 添加
        tigrNum.append(line[tigrLoc + 4: tigrLoc + 8])
        waterVpor.append(line[tcwvLoc + 5: tcwvLoc + 9])
        
        start += 1
    
    # 合并列表
    excelMatrix = list(zip(tigrNum, waterVpor))
    # 排序
    # excelMatrix = sorted(excelMatrix, key=(lambda x : x[1]))
            
    # 保存到excel表中
    excelFileName = '廓线-水汽对照表'
    writer = pd.ExcelWriter("./result/" + excelFileName + ".xlsx")
    data_df = pd.DataFrame(data=excelMatrix, columns=excelHeaders)
    data_df.to_excel(writer, 'page_1', float_format='%.5f', index=False, header=True)
    writer.save()
    
def lineCount(fileName):
    count = 0
    f = open(fileName, 'r')
    while True:
        buffer = f.read(1024 * 8192)
        if not buffer:
            break
        count += buffer.count('\n')
    f.close()
    return count