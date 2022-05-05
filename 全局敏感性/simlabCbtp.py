import numpy as np
import os
import linecache

'''
simlab全局敏感性分析
读取生成的sam文件数据
导入将要分析的模型中
然后根据其生成所需文件
cbtp模型计算得到的cbtpev cbtpes 有效发射率
'''
def readSam(fileLoc, saveLoc):
    # 行数
    countLine = int(linecache.getline(fileLoc, 2))
    
    # 结果列表
    cbtpevLst = []
    cbtpesLst = []
    
    # 从第五行开始
    for i in range(countLine):
        nowLine = i + 5
        lst = linecache.getline(fileLoc, nowLine).split()
        cbtpev, cbtpes = calCbtpModel(float(lst[0]), float(lst[1]), float(lst[2]), float(lst[3]))
        cbtpevLst.append(cbtpev)
        cbtpesLst.append(cbtpes)
    
    # 写入
    f = open(saveLoc, 'w')
    f.write('2\n')
    f.write('cbtpev\n')
    f.write('cbtpes\n')
    f.write('time = no\n')
    f.write(str(countLine) + '\n')
    for i in range(countLine):
        f.write(str(cbtpevLst[i]) + '\t' + str(cbtpesLst[i]) + '\n')

def calCbtpModel(lai, es, ev, vza):
    # 引入必要的参数
    c1 = 3.7404e8
    c2 = 14387
    
    # 角度转弧度
    vza_rad = vza / 180 * np.pi
    
    # 拦截率
    i0 = 1 - np.exp(-0.5 * lai / np.cos(vza_rad))
    
    # i0hemi 半球平均拦截率，经验公式
    i0hemi = 1 - np.exp(-0.825 * lai)
    
    # start数值积分法计算TIR-P
    fun_up = 0
    fun_down = 0
    
    for j in np.arange(0, float(lai), 0.01):
        fun_up = fun_up + np.exp(-0.5 * j / np.cos(vza_rad)) * 0.5 / np.cos(vza_rad) * 0.5 * np.exp(-0.8 * pow(j, 0.9)) / (1 - np.exp(-0.5 * lai / np.cos(vza_rad))) * 0.01
        fun_down = fun_down + np.exp(-0.5 * j / np.cos(vza_rad)) * 0.5 / np.cos(vza_rad) * 0.5 * np.exp(-0.8 * pow(lai - j, 0.9)) / (1 - np.exp(-0.5 * lai / np.cos(vza_rad))) * 0.01
    # elup 向上逃逸概率 eldown向下逃逸概率
    elup = round(fun_up, 4)
    eldown = round(fun_down, 4)
    
    # p 再碰撞概率
    p = 1 - elup - eldown
    
    # 计算DCE123简化模型
    # c1 光子被叶片反射后碰到叶片的概率
    c1 = (1 - ev) * p                         
    # c2 光子被土壤反射后碰撞叶片的概率
    c2 = (1 - es) * i0hemi
    # c3 光子被叶片反射后碰撞土壤的概率
    c3 = (1 - ev) * eldown
    
    # CBT-P模型 植被和土壤有效发射率
    cbtpev = i0 * ev * (1 + c1 + c1 * c1 + c3 * c2) + (1 - i0) * ev * (c2 + c2 * c1)
    cbtpes = i0 * es * (c3 + c1 * c3) + (1 - i0) * es * (1 + c2 * c3)
    
    return cbtpev, cbtpes

fileLoc = r'D:/simlab/soft/data/cbtp.sam'
saveLoc = r'D:/simlab/soft/data/cbtpRes.txt'
readSam(fileLoc, saveLoc)