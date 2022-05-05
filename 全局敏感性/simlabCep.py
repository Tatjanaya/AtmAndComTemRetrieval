import numpy as np
import os
import linecache

'''
cep模型全局敏感性分析
'''
def readSam(fileLoc, saveLoc):
    # 行数
    countLine = int(linecache.getline(fileLoc, 2))
    
    # 结果列表
    ecaoP4sailLst = []
    
    # 从第五行开始
    for i in range(countLine):
        nowLine = i + 5
        lst = linecache.getline(fileLoc, nowLine).split()
        ecaoP4sail = calCep(float(lst[0]), float(lst[1]), float(lst[2]), float(lst[3]))
        ecaoP4sailLst.append(ecaoP4sail)
    
    # 写入
    f = open(saveLoc, 'w')
    f.write('1\n')
    f.write('ecaoP4sail\n')
    f.write('time = no\n')
    f.write(str(countLine) + '\n')
    for i in range(countLine):
        f.write(str(ecaoP4sailLst[i]) + '\n')

def calCep(lai, vza, ev, es):
    
    # 角度转弧度
    vza_rad = vza / 180 * np.pi
    
    # i0 拦截率
    i0 = 1 - np.exp(-0.5 * lai / np.cos(vza_rad))
    
    # start数值积分法计算TIR-P
    fun_up = 0
    fun_down = 0
    
    for j in np.arange(0, float(lai), 0.01):
        fun_up = fun_up + np.exp(-0.5 * j / np.cos(vza_rad)) * 0.5 / np.cos(vza_rad) * 0.5 * np.exp(-0.8 * np.power(j, 0.9)) / (1 - np.exp(-0.5 * lai / np.cos(vza_rad))) * 0.01
        fun_down = fun_down + np.exp(-0.5 * j / np.cos(vza_rad)) * 0.5 / np.cos(vza_rad) * 0.5 * np.exp(-0.8 * np.power(lai - j, 0.9)) / (1 - np.exp(-0.5 * lai / np.cos(vza_rad))) * 0.01
    # elup 向上逃逸概率 eldown向下逃逸概率
    elup = round(fun_up, 4)
    eldown = round(fun_down, 4)
    
    # p 再碰撞概率
    p = 1 - elup - eldown
    
    # cavity 反射后未被植被吸收的部分，可能前向漫反射，也可能后向漫反射
    cavity = 1 - p * (1 - ev)
    
    # i0hemi 半球平均拦截率，经验公式
    i0hemi = 1 - np.exp(-0.825 * lai)
    
    # rc1 冠层前向漫反射率
    # rc2 冠层后向漫反射率
    rc1 = (1 - ev) * eldown / cavity
    rc2 = (1 - ev) * elup / cavity

    # 计算e1 e2 e3 e4 e5
    mulsv = 1 - rc2 * (1 - es) * i0hemi
    
    e1 = i0 * ev / cavity
    e2 = (1 - i0) * (1 - es) * i0hemi * ev / cavity / mulsv
    e3 = i0 * rc1 * (1 - es) * i0hemi * ev / cavity / mulsv
    e4 = (1 - i0) * es / mulsv
    e5 = i0 * rc1 * es / mulsv

    # ecaoP4sail 冠层方向发射率
    ecaoP4sail = e1 + e2 + e3 + e4 + e5
    
    return ecaoP4sail

fileLoc = r'D:/simlab/soft/data/ce-p/cep.sam'
saveLoc = r'D:/simlab/soft/data/ce-p/cepRes.txt'
readSam(fileLoc, saveLoc)