import numpy as np
from SALib.analyze import sobol
from SALib.plotting.bar import plot as barplot
import matplotlib.pyplot as plt
from SALib.sample import saltelli

def genResTxt(paraLst, vza, niheFile, bandMark, tarLoc):
    problem = {
        'num_vars': 6,
        'names': ['lai', 'es', 'ev', 'Tv', 'Ts', 'wv'],
        'bounds': [[0.5, 7.0],
                [0.94, 0.96],
                [0.96, 0.99],
                [280, 290],
                [290, 300],
                [0.1, 2.5]]
    }
    res_lst, lenp = calResLst(paraLst, vza, niheFile, bandMark)
    param_values = saltelli.sample(problem, 100)
    r_lst = []
    for p in param_values:
        lai = p[0]
        es = p[1]
        ev = p[2]
        tv = p[3]
        ts = p[4]
        wv = p[5]
        r_lst.append(calToaRes(lai, ev, es, tv, ts, wv, niheFile, vza, bandMark))
    
    # f = open(tarLoc, 'w')
    # f.write('1\n')
    # f.write('ltoa\n')
    # f.write('time = no\n')
    # f.write(str(lenp) + '\n')
    # for i in range(lenp):
    #     f.write(str(res_lst[i]) + '\n')
    
    Si = sobol.analyze(problem, np.array(r_lst), calc_second_order=False)
    
    print('S1:', Si['S1'])
    # Si_df = Si.to_df()
    # barplot(Si_df[0])
    # plt.show()
    # plt.savefig('./res/' + str(vza) + '_' + str(bandMark) + '.png')
    print('ST', Si['ST'])
    
    f = open('./res/salib_res_salib.txt', 'a')
    f.write('vza: ' + str(vza) + ' ' + 'band: ' + str(bandMark))
    f.write('\n')
    f.write('lai es ev Tv Ts wv')
    f.write('\n')
    f.write('S1: ' + str(Si['S1']))
    f.write('\n')
    f.write('S1_conf' + str(Si['S1_conf']))
    f.write('\n')
    f.write('ST: ' + str(Si['ST']))
    f.write('\n')
    f.write('ST_conf: ' + str(Si['ST_conf']))
    f.write('\n')

# 假定了就是0 lai es ev tv ts wv
def calResLst(samFile, vza, niheFile, bandMark):
    paraLst = readSam(samFile)
    res_lst = []
    for para in paraLst:
        lai = para[0]
        es = para[1]
        ev = para[2]
        tv = para[3]
        ts = para[4]
        wv = para[5]
        l_toa = calToaRes(lai, ev, es, tv, ts, wv, niheFile, vza, bandMark)
        res_lst.append(l_toa)
    return res_lst, len(paraLst)
    
# 获取四个TOA
def calToaRes(lai, ev, es, tv, ts, wv, niheFile, vza, bandMark):
    # 必要的参数
    c1 = 3.7404e8
    c2 = 14387

    if bandMark == 'S8':
        band = 10.852
    else:
        band = 12
    # 首先计算cep cbtp
    cep = calCep(lai, vza, ev, es)
    
    cbtpev, cbtpes = calCbtpModel(lai, es, ev, vza)
    
    # 接着计算大气参数
    trans_lst, up_lst, down_lst = readNiheFile(niheFile, vza, bandMark)
    
    t = trans_lst[0] * wv ** 2 + trans_lst[1] * wv + trans_lst[2] 
    u = up_lst[0] * wv ** 2 + up_lst[1] * wv + up_lst[2]
    d = down_lst[0] * wv ** 2 + down_lst[1] * wv + down_lst[2]
    
    # T -> L
    l_tv = c1 / (np.pi * pow(band, 5) * (np.exp(c2 / band / tv) - 1))
    l_ts = c1 / (np.pi * pow(band, 5) * (np.exp(c2 / band / ts) - 1))
    
    # 计算toa
    l_toa = t * ((1 - cep) * d + cbtpes * l_ts + cbtpev * l_tv) + u
    
    return l_toa

def readNiheFile(niheFile, angle, bandMark):
    f = open(niheFile, 'r')
    lines = f.readlines()
    # 找到waterVapor_trans_angle_bandMark
    tranStr = 'waterVapor_trans_' + str(angle) + '_' + str(bandMark)
    # 找到waterVapor_up_angle_bandMark
    upStr = 'waterVapor_up_' + str(angle) + '_' + str(bandMark)
    # 找到waterVapor_down_angle_bandMark
    downStr = 'waterVapor_down_' + str(angle) + '_' + str(bandMark)
    for i in range(len(lines)):
        temp_line = lines[i].strip()
        if tranStr == temp_line:
            resTran = lines[i + 2].split(':')[1].strip()
            tran_lst = resTran.strip().split()
            tran_lst = [float(x) for x in tran_lst]
        if upStr == temp_line:
            resUp = lines[i + 2].split(':')[1].strip() 
            up_lst = resUp.strip().split()
            up_lst = [float(x) for x in up_lst]
        if downStr == temp_line:
            resDown = lines[i + 2].split(':')[1].strip()
            down_lst = resDown.strip().split()
            down_lst = [float(x) for x in down_lst]
    # print(tran_lst)
    # print(up_lst)
    # print(down_lst)
    return tran_lst, up_lst, down_lst

def readSam(samFile):
    f = open(samFile, 'r')
    lines = f.readlines()
    # 读取行数
    count = int(lines[1].strip())
    # 读取后文
    paraLst = []
    for i in range(count):
        tempLst = lines[4 + i].strip().split()
        tempLst = [float(x) for x in tempLst]
        paraLst.append(tempLst)
    return paraLst

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

samFile = './newdata/total1.sam'
# vza = 0
# bandMark = 'S8'
tarLoc = './newdata/res.txt'
niheFile = './newdata/nihe_wv.txt'

vza_lst = range(0, 60, 5)
for vza in vza_lst:
    genResTxt(samFile, vza, niheFile, 'S8', tarLoc)
    genResTxt(samFile, vza, niheFile, 'S9', tarLoc)

# lai = 1
# es = 0.94
# ev = 0.98
# cbtpev, cbtpes = calCbtpModel(lai, es, ev, vza)
# print(cbtpev)
# print(cbtpes)