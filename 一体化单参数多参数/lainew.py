import numpy as np
import os
import random
from sklearn.metrics import mean_squared_error
import seaborn as sns
import matplotlib.pyplot as plt

lai_lst = np.arange(0.5, 7.1, 0.1)
vza = 0
ev = 0.98
es = 0.94
niheFile = './data/nihe_wv.txt'
wv = 2.0
tv = 300
ts = 320
bandMark = 'S8'

def main(lai_lst, vza, ev, es, niheFile, wv, tv, ts, bandMark):
    # 
    rmse_lst = []
    for lai in lai_lst:
        # 1000
        right_toa = []
        err_toa = []
        for _ in range(1000):
            toa, toa_err = calToa(lai, vza, ev, es, niheFile, wv, tv, ts, bandMark)
            right_toa.append(toa)
            err_toa.append(toa_err)
        # 计算rmse
        mse = mean_squared_error(right_toa, err_toa)
        rmse = np.sqrt(mse)
        rmse_lst.append(rmse)
        
        print('lai: ' + str(lai) + ' finish')
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 黑体
    plt.rcParams['axes.unicode_minus'] = False    # 解决无法显示符号的问题
    # 坐标轴的刻度设置向内(in)或向外(out)
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    sns.set(style="white", font='SimHei')
    
    plt.xlabel("LAI")
    plt.ylabel("大气层顶温度RMSE")
    plt.grid()
    plt.plot(lai_lst, rmse_lst, color="violet", linestyle='--')
    # plt.savefig("./res/" + "组分温度-模型lai.png", format='png', dpi=300)
    plt.show()
    # plt.clf()

def calToa(lai, vza, ev, es, niheFile, wv, tv, ts, bandMark):
    # 必要的参数
    c1 = 3.7404e8
    c2 = 14387
    
    if bandMark == 'S8':
        band = 10.852
    else:
        band = 12
        
    cep = calCep(lai, vza, ev, es)
    cbtpev, cbtpes = calCbtpModel(lai, es, ev, vza)
    
    trans_wvlst, up_wvlst, down_wvlst = readNiheFile(niheFile, vza, bandMark)
    
    trans = trans_wvlst[0] * wv ** 2 + trans_wvlst[1] * wv + trans_wvlst[2]
    up = up_wvlst[0] * wv ** 2 + up_wvlst[1] * wv + up_wvlst[2]
    down = down_wvlst[0] * wv ** 2 + down_wvlst[1] * wv + down_wvlst[2]
    
    l_tv = c1 / (np.pi * pow(band, 5) * (np.exp(c2 / band / tv) - 1))
    l_ts = c1 / (np.pi * pow(band, 5) * (np.exp(c2 / band / ts) - 1))
    
    # 正确的ltoa
    ltoa = trans * ((1 - cep) * down + cbtpev * l_ts + cbtpev * l_tv) + up
    
    # 叠加误差的lai
    if lai <= 0.5:
        lai_err = lai - 0.4 + random.random()
    else:
        lai_err = lai - 0.5 + random.random()
    
    # 计算cep cbtp
    cep_err = calCep(lai_err, vza, ev, es)
    cbtpev_err, cbtpes_err = calCbtpModel(lai_err, es, ev, vza)
    
    ltoa_err = trans * ((1 - cep_err) * down + cbtpev_err * l_ts + cbtpev_err * l_tv) + up
    
    toa = c2 / band / np.log((c1 / (np.pi * pow(band, 5) * ltoa)) + 1)
    toa_err = c2 / band / np.log((c1 / (np.pi * pow(band, 5) * ltoa_err)) + 1)
    
    return toa, toa_err
    
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

main(lai_lst, vza, ev, es, niheFile, wv, tv, ts, bandMark)