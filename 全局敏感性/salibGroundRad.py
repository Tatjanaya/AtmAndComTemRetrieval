import numpy as np
from SALib.analyze import sobol
from SALib.plotting.bar import plot as barplot
import matplotlib.pyplot as plt
from SALib.sample import saltelli

def genResTxt(vza, bandMark):
    problem = {
        'num_vars': 5,
        'names': ['lai', 'es', 'ev', 'Tv', 'Ts'],
        'bounds': [[0.5, 7.0],
                [0.94, 0.96],
                [0.96, 0.99],
                [280, 290],
                [290, 300]]
    }
    
    param_values = saltelli.sample(problem, 100)
    
    r_lst = []
    for p in param_values:
        lai = p[0]
        es = p[1]
        ev = p[2]
        tv = p[3]
        ts = p[4]
        r_lst.append(calGroundRad(lai, es, ev, vza, tv, ts, bandMark))
    
    Si = sobol.analyze(problem, np.array(r_lst))
    
    print('S1:', Si['S1'])
    
def calGroundRad(lai, es, ev, vza, tv, ts, bandMark):
    c1 = 3.7404e8
    c2 = 14387
    
    if bandMark == 'S8':
        band = 10.852
    else:
        band = 12
        
    cbtpev, cbtpes = calCbtpModel(lai, es, ev, vza)
    
    l_tv = c1 / (np.pi * pow(band, 5) * (np.exp(c2 / band / tv) - 1))
    l_ts = c1 / (np.pi * pow(band, 5) * (np.exp(c2 / band / ts) - 1))
    
    return cbtpes * l_ts + cbtpev * l_tv
    
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

vza = 0
bandMark = 'S8'
genResTxt(vza, bandMark)