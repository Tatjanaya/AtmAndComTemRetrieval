import re
from symtable import Symbol
import numpy as np
from sklearn.metrics import mean_squared_error
import seaborn as sns
import os
import sympy
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# 劈窗算法的相对误差
# 思路 计算平均相对偏差
# 平均相对偏差是地表辐射亮温反演值和真值换算为地表辐亮度 计算二者绝对误差占地表辐射亮度真值比例的均值
# 以0 55为例 主要探讨不同水汽含量分组的RMSE随Nadir方向观测天顶角变化情况
# 以0-2.5水汽含量为例
# LAI 0.1-6.0
# es 0.94 
# ev 0.98
# tv ts 300 305 / 300 310
# nadir oblique 0 55

# 首先根据上述条件推断出正确的TOA 然后根据错误的劈窗参数反推得出错误的地面热辐射 再解算出Ts Tv
# 劈窗算法分为 可能 Nadir角度出错 可能 水汽含量分组 出错
# 这里先计算Nadir角度出错
# 从 5 - 30

# 拟合参数结果为(结构为a * x + b * (x - y) + c * (x - y)**2 + d): 
# 参数为: nadir: 0 oblique: 55 lowVapor: 0 highVapor: 2.5
# S8波段0°: [ 1.0174  0.6045  0.2652 -3.7037]
# S8波段55°: [ 1.0223  0.9217  0.2249 -4.3483]
# S9波段0°: [ 1.0292 -0.0667  0.4339 -6.7047]
# S9波段55°: [ 1.035   0.3792  0.3414 -7.5637]

# 拟合参数结果为(结构为a * x + b * (x - y) + c * (x - y)**2 + d): 
# 参数为: nadir: 5 oblique: 55 lowVapor: 0 highVapor: 2.5
# S8波段5°: [ 1.0175  0.6063  0.2649 -3.7083]
# S8波段55°: [ 1.0223  0.9217  0.2249 -4.3483]
# S9波段5°: [ 1.0292 -0.0641  0.4333 -6.7113]
# S9波段55°: [ 1.035   0.3792  0.3414 -7.5637]

# 拟合参数结果为(结构为a * x + b * (x - y) + c * (x - y)**2 + d): 
# 参数为: nadir: 10 oblique: 55 lowVapor: 0 highVapor: 2.5
# S8波段10°: [ 1.0176  0.6118  0.2641 -3.7216]
# S8波段55°: [ 1.0223  0.9217  0.2249 -4.3483]
# S9波段10°: [ 1.0294 -0.0562  0.4314 -6.7311]
# S9波段55°: [ 1.035   0.3792  0.3414 -7.5637]

# 拟合参数结果为(结构为a * x + b * (x - y) + c * (x - y)**2 + d): 
# 参数为: nadir: 15 oblique: 55 lowVapor: 0 highVapor: 2.5
# S8波段15°: [ 1.0177  0.6211  0.2628 -3.7451]
# S8波段55°: [ 1.0223  0.9217  0.2249 -4.3483]
# S9波段15°: [ 1.0296 -0.043   0.4283 -6.7649]
# S9波段55°: [ 1.035   0.3792  0.3414 -7.5637]

# 拟合参数结果为(结构为a * x + b * (x - y) + c * (x - y)**2 + d): 
# 参数为: nadir: 20 oblique: 55 lowVapor: 0 highVapor: 2.5
# S8波段20°: [ 1.0179  0.6346  0.2608 -3.7782]
# S8波段55°: [ 1.0223  0.9217  0.2249 -4.3483]
# S9波段20°: [ 1.0298 -0.0237  0.4237 -6.8122]
# S9波段55°: [ 1.035   0.3792  0.3414 -7.5637]

# 拟合参数结果为(结构为a * x + b * (x - y) + c * (x - y)**2 + d): 
# 参数为: nadir: 25 oblique: 55 lowVapor: 0 highVapor: 2.5
# S8波段25°: [ 1.0182  0.6526  0.2581 -3.8209]
# S8波段55°: [ 1.0223  0.9217  0.2249 -4.3483]
# S9波段25°: [ 1.0302e+00  1.9000e-03  4.1770e-01 -6.8736e+00]
# S9波段55°: [ 1.035   0.3792  0.3414 -7.5637]

# 拟合参数结果为(结构为a * x + b * (x - y) + c * (x - y)**2 + d): 
# 参数为: nadir: 30 oblique: 55 lowVapor: 0 highVapor: 2.5
# S8波段30°: [ 1.0186  0.6759  0.2548 -3.8747]
# S8波段55°: [ 1.0223  0.9217  0.2249 -4.3483]
# S9波段30°: [ 1.0307  0.0351  0.41   -6.9495]
# S9波段55°: [ 1.035   0.3792  0.3414 -7.5637]
def calSplitWindow(lai, ev_real, es_real, nadir_angle, oblique_angle, Tv_real, Ts_real, para_real, para_err):
    # 没文件夹创建文件夹
    dirs = './res/'
    if not os.path.exists(dirs):
        os.makedirs(dirs)
        
    # 引入必要的参数
    c1 = 3.7404e8   
    c2 = 14387

    band1 = 10.852
    band2 = 12
    
    es_8 = es_real
    es_9 = es_real
    ev_8 = ev_real
    ev_9 = ev_real
    
    # 两个角度的弧度值
    sita = [nadir_angle / 180 * np.pi, oblique_angle / 180 * np.pi]
    
    # 真实Tv Ts 辐亮度
    L_v_real_8 = c1 / (np.pi * pow(band1, 5) * (np.exp(c2 / band1 / Tv_real) - 1))
    L_s_real_8 = c1 / (np.pi * pow(band1, 5) * (np.exp(c2 / band1 / Ts_real) - 1))
    L_v_real_9 = c1 / (np.pi * pow(band2, 5) * (np.exp(c2 / band2 / Tv_real) - 1))
    L_s_real_9 = c1 / (np.pi * pow(band2, 5) * (np.exp(c2 / band2 / Ts_real) - 1))
    
    # 先根据正确的lai和正确的ts tv计算正确的地面热辐射
    cbtpev_8_ang_1, cbtpev_8_ang_2, cbtpev_9_ang_1, cbtpev_9_ang_2, cbtpes_8_ang_1, cbtpes_8_ang_2, cbtpes_9_ang_1, cbtpes_9_ang_2 = \
            calCbtp(sita, lai, es_8, es_9, ev_8, ev_9)
    # 反推地面热辐射 -> 理论上这是我们从卫星影像上根据劈窗法得到的结果
    ground_8_ang_1 = cbtpev_8_ang_1 * L_v_real_8 + cbtpes_8_ang_1 * L_s_real_8
    ground_8_ang_2 = cbtpev_8_ang_2 * L_v_real_8 + cbtpes_8_ang_2 * L_s_real_8
    ground_9_ang_1 = cbtpev_9_ang_1 * L_v_real_9 + cbtpes_9_ang_1 * L_s_real_9
    ground_9_ang_2 = cbtpev_9_ang_2 * L_v_real_9 + cbtpes_9_ang_2 * L_s_real_9
    
    # 普朗克公式得到地表亮温
    bt_8_n = c2 / band1 / np.log((c1 / (np.pi * pow(band1, 5) * ground_8_ang_1)) + 1)
    bt_8_o = c2 / band1 / np.log((c1 / (np.pi * pow(band1, 5) * ground_8_ang_2)) + 1)
    bt_9_n = c2 / band2 / np.log((c1 / (np.pi * pow(band2, 5) * ground_9_ang_1)) + 1)
    bt_9_o = c2 / band2 / np.log((c1 / (np.pi * pow(band2, 5) * ground_9_ang_2)) + 1)
    
    # 根据计算参数得到正确的TOA
    real_1_lst = para_real[0]
    real_2_lst = para_real[1]
    real_3_lst = para_real[2]
    real_4_lst = para_real[3]
    
    # TOA四个温度值 这四个温度值是经由正确条件推导出来的
    # print(calFunByFsolve(real_1_lst, real_3_lst, bt_8_n, bt_9_n))
    # print(calFunByFsolve(real_2_lst, real_4_lst, bt_8_o, bt_9_o))
    t_8_n = float(calFunByFsolve(real_1_lst, real_3_lst, bt_8_n, bt_9_n)[0])
    t_9_n = float(calFunByFsolve(real_1_lst, real_3_lst, bt_8_n, bt_9_n)[1])
    t_8_o = float(calFunByFsolve(real_2_lst, real_4_lst, bt_8_o, bt_9_o)[0])
    t_9_o = float(calFunByFsolve(real_2_lst, real_4_lst, bt_8_o, bt_9_o)[1])
    
    # 接下来由错误条件推导出错误的地面热辐射 从而得到错误的Tv Ts
    err_1_lst = para_err[0]
    err_2_lst = para_err[1]
    err_3_lst = para_err[2]
    err_4_lst = para_err[3]
    
    # 推导出错误的地面热辐射
    err_t_n_8 = err_1_lst[0] * t_8_n + err_1_lst[1] * (t_8_n - t_9_n) + err_1_lst[2] * (t_8_n - t_9_n)**2 + err_1_lst[3]
    err_t_o_8 = err_2_lst[0] * t_8_o + err_2_lst[1] * (t_8_o - t_9_o) + err_2_lst[2] * (t_8_o - t_9_o)**2 + err_2_lst[3]
    err_t_n_9 = err_3_lst[0] * t_8_n + err_3_lst[1] * (t_8_n - t_9_n) + err_3_lst[2] * (t_8_n - t_9_n)**2 + err_3_lst[3]
    err_t_o_9 = err_4_lst[0] * t_8_o + err_4_lst[1] * (t_8_o - t_9_o) + err_4_lst[2] * (t_8_o - t_9_o)**2 + err_4_lst[3]
    
    # 普朗克转换为辐亮度
    err_8_ang_1 = c1 / (np.pi * pow(band1, 5) * (np.exp(c2 / band1 / err_t_n_8) - 1))
    err_8_ang_2 = c1 / (np.pi * pow(band1, 5) * (np.exp(c2 / band1 / err_t_o_8) - 1))
    err_9_ang_1 = c1 / (np.pi * pow(band2, 5) * (np.exp(c2 / band2 / err_t_n_9) - 1))
    err_9_ang_2 = c1 / (np.pi * pow(band2, 5) * (np.exp(c2 / band2 / err_t_o_9) - 1))
    
    # 解算出Ts Tv
    cbtpev_8_ang_1, cbtpev_8_ang_2, cbtpev_9_ang_1, cbtpev_9_ang_2, cbtpes_8_ang_1, cbtpes_8_ang_2, cbtpes_9_ang_1, cbtpes_9_ang_2 = \
        calCbtp(sita, lai, es_real, es_real, ev_real, ev_real)
    
    L_err_v_8 = (cbtpes_8_ang_1 * err_8_ang_2 - cbtpes_8_ang_2 * err_8_ang_1) / (cbtpev_8_ang_2 * cbtpes_8_ang_1 - cbtpes_8_ang_2 * cbtpev_8_ang_1)
    L_err_s_8 = (err_8_ang_1 - cbtpev_8_ang_1 * L_err_v_8) / cbtpes_8_ang_1

    L_err_v_9 = (cbtpes_9_ang_1 * err_9_ang_2 - cbtpes_9_ang_2 * err_9_ang_1) / (cbtpev_9_ang_2 * cbtpes_9_ang_1 - cbtpes_9_ang_2 * cbtpev_9_ang_1)
    L_err_s_9 = (err_9_ang_1 - cbtpev_9_ang_1 * L_err_v_9) / cbtpes_9_ang_1
    
    # 普朗克反算Tv Ts
    Tv_err_8 = c2 / band1 / np.log((c1 / (np.pi * pow(band1, 5) * L_err_v_8)) + 1)
    Ts_err_8 = c2 / band1 / np.log((c1 / (np.pi * pow(band1, 5) * L_err_s_8)) + 1)

    Tv_err_9 = c2 / band2 / np.log((c1 / (np.pi * pow(band2, 5) * L_err_v_9)) + 1)
    Ts_err_9 = c2 / band2 / np.log((c1 / (np.pi * pow(band2, 5) * L_err_s_9)) + 1)
    
    # 
    Tv_err = (Tv_err_8 + Tv_err_9) / 2
    Ts_err = (Ts_err_8 + Ts_err_9) / 2
    
    print(Tv_err)
    print(Ts_err)
    return Tv_err, Ts_err
    
# cbtp模型计算
def calCbtp(sita, LAI, es_8, es_9, ev_8, ev_9):
    # i0 拦截率
    i0_ang_1 = 1 - np.exp(-0.5 * LAI / np.cos(sita[0]))
    i0_ang_2 = 1 - np.exp(-0.5 * LAI / np.cos(sita[1]))

    # i0hemi 半球平均拦截率，经验公式
    i0hemi = 1 - np.exp(-0.825 * LAI)

    # start数值积分法计算TIR-P
    fun_ang_1_up = 0
    fun_ang_2_up = 0
    fun_ang_1_down = 0
    fun_ang_2_down = 0

    for j in np.arange(0, float(LAI), 0.01):
        fun_ang_1_up = fun_ang_1_up + np.exp(-0.5 * j / np.cos(sita[0])) * 0.5 / np.cos(sita[0]) * 0.5 * np.exp(-0.8 * pow(j, 0.9)) / (1 - np.exp(-0.5 * LAI / np.cos(sita[0]))) * 0.01
        fun_ang_2_up = fun_ang_2_up + np.exp(-0.5 * j / np.cos(sita[1])) * 0.5 / np.cos(sita[1]) * 0.5 * np.exp(-0.8 * pow(j, 0.9)) / (1 - np.exp(-0.5 * LAI / np.cos(sita[1]))) * 0.01
        fun_ang_1_down = fun_ang_1_down + np.exp(-0.5 * j / np.cos(sita[0])) * 0.5 / np.cos(sita[0]) * 0.5 * np.exp(-0.8 * pow(LAI - j, 0.9)) / (1 - np.exp(-0.5 * LAI / np.cos(sita[0]))) * 0.01
        fun_ang_2_down = fun_ang_2_down + np.exp(-0.5 * j / np.cos(sita[1])) * 0.5 / np.cos(sita[1]) * 0.5 * np.exp(-0.8 * pow(LAI - j, 0.9)) / (1 - np.exp(-0.5 * LAI / np.cos(sita[1]))) * 0.01
    # elup 向上逃逸概率 eldown向下逃逸概率
    elup_ang_1 = round(fun_ang_1_up, 4)
    elup_ang_2 = round(fun_ang_2_up, 4)
    eldown_ang_1 = round(fun_ang_1_down, 4)
    eldown_ang_2 = round(fun_ang_2_down, 4) 

    # p 再碰撞概率
    p_ang_1 = 1 - elup_ang_1 - eldown_ang_1
    p_ang_2 = 1 - elup_ang_2 - eldown_ang_2

    # 计算DCE123简化模型
    # c1 光子被叶片反射后碰到叶片的概率
    c1_8_ang_1 = (1 - ev_8) * p_ang_1                           
    c1_8_ang_2 = (1 - ev_8) * p_ang_2
    c1_9_ang_1 = (1 - ev_9) * p_ang_1
    c1_9_ang_2 = (1 - ev_9) * p_ang_2
    # c2 光子被土壤反射后碰撞叶片的概率
    c2_8_ang_1 = (1 - es_8) * i0hemi
    c2_8_ang_2 = (1 - es_8) * i0hemi
    c2_9_ang_1 = (1 - es_9) * i0hemi
    c2_9_ang_2 = (1 - es_9) * i0hemi
    # c3 光子被叶片反射后碰撞土壤的概率
    c3_8_ang_1 = (1 - ev_8) * eldown_ang_1
    c3_8_ang_2 = (1 - ev_8) * eldown_ang_2
    c3_9_ang_1 = (1 - ev_9) * eldown_ang_1
    c3_9_ang_2 = (1 - ev_9) * eldown_ang_2

    # CBT-P模型 植被和土壤有效发射率
    cbtpev_8_ang_1 = i0_ang_1 * ev_8 * (1 + c1_8_ang_1 + c1_8_ang_1 * c1_8_ang_1 + c3_8_ang_1 * c2_8_ang_1) + (1 - i0_ang_1) * ev_8 * (c2_8_ang_1 + c2_8_ang_1 * c1_8_ang_1)
    cbtpev_8_ang_2 = i0_ang_2 * ev_8 * (1 + c1_8_ang_2 + c1_8_ang_2 * c1_8_ang_2 + c3_8_ang_2 * c2_8_ang_2) + (1 - i0_ang_2) * ev_8 * (c2_8_ang_2 + c2_8_ang_2 * c1_8_ang_2)
    cbtpev_9_ang_1 = i0_ang_1 * ev_9 * (1 + c1_9_ang_1 + c1_9_ang_1 * c1_9_ang_1 + c3_9_ang_1 * c2_9_ang_1) + (1 - i0_ang_1) * ev_9 * (c2_9_ang_1 + c2_9_ang_1 * c1_9_ang_1)
    cbtpev_9_ang_2 = i0_ang_2 * ev_9 * (1 + c1_9_ang_2 + c1_9_ang_2 * c1_9_ang_2 + c3_9_ang_2 * c2_9_ang_2) + (1 - i0_ang_2) * ev_9 * (c2_9_ang_2 + c2_9_ang_2 * c1_9_ang_2)

    cbtpes_8_ang_1 = i0_ang_1 * es_8 * (c3_8_ang_1 + c1_8_ang_1 * c3_8_ang_1) + (1 - i0_ang_1) * es_8 * (1 + c2_8_ang_1 * c3_8_ang_1)
    cbtpes_8_ang_2 = i0_ang_2 * es_8 * (c3_8_ang_2 + c1_8_ang_2 * c3_8_ang_2) + (1 - i0_ang_2) * es_8 * (1 + c2_8_ang_2 * c3_8_ang_2)
    cbtpes_9_ang_1 = i0_ang_1 * es_9 * (c3_9_ang_1 + c1_9_ang_1 * c3_9_ang_1) + (1 - i0_ang_1) * es_9 * (1 + c2_9_ang_1 * c3_9_ang_1)
    cbtpes_9_ang_2 = i0_ang_2 * es_9 * (c3_9_ang_2 + c1_9_ang_2 * c3_9_ang_2) + (1 - i0_ang_2) * es_9 * (1 + c2_9_ang_2 * c3_9_ang_2)
    
    return cbtpev_8_ang_1, cbtpev_8_ang_2, cbtpev_9_ang_1, cbtpev_9_ang_2, cbtpes_8_ang_1, cbtpes_8_ang_2, cbtpes_9_ang_1, cbtpes_9_ang_2

# 求解非线性二元二次方程组
def calFunBySympy(fun_1_para_lst, fun_2_para_lst, BT_1, BT_2):
    a1 = fun_1_para_lst[0]
    b1 = fun_1_para_lst[1]
    c1 = fun_1_para_lst[2]
    d1 = fun_1_para_lst[3]
    
    a2 = fun_2_para_lst[0]
    b2 = fun_2_para_lst[1]
    c2 = fun_2_para_lst[2]
    d2 = fun_2_para_lst[3]
    x = sympy.Symbol('x', real=True)
    y = sympy.Symbol('y', real=True)
    res = sympy.solve([a1 * x + b1 * (x - y) + c1 * (x - y)**2 + d1 - BT_1, a2 * x + b2* (x - y) + c2 * (x - y)**2 + d2 - BT_2], [x, y])
    return res
    # print(res)

def calFunByFsolve(fun_1_para_lst, fun_2_para_lst, BT_1, BT_2):
    a1 = fun_1_para_lst[0]
    b1 = fun_1_para_lst[1]
    c1 = fun_1_para_lst[2]
    d1 = fun_1_para_lst[3]
    
    a2 = fun_2_para_lst[0]
    b2 = fun_2_para_lst[1]
    c2 = fun_2_para_lst[2]
    d2 = fun_2_para_lst[3]
    def func(i):
        x, y = i[0], i[1]
        return [a1 * x + b1 * (x - y) + c1 * (x - y)**2 + d1 - BT_1, a2 * x + b2* (x - y) + c2 * (x - y)**2 + d2 - BT_2]

    r = fsolve(func, [300, 300])
    # print(r)
    return r
    
def drawPic(nadir_angle_lst, rmse_tv_lst, rmse_ts_lst, lai):
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 黑体
    plt.rcParams['axes.unicode_minus'] = False    # 解决无法显示符号的问题
    # 坐标轴的刻度设置向内(in)或向外(out)
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    sns.set(style="white", font='SimHei')

    plt.xlabel("Nadir方向观测天顶角")
    plt.ylabel("组分温度RMSE")
    plt.grid()
    plt.plot(nadir_angle_lst, rmse_tv_lst, color="green", linestyle='--', label='植被')
    plt.plot(nadir_angle_lst, rmse_ts_lst, color="red", linestyle='--', label='土壤')
    plt.legend(loc='upper right')
    plt.savefig("./res/" + str(lai) + "-劈窗-随角度.png", format='png', dpi=300)
    plt.clf()


ev_real = 0.98
es_real = 0.94
nadir_angle = 0
oblique_angle = 55
Tv_real_lst = [300, 300]
Ts_real_lst = [305, 310]
# 0 55
para_real = [[1.0174, 0.6045, 0.2652, -3.7037], 
             [1.0223, 0.9217, 0.2249, -4.3483], 
             [1.0292, -0.0667, 0.4339, -6.7047], 
             [1.035, 0.3792, 0.3414, -7.5637]]
# 5 55
para_err_5 = [[1.0175, 0.6063, 0.2649, -3.7083], 
             [1.0223, 0.9217, 0.2249, -4.3483], 
             [1.0292, -0.0641, 0.4333, -6.7113], 
             [1.035, 0.3792, 0.3414, -7.5637]]
# 10 55
para_err_10 = [[1.0176, 0.6118, 0.2641, -3.7216], 
               [1.0223, 0.9217, 0.2249, -4.3483], 
               [1.0294, -0.0562, 0.4314, -6.7311], 
               [1.035, 0.3792, 0.3414, -7.5637]]
# 15 55
para_err_15 = [[1.0177, 0.6211, 0.2628, -3.7451], 
               [1.0223, 0.9217, 0.2249, -4.3483], 
               [1.0296, -0.043, 0.4283, -6.7649], 
               [1.035, 0.3792, 0.3414, -7.5637]]
# 20 55
para_err_20 = [[1.0179, 0.6346, 0.2608, -3.7782], 
               [1.0223, 0.9217, 0.2249, -4.3483], 
               [1.0298, -0.0237, 0.4237, -6.8122], 
               [1.035, 0.3792, 0.3414, -7.5637]]
# 25 55
para_err_25 = [[1.0182, 0.6526, 0.2581, -3.8209], 
               [1.0223, 0.9217, 0.2249, -4.3483], 
               [1.0302, 0.0019, 0.4177, -6.8736], 
               [1.035, 0.3792, 0.3414, -7.5637]]
# 30 55
para_err_30 = [[1.0186, 0.6759, 0.2548, -3.8747],
               [1.0223, 0.9217, 0.2249, -4.3483], 
               [1.0307, 0.0351, 0.41, -6.9495], 
               [1.035, 0.3792, 0.3414, -7.5637]]

para_err_lst = [para_err_5, para_err_10, para_err_15, para_err_20, para_err_25, para_err_30]

# calSplitWindow(lai, ev_real, es_real, nadir_angle, oblique_angle, Tv_real, Ts_real, para_real, para_err_5)
lai_05 = []
lai_1 = []
lai_2 = []
lai_4 = []
lai_6 = []
lai_lst = [0.5, 1.0, 2.0, 4.0, 6.0]
nadir_angle_lst = [5, 10, 15, 20, 25, 30]

rmse_tv_lst = []
rmse_ts_lst = []

# 每个LAI画张图
# for lai in lai_lst:
#     tv_standard = []
#     ts_standard = []
#     tv_err_lst = []
#     ts_err_lst = []
#     for para_err in para_err_lst:
#         for i in range(len(Tv_real_lst)):
#             Tv_real = Tv_real_lst[i]
#             Ts_real = Ts_real_lst[i]
#             tv_standard.append(Tv_real)
#             ts_standard.append(Ts_real)
#             Tv_err, Ts_err = calSplitWindow(lai, ev_real, es_real, nadir_angle, oblique_angle, Tv_real, Ts_real, para_real, para_err)
#             tv_err_lst.append(Tv_err)
#             ts_err_lst.append(Ts_err)
#         # 计算每个para_err的rmse
#         mse_tv = mean_squared_error(tv_standard, tv_err_lst)
#         rmse_tv = np.sqrt(mse_tv)
#         rmse_tv_lst.append(rmse_tv)

#         mse_ts = mean_squared_error(ts_standard, ts_err_lst)
#         rmse_ts = np.sqrt(mse_ts)
#         rmse_ts_lst.append(rmse_ts)
#     # 绘制
#     drawPic(nadir_angle_lst, rmse_tv_lst, rmse_ts_lst, lai)

# for lai in lai_lst:
#     for para_err in para_err_lst:
#         for i in range(len(Tv_real_lst)):
#             Tv_real = Tv_real_lst[i]
#             Ts_real = Ts_real_lst[i]
#             calSplitWindow(lai, ev_real, es_real, nadir_angle, oblique_angle, Tv_real, Ts_real, para_real, para_err)
# calSplitWindow(4.0, ev_real, es_real, nadir_angle, oblique_angle, 300, 305, para_real, para_err_5)
lai = 1.0
Tv_real = 300
Ts_real = 305
calSplitWindow(lai, ev_real, es_real, nadir_angle, oblique_angle, Tv_real, Ts_real, para_real, para_err_5)
calSplitWindow(lai, ev_real, es_real, nadir_angle, oblique_angle, Tv_real, Ts_real, para_real, para_err_10)
calSplitWindow(lai, ev_real, es_real, nadir_angle, oblique_angle, Tv_real, Ts_real, para_real, para_err_15)
calSplitWindow(lai, ev_real, es_real, nadir_angle, oblique_angle, Tv_real, Ts_real, para_real, para_err_20)
calSplitWindow(lai, ev_real, es_real, nadir_angle, oblique_angle, Tv_real, Ts_real, para_real, para_err_25)
calSplitWindow(lai, ev_real, es_real, nadir_angle, oblique_angle, Tv_real, Ts_real, para_real, para_err_30)