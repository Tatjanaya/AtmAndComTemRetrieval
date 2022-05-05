import re
from symtable import Symbol
import numpy as np
from sklearn.metrics import mean_squared_error
import seaborn as sns
import os
import sympy
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# 上一步是Nadir角度发生偏差
# 这步是相同角度下水汽含量分组发生偏差
# 本质上是相同的 都是由劈窗参数不一致所引起的误差
# 拟合参数结果为(结构为a * x + b * (x - y) + c * (x - y)**2 + d): 
# 参数为: nadir: 0 oblique: 55 lowVapor: 0 highVapor: 2.5
# S8波段0°: [ 1.0174  0.6045  0.2652 -3.7037]
# S8波段55°: [ 1.0223  0.9217  0.2249 -4.3483]
# S9波段0°: [ 1.0292 -0.0667  0.4339 -6.7047]
# S9波段55°: [ 1.035   0.3792  0.3414 -7.5637]
# RMSE: 0.4768991302179044
# RMSE S8波段0°: 0.2962
# RMSE S8波段55°: 0.4097
# RMSE S9波段0°: 0.5136
# RMSE S9波段55°: 0.6247

# 拟合参数结果为(结构为a * x + b * (x - y) + c * (x - y)**2 + d): 
# 参数为: nadir: 0 oblique: 55 lowVapor: 2 highVapor: 3.5
# S8波段0°: [ 0.9676  2.3078 -0.0222  9.3004]
# S8波段55°: [ 9.04400e-01  3.22900e+00 -2.64000e-02  2.64505e+01]
# S9波段0°: [ 0.9948  2.0776 -0.025   1.7955]
# S9波段55°: [ 0.9196  3.1433 -0.0296 22.3866]
# RMSE: 0.8704631875457502
# RMSE S8波段0°: 0.5577
# RMSE S8波段55°: 0.8600
# RMSE S9波段0°: 0.9048
# RMSE S9波段55°: 1.0778

# 拟合参数结果为(结构为a * x + b * (x - y) + c * (x - y)**2 + d): 
# 参数为: nadir: 0 oblique: 55 lowVapor: 3 highVapor: 4.5
# S8波段0°: [ 0.9504  2.9501 -0.0456 12.8949]
# S8波段55°: [ 0.913   4.2432 -0.0844 21.127 ]
# S9波段0°: [ 0.9658  2.8903 -0.054   8.693 ]
# S9波段55°: [ 0.9194  4.248  -0.0903 19.4662]
# RMSE: 1.1702126785914992
# RMSE S8波段0°: 0.6549
# RMSE S8波段55°: 1.4186
# RMSE S9波段0°: 0.9117
# RMSE S9波段55°: 1.4849

# 拟合参数结果为(结构为a * x + b * (x - y) + c * (x - y)**2 + d): 
# 参数为: nadir: 0 oblique: 55 lowVapor: 4 highVapor: 6.5
# S8波段0°: [ 0.9499  3.1099  0.0171 11.7235]
# S8波段55°: [ 0.9681  4.5873 -0.0065  2.6306]
# S9波段0°: [0.9626 3.0925 0.0104 8.2628]
# S9波段55°: [ 0.9722  4.6139 -0.0121  1.5902]
# RMSE: 1.8318782390116408
# RMSE S8波段0°: 1.0504
# RMSE S8波段55°: 2.3354
# RMSE S9波段0°: 1.1735
# RMSE S9波段55°: 2.3428
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
    t_8_n = float(calFunBySympy(real_1_lst, real_3_lst, bt_8_n, bt_9_n)[0][0])
    t_9_n = float(calFunBySympy(real_1_lst, real_3_lst, bt_8_n, bt_9_n)[0][1])
    t_8_o = float(calFunBySympy(real_2_lst, real_4_lst, bt_8_o, bt_9_o)[0][0])
    t_9_o = float(calFunBySympy(real_2_lst, real_4_lst, bt_8_o, bt_9_o)[0][1])
    
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
    x = sympy.Symbol('x')
    y = sympy.Symbol('y')
    res = sympy.solve([a1 * x + b1 * (x - y) + c1 * (x - y)**2 + d1 - BT_1, a2 * x + b2* (x - y) + c2 * (x - y)**2 + d2 - BT_2], [x, y])
    return res
    # print(res)

lai = 4
ev_real = 0.98
es_real = 0.94
nadir_angle = 0
oblique_angle = 55
Tv_real = 300
Ts_real = 305
# 0 - 2.5
para_real = [[1.0174, 0.6045, 0.2652, -3.7037], 
             [1.0223, 0.9217, 0.2249, -4.3483], 
             [1.0292, -0.0667, 0.4339, -6.7047], 
             [1.035, 0.3792, 0.3414, -7.5637]]
# 2 - 3.5
para_err_235 = [[0.9676, 2.3078, -0.0222, 9.3004], 
             [0.9044, 3.2290, -0.0264, 26.4505], 
             [0.9948, 2.0776, -0.0250, 1.7955], 
             [0.9196, 3.1433, -0.0296, 22.3866]]
# 3 - 4.5
para_err_345 = [[0.9504, 2.9501, -0.0456, 12.8949], 
                [0.913, 4.2432, -0.0844, 21.127], 
                [0.9658, 2.8903, -0.054, 8.693], 
                [0.9194, 4.248, -0.0903, 19.4662]]
# 4 - 6.5
para_err_465 = [[0.9499, 3.1099, 0.0171, 11.7235], 
                [0.9681, 4.5873, -0.0065, 2.6306], 
                [0.9626, 3.0925, 0.0104, 8.2628], 
                [0.9722, 4.6139, -0.0121, 1.5902]]
calSplitWindow(lai, ev_real, es_real, nadir_angle, oblique_angle, Tv_real, Ts_real, para_real, para_err_235)
calSplitWindow(lai, ev_real, es_real, nadir_angle, oblique_angle, Tv_real, Ts_real, para_real, para_err_345)
calSplitWindow(lai, ev_real, es_real, nadir_angle, oblique_angle, Tv_real, Ts_real, para_real, para_err_465)
