from time import time
import numpy as np
import os
import random
from sklearn.metrics import mean_squared_error
import seaborn as sns
import matplotlib.pyplot as plt

# lai 从 0.5 - 6.0
# 正负0.5的随机误差
# es 0.94 ev 0.98
# Tv Ts [300 305] [300 310] [300 320]
# 做10000次
def calLaiLstRmse(laiLst, ev_real, es_real, nadir_angle, oblique_angle, Tv_real, Ts_real_lst, times):
    # 没文件夹创建文件夹
    dirs = './res/'
    if not os.path.exists(dirs):
        os.makedirs(dirs)
    
    # 每个值对应一个LAI的rmse
    rmseLst_tv = []
    rmseLst_ts = []
    
    for lai in laiLst:
        rmse_tv, rmse_ts = calLaiRmse(lai, ev_real, es_real, nadir_angle, oblique_angle, Tv_real, Ts_real_lst, times)
        rmseLst_tv.append(rmse_tv)
        rmseLst_ts.append(rmse_ts)
    
    # 绘图
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 黑体
    plt.rcParams['axes.unicode_minus'] = False    # 解决无法显示符号的问题
    # 坐标轴的刻度设置向内(in)或向外(out)
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    sns.set(style="white", font='SimHei')
    
    # x轴是laiLst y轴是rmse
    drawVegPic(laiLst, rmseLst_tv)
    drawSoilPic(laiLst, rmseLst_ts)

# 一个lai下的rmse计算
# 总共30000组 10000组305 10000组310 10000组320
# return 这个lai下的rmse
def calLaiRmse(lai, ev_real, es_real, nadir_angle, oblique_angle, Tv_real, Ts_real_lst, times):
    # 引入必要的参数
    c1 = 3.7404e8   
    c2 = 14387

    band1 = 10.852
    band2 = 12
    # sita
    sita = [nadir_angle / 180 * np.pi, oblique_angle / 180 * np.pi]
    
    # 标准tv ts
    tv_standard = [300] * times * 3
    ts_standard = []
    for Ts_real in Ts_real_lst:
        ts_standard += [Ts_real] * times
    
    # 得到的tv ts
    tv_err_lst = []
    ts_err_lst = []
    
    # 给LAI叠加±0.5的随机误差
    # 10000次循环 305
    # 先计算合理的地面热辐射
    for Ts_real in Ts_real_lst:
        cbtpev_8_ang_1, cbtpev_8_ang_2, cbtpev_9_ang_1, cbtpev_9_ang_2, cbtpes_8_ang_1, cbtpes_8_ang_2, cbtpes_9_ang_1, cbtpes_9_ang_2 = \
            calCbtp(sita, lai, es_real, es_real, ev_real, ev_real)
        # 真实Tv Ts 辐亮度
        L_v_real_8 = c1 / (np.pi * pow(band1, 5) * (np.exp(c2 / band1 / Tv_real) - 1))
        L_s_real_8 = c1 / (np.pi * pow(band1, 5) * (np.exp(c2 / band1 / Ts_real) - 1))
        L_v_real_9 = c1 / (np.pi * pow(band2, 5) * (np.exp(c2 / band2 / Tv_real) - 1))
        L_s_real_9 = c1 / (np.pi * pow(band2, 5) * (np.exp(c2 / band2 / Ts_real) - 1))
        # 反推地面热辐射 -> 理论上这是我们从卫星影像上根据劈窗法得到的结果
        ground_8_ang_1 = cbtpev_8_ang_1 * L_v_real_8 + cbtpes_8_ang_1 * L_s_real_8
        ground_8_ang_2 = cbtpev_8_ang_2 * L_v_real_8 + cbtpes_8_ang_2 * L_s_real_8
        ground_9_ang_1 = cbtpev_9_ang_1 * L_v_real_9 + cbtpes_9_ang_1 * L_s_real_9
        ground_9_ang_2 = cbtpev_9_ang_2 * L_v_real_9 + cbtpes_9_ang_2 * L_s_real_9
        
        # temp_count = 0
        for _ in range(times):
            lai_err = (lai - 0.5) + random.random()
            if lai_err <= 0.1:
                lai_err = 0.1
            # 根据错误lai计算cbtp
            cbtpev_8_ang_1, cbtpev_8_ang_2, cbtpev_9_ang_1, cbtpev_9_ang_2, cbtpes_8_ang_1, cbtpes_8_ang_2, cbtpes_9_ang_1, cbtpes_9_ang_2 = \
                calCbtp(sita, lai_err, es_real, es_real, ev_real, ev_real)
            # 根据错误的lai推导出来的错误的cbtp反推计算ts tv
            L_err_v_8 = (cbtpes_8_ang_1 * ground_8_ang_2 - cbtpes_8_ang_2 * ground_8_ang_1) / (cbtpev_8_ang_2 * cbtpes_8_ang_1 - cbtpes_8_ang_2 * cbtpev_8_ang_1)
            L_err_s_8 = (ground_8_ang_1 - cbtpev_8_ang_1 * L_err_v_8) / cbtpes_8_ang_1

            L_err_v_9 = (cbtpes_9_ang_1 * ground_9_ang_2 - cbtpes_9_ang_2 * ground_9_ang_1) / (cbtpev_9_ang_2 * cbtpes_9_ang_1 - cbtpes_9_ang_2 * cbtpev_9_ang_1)
            L_err_s_9 = (ground_9_ang_1 - cbtpev_9_ang_1 * L_err_v_9) / cbtpes_9_ang_1
            
            # 普朗克
            if (c1 / (np.pi * pow(band1, 5) * L_err_v_8)) + 1 <= 0:
                print('lai_err: ' + str(lai_err))
            Tv_err_8 = c2 / band1 / np.log((c1 / (np.pi * pow(band1, 5) * L_err_v_8)) + 1)
            Ts_err_8 = c2 / band1 / np.log((c1 / (np.pi * pow(band1, 5) * L_err_s_8)) + 1)

            Tv_err_9 = c2 / band2 / np.log((c1 / (np.pi * pow(band2, 5) * L_err_v_9)) + 1)
            Ts_err_9 = c2 / band2 / np.log((c1 / (np.pi * pow(band2, 5) * L_err_s_9)) + 1)
            
            # 
            Tv_err = (Tv_err_8 + Tv_err_9) / 2
            Ts_err = (Ts_err_8 + Ts_err_9) / 2
            
            # append
            tv_err_lst.append(Tv_err)
            ts_err_lst.append(Ts_err)
            # print('Finish: ' + str(temp_count))
            # temp_count += 1
    # 计算rmse
    mse_tv = mean_squared_error(tv_standard, tv_err_lst)
    rmse_tv = np.sqrt(mse_tv)
    
    mse_ts = mean_squared_error(ts_standard, ts_err_lst)
    rmse_ts = np.sqrt(mse_ts)
    
    return rmse_tv, rmse_ts
        
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

def drawVegPic(laiLst, tv_rmseLst):
    plt.xlabel("LAI")
    plt.ylabel("植被组分温度RMSE")
    plt.grid()
    plt.plot(laiLst, tv_rmseLst, color="black", linestyle='--')
    plt.savefig("./res/" + "植被组分温度RMSE.png", format='png', dpi=300)
    plt.clf()
    
def drawSoilPic(laiLst, ts_rmseLst):
    plt.xlabel("LAI")
    plt.ylabel("土壤组分温度RMSE")
    plt.grid()
    plt.plot(laiLst, ts_rmseLst, color="black", linestyle='--')
    plt.savefig("./res/" + "土壤组分温度RMSE.png", format='png', dpi=300)
    plt.clf()
        
Tv_real = 300
Ts_real_lst = [305, 310, 320]
nadir_angle = 0
oblique_angle = 55
ev_real = 0.98
es_real = 0.94
laiLst = np.arange(0.5, 6.1, 0.1)
times = 10000

calLaiLstRmse(laiLst, ev_real, es_real, nadir_angle, oblique_angle, Tv_real, Ts_real_lst, times)