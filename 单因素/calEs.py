from cProfile import label
import numpy as np
from sklearn.metrics import mean_squared_error
import seaborn as sns
import os
import matplotlib.pyplot as plt

# 思路 0.92-0.96 0.01为间隔
# 先推导出正确的地面热辐射 然后根据错误的结果推导出组分温度 查看和正确的组分温度之间的差距
def calEsErrLstRes(lai, ev_real, es_real, nadir_angle, oblique_angle, Tv_real_lst, Ts_real_lst, es_err_lst):
    tv_rmse_res = []
    ts_rmse_res = []
    for es_err in es_err_lst:
        rmse_tv, rmse_ts = calEs(lai, ev_real, es_real, nadir_angle, oblique_angle, Tv_real_lst, Ts_real_lst, es_err)
        tv_rmse_res.append(rmse_tv)
        ts_rmse_res.append(rmse_ts)
        
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 黑体
    plt.rcParams['axes.unicode_minus'] = False    # 解决无法显示符号的问题
    # 坐标轴的刻度设置向内(in)或向外(out)
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    sns.set(style="white", font='SimHei')
    
    drawPic(es_err_lst, tv_rmse_res, ts_rmse_res, lai)

# 单个lai 单个es_err 一整个Tv Ts lst导出的情况
# 每个lai独立绘制一张图 该函数计算一个lai 一个es_err下的rmse（因为有多组lst）
def calEs(lai, ev_real, es_real, nadir_angle, oblique_angle, Tv_real_lst, Ts_real_lst, es_err):
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
    
    # 储备lst
    tv_err_lst = []
    ts_err_lst = []
    tv_standard = []
    ts_standard = []
    for i in range(len(Tv_real_lst)):
        tv_standard.append(Tv_real_lst[i])
        ts_standard.append(Ts_real_lst[i])
    
    for i in range(len(Tv_real_lst)):
        Tv_real = Tv_real_lst[i]
        Ts_real = Ts_real_lst[i]
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
        
        # 接着从错误的es推导出cbtp
        cbtpev_8_ang_1, cbtpev_8_ang_2, cbtpev_9_ang_1, cbtpev_9_ang_2, cbtpes_8_ang_1, cbtpes_8_ang_2, cbtpes_9_ang_1, cbtpes_9_ang_2 = \
            calCbtp(sita, lai, es_err, es_err, ev_real, ev_real)
        
        # 推出后反推计算ts tv
        L_err_v_8 = (cbtpes_8_ang_1 * ground_8_ang_2 - cbtpes_8_ang_2 * ground_8_ang_1) / (cbtpev_8_ang_2 * cbtpes_8_ang_1 - cbtpes_8_ang_2 * cbtpev_8_ang_1)
        L_err_s_8 = (ground_8_ang_1 - cbtpev_8_ang_1 * L_err_v_8) / cbtpes_8_ang_1

        L_err_v_9 = (cbtpes_9_ang_1 * ground_9_ang_2 - cbtpes_9_ang_2 * ground_9_ang_1) / (cbtpev_9_ang_2 * cbtpes_9_ang_1 - cbtpes_9_ang_2 * cbtpev_9_ang_1)
        L_err_s_9 = (ground_9_ang_1 - cbtpev_9_ang_1 * L_err_v_9) / cbtpes_9_ang_1
        
        # 普朗克
        Tv_err_8 = c2 / band1 / np.log((c1 / (np.pi * pow(band1, 5) * L_err_v_8)) + 1)
        Ts_err_8 = c2 / band1 / np.log((c1 / (np.pi * pow(band1, 5) * L_err_s_8)) + 1)

        Tv_err_9 = c2 / band2 / np.log((c1 / (np.pi * pow(band2, 5) * L_err_v_9)) + 1)
        Ts_err_9 = c2 / band2 / np.log((c1 / (np.pi * pow(band2, 5) * L_err_s_9)) + 1)
        
        # 
        Tv_err = (Tv_err_8 + Tv_err_9) / 2
        Ts_err = (Ts_err_8 + Ts_err_9) / 2
        
        tv_err_lst.append(Tv_err)
        ts_err_lst.append(Ts_err)
    
    mse_tv = mean_squared_error(tv_standard, tv_err_lst)
    rmse_tv = np.sqrt(mse_tv)
    
    mse_ts = mean_squared_error(ts_standard, ts_err_lst)
    rmse_ts = np.sqrt(mse_ts)
    
    # 返回单个err的rmse
    
    # 不算rmse 算温差
    tv_standard = np.array(tv_standard)
    ts_standard = np.array(ts_standard)
    tv_err_lst = np.array(tv_err_lst)
    ts_err_lst = np.array(ts_err_lst)
    
    chaju_tv_lst = tv_err_lst - tv_standard
    chaju_ts_lst = ts_err_lst - ts_standard
    # 计算平均
    cha_tv = sum(chaju_tv_lst) / len(chaju_tv_lst)
    cha_ts = sum(chaju_ts_lst) / len(chaju_ts_lst)
    
    print(tv_err_lst)
    print(ts_err_lst)
    
    return cha_tv, cha_ts
    # return rmse_tv, rmse_ts

# 绘制
def drawPic(es_err_lst, tv_rmse_lst, ts_rmse_lst, lai):
    plt.xlabel("给定的土壤组分发射率")
    plt.ylabel("组分温度差值")
    plt.grid()
    plt.plot(es_err_lst, tv_rmse_lst, color="green", linestyle='--', label='植被')
    plt.plot(es_err_lst, ts_rmse_lst, color="red", linestyle='--', label='土壤')
    plt.legend(loc='upper right')
    plt.savefig("./res/" + "es_err" + str(lai) + "组分温度差值.png", format='png', dpi=300)
    plt.clf()

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

# 0.5 1.0 2.0 4.0
ev_real = 0.98
es_real = 0.94
nadir_angle = 0
oblique_angle = 55
Tv_real = 300
Ts_real = 310
es_err = 0.92
lai = 2
# calEs(lai, ev_real, es_real, nadir_angle, oblique_angle, Tv_real, Ts_real, es_err)
Tv_real_lst = [300, 300]
Ts_real_lst = [305, 310]
es_err_lst = [0.92, 0.93, 0.94, 0.95, 0.96]
# calEsErrLstRes(0.5, ev_real, es_real, nadir_angle, oblique_angle, Tv_real_lst, Ts_real_lst, es_err_lst)
# calEsErrLstRes(1.0, ev_real, es_real, nadir_angle, oblique_angle, Tv_real_lst, Ts_real_lst, es_err_lst)
# calEsErrLstRes(2.0, ev_real, es_real, nadir_angle, oblique_angle, Tv_real_lst, Ts_real_lst, es_err_lst)
# calEsErrLstRes(4.0, ev_real, es_real, nadir_angle, oblique_angle, Tv_real_lst, Ts_real_lst, es_err_lst)
calEsErrLstRes(6.0, ev_real, es_real, nadir_angle, oblique_angle, Tv_real_lst, Ts_real_lst, es_err_lst)