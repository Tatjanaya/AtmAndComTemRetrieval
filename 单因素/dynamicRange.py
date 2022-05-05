import numpy as np
import os

def calTvTsWvFirVal(lai, Tv_real, Ts_real, wv_real, ev_real, es_real, nadir_angle, oblique_angle, Tv_err, Ts_err, wv_err, tranRelaLst, upRelaLst, downRelaLst):
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
    cbtpev_8_ang_1, cbtpev_8_ang_2, cbtpev_9_ang_1, cbtpev_9_ang_2, cbtpes_8_ang_1, cbtpes_8_ang_2, cbtpes_9_ang_1, cbtpes_9_ang_2, \
        ecaoP4sail_8_ang_1, ecaoP4sail_8_ang_2, ecaoP4sail_9_ang_1, ecaoP4sail_9_ang_2 = \
            calCbtp(sita, lai, es_8, es_9, ev_8, ev_9)
    # 反推地面热辐射 -> 理论上这是我们从卫星影像上根据劈窗法得到的结果
    ground_8_ang_1 = cbtpev_8_ang_1 * L_v_real_8 + cbtpes_8_ang_1 * L_s_real_8
    ground_8_ang_2 = cbtpev_8_ang_2 * L_v_real_8 + cbtpes_8_ang_2 * L_s_real_8
    ground_9_ang_1 = cbtpev_9_ang_1 * L_v_real_9 + cbtpes_9_ang_1 * L_s_real_9
    ground_9_ang_2 = cbtpev_9_ang_2 * L_v_real_9 + cbtpes_9_ang_2 * L_s_real_9
    
    # 算出8波段0度lup和tran
    tran_8_n, lup_8_n = wv2tranAndLup(wv_real)
    
    # 算出其他大气参数
    tran_8_o, tran_9_n, tran_9_o = calOtherTran(tran_8_n)
    lup_8_o, lup_9_n, lup_9_o = calOtherUp(lup_8_n)
    down_8_n, down_8_o, down_9_n, down_9_o = calOtherDown(lup_8_n)
    
    # 解算出四个LToa
    l_toa_8_n = tran_8_n * (ground_8_ang_1 + (1 - ecaoP4sail_8_ang_1) * down_8_n) + lup_8_n
    l_toa_8_o = tran_8_o * (ground_8_ang_2 + (1 - ecaoP4sail_8_ang_2) * down_8_o) + lup_8_o
    l_toa_9_n = tran_9_n * (ground_9_ang_1 + (1 - ecaoP4sail_9_ang_1) * down_9_n) + lup_9_n
    l_toa_9_o = tran_9_o * (ground_9_ang_2 + (1 - ecaoP4sail_9_ang_2) * down_9_o) + lup_9_o
    
    # 在有了四个TOA之后 可以进行牛顿迭代解算
    # 目前将采用四个初值 Tv_err Ts_real wv_real -> tran_8_n lup_8_n
    # x是初值 Tv采用错误的值进行代替
    # 错误的大气参数
    tran_8_n_err, lup_8_n_err = wv2tranAndLup(wv_err)
    x = [Tv_err, Ts_err, tran_8_n_err, lup_8_n_err]
    num = 4
    ecaoP4sailLst = [ecaoP4sail_8_ang_1, ecaoP4sail_8_ang_2, ecaoP4sail_9_ang_1, ecaoP4sail_9_ang_2]
    cbtpevLst = [cbtpev_8_ang_1, cbtpev_8_ang_2, cbtpev_9_ang_1, cbtpev_9_ang_2]
    cbtpesLst = [cbtpes_8_ang_1, cbtpes_8_ang_2, cbtpes_9_ang_1, cbtpes_9_ang_2]
    L_bt_lst = [l_toa_8_n, l_toa_8_o, l_toa_9_n, l_toa_9_o]
    return Newton(x, num, ecaoP4sailLst, band1, band2, c1, c2, cbtpevLst, cbtpesLst, L_bt_lst, tranRelaLst, upRelaLst, downRelaLst)

# 下面计算的都仅限于0 55角度的计算
# 水汽与8波段0度 透过率 上行辐射转换
# wv单位 g/cm2
# -0.005x*x+-0.083x+0.981	0.060x*x+0.640x+-0.064
def wv2tranAndLup(wv):
    tran = -0.005 * wv * wv -0.083 * wv + 0.981
    lup = 0.060 * wv * wv + 0.640 * wv -0.064
    return tran, lup

# 算出其他透过率
# trans_0_S8 trans_55_S8	trans_0_S8 trans_0_S9		trans_0_S8 trans_55_S9
# 0.645x*x+0.429x+-0.076	0.789x*x+0.156x+0.071		1.928x*x+-1.184x+0.273
def calOtherTran(tran_8_n):
    tran_8_o = 0.645 * tran_8_n * tran_8_n + 0.429 * tran_8_n - 0.076
    tran_9_n = 0.789 * tran_8_n * tran_8_n + 0.156 * tran_8_n + 0.071
    tran_9_o = 1.928 * tran_8_n * tran_8_n - 1.184 * tran_8_n + 0.273
    return tran_8_o, tran_9_n, tran_9_o

# up_ang_1_S8 up_ang_2_S8			up_ang_1_S8 up_ang_1_S9			up_ang_1_S8 up_ang_2_S9
# -0.072x*x+1.679x+0.008	-0.075x*x+1.526x+0.014	-0.187x*x+2.298x+0.052
def calOtherUp(lup_8_n):
    lup_8_o = -0.072 * lup_8_n * lup_8_n + 1.679 * lup_8_n + 0.008
    lup_9_n = -0.075 * lup_8_n * lup_8_n + 1.526 * lup_8_n + 0.014
    lup_9_o = -0.187 * lup_8_n * lup_8_n + 2.298 * lup_8_n + 0.052
    return lup_8_o, lup_9_n, lup_9_o

# up_ang_1_S8 down_0_S8		up_ang_1_S8 down_55_S8		up_ang_1_S8 down_0_S9		up_ang_1_S8 down_55_S9
# 0.006x*x+1.033x+0.001		-0.055x*x+1.708x+0.012	-0.073x*x+1.622x+0.009	-0.177x*x+2.429x+0.043
def calOtherDown(lup_8_n):
    down_8_n = 0.006 * lup_8_n * lup_8_n + 1.033 * lup_8_n + 0.001
    down_8_o = -0.055 * lup_8_n * lup_8_n + 1.708 * lup_8_n + 0.012
    down_9_n = -0.073 * lup_8_n * lup_8_n + 1.622 * lup_8_n + 0.00
    down_9_o = -0.177 * lup_8_n * lup_8_n + 2.429 * lup_8_n + 0.043
    return down_8_n, down_8_o, down_9_n, down_9_o

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
    
    # cavity 反射后未被植被吸收的部分，可能前向漫反射，也可能后向漫反射
    cavity_8_ang_1 = 1 - p_ang_1 * (1 - ev_8)
    cavity_8_ang_2 = 1 - p_ang_2 * (1 - ev_8)
    cavity_9_ang_1 = 1 - p_ang_1 * (1 - ev_9)
    cavity_9_ang_2 = 1 - p_ang_2 * (1 - ev_9)

    # rc1 冠层前向漫反射率
    # rc2 冠层后向漫反射率
    rc1_8_ang_1 = (1 - ev_8) * eldown_ang_1 / cavity_8_ang_1
    rc1_8_ang_2 = (1 - ev_8) * eldown_ang_2 / cavity_8_ang_2
    rc1_9_ang_1 = (1 - ev_9) * eldown_ang_1 / cavity_9_ang_1
    rc1_9_ang_2 = (1 - ev_9) * eldown_ang_2 / cavity_9_ang_2

    rc2_8_ang_1 = (1 - ev_8) * elup_ang_1 / cavity_8_ang_1
    rc2_8_ang_2 = (1 - ev_8) * elup_ang_2 / cavity_8_ang_2
    rc2_9_ang_1 = (1 - ev_9) * elup_ang_1 / cavity_9_ang_1
    rc2_9_ang_2 = (1 - ev_9) * elup_ang_2 / cavity_9_ang_2

    # 计算e1 e2 e3 e4 e5
    mulsv_8_ang_1 = 1 - rc2_8_ang_1 * (1 - es_8) * i0hemi
    mulsv_8_ang_2 = 1 - rc2_8_ang_2 * (1 - es_8) * i0hemi
    mulsv_9_ang_1 = 1 - rc2_9_ang_1 * (1 - es_9) * i0hemi
    mulsv_9_ang_2 = 1 - rc2_9_ang_2 * (1 - es_9) * i0hemi        

    e1_8_ang_1 = i0_ang_1 * ev_8 / cavity_8_ang_1
    e1_8_ang_2 = i0_ang_2 * ev_8 / cavity_8_ang_2
    e1_9_ang_1 = i0_ang_1 * ev_9 / cavity_9_ang_1
    e1_9_ang_2 = i0_ang_2 * ev_9 / cavity_9_ang_2

    e2_8_ang_1 = (1 - i0_ang_1) * (1 - es_8) * i0hemi * ev_8 / cavity_8_ang_1 / mulsv_8_ang_1
    e2_8_ang_2 = (1 - i0_ang_2) * (1 - es_8) * i0hemi * ev_8 / cavity_8_ang_2 / mulsv_8_ang_2
    e2_9_ang_1 = (1 - i0_ang_1) * (1 - es_9) * i0hemi * ev_9 / cavity_9_ang_1 / mulsv_9_ang_1
    e2_9_ang_2 = (1 - i0_ang_2) * (1 - es_9) * i0hemi * ev_9 / cavity_9_ang_2 / mulsv_9_ang_2

    e3_8_ang_1 = i0_ang_1 * rc1_8_ang_1 * (1 - es_8) * i0hemi * ev_8 / cavity_8_ang_1 / mulsv_8_ang_1
    e3_8_ang_2 = i0_ang_2 * rc1_8_ang_2 * (1 - es_8) * i0hemi * ev_8 / cavity_8_ang_2 / mulsv_8_ang_2
    e3_9_ang_1 = i0_ang_1 * rc1_9_ang_1 * (1 - es_9) * i0hemi * ev_9 / cavity_9_ang_1 / mulsv_9_ang_1
    e3_9_ang_2 = i0_ang_2 * rc1_9_ang_2 * (1 - es_9) * i0hemi * ev_9 / cavity_9_ang_2 / mulsv_9_ang_2

    e4_8_ang_1 = (1 - i0_ang_1) * es_8 / mulsv_8_ang_1
    e4_8_ang_2 = (1 - i0_ang_2) * es_8 / mulsv_8_ang_2 
    e4_9_ang_1 = (1 - i0_ang_1) * es_9 / mulsv_9_ang_1
    e4_9_ang_2 = (1 - i0_ang_2) * es_9 / mulsv_9_ang_2 

    e5_8_ang_1 = i0_ang_1 * rc1_8_ang_1 * es_8 / mulsv_8_ang_1
    e5_8_ang_2 = i0_ang_2 * rc1_8_ang_2 * es_8 / mulsv_8_ang_2
    e5_9_ang_1 = i0_ang_1 * rc1_9_ang_1 * es_9 / mulsv_9_ang_1
    e5_9_ang_2 = i0_ang_2 * rc1_9_ang_2 * es_9 / mulsv_9_ang_2

    # ecaoP4sail 冠层方向发射率
    ecaoP4sail_8_ang_1 = e1_8_ang_1 + e2_8_ang_1 + e3_8_ang_1 + e4_8_ang_1 + e5_8_ang_1
    ecaoP4sail_8_ang_2 = e1_8_ang_2 + e2_8_ang_2 + e3_8_ang_2 + e4_8_ang_2 + e5_8_ang_2
    ecaoP4sail_9_ang_1 = e1_9_ang_1 + e2_9_ang_1 + e3_9_ang_1 + e4_9_ang_1 + e5_9_ang_1
    ecaoP4sail_9_ang_2 = e1_9_ang_2 + e2_9_ang_2 + e3_9_ang_2 + e4_9_ang_2 + e5_9_ang_2

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
    
    return cbtpev_8_ang_1, cbtpev_8_ang_2, cbtpev_9_ang_1, cbtpev_9_ang_2, cbtpes_8_ang_1, cbtpes_8_ang_2, cbtpes_9_ang_1, cbtpes_9_ang_2, \
        ecaoP4sail_8_ang_1, ecaoP4sail_8_ang_2, ecaoP4sail_9_ang_1, ecaoP4sail_9_ang_2

# 迭代计算
def Fun(x, num, ecaoP4sailLst, band1, band2, c1, c2, cbtpevLst, cbtpesLst, L_bt_lst, tranRelaLst, upRelaLst, downRelaLst): 
    i = num
    f = np.zeros((i),dtype=float)
    f[0] = x[2] * ((1 - ecaoP4sailLst[0]) * (downRelaLst[0][0] * x[3] * x[3] + downRelaLst[0][1] * x[3] + downRelaLst[0][2]) + (cbtpevLst[0] * c1 / (np.pi * pow(band1, 5) * (np.exp(c2 / band1 / x[0]) - 1))) + (cbtpesLst[0] * c1 / (np.pi * pow(band1, 5) * (np.exp(c2 / band1 / x[1]) - 1)))) + x[3] - L_bt_lst[0]
    # 8 55
    f[1] = (tranRelaLst[0][0] * x[2] * x[2] + tranRelaLst[0][1] * x[2] + tranRelaLst[0][2]) * ((1 - ecaoP4sailLst[1]) * (downRelaLst[1][0] * x[3] * x[3] + downRelaLst[1][1] * x[3] + downRelaLst[1][2]) + (cbtpevLst[1] * c1 / (np.pi * pow(band1, 5) * (np.exp(c2 / band1 / x[0]) - 1))) + (cbtpesLst[1] * c1 / (np.pi * pow(band1, 5) * (np.exp(c2 / band1 / x[1]) - 1)))) + (upRelaLst[0][0] * x[3] * x[3] + upRelaLst[0][1] * x[3] + upRelaLst[0][2]) - L_bt_lst[1]
    # 9 0 
    f[2] = (tranRelaLst[1][0] * x[2] * x[2] + tranRelaLst[1][1] * x[2] + tranRelaLst[1][2]) * ((1 - ecaoP4sailLst[2]) * (downRelaLst[2][0] * x[3] * x[3] + downRelaLst[2][1] * x[3] + downRelaLst[2][2]) + (cbtpevLst[2] * c1 / (np.pi * pow(band2, 5) * (np.exp(c2 / band2 / x[0]) - 1))) + (cbtpesLst[2] * c1 / (np.pi * pow(band2, 5) * (np.exp(c2 / band2 / x[1]) - 1)))) + (upRelaLst[1][0] * x[3] * x[3] + upRelaLst[1][1] * x[3] + upRelaLst[1][2]) - L_bt_lst[2]
    # 9 55
    f[3] = (tranRelaLst[2][0] * x[2] * x[2] + tranRelaLst[2][1] * x[2] + tranRelaLst[2][2]) * ((1 - ecaoP4sailLst[3]) * (downRelaLst[3][0] * x[3] * x[3] + downRelaLst[3][1] * x[3] + downRelaLst[3][2]) + (cbtpevLst[3] * c1 / (np.pi * pow(band2, 5) * (np.exp(c2 / band2 / x[0]) - 1))) + (cbtpesLst[3] * c1 / (np.pi * pow(band2, 5) * (np.exp(c2 / band2 / x[1]) - 1)))) + (upRelaLst[2][0] * x[3] * x[3] + upRelaLst[2][1] * x[3] + upRelaLst[2][2]) - L_bt_lst[3]
    return f

def dfun(x, num, ecaoP4sailLst, band1, band2, c1, c2, cbtpevLst, cbtpesLst, L_bt_lst, tranRelaLst, upRelaLst, downRelaLst):                         #计算雅可比矩阵的逆矩阵
    df = np.zeros((num, num), dtype=float)
    dx = 0.00001                           # 
    x1 = np.copy(x)
    for i in range(0, num):              # 求导数，i是列，j是行
        for j in range(0, num):
            x1 = np.copy(x)
            x1[j] = x1[j] + dx           #x+dx
            df[i,j] = (Fun(x1, num, ecaoP4sailLst, band1, band2, c1, c2, cbtpevLst, cbtpesLst, L_bt_lst, tranRelaLst, upRelaLst, downRelaLst)[i] - Fun(x, num, ecaoP4sailLst, band1, band2, c1, c2, cbtpevLst, cbtpesLst, L_bt_lst, tranRelaLst, upRelaLst, downRelaLst)[i]) / dx   #f(x+dx)-f（x）/dx
    df_1 = np.linalg.inv(df)                              #计算逆矩阵
    return df_1

def Newton(x, num, ecaoP4sailLst, band1, band2, c1, c2, cbtpevLst, cbtpesLst, L_bt_lst, tranRelaLst, upRelaLst, downRelaLst):
    x1 = np.copy(x)
    i = 0
    delta = np.copy(x)
    try:
        while(np.sum(abs(delta)) > 1.e-8 and i < 20):  #控制循环次数
            x1 = x - np.dot(dfun(x, num, ecaoP4sailLst, band1, band2, c1, c2, cbtpevLst, cbtpesLst, L_bt_lst, tranRelaLst, upRelaLst, downRelaLst), Fun(x, num, ecaoP4sailLst, band1, band2, c1, c2, cbtpevLst, cbtpesLst, L_bt_lst, tranRelaLst, upRelaLst, downRelaLst))  #公式
            delta = x1 - x                     #比较x的变化
            x = x1
            i = i + 1
            # print(x)
        return x
    except np.linalg.LinAlgError as err:
        pass
    return [np.nan, np.nan, np.nan, np.nan]

lai_lst = [0.5, 1.0, 1.5, 2.0, 4.0, 6.0]
Tv_real = 300
Ts_real = 310
wv_real = 2
ev_real = 0.98
es_real = 0.94
nadir_angle = 0
oblique_angle = 55
tranRelaLst = [[0.645, 0.429, -0.076], 
               [0.789, 0.156, 0.071], 
               [1.928, -1.184, 0.273]]
upRelaLst = [[-0.072, 1.679, 0.008], 
             [-0.075, 1.526, 0.014], 
             [-0.187, 2.298, 0.052]]
downRelaLst = [[0.006, 1.033, 0.001], 
               [-0.055, 1.708, 0.012], 
               [-0.073, 1.622, 0.009], 
               [-0.177, 2.429, 0.043]]
f = open('./res/resDyRan.txt','w')
for lai in lai_lst:
    for Tv_err in range(280, 331):
        for Ts_err in range(Tv_err, Tv_err + 21):
            for wv_err in np.arange(0, 4, 0.1):
                x = calTvTsWvFirVal(lai, Tv_real, Ts_real, wv_real, ev_real, es_real, nadir_angle, oblique_angle, Tv_err, Ts_err, wv_err, tranRelaLst, upRelaLst, downRelaLst)
                if abs(x[0] - Tv_real) < 0.5 and abs(x[1] - Ts_real) < 0.5 and abs(x[2] - 0.795) < 0.01 and abs(x[3] - 1.456) < 0.02:
                    f.write('LAI: ' + str(lai) + ' Tv_err: ' + str(Tv_err) + ' Ts_err: ' + str(Ts_err) + ' wv_err: ' + str(wv_err))
                    f.write('\n')
                    
f.close()