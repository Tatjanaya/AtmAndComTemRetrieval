import numpy as np
import linecache
import simlabCbtp
import simlabCep

def calTotalModel(lai, es_8, ev_8, es_9, ev_9, Tv_fir, Ts_fir, wv_fir, L_toa_8_nadir, L_toa_8_oblique, L_toa_9_nadir, L_toa_9_oblique):
    # 引入必要的参数
    c1 = 3.7404e8
    c2 = 14387
    
    # 波段
    band1 = 10.852
    band2 = 12
    
    # vza
    vza_nadir = 0
    vza_oblique = 55
    
    # T转L
    
    
    # 水汽转透过率和上行
    # -0.005x*x+-0.083x+0.981	0.060x*x+0.640x+-0.064
    tran_fir = -0.005 * wv_fir * wv_fir - 0.083 * wv_fir + 0.981
    lup_fir = 0.060 * wv_fir * wv_fir + 0.640 * wv_fir - 0.064
    
    x = [Tv_fir, Ts_fir, tran_fir, lup_fir]
    
    # calCep(lai, vza, ev, es)
    ecaoP4sail_8_nadir = simlabCep.calCep(lai, vza_nadir, ev_8, es_8)
    ecaoP4sail_8_oblique = simlabCep.calCep(lai, vza_oblique, ev_8, es_8)
    ecaoP4sail_9_nadir = simlabCep.calCep(lai, vza_nadir, ev_9, es_9)
    ecaoP4sail_9_oblique = simlabCep.calCep(lai, vza_oblique, ev_9, es_9)
    
    # calCbtpModel(lai, es, ev, vza)
    cbtpev_8_nadir, cbtpes_8_nadir = simlabCbtp.calCbtpModel(lai, es_8, ev_8, 0)
    cbtpev_8_oblique, cbtpes_8_oblique = simlabCbtp.calCbtpModel(lai, es_8, ev_8, 55)
    cbtpev_9_nadir, cbtpes_9_nadir = simlabCbtp.calCbtpModel(lai, es_9, ev_9, 0)
    cbtpev_9_oblique, cbtpes_9_oblique = simlabCbtp.calCbtpModel(lai, es_9, ev_9, 55)
    
    tran_8_n = tran_fir
    lup_8_n = lup_fir
    
    # 算出其他大气参数
    tran_8_o, tran_9_n, tran_9_o = calOtherTran(tran_8_n)
    lup_8_o, lup_9_n, lup_9_o = calOtherUp(lup_8_n)
    down_8_n, down_8_o, down_9_n, down_9_o = calOtherDown(lup_8_n)
    
    ecaoP4sailLst = [ecaoP4sail_8_nadir, ecaoP4sail_8_oblique, ecaoP4sail_9_nadir, ecaoP4sail_9_oblique]
    cbtpevLst = [cbtpev_8_nadir, cbtpev_8_oblique, cbtpev_9_nadir, cbtpev_9_oblique]
    cbtpesLst = [cbtpes_8_nadir, cbtpes_8_oblique, cbtpes_9_nadir, cbtpes_9_oblique]
    L_bt_lst = [L_toa_8_nadir, L_toa_8_oblique, L_toa_9_nadir, L_toa_9_oblique]
    tranRelaLst = [tran_8_n, tran_8_o, tran_9_n, tran_9_o]
    upRelaLst = [lup_8_n, lup_8_o, lup_9_n, lup_9_o]
    downRelaLst = [down_8_n, down_8_o, down_9_n, down_9_o]
    
    return Newton(x, 4, ecaoP4sailLst, band1, band2, c1, c2, cbtpevLst, cbtpesLst, L_bt_lst, tranRelaLst, upRelaLst, downRelaLst)
    
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
            print(x)
        return x
    except np.linalg.LinAlgError as err:
        pass
    return [np.nan, np.nan, np.nan, np.nan]

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