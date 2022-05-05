import xlrd
import sys
import os
import numpy as np
from scipy import optimize
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
import linecache
""" 劈窗算法获取地面亮温

@splitWindow: 劈窗算法
excelFile: MODTRAN拟合结果excel表
lutChart: 廓线和水汽含量查找表
lowVapor: 水汽含量左边界
highVapor: 水汽含量右边界
atmFile: TIGR大气廓线数据
ang_1: 角度1
ang_2: 角度2

@extractAtmData: 从MODTRAN拟合结果中找到所需的一列数据
excelFile: res.excel文件位置
needLst: 符合条件的水汽含量相对位置列表
matchCol: 匹配的某一列

@drawPic: 绘制图像
BT: 亮温
T: 温度
ang_1: 角度1
ang_2: 角度2
lowWV: 最低水汽
highWV: 最高水汽
nowBand: 当前波段
nowAng: 当前角度
color: 选取颜色
"""
def splitWindow(excelFile, lutChart, atmFile, lowVapor, highVapor, ang_1, ang_2):
    # 如果不存在文件夹 则创建
    dirs = './result-new2/'
    if not os.path.exists(dirs):
        os.makedirs(dirs)
    
    if highVapor < lowVapor:
        return
    
    # 读取excel文件中水汽含量一列 TIGR编号是没有定位能力的 这是因为TIGR廓线选取避开气候恶劣区域（如暴雨等） 因此例编号100的廓线并不是第100条廓线
    content = xlrd.open_workbook(lutChart)
    sh = content.sheet_by_name('page_1')
    waterVapor = []
    for i in range(1, sh.nrows):
        waterVapor.append(float(sh.cell(i, 1).value))
    
    # 列表的长度是nrows - 1 1是表头
    lenWV = sh.nrows - 1
    
    # 存储符合条件的廓线位置 位置是相对第一个的偏差 而不是廓线本身的编号
    # 廓线编号不能定位chn绝对位置 例TIGR0100并不是出现的第100条廓线 因为去除了恶劣天气（如雷暴等）的廓线
    needLst = []
    
    # 如果下限高于上限 那么一个空列表没有计算意义
    if lowVapor > highVapor:
        print('Something wrong happens')
        sys.exit()
        
    for i in range(lenWV):
        if waterVapor[i] >= lowVapor and waterVapor[i] <= highVapor:
            needLst.append(i)
    
    # 空列表无意义
    if len(needLst) < 1:
        print('Something wrong happens')
        sys.exit()
    
    # 参与计算的廓线数量
    rnk = len(needLst)
    
    # 提取ang_1 - ang_2 透过率 上行辐射 下行辐射（不变）
    # 第一行为8波段 第二行为9波段
    atmDownRad_ang_1 = np.zeros((2, rnk), dtype=float)
    atmDownRad_ang_2 = np.zeros((2, rnk), dtype=float)
    
    atmUpRad_ang_1 = np.zeros((2, rnk), dtype=float)
    atmUpRad_ang_2 = np.zeros((2, rnk), dtype=float)
    
    atmUpTran_ang_1 = np.zeros((2, rnk), dtype=float)
    atmUpTran_ang_2 = np.zeros((2, rnk), dtype=float)
    
    atmDownRad_ang_1[0] = extractAtmData(excelFile, needLst, 'down_' + str(ang_1) + '_S8')
    atmDownRad_ang_1[1] = extractAtmData(excelFile, needLst, 'down_' + str(ang_1) + '_S9')
    
    atmDownRad_ang_2[0] = extractAtmData(excelFile, needLst, 'down_' + str(ang_2) + '_S8')
    atmDownRad_ang_2[1] = extractAtmData(excelFile, needLst, 'down_' + str(ang_2) + '_S9')
    
    atmUpRad_ang_1[0] = extractAtmData(excelFile, needLst, 'up_' + str(ang_1) + '_S8')
    atmUpRad_ang_1[1] = extractAtmData(excelFile, needLst, 'up_' + str(ang_1) + '_S9')
    
    atmUpRad_ang_2[0] = extractAtmData(excelFile, needLst, 'up_' + str(ang_2) + '_S8')
    atmUpRad_ang_2[1] = extractAtmData(excelFile, needLst, 'up_' + str(ang_2) + '_S9')
    
    atmUpTran_ang_1[0] = extractAtmData(excelFile, needLst, 'trans_' + str(ang_1) + '_S8')
    atmUpTran_ang_1[1] = extractAtmData(excelFile, needLst, 'trans_' + str(ang_1) + '_S9')
    
    atmUpTran_ang_2[0] = extractAtmData(excelFile, needLst, 'trans_' + str(ang_2) + '_S8')
    atmUpTran_ang_2[1] = extractAtmData(excelFile, needLst, 'trans_' + str(ang_2) + '_S9')
    
    # 引入必要的参数
    c1 = 3.7404e8   
    c2 = 14387

    band1 = 10.852
    band2 = 12

    nums = 6
    kind = 10
    
    # T 现在是个rnk * 6的矩阵
    T_ang_1 = np.zeros((rnk, 6), dtype=float)
    T_ang_2 = np.zeros((rnk, 6), dtype=float)
    
    loc = 0
    
    # T 需要提取廓线文件中底层温度
    for start in needLst:
        # 提取廓线第一行数据 气压1013mb
        line = linecache.getline(atmFile, start * 43 + 3).strip()
        
        # 真值
        real_T_ang_1 = float(line.split()[2])
        real_T_ang_2 = float(line.split()[2])
        
        if real_T_ang_1 >= 280:
            T_ang_1[loc][0] = real_T_ang_1 - 5
            T_ang_1[loc][1] = real_T_ang_1
            T_ang_1[loc][2] = real_T_ang_1 + 5
            T_ang_1[loc][3] = real_T_ang_1 + 10
            T_ang_1[loc][4] = real_T_ang_1 + 15
            T_ang_1[loc][5] = real_T_ang_1 + 20
        else:
            T_ang_1[loc][0] = real_T_ang_1 - 5
            T_ang_1[loc][1] = real_T_ang_1 - 5
            T_ang_1[loc][2] = real_T_ang_1
            T_ang_1[loc][3] = real_T_ang_1
            T_ang_1[loc][4] = real_T_ang_1 + 5
            T_ang_1[loc][5] = real_T_ang_1 + 5
        if real_T_ang_2 >= 280:
            T_ang_2[loc][0] = real_T_ang_2 - 5
            T_ang_2[loc][1] = real_T_ang_2
            T_ang_2[loc][2] = real_T_ang_2 + 5
            T_ang_2[loc][3] = real_T_ang_2 + 10
            T_ang_2[loc][4] = real_T_ang_2 + 15
            T_ang_2[loc][5] = real_T_ang_2 + 20
        else:
            T_ang_2[loc][0] = real_T_ang_2 - 5
            T_ang_2[loc][1] = real_T_ang_2 - 5
            T_ang_2[loc][2] = real_T_ang_2
            T_ang_2[loc][3] = real_T_ang_2
            T_ang_2[loc][4] = real_T_ang_2 + 5
            T_ang_2[loc][5] = real_T_ang_2 + 5
        
        loc += 1
    
    # 10种地物在11和12um处的发射率
    emm = np.array([[0.992692, 0.99328, 0.974837, 0.9653, 0.967658, 0.94871, 0.951261, 0.970933, 0.954495, 0.969078],
                [0.987771, 0.986645, 0.986192, 0.978892, 0.971629, 0.953049, 0.961569, 0.973487, 0.967363, 0.972328]])
    
    # 分地物 建
    B_sensor_s8_ang_1 = np.zeros((kind, rnk, nums), dtype=float)
    B_sensor_s9_ang_1 = np.zeros((kind, rnk, nums), dtype=float)

    B_sensor_s8_ang_2 = np.zeros((kind, rnk, nums), dtype=float)
    B_sensor_s9_ang_2 = np.zeros((kind, rnk, nums), dtype=float)

    # 地物辐射值
    B_sur8_ang_1 = np.zeros((kind, rnk, nums), dtype=float)
    B_sur8_ang_2 = np.zeros((kind, rnk, nums), dtype=float)
    B_sur9_ang_1 = np.zeros((kind, rnk, nums), dtype=float)
    B_sur9_ang_2 = np.zeros((kind, rnk, nums), dtype=float)
    
    # 引入大气纠正方程，得到不同地物发射率的传感器观测辐射值
    # rnk * 6
    L_s8_0 = c1 / (np.pi * pow(band1, 5) * (np.exp(c2 / band1 / T_ang_1) - 1))
    L_s8_55 = c1 / (np.pi * pow(band1, 5) * (np.exp(c2 / band1 / T_ang_2) - 1))
    L_s9_0 = c1 / (np.pi * pow(band2, 5) * (np.exp(c2 / band2 / T_ang_1) - 1))
    L_s9_55 = c1 / (np.pi * pow(band2, 5) * (np.exp(c2 / band2 / T_ang_2) - 1))

    for k in range(kind):
        for j in range(nums):
            for i in range(rnk):
                B_sensor_s8_ang_1[k][i][j] = atmUpTran_ang_1[0][i] * (L_s8_0[i][j] * emm[0][k] + atmDownRad_ang_1[0][i] * (1 - emm[0][k])) + atmUpRad_ang_1[0][i]
                B_sensor_s9_ang_1[k][i][j] = atmUpTran_ang_1[1][i] * (L_s9_0[i][j] * emm[1][k] + atmDownRad_ang_1[1][i] * (1 - emm[1][k])) + atmUpRad_ang_1[1][i]

                B_sensor_s8_ang_2[k][i][j] = atmUpTran_ang_2[0][i] * (L_s8_55[i][j] * emm[0][k] + atmDownRad_ang_2[0][i] * (1 - emm[0][k])) + atmUpRad_ang_2[0][i]
                B_sensor_s9_ang_2[k][i][j] = atmUpTran_ang_2[1][i] * (L_s9_55[i][j] * emm[1][k] + atmDownRad_ang_2[1][i] * (1 - emm[1][k])) + atmUpRad_ang_2[1][i]

                B_sur8_ang_1[k][i][j] = emm[0][k] * L_s8_0[i][j]
                B_sur8_ang_2[k][i][j] = emm[0][k] * L_s8_55[i][j]

                B_sur9_ang_1[k][i][j] = emm[1][k] * L_s9_0[i][j]
                B_sur9_ang_2[k][i][j] = emm[1][k] * L_s9_55[i][j]

    # 地表辐射亮温求出
    BT_s8_ang_1 = c2 / band1 / np.log((c1 / (np.pi * pow(band1, 5) * B_sur8_ang_1))+1)
    BT_s8_ang_2 = c2 / band1 / np.log((c1 / (np.pi * pow(band1, 5) * B_sur8_ang_2))+1)
    BT_s9_ang_1 = c2 / band2 / np.log((c1 / (np.pi * pow(band2, 5) * B_sur9_ang_1))+1)
    BT_s9_ang_2 = c2 / band2 / np.log((c1 / (np.pi * pow(band2, 5) * B_sur9_ang_2))+1)

    # 传感器观测到的8 9波段亮温值
    T_s8_ang_1 = c2 / band1 / np.log((c1 / (np.pi * pow(band1, 5) * B_sensor_s8_ang_1))+1)
    T_s8_ang_2 = c2 / band1 / np.log((c1 / (np.pi * pow(band1, 5) * B_sensor_s8_ang_2))+1)

    T_s9_ang_1 = c2 / band2 / np.log((c1 / (np.pi * pow(band2, 5) * B_sensor_s9_ang_1))+1)
    T_s9_ang_2 = c2 / band2 / np.log((c1 / (np.pi * pow(band2, 5) * B_sensor_s9_ang_2))+1)

    BT_s8_ang_1 = BT_s8_ang_1.ravel()
    BT_s8_ang_2 = BT_s8_ang_2.ravel()
    BT_s9_ang_1 = BT_s9_ang_1.ravel()
    BT_s9_ang_2 = BT_s9_ang_2.ravel()

    T_s8_ang_1 = T_s8_ang_1.ravel()
    T_s8_ang_2 = T_s8_ang_2.ravel()
    T_s9_ang_1 = T_s9_ang_1.ravel()
    T_s9_ang_2 = T_s9_ang_2.ravel()
    
    # 做拟合
    def func(x, y, p):
        a, b, c, d = p
        return a * x + b * (x - y) + c * (x - y)**2 + d

    def residuals(p, z, x, y):
        return z - func(x, y, p)
    
    plsq = optimize.leastsq(residuals, np.array([0, 0, 0, 0]), args=(BT_s8_ang_1, T_s8_ang_1, T_s9_ang_1))
    # plsq = optimize.leastsq(residuals, np.array([0, 0, 0, 0]), args=(BT_s8_ang_1, T_s8_ang_2, T_s9_ang_2))
    a_res, b_res, c_res, d_res = plsq[0] # 获得拟合结果

    plsq1 = optimize.leastsq(residuals, np.array([0, 0, 0, 0]), args=(BT_s8_ang_2, T_s8_ang_2, T_s9_ang_2))
    # plsq1 = optimize.leastsq(residuals, np.array([0, 0, 0, 0]), args=(BT_s8_ang_2, T_s8_ang_1, T_s9_ang_1))
    a_res1, b_res1, c_res1, d_res1 = plsq1[0] # 获得拟合结果

    plsq2 = optimize.leastsq(residuals, np.array([0, 0, 0, 0]), args=(BT_s9_ang_1, T_s8_ang_1, T_s9_ang_1))
    # plsq2 = optimize.leastsq(residuals, np.array([0, 0, 0, 0]), args=(BT_s9_ang_1, T_s8_ang_2, T_s9_ang_2))
    a_res2, b_res2, c_res2, d_res2 = plsq2[0] # 获得拟合结果

    plsq3 = optimize.leastsq(residuals, np.array([0, 0, 0, 0]), args=(BT_s9_ang_2, T_s8_ang_2, T_s9_ang_2))
    # plsq3 = optimize.leastsq(residuals, np.array([0, 0, 0, 0]), args=(BT_s9_ang_2, T_s8_ang_1, T_s9_ang_1))
    a_res3, b_res3, c_res3, d_res3 = plsq3[0] # 获得拟合结果
    
    print('----参数----')
    print(plsq[0])
    print(plsq1[0])
    print(plsq2[0])
    print(plsq3[0])
    print('------------')
    
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 黑体
    plt.rcParams['axes.unicode_minus'] = False    # 解决无法显示符号的问题
    # 坐标轴的刻度设置向内(in)或向外(out)
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    sns.set(style="whitegrid", font='SimHei')

    T_test_8_ang_1 = a_res * T_s8_ang_1 + b_res * (T_s8_ang_1 - T_s9_ang_1) + c_res * (T_s8_ang_1 - T_s9_ang_1)**2 + d_res
    T_test_8_ang_2 = a_res1 * T_s8_ang_2 + b_res1 * (T_s8_ang_2 - T_s9_ang_2) + c_res1 * (T_s8_ang_2 - T_s9_ang_2)**2 + d_res1
    T_test_9_ang_1 = a_res2 * T_s8_ang_1 + b_res2 * (T_s8_ang_1 - T_s9_ang_1) + c_res2 * (T_s8_ang_1 - T_s9_ang_1)**2 + d_res2
    T_test_9_ang_2 = a_res3 * T_s8_ang_2 + b_res3 * (T_s8_ang_2 - T_s9_ang_2) + c_res3 * (T_s8_ang_2 - T_s9_ang_2)**2 + d_res3

    T_total = np.hstack((T_test_8_ang_1, T_test_8_ang_2, T_test_9_ang_1, T_test_9_ang_2))
    BT_total = np.hstack((BT_s8_ang_1, BT_s8_ang_2, BT_s9_ang_1, BT_s9_ang_2))
    mse = mean_squared_error(BT_total, T_total)
    rmse = np.sqrt(mse)

    mse_8_ang_1 = mean_squared_error(BT_s8_ang_1, T_test_8_ang_1)
    rmse_8_ang_1 = np.sqrt(mse_8_ang_1)

    mse_8_ang_2 = mean_squared_error(BT_s8_ang_2, T_test_8_ang_2)
    rmse_8_ang_2 = np.sqrt(mse_8_ang_2)

    mse_9_ang_1 = mean_squared_error(BT_s9_ang_1, T_test_9_ang_1)
    rmse_9_ang_1 = np.sqrt(mse_9_ang_1)

    mse_9_ang_2 = mean_squared_error(BT_s9_ang_2, T_test_9_ang_2)
    rmse_9_ang_2 = np.sqrt(mse_9_ang_2)

    # 平均误差系数
    L_T_test_8_ang_1 = 3.7404e8 / (np.pi * pow(band1, 5) * (np.exp(14387 / band1 / np.array(T_test_8_ang_1)) - 1))
    L_T_test_8_ang_2 = 3.7404e8 / (np.pi * pow(band1, 5) * (np.exp(14387 / band1 / np.array(T_test_8_ang_2)) - 1))
    L_T_test_9_ang_1 = 3.7404e8 / (np.pi * pow(band2, 5) * (np.exp(14387 / band2 / np.array(T_test_9_ang_1)) - 1))
    L_T_test_9_ang_2 = 3.7404e8 / (np.pi * pow(band2, 5) * (np.exp(14387 / band2 / np.array(T_test_9_ang_2)) - 1))

    L_BT_s8_ang_1 = 3.7404e8 / (np.pi * pow(band1, 5) * (np.exp(14387 / band1 / np.array(BT_s8_ang_1)) - 1))
    L_BT_s8_ang_2 = 3.7404e8 / (np.pi * pow(band1, 5) * (np.exp(14387 / band1 / np.array(BT_s8_ang_2)) - 1))
    L_BT_s9_ang_1 = 3.7404e8 / (np.pi * pow(band2, 5) * (np.exp(14387 / band2 / np.array(BT_s9_ang_1)) - 1))
    L_BT_s9_ang_2 = 3.7404e8 / (np.pi * pow(band2, 5) * (np.exp(14387 / band2 / np.array(BT_s9_ang_2)) - 1))

    x1 = np.sum(np.abs((L_T_test_8_ang_1 - L_BT_s8_ang_1) / L_BT_s8_ang_1)) / len(L_BT_s8_ang_1)
    print("x1 = " + str(format(x1, '.4f')))
    x2 = np.sum(np.abs((L_T_test_8_ang_2 - L_BT_s8_ang_2) / L_BT_s8_ang_2)) / len(L_BT_s8_ang_2)
    print("x2 = " + str(format(x2, '.4f')))
    x3 = np.sum(np.abs((L_T_test_9_ang_1 - L_BT_s9_ang_1) / L_BT_s9_ang_1)) / len(L_BT_s9_ang_1)
    print("x3 = " + str(format(x3, '.4f')))
    x4 = np.sum(np.abs((L_T_test_9_ang_2 - L_BT_s9_ang_2) / L_BT_s9_ang_2)) / len(L_BT_s9_ang_2)
    print("x4 = " + str(format(x4, '.4f')))
    
    drawPic(BT_s8_ang_1, T_test_8_ang_1, ang_1, ang_2, lowVapor, highVapor, 'S8', ang_1, 'r')
    drawPic(BT_s8_ang_2, T_test_8_ang_2, ang_1, ang_2, lowVapor, highVapor, 'S8', ang_2, 'slateblue')
    drawPic(BT_s9_ang_1, T_test_9_ang_1, ang_1, ang_2, lowVapor, highVapor, 'S9', ang_1, 'violet')
    drawPic(BT_s9_ang_2, T_test_9_ang_2, ang_1, ang_2, lowVapor, highVapor, 'S9', ang_2, 'seagreen')
    
    # a * x + b * (x - y) + c * (x - y)**2 + d
    f = open("./result-new2/splitWindowResNew.txt", "a+")
    f.write("拟合参数结果为(结构为a * x + b * (x - y) + c * (x - y)**2 + d): " + "\n" \
                + "参数为: " + "nadir: " + str(ang_1) + " oblique: " + str(ang_2) + " lowVapor: " + str(lowVapor) + " highVapor: " + str(highVapor) + "\n" \
                    + "S8波段" + str(ang_1) + "°: " + str(formatTrans(plsq[0])) + "\n" \
                        + "S8波段" + str(ang_2) + "°: " + str(formatTrans(plsq1[0])) + "\n" \
                            + "S9波段" + str(ang_1) + "°: " + str(formatTrans(plsq2[0])) + "\n" \
                                + "S9波段" + str(ang_2) + "°: " + str(formatTrans(plsq3[0])) + "\n" \
                                    + "RMSE: " + str(rmse) + "\n" \
                                        + "RMSE S8波段" + str(ang_1) + "°: " + str(format(rmse_8_ang_1, '.4f')) + "\n" \
                                            + "RMSE S8波段" + str(ang_2) + "°: " + str(format(rmse_8_ang_2, '.4f')) + "\n" \
                                                + "RMSE S9波段" + str(ang_1) + "°: " + str(format(rmse_9_ang_1, '.4f')) + "\n" \
                                                    + "RMSE S9波段" + str(ang_2) + "°: " + str(format(rmse_9_ang_2, '.4f')) + "\n" \
                                                        + "平均相对偏差: " + "8_n(" + str(format(x1, '.4f')) + ") " + "8_o(" + str(format(x2, '.4f')) + ") " + "9_n(" + str(format(x3, '.4f')) + ") " + "9_o(" + str(format(x4, '.4f')) + ") " + "\n")

def drawPic(BT, T, ang_1, ang_2, lowWV, highWV, nowBand, nowAng, color):
    plt.title("水汽含量" + str(lowWV) + "-" + str(highWV) + "gm/cm2  " + "角度组合: " + str(ang_1) + "° " + str(ang_2) + "°")
    plt.xlabel("地表辐射真值(" + str(nowBand) + "波段" + str(nowAng) + "°" + ")/K")
    plt.ylabel("地表辐射模拟值(" + str(nowBand) + "波段" + str(nowAng) + "°" + ")/K")
    plt.scatter(BT, T, s=0.01, c=color, marker='o')
    plt.plot(np.linspace(240, 330, 60), np.linspace(240, 330, 60), color="lemonchiffon")
    plt.savefig("./result-new2/" + str(lowWV) + "-" + str(highWV) + "_" + str(nowBand) + "_" + str(nowAng) + "°.png", format='png', dpi=300)
    plt.clf()
    
def formatTrans(lst):
    for i in range(len(lst)):
        lst[i] = format(lst[i], '.4f')
    return lst
    
def extractAtmData(excelFile, needLst, matchCol):
    content = xlrd.open_workbook(excelFile)
    sh = content.sheet_by_name('page_1')
    
    res = []
    
    # 找到匹配的列
    tempCol = 0
    for i in range(sh.ncols):
        if matchCol == sh.cell(0, i).value:
            break
        tempCol += 1
    
    # 如果匹配不到 说明matchCol有问题 返回空列表
    if tempCol >= sh.ncols:
        return np.zeros((1, len(needLst)))
    
    # 找到指定位置的列后 根据needLst存储
    for need in needLst:
        res.append(sh.cell(need + 1, tempCol).value)
    
    return res