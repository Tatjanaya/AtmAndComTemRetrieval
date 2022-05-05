from cProfile import label
from turtle import color
import numpy as np
import toaCepCbtp
from sklearn.metrics import mean_squared_error
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import matplotlib.ticker as ticker

'''
300 305 / 300 310 / 300 320
0 55
lai 0.5 - 4
wv 0.5 1.0 2.0
'''
def calMulti(lai_lst, ev, es, niheFile, wv, tv, ts):
    toa_lst_s8 = []
    toa_lst_s9 = []
    for lai in lai_lst:
        toa_lst_s8.append(toaCepCbtp.calToa(lai, 0, ev, es, niheFile, wv, tv, ts, 'S8') - toaCepCbtp.calToa(lai, 55, ev, es, niheFile, wv, tv, ts, 'S8'))
        toa_lst_s9.append(toaCepCbtp.calToa(lai, 0, ev, es, niheFile, wv, tv, ts, 'S9') - toaCepCbtp.calToa(lai, 55, ev, es, niheFile, wv, tv, ts, 'S9'))
        
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 黑体
    plt.rcParams['axes.unicode_minus'] = False    # 解决无法显示符号的问题
    # 坐标轴的刻度设置向内(in)或向外(out)
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    sns.set(style="whitegrid", font='SimHei')
    
    plt.xlabel("LAI")
    plt.ylabel("观测角0°和55°大气层顶亮温温差/K")
    plt.ylim(ymin=1.5, ymax=6.5)
    plt.plot(lai_lst, toa_lst_s8, color="blueviolet", linestyle='solid', label='S8')
    plt.plot(lai_lst, toa_lst_s9, color="dodgerblue", linestyle='solid', label='S9')
    plt.legend(loc='upper left')
    plt.savefig("./res-new/多参数/" + "多参数-" + str(wv) + "-" + str(ts - tv) + ".png", format='png', dpi=300)
    plt.clf()
    
ev = 0.98
es = 0.94
niheFile = './data/nihe_wv.txt'
tv = 300
ts = 310
wv_lst = [0.1, 0.5, 1.0, 2.0]
ts_lst = [305, 310, 320]
lai_lst = np.arange(0.5, 4.1, 0.1)
for wv in wv_lst:
    for ts in ts_lst:
        calMulti(lai_lst, ev, es, niheFile, wv, tv, ts)