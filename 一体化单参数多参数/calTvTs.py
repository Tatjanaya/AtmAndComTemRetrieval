from turtle import color
import numpy as np
import toaCepCbtp
from sklearn.metrics import mean_squared_error
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import matplotlib.ticker as ticker

'''
计算大气水汽含量变化 随 LAI 对观测值的影响
tv 300
ts 320
ev 0.95-0.99
es 0.94
band vza S8 0 画个图 S8 55 画个图
lai 2
'''
def calTv(lai, vza, ev, es, niheFile, wv, tv_lst, ts, bandMark):
    toa_lst = []
    for tv in tv_lst:
        toa_lst.append(toaCepCbtp.calToa(lai, vza, ev, es, niheFile, wv, tv, ts, bandMark))
        
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 黑体
    plt.rcParams['axes.unicode_minus'] = False    # 解决无法显示符号的问题
    # 坐标轴的刻度设置向内(in)或向外(out)
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    sns.set(style="white", font='SimHei')
    
    plt.xlabel("植被组分温度Tv/K")
    plt.ylabel("大气层顶温度/K")
    plt.grid()
    plt.plot(tv_lst, toa_lst, color="palegreen", linestyle='solid')
    plt.savefig("./res/" + "单因素tv" + str(lai) + ".png", format='png', dpi=300)
    
def calTs(lai, vza, ev, es, niheFile, wv, tv, ts_lst, bandMark):
    toa_lst = []
    for ts in ts_lst:
        toa_lst.append(toaCepCbtp.calToa(lai, vza, ev, es, niheFile, wv, tv, ts, bandMark))
        
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 黑体
    plt.rcParams['axes.unicode_minus'] = False    # 解决无法显示符号的问题
    # 坐标轴的刻度设置向内(in)或向外(out)
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    sns.set(style="white", font='SimHei')
    
    plt.xlabel("土壤组分温度Tv/K")
    plt.ylabel("大气层顶温度/K")
    plt.grid()
    plt.plot(ts_lst, toa_lst, color="palegreen", linestyle='solid')
    plt.savefig("./res/" + "单因素ts" + str(lai) + ".png", format='png', dpi=300)
    
lai = 4
vza = 0
bandMark = 'S8'
ev = 0.98
es = 0.94
niheFile = './data/nihe_wv.txt'
wv = 2.0
tv_lst = np.arange(290, 310.01, 0.01)
ts = 320
wv = 2.0
vza = 0
calTv(lai, vza, ev, es, niheFile, wv, tv_lst, ts, bandMark)