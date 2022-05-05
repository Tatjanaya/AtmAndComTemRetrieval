from turtle import color
import numpy as np
import toaCepCbtp
from sklearn.metrics import mean_squared_error
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import matplotlib.ticker as ticker

'''
计算组分发射率变化 随 LAI 对观测值的影响
tv 300
ts 320
ev 0.98
es 0.94
band vza S8 0 画个图 S8 55 画个图
lai 0.5-7
'''
def calLai(lai_lst, vza, ev, es, niheFile, wv, tv, ts, bandMark):
    toa_lst = []
    for lai in lai_lst:
        toa_lst.append(toaCepCbtp.calToa(lai, vza, ev, es, niheFile, wv, tv, ts, bandMark))
    
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 黑体
    plt.rcParams['axes.unicode_minus'] = False    # 解决无法显示符号的问题
    # 坐标轴的刻度设置向内(in)或向外(out)
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    sns.set(style="white", font='SimHei')
    
    plt.xlabel("植被发射率εv")
    plt.ylabel("大气层顶温度K")
    plt.grid()
    plt.plot(lai_lst, toa_lst, color="deepskyblue", linestyle='solid')
    plt.savefig("./res/" + "单因素lai.png", format='png', dpi=300)

lai_lst = np.arange(0.5, 7.1, 0.1)
vza = 0
bandMark = 'S8'
ev = 0.98
es = 0.94
niheFile = './data/nihe_wv.txt'
wv = 2.0
tv = 300
ts = 320
calLai(lai_lst, vza, ev, es, niheFile, wv, tv, ts, bandMark)