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
def calWv(lai, vza, ev, es, niheFile, wv_lst, tv, ts, bandMark):
    toa_lst = []
    for wv in wv_lst:
        toa_lst.append(toaCepCbtp.calToa(lai, vza, ev, es, niheFile, wv, tv, ts, bandMark))
        
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 黑体
    plt.rcParams['axes.unicode_minus'] = False    # 解决无法显示符号的问题
    # 坐标轴的刻度设置向内(in)或向外(out)
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    sns.set(style="white", font='SimHei')
    
    plt.xlabel("大气水汽含量wv(g/cm2)")
    plt.ylabel("大气层顶温度K")
    plt.grid()
    plt.plot(wv_lst, toa_lst, color="navy", linestyle='solid')
    plt.savefig("./res/" + "单因素wv.png", format='png', dpi=300)
    
lai = 2
vza = 0
bandMark = 'S8'
ev = 0.98
es = 0.94
niheFile = './data/nihe_wv.txt'
wv = 2.0
tv = 300
ts = 320
wv_lst = np.arange(0.1, 4.01, 0.01)
calWv(lai, vza, ev, es, niheFile, wv_lst, tv, ts, bandMark)

wv = 0.5
print(toaCepCbtp.calToa(lai, vza, ev, es, niheFile, wv, tv, ts, bandMark))
wv = 1.5
print(toaCepCbtp.calToa(lai, vza, ev, es, niheFile, wv, tv, ts, bandMark))