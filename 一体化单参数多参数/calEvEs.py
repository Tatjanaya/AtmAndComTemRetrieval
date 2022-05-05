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
ev 0.95-0.99
es 0.94
band vza S8 0 画个图 S8 55 画个图
lai 2
'''
def calEv(lai, vza, ev_lst, es, niheFile, wv, tv, ts, bandMark):
    toa_lst = []
    for ev in ev_lst:
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
    plt.plot(ev_lst, toa_lst, color="green", linestyle='solid')
    plt.savefig("./res/" + "单因素ev-" + str(lai) + ".png", format='png', dpi=300)
    plt.clf()
    # plt.show()
    
def calEs(lai, vza, ev, es_lst, niheFile, wv, tv, ts, bandMark):
    toa_lst = []
    for es in es_lst:
        toa_lst.append(toaCepCbtp.calToa(lai, vza, ev, es, niheFile, wv, tv, ts, bandMark))
     
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 黑体
    plt.rcParams['axes.unicode_minus'] = False    # 解决无法显示符号的问题
    # 坐标轴的刻度设置向内(in)或向外(out)
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    sns.set(style="white", font='SimHei')
    plt.ticklabel_format(style='plain')
    plt.gca().xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    plt.xlabel("土壤发射率εs")
    plt.ylabel("大气层顶温度K")
    plt.grid()
    plt.plot(es_lst, toa_lst, color="brown", linestyle='solid')
    plt.yticks([307.11, 307.12, 307.13, 307.14, 307.15, 307.16, 307.17, 307.18], ['307.11', '307.12', '307.13', '307.14', '307.15', '307.16', '307.17', '307.18'])
    plt.savefig("./res/" + "单因素es-" + str(lai) + ".png", format='png', dpi=300)
    plt.clf()
    # plt.show()

lai = 2
vza = 0
bandMark = 'S8'
ev = 0.98
es = 0.94
niheFile = './data/nihe_wv.txt'
wv = 2.0
tv = 300
ts = 320
ev_lst = np.arange(0.95, 0.99, 0.001)
es_lst = np.arange(0.92, 0.96, 0.001)
calEv(0.5, vza, ev_lst, es, niheFile, wv, tv, ts, bandMark)
calEv(1.0, vza, ev_lst, es, niheFile, wv, tv, ts, bandMark)
calEv(2.0, vza, ev_lst, es, niheFile, wv, tv, ts, bandMark)
calEv(4.0, vza, ev_lst, es, niheFile, wv, tv, ts, bandMark)
# calEs(0.5, vza, ev, es_lst, niheFile, wv, tv, ts, bandMark)
# calEs(1.0, vza, ev, es_lst, niheFile, wv, tv, ts, bandMark)
# calEs(2.0, vza, ev, es_lst, niheFile, wv, tv, ts, bandMark)
# calEs(4.0, vza, ev, es_lst, niheFile, wv, tv, ts, bandMark)