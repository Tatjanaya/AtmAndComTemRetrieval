from cProfile import label
from turtle import color
import numpy as np
import toaCepCbtp
from sklearn.metrics import mean_squared_error
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import matplotlib.ticker as ticker

lai = 4.0
vza = 0
# vza = 55
ev = 0.98
es = 0.94
niheFile = './data/nihe_wv.txt'
wv = 2.0
tv = 300
ts = 320
print(toaCepCbtp.calToa(lai, vza, ev, es, niheFile, wv, tv, ts, 'S8'))
print(toaCepCbtp.calToa(lai, vza, ev, es, niheFile, wv, tv, ts, 'S9'))

plt.rcParams['font.sans-serif'] = ['SimHei']  # 黑体
plt.rcParams['axes.unicode_minus'] = False    # 解决无法显示符号的问题
# 坐标轴的刻度设置向内(in)或向外(out)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
sns.set(style="white", font='SimHei')

s8_lst = [toaCepCbtp.calToa(0.5, vza, ev, es, niheFile, wv, tv, ts, 'S8'), toaCepCbtp.calToa(1.0, vza, ev, es, niheFile, wv, tv, ts, 'S8'), \
            toaCepCbtp.calToa(2.0, vza, ev, es, niheFile, wv, tv, ts, 'S8'), toaCepCbtp.calToa(4.0, vza, ev, es, niheFile, wv, tv, ts, 'S8')]
s9_lst = [toaCepCbtp.calToa(0.5, vza, ev, es, niheFile, wv, tv, ts, 'S9'), toaCepCbtp.calToa(1.0, vza, ev, es, niheFile, wv, tv, ts, 'S9'), \
            toaCepCbtp.calToa(2.0, vza, ev, es, niheFile, wv, tv, ts, 'S9'), toaCepCbtp.calToa(4.0, vza, ev, es, niheFile, wv, tv, ts, 'S9')]

lai_lst = [0.5, 1.0, 2.0, 4.0]

# plt.xlabel("LAI")
# plt.ylabel("大气层顶温度/K")
# plt.grid()
# plt.plot(lai_lst, s8_lst, color="blueviolet", linestyle='solid', label='S8')
# plt.plot(lai_lst, s9_lst, color="dodgerblue", linestyle='solid', label='S9')
# plt.legend(loc='upper right')
# # plt.show()
# plt.savefig("./res/" + "单因素Band.png", format='png', dpi=300)