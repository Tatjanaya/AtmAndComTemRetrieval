import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
import seaborn as sns 

# LAI = 0.5
# tv_5 = [299.9385526173246, 299.93531386738357]
# ts_5 = [305.0340702083756, 310.03411877448013]

# tv_10 = [299.7557312483843, 299.7444166125438]
# ts_10 = [305.1379323695733, 310.1387783612861]

# tv_15 = [299.6047594873387, 299.58624861212513]
# ts_15 = [305.223522392102, 310.2253127290643]

# tv_20 = [299.4062634578588, 299.3779730043855]
# ts_20 = [305.3355719861562, 310.3386024516793]

# tv_25 = [299.0065271601509, 298.9597645891752]
# ts_25 = [305.56066602384067, 310.5656718290775]

# tv_30 = [298.4753939721106, 298.40420046445274]
# ts_30 = [305.8580406410705, 310.8655322726255]

# LAI = 1.0
# tv_5 = [299.9728038288179, 299.9710813969465]
# ts_5 = [305.0384679950073, 310.0388266896212]

# tv_10 = [299.89081130401223, 299.88501298649066]
# ts_10 = [305.15542231371694, 310.15731091121927]

# tv_15 = [299.82453059437273, 299.8138016114399]
# ts_15 = [305.2498362009916, 310.2552047947406]

# tv_20 = [299.7375767653426, 299.7201690190933]
# ts_20 = [305.3733009817138, 310.3832997112298]

# tv_25 = [299.5604835496453, 299.532264828494]
# ts_25 = [305.62436985067876, 310.64016923492534]

# tv_30 = [299.32599975102255, 299.2831660316431]
# ts_30 = [305.95551676152047, 310.979315460891]

# LAI = 2.0
# tv_5 = [299.98926049458856, 299.98888585340626]
# ts_5 = [305.05064024903464, 310.0503282068914]

# tv_10 = [299.95673879196863, 299.95528957714725]
# ts_10 = [305.2047732276608, 310.2040111091849]

# tv_15 = [299.931460578966, 299.92873904705635]
# ts_15 = [305.32443409493, 310.325334666863]

# tv_20 = [299.89849081379214, 299.89400444689016]
# ts_20 = [305.4800303033954, 310.4833272663849]

# tv_25 = [299.8295803277709, 299.8222681108855]
# ts_25 = [305.8048165208113, 310.80945438544745]

# tv_30 = [299.7385041129836, 299.7274529151615]
# ts_30 = [306.2325685484486, 311.2389788730509]

# LAI = 4.0
# tv_5 = [299.9965475237279, 299.9965013531179]
# ts_5 = [305.1057572977112, 310.1029130242272]

# tv_10 = [299.98608373966846, 299.98590531336214]
# ts_10 = [305.42788450831324, 310.4177955911324]

# tv_15 = [299.9780910881248, 299.97775683246755]
# ts_15 = [305.673403040137, 310.65948517351416]

# tv_20 = [299.9677011826932, 299.9671499368168]
# ts_20 = [305.9912711116345, 310.9723340568888]

# tv_25 = [299.9457165649991, 299.9448184363325]
# ts_25 = [306.6619475896114, 311.6297761462645]

# tv_30 = [299.9166733381534, 299.9153180829237]
# ts_30 = [307.54266613822614, 312.49327918494225]

# LAI = 6.0
tv_5 = [299.9985324933609, 299.99852511794995]
ts_5 = [305.2519716097785, 310.2431972908447]

tv_10 = [299.9940828936292, 299.9940543855624]
ts_10 = [306.01930190976555, 310.98757490556767]

tv_15 = [299.990697193409, 299.99064373018143]
ts_15 = [306.6001198946319, 311.5527865289197]

tv_20 = [299.9862996641693, 299.9862114182289]
ts_20 = [307.34930370316306, 312.28123431464275]

tv_25 = [299.97696872897313, 299.9768249904086]
ts_25 = [308.92706142004647, 313.81442896608735]

tv_30 = [299.9646429395504, 299.9644261275986]
ts_30 = [310.9838171941356, 315.8145855140403]

tv_standard = [300, 300]
ts_standard = [305, 310]

nadir_angle_lst = [5, 10, 15, 20, 25, 30]
rmse_tv_lst = []
rmse_ts_lst = []

mse_tv = mean_squared_error(tv_standard, tv_5)
rmse_tv = np.sqrt(mse_tv)
rmse_tv_lst.append(rmse_tv)

mse_ts = mean_squared_error(ts_standard, ts_5)
rmse_ts = np.sqrt(mse_ts)
rmse_ts_lst.append(rmse_ts)

mse_tv = mean_squared_error(tv_standard, tv_10)
rmse_tv = np.sqrt(mse_tv)
rmse_tv_lst.append(rmse_tv)

mse_ts = mean_squared_error(ts_standard, ts_10)
rmse_ts = np.sqrt(mse_ts)
rmse_ts_lst.append(rmse_ts)

mse_tv = mean_squared_error(tv_standard, tv_15)
rmse_tv = np.sqrt(mse_tv)
rmse_tv_lst.append(rmse_tv)

mse_ts = mean_squared_error(ts_standard, ts_15)
rmse_ts = np.sqrt(mse_ts)
rmse_ts_lst.append(rmse_ts)

mse_tv = mean_squared_error(tv_standard, tv_20)
rmse_tv = np.sqrt(mse_tv)
rmse_tv_lst.append(rmse_tv)

mse_ts = mean_squared_error(ts_standard, ts_20)
rmse_ts = np.sqrt(mse_ts)
rmse_ts_lst.append(rmse_ts)

mse_tv = mean_squared_error(tv_standard, tv_25)
rmse_tv = np.sqrt(mse_tv)
rmse_tv_lst.append(rmse_tv)

mse_ts = mean_squared_error(ts_standard, ts_25)
rmse_ts = np.sqrt(mse_ts)
rmse_ts_lst.append(rmse_ts)

mse_tv = mean_squared_error(tv_standard, tv_30)
rmse_tv = np.sqrt(mse_tv)
rmse_tv_lst.append(rmse_tv)

mse_ts = mean_squared_error(ts_standard, ts_30)
rmse_ts = np.sqrt(mse_ts)
rmse_ts_lst.append(rmse_ts)

plt.rcParams['font.sans-serif'] = ['SimHei']  # 黑体
plt.rcParams['axes.unicode_minus'] = False    # 解决无法显示符号的问题
# 坐标轴的刻度设置向内(in)或向外(out)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
sns.set(style="white", font='SimHei')

plt.xlabel("Nadir方向观测天顶角")
plt.ylabel("组分温度RMSE")
plt.grid()
plt.plot(nadir_angle_lst, rmse_tv_lst, color="green", linestyle='--', label='植被')
plt.plot(nadir_angle_lst, rmse_ts_lst, color="red", linestyle='--', label='土壤')
plt.legend(loc='upper right')
plt.savefig("./res/" + "lai6.0-劈窗-随角度.png", format='png', dpi=300)
plt.clf()