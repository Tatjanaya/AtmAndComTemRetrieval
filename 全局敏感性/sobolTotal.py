import numpy as np
import simlabCbtp
import simlabCep
from SALib.sample import saltelli
from SALib.analyze import sobol
import linecache

# S8 0°为例
def calTotalModel(paraLst):
    
    lai = paraLst[0]
    es = paraLst[1]
    ev = paraLst[2]
    Tv = paraLst[3]
    Ts = paraLst[4]
    wv = paraLst[5]
    
    # 引入必要的参数
    c1 = 3.7404e8
    c2 = 14387
    
    # S8
    band1 = 10.852
    
    # S8 0°
    tran = -0.005 * wv * wv - 0.083 * wv + 0.981
    lup = 0.060 * wv * wv + 0.640 * wv - 0.064
    ldown = 0.006 * lup * lup + 1.033 * lup + 0.001
    
    # T转L
    L_Tv = c1 / (np.pi * pow(band1, 5) * (np.exp(c2 / band1 / Tv) - 1))
    L_Ts = c1 / (np.pi * pow(band1, 5) * (np.exp(c2 / band1 / Ts) - 1))
    
    # 计算cep +++vza+++
    cep = simlabCep.calCep(lai, 0, ev, es)
    
    # 计算cbtp +++vza+++
    cbtpev, cbtpes = simlabCbtp.calCbtpModel(lai, es, ev, 0)
    
    # ltoa
    L_toa = tran * ((1 - cep) * ldown + cbtpes * L_Ts + cbtpev * L_Tv) + lup
    return L_toa

def getEffData(samFile):
    # 行数
    countLine = int(linecache.getline(samFile, 2))
    # 存储结果表
    resLst = []
    
    # 从第五行开始
    for i in range(10):
        nowLine = i + 5
        lst = linecache.getline(samFile, nowLine).split()
        # 判断 Tv <= Ts Tv + 20 > Ts
        if float(lst[3]) <= float(lst[4]) and float(lst[3]) + 20 > float(lst[4]):
            tempLst = []
            tempLst.append(float(lst[0]))
            tempLst.append(float(lst[1]))
            tempLst.append(float(lst[2]))
            tempLst.append(float(lst[3]))
            tempLst.append(float(lst[4]))
            tempLst.append(float(lst[5]))
            resLst.append(tempLst)
    
    # 转置一下
    resLst = np.array(resLst)
    resLst = resLst.T
    return resLst
    
samFile = r'D:/simlab/soft/data/total/total.sam'

problem = {
    'num_vars': 6,
    'names': ['lai', 'es', 'ev', 'Tv', 'Ts', 'wv'],
    'bounds': [[0.5, 6.0],
               [0.92, 0.96],
               [0.95, 0.99],
               [280, 310],
               [280, 310],
               [0.1, 2.5]]
}

param_values = getEffData(samFile)
# 保存出来看看
np.savetxt('./data/param_values.txt', param_values)
Y = calTotalModel(param_values)
print(param_values.shape, Y.shape)
Si = sobol.analyze(problem, Y, print_to_console=True)
print()