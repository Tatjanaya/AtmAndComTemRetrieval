import numpy as np
import linecache
import os
import simlabCbtp

'''
分析Tv Ts初值全局敏感性
'''
def readSam(fileLoc, saveLoc):
    # 行数
    countLine = int(linecache.getline(fileLoc, 2))
    
    # 结果列表
    TvLst = []
    TsLst = []
    
    # 从第五行开始
    for i in range(countLine):
        nowLine = i + 5
        lst = linecache.getline(fileLoc, nowLine).split()
        Tv, Ts = calTvTs(float(lst[0]), float(lst[1]), float(lst[2]), float(lst[3]), float(lst[4]), float(lst[5]), float(lst[6]))
        TvLst.append(Tv)
        TsLst.append(Ts)
    
    # 写入
    f = open(saveLoc, 'w')
    f.write('2\n')
    f.write('cbtpev\n')
    f.write('cbtpes\n')
    f.write('time = no\n')
    f.write(str(countLine) + '\n')
    for i in range(countLine):
        f.write(str(TvLst[i]) + '\t' + str(TsLst[i]) + '\n')

def calTvTs(lai, es, ev, L_ground_8_nadir, L_ground_8_oblique, L_ground_9_nadir, L_ground_9_oblique):
    # 引入必要的参数
    c1 = 3.7404e8
    c2 = 14387
    
    # 波段
    band1 = 10.852
    band2 = 12
    
    # vza
    vza_nadir = 0
    vza_oblique = 55
    
    cbtpev_nadir, cbtpes_nadir = simlabCbtp.calCbtpModel(lai, es, ev, vza_nadir)
    cbtpev_oblique, cbtpes_oblique = simlabCbtp.calCbtpModel(lai, es, ev, vza_oblique)
    
    L_v_8 = (cbtpes_nadir * L_ground_8_oblique - cbtpes_oblique * L_ground_8_nadir) / (cbtpev_oblique * cbtpes_nadir - cbtpes_oblique * cbtpev_nadir)
    L_s_8 = (L_ground_8_nadir - cbtpev_nadir * L_v_8) / cbtpes_nadir

    L_v_9 = (cbtpes_nadir * L_ground_9_oblique - cbtpes_oblique * L_ground_9_nadir) / (cbtpev_oblique * cbtpes_nadir - cbtpes_oblique * cbtpev_nadir)
    L_s_9 = (L_ground_9_nadir - cbtpev_nadir * L_v_9) / cbtpes_nadir
    
    Tv_8 = c2 / band1 / np.log((c1 / (np.pi * pow(band1, 5) * L_v_8))+1)
    Ts_8 = c2 / band1 / np.log((c1 / (np.pi * pow(band1, 5) * L_s_8))+1)

    Tv_9 = c2 / band2 / np.log((c1 / (np.pi * pow(band2, 5) * L_v_9))+1)
    Ts_9 = c2 / band2 / np.log((c1 / (np.pi * pow(band2, 5) * L_s_9))+1)

    Tv = (Tv_8 + Tv_9) / 2
    Ts = (Ts_8 + Ts_9) / 2
    
    return Tv, Ts