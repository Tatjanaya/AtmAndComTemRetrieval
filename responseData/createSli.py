import numpy as np
""" 生成SLSTR光谱响应函数Sli文件

"""
def createSli(S1txt, S2txt, S3txt, locTxt):
    f1 = open(S1txt, 'r')
    f2 = open(S2txt, 'r')
    f3 = open(S3txt, 'r')
    
    # 读取文件内容
    f1_content = f1.readlines()[4:]
    f2_content = f2.readlines()[4:]
    f3_content = f3.readlines()[4:]
    
    waveLen = range(520, 906)
    
    # print(f1_content[0].strip().split())
    
    # 最后存放到txt的内容文本 最后倒置
    res = np.zeros((4, len(waveLen)))
    
    # 先把waveLen写入第一列
    res[0] = waveLen
        
    # f1_content
    for i in range(len(f1_content)):
        line = f1_content[i].strip().split()
        # 确定应该放哪里
        loc = waveNumber2waveLength(float(line[0]))
        # 比如loc = 520 就是第一个
        res[1][loc - 520] = float(line[1])
        
    # f2_content
    for i in range(len(f2_content)):
        line = f2_content[i].strip().split()
        # 确定应该放哪里
        loc = waveNumber2waveLength(float(line[0]))
        # 比如loc = 520 就是第一个
        res[2][loc - 520] = float(line[1])
        
    # f3_content
    for i in range(len(f3_content)):
        line = f3_content[i].strip().split()
        # 确定应该放哪里
        loc = waveNumber2waveLength(float(line[0]))
        # 比如loc = 520 就是第一个
        res[3][loc - 520] = float(line[1])
        
    # 转置
    res = res.T
    
    # 写入txt中
    np.savetxt(locTxt, res, '%d %.4f %.4f %.4f')

def waveNumber2waveLength(waveNumber):
    return round(1e7 / waveNumber)
    
createSli('./data/rtcoef_sentinel3_1_slstr_srf_ch01.txt', './data/rtcoef_sentinel3_1_slstr_srf_ch02.txt', './data/rtcoef_sentinel3_1_slstr_srf_ch03.txt', 'slstr.txt')