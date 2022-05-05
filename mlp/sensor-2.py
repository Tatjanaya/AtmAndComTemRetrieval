from sklearn import preprocessing
import numpy as np
from sklearn.neural_network import MLPRegressor
import toaCepCbtp
from sklearn.metrics import mean_squared_error
import random
import re
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
# 生成数据
# x_data (5000, 5) -> y_data (5000, 3)
# 待会也28分流
def genData(lai_lst, wv_lst, tv_lst, es_8, es_9, ev_8, ev_9, nadir, oblique, wvFile, niheFile):
    c1 = 3.7404e8   
    c2 = 14387

    band1 = 10.852
    band2 = 12
    # x_data -> 4 toa 1 lai
    x_data = np.zeros((5000, 5))
    y_data = np.zeros((5000, 3))
    # 水汽推导
    f = open(wvFile, 'r')
    lines = f.readlines()
    f.close()
    # 读取 waterVapor trans_0_S8
    tarTranStr = 'waterVapor trans_' + str(nadir) + '_S8'
    for i in range(len(lines)):
        if lines[i].strip() == tarTranStr:
            tran_8n_lst = re.findall('-?\d+\.?\d*', lines[i + 1])
            tran_8n_lst = [float(x) for x in tran_8n_lst]
            break
    tarTranStr = 'waterVapor up_' + str(nadir) + '_S8'
    for i in range(len(lines)):
        if lines[i].strip() == tarTranStr:
            up_8n_lst = re.findall('-?\d+\.?\d*', lines[i + 1])
            up_8n_lst = [float(x) for x in up_8n_lst]
            break
    # 读取关系式
    f = open(niheFile, 'r')
    lines = f.readlines()
    tran_9n_lst_str = 'trans_' + str(nadir) + '_S8 trans_' + str(nadir) + '_S9'
    tran_8o_lst_str = 'trans_' + str(nadir) + '_S8 trans_' + str(oblique) + '_S8'
    tran_9o_lst_str = 'trans_' + str(nadir) + '_S8 trans_' + str(oblique) + '_S9'
    
    up_9n_lst_str = 'up_' + str(nadir) + '_S8 up_' + str(nadir) + '_S9'
    up_8o_lst_str = 'up_' + str(nadir) + '_S8 up_' + str(oblique) + '_S8'
    up_9o_lst_str = 'up_' + str(nadir) + '_S8 up_' + str(oblique) + '_S9'
    
    down_8n_lst_str = 'up_' + str(nadir) + '_S8 down_' + str(nadir) + '_S8'
    down_9n_lst_str = 'up_' + str(nadir) + '_S8 down_' + str(nadir) + '_S9'
    down_8o_lst_str = 'up_' + str(nadir) + '_S8 down_' + str(oblique) + '_S8'
    down_9o_lst_str = 'up_' + str(nadir) + '_S8 down_' + str(oblique) + '_S9'
    for i in range(len(lines)):
        if lines[i].strip() == tran_9n_lst_str:
            tran_9n_lst = re.findall('-?\d+\.?\d*', lines[i + 1])
            tran_9n_lst = [float(x) for x in tran_9n_lst]
        if lines[i].strip() == tran_8o_lst_str:
            tran_8o_lst = re.findall('-?\d+\.?\d*', lines[i + 1])
            tran_8o_lst = [float(x) for x in tran_8o_lst]
        if lines[i].strip() == tran_9o_lst_str:
            tran_9o_lst = re.findall('-?\d+\.?\d*', lines[i + 1])
            tran_9o_lst = [float(x) for x in tran_9o_lst]
            
        if lines[i].strip() == up_9n_lst_str:
            up_9n_lst = re.findall('-?\d+\.?\d*', lines[i + 1])
            up_9n_lst = [float(x) for x in up_9n_lst]
        if lines[i].strip() == up_8o_lst_str:
            up_8o_lst = re.findall('-?\d+\.?\d*', lines[i + 1])
            up_8o_lst = [float(x) for x in up_8o_lst]
        if lines[i].strip() == up_9o_lst_str:
            up_9o_lst = re.findall('-?\d+\.?\d*', lines[i + 1])
            up_9o_lst = [float(x) for x in up_9o_lst]
            
        if lines[i].strip() == down_8n_lst_str:
            down_8n_lst = re.findall('-?\d+\.?\d*', lines[i + 1])
            down_8n_lst = [float(x) for x in down_8n_lst]
        if lines[i].strip() == down_9n_lst_str:
            down_9n_lst = re.findall('-?\d+\.?\d*', lines[i + 1])
            down_9n_lst = [float(x) for x in down_9n_lst]
        if lines[i].strip() == down_8o_lst_str:
            down_8o_lst = re.findall('-?\d+\.?\d*', lines[i + 1])
            down_8o_lst = [float(x) for x in down_8o_lst]
        if lines[i].strip() == down_9o_lst_str:
            down_9o_lst = re.findall('-?\d+\.?\d*', lines[i + 1])
            down_9o_lst = [float(x) for x in down_9o_lst]
    f.close()
    
    for locc in range(5000):
        # 随机lai
        lai = lai_lst[0] + (lai_lst[-1] - lai_lst[0]) * random.random()
        # 随机水汽
        wv = wv_lst[0] + (wv_lst[-1] - wv_lst[0]) * random.random()
        # 随机tv ts
        tv = tv_lst[0] + (tv_lst[-1] - tv_lst[0]) * random.random()
        ts = tv + 20 * random.random()
        # 计算8n
        tran_8n = tran_8n_lst[0] * wv ** 2 + tran_8n_lst[1] * wv + tran_8n_lst[2]
        up_8n = up_8n_lst[0] * wv ** 2 + up_8n_lst[1] * wv + up_8n_lst[2]
        
        # 其他
        tran_9n = tran_9n_lst[0] * tran_8n ** 2 + tran_9n_lst[1] * tran_8n + tran_9n_lst[2]
        tran_8o = tran_8o_lst[0] * tran_8n ** 2 + tran_8o_lst[1] * tran_8n + tran_8o_lst[2]
        tran_9o = tran_9o_lst[0] * tran_8n ** 2 + tran_9o_lst[1] * tran_8n + tran_9o_lst[2]
        
        up_8o = up_8o_lst[0] * up_8n ** 2 + up_8o_lst[1] * up_8n + up_8o_lst[2]
        up_9n = up_9n_lst[0] * up_8n ** 2 + up_9n_lst[1] * up_8n + up_9n_lst[2]
        up_9o = up_9o_lst[0] * up_8n ** 2 + up_9o_lst[1] * up_8n + up_9o_lst[2]
        
        down_8n = down_8n_lst[0] * up_8n ** 2 + down_8n_lst[1] * up_8n + down_8n_lst[2]
        down_8o = down_8o_lst[0] * up_8n ** 2 + down_8o_lst[1] * up_8n + down_8o_lst[2]
        down_9n = down_9n_lst[0] * up_8n ** 2 + down_9n_lst[1] * up_8n + down_9n_lst[2]
        down_9o = down_9o_lst[0] * up_8n ** 2 + down_9o_lst[1] * up_8n + down_9o_lst[2]
        # cep
        cep_8_n = toaCepCbtp.calCep(lai, nadir, ev_8, es_8)
        cep_9_n = toaCepCbtp.calCep(lai, nadir, ev_9, es_9)
        
        cep_8_o = toaCepCbtp.calCep(lai, oblique, ev_8, es_8)
        cep_9_o = toaCepCbtp.calCep(lai, oblique, ev_9, es_9)
        
        # cbtp
        cbtpev_8_n, cbtpes_8_n = toaCepCbtp.calCbtpModel(lai, es_8, ev_8, nadir)
        cbtpev_9_n, cbtpes_9_n = toaCepCbtp.calCbtpModel(lai, es_9, ev_9, nadir)
        
        cbtpev_8_o, cbtpes_8_o = toaCepCbtp.calCbtpModel(lai, es_8, ev_8, oblique)
        cbtpev_9_o, cbtpes_9_o = toaCepCbtp.calCbtpModel(lai, es_9, ev_9, oblique)
        
        # T -> L
        ls8 = c1 / (np.pi * pow(band1, 5) * (np.exp(c2 / band1 / ts) - 1))
        ls9 = c1 / (np.pi * pow(band2, 5) * (np.exp(c2 / band2 / ts) - 1))
        
        lv8 = c1 / (np.pi * pow(band1, 5) * (np.exp(c2 / band1 / tv) - 1))
        lv9 = c1 / (np.pi * pow(band2, 5) * (np.exp(c2 / band2 / tv) - 1))
        
        # toa
        # 水汽推导 8n 透过率 上行辐射
        l_toa_8n = tran_8n * ((1 - cep_8_n) * down_8n + cbtpev_8_n * lv8 + cbtpes_8_n * ls8) + up_8n
        l_toa_8o = tran_8o * ((1 - cep_8_o) * down_8o + cbtpev_8_o * lv8 + cbtpes_8_o * ls8) + up_8o
        l_toa_9n = tran_9n * ((1 - cep_9_n) * down_9n + cbtpev_9_n * lv9 + cbtpes_9_n * ls9) + up_9n 
        l_toa_9o = tran_9o * ((1 - cep_9_o) * down_9o + cbtpev_9_o * lv9 + cbtpes_9_o * ls9) + up_9o
        
        # l -> t
        toa_8n = c2 / band1 / np.log((c1 / (np.pi * pow(band1, 5) * l_toa_8n))+1)
        toa_8o = c2 / band1 / np.log((c1 / (np.pi * pow(band1, 5) * l_toa_8o))+1)
        toa_9n = c2 / band2 / np.log((c1 / (np.pi * pow(band2, 5) * l_toa_9n))+1)
        toa_9o = c2 / band2 / np.log((c1 / (np.pi * pow(band2, 5) * l_toa_9o))+1)
        
        x_data[locc][0] = toa_8n
        x_data[locc][1] = toa_8o
        x_data[locc][2] = toa_9n
        x_data[locc][3] = toa_9o
        x_data[locc][4] = lai
        
        y_data[locc][0] = tv
        y_data[locc][1] = ts
        y_data[locc][2] = wv
    
    return x_data, y_data

def ml(x_data, y_data):
    # 标准化函数
    x_minMax = preprocessing.MinMaxScaler()
    y_minMax = preprocessing.MinMaxScaler()
    # 标准化
    # x = x_minMax.fit_transform(x_data)
    # y = y_minMax.fit_transform(y_data)
    x = x_data
    y = y_data

    # 28 
    # x_train, x_test, y_train, y_test = train_test_split(x,y,test_size = 0.2)
    
    x_train = x 
    y_train = y 
    
    # 2016-06-06
    # x_test = np.array([286.97, 285.71, 286.88, 285.52, 1.5]).reshape(1, 5)
    # y_test = np.array([288.5, 289.7, 0.5]).reshape(1, 3)
    
    # 2016-05-06
    x_test = np.array([298.31, 296.83, 297.20, 294.75, 1.5]).reshape(1, 5)
    y_test = np.array([297.0, 302.4, 0.5]).reshape(1, 3)
    
    # 标准化
    scaler = preprocessing.StandardScaler().fit(x_train)
    x_train = scaler.transform(x_train)
    x_test = scaler.transform(x_test)
    # 模型构建
    #fit1 = MLPRegressor(
     #   hidden_layer_sizes=(100,50), activation='relu',solver='adam', alpha=0.01,max_iter=200)
    fit1 = MLPRegressor(hidden_layer_sizes=(120, 100, 150), activation='relu', solver='adam', alpha=0.0001,
                        batch_size='auto',
                         learning_rate_init=0.001, power_t=0.5, max_iter=5000, shuffle=True,
                        random_state=1, tol=0.0001, verbose=False, warm_start=False, momentum=0.9,
                        nesterovs_momentum=True,
                        early_stopping=False, beta_1=0.9, beta_2=0.999, epsilon=1e-08)
    print("fitting model right now")
    print(x_train.shape,y_train.shape)
    fit1.fit(x_train,y_train)
    pred1_train = fit1.predict(x_train)
    
    '''计算训练集 MSE'''
    mse_1 = mean_squared_error(pred1_train,y_train)
    print ("Train ERROR = ", mse_1)
    #print(mean_squared_error(pred1_train[:,1],y_train[:,1]))
    #print(mean_squared_error(pred1_train[:,2],y_train[:,2]))
    '''计算测试集mse'''
    pred1_test = fit1.predict(x_test)
    mse_2 = mean_squared_error(pred1_test,y_test)
    print ("Test ERROR = ", mse_2)
    
    print(y_test[0])
    print(pred1_test[0])

lai_lst = [0.5, 3.0]
# wv_lst = [0.45, 0.55]
wv_lst = [0.1, 4.0]
tv_lst = [285, 310]
# 11.3um-2018.12.17
# 获取9波段植被 土壤发射率
es_9=0.955025
ev_9=0.967125
# 10.6um-2018.12.17
# 获取8波段植被 土壤发射率
es_8=0.95273
ev_8=0.967765
nadir = 0
oblique = 55
wvFile = './data/WaterVaporResult.txt'
niheFile = './data/MODTRANResult.txt'
x_data, y_data = genData(lai_lst, wv_lst, tv_lst, es_8, es_9, ev_8, ev_9, nadir, oblique, wvFile, niheFile)
ml(x_data, y_data)