import numpy as np
import gdal
import os

# 如果不存在文件夹 则创建
dirs = './res/'
if not os.path.exists(dirs):
    os.makedirs(dirs)

def tifPara(tifFile):
    ds = gdal.Open(tifFile)
    col = ds.RasterXSize
    row = ds.RasterYSize
    rasBand = ds.GetRasterBand(1)
    img = ds.ReadAsArray(0, 0, col, row)
    geoTransform = ds.GetGeoTransform()
    proj = ds.GetProjection()
    return img, row, col, rasBand.DataType, geoTransform, proj

def saveTif(nArray, nInfo, targetLoc, geoTransform, proj):
    driver = gdal.GetDriverByName('GTiff')
    op = driver.Create(targetLoc, nArray.shape[1], nArray.shape[0], 1, nInfo)
    op.SetGeoTransform(geoTransform)
    op.SetProjection(proj)
    op.GetRasterBand(1).WriteArray(nArray)

# 加载初值
# Ts_fir, row, col, rasDataType, geoTransform, proj = tifPara('./fir/tsFir.tif')
# Tv_fir, row, col, rasDataType, geoTransform, proj = tifPara('./fir/tvFir.tif')

# wv_fir, row, col, rasDataType, geoTransform, proj = tifPara('./wv/wv_cal.tif')

wv_fir = 0.5

lup_fir = 0.127 * wv_fir * wv_fir + 0.439 * wv_fir +0.001
t_fir = -0.009 * wv_fir * wv_fir +-0.070 * wv_fir + 0.977

# 加载左值
Ts8_0, row, col, rasDataType, geoTransform, proj = tifPara('./toatif/s8ntoa.tif')
Ts9_0, row, col, rasDataType, geoTransform, proj = tifPara('./toatif/s9ntoa.tif')
Ts8_55, row, col, rasDataType, geoTransform, proj = tifPara('./toatif/s8otoa.tif')
Ts9_55, row, col, rasDataType, geoTransform, proj = tifPara('./toatif/s9otoa.tif')
# lc_lai, row, col, rasDataType, geoTransform, proj = tifPara('./lai/laitengchong11.tif')

# print(Tv_fir.shape)
# print(Ts_fir.shape)
# print(wv_fir.shape)
# print(lup_fir.shape)
# print(t_fir.shape)
# print(Ts8_0.shape)
# print(Ts9_0.shape)
# print(Ts8_55.shape)
# print(Ts9_55.shape)
# print(lc_lai.shape)

count = 0

num = 4

c1 = 3.7404e8
c2 = 14387

band1 = 10.852
band2 = 12

# 11.3um-2018.12.17
# 获取9波段植被 土壤发射率
es_9=0.955025
ev_9=0.967125
# 10.6um-2018.12.17
# 获取8波段植被 土壤发射率
es_8=0.95273
ev_8=0.967765

# 算左值亮温
L_bt_s8_0 = c1 / (np.pi * pow(band1, 5) * (np.exp(c2 / band1 / Ts8_0) - 1))
L_bt_s8_55 = c1 / (np.pi * pow(band1, 5) * (np.exp(c2 / band1 / Ts8_55) - 1))
L_bt_s9_0 = c1 / (np.pi * pow(band2, 5) * (np.exp(c2 / band2 / Ts9_0) - 1))
L_bt_s9_55 = c1 / (np.pi * pow(band2, 5) * (np.exp(c2 / band2 / Ts9_55) - 1))

# 两个角度的弧度值 0 55
sita = [10 / 180 * np.pi, 55 / 180 * np.pi]

tv_fin = np.zeros((row, col), dtype=float)
ts_fin = np.zeros((row, col), dtype=float)
t_fin = np.zeros((row, col), dtype=float)
lup_fin = np.zeros((row, col), dtype=float)

# # 方程组容器
# f = np.zeros((4),dtype=float)
# # 初值
# x_fir = np.ones((4),dtype=float)

for xx in range(row):
    for yy in range(col):

        # 方程组容器
        f = np.zeros((4),dtype=float)
        # 初值
        x_fir = np.ones((4),dtype=float)

        x_fir[0] = 290#Tv_fir[xx][yy] # Tv
        x_fir[1] = 300#Ts_fir[xx][yy] # Ts
        x_fir[2] = t_fir
        x_fir[3] = lup_fir

        # LAI = lc_lai[xx][yy]
        LAI = 2.222

        if (np.isnan(LAI) == False and LAI > 0):

            # i0 拦截率
            i0_0 = 1 - np.exp(-0.5 * LAI / np.cos(sita[0]))
            i0_55 = 1 - np.exp(-0.5 * LAI / np.cos(sita[1]))

            # i0hemi 半球平均拦截率，经验公式
            i0hemi = 1 - np.exp(-0.825 * LAI)

            # start数值积分法计算TIR-P
            fun_0_up = 0
            fun_55_up = 0
            fun_0_down = 0
            fun_55_down = 0
            for j in np.arange(0, float(LAI), 0.01):
                fun_0_up = fun_0_up + np.exp(-0.5 * j / np.cos(sita[0])) * 0.5 / np.cos(sita[0]) * 0.5 * np.exp(-0.8 * pow(j, 0.9)) / (1 - np.exp(-0.5 * LAI / np.cos(sita[0]))) * 0.01
                fun_55_up = fun_55_up + np.exp(-0.5 * j / np.cos(sita[1])) * 0.5 / np.cos(sita[1]) * 0.5 * np.exp(-0.8 * pow(j, 0.9)) / (1 - np.exp(-0.5 * LAI / np.cos(sita[1]))) * 0.01
                fun_0_down = fun_0_down + np.exp(-0.5 * j / np.cos(sita[0])) * 0.5 / np.cos(sita[0]) * 0.5 * np.exp(-0.8 * pow(LAI - j, 0.9)) / (1 - np.exp(-0.5 * LAI / np.cos(sita[0]))) * 0.01
                fun_55_down = fun_55_down + np.exp(-0.5 * j / np.cos(sita[1])) * 0.5 / np.cos(sita[1]) * 0.5 * np.exp(-0.8 * pow(LAI - j, 0.9)) / (1 - np.exp(-0.5 * LAI / np.cos(sita[1]))) * 0.01
            # elup 向上逃逸概率 eldown向下逃逸概率
            elup_0 = round(fun_0_up, 4)
            elup_55 = round(fun_55_up, 4)
            eldown_0 = round(fun_0_down, 4)
            eldown_55 = round(fun_55_down, 4) 
            # p 再碰撞概率
            p_0 = 1 - elup_0 - eldown_0
            p_55 = 1 - elup_55 - eldown_55

            # cavity 反射后未被植被吸收的部分，可能前向漫反射，也可能后向漫反射
            cavity_8_0 = 1 - p_0 * (1 - ev_8)
            cavity_8_55 = 1 - p_55 * (1 - ev_8)
            cavity_9_0 = 1 - p_0 * (1 - ev_9)
            cavity_9_55 = 1 - p_55 * (1 - ev_9)

            # rc1 冠层前向漫反射率
            # rc2 冠层后向漫反射率
            rc1_8_0 = (1 - ev_8) * eldown_0 / cavity_8_0
            rc1_8_55 = (1 - ev_8) * eldown_55 / cavity_8_55
            rc1_9_0 = (1 - ev_9) * eldown_0 / cavity_9_0
            rc1_9_55 = (1 - ev_9) * eldown_55 / cavity_9_55

            rc2_8_0 = (1 - ev_8) * elup_0 / cavity_8_0
            rc2_8_55 = (1 - ev_8) * elup_55 / cavity_8_55
            rc2_9_0 = (1 - ev_9) * elup_0 / cavity_9_0
            rc2_9_55 = (1 - ev_9) * elup_55 / cavity_9_55

            # 计算e1 e2 e3 e4 e5
            mulsv_8_0 = 1 - rc2_8_0 * (1 - es_8) * i0hemi
            mulsv_8_55 = 1 - rc2_8_55 * (1 - es_8) * i0hemi
            mulsv_9_0 = 1 - rc2_9_0 * (1 - es_9) * i0hemi
            mulsv_9_55 = 1 - rc2_9_55 * (1 - es_9) * i0hemi        

            e1_8_0 = i0_0 * ev_8 / cavity_8_0
            e1_8_55 = i0_55 * ev_8 / cavity_8_55
            e1_9_0 = i0_0 * ev_9 / cavity_9_0
            e1_9_55 = i0_55 * ev_9 / cavity_9_55

            e2_8_0 = (1 - i0_0) * (1 - es_8) * i0hemi * ev_8 / cavity_8_0 / mulsv_8_0
            e2_8_55 = (1 - i0_55) * (1 - es_8) * i0hemi * ev_8 / cavity_8_55 / mulsv_8_55
            e2_9_0 = (1 - i0_0) * (1 - es_9) * i0hemi * ev_9 / cavity_9_0 / mulsv_9_0
            e2_9_55 = (1 - i0_55) * (1 - es_9) * i0hemi * ev_9 / cavity_9_55 / mulsv_9_55

            e3_8_0 = i0_0 * rc1_8_0 * (1 - es_8) * i0hemi * ev_8 / cavity_8_0 / mulsv_8_0
            e3_8_55 = i0_55 * rc1_8_55 * (1 - es_8) * i0hemi * ev_8 / cavity_8_55 / mulsv_8_55
            e3_9_0 = i0_0 * rc1_9_0 * (1 - es_9) * i0hemi * ev_9 / cavity_9_0 / mulsv_9_0
            e3_9_55 = i0_55 * rc1_9_55 * (1 - es_9) * i0hemi * ev_9 / cavity_9_55 / mulsv_9_55

            e4_8_0 = (1 - i0_0) * es_8 / mulsv_8_0
            e4_8_55 = (1 - i0_55) * es_8 / mulsv_8_55 
            e4_9_0 = (1 - i0_0) * es_9 / mulsv_9_0
            e4_9_55 = (1 - i0_55) * es_9 / mulsv_9_55 
        
            e5_8_0 = i0_0 * rc1_8_0 * es_8 / mulsv_8_0
            e5_8_55 = i0_55 * rc1_8_55 * es_8 / mulsv_8_55
            e5_9_0 = i0_0 * rc1_9_0 * es_9 / mulsv_9_0
            e5_9_55 = i0_55 * rc1_9_55 * es_9 / mulsv_9_55

            # ecaoP4sail 冠层方向发射率
            ecaoP4sail_8_0 = e1_8_0 + e2_8_0 + e3_8_0 + e4_8_0 + e5_8_0
            ecaoP4sail_8_55 = e1_8_55 + e2_8_55 + e3_8_55 + e4_8_55 + e5_8_55
            ecaoP4sail_9_0 = e1_9_0 + e2_9_0 + e3_9_0 + e4_9_0 + e5_9_0
            ecaoP4sail_9_55 = e1_9_55 + e2_9_55 + e3_9_55 + e4_9_55 + e5_9_55

            # 计算DCE123简化模型
            # c1 光子被叶片反射后碰到叶片的概率
            c1_8_0 = (1 - ev_8) * p_0                           
            c1_8_55 = (1 - ev_8) * p_55
            c1_9_0 = (1 - ev_9) * p_0
            c1_9_55 = (1 - ev_9) * p_55
            # c2 光子被土壤反射后碰撞叶片的概率
            c2_8_0 = (1 - es_8) * i0hemi
            c2_8_55 = (1 - es_8) * i0hemi
            c2_9_0 = (1 - es_9) * i0hemi
            c2_9_55 = (1 - es_9) * i0hemi
            # c3 光子被叶片反射后碰撞土壤的概率
            c3_8_0 = (1 - ev_8) * eldown_0
            c3_8_55 = (1 - ev_8) * eldown_55
            c3_9_0 = (1 - ev_9) * eldown_0
            c3_9_55 = (1 - ev_9) * eldown_55

            # CBT-P模型 植被和土壤有效发射率
            cbtpev_8_0 = i0_0 * ev_8 * (1 + c1_8_0 + c1_8_0 * c1_8_0 + c3_8_0 * c2_8_0) + (1 - i0_0) * ev_8 * (c2_8_0 + c2_8_0 * c1_8_0)
            cbtpev_8_55 = i0_55 * ev_8 * (1 + c1_8_55 + c1_8_55 * c1_8_55 + c3_8_55 * c2_8_55) + (1 - i0_55) * ev_8 * (c2_8_55 + c2_8_55 * c1_8_55)
            cbtpev_9_0 = i0_0 * ev_9 * (1 + c1_9_0 + c1_9_0 * c1_9_0 + c3_9_0 * c2_9_0) + (1 - i0_0) * ev_9 * (c2_9_0 + c2_9_0 * c1_9_0)
            cbtpev_9_55 = i0_55 * ev_9 * (1 + c1_9_55 + c1_9_55 * c1_9_55 + c3_9_55 * c2_9_55) + (1 - i0_55) * ev_9 * (c2_9_55 + c2_9_55 * c1_9_55)
        
            cbtpes_8_0 = i0_0 * es_8 * (c3_8_0 + c1_8_0 * c3_8_0) + (1 - i0_0) * es_8 * (1 + c2_8_0 * c3_8_0)
            cbtpes_8_55 = i0_55 * es_8 * (c3_8_55 + c1_8_55 * c3_8_55) + (1 - i0_55) * es_8 * (1 + c2_8_55 * c3_8_55)
            cbtpes_9_0 = i0_0 * es_9 * (c3_9_0 + c1_9_0 * c3_9_0) + (1 - i0_0) * es_9 * (1 + c2_9_0 * c3_9_0)
            cbtpes_9_55 = i0_55 * es_9 * (c3_9_55 + c1_9_55 * c3_9_55) + (1 - i0_55) * es_9 * (1 + c2_9_55 * c3_9_55)
            
            # print(cbtpev_8_0)
            # print(cbtpev_8_55)
            # print(cbtpev_9_0)
            # print(cbtpev_9_55)
            
            # print(cbtpes_8_0)
            # print(cbtpes_8_55)
            # print(cbtpes_9_0)
            # print(cbtpes_9_55)

            # 0.298x*x+0.876x+-0.176 1.334x*x+-0.787x+0.476 2.277x*x+-1.947x+0.697
            # -0.044x*x+1.467x+0.008 -0.164x*x+1.684x+-0.009 -0.288x*x+2.324x+0.006
            # 0.007x*x+1.026x+0.002 -0.030x*x+1.485x+0.014 -0.159x*x+1.770x+-0.016 -0.270x*x+2.416x+-0.002
            
            def Fun(x,num): 
                i = num
                f = np.zeros((i),dtype=float)
                f[0] = x[2] * ((1 - ecaoP4sail_8_0) * (0.008 * x[3] * x[3] + 1.026 * x[3] + 0.002) + (cbtpev_8_0 * c1 / (np.pi * pow(band1, 5) * (np.exp(c2 / band1 / x[0]) - 1))) + (cbtpes_8_0 * c1 / (np.pi * pow(band1, 5) * (np.exp(c2 / band1 / x[1]) - 1)))) + x[3] - L_bt_s8_0[xx][yy]
                f[1] = (0.472 * x[2] * x[2] + 0.714 * x[2] + -0.190) * ((1 - ecaoP4sail_8_55) * (-0.051 * x[3] * x[3] + 1.670 * x[3] + 0.017) + (cbtpev_8_55 * c1 / (np.pi * pow(band1, 5) * (np.exp(c2 / band1 / x[0]) - 1))) + (cbtpes_8_55 * c1 / (np.pi * pow(band1, 5) * (np.exp(c2 / band1 / x[1]) - 1)))) + (-0.069 * x[3] * x[3] + 1.652 * x[3] + 0.010) - L_bt_s8_55[xx][yy]
                # f[2] = (1.593 * x[2] * x[2] + -1.250 * x[2] + 0.681) * ((1 - ecaoP4sail_9_0) * (-0.187 * x[3] * x[3] + 1.810 * x[3] + -0.017) + (cbtpev_9_0 * c1 / (np.pi * pow(band2, 5) * (np.exp(c2 / band2 / x[0]) - 1))) + (cbtpes_9_0 * c1 / (np.pi * pow(band2, 5) * (np.exp(c2 / band2 / x[1]) - 1)))) + (-0.194 * x[3] * x[3] + 1.727 * x[3] + -0.013) - L_bt_s9_0[xx][yy]
                f[2] = 2.587 * x[2] * x[2] +-12.181 * x[2] + 9.406 - x[3]
                f[3] = (3.133 * x[2] * x[2] - 3.295 * x[2] + 1.190) * ((1 - ecaoP4sail_9_55) * (-0.369 * x[3] * x[3] + 2.733 * x[3] + -0.000) + (cbtpev_9_55 * c1 / (np.pi * pow(band2, 5) * (np.exp(c2 / band2 / x[0]) - 1))) + (cbtpes_9_55 * c1 / (np.pi * pow(band2, 5) * (np.exp(c2 / band2 / x[1]) - 1)))) + (-0.393 * x[3] * x[3] + 2.638 * x[3] + 0.004) - L_bt_s9_55[xx][yy]
                return f

            def dfun(x,num):                         #计算雅可比矩阵的逆矩阵
                df = np.zeros((num,num),dtype=float)
                dx = 0.00001                           # 
                x1 = np.copy(x)
                for i in range(0,num):              # 求导数，i是列，j是行
                    for j in range(0,num):
                        x1 = np.copy(x)
                        x1[j] = x1[j]+dx           #x+dx
                        df[i,j] = (Fun(x1,num)[i]-Fun(x,num)[i])/dx   #f(x+dx)-f（x）/dx
                df_1 = np.linalg.inv(df)                              #计算逆矩阵
                return df_1

            def Newton(x,num):
                x1 = np.copy(x)
                i = 0
                delta = np.copy(x)
            #    dfun0=dfun(x,num)          #也可以使用简化牛顿法
                try:
                    while(np.sum(abs(delta)) > 1.e-8 and i < 20):  #控制循环次数
                        x1 = x-np.dot(dfun(x,num),Fun(x,num))  #公式
                        delta = x1-x                     #比较x的变化
                        x = x1
                        i = i+1
                        # print(x)
                    return x
                except np.linalg.LinAlgError as err:
                    pass
                return [np.nan, np.nan, np.nan, np.nan]
        
            count += 1
            # print("第" + str(count) + "次循环")
            # print('x_fir: ' + str(x_fir))
            a = Newton(x_fir,num)
            print(a)     

            # print("--------")
            # print('----差别----')
            # print(Tv_fir[xx][yy] - a[0])
            # print(Ts_fir[xx][yy] - a[1])
            # print('---下一个---')
            if a[0] > 200 and a[0] < 350 and a[1] > 200 and a[1] < 350:
                tv_fin[xx][yy] = a[0]
                ts_fin[xx][yy] = a[1]
                t_fin[xx][yy] = a[2]
                lup_fin[xx][yy] = a[3]
            else:
                tv_fin[xx][yy] = np.nan
                ts_fin[xx][yy] = np.nan
                t_fin[xx][yy] = np.nan
                lup_fin[xx][yy] = np.nan
        else:
            count += 1
            print("第" + str(count) + "次循环")
            a = [np.nan, np.nan, np.nan, np.nan]
            # print(a)     

            # print("--------")
            # print('----差别----')
            # print(Tv_fir[xx][yy] - a[0])
            # print(Ts_fir[xx][yy] - a[1])
            # print('---下一个---')
            tv_fin[xx][yy] = a[0]
            ts_fin[xx][yy] = a[1]
            t_fin[xx][yy] = a[2]
            lup_fin[xx][yy] = a[3]

# tv_fin[np.isnan(tv_fin)] = 0
# ts_fin[np.isnan(ts_fin)] = 0

# t_fin[t_fin == 0.9] = 0
# lup_fin[lup_fin == 0.5] = 0

# np.savetxt("./res/tv_fin.txt", tv_fin)
# np.savetxt("./res/ts_fin.txt", ts_fin)
# np.savetxt("./res/t_fin.txt", t_fin)
# np.savetxt("./res/lup_fin.txt", lup_fin) img, row, col, rasDataType, geoTransform, proj

saveTif(tv_fin, rasDataType, './res/tv_fin.tif', geoTransform, proj)
saveTif(ts_fin, rasDataType, './res/ts_fin.tif', geoTransform, proj)
saveTif(t_fin, rasDataType, './res/t_fin.tif', geoTransform, proj)
saveTif(lup_fin, rasDataType, './res/lup_fin.tif', geoTransform, proj)