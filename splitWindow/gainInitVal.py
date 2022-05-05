import gdal
import sys
import linecache
import numpy as np
""" 从一张影像上获取初值

@gainInitVal: 从原始影像获取植被和土壤组分温度初值
initTif: 原始影像
LAITif: LAI反演结果影像
splitWindowFile: 劈窗结果文件
ang: 两个角度
"""
def gainInitVal(initTif_8_nadir, initTif_8_oblique, initTif_9_nadir, initTif_9_oblique, LAITif, splitWindowFile, ang_1, ang_2):
    # 引入必要的参数
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
    
    # 两个角度的弧度值
    sita = [ang_1 / 180 * np.pi, ang_2 / 180 * np.pi]
    
    # 读取预处理后双角度双通道影像
    ds1 = gdal.Open(initTif_8_nadir)
    ds2 = gdal.Open(initTif_8_oblique)
    ds3 = gdal.Open(initTif_9_nadir)
    ds4 = gdal.Open(initTif_9_oblique)
    
    row = ds1.RasterXSize
    col = ds1.RasterYSize
    
    bt_8_nadir = ds1.ReadAsArray(0, 0, row, col)
    bt_8_oblique = ds3.ReadAsArray(0, 0, row, col)
    bt_9_nadir = ds2.ReadAsArray(0, 0, row, col)
    bt_9_oblique = ds4.ReadAsArray(0, 0, row, col)
    
    # 读取劈窗结果
    line_8_nadir = linecache.getline(splitWindowFile, 2).strip()
    line_8_oblique = linecache.getline(splitWindowFile, 3).strip()
    line_9_nadir = linecache.getline(splitWindowFile, 4).strip()
    line_9_oblique = linecache.getline(splitWindowFile, 5).strip()
    
    # 将数据裁出
    nadir_8_para = list(map(float, line_8_nadir[line_8_nadir.find("[") + 1: -1].strip().split()))
    oblique_8_para = list(map(float, line_8_oblique[line_8_oblique.find("[") + 1: -1].strip().split()))
    nadir_9_para = list(map(float, line_9_nadir[line_9_nadir.find("[") + 1: -1].strip().split()))
    oblique_9_para = list(map(float, line_9_oblique[line_9_oblique.find("[") + 1: -1].strip().split()))
    
    # 劈窗处理 a * x + b * (x - y) + c * (x - y)**2 + d)
    ground_nadir_8 = nadir_8_para[0] * bt_8_nadir + nadir_8_para[1] * (bt_8_nadir - bt_9_nadir) + nadir_8_para[2] * (bt_8_nadir - bt_9_nadir)**2 + nadir_8_para[3]
    ground_oblique_8 = oblique_8_para[0] * bt_8_oblique + oblique_8_para[1] * (bt_8_oblique - bt_9_oblique) + oblique_8_para[2] * (bt_8_oblique - bt_9_oblique)**2 + oblique_8_para[3]
    gournd_nadir_9 = nadir_9_para[0] * bt_8_nadir + nadir_9_para[1] * (bt_8_nadir - bt_9_nadir) + nadir_9_para[2] * (bt_8_nadir - bt_9_nadir)**2 + nadir_9_para[3]
    ground_oblique_9 = oblique_9_para[0] * bt_8_oblique + oblique_9_para[1] * (bt_8_oblique - bt_9_oblique) + oblique_9_para[2] * (bt_8_oblique - bt_9_oblique)**2 + oblique_9_para[3]
    
    # 得到地面亮温结果可以先保存
    np.savetxt('./result/ground_nadir_8.txt', ground_nadir_8)
    np.savetxt('./result/ground_oblique_8.txt', ground_oblique_8)
    np.savetxt('./result/gournd_nadir_9.txt', gournd_nadir_9)
    np.savetxt('./result/ground_oblique_9.txt', ground_oblique_9)
    
    # 读取LAI
    dsLAI = gdal.Open(LAITif)
    rowLAI = dsLAI.RasterXSize
    colLAI = dsLAI.RasterYSize
    # 如果LAI的行列和原始影像不一致 那么退出程序
    if rowLAI != row or colLAI != col:
        print('Something wrong happens')
        sys.exit()
        
    LAIMatrix = dsLAI.ReadAsArray(0, 0, row, col)
    
    # 准备结果容器
    L_v_8 = np.zeros((row, col), dtype=float)
    L_v_9 = np.zeros((row, col), dtype=float)
    L_s_8 = np.zeros((row, col), dtype=float)
    L_s_9 = np.zeros((row, col), dtype=float)
    
    for xx in range(row):
        for yy in range(col):
            L_8_ang_1 = ground_nadir_8[xx][yy]
            L_8_ang_2 = ground_oblique_8[xx][yy]
            L_9_ang_1 = gournd_nadir_9[xx][yy]
            L_9_ang_2 = ground_oblique_9[xx][yy]
            
            LAI = LAIMatrix[xx][yy]
            
            if (np.isnan(LAI) == False):

                # i0 拦截率
                i0_ang_1 = 1 - np.exp(-0.5 * LAI / np.cos(sita[0]))
                i0_ang_2 = 1 - np.exp(-0.5 * LAI / np.cos(sita[1]))

                # i0hemi 半球平均拦截率，经验公式
                i0hemi = 1 - np.exp(-0.825 * LAI)

                # start数值积分法计算TIR-P
                fun_ang_1_up = 0
                fun_ang_2_up = 0
                fun_ang_1_down = 0
                fun_ang_2_down = 0

                for j in np.arange(0, float(LAI), 0.01):
                    fun_ang_1_up = fun_ang_1_up + np.exp(-0.5 * j / np.cos(sita[0])) * 0.5 / np.cos(sita[0]) * 0.5 * np.exp(-0.8 * pow(j, 0.9)) / (1 - np.exp(-0.5 * LAI / np.cos(sita[0]))) * 0.01
                    fun_ang_2_up = fun_ang_2_up + np.exp(-0.5 * j / np.cos(sita[1])) * 0.5 / np.cos(sita[1]) * 0.5 * np.exp(-0.8 * pow(j, 0.9)) / (1 - np.exp(-0.5 * LAI / np.cos(sita[1]))) * 0.01
                    fun_ang_1_down = fun_ang_1_down + np.exp(-0.5 * j / np.cos(sita[0])) * 0.5 / np.cos(sita[0]) * 0.5 * np.exp(-0.8 * pow(LAI - j, 0.9)) / (1 - np.exp(-0.5 * LAI / np.cos(sita[0]))) * 0.01
                    fun_ang_2_down = fun_ang_2_down + np.exp(-0.5 * j / np.cos(sita[1])) * 0.5 / np.cos(sita[1]) * 0.5 * np.exp(-0.8 * pow(LAI - j, 0.9)) / (1 - np.exp(-0.5 * LAI / np.cos(sita[1]))) * 0.01
                # elup 向上逃逸概率 eldown向下逃逸概率
                elup_ang_1 = round(fun_ang_1_up, 4)
                elup_ang_2 = round(fun_ang_2_up, 4)
                eldown_ang_1 = round(fun_ang_1_down, 4)
                eldown_ang_2 = round(fun_ang_2_down, 4) 

                # p 再碰撞概率
                p_ang_1 = 1 - elup_ang_1 - eldown_ang_1
                p_ang_2 = 1 - elup_ang_2 - eldown_ang_2

                # 计算DCE123简化模型
                # c1 光子被叶片反射后碰到叶片的概率
                c1_8_ang_1 = (1 - ev_8) * p_ang_1                           
                c1_8_ang_2 = (1 - ev_8) * p_ang_2
                c1_9_ang_1 = (1 - ev_9) * p_ang_1
                c1_9_ang_2 = (1 - ev_9) * p_ang_2
                # c2 光子被土壤反射后碰撞叶片的概率
                c2_8_ang_1 = (1 - es_8) * i0hemi
                c2_8_ang_2 = (1 - es_8) * i0hemi
                c2_9_ang_1 = (1 - es_9) * i0hemi
                c2_9_ang_2 = (1 - es_9) * i0hemi
                # c3 光子被叶片反射后碰撞土壤的概率
                c3_8_ang_1 = (1 - ev_8) * eldown_ang_1
                c3_8_ang_2 = (1 - ev_8) * eldown_ang_2
                c3_9_ang_1 = (1 - ev_9) * eldown_ang_1
                c3_9_ang_2 = (1 - ev_9) * eldown_ang_2

                # CBT-P模型 植被和土壤有效发射率
                cbtpev_8_ang_1 = i0_ang_1 * ev_8 * (1 + c1_8_ang_1 + c1_8_ang_1 * c1_8_ang_1 + c3_8_ang_1 * c2_8_ang_1) + (1 - i0_ang_1) * ev_8 * (c2_8_ang_1 + c2_8_ang_1 * c1_8_ang_1)
                cbtpev_8_ang_2 = i0_ang_2 * ev_8 * (1 + c1_8_ang_2 + c1_8_ang_2 * c1_8_ang_2 + c3_8_ang_2 * c2_8_ang_2) + (1 - i0_ang_2) * ev_8 * (c2_8_ang_2 + c2_8_ang_2 * c1_8_ang_2)
                cbtpev_9_ang_1 = i0_ang_1 * ev_9 * (1 + c1_9_ang_1 + c1_9_ang_1 * c1_9_ang_1 + c3_9_ang_1 * c2_9_ang_1) + (1 - i0_ang_1) * ev_9 * (c2_9_ang_1 + c2_9_ang_1 * c1_9_ang_1)
                cbtpev_9_ang_2 = i0_ang_2 * ev_9 * (1 + c1_9_ang_2 + c1_9_ang_2 * c1_9_ang_2 + c3_9_ang_2 * c2_9_ang_2) + (1 - i0_ang_2) * ev_9 * (c2_9_ang_2 + c2_9_ang_2 * c1_9_ang_2)
            
                cbtpes_8_ang_1 = i0_ang_1 * es_8 * (c3_8_ang_1 + c1_8_ang_1 * c3_8_ang_1) + (1 - i0_ang_1) * es_8 * (1 + c2_8_ang_1 * c3_8_ang_1)
                cbtpes_8_ang_2 = i0_ang_2 * es_8 * (c3_8_ang_2 + c1_8_ang_2 * c3_8_ang_2) + (1 - i0_ang_2) * es_8 * (1 + c2_8_ang_2 * c3_8_ang_2)
                cbtpes_9_ang_1 = i0_ang_1 * es_9 * (c3_9_ang_1 + c1_9_ang_1 * c3_9_ang_1) + (1 - i0_ang_1) * es_9 * (1 + c2_9_ang_1 * c3_9_ang_1)
                cbtpes_9_ang_2 = i0_ang_2 * es_9 * (c3_9_ang_2 + c1_9_ang_2 * c3_9_ang_2) + (1 - i0_ang_2) * es_9 * (1 + c2_9_ang_2 * c3_9_ang_2)

                L_v_8[xx][yy] = (cbtpes_8_ang_1 * L_8_ang_2 - cbtpes_8_ang_2 * L_8_ang_1) / (cbtpev_8_ang_2 * cbtpes_8_ang_1 - cbtpes_8_ang_2 * cbtpev_8_ang_1)
                L_s_8[xx][yy] = (L_8_ang_1 - cbtpev_8_ang_1 * L_v_8[xx][yy]) / cbtpes_8_ang_1

                L_v_9[xx][yy] = (cbtpes_9_ang_1 * L_9_ang_2 - cbtpes_9_ang_2 * L_9_ang_1) / (cbtpev_9_ang_2 * cbtpes_9_ang_1 - cbtpes_9_ang_2 * cbtpev_9_ang_1)
                L_s_9[xx][yy] = (L_9_ang_1 - cbtpev_9_ang_1 * L_v_9[xx][yy]) / cbtpes_9_ang_1
            else:
                L_v_8[xx][yy] = np.nan
                L_s_8[xx][yy] = np.nan

                L_v_9[xx][yy] = np.nan
                L_s_9[xx][yy] = np.nan
                
    Tv_8 = c2 / band1 / np.log((c1 / (np.pi * pow(band1, 5) * L_v_8))+1)
    Ts_8 = c2 / band1 / np.log((c1 / (np.pi * pow(band1, 5) * L_s_8))+1)

    Tv_9 = c2 / band2 / np.log((c1 / (np.pi * pow(band2, 5) * L_v_9))+1)
    Ts_9 = c2 / band2 / np.log((c1 / (np.pi * pow(band2, 5) * L_s_9))+1)

    Tv = (Tv_8 + Tv_9) / 2
    Ts = (Ts_8 + Ts_9) / 2
    
    np.savetxt('./result/Tv_res.txt', Tv)
    np.savetxt('./result/Ts_res.txt', Ts)