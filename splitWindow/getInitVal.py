import gdal
import linecache
import numpy as np
from numpy.lib.npyio import save
""" 获取初值
@getInitTsAndTv: 获取Ts和Tv初值
s8_bt_nadir: S8 nadir TOA亮温
s8_bt_oblique: S8 oblique TOA亮温
s9_bt_nadir: S9 nadir TOA亮温
s9_bt_oblique: S9 oblique TOA亮温
LAITif: LAI图像
sat_vza_nadir: VZA nadir 角度
sat_vza_oblique: VZA oblique 角度
splitWindowFile: 劈窗结果文件

@tifPara: 获取tif参数
tifFile: tif文件

@saveTif: 保存tif
nArray: ndarray格式的二维数组
nInfo: RasterBand
targetLoc: 目标保存位置
"""
def getInitTsAndTv(s8_bt_nadir, s8_bt_oblique, s9_bt_nadir, s9_bt_oblique, LAITif, sat_vza_nadir, sat_vza_oblique, splitWindowFile):
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
    
    # 读取 s8 s9 nadir oblique 亮温数据
    s8_nadir, row, col, rasDataType = tifPara(s8_bt_nadir)
    s8_oblique, _, _, _ = tifPara(s8_bt_oblique)
    s9_nadir, _, _, _ = tifPara(s9_bt_nadir)
    s9_oblique, _, _, _ = tifPara(s9_bt_oblique)
    
    # 读取LAI数据
    LAIMatrix, _, _, _ = tifPara(LAITif)
    
    # 读取vza 数据
    vza_nadir, _, _, _ = tifPara(sat_vza_nadir)
    vza_oblique, _, _, _ = tifPara(sat_vza_oblique)
    
    # 读取劈窗系数
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
    ground_nadir_8 = nadir_8_para[0] * s8_nadir + nadir_8_para[1] * (s8_nadir - s9_nadir) + nadir_8_para[2] * (s8_nadir - s9_nadir)**2 + nadir_8_para[3]
    ground_oblique_8 = oblique_8_para[0] * s8_oblique + oblique_8_para[1] * (s8_oblique - s9_oblique) + oblique_8_para[2] * (s8_oblique - s9_oblique)**2 + oblique_8_para[3]
    gournd_nadir_9 = nadir_9_para[0] * s8_nadir + nadir_9_para[1] * (s8_nadir - s9_nadir) + nadir_9_para[2] * (s8_nadir - s9_nadir)**2 + nadir_9_para[3]
    ground_oblique_9 = oblique_9_para[0] * s8_oblique + oblique_9_para[1] * (s8_oblique - s9_oblique) + oblique_9_para[2] * (s8_oblique - s9_oblique)**2 + oblique_9_para[3]
    
    # 角度弧度化
    vza_nadir = vza_nadir / 180 * np.pi
    vza_oblique = vza_oblique / 180 * np.pi
    
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
                i0_ang_1 = 1 - np.exp(-0.5 * LAI / np.cos(vza_nadir[xx][yy]))
                i0_ang_2 = 1 - np.exp(-0.5 * LAI / np.cos(vza_oblique[xx][yy]))

                # i0hemi 半球平均拦截率，经验公式
                i0hemi = 1 - np.exp(-0.825 * LAI)

                # start数值积分法计算TIR-P
                fun_ang_1_up = 0
                fun_ang_2_up = 0
                fun_ang_1_down = 0
                fun_ang_2_down = 0

                for j in np.arange(0, float(LAI), 0.01):
                    fun_ang_1_up = fun_ang_1_up + np.exp(-0.5 * j / np.cos(vza_nadir[xx][yy])) * 0.5 / np.cos(vza_nadir[xx][yy]) * 0.5 * np.exp(-0.8 * pow(j, 0.9)) / (1 - np.exp(-0.5 * LAI / np.cos(vza_nadir[xx][yy]))) * 0.01
                    fun_ang_2_up = fun_ang_2_up + np.exp(-0.5 * j / np.cos(vza_oblique[xx][yy])) * 0.5 / np.cos(vza_oblique[xx][yy]) * 0.5 * np.exp(-0.8 * pow(j, 0.9)) / (1 - np.exp(-0.5 * LAI / np.cos(vza_oblique[xx][yy]))) * 0.01
                    fun_ang_1_down = fun_ang_1_down + np.exp(-0.5 * j / np.cos(vza_nadir[xx][yy])) * 0.5 / np.cos(vza_nadir[xx][yy]) * 0.5 * np.exp(-0.8 * pow(LAI - j, 0.9)) / (1 - np.exp(-0.5 * LAI / np.cos(vza_nadir[xx][yy]))) * 0.01
                    fun_ang_2_down = fun_ang_2_down + np.exp(-0.5 * j / np.cos(vza_oblique[xx][yy])) * 0.5 / np.cos(vza_oblique[xx][yy]) * 0.5 * np.exp(-0.8 * pow(LAI - j, 0.9)) / (1 - np.exp(-0.5 * LAI / np.cos(vza_oblique[xx][yy]))) * 0.01
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
    
    # 保存Tv和Ts初值
    saveTif(Tv, rasDataType, './result/Tv.tif')
    saveTif(Ts, rasDataType, './result/Ts.tif')
    
def tifPara(tifFile):
    ds = gdal.Open(tifFile)
    col = ds.RasterXSize
    row = ds.RasterYSize
    rasBand = ds.GetRasterBand(1)
    img = ds.ReadAsArray(0, 0, col, row)
    return img, row, col, rasBand.DataType

def saveTif(nArray, nInfo, targetLoc):
    driver = gdal.GetDriverByName('GTiff')
    op = driver.Create(targetLoc, nArray.shape[1], nArray.shape[0], 1, nInfo.DataType)
    op.GetRasterBand(1).WriteArray(nArray)