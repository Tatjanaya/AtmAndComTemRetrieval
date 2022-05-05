import numpy as np
import belta as be
import sympy
from scipy import integrate
import pandas as pd
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import os
import readtxt

plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

def simulss(Leafpath, Soilpath, Bands=3, v_angle_min=0, v_angle_max=0, v_interval=5, s_angle_min=0, s_angle_max=0, s_interval=5, LAI_min=1, LAI_max=1, LAI_interval=1, ss_min=0, ss_max=0, ss_interval=5, vv_min=0, vv_max=0, vv_interval=5, Ha=0, Hb=1.6, nsss=0.5, Rsss=2, flag=0):
    # 如果不存在文件夹 则创建
    dirs = './res/'
    if not os.path.exists(dirs):
        os.makedirs(dirs)
    
    if v_angle_max < v_angle_min:
        print("输入数据不正确")
        os._exit(0)

    if s_angle_max < s_angle_min:
        print("输入数据不正确")
        os._exit(0)

    if LAI_max < LAI_min:
        print("输入数据不正确")

    if v_interval == 0:
        v_interval = 5

    if s_interval == 0:
        s_interval = 5

    v_number = (int)((v_angle_max + v_interval - v_angle_min) / v_interval)
    s_number = (int)((s_angle_max + s_interval - s_angle_min) / s_interval)
    vv_number = (int)((vv_max + vv_interval - vv_min) / vv_interval)
    ss_number = (int)((ss_max + ss_interval - ss_min) / ss_interval)

    Ha = Ha + 0.001  # 防止为0时发生异常
    Hb = Hb + 0.001

    G = 0.5
    hsss = Hb

    length_min = 0.00000001
    length_max = 1

    leaf_ratio0 = 0.99999999
    nonleaf_ratio0 = 0.00000001

    total = vv_number * ss_number * s_number * v_number
    k_ratio = np.zeros((total, 10 + Bands))

    # 植被反射率和波段信息
    Rcs = [None] * Bands
    Ban = [None] * Bands
    Rct = [None] * Bands
    Rcs, Rct = readtxt.read2txt(Leafpath, Rcs, Rct, Bands, flag)
    Ban = readtxt.readtxt(Leafpath, Ban, Bands, flag)

    # 土壤反射率
    Rgs = [None] * Bands
    Rgs = readtxt.read1txt(Soilpath, Rgs, Bands, flag)  # 用户输入路径

    # 非叶反射率,默认
    Rncs = [0.001] * Bands

    # 各波段天空散射光比例
    belta = np.zeros((s_number, Bands))
    belta = be.Bel(belta, s_angle_min, s_interval, Ban)

    # 单叶片反照率
    w1 = [0.0] * Bands
    for i in range(Bands):
        w1[i] = Rcs[i] + Rct[i]

    Rm = np.zeros((s_number, Bands))

    for LAIa in np.arange(LAI_min, LAI_max + LAI_interval, LAI_interval):
        if LAIa > LAI_max:
            break
        else:
            count = 0

            for vv_i in range(vv_number):
                for ss_j in range(ss_number):
                    for i in range(v_number):
                        for j in range(s_number):
                            vv_need = vv_i * vv_interval + vv_min
                            ss_need = ss_j * ss_interval + ss_min
                            vv_ss = abs(vv_need - ss_need)
                            m_s = s_angle_min + j * s_interval
                            m_v = v_angle_min + i * v_interval
                            lamtai = 1
                            if Rsss > 0:
                                usss = LAIa * 3 / 4 / nsss / np.pi / hsss / Rsss / Rsss
                                ksss = hsss / 2 / Rsss
                                lamtai = min(1, 3 * np.sqrt((np.cos(m_v / 180 * np.pi) ** 2) + (ksss ** 2) * (
                                            np.sin(m_v / 180 * np.pi) ** 2)) / 8 / G / Rsss / ksss / usss)
                            elif Rsss < 0:
                                vv_tocal = vv_need
                                if vv_need < 5:
                                    vv_tocal = 5
                                A1 = Ha
                                A2 = Hb
                                H = nsss
                                A1 = A1 / np.sin(vv_tocal / 180 * np.pi)
                                A2 = A2 / np.sin(vv_tocal / 180 * np.pi)
                                # 求临界角
                                c = np.arctan(A2 / H) * 180 / np.pi
                                ci = 0
                                if m_v >= c:
                                    ci = 1
                                else:
                                    p1 = np.exp(-G * LAIa / np.cos(m_v / 180 * np.pi))
                                    p = p1 + (1 - p1) * (A2 - H * np.tan(m_v / 180 * np.pi)) / (A1 + A2)
                                    ci = -np.log(p) * np.cos(m_v / 180 * np.pi) / G / LAIa
                                lamtai = ci
                            print(lamtai)

                            p = 0.7 * np.exp(0.01 * lamtai * LAIa) - 0.66 * np.exp(-0.8 * lamtai * LAIa)
                            for ii in range(Bands):
                                for jj in range(s_number):
                                    # ceita = 60.21  # 太阳入射天顶角
                                    ui = np.cos(np.pi / 180 * (s_angle_min + jj * s_interval))
                                    i0 = 1 - np.exp(-lamtai * G * LAIa / ui)  # 直射光拦截概率
                                    ia = 1 - np.exp(-math.pow((lamtai * LAIa), 0.9) * 0.8)
                                    m1 = 0.5 * (1 - belta[jj][ii]) * i0 * math.pow((w1[ii]), 2) * p * (1 - p) / (
                                                1 - w1[ii] * p)
                                    m1 = round(m1, 5)

                                    # S1 = np.exp(-0.5867 * LAIa)
                                    S1 = 0.46279 * np.exp(-0.46779 * lamtai * LAIa) + 0.01127
                                    m2 = belta[jj][ii] * (1 - S1) * ia / 2 * math.pow((w1[ii]), 2) * p * (1 - p) / (
                                                1 - w1[ii] * p)
                                    m2 = round(m2, 5)

                                    Sbs = i0 * math.pow(w1[ii], 2) * m1 * (1 - m1) / (1 - w1[ii] * m1)
                                    rc = (0.5 / ia) * Sbs
                                    m3 = (1 - belta[jj][ii]) * (1 - i0) * Rgs[ii] * ia * rc / (
                                                1 - Rgs[ii] * ia * rc) * (1 + Rgs[ii] * (1 - ia))
                                    m3 = round(m3, 5)

                                    m4 = (1 - belta[jj][ii]) * i0 * rc * Rgs[ii] / (1 - Rgs[ii] * ia * rc) * (
                                                1 - ia + ia * rc)
                                    m4 = round(m4, 5)

                                    m5 = belta[jj][ii] * S1 * (1 - ia) * Rgs[ii] * ia * rc / (1 - Rgs[ii] * ia * rc) * (
                                                1 + Rgs[ii] * (1 - ia))
                                    m5 = round(m5, 5)

                                    m6 = belta[jj][ii] * (1 - S1) * Rgs[ii] * ia * rc / (1 - Rgs[ii] * ia * rc) * (
                                                1 - ia + ia * rc)
                                    m6 = round(m6, 5)

                                    Tm = (m1 + m2 + m3 + m4 + m5 + m6) / np.pi  # 总的多次散射项
                                    Rm[jj][ii] = round(Tm, 6)

                            lc = 0  # 待求lc
                            lci = 0  # 待求lci
                            lcv = 0  # 待求lcv

                            # 地表天空散射光占比
                            def fx(thx):
                                return np.cos(thx) * np.exp(-G * lamtai * 2 / np.cos(thx)) * np.sin(thx)

                            S, _ = integrate.quad(fx, 0, np.pi / 2)

                            # 确定俩椭圆位置 形状
                            r = np.sqrt(lamtai * G * LAIa / np.pi)  # 假设天顶角为0时，为圆形，半径r
                            lv = G * lamtai * LAIa / np.cos(m_v / 180 * np.pi)  # 观测方向无因次拦截因子
                            li = G * lamtai * LAIa / np.cos(m_s / 180 * np.pi)  # 入射方向无因次拦截因子
                            ui = np.cos(m_s / 180 * np.pi)
                            uv = np.cos(m_v / 180 * np.pi)
                            ul = LAIa / Hb  # 体密度函数
                            # 以长半轴为x轴方向
                            ci = 0.5 * (Ha + Hb) * ul * np.tan(m_s / 180 * np.pi)  # 入射方向椭圆中心在x轴位置
                            cv = 0.5 * (Ha + Hb) * ul * np.tan(m_v / 180 * np.pi)  # 观测方向椭圆中心在x轴位置

                            # 注意这里的坐标系不一定是同一个，都是相对各自椭圆长半轴方向的坐标系
                            # 短半轴仍为r
                            rli = r / ui  # 入射方向椭圆长半轴
                            rlv = r / uv  # 观测方向椭圆长半轴
                            x = sympy.Symbol('x', real=True)
                            y = sympy.Symbol('y', real=True)
                            solved_value = sympy.solve([(x - ci) ** 2 / rli ** 2 + y ** 2 / r ** 2 - 1,
                                                        (x * sympy.cos(vv_ss / 180 * np.pi) + y * sympy.sin(
                                                            vv_ss / 180 * np.pi) - cv) ** 2 / rlv ** 2 + (
                                                                    y * sympy.cos(vv_ss / 180 * np.pi) - x * sympy.sin(
                                                                vv_ss / 180 * np.pi)) ** 2 / r ** 2 - 1],
                                                       [x, y])

                            # 计算相交面积lc
                            # 如果两个椭圆完全重合，那么数据类型是字典
                            # 完全不重合，返回空列表，判断A中心是否在B内 || B中心是否在A内 判断是分离还是包含关系
                            # 普遍情况下，有两个交点，利用海伦公式和椭圆扇形面积公式求取lc
                            # 两个交点的坐标系是i椭圆坐标系 俩交点 (ci, 0) (vtoix, vtoiy)
                            vtoix = cv * np.cos(vv_ss / 180 * np.pi)
                            vtoiy = -cv * np.sin(vv_ss / 180 * np.pi)  # 将v椭圆中心转换为i椭圆坐标系
                            itovx = ci * np.cos(vv_ss / 180 * np.pi)
                            itovy = -ci * np.sin(vv_ss / 180 * np.pi)  # 将i椭圆中心转换为v椭圆坐标系

                            if len(solved_value) == 0:
                                if math.pow((vtoix - ci), 2) / math.pow(rli, 2) + math.pow(vtoiy, 2) / math.pow(r,
                                                                                                                2) - 1 < 0 or math.pow(
                                        (itovx - cv), 2) / math.pow(rlv, 2) + math.pow(itovy, 2) / math.pow(r,
                                                                                                            2) - 1 < 0:
                                    # lc = (lv * np.cos(m_v / 180 * np.pi)) if min(lv, li)==lv else (li * np.cos(m_s / 180 * np.pi))
                                    lc = min(li, lv)
                                    # print('情况1')
                                else:
                                    lc = 0
                                    # print('情况2')
                            elif (str(type(
                                    solved_value[0][0].evalf())) == '<class \'sympy.core.numbers.Float\'>') == False:
                                lc = lv
                                # print('情况3')
                            else:
                                # print('情况4')
                                x1 = float(str((solved_value[0][0].evalf())))
                                y1 = float(str((solved_value[0][1].evalf())))
                                x2 = float(str((solved_value[1][0].evalf())))
                                y2 = float(str((solved_value[1][1].evalf())))
                                # 求i椭圆扇形面积
                                lci_sector = 0.5 * rli * r * np.arccos(
                                    (r * r * (x1 - ci) * (x2 - ci) + rli * rli * y1 * y2) / (r * r * rli * rli))
                                along = np.sqrt(math.pow(x1 - x2, 2) + math.pow(y1 - y2, 2))
                                blong = np.sqrt(math.pow(x1 - ci, 2) + math.pow(y1, 2))
                                clong = np.sqrt(math.pow(x2 - ci, 2) + math.pow(y2, 2))
                                p_helen = (along + blong + clong) / 2
                                s_helen = np.sqrt(p_helen * (p_helen - along) * (p_helen - blong) * (p_helen - clong))
                                lci = lci_sector - s_helen
                                # 求v椭圆扇形面积
                                # (cv, 0)
                                # 将两个交点转为v坐标
                                x1tov = x1 * np.cos(vv_ss / 180 * np.pi) + y1 * np.sin(vv_ss / 180 * np.pi)
                                y1tov = y1 * np.cos(vv_ss / 180 * np.pi) - x1 * np.sin(vv_ss / 180 * np.pi)
                                x2tov = x2 * np.cos(vv_ss / 180 * np.pi) + y2 * np.sin(vv_ss / 180 * np.pi)
                                y2tov = y2 * np.cos(vv_ss / 180 * np.pi) - x2 * np.sin(vv_ss / 180 * np.pi)
                                lcv_sector = 0.5 * rlv * r * np.arccos(
                                    (r * r * (x1tov - cv) * (x2tov - cv) + rlv * rlv * y1tov * y2tov) / (
                                                r * r * rlv * rlv))
                                alongv = np.sqrt(math.pow(x1tov - x2tov, 2) + math.pow(y1tov - y2tov, 2))
                                blongv = np.sqrt(math.pow(x1tov - cv, 2) + math.pow(y1tov, 2))
                                clongv = np.sqrt(math.pow(x2tov - cv, 2) + math.pow(y2tov, 2))
                                p_helenv = (alongv + blongv + clongv) / 2
                                s_helenv = np.sqrt(
                                    p_helenv * (p_helenv - alongv) * (p_helenv - blongv) * (p_helenv - clongv))
                                lcv = lcv_sector - s_helenv
                                # lc = lci * np.cos(m_s / 180 * np.pi) + lcv * np.cos(m_v / 180 * np.pi)  # lc面积
                                lc = lci + lcv
                                # lc = lc / (np.pi * r * rlv)

                            # 求四分量kg, kc, kt, kz
                            kg = np.exp(-li - lv + lc)
                            kz = np.exp(-lv) - kg
                            kc = 1 - np.exp(-lv) - np.exp(-li) + np.exp(-li - lv + lc)
                            kt = 1 - np.exp(-lv) - kc
                            kcl = kc * leaf_ratio0
                            kcnl = kc * nonleaf_ratio0
                            kcnl_1 = kc - kcl - kcnl + (1 - kc - kg - kz) * (1 - leaf_ratio0)

                            # 装盘
                            k_ratio[count][0] = LAIa
                            k_ratio[count][1] = s_angle_min + j * s_interval
                            k_ratio[count][2] = ss_j * ss_interval + ss_min
                            k_ratio[count][3] = v_angle_min + i * v_interval
                            k_ratio[count][4] = vv_i * vv_interval + vv_min
                            k_ratio[count][5] = kc
                            k_ratio[count][6] = kt
                            k_ratio[count][7] = kg
                            k_ratio[count][8] = kz
                            tem1 = (s_number * i + j) % s_number
                            for k in range(Bands):
                                k_ratio[count][9 + k] = (1 - belta[tem1][k]) * (kcl * Rcs[k] + kg * Rgs[k]) + \
                                                        (1 - S) * belta[tem1][k] * kt * Rcs[k] + S * belta[tem1][
                                                            k] * kz * Rgs[k] + \
                                                        (1 - belta[tem1][k]) * kcnl * Rncs[k] + (1 - S) * belta[tem1][
                                                            k] * kcnl_1 * Rncs[k] + Rm[tem1][k]
                            k_ratio[count][12] = lamtai
                            count = count + 1
                            print('Finish! ' + str(count))

            writer = pd.ExcelWriter('./res/' + str(round(LAIa, 2)) + '.xlsx')
            data_df = pd.DataFrame(k_ratio)
            data_df.to_excel(writer, 'page_1', float_format='%.5f', index=False, header=False)
            writer.save()