from cProfile import label
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import xlrd

'''
绘制3.1节所需图片
S8 0-2.5水汽 0°观测天顶角 透过率 上行辐射 下行辐射 散点图
'''
excelFile = '../data/res.xlsx'
# 读excel
content = xlrd.open_workbook(excelFile)
sh = content.sheet_by_name('page_1')

# 找到水汽和S8 0 透过率 上行辐射 下行辐射 
# trans_0_S8    up_0_S8    down_0_S8
tranCol = 'trans_0_S8'
upCol = 'up_0_S8'
downCol = 'down_0_S8'
wvCol = 'waterVapor'
for i in range(sh.ncols):
    if tranCol == sh.cell(0, i).value:
        tranLine = i
    if upCol == sh.cell(0, i).value:
        upLine = i
    if downCol == sh.cell(0, i).value:
        downLine = i
    if wvCol == sh.cell(0, i).value:
        wvLine = i
        
# 现在知道具体第几列了
# 提取数据
tranLst = []
upLst = []
downLst = []
wvLst = []
for i in range(1, sh.nrows):
    if float(sh.cell(i, wvLine).value) < 2.5:
        tranLst.append(float(sh.cell(i, tranLine).value))
        upLst.append(float(sh.cell(i, upLine).value))
        downLst.append(float(sh.cell(i, downLine).value))
        wvLst.append(float(sh.cell(i, wvLine).value))

# 
plt.rcParams['font.sans-serif'] = ['SimHei']  # 黑体
plt.rcParams['axes.unicode_minus'] = False    # 解决无法显示符号的问题
# 坐标轴的刻度设置向内(in)或向外(out)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
sns.set(style="white", font='SimHei')

plt.grid()
# plt.xlim(xmin=0, xmax=2.5)
# plt.ylim(ymin=0.6, ymax=1.0)
# plt.xlabel('水汽含量(g/cm2)')
# plt.ylabel('S8波段0°透过率')
# plt.scatter(wvLst, tranLst, s=1, c='r', marker='o')
# plt.show()

plt.xlim(xmin=0, xmax=2.5)
plt.ylim(ymin=0, ymax=2.5)
plt.xlabel('水汽含量(g/cm2)')
plt.ylabel('S8波段0°下行辐射(W/m2 sr)')
plt.scatter(wvLst, downLst, s=1, c='b', marker='o', label='下行辐射')
plt.show()