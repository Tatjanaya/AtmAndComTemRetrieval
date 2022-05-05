from csv import excel_tab
import createTape5
import runMODTRAN
import extractChnData
import drawLineChartWvAndMid
import drawLineChart

# angleLst = [0, 10, 20, 30, 40, 50, 55]
angleLst = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55]

fileNameInput = 'TIGR_atmProfilesSelection.txt'
fileNameOutput = './tape5/trans_10.tp5'
fltPath = 'D:\Sentinel3_SLSTR_S8_S9.flt'

file_dir = r'E:\integrate\modtran\tape5'
mod5_dir = r'E:\1209\modtran5.2\MODTRAN'

excelDir = r'E:\integrate\modtran\result'
chnDir = r'E:\integrate\modtran\tape5'

excelFile = r'E:\integrate\modtran\result\res.xlsx'
duibiexcel = r'E:/integrate/data/廓线-水汽对照表.xlsx'

weidu = 'mid'
wv = [0, 2.5]

# kind = 1
# angle = 10
# 创建Tape5 文件 指定透过率前缀tran 上行前缀up 下行前缀down
# for angle in angleLst:
#     createTape5.tape5Generate(fileNameInput, './tape5/trans_' + str(angle) + '.tp5', fltPath, 1, angle)
#     createTape5.tape5Generate(fileNameInput, './tape5/up_' + str(angle) + '.tp5', fltPath, 2, angle)
#     createTape5.tape5Generate(fileNameInput, './tape5/down_' + str(angle) + '.tp5', fltPath, 3, angle)

# Run MODTRAN
# runMODTRAN.runMODTRAN(file_dir, mod5_dir)

# 提取chn里必要数据
# extractChnData.extractChnData(excelDir, chnDir, fileNameInput)
# waterVapor
# 绘图
trans_up_down_lst = ['_0_S8', '_0_S9', '_5_S8', '_5_S9', '_10_S8', '_10_S9', '_15_S8', '_15_S9', '_20_S8', '_20_S9', '_25_S8', '_25_S9', \
    '_30_S8', '_30_S9', '_35_S8', '_35_S9', '_40_S8', '_40_S9', '_45_S8', '_45_S9', '_50_S8', '_50_S9', '_55_S8', '_55_S9']

ang_lst = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55]

# for ang in ang_lst:
#     for tud in trans_up_down_lst:
#         drawLineChartWvAndMid.drawLineChart(excelFile, 'trans_' + str(ang) + '_S8', 'trans' + tud, weidu, wv, duibiexcel)
#         print(str(ang) + '----' + 'trans' + tud)
# for ang in ang_lst: 
#     for tud in trans_up_down_lst:   
#         drawLineChartWvAndMid.drawLineChart(excelFile, 'up_' + str(ang) + '_S8', 'up' + tud, weidu, wv, duibiexcel)
#         print(str(ang) + '----' + 'up' + tud)
# for ang in ang_lst:
#     drawLineChartWvAndMid.drawLineChart(excelFile, 'trans_' + str(ang) + '_S8', 'up_' + str(ang) + '_S8', weidu, wv, duibiexcel)
    # print(str(ang) + '----' + 'down' + tud)

# drawLineChart.drawLineChart(excelFile, 'trans_0_S8', 'trans_55_S8')
# drawLineChart.drawLineChart(excelFile, 'trans_0_S8', 'trans_0_S9')
# drawLineChart.drawLineChart(excelFile, 'trans_0_S8', 'trans_55_S9')
drawLineChart.drawLineChart(excelFile, 'up_0_S8', 'down_0_S8')
drawLineChart.drawLineChart(excelFile, 'up_0_S8', 'down_55_S8')
drawLineChart.drawLineChart(excelFile, 'up_0_S8', 'down_0_S9')
drawLineChart.drawLineChart(excelFile, 'up_0_S8', 'down_55_S9')

# for ang in ang_lst:
#     drawLineChartWvAndMid.waterVaporDraw(excelFile, 'trans_' + str(ang) + '_S8', 0, 2.5)
#     drawLineChartWvAndMid.waterVaporDraw(excelFile, 'trans_' + str(ang) + '_S9', 0, 2.5)
# for ang in ang_lst:
#     drawLineChartWvAndMid.waterVaporDraw(excelFile, 'up_' + str(ang) + '_S8', 0, 2.5)
#     drawLineChartWvAndMid.waterVaporDraw(excelFile, 'up_' + str(ang) + '_S9', 0, 2.5)
# for ang in ang_lst:
#     drawLineChartWvAndMid.waterVaporDraw(excelFile, 'down_' + str(ang) + '_S8', 0, 2.5)
#     drawLineChartWvAndMid.waterVaporDraw(excelFile, 'down_' + str(ang) + '_S9', 0, 2.5)