import extractWaterVapor
import splitWindowNew
import splitWindow

tigrFile = '../data/TIGR_atmProfilesSelection.txt'
excelFile = '../data/res.xlsx'
lutChart = './result/廓线-水汽对照表.xlsx'

# lowVapor = 0
# highVapor = 2.5

# 提取各条廓线编号和水汽含量
# extractWaterVapor.extractWaterVapor(tigrFile)

# 根据输入的水汽含量范围进行劈窗算法
# splitWindow.splitWindow(excelFile, lutChart, tigrFile, lowVapor, highVapor, 0, 55)

oblique = 55
nadirLst = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
nadirLst2 = [0]
weiduLst = ['trop', 'mid', 'polar']

vaporLst = [(0, 2.5), (2, 3.5), (3, 4.5), (4, 6.5)]
vaporLst2 = [(0, 2.5)]
for nadir in nadirLst:
    for vaporDou in vaporLst:
        lowVapor = vaporDou[0]
        highVapor = vaporDou[1]
        splitWindow.splitWindow(excelFile, lutChart, tigrFile, lowVapor, highVapor, nadir, oblique)

# for nadir in nadirLst:
#     for vaporDou in vaporLst:
#         for weidu in weiduLst:
#             lowVapor = vaporDou[0]
#             highVapor = vaporDou[1]
#             splitWindowNew.splitWindow(excelFile, lutChart, tigrFile, lowVapor, highVapor, nadir, oblique, weidu)
        
# 组分温度初值获取
