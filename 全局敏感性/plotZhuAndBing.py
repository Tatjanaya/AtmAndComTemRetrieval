import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 
import re

def main(salibFile):
    markLst, dataLst = readFile(salibFile)
    drawBing(dataLst, markLst)
    drawZhu(dataLst)
    
def drawZhu(dataLst):
    # 设置中文画图
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 中文字体设置-黑体
    plt.rcParams['axes.unicode_minus'] = False  # 解决保存图像是负号'-'显示为方块的问题
    sns.set(style="white", font='SimHei')
    
    labels = ['0°', '5°', '10°', '15°', '20°', '25°', '30°', '35°', '40°', '45°', '50°', '55°']
    
    s8_lai = []
    s8_es = []
    s8_ev = []
    s8_Tv = []
    s8_Ts = []
    s8_wv = []
    
    s9_lai = []
    s9_es = []
    s9_ev = []
    s9_Tv = []
    s9_Ts = []
    s9_wv = []
    
    flag = True
    for data in dataLst:
        if flag:
            s8_lai.append(data[0])
            s8_es.append(data[1])
            s8_ev.append(data[2])
            s8_Tv.append(data[3])
            s8_Ts.append(data[4])
            s8_wv.append(data[5])
        else:
            s9_lai.append(data[0])
            s9_es.append(data[1])
            s9_ev.append(data[2])
            s9_Tv.append(data[3])
            s9_Ts.append(data[4])
            s9_wv.append(data[5])
        flag = not flag
    
    s8_data = [s8_lai, s8_es, s8_ev, s8_Tv, s8_Ts, s8_wv]
    s9_data = [s9_lai, s9_es, s9_ev, s9_Tv, s9_Ts, s9_wv]
    
    create_multi_bars(labels, s8_data, 'S8')
    create_multi_bars(labels, s9_data, 'S9')
    
def create_multi_bars(labels, datas, bandMark, tick_step=10, group_gap=2, bar_gap=0):
    '''
    labels : x轴坐标标签序列
    datas :数据集,二维列表,要求列表每个元素的长度必须与labels的长度一致
    tick_step :默认x轴刻度步长为1,通过tick_step可调整x轴刻度步长。
    group_gap : 柱子组与组之间的间隙,最好为正值,否则组与组之间重叠
    bar_gap :每组柱子之间的空隙,默认为0,每组柱子紧挨,正值每组柱子之间有间隙,负值每组柱子之间重叠
    '''
    # 画布大小
    # plt.figure(dpi=300)
    plt.subplots_adjust(left=0.12, bottom=0.11, right=0.90, top=0.88, wspace=0.20, hspace=0.20)
    # ticks为x轴刻度
    ticks = np.arange(len(labels)) * tick_step
    # group_num为数据的组数,即每组柱子的柱子个数
    group_num = len(datas)
    # group_width为每组柱子的总宽度,group_gap 为柱子组与组之间的间隙。
    group_width = tick_step - group_gap
    # bar_span为每组柱子之间在x轴上的距离,即柱子宽度和间隙的总和
    bar_span = group_width / group_num
    # bar_width为每个柱子的实际宽度
    bar_width = bar_span - bar_gap
    # baseline_x为每组柱子第一个柱子的基准x轴位置,随后的柱子依次递增bar_span即可
    baseline_x = ticks - (group_width - bar_span) / 2
    data_name = ['LAI', 'εs', 'εv', 'Tv', 'Ts', 'wv']
    temp_count = 0
    for index, y in enumerate(datas):
        plt.bar(baseline_x + index * bar_span, y, bar_width, label=data_name[temp_count])
        temp_count += 1
    plt.ylabel('一阶敏感度')
    plt.xlabel('vza')
    plt.legend(bbox_to_anchor=(1, 1),)
    # plt.title('S8波段')
    # x轴刻度标签位置与x轴刻度一致
    plt.xticks(ticks, labels, rotation=70)
    plt.subplots_adjust(right=0.8)
    # plt.show()
    plt.savefig('./res/zhu/' + bandMark + '.png')
    plt.clf()

# 饼图
def drawBing(dataLst, markLst):
    for i in range(len(dataLst)):
        data = dataLst[i]
        mark = markLst[i]
        s1_abs = [np.abs(x) for x in data]
        s1_total = sum(s1_abs)
        s1_abs /= s1_total
        plt.pie(s1_abs,
            labels=['lai', 'εs', 'εv', 'Tv', 'Ts', 'wv'],
            colors=['#e80e09', '#e8a509', '#41e809', '#7c09e8', '#ca09e8', '#09e8e8'],
            explode=(0, 0, 0, 0, 0, 0), # 第几部分突出显示 值越大 离中心越远
            autopct='%.2f%%', # 格式化输出百分比
            )
        plt.savefig('./res/bing/bing_' + mark + '.png')
        plt.clf()

def readFile(salibFile):
    f = open(salibFile, 'r')
    lines = f.readlines()
    f.close()
    
    # markLst 存 str -> vza_band -> 0_S8
    markLst = []
    # dataLst 存 [] 
    dataLst = []
    
    # 0 5 10 15 20 25 30 35 40 45 50 55
    for i in range(24):
        markLine = lines[i * 6].strip()
        vza = markLine.split()[1].strip()
        band = markLine.split()[-1].strip()
        markLst.append(str(vza) + '_' + str(band))
    
        # 正则取数字
        resLst = re.findall(r'-?\d+\.?\d*', lines[i * 6 + 2].strip())
        # 会取到S1的1
        del resLst[0]
        resLst = [float(x) for x in resLst]
        
        dataLst.append(resLst)
    
    return markLst, dataLst

salibFile = './res/salib_res_salib.txt'  
main(salibFile)