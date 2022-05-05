import os
import sys

def updateFile(file,row,num):
    """
    修改文件中第几行的字符串
    :param file:文件名
    :param row:行数
    :param num:新数字
    :return:
    """

    with open(file, "r", encoding="utf-8") as f:
        list = f.readlines()

    for row_num,line in enumerate(list):
        line = line.rstrip()
        if(row_num == int(row)-1):
            line = str(num)
        list[row_num] = line + "\n"

    with open(file, "w", encoding="utf-8") as f:
        f.writelines(list)

    f.close()

def takeValue(file1):
    """
    读取文件1中天空散射光比例的数值
    :param file1:文件1
    :return: p:天空直射光比例
    """
    s = []
    with open(file1, "r", encoding="utf-8") as f:
        list = f.readlines()

    for row_num, line in enumerate(list):
        # line = line.rstrip()
        lsp = []
        if (row_num == 65):
            lsp = line.strip('\n').split(' ')
        for i in lsp:
            s.append(i)
    f.close()
    # print(s[16])
    p = '0'+ s[16]
    return p

def bel(array, s_angle_min, s_interval, band):
    a = array.shape[0]
    b = array.shape[1]
    for i in range(a):
        s_angle = s_angle_min + i * s_interval
        tempTable = runTempTable(s_angle)
        for j in range(b):
            if band[j] <= 533 and band[j] > 439:
                array[i][j] = 1 - float(tempTable[0])
            elif band[j] <= 583 and band[j] > 538:
                array[i][j] = 1 - float(tempTable[1])
            elif band[j] <= 684 and band[j] > 640:
                array[i][j] = 1 - float(tempTable[2])
            elif band[j] <= 714 and band[j] > 695:
                array[i][j] = 1 - float(tempTable[3])
            elif band[j] <= 749 and band[j] > 731:
                array[i][j] = 1 - float(tempTable[4])
            elif band[j] <= 797 and band[j] > 769:
                array[i][j] = 1 - float(tempTable[5])
            elif band[j] <= 881 and band[j] > 847:
                array[i][j] = 1 - float(tempTable[6])
            elif band[j] <= 1682 and band[j] > 1539:
                array[i][j] = 1 - float(tempTable[7])
            elif band[j] <= 2320 and band[j] > 2078:
                array[i][j] = 1 - float(tempTable[8])
            else:
                print('band value is not in interval ' + str(band[j]))
                sys.exit()
    
    return array
                
                

# 返回的是一个太阳天顶角的结果值 res_lst = [各个波段]    
def runTempTable(s_angle):
    # 波段值和范围
    band_lst = [[0.439,0.533,'B2.txt'],[0.538,0.583,'B3.txt'],[0.646,0.684,'B4.txt'],[0.695,0.714,'B5.txt'],[0.731,0.749,'B6.txt'],[0.769,0.797,'B7.txt'],[0.847,0.881,'B8.txt'],[1.539,1.682,'B11.txt'],[2.078,2.320,'B12.txt']]
    result = []
    # step1 更新参数文件
    updateFile(r'B_input.txt', 2, s_angle)
    for n in band_lst:
        updateFile(r'B_input.txt', 14, n[0])
        updateFile(r'B_input.txt', 15, n[1])
        updateFile(r'B_input.txt', 16, n[2])
        for i in range(0, 1):
            if i == 0:
                i = 30
            updateFile(r'B_input.txt', 10, i)
            # step2 运行exe程序
            os.system(r'6S_ATM_corr.EXE<B_input.txt')
            # step3 读取输出文件中天空散射光比例的值
            result.append(takeValue(r'sixs.out'))  
    # 返回result列表
    return result