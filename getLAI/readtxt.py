import os


def readtxt(path, array, band, flag):  # 以空格为间隔
    f1 = open(path, 'r')
    lines1 = f1.readlines()
    c_row = 0
    for line in lines1:
        if flag == 0:
            clist = line.strip('\n').split(' ')
        else:
            clist = line.strip('\n').split('\t')
        if c_row < band:
            array[c_row] = float(clist[0])
            c_row += 1
    f1.close()
    if c_row != band or c_row == 0:
        print('输入文件有误，请检查后重新输入再运行程序')  # 警告用户
        os._exit(0)
    return array


def read1txt(path, array, band, flag):  # 以空格为间隔
    f1 = open(path, 'r')
    lines1 = f1.readlines()
    c_row = 0
    for line in lines1:
        if flag == 0:
            clist = line.strip('\n').split(' ')
        else:
            clist = line.strip('\n').split('\t')
        if c_row < band:
            array[c_row] = float(clist[1])
            c_row += 1
    f1.close()
    if c_row != band or c_row == 0:
        print('输入文件有误，请检查后重新输入再运行程序')  # 警告用户
        os._exit(0)
    return array

def read2txt(path, array1, array2, band, flag):  # 以空格为间隔
    f1 = open(path, 'r')
    lines1 = f1.readlines()
    c_row = 0
    for line in lines1:
        if flag == 0:
            clist = line.strip('\n').split(' ')
        else:
            clist = line.strip('\n').split('\t')
        if c_row < band:
            array1[c_row] = float(clist[1])
            array2[c_row] = float(clist[2])
            c_row += 1
    f1.close()
    if c_row != band or c_row == 0:
        print('输入文件有误，请检查后重新输入再运行程序')  # 警告用户
        os._exit(0)
    return array1, array2

