import os
"""运行一个文件夹下所有tape5文件

1. 将resource文件夹下的所有文件路径写入mod5root.in
2. 调用MODTRAN.exe

@ writeMod5root:
    file_dir: tape5文件夹位置
    mod5_dir: MODTRAN文件夹位置

@ runMODTRAN
    mod5_dir: MODTRAN文件夹位置
"""
def writeMod5root(file_dir, mod5_dir):
    # 文件夹中所有文件
    # 取出这个文件夹下所有'.tp5'文件
    # 如果后缀不为tp5，则从列表中全部删除
    fileLst = []
    for _, _, k in os.walk(file_dir):
        for fileName in k:
            if fileName.endswith('.tp5'):
                fileLst.append(fileName)

    # 写入mod5root.in
    with open(mod5_dir + r'\mod5root.in', 'w', encoding='utf-8') as f:
        for fileName in fileLst:
            tempStr = file_dir + "\\" + fileName
            tempStr.replace('\\\\', '\\')
            f.write(tempStr + "\n")

def runMODTRAN(file_dir, mod5_dir):
    writeMod5root(file_dir, mod5_dir)
    # 切换工作目录
    os.chdir(mod5_dir)
    runExe = "MOD5_win.exe"
    # 执行MODTRAN exe文件
    os.system(runExe)
