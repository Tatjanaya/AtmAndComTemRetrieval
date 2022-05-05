""" 定义生成Tape5文件

@ 参数设置:
fileNameInput: TIGR大气廓线文件
fileNameOutput: 输出Tape5文件位置
fltPath: 光谱响应函数文件绝对路径
kind: 大气参数类型
    kind = 1: 透过率
    kind = 2: 大气上行辐射
    kind = 3: 大气下行辐射
angle: 天顶角

@ TIGR大气廓线文件概览:
1-236-热带 
237-375-中纬度夏季 
376-494-中纬度冬季 
495-573-极地夏季
574-946-极地冬季
"""
def tape5Generate(fileNameInput, fileNameOutput, fltPath, kind, angle):
    for i in range(946):
        # card 5 IRPT
        irpt = 1
        if i == 945:
            irpt = 0
        else:
            irpt = 1
        # 文件读取位置
        start = i * 43
        # 高程 压强 气温 水汽含量 二氧化碳（0）臭氧 
        altitude = []
        pressure = []
        temperature = []
        waterVapor = []
        ozone = []
        
        with open(fileNameInput, 'r') as f:
            # 判断季节，默认热带
            season = 1
            if i < 235:
                season = 1
            elif i < 375:
                season = 2
            elif i < 494:
                season = 3
            elif i < 573:
                season = 4
            else:
                season = 5
            
            # 季节字符串
            seasonStr = str(season) + '    '
            seasonStrLink = seasonStr + seasonStr + seasonStr + seasonStr
            
            # AAH Cxxx
            AAHSeason = str(season + 1)
            AAHSeasonStr = 'AAH C' + AAHSeason + AAHSeason + AAHSeason
            
            # 依次读取各条廓线高程 压强 气温 水汽含量 二氧化碳（0）臭氧 
            for line in f.readlines()[start + 2: start + 42]:
                curline = line.strip().split()
                floatline = list(map(float, curline))
                altitude.append(floatline[0])
                pressure.append(floatline[1])
                temperature.append(floatline[2])
                waterVapor.append(floatline[3])
                ozone.append(floatline[4])
        
        with open(fileNameOutput, 'a') as f:
            if kind == 1:
                f.write('TMF 7    2    0   -1    0    0    ' + seasonStrLink + '0    1    0   0.001    0.0   !card1' + '\n' \
                    + 'FFF  8 0.0   365.000                    01 F T T       0.000      0.00     0.000     0.000     0.000         0   !card1a' + '\n' \
                        + fltPath + '\n' \
                            + '    1    0    0    0    0    0    23.000     0.000     0.000     0.000     0.500   !card2' + '\n' \
                                '   40    0    0                           0.0    0     1.000    28.964  !card2c' + '\n')
                for j in range(len(altitude)):
                    f.write("   " + str(format(altitude[j], '.3f')).rjust(7) + " " + str(format(pressure[j], '.3e')).lower() + " " \
                        + str(format(temperature[j], '.3e')).lower() + " " + str(format(waterVapor[j], '.3e')).lower() + " " + str(format(0.000, '.3e')).lower() + " " + str(format(ozone[j], '.3e')).lower() + AAHSeasonStr + "\n")
                f.write("   100.000     0.002   " + str(format(180 - angle, '.3f')) + "     0.000     0.000     0.000    0          0.000 !card3" + "\n" \
                    + "   737.000  1000.000       1.0       2.0RM              F  1             !card4" + "\n" \
                        + "    " + str(irpt) + " !card5" + "\n")
            elif kind == 2:
                f.write('TMF 7    2    1   -1    0    0    ' + seasonStrLink + '0    1    0   0.001    0.0   !card1' + '\n' \
                    + 'FFF  8 0.0   365.000                    01 F T T       0.000      0.00     0.000     0.000     0.000         0   !card1a' + '\n' \
                        + fltPath + '\n' \
                            + '    1    0    0    0    0    0    23.000     0.000     0.000     0.000     0.500   !card2' + '\n' \
                                '   40    0    0                           0.0    0     1.000    28.964  !card2c' + '\n')
                for j in range(len(altitude)):
                    f.write("   " + str(format(altitude[j], '.3f')).rjust(7) + " " + str(format(pressure[j], '.3e')).lower() + " " \
                        + str(format(temperature[j], '.3e')).lower() + " " + str(format(waterVapor[j], '.3e')).lower() + " " + str(format(0.000, '.3e')).lower() + " " + str(format(ozone[j], '.3e')).lower() + AAHSeasonStr + "\n")
                f.write("   100.000     0.002   " + str(format(180 - angle, '.3f')) + "     0.000     0.000     0.000    0          0.000 !card3" + "\n" \
                    + "   737.000  1000.000       1.0       2.0RM              F  1             !card4" + "\n" \
                        + "    " + str(irpt) + " !card5" + "\n")
            else:
                f.write('TMF 7    2    1   -1    0    0    ' + seasonStrLink + '0    1    0   0.001    0.0   !card1' + '\n' \
                    + 'FFF  8 0.0   365.000                    01 F T T       0.000      0.00     0.000     0.000     0.000         0   !card1a' + '\n' \
                        + fltPath + '\n' \
                            + '    1    0    0    0    0    0    23.000     0.000     0.000     0.000     0.500   !card2' + '\n' \
                                '   40    0    0                           0.0    0     1.000    28.964  !card2c' + '\n')
                for j in range(len(altitude)):
                    f.write("   " + str(format(altitude[j], '.3f')).rjust(7) + " " + str(format(pressure[j], '.3e')).lower() + " " \
                        + str(format(temperature[j], '.3e')).lower() + " " + str(format(waterVapor[j], '.3e')).lower() + " " + str(format(0.000, '.3e')).lower() + " " + str(format(ozone[j], '.3e')).lower() + AAHSeasonStr + "\n")
                f.write("     0.002   100.000    " + str(format(angle, '.3f')) + "     0.000     0.000     0.000    0          0.000 !card3" + "\n" \
                    + "   737.000  1000.000       1.0       2.0RM              F  1             !card4" + "\n" \
                        + "    " + str(irpt) + " !card5" + "\n")