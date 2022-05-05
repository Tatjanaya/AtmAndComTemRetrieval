import os
import gdal
import ndviCal
import simuls
import lutCalLAI

# simuls.simulss('leaf_test.txt', 'soil_test.txt', Bands=3, v_angle_min=30, v_angle_max=30, v_interval=5, s_angle_min=40, s_angle_max=40, s_interval=5, LAI_min=1, LAI_max=6.1, LAI_interval=0.1, ss_min=13, ss_max=13, ss_interval=5, vv_min=100, vv_max=100, vv_interval=5, Ha=0, Hb=1.6, nsss=0.5, Rsss=2, flag=0)
# tifLoc = 'res.tif'
# saihanba = 'E:\integrate\getLAI\data\saihanba_test.tif'
# saihanba_ndvi = 'E:\integrate\getLAI\data\saihanba_ndvi.tif'
# redBandLoc = 2
# nirBandLoc = 3
# ndviLoc = 'ndvi.tif'
# ndviCal.ndviCal(tifLoc, redBandLoc, nirBandLoc, saihanba_ndvi)

# vegkind = 'xiaomai'
# simuls.simulss('./saihanba/reflect/xiaomai.txt', './saihanba/reflect/soil.txt', Bands=3, v_angle_min=10, v_angle_max=10, v_interval=5, s_angle_min=30, s_angle_max=30, s_interval=5, LAI_min=0.1, LAI_max=6.1, LAI_interval=0.1, ss_min=144, ss_max=144, ss_interval=5, vv_min=105, vv_max=105, vv_interval=5, Ha=0, Hb=1.6, nsss=0.5, Rsss=2, flag=0, vegkind=vegkind)

# vegkind = 'zhangzisong'
# simuls.simulss('./saihanba/reflect/zhangzisong.txt', './saihanba/reflect/soil.txt', Bands=3, v_angle_min=10, v_angle_max=10, v_interval=5, s_angle_min=30, s_angle_max=30, s_interval=5, LAI_min=0.1, LAI_max=6.1, LAI_interval=0.1, ss_min=144, ss_max=144, ss_interval=5, vv_min=105, vv_max=105, vv_interval=5, Ha=0, Hb=1.6, nsss=0.5, Rsss=2, flag=0, vegkind=vegkind)

# vegkind = 'yumi'
# simuls.simulss('./saihanba/reflect/yumi.txt', './saihanba/reflect/soil.txt', Bands=3, v_angle_min=10, v_angle_max=10, v_interval=5, s_angle_min=30, s_angle_max=30, s_interval=5, LAI_min=0.1, LAI_max=6.1, LAI_interval=0.1, ss_min=144, ss_max=144, ss_interval=5, vv_min=105, vv_max=105, vv_interval=5, Ha=0, Hb=1.6, nsss=0.5, Rsss=2, flag=0, vegkind=vegkind)

# vegkind = 'zhangzisong'
# simuls.simulss('./saihanba/reflect/zhangzisong.txt', './saihanba/reflect/soil.txt', Bands=3, v_angle_min=10, v_angle_max=10, v_interval=5, s_angle_min=30, s_angle_max=30, s_interval=5, LAI_min=0.1, LAI_max=6.1, LAI_interval=0.1, ss_min=144, ss_max=144, ss_interval=5, vv_min=105, vv_max=105, vv_interval=5, Ha=0, Hb=1.6, nsss=0.5, Rsss=2, flag=0, vegkind=vegkind)

# vegkind = 'baihua'
# sim_new.simulss('./saihanba/reflect/baihua.txt', './saihanba/reflect/soil.txt', Bands=3, v_angle_min=10, v_angle_max=10, v_interval=5, s_angle_min=30, s_angle_max=30, s_interval=5, LAI_min=6.0, LAI_max=6.1, LAI_interval=0.1, ss_min=144, ss_max=144, ss_interval=5, vv_min=105, vv_max=105, vv_interval=5, Ha=0, Hb=1.6, nsss=0.5, Rsss=2, flag=0, vegkind=vegkind)


# vegkind = 'ejina'
# simuls.simulss('./leaf_tc.txt', './soil_tc.txt', Bands=3, v_angle_min=11, v_angle_max=11, v_interval=5, s_angle_min=36, s_angle_max=36, s_interval=5, LAI_min=0.1, LAI_max=7.1, LAI_interval=0.1, ss_min=131, ss_max=131, ss_interval=5, vv_min=314, vv_max=314, vv_interval=5, Ha=0, Hb=1.6, nsss=0.5, Rsss=2, flag=0, vegkind=vegkind)

# vegkind = 'yumi'
# simuls.simulss('./224run/yumi.txt', './224run/soil.txt', Bands=3, v_angle_min=0, v_angle_max=0, v_interval=5, s_angle_min=0, s_angle_max=0, s_interval=10, LAI_min=3.0, LAI_max=3.0, LAI_interval=0.1, ss_min=144, ss_max=144, ss_interval=5, vv_min=105, vv_max=105, vv_interval=5, Ha=0, Hb=1.6, nsss=0.5, Rsss=2, flag=0, vegkind=vegkind)
# simuls.simulss('./224run/yumi.txt', './224run/soil.txt', Bands=3, v_angle_min=0, v_angle_max=0, v_interval=5, s_angle_min=10, s_angle_max=10, s_interval=10, LAI_min=3.0, LAI_max=3.0, LAI_interval=0.1, ss_min=144, ss_max=144, ss_interval=5, vv_min=105, vv_max=105, vv_interval=5, Ha=0, Hb=1.6, nsss=0.5, Rsss=2, flag=0, vegkind=vegkind)
# simuls.simulss('./224run/yumi.txt', './224run/soil.txt', Bands=3, v_angle_min=0, v_angle_max=0, v_interval=5, s_angle_min=20, s_angle_max=20, s_interval=10, LAI_min=3.0, LAI_max=3.0, LAI_interval=0.1, ss_min=144, ss_max=144, ss_interval=5, vv_min=105, vv_max=105, vv_interval=5, Ha=0, Hb=1.6, nsss=0.5, Rsss=2, flag=0, vegkind=vegkind)
# simuls.simulss('./224run/yumi.txt', './224run/soil.txt', Bands=3, v_angle_min=0, v_angle_max=0, v_interval=5, s_angle_min=30, s_angle_max=30, s_interval=10, LAI_min=3.0, LAI_max=3.0, LAI_interval=0.1, ss_min=144, ss_max=144, ss_interval=5, vv_min=105, vv_max=105, vv_interval=5, Ha=0, Hb=1.6, nsss=0.5, Rsss=2, flag=0, vegkind=vegkind)
# simuls.simulss('./224run/yumi.txt', './224run/soil.txt', Bands=3, v_angle_min=0, v_angle_max=0, v_interval=5, s_angle_min=40, s_angle_max=40, s_interval=10, LAI_min=3.0, LAI_max=3.0, LAI_interval=0.1, ss_min=144, ss_max=144, ss_interval=5, vv_min=105, vv_max=105, vv_interval=5, Ha=0, Hb=1.6, nsss=0.5, Rsss=2, flag=0, vegkind=vegkind)
# simuls.simulss('./224run/yumi.txt', './224run/soil.txt', Bands=3, v_angle_min=0, v_angle_max=0, v_interval=5, s_angle_min=50, s_angle_max=50, s_interval=10, LAI_min=3.0, LAI_max=3.0, LAI_interval=0.1, ss_min=144, ss_max=144, ss_interval=5, vv_min=105, vv_max=105, vv_interval=5, Ha=0, Hb=1.6, nsss=0.5, Rsss=2, flag=0, vegkind=vegkind)
# simuls.simulss('./224run/yumi.txt', './224run/soil.txt', Bands=3, v_angle_min=10, v_angle_max=10, v_interval=5, s_angle_min=50, s_angle_max=50, s_interval=10, LAI_min=6.0, LAI_max=6.0, LAI_interval=0.1, ss_min=144, ss_max=144, ss_interval=5, vv_min=105, vv_max=105, vv_interval=5, Ha=0, Hb=1.6, nsss=0.5, Rsss=2, flag=0, vegkind=vegkind)


# tifLoc = r"E:\integrate\getLAI\FPan\S2B_band6_20210807.tif"
# tarLoc = r'E:\integrate\getLAI\FPan\res\LAI\0807.tif'
# detector = 10000
# lutCalLAI.lutCalLAI(tifLoc, './res/LAI.xlsx', './FPan/LAI_test.tif', 2, 3, 8, detector)

# 获得target_dir下所有tif文件
# target_dir = './FPan/res'
# ex_dir = 'E:\integrate\getLAI\res\LAI_sen2.xlsx'
# fpan_dir = './FPan/res/LAI'
# lst = os.listdir(target_dir)
# fileLst = []
# for i in lst:
#     if os.path.splitext(i)[1] == '.tif':
#         fileLst.append(i)
        
# for file in fileLst:
#     lutCalLAI.lutCalLAI(file, ex_dir, fpan_dir + '/' + file, 2, 3, 8)


vegkind = 'luohenew'
simuls.simulss('./leaf_tc.txt', './soil_tc.txt', Bands=3, v_angle_min=11, v_angle_max=11, v_interval=5, s_angle_min=36, s_angle_max=36, s_interval=5, LAI_min=0.1, LAI_max=7.1, LAI_interval=0.1, ss_min=131, ss_max=131, ss_interval=5, vv_min=314, vv_max=314, vv_interval=5, Ha=0, Hb=1.6, nsss=0.5, Rsss=2, flag=0, vegkind=vegkind)

# tifLoc = r'./额济纳数据/ejina123.tif'
# tarLoc = r'./额济纳数据/res/laiejina.tif'

lutCalLAI.lutCalLAI(r'E:/integrate/黑水/sen2/process/lss.tif', r'E:/integrate/黑水/sen2/ejina.xlsx', r'E:/integrate/黑水/sen2/process/reslai.tif', 0, 1, 2)