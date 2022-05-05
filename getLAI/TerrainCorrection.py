from math import *
from PIL import Image
import numpy as np
import numba as nb
import matplotlib.pyplot as pyplot
from osgeo import gdal
import writeoneRaster
import readtiff

ColNum=3666
RowNum=2923
Left=507
Right=3663
Up=0
Down=2923
pix_size=20



rou_x_green=[[0.00001 for i in range(Right-Left)] for j in range(RowNum)] #纠正后的反射率
rou_x_red=[[0.00001 for i in range(Right-Left)] for j in range(RowNum)]
rou_x_nir=[[0.00001 for i in range(Right-Left)] for j in range(RowNum)]
rou_tol=[rou_x_green, rou_x_red, rou_x_nir] #大气层顶反射率
L_top=[[[0.00001 for i in range(Right-Left)] for j in range(RowNum)] for k in range(3)]
V=[] #天空可视因子
Ea=[[[0 for i in range(ColNum)] for j in range(RowNum)] for k in range(3)] #Ea
Ed=[]
Ef=[]
DEM=[]
Asp_mat=[[0.00001 for i in range(ColNum)] for j in range(RowNum)]
slope_mat=[[0.00001 for i in range(ColNum)] for j in range(RowNum)]
T_mat = [[[0.001,0.001,0.001] for i in range(ColNum)] for j in range(RowNum)]
ShadowFactor=[] #遮蔽因子
cosi_mat=[] #有效太阳入射角
rou_res=[[[1 for i in range(Right-Left)] for j in range(RowNum)] for k in range(3)]#计算结果


S0=1367 #太阳常数
Day=219

ESUN = [1845.7, 1528.5, 970.4]
k0=[0.767, 0.821, 0.883] #天空直射辐射比例
k=[1/k0[0]-1, 1/k0[1]-1, 1/k0[2]-1] #散射与直射之比
tau_green=0 #气溶胶厚度
tau_red=0
tau_nir=0
tau_rayleigh=[0.07,0.04,0.02]
tau_aero = 0.2347
tau_tot=[tau_aero + tau_rayleigh[0], tau_aero + tau_rayleigh[1], tau_aero + tau_rayleigh[2]]
tau_tot1=[0.084,0.0399,0.0023]
tau_tot=[0.063,0.03,0.00174]
trans_tot=[0.93, 0.966, 0.998]
theta_s = 29.29 #太阳入射角
theta_s_arc = theta_s*pi/180
cos_theta_s = cos(theta_s_arc)
tg_green = -0.590*cos_theta_s*cos_theta_s + 0.974*cos_theta_s + 0.527
tg_red = -0.411*cos_theta_s*cos_theta_s + 0.677*cos_theta_s + 0.668
tg_nir = -0.065*cos_theta_s*cos_theta_s + 0.139*cos_theta_s + 0.894
tg_list=[tg_green, tg_red, tg_nir]
Lp=[19.3,10.446,3.327] #程辐射
#T=1 #T(theta_s,0)
#E0= S0*(1+0.034*cos(2*pi*Day/365.25)) #垂直于大气层顶的太阳辐射
#E0=[ESUN[0]*exp(-tau_tot[0]/cos_theta_s), ESUN[1]*exp(-tau_tot[1]/cos_theta_s), ESUN[2]*exp(-tau_tot[2]/cos_theta_s)]
E0=[ESUN[0]*(1+0.034*cos(2*pi*Day/365.25)),ESUN[1]*(1+0.034*cos(2*pi*Day/365.25)),ESUN[2]*(1+0.034*cos(2*pi*Day/365.25))]
Ei=[E0[0]*exp(-tau_tot[0]/cos_theta_s), E0[1]*exp(-tau_tot[1]/cos_theta_s), E0[2]*exp(-tau_tot[2]/cos_theta_s)] #Ei
print('Ei:',Ei)
D = 1.0139 * 1.496e11 #日地距离 单位：米

def cosi_minus_coss(cosi_mat):
    pyplot.imshow(np.array(cosi_mat))
    pyplot.show()
    dif_mat=[[0 for i in range(ColNum)] for j in range(RowNum)]
    sum = 0
    for i in range(RowNum):
        for j in  range(ColNum):
            ele = round(10000*(cosi_mat[i][j]-cos_theta_s))
            dif_mat[i][j]+=ele
            sum += ele
    sum = sum/RowNum/ColNum
    print("SUM=",sum)
    pyplot.imshow(np.array(dif_mat))
    pyplot.show()
    return


def save_rou(path,band):
    with open(path, 'w') as rou_file:
        for i in range(RowNum):
            for j in range(Right-Left-1):
                rou_file.write(str(rou_res[band][i][j])+',')
            rou_file.write(str(rou_res[band][i][-1])+'\n')

def save_Ea(band):
    path="D:\\Ea_"+str(band)+".txt"
    with open(path, 'w') as Ea_file:
        for i in range(RowNum):
            for j in range(ColNum-1):
                Ea_file.write(str(round(Ea[band][i][j], 8)) + ',')
            Ea_file.write(str(round(Ea[band][i][-1], 8)) + '\n')

def read_Ea(band):
    path = "D:\\Ea_" + str(band) + ".txt"
    with open(path,'r') as Ea_file:
        for i in range(RowNum):
            input_lst = Ea_file.readline().rstrip().split(',')
            for j in range(ColNum):
                Ea[0][i][j] = float(input_lst[j])


def save_T(T_matrix):
    with open('D:\\T1.txt', 'w') as T1_file:
        for i in range(RowNum):
            for j in range(ColNum-1):
                T1_file.write(str(round(T_matrix[i][j][0], 8)) + ',')
            T1_file.write(str(round(T_matrix[i][-1][0], 8)) + '\n')
    with open('D:\\T2.txt', 'w') as T2_file:
        for i in range(RowNum):
            for j in range(ColNum-1):
                T2_file.write(str(round(T_matrix[i][j][1], 8)) + ',')
            T2_file.write(str(round(T_matrix[i][-1][1], 8)) + '\n')
    with open('D:\\T3.txt', 'w') as T3_file:
        for i in range(RowNum):
            for j in range(ColNum-1):
                T3_file.write(str(round(T_matrix[i][j][2], 8)) + ',')
            T3_file.write(str(round(T_matrix[i][-1][2], 8)) + '\n')

def read_T():
    with open('D:\\T1.txt' ,'r') as T1_file:
        for i in range(RowNum):
            input_lst = T1_file.readline().rstrip().split(',')
            for j in range(ColNum):
                T_mat[i][j][0] = float(input_lst[j])
    with open('D:\\T2.txt' ,'r') as T2_file:
        for i in range(RowNum):
            input_lst = T2_file.readline().rstrip().split(',')
            for j in range(ColNum):
                T_mat[i][j][1] = float(input_lst[j])
    with open('D:\\T3.txt' ,'r') as T3_file:
        for i in range(RowNum):
            input_lst = T3_file.readline().rstrip().split(',')
            for j in range(ColNum):
                T_mat[i][j][2] = float(input_lst[j])


def read_V(V):
    with open('D:\\V_res.txt' ,'r') as V_file:
        for i in range(RowNum):
            input_lst = V_file.readline().rstrip().split(',')
            V.append([])
            for j in range(ColNum):
                V[i].append(0)
                V[i][-1] += (float(input_lst[j]))

def read_ShadowFactor(ShadowFactor):
    with open('D:\\Theta_res.txt', 'r') as SF_file:
        for i in range(RowNum):
            input_lst = SF_file.readline().rstrip().split(',')
            ShadowFactor.append([])
            for j in range(ColNum):
                ShadowFactor[i].append(0)
                ShadowFactor[i][-1] += (int(input_lst[j]))


def get_cos(vec1,vec2):
    ans=((vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2])/
           (np.sqrt(vec1[0]*vec1[0]+vec1[1]*vec1[1]+vec1[2]*vec1[2])*
               np.sqrt(vec2[0]*vec2[0]+vec2[1]*vec2[1]+vec2[2]*vec2[2])))
    if ans>1:
        ans=1
    else:
        if ans<-1:
            ans=-1
    return (ans)


def Calc_T_mat(T_mat):
    for i in range(RowNum):
        for j in range(ColNum):
            if Asp_mat[i][j]<90:
                asp_temp=(Asp_mat[i][j]+270)*pi/180
            else:
                asp_temp=(Asp_mat[i][j]-90)*pi/180
            slp_temp = slope_mat[i][j]*pi/180
            #asp_temp = Asp_mat[i][j]*pi/180
            T_mat[i][j] = [sin(slp_temp) * cos(asp_temp), sin(slp_temp) * sin(asp_temp), cos(slp_temp)]

def Calc_Ed(Ed,band):
    for i in range(RowNum):
        for j in range(ColNum):
            Ed[i][j] = ShadowFactor[i][j]*cosi_mat[i][j]*Ei[band]

def Calc_Ltop(band):
    for i in range(RowNum):
        for j in range(Left,Right):
            L_top[band][i][j-Left] = rou_tol[band][i][j-Left] * cos_theta_s * ESUN[band] / (pi  * 10000)
            #print('rou_top:',rou_tol[band][i][j-Left])
            #print('L_top:',L_top[band][i][j-Left])

def Calc_Ef(Ef, band):
    for i in range(RowNum):
        for j in range(ColNum):
            Ef[i][j] = V[i][j]*k[band]*Ei[band]

def Calc_Ea(band):
    exp_tau_band=exp(tau_tot[band])
    for i in range(RowNum):
        print(i)
        for j in range(ColNum):
            if j<Left or j>=Right:  #不在区域内
                Ea[band][i][j] = 0
                continue
            #遍历周边八个像元
            #print("x,y:",i,j)
            temp = 0
            for row in range(max(0,i-1),min(i+2,RowNum)):
                for col in range(max(Left,j-1),min(j+2,Right)):
                    if (row==i and col==j):#跳过自己
                        continue
                    vec_MN=[(i-row)*pix_size, (j-col)*pix_size, DEM[row][col]-DEM[i][j]]
                    cos_TM=get_cos(vec_MN,T_mat[row][col])
                    cos_TM=max(0,cos_TM)
                    cos_TN=-(get_cos(vec_MN,T_mat[i][j])) #vec_NM
                    cos_TN=max(0,cos_TN)
                    dS=pix_size*pix_size*sqrt(1 + (T_mat[i][j][0]/T_mat[i][j][2])**2 + (T_mat[i][j][1]/T_mat[i][j][2])**2)
                    temp += exp_tau_band*(L_top[band][row][col-Left]-Lp[band])*cos_TM*cos_TN*dS/\
                            (vec_MN[0]*vec_MN[0] + vec_MN[1]*vec_MN[1] + vec_MN[2]*vec_MN[2])
                    print('vector('+str(i)+','+str(j)+'):','MN',vec_MN,'TM',T_mat[row][col], 'TN',T_mat[i][j])
                    print('result('+str(i)+','+str(j)+'):',L_top[band][i][j-Left],Lp[band],cos_TM,cos_TN)
            Ea[band][i][j] = temp
        #print(Ea[0][i])
    return

def Calc_rou_x(rou_x, band):
    sum_L=0
    sum_rou=0
    path="F:\\s2-20210807\\MakeRas_tif11.tif"
    temp_data = Image.open(path)
    #L2A_mat = np.array(temp_data)
    #print(np.shape(L2A_mat))
    rou_mat=[[0 for i in range(Right-Left)] for j in range(RowNum)]
    for i in range(RowNum):
        for j in range(Left,Right):
            #L_ij = rou_x[i][j]*(Ea[i][j]+Ed[i][j]+Ef[i][j])
            L_ij = L_top[band][i][j-Left]
            #L_x = (L_ij-tg_list[band]*Lp[band])*(cos(theta_s)+k[band])/(ShadowFactor[i][j]*cosi_mat[i][j]+V[i][j]*k[band]+Ea[band][i][j]/Ei[band])\
            #      +tg_list[band]*Lp[band]
            L_x = (L_ij-tg_list[band]*Lp[band])*(cos_theta_s+k[band])/(cosi_mat[i][j]+k[band])+tg_list[band]*Lp[band]
            #Lx_mat[i][j-Left] += L_ij - L_x
            Ed_adj = cos_theta_s*E0[band]*exp(-tau_tot[band]/cos_theta_s)
            Ef_adj= k[band]*Ei[band]
            rou_x[band][i][j-Left]=round(10000*pi*L_x/(Ed_adj+Ef_adj))
            #rou_mat[i][j-Left] += (round(10000*pi*L_x/(Ed_adj+Ef_adj))-L2A_mat[i][j])
            #sum_L+=L_x
            #sum_rou+=10000*pi*L_x/(Ed_adj+Ef_adj)
            #print('x,y:',i,j,'L_top:',L_ij,'L_adj:',L_x,'rou_adj:',round(10000*pi*L_x/(Ed_adj+Ef_adj)))
            #rou_x[band][i][j-Left] = round(10000*L_x/((cos_theta_s+k[band])*Ei[band]))
        print('line',i)
    #print('sum_L=',sum_L/RowNum/(Right-Left))
    #print('sum_rou=',sum_rou/RowNum/(Right-Left))
    #pyplot.imshow(rou_mat)
    #pyplot.show()
    pyplot.imshow(rou_x[band])
    pyplot.show()


#read_Ea(0)
#pyplot.imshow(np.array(Ea[0]))
#pyplot.show()
#exit(1)

bandnum = 3
#读入影像
ds = gdal.Open("E:\\cliped\\b3_ex.tif")
geoTransform = ds.GetGeoTransform()
proj = ds.GetProjection()  # 地图投影信息
temp_data = Image.open("E:\\cliped\\b3_ex.tif")
rou_x_green = np.array(temp_data)
temp_data = Image.open("E:\\cliped\\b4_ex.tif")
rou_x_red = np.array(temp_data)
temp_data = Image.open("E:\\cliped\\b8a_ex.tif")
rou_x_nir = np.array(temp_data)
rou_tol=[rou_x_green, rou_x_red, rou_x_nir]
#读入数据
#DEM,slope,asp,ShadowFactor,cosi,V_mat
DEM_data = Image.open("D:\\BaiduNetdiskDownload\\N040E115_N045E120\\use_match\\ele_match2.tif")
DEM = np.array(DEM_data)
print('DEM:',len(DEM),len(DEM[0]))
print(DEM[0])
slp_data = Image.open("D:\\BaiduNetdiskDownload\\N040E115_N045E120\\use_match\\slp_match2.tif")
#slope_mat = np.array(slp_data)
#print('slp:',len(slope_mat),len(slope_mat[0]))
#print(slope_mat[0])
#Asp_data = Image.open("D:\\BaiduNetdiskDownload\\N040E115_N045E120\\use_match\\asp_match2.tif")
#Asp_mat = np.array(Asp_data)
#print('Asp:',len(Asp_mat),len(Asp_mat[0]))
#print(Asp_mat[0])
cosi_data = Image.open("D:\\BaiduNetdiskDownload\\N040E115_N045E120\\use_match\\cosi_match2.tif")
cosi_mat = np.array(cosi_data)
pyplot.imshow(cosi_mat)
pyplot.show()
print('cosi:',len(cosi_mat),len(cosi_mat[0]))
print(cosi_mat[0])
cosi_minus_coss(cosi_mat)
#read_V(V)
#print('V:',len(V),len(V[0]))
#print(V[0])
#read_ShadowFactor(ShadowFactor)
#print('SF:',len(ShadowFactor),len(ShadowFactor[0]))
#print(ShadowFactor[0])
#Calc_T_mat(T_mat)
#read_T()
#print('T_mat:',len(T_mat),len(T_mat[0]))
#print(T_mat[0])
#save_T(T_mat)
print('terrain data prepared')


#Calc_Ltop(0)
#Calc_Ea(0)
#print('Ea_0 calculating finished')
#exit(1)
#save_Ea(0)
#read_Ea(0)

print('sensor data prepared')
for i in range(0,2):
    Calc_Ltop(i)
    pyplot.imshow(np.array(L_top[i]))
    pyplot.show()
    Calc_rou_x(rou_res,i)
    save_path="E:\\rou_res"+str(i)+".txt"
    outpath="E:\\rou_res"+str(i)+".tif"
    #save_rou(save_path,i)
    writeoneRaster.writeOneRaster(outpath,np.array(rou_res[i]),geoTransform,proj,RowNum,Right-Left)
    pyplot.imshow(np.array(rou_res[i]))
    pyplot.show()
