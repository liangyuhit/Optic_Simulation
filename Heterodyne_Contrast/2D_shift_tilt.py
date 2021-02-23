# -*- coding: utf-8 -*-
'''
Created on Feb 22, 2021

@author: yl
'''
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
plt.rcParams['font.sans-serif']=['SimSun'] #用来正常显示中文标签

j = complex(0, 1)
c = 3e8 # 光速 [m/s]
Lamda = 633e-9 # 光波长 [m]
f_1 = c / Lamda # 光频率 [Hz]
f_2 = f_1 + 1e7
w_1 = 2*np.pi*f_1
w_2 = 2*np.pi*f_2

Fs = 100 * f_1 # 数据点频率 [Hz]
T0 = 1 / Fs  # 数据点间隔 [s]
N = int(1e6) # 数据点数
T = N * T0 # 数据时长 [s]
timeline = np.arange(N) * T0 # 时间序列

w_0 = 3e-3
Z_R = np.pi * w_0**2 / Lamda

def w_z(w_0, Z_R, Z):
    return w_0*np.sqrt(1 + (Z/Z_R)**2)

def R_z(Z_R, Z):
    return Z*(1+(Z_R/Z)**2)

def Phi_Gouy(Z_R, Z):
    return np.arctan(Z/Z_R)    

def f(X, Y, L=0.5, V_m_x=0e-5, V_m_y=0, V_r_x=0, V_r_y=0, M_m=0e-3, N_m=0, M_r=0, N_r=0, w_0=w_0, Z_R=Z_R,A_mea=1,A_ref=1,t=0):
    
    k = 2 * np.pi / Lamda
    I_0 = 1
        
    d_mea = (X*V_m_x + Y*V_m_y - L*V_m_x**2 - L*V_m_y**2 - M_m*V_m_x - N_m*V_m_y) / np.sqrt(V_m_x**2 + V_m_y**2 + 1)
    d_ref = (X*V_r_x + Y*V_r_y - L*V_r_x**2 - L*V_r_y**2 - M_r*V_r_x - N_r*V_r_y) / np.sqrt(V_r_x**2 + V_r_y**2 + 1)
    
    r_mea = np.sqrt( (X-M_m)**2 + (Y-N_m)**2 + (L)**2 - ((X-M_m)*V_m_x + (Y-N_m)*V_m_y + L)**2 / (V_m_x**2 + V_m_y**2 +1) )
    r_ref = np.sqrt( (X-M_r)**2 + (Y-N_r)**2 + (L)**2 - ((X-M_r)*V_r_x + (Y-N_r)*V_r_y + L)**2 / (V_r_x**2 + V_r_y**2 +1) )
    
    Z_r = L * np.sqrt(1+V_r_x**2+V_r_y**2) + d_ref
    Z_m = L * np.sqrt(1+V_m_x**2+V_m_y**2) + d_mea
    
    R_ref = R_z(Z_R, Z_r)
    R_mea = R_z(Z_R, Z_m)
    
    Phi_gouy_ref = Phi_Gouy(Z_R, Z_r)
    Phi_gouy_mea = Phi_Gouy(Z_R, Z_m)
    
    E_m_0 = A_mea*np.cos(w_1*t)
    E_r_0 = A_ref*np.cos(w_2*t)
    w_z_m = w_z(w_0, Z_R, Z_m)
    w_z_r = w_z(w_0, Z_R, Z_r)
    
    E_m = A_mea * w_0/w_z_m * np.exp(-(r_mea**2)/(w_z_m**2)) * np.cos(w_1*t + k*Z_m + k*r_mea**2/2/R_mea + Phi_gouy_mea)
    E_r = A_ref * w_0/w_z_r * np.exp(-(r_ref**2)/(w_z_r**2)) * np.cos(w_2*t + k*Z_r + k*r_ref**2/2/R_ref + Phi_gouy_ref)
    
    A_m = A_mea * w_0/w_z_m * np.exp(-(r_mea**2)/(w_z_m**2))
    A_r = A_ref * w_0/w_z_r * np.exp(-(r_ref**2)/(w_z_r**2))
    
    I = A_m**2 + A_r**2 + 2*A_m*A_r*np.cos(2*np.pi*5e6*t*1e8 + k*(Z_m-Z_r) + k*(r_mea**2/2/R_mea-r_ref**2/2/R_ref) + (Phi_gouy_mea-Phi_gouy_ref))
    return np.sum(I)


screen_diameter = 10e-3 # in m
dx = np.linspace(-screen_diameter/2, screen_diameter/2, num=101)
dy = dx
X, Y = np.meshgrid(dx, dy)

'''
    光斑平移的影响，平行光，不同光斑尺寸
'''

M_m_set = np.linspace(start=-10e-3, stop=10e-3, num=301)
V_m_x_set = np.linspace(start=-200e-6, stop=200e-6, num=301)
w_0_set = [1.5e-3, 3e-3, 4.5e-3]
Z_R_set = [np.pi * i**2 / Lamda for i in w_0_set]

contrast_set = np.empty(shape=(len(w_0_set),len(M_m_set),len(V_m_x_set)))

for i in range(len(w_0_set)):
    w_0 = w_0_set[i]
    Z_R = Z_R_set[i]
    for j in tqdm(range(len(M_m_set))):
        for p in range(len(V_m_x_set)):
            result = []
            for k in range(100):
                a = f(X,Y,M_m=M_m_set[j],V_m_x=V_m_x_set[p],w_0=w_0,Z_R=Z_R,t=timeline[k])
                result.append(a)
            con = (np.max(result) - np.min(result))/np.max(result)
            contrast_set[i,j,p] = con

'''
    Writting
''' 
with open('2D_shift_tilt.npy', 'wb') as f:
    np.savez(f,M_m_set=M_m_set,V_m_x_set=V_m_x_set, w_0_set=w_0_set, contrast_set=contrast_set)
# '''
#     Plot
# '''
# colorset = ['black','blue','red','green']
# for i in range(len(w_0_set)):
#     plt.plot(M_m_set*1e3, contrast_set[i], color=colorset[i],label='光斑直径%i mm'%(w_0_set[i]*1e3*2))    
# plt.ylabel('对比度')
# plt.xlabel('平移量/mm')
# plt.legend(loc='upper left')
# plt.show()
