# -*- coding: utf-8 -*-
'''
Created on Feb 20, 2021

@author: yl
'''
import numpy as np
import matplotlib.pyplot as plt

j = complex(0, 1)
c = 3e8 # 光速 [m/s]
Lamda = 633e-9 # 光波长 [m]
Fc = c / Lamda # 光频率 [Hz]

w_0 = 3e-3
Z_R = np.pi * w_0**2 / Lamda

def fit_func(x, a, b):
    return a*x+b

def w_z(w_0, Z_R, Z):
    return w_0*np.sqrt(1 + (Z/Z_R)**2)

def R_z(Z_R, Z):
    return Z*(1+(Z_R/Z)**2)

def Phi_Gouy(Z_R, Z):
    return np.arctan(Z/Z_R)    

def f(X, Y, V_r_x=0, V_r_y=0, V_m_x=0, V_m_y=0, M=0, N=0, L=0.12, D=0, w_0=w_0, Z_R=Z_R,):
    
    k = 2 * np.pi / Lamda
    I_0 = 1
    
    Z_p_r = M + N + L * (2/(1-V_r_x**2-V_r_y**2) - 1)
    Z_p_m = M + N + (L+D) * (2/(1-V_m_x**2-V_m_y**2) - 1)
    
    diff_Z_p = D * 2 / (1-V_m_x**2-V_m_y**2) + L * 2 * ((V_m_x**2+V_m_y**2)-(V_r_x**2+V_r_y**2)) / (1-V_m_x**2-V_m_y**2) / (1-V_r_x**2-V_r_y**2)
    
    d_mea = (X*2*V_m_x+Y*2*V_m_y)/(1+V_m_x**2+V_m_y**2) - (L+D)*4*(V_m_x**2+V_m_y**2)/(1-(V_m_x**2+V_m_y**2)**2)
    d_ref = (X*2*V_r_x+Y*2*V_r_y)/(1+V_r_x**2+V_r_y**2) - L*4*(V_r_x**2+V_r_y**2)/(1-(V_r_x**2+V_r_y**2)**2)
    
    r_mea = np.sqrt( X**2 + Y**2 + (L+D)**2 - ((2*X*V_m_x+2*Y*V_m_y+(L+D)*(1-V_m_x**2-V_m_y**2))/(1+V_m_x**2+V_m_y**2))**2 )
    r_ref = np.sqrt( X**2 + Y**2 + (L)**2 - ((2*X*V_r_x+2*Y*V_r_y+L*(1-V_r_x**2-V_r_y**2))/(1+V_r_x**2+V_r_y**2))**2 )
    
    R_ref = R_z(Z_R, Z_p_r+d_ref)
    R_mea = R_z(Z_R, Z_p_m+d_mea)
    
    Phi_gouy_ref = Phi_Gouy(Z_R, Z_p_r+d_ref)
    Phi_gouy_mea = Phi_Gouy(Z_R, Z_p_m+d_mea)
    
    diff_phi = k * (diff_Z_p + (d_mea-d_ref) + r_mea**2/2/R_mea - r_ref**2/2/R_ref) + Phi_gouy_mea - Phi_gouy_ref
    #     A_beat = 0.5 * I_0 * w_0**2 / w_z(w_0, Z_R, Z_p_r+d_ref) / w_z(w_0, Z_R, Z_p_m+d_mea) * np.exp((-X**2-Y**2)/w_z(w_0, Z_R, Z_p_r+d_ref) / w_z(w_0, Z_R, Z_p_m+d_mea))
    I_beat = 1/w_z(w_0, Z_R, Z_p_r+d_ref)**2*np.exp(-2*r_ref**2/w_z(w_0, Z_R, Z_p_r+d_ref)**2) + 1/w_z(w_0, Z_R, Z_p_m+d_mea)**2*np.exp(-2*r_mea**2/w_z(w_0, Z_R, Z_p_m+d_mea)**2) + 2/w_z(w_0, Z_R, Z_p_r+d_ref)/w_z(w_0, Z_R, Z_p_m+d_mea)*np.exp(-r_ref**2/w_z(w_0, Z_R, Z_p_r+d_ref)**2-r_mea**2/w_z(w_0, Z_R, Z_p_m+d_mea)**2)*np.cos(diff_phi)
    I_beat = 0.5*I_0*w_0**2*I_beat
#     I_beat = I_beat.astype(np.int)
    return I_beat

def center_S(V_r_x=0, V_r_y=0, V_m_x=0, V_m_y=0, L=0.12,D=0):
    center_ref = [2*L*V_r_x/(1-V_r_x**2-V_r_y**2),2*L*V_r_y/(1-V_r_x**2-V_r_y**2)]
    center_mea = [2*(L+D)*V_m_x/(1-V_m_x**2-V_m_y**2),2*(L+D)*V_m_y/(1-V_m_x**2-V_m_y**2)]
    print(center_ref, center_mea)
    return center_ref, center_mea

screen_diameter = 10e-3 # in m
dx = np.linspace(-screen_diameter/2, screen_diameter/2, num=1001)
dy = dx
X, Y = np.meshgrid(dx, dy)

img = f(X, Y,V_m_x=1e-3, D=Lamda*0/2)
center = center_S(V_m_x=1e-3,D=Lamda*0/2)

extent=[-screen_diameter/2*1e3,screen_diameter/2*1e3,-screen_diameter/2*1e3,screen_diameter/2*1e3] #in mm
plt.plot(center[0][0]*1000,center[0][1]*1000,marker='o',markersize=6,color='b')
plt.plot(center[1][0]*1000,center[1][1]*1000,marker='x',markersize=6,color='r',mew=2)
plt.imshow(img,origin='lower',extent=extent,cmap='gray', vmin=0, vmax=2)
plt.show()
