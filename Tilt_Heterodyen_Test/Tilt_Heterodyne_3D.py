# -*- coding: utf-8 -*-
 
import matplotlib.pyplot as plt
import numpy as np


j = complex(0, 1)
c = 3e8 # 光速 [m/s]
Lamda = 633e-9 # 光波长 [m]
Fc = c / Lamda # 光频率 [Hz]

screen_diameter = 2e-3
I_0 = 1
A1 = 1 #幅值
A2 = 1
Phase1 = 0 #相位
Phase2 = 0

# dx = np.linspace(screen_diameter/2*(-1), screen_diameter/2, num=301)
dx = np.linspace(-1e-3, 1e-3, num=301)
dy = dx

# V_x, V_y, Vz = 0.002, 0.000, 1 # Reference
V_x, V_y, Vz = np.arcsin(1*Lamda/0.75e-3/2/3), -1*np.arcsin(1*Lamda/0.75e-3/2/3/10),  1 # Reference
V_x_2, V_y_2, Vz_2 = 0, 0.00, 1 # Measurement
parameter_note = 'Ref. (%f,'%V_x + ' %f,'%V_y + '%f)\n'%Vz + 'Mea. (%f,'%V_x_2 + ' %f,'%V_y_2 + '%f)\n'%Vz_2

px = [0, 0.25e-3, 0.5e-3, 0.75e-3]
py = [0, -0.5e-3]

L = 0
D = 0

for Y in py:
    for X in px:
        D_d = D - X * np.tan(V_x_2 - V_x) + Y * np.tan(V_y_2 - V_y)
        L_d = L + X * V_x + Y * V_y
        diff_L_d_2 = D_d + (L_d + D_d) * (1 - (V_x_2**2 + V_y_2**2)) / (1 + (V_x_2**2 + V_y_2**2)) - L_d * (1 - (V_x**2 + V_y**2)) / (1 + (V_x**2 + V_y**2))
        Phase_change = 2*np.pi/Lamda*diff_L_d_2/np.pi*180
        print(Phase_change)

X, Y = np.meshgrid(dx, dy)
D_d = D - X * np.tan(V_x_2 - V_x) + Y * np.tan(V_y_2 - V_y)
L_d = L + X * V_x + Y * V_y
diff_L_d_2 = D_d + (L_d + D_d) * (1 - (V_x_2**2 + V_y_2**2)) / (1 + (V_x_2**2 + V_y_2**2)) - L_d * (1 - (V_x**2 + V_y**2)) / (1 + (V_x**2 + V_y**2))
# Phase_change = 2*np.pi/Lamda*diff_L
Z = 2*np.pi/Lamda*diff_L_d_2  /np.pi*180

if 1:
    plt.figure('Heterodyne_3D')
    c = plt.gca().pcolor(X*1000, Y*1000, Z, cmap='RdBu')
    plt.plot([i*1000 for i in px], [i*1000 for i in [0,0,0,0]], color='black', linestyle='solid', linewidth=0, marker='o', fillstyle='none', markeredgecolor='cyan', markeredgewidth=2)
    plt.plot([i*1000 for i in px], [i*1000 for i in [-0.5e-3,-0.5e-3,-0.5e-3,-0.5e-3]], color='black', linestyle='solid', linewidth=0, marker='o', fillstyle='none', markeredgecolor='cyan', markeredgewidth=2)
    plt.text(0, 0, parameter_note, fontsize=10, color = 'cyan', style='oblique')
    plt.xlabel('dx [mm]')
    plt.ylabel('dy [mm]')
    plt.colorbar(c, ax=plt.gca())

    plt.show()