# -*- coding: utf-8 -*-
 
import matplotlib.pyplot as plt
import numpy as np

j = complex(0, 1)
c = 3e8 # 光速 [m/s]
timeline = np.linspace(0,1,1001)

Lamda = 633e-9 # 光波长 [m]
Fc = c / Lamda # 光频率 [Hz]

I_0 = 1

V_x, V_y, Vz = 0.0002*633/600, 0.000, 1
V_x_2, V_y_2, Vz_2 =0.00, 0.00,  1
print('V_x = %f'%V_x)

L = 0.2 
D = 0 + timeline * 0.5 * Lamda
R = 0.5
T = 1-R


fig_color = ['black', 'red', 'blue', 'green', 'magenta', 'yellow', 'cyan', 'grey', 'orange', 'brown']


dx = [-0.25e-3, 0, 0.25e-3, 0.5e-3]
dy = [0]
# dy = [0, -0.5e-3]
I_origin = []

### Origin Values
for i in range(len(dy)):
    for k in range(len(dx)):
        D_d = D - dx[k] * np.tan(V_x_2 - V_x) + dy[i] * np.tan(V_y_2 - V_y)
        L_d = L + dx[k] * V_x + dy[i] * V_y
        diff_L_d_2 = D_d + (L_d + D_d) * (1 - (V_x_2**2 + V_y_2**2)) / (1 + (V_x_2**2 + V_y_2**2)) - L_d * (1 - (V_x**2 + V_y**2)) / (1 + (V_x**2 + V_y**2))
        a = np.array(2 * I_0 *2*(R*T)**0.5 * (1 + np.cos(2 * np.pi * diff_L_d_2 / Lamda)))
        I_origin.append(a)



##### Parameter changed

I_changed = []
# V_x_2 = 0.0005 
# I_0 = I_0 *1.01
# R = R * 1.1
# T = 1-R

for i in range(len(dy)):
    for k in range(len(dx)):
        D_d = D - dx[k] * np.tan(V_x_2 - V_x) + dy[i] * np.tan(V_y_2 - V_y)
        L_d = L + dx[k] * V_x + dy[i] * V_y
        diff_L_d_2 = D_d + (L_d + D_d) * (1 - (V_x_2**2 + V_y_2**2)) / (1 + (V_x_2**2 + V_y_2**2)) - L_d * (1 - (V_x**2 + V_y**2)) / (1 + (V_x**2 + V_y**2))
        a = np.array(2 * I_0 * 2*(R*T)**0.5 *(1 + np.cos(2 * np.pi * diff_L_d_2 / Lamda)))
        I_changed.append(a)

if 1:

    plt.figure(1)
    for k in range(len(dx)):
        plt.subplot(2,4,k+1)
        plt.plot(D/Lamda, I_origin[k], color=fig_color[k], linestyle='dashed', linewidth=1, marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
        plt.plot(D/Lamda, I_changed[k], color=fig_color[k], linestyle='solid', linewidth=1, marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
        plt.xlabel('D/Lamda')
#     plt.title('Phase Delay = %f [deg]'%(k * np.pi/4 /(2*np.pi) * 360))

    for k in range(len(dx)-1):
        plt.subplot(2,4,k+6)    
        plt.plot(I_origin[0], I_origin[k+1], color=fig_color[k+1], linestyle='dashed', linewidth=1, marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
        plt.plot(I_changed[0], I_changed[k+1], color=fig_color[k+1], linestyle='solid', linewidth=1, marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.5, hspace=0.5)
    plt.get_current_fig_manager().window.setGeometry(20, 50, 1600, 700)

            
    plt.show()
    
 