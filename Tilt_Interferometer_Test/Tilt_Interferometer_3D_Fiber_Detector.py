# -*- coding: utf-8 -*-
 
import matplotlib.pyplot as plt
import numpy as np

j = complex(0, 1)
c = 3e8 # 光速 [m/s]

Lamda = 633e-9 # 光波长 [m]
Fc = c / Lamda # 光频率 [Hz]

I_0 = 1

V_x, V_y, Vz = 0.0002*633/600, 0.000, 1
V_x_2, V_y_2, Vz_2 = 0.00, 0.00, 1

timeline = np.linspace(0,1,1001)
L = 0.2
D = timeline * 2 * Lamda


dx = [-0.25e-3, 0, 0.25e-3, 0.5e-3]
dy = [0, -0.5e-3]
I = []

for i in range(len(dy)):
    for k in range(len(dx)):
        D_d = D - dx[k] * np.tan(V_x_2 - V_x) + dy[i] * np.tan(V_y_2 - V_y)
        L_d = L + dx[k] * V_x + dy[i] * V_y
        diff_L_d_2 = D_d + (L_d + D_d) * (1 - (V_x_2**2 + V_y_2**2)) / (1 + (V_x_2**2 + V_y_2**2)) - L_d * (1 - (V_x**2 + V_y**2)) / (1 + (V_x**2 + V_y**2))
        a = np.array(2 * I_0 * (1 + np.cos(2 * np.pi * diff_L_d_2 / Lamda)))
        I.append(a)


if 0:
    plt.figure(1)
    
    for m in range(len(dx)*len(dy)):
        plt.subplot(2,4,m+1)
        plt.plot(D/Lamda, I[m])
        plt.xlabel('Movement/Lamda')
        plt.ylabel('I @d%i'%(m+1))
            
    plt.show()
    
    
if 1:
    plt.figure(2)
    plt.plot(I[0], I[1], color='blue', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    plt.plot(I[0], I[2], color='red', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    plt.plot(I[0], I[3], color='black', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)

    plt.show()
