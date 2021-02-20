# -*- coding: utf-8 -*-
 
import matplotlib.pyplot as plt
import numpy as np

j = complex(0, 1)
c = 3e8 # 光速 [m/s]
Lamda = 633e-9 # 光波长 [m]
Fc = c / Lamda # 光频率 [Hz]

# Fs = 8 * Fc # 数据点频率 [Hz]
# T0 = 1 / Fs  # 数据点间隔 [s]
# N = 2**12 # 数据点数
# T = N * T0 # 数据时长 [s]
# timeline = np.arange(N) * T0 # 时间序列
 
screen_diameter = 2e-3
I_0 = 1
angle_beta = np.arcsin(10*Lamda/screen_diameter)
angle_beta = np.pi/2
angle_alpha = angle_beta / 2
d = np.linspace(screen_diameter/2*(-1), screen_diameter/2, num=100001)

L_0 = 0.2
D_0 = 0.5
# D_0 = 0.5 * 2 * Lamda * d / (screen_diameter/2)
print(D_0)
# d_repeat = Lamda * (1 + np.tan(angle_alpha)*np.tan(angle_beta)) / np.tan(angle_alpha) / (1 + 1/np.cos(angle_beta))
# print('d_repeat = %e [m]'%d_repeat)

D_d = D_0 - d * np.tan(angle_alpha)
L_d = L_0 + d * np.tan(angle_alpha)

diff_L_d = 2 * D_d + L_d * (1 - 1/np.cos(angle_beta) + 2*np.tan(angle_alpha)*np.tan(angle_beta)) / (1 + np.tan(angle_alpha)*np.tan(angle_beta))
I_beat_d = 2 * I_0 * (1 + np.cos(2 * np.pi * diff_L_d / Lamda))
print('Angle_Beta = %e'%angle_beta)

#### For Two Tilted Mirrors:
angle_beta_2 = angle_beta * (1)
angle_alpha_2 = angle_beta_2 / 2
L_0_2 = L_0 + D_0
D_d = D_0 - d * np.tan(angle_alpha) + d * np.tan(angle_alpha_2)
L_d = L_0 + d * np.tan(angle_alpha)
L_d_2 = L_0_2 + d * np.tan(angle_alpha_2)


 
 
diff_L_d_2 = D_d + L_d_2 * (1/np.cos(angle_beta_2) - np.tan(angle_alpha_2)*np.tan(angle_beta_2)) / (1 + np.tan(angle_alpha_2)*np.tan(angle_beta_2)) - L_d * (1/np.cos(angle_beta) - np.tan(angle_alpha)*np.tan(angle_beta)) / (1 + np.tan(angle_alpha)*np.tan(angle_beta))
I_beat_d_2 = 2 * I_0 * (1 + np.cos(2 * np.pi * diff_L_d_2 / Lamda))
# print(I_beat_d_2)
print('Angle_Beta_2 = %e'%angle_beta_2)


if 1:
    plt.figure(1)
#     plt.subplot(211)
    plt.plot(d*1000, I_beat_d, color='blue', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
#     plt.plot(d*1000, I_beat_d_2, color='red', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)

    plt.xlabel('Screen Position [mm]')
    plt.ylabel('Intensity [A.U.]')
    plt.ylim(-0.2,4.2)
    plt.title('One tilted mirror')
    plt.grid(which = 'both')
    
#     plt.subplot(212)
#     plt.plot(d*1000, I_beat_d_2, color='blue', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
#     plt.xlabel('Screen Position [mm]')
#     plt.ylabel('Intensity [A.U.]')
#     plt.ylim(-0.5,4.5)
#     plt.title('Two tilted mirrors')
#     plt.grid(which = 'both')
    plt.show()
