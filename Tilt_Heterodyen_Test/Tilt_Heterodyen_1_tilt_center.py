# -*- coding: utf-8 -*-
 
import matplotlib.pyplot as plt
import numpy as np

j = complex(0, 1)
c = 3e8 # 光速 [m/s]
Lamda = 633e-9 # 光波长 [m]
Fc = c / Lamda # 光频率 [Hz]
Dif_F = 1e6 # 外差频率
Fs = 100 * Dif_F # 数据点频率 [Hz]
T0 = 1 / Fs  # 数据点间隔 [s]
N = 0.5e3 # 数据点数
T = N * T0 # 数据时长 [s]
timeline = np.linspace(0, N, N+1) * T0 # 时间序列
print(T)


screen_diameter = 2e-3
I_0 = 1
A1 = 1 #幅值
A2 = 1
Phase1 = 0 #相位
Phase2 = 0

angle_beta = np.arcsin(10*Lamda/screen_diameter)
angle_alpha = angle_beta / 2
print('Angle_Alpha = %e [rad]'%angle_alpha)

L = 0.2
D = 0
Displacement = np.linspace(0, 0.5, num=N+1) * Lamda
D = D + Displacement

diff_L = 2 * D + 2 * L * np.tan(angle_alpha)**2 / (1 + np.tan(angle_alpha)**2)
Phase_change = 2*np.pi/Lamda*diff_L
I_beat_mea = (A1/np.sqrt(2))**2 + (A2/np.sqrt(2))**2 + 2*(A1/np.sqrt(2))*(A2/np.sqrt(2))*np.cos(2*np.pi*Dif_F*timeline + Phase_change)
I_beat_ref = (A1/np.sqrt(2))**2 + (A2/np.sqrt(2))**2 + 2*(A1/np.sqrt(2))*(A2/np.sqrt(2))*np.cos(2*np.pi*Dif_F*timeline)

if 1: # Center Point
    plt.figure(1)
    plt.title('Center Point Intensity(One tilted mirror)')
#     plt.plot(timeline, I_beat_ref, color='blue', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
#     plt.plot(timeline, I_beat_mea, color='red', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    plt.plot(Displacement/Lamda, I_beat_ref, color='blue', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    plt.plot(Displacement/Lamda, I_beat_mea, color='red', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
#     plt.plot(Displacement/Lamda, I_beat_final, color='red', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    plt.xlabel('Displacement/Lamda')
    plt.ylabel('Intensity [A.U.]').set_color('blue')
    [i.set_color("blue") for i in plt.gca().get_yticklabels()]
#     plt.ylim(-0.2,4.2)
    plt.grid(which = 'both')
    
    plt.twinx().plot(Displacement/Lamda, Phase_change, color='black', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
#     plt.ylim(-0.2,4.2)
    plt.ylabel('Phase_change').set_color('black')
    [i.set_color("black") for i in plt.gca().get_yticklabels()]
    plt.grid(which = 'both')
    
    plt.show()
