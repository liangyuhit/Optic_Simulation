# -*- coding: utf-8 -*-
 
import matplotlib.pyplot as plt
import numpy as np

j = complex(0, 1)
c = 3e8 # 光速 [m/s]
Lamda = 1550e-9 # 光波长 [m]
Fc = c / Lamda # 光频率 [Hz]

  
Fs = 100 * Fc # 数据点频率 [Hz]
T0 = 1 / Fs  # 数据点间隔 [s]
N = 1e6 # 数据点数
T = N * T0 # 数据时长 [s]
timeline = np.arange(N) * T0 # 时间序列
  
A1 = 1#幅值
A2 = A1
Phase1 = 0 #相位
# Phase2 =  np.pi * 1.8
L = Lamda * 4
Phase2 = L / Lamda * np.pi * 2 * np.arange(N)/N
E1 = A1 * np.exp(j * (2 * np.pi * Fc * timeline + Phase1)) #电场式
E2 = A2 * np.exp(j * (2 * np.pi * Fc * timeline + Phase2)) #电场式
E = E1 + E2
I = E * np.conj(E)
print(T)

if 1:
    plt.figure(1)
     
    plt.subplot(211)
    plt.plot(timeline, E.real, color='blue', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    plt.xlabel('Time [s]')
    plt.ylabel('E [V]')
#     plt.ylim(-2,2)
    plt.grid(which = 'both')
     
    plt.subplot(212)
    plt.plot(timeline, I.real, color='blue', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    plt.xlabel('Time [s]')
    plt.ylabel('I [V]')
#     plt.ylim(-2,2)
    plt.grid(which = 'both')
     
     
    plt.show()
 





