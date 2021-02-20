# -*- coding: utf-8 -*-
 
import matplotlib.pyplot as plt
import numpy as np

j = complex(0, 1)
c = 3e8 # 光速 [m/s]
Lamda = 1550e-9 # 光波长 [m]
Fc = c / Lamda # 光频率 [Hz]
Dif_F = 1e11 # 外差频率

  
Fs = 100 * Fc # 数据点频率 [Hz]
T0 = 1 / Fs  # 数据点间隔 [s]
N = 1e6 # 数据点数
T = N * T0 # 数据时长 [s]
timeline = np.arange(N) * T0 # 时间序列
print(T)

A1 = 1 #幅值
A2 = 1
Phase_REF = 0 #参考信号相位
Phase_r = 0# 参考光
# Phase3 = np.pi * 1
L = Lamda * 0.5
Phase_m = L / Lamda * np.pi * 2 * np.arange(N)/N #测量光相位

E_r1 = A1/np.sqrt(2) * np.exp(j * (2 * np.pi * Fc * timeline + Phase_REF)) #参考信号1
E_r2 = A2/np.sqrt(2) * np.exp(j * (2 * np.pi * (Fc + Dif_F) * timeline + Phase_REF)) #参考信号2
E_m1 = A2/np.sqrt(2) * np.exp(j * (2 * np.pi * Fc * timeline + Phase_r)) #参考光
E_m2 = A2/np.sqrt(2) * np.exp(j * (2 * np.pi * (Fc + Dif_F) * timeline + Phase_m)) #测量光


Eref = E_r1 + E_r2
Iref = Eref * np.conj(Eref)
Iref_calc = (A1/np.sqrt(2))**2 + (A2/np.sqrt(2))**2 + 2*(A1/np.sqrt(2))*(A2/np.sqrt(2))*np.cos(2*np.pi*Dif_F*timeline)

Emeas = E_m1 + E_m2
Imeas = Emeas * np.conj(Emeas)
Imeas_calc = (A1/np.sqrt(2))**2 + (A2/np.sqrt(2))**2 + 2*(A1/np.sqrt(2))*(A2/np.sqrt(2))*np.cos(2*np.pi*Dif_F*timeline + Phase_m - Phase_r)


if 1:
    plt.figure(1)
     
    plt.subplot(211)
    plt.plot(timeline, Eref.real, color='blue', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    plt.plot(timeline, Emeas.real, color='red', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    
    plt.xlabel('Time [s]')
    plt.ylabel('E [V]')
#     plt.ylim(-2,2)
    plt.grid(which = 'both')
     
    plt.subplot(212)
    plt.plot(timeline, Iref.real, color='blue', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    plt.plot(timeline, Imeas.real, color='red', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    plt.plot(timeline, Iref_calc, color='blue', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    plt.plot(timeline, Imeas_calc, color='red', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)

    plt.xlabel('Time [s]')
    plt.ylabel('I [V]')
#     plt.ylim(-2,2)
    plt.grid(which = 'both')
     
     
    plt.show()
 






