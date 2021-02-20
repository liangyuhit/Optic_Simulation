# -*- coding: utf-8 -*-
 
import matplotlib.pyplot as plt
import numpy as np

j = complex(0, 1)
c = 3e8 # 光速 [m/s]
Lamda = 1550e-9 # 光波长 [m]
Fc = c / Lamda # 光频率 [Hz]
print(Fc)

T = 100 / Fc # 数据时长 [s]
N = 1e4
T0 = T / N
timeline = np.arange(N) * T0 # 时间序列

A = 1 #幅值
Phase1 = 0
L1 = Lamda
L2 = L1 + 4 * Lamda / N * np.arange(N)

E0 = A * np.exp(j * (2 * np.pi * Fc * timeline + Phase1)) #电场式
# I0 = E_ref * np.conj(E_ref)

E1 = E0/(2**0.5)
E2 = E0/(2**0.5) * np.exp(j * np.pi/2)
E3 = E1 * np.exp(j * 2 * np.pi * L1 / Lamda)
E4 = E2 * np.exp(j * 2 * np.pi * L2 / Lamda)
E5 = E3/(2**0.5) * np.exp(j * np.pi/2)
E6 = E4/(2**0.5) * np.exp(j * np.pi/2)
E7 = E3/(2**0.5)
E8 = E4/(2**0.5) 
 
 
E_A = E5 + E8
E_B = E6 + E7
 
I_A = E_A * np.conj(E_A)
I_B = E_B * np.conj(E_B)

E_beat = E3 + E4
I_beat = E_beat * np.conj(E_beat)

print(T)


if 1:
    plt.figure(1)
     
    plt.subplot(311)
    plt.plot(timeline, E_A.real, color='blue', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    plt.plot(timeline, E_B.real, color='red', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)

    plt.xlabel('Time [s]')
    plt.ylabel('E [V]')
    plt.grid(which = 'both')
     
    plt.subplot(312)
    plt.plot(timeline, I_A.real, color='blue', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    plt.plot(timeline, I_B.real, color='red', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)


    plt.xlabel('Time [s]')
    plt.ylabel('I [V]')
    plt.grid(which = 'both')
    
    plt.subplot(313)
#     plt.plot(timeline, I_A.real+I_B.real, color='blue', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
#     plt.plot(timeline, I_A.real-I_B.real, color='red', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    plt.plot(timeline, I_beat.real, color='red', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    plt.xlabel('Time [s]')
    plt.ylabel('I [V]')
    plt.grid(which = 'both')
    
     
    plt.show()
