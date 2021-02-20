# -*- coding: utf-8 -*-
 
import matplotlib.pyplot as plt
import numpy as np

j = complex(0, 1)
c = 3e8 # 光速 [m/s]
Lamda = 1550e-9 # 光波长 [m]
Fc = c / Lamda # 光频率 [Hz]
print(Fc)
Fm = 5e13 # 调制频率

Fs = 8 * Fc # 数据点频率 [Hz]
T0 = 1 / Fs  # 数据点间隔 [s]
N = 2**12 # 数据点数
T = N * T0 # 数据时长 [s]
timeline = np.arange(N) * T0 # 时间序列

Mod_Depth = 0.8
A1 = 1 #幅值
Phase1 = 0
L = Lamda * 0.06
Phase_pm = L / Lamda * np.pi * 2 * np.arange(N)/N
Phase_pm = Mod_Depth * np.sin(2 * np.pi * Fm * timeline) + Phase_pm
# Phase_pm = Mod_Depth * np.sin(2 * np.pi * Fm * timeline)

E_ref = A1 * np.exp(j * (2 * np.pi * Fc * timeline + Phase1)) #电场式
I_ref = E_ref * np.conj(E_ref)

E_pm = A1 * np.exp(j * (2 * np.pi * Fc * timeline + Phase_pm)) #电场式
I_pm = E_pm * np.conj(E_pm)

# E_2nd = A1 * np.exp(j * (2 * np.pi * Fc * timeline + Phase2)) #电场式

# I_beat = np.cos(np.sin(2 * np.pi * Fm * timeline))
E_beat = E_pm + E_ref
I_beat = E_beat * np.conj(E_beat)

# E_beat_delay = E_beat * np.exp(j * np.pi)
# I_beat_delay = E_beat_delay * np.conj(E_beat_delay)

print(T)

if 1:
    freqline = np.fft.fftfreq(N, d=T0)[1:N//2]
#     Data_FFT= np.array(abs(np.fft.fft(E_pm))[1:N//2])*2/N
    Data_FFT= np.array(abs(np.fft.fft(I_beat))[1:N//2])*2/N

if 1:
    plt.figure(1)
     
    plt.subplot(311)
    plt.plot(timeline, E_beat.real, color='blue', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
#     plt.plot(timeline, E_beat_delay.real, color='red', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    plt.xlabel('Time [s]')
    plt.ylabel('E [V]')
#     plt.ylim(-2,2)
    plt.grid(which = 'both')
     
    plt.subplot(312)
    plt.plot(timeline, I_beat.real, color='blue', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
#     plt.plot(timeline, I_beat_delay.real, color='red', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)

    plt.xlabel('Time [s]')
    plt.ylabel('I [V]')
#     plt.ylim(-2,2)
    plt.grid(which = 'both')
    
    plt.subplot(313)
#     plt.plot(freqline, Data_FFT, color='blue', marker='o', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    plt.semilogy(freqline, Data_FFT, color='blue', marker='o', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)

    plt.xlabel('Freq [Hz]')
    plt.ylabel('FFT [V]')
#     plt.ylim(-2,2)
    plt.grid(which = 'both')
     
    plt.show()

