# -*- coding: utf-8 -*-
 
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

def FFT_cal(Data, tau0): ### tau0 = (screen_diameter/len(dx))
    N_Data = len(Data)
    freqline = np.fft.fftfreq(len(Data), d=tau0)
    Data_FFT = np.fft.fft(Data) #two sides FFT
    
    Data_Phase = np.angle(Data_FFT)
    Data_FFT= abs(Data_FFT)
    
    freqline_half = freqline[1:len(Data)//2] # One side
    Data_FFT = Data_FFT[1:len(Data)//2]
    Data_Phase = Data_Phase[1:len(Data)//2]
    Data_FFT = np.array(Data_FFT)   
    Data_FFT = Data_FFT*2/len(Data)
    
    return freqline_half, Data_FFT, Data_Phase


j = complex(0, 1)
c = 3e8 # 光速 [m/s]
Lamda = 633e-9 # 光波长 [m]
Fc = c / Lamda # 光频率 [Hz]

screen_diameter = 1.5e-3
I_0 = 1
dx = np.linspace(screen_diameter/2*(-1), screen_diameter/2, num=1000)
dy = dx

slope_cal_set = []
slope_sim_set = []
Vx_cal_set = []
Vy_cal_set = []
V_x_2_set = []
V_y_2_set = []
for m in range(10):
    V_x, V_y, Vz = 0.001, 0.001, 1 # Reference
    V_x_2, V_y_2, Vz_2 = 0.0000001*m, 0.0000005*m, 1 # Measurement
    V_x_2_set.append(V_x_2)
    V_y_2_set.append(V_y_2)
    
    tana = (V_x**2 + V_y**2)**0.5
    tana_2 = (V_x_2**2 + V_y_2**2)**0.5
    
    L = 0.2
    D = 0.1
     
    fringe_slope = (V_x_2 * tana**2 - V_x * tana_2**2) / (V_y_2 * tana**2 - V_y * tana_2**2)
#     print(fringe_slope)
    slope_cal_set.append(fringe_slope)
    
    X, Y = np.meshgrid(dx, dy)
    D_d = D - X * np.tan(V_x_2 - V_x) + Y * np.tan(V_y_2 - V_y)
    L_d = L + X * V_x + Y * V_y
    diff_L_d_2 = D_d + (L_d + D_d) * (1 - (V_x_2**2 + V_y_2**2)) / (1 + (V_x_2**2 + V_y_2**2)) - L_d * (1 - (V_x**2 + V_y**2)) / (1 + (V_x**2 + V_y**2))
    Z = 2 * I_0 * (1 + np.cos(2 * np.pi * diff_L_d_2 / Lamda))
    # print(center_x, center_y)
    # center_horizontal = Z[center_x]
    horizontal_set = []
    vertical_set = []
    for i in range(int(len(Z)/10)):
        horizontal = []
        for j in range(len(Z[i*10])):
            horizontal.append(Z[i*10][j])
    #     print(len(horizontal))
    #     print(horizontal)
        for i in range(99*len(horizontal)):
            horizontal.append(0)
        FFT_horizontal = FFT_cal(horizontal, screen_diameter/len(dx))
        horizontal_freq = FFT_horizontal[0][FFT_horizontal[1][100:].argmax()+100]
    #     print(horizontal_freq)
        horizontal_set.append(horizontal_freq)
    
    for j in range(int(len(Z[0])/10)):
        vertical = []
        for i in range(len(Z)):
            vertical.append(Z[i][j*10])
        for i in range(99*len(vertical)):
            vertical.append(0)
        FFT_vertical = FFT_cal(vertical, screen_diameter/len(dy))
        vertical_freq = FFT_vertical[0][FFT_vertical[1][100:].argmax()+100]
        vertical_set.append(vertical_freq)
    
#     print(FFT_horizontal[0])
    print('horizontal freq:%f'%np.average(horizontal_set))
    print('vertical freq:%f'%np.average(vertical_set))
    period_x = 1/np.average(horizontal_set)
    period_y = 1/np.average(vertical_set)
    print('horizontal period:%f'%period_x)
    print('vertical period:%f'%period_y)
    V_x_cal = -Lamda/2/period_x + V_x
    V_y_cal = -Lamda/2/period_y + V_y
    Vx_cal_set.append(V_x_cal)
    Vy_cal_set.append(V_y_cal)
    print('horizontal angle:%f'%V_x_cal, 'vertical angle:%f'%V_y_cal)
    slope_sim = np.average(horizontal_set)/np.average(vertical_set)
    slope_sim_set.append(slope_sim)
    print(fringe_slope, slope_sim)
    
if 1:
    plt.figure('x')
    plt.plot(np.array(V_x_2_set)*1e6, 'k', marker='o',fillstyle='full', markeredgecolor='k')
    plt.plot(np.array(Vx_cal_set)*1e6, 'b', marker='o',fillstyle='full', markeredgecolor='b')
    
    plt.figure('y')
    plt.plot(np.array(V_y_2_set)*1e6, 'k', marker='o',fillstyle='full', markeredgecolor='k')
    plt.plot(np.array(Vy_cal_set)*1e6, 'b', marker='o',fillstyle='full', markeredgecolor='b')

    plt.show()
    
if 0:
    plt.figure(1)
#     plt.plot(horizontal)
#     plt.plot(FFT_horizontal[0], FFT_horizontal[1], color='blue', marker='o', fillstyle='full', markeredgecolor='blue', markeredgewidth=0.0)
    plt.plot(slope_cal_set, color='blue', marker='o', fillstyle='full', markeredgecolor='blue', markeredgewidth=0.0)
    plt.plot(slope_sim_set, color='red', marker='o', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    
    plt.figure(2)
    plt.plot(np.array(slope_cal_set) - np.array(slope_sim_set))
    
    plt.show()

if 0:
    plt.figure(1)
    c = plt.gca().pcolor(X*1000, Y*1000, Z, cmap='gray')
    plt.xlabel('dx [mm]')
    plt.ylabel('dy [mm]')
    plt.colorbar(c, ax=plt.gca())

#     plt.figure(2)
#     plt.subplot(2,1,1)
# #     plt.plot(center_horizontal)
# #     plt.subplot(3,1,2)
#     plt.plot(FFT_horizontal[0], FFT_horizontal[1])
# #     plt.subplot(3,1,3)
# #     plt.plot(FFT_horizontal[0], FFT_horizontal[2])
#   
# #     plt.figure(3)
#     plt.subplot(2,1,2)
# #     plt.plot(center_vertical)
# #     plt.subplot(3,1,2)
#     plt.plot(FFT_vertical[0], FFT_vertical[1])
# #     plt.subplot(3,1,3)
# #     plt.plot(FFT_vertical[0], FFT_vertical[2])
#      
    plt.show()
