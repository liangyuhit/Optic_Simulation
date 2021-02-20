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
dx = np.linspace(screen_diameter/2*(-1), screen_diameter/2, num=1001)
dy = dx

V_x, V_y, Vz = 0.002, 0.002, 1 # Reference
V_x_2, V_y_2, Vz_2 = 0.000, 0.000, 1 # Measurement

tana = (V_x**2 + V_y**2 + 1)**0.5
tana_2 = (V_x_2**2 + V_y_2**2 + 1)**0.5

L = 0.2
D = 0.1
 
fringe_slope = (V_x_2 * tana**2 - V_x * tana_2**2) / (V_y_2 * tana**2 - V_y * tana_2**2)
print(fringe_slope)

X, Y = np.meshgrid(dx, dy)
D_d = D - X * np.tan(V_x_2 - V_x) + Y * np.tan(V_y_2 - V_y)
L_d = L + X * V_x + Y * V_y
diff_L_d_2 = D_d + (L_d + D_d) * (1 - (V_x_2**2 + V_y_2**2)) / (1 + (V_x_2**2 + V_y_2**2)) - L_d * (1 - (V_x**2 + V_y**2)) / (1 + (V_x**2 + V_y**2))
Z = 2 * I_0 * (1 + np.cos(2 * np.pi * diff_L_d_2 / Lamda))
center_x = int(len(dy)/2)+1
center_y = int(len(dx)/2)+1
# print(center_x, center_y)
# center_horizontal = Z[center_x]
center_horizontal = []
center_vertical = []
for i in range(len(Z[center_x])):
    center_horizontal.append(Z[center_x][i])
for i in range(len(Z)):
    center_vertical.append(Z[i][center_y])
# for i in range(2**15-len(Z)):
#     center_horizontal.append(0)
#     center_vertical.append(0)


print(len(center_horizontal))
FFT_horizontal = FFT_cal(center_horizontal, screen_diameter/len(dx))
FFT_vertical = FFT_cal(center_vertical, screen_diameter/len(dy))

center_horizontal_freq = FFT_horizontal[0][FFT_horizontal[1][100:].argmax()+100]
center_vertical_freq = FFT_vertical[0][FFT_vertical[1][100:].argmax()+100]
print(center_horizontal_freq, center_vertical_freq, center_horizontal_freq/center_vertical_freq)
print(FFT_horizontal[2])

if 1:
#     plt.figure(1)
#     c = plt.gca().pcolor(X*1000, Y*1000, Z, cmap='gray')
#     plt.xlabel('dx [mm]')
#     plt.ylabel('dy [mm]')
#     plt.colorbar(c, ax=plt.gca())

    plt.figure(2)
    plt.subplot(3,1,1)
    plt.plot(center_horizontal)
    plt.subplot(3,1,2)
    plt.plot(FFT_horizontal[0], FFT_horizontal[1])
    plt.subplot(3,1,3)
    plt.plot(FFT_horizontal[0], FFT_horizontal[2])
 
    plt.figure(3)
    plt.subplot(3,1,1)
    plt.plot(center_vertical)
    plt.subplot(3,1,2)
    plt.plot(FFT_vertical[0], FFT_vertical[1])
    plt.subplot(3,1,3)
    plt.plot(FFT_vertical[0], FFT_vertical[2])
    
    plt.show()