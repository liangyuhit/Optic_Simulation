# -*- coding: utf-8 -*-
 
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


j = complex(0, 1)
c = 3e8 # 光速 [m/s]
Lamda = 633e-9 # 光波长 [m]
Fc = c / Lamda # 光频率 [Hz]

screen_diameter = 2e-3
I_0 = 1
angle_beta = np.arcsin(10*Lamda/screen_diameter)
angle_alpha = angle_beta / 2
print('Angle_Alpha = %e [rad]'%angle_alpha)

L = 0.2
D = 0.1

# d_repeat = Lamda/2 * (np.tan(angle_alpha) + 1/np.tan(angle_alpha))
d = np.linspace(screen_diameter/2*(-1), screen_diameter/2, num=301)

D_d = D - d * np.tan(angle_alpha)
L_d = L + d * np.tan(angle_alpha)

diff_L_d = 2 * D_d + 2 * L_d * np.tan(angle_alpha)**2 / (1 + np.tan(angle_alpha)**2)
I_beat_d = 2 * I_0 * (1 + np.cos(2 * np.pi * diff_L_d / Lamda))

if 1: # Center Point
    plt.figure(1)
    plt.title('Intensity Distribution(One tilted mirror)')
    plt.plot(d*1000, I_beat_d, color='blue', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    plt.xlabel('Screen [mm]')
    plt.ylabel('Intensity [A.U.]')
#     [i.set_color("blue") for i in plt.gca().get_yticklabels()]
    plt.ylim(-0.2,4.2)
    plt.grid(which = 'both')    
    plt.show()
