# -*- coding: utf-8 -*-
 
import matplotlib.pyplot as plt
import numpy as np

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
Displacement = np.linspace(0, 5, num=100001) * Lamda
D = D + Displacement
# diff_L = 2 * D + L * (1 - 1/np.cos(angle_beta) + 2*np.tan(angle_alpha)*np.tan(angle_beta)) / (1 + np.tan(angle_alpha)*np.tan(angle_beta))
# I_beat = 2 * I_0 * (1 + np.cos(2 * np.pi * diff_L / Lamda))

diff_L = 2 * D + 2 * L * np.tan(angle_alpha)**2 / (1 + np.tan(angle_alpha)**2)
I_beat = 2 * I_0 * (1 + np.cos(2 * np.pi * diff_L / Lamda))


if 1: # Center Point
    plt.figure(1)
    plt.title('Center Point Intensity(One tilted mirror)')
    plt.plot(Displacement/Lamda, I_beat, color='blue', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
#     plt.plot(Displacement/Lamda, I_beat_final, color='red', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    plt.xlabel('Displacement/Lamda')
    plt.ylabel('Intensity [A.U.]').set_color('blue')
    [i.set_color("blue") for i in plt.gca().get_yticklabels()]
    plt.ylim(-0.2,4.2)
    plt.grid(which = 'both')
    
    plt.twinx().plot(Displacement/Lamda, Displacement/Lamda, color='red', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
#     plt.ylim(-0.2,4.2)
    plt.ylabel('Displacement/Lamda').set_color('red')
    [i.set_color("red") for i in plt.gca().get_yticklabels()]
    plt.grid(which = 'both')
    
    plt.show()
