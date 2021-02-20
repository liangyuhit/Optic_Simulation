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

angle_beta_2 = angle_beta * (0.5)
angle_alpha_2 = angle_beta_2 / 2
print('Angle_Alpha_2 = %e [rad]'%angle_alpha_2)


L = 0.2
D = 0.1
Displacement = np.linspace(0, 5, num=100001) * Lamda
D = D + Displacement

diff_L_2 = D + (L + D) * np.cos(2*angle_alpha_2) - L * np.cos(2*angle_alpha)
I_beat_2 = 2 * I_0 * (1 + np.cos(2 * np.pi * diff_L_2 / Lamda))


if 1: # Center Point
    plt.figure(1)
    plt.title('Center Point Intensity(Two tilted mirror)')
    plt.plot(Displacement/Lamda, I_beat_2, color='blue', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
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
