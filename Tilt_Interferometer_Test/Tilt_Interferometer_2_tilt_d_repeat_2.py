
# -*- coding: utf-8 -*-
 
import matplotlib.pyplot as plt
import numpy as np


c = 3e8 # 光速 [m/s]
Lamda = 633e-9 # 光波长 [m]
Fc = c / Lamda # 光频率 [Hz]



# angle_beta = np.linspace(0.0005, 0.05, 10001)
angle_beta = 3.165e-3
angle_alpha = angle_beta / 2
print('Angle_Alpha = %e [rad]'%angle_alpha)

# angle_beta_2 = angle_beta * (1 + np.hstack((np.linspace(-2, -0.02, 1000),np.linspace(0.02, 2, 1000))))
angle_beta_2 = angle_beta * (1 + np.linspace(0.02, 4, 5000))
angle_alpha_2 = angle_beta_2 / 2
# print('Angle_Beta_2 = %e [rad]'%angle_beta_2)

d_repeat_2 = Lamda / (np.sin(2*angle_alpha_2) - np.sin(2*angle_alpha))


if 1: # Center Point
    plt.figure(1)
    plt.title('fringe distance')
    plt.plot((angle_beta_2)/angle_beta, d_repeat_2*1000, color='blue', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    
    angle_beta_2 = angle_beta * (1 + np.linspace(-4, -0.02, 5000))
    angle_alpha_2 = angle_beta_2 / 2
    # print('Angle_Beta_2 = %e [rad]'%angle_beta_2)
    
    d_repeat_2 = Lamda / (np.sin(2*angle_alpha_2) - np.sin(2*angle_alpha))

    plt.plot((angle_beta_2)/angle_beta, d_repeat_2*1000, color='blue', marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)

    
    
    plt.xlabel('beta_2/beta')
    plt.ylabel('d_repeat_2 [mm]')
#     [i.set_color("blue") for i in plt.gca().get_yticklabels()]
    plt.ylim(-10, 10)
    plt.grid(which = 'both')
    
    plt.show()
