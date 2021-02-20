# -*- coding: utf-8 -*-
 
import matplotlib.pyplot as plt
import numpy as np
from Ellipse_Fit_Func import Ellipse_fit

export_path = 'C:/Users/yu03/Desktop/Gabor/Interferometer Sim/Simulation_Result/Test_Result/New_Configuration/'

V_x_2_set = []
V_y_2_set = []
    
for m in range(8):
    j = complex(0, 1)
    c = 3e8 # ���� [m/s]
    timeline = np.linspace(0,1,501)
    
    Lamda = 633e-9 # �Ⲩ�� [m]
    Fc = c / Lamda # ��Ƶ�� [Hz]
    
    I_0 = 1
    A1 = 1 #幅值
    A2 = 1
    Phase1 = 0 #相位
    Phase2 = 0
    
    V_x, V_y, Vz = np.arcsin(1*Lamda/0.75e-3/2/3),-1*np.arcsin(1*Lamda/0.75e-3/2/3/10),  1 # Reference
    V_x_2, V_y_2, Vz_2 = 0.00, 0.00,  1
    
    L = 0
    D = 0 + timeline * 0.5 * Lamda
    R = 0.5
    T = 1-R
    
    fig_color = ['black', 'red', 'blue', 'green', 'magenta', 'cyan', 'orange', 'brown', 'grey', 'yellow']
    
    dx = [0, 0.25e-3, 0.5e-3, 0.75e-3]
    dy = [0, -0.5e-3]

    if 1:    ###modulation
        V_x_2 = V_x/20*m
        V_y_2 = V_x/20*m
        
        V_x_2_set.append(V_x_2)
        V_y_2_set.append(V_y_2)
        print('V_x_2 = %f; '%V_x_2 + 'V_y_2 = %f'%V_y_2)        
        parameter_note = 'Ref. (%f,'%V_x + ' %f,'%V_y + '%f)\n'%Vz + 'Mea. (%f,'%V_x_2 + ' %f,'%V_y_2 + '%f)\n'%Vz_2
    
    ### Origin Values
    Phase_set = []
    
    for i in range(len(dy)):
        for j in range(len(dx)):
            D_d = D - dx[j] * np.tan(V_x_2 - V_x) + dy[i] * np.tan(V_y_2 - V_y)
            L_d = L + dx[j] * V_x + dy[i] * V_y
            diff_L_d_2 = D_d + (L_d + D_d) * (1 - (V_x_2**2 + V_y_2**2)) / (1 + (V_x_2**2 + V_y_2**2)) - L_d * (1 - (V_x**2 + V_y**2)) / (1 + (V_x**2 + V_y**2))
            Phase = 2*np.pi/Lamda*diff_L_d_2  /np.pi*180
            Phase_set.append(Phase)
            
    slope_set = []
    start_phase_set = []
    for i in Phase_set:
        slope = i[-1] - i[0]
        start_phase = i[0]
        slope_set.append(slope)
        start_phase_set.append(start_phase)
    if 1:
        plt.figure('Calc.Intensity_%i'%m)
        plt.gcf().set_size_inches(16,14)
                
        for j in range(8): ### Line 2
            plt.subplot(2,4,j+1)
            plt.plot(D/Lamda, Phase_set[j], label='Phase_%i\n'%j + '%f\n'%start_phase_set[j] + '%f'%slope_set[j], color=fig_color[j], linestyle='solid', linewidth=1, marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
            legend = plt.legend(loc='upper left')
            for text in legend.get_texts():
                plt.setp(text, color = fig_color[j])
#         plt.subplots_adjust(left=None, bottom=None, right=None, top=0.95, wspace=None, hspace=0.2)
#         plt.get_current_fig_manager().window.setGeometry(20, 50, 1200, 1000)
#         plt.savefig(export_path + 'Calc.Intensity_%i.jpg'%m, dpi=300)
#         plt.close()
        plt.show()