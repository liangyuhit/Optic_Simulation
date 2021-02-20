# -*- coding: utf-8 -*-
 
import matplotlib.pyplot as plt
import numpy as np
from Ellipse_Fit_Func import Ellipse_fit

export_path = 'C:/Users/yu03/Desktop/Gabor/Interferometer Sim/Simulation_Result/Test_Result/0_degree/'

V_x_2_set = []
V_y_2_set = []

a_set = []
b_set = []
x0_set = []
y0_set = []
Delta_set = []

for i in range(10):
    a_set.append([])
    b_set.append([])
    x0_set.append([])
    y0_set.append([])
    Delta_set.append([])
    
for m in range(10):
    j = complex(0, 1)
    c = 3e8 # ���� [m/s]
    timeline = np.linspace(0,1,501)
    
    Lamda = 633e-9 # �Ⲩ�� [m]
    Fc = c / Lamda # ��Ƶ�� [Hz]
    
    I_0 = 1
    
    V_x, V_y, Vz = np.arcsin(1*Lamda/0.75e-3/2/3),-1*np.arcsin(1*Lamda/0.75e-3/2/3/10),  1 # Reference
    V_x_2, V_y_2, Vz_2 = 0.00, 0.00,  1
    
    L = 0.2 
    D = 0 + timeline * 0.5 * Lamda
    R = 0.5
    T = 1-R
    
    
    fig_color = ['black', 'red', 'blue', 'green', 'magenta', 'cyan', 'orange', 'brown', 'grey', 'yellow']
    
    
    dx = [0, 0.25e-3, 0.5e-3, 0.75e-3]
    # dy = [0]
    dy = [0, -0.5e-3]
    I_origin = []


    if 1:    ###modulation
        V_x_2 = V_x/20*m
        V_y_2 = V_x/20*m
        
        V_x_2_set.append(V_x_2)
        V_y_2_set.append(V_y_2)
        print('V_x_2 = %f; '%V_x_2 + 'V_y_2 = %f'%V_y_2)        
        parameter_note = 'Ref. (%f,'%V_x + ' %f,'%V_y + '%f)\n'%Vz + 'Mea. (%f,'%V_x_2 + ' %f,'%V_y_2 + '%f)\n'%Vz_2
    
    ### Origin Values
    for i in range(len(dy)):
        for j in range(len(dx)):
            D_d = D - dx[j] * np.tan(V_x_2 - V_x) + dy[i] * np.tan(V_y_2 - V_y)
            L_d = L + dx[j] * V_x + dy[i] * V_y
            diff_L_d_2 = D_d + (L_d + D_d) * (1 - (V_x_2**2 + V_y_2**2)) / (1 + (V_x_2**2 + V_y_2**2)) - L_d * (1 - (V_x**2 + V_y**2)) / (1 + (V_x**2 + V_y**2))
            a = np.array(2 * I_0 *2*(R*T)**0.5 * (1 + np.cos(2 * np.pi * diff_L_d_2 / Lamda)))
            I_origin.append(a)
    
    
    ##### Ecllipse Fitting
    
    fit_xaxis = fit_yaxis = np.linspace(0, 4, 501)
    fit_xaxis, fit_yaxis = np.meshgrid(fit_xaxis, fit_yaxis)
    fit_ellipse = []
    text = []
     
    for i in range(len(dx)*len(dy)-1):
        fit = Ellipse_fit(I_origin[0], I_origin[i+1])
        parameters = np.array(fit.fitEllipse())
        fit_A, fit_B, fit_C, fit_D, fit_F, fit_G = parameters
        fit_A, fit_B, fit_C, fit_D, fit_E, fit_F = parameters/fit_A
        
        fit_result = (fit_xaxis**2) + fit_B*fit_xaxis*fit_yaxis + fit_C*(fit_yaxis**2) +fit_D*fit_xaxis + fit_E*fit_yaxis + fit_F
        x0 = (2*fit_C*fit_D - fit_B*fit_E)/(fit_B**2-4*fit_C)
        y0 = (2*fit_E - fit_B*fit_D)/(fit_B**2-4*fit_C)
        a = ((x0**2 + (y0**2)*fit_C + x0*y0*fit_B - fit_F)/(1-(fit_B**2)/(4*fit_C)))**0.5
        b = ((a**2)/fit_C)**0.5
        Delta = np.arccos(-1*fit_B/2/(fit_C**0.5))
         
        fit_ellipse.append(fit_result)
        x0_set[i].append(x0)
        y0_set[i].append(y0)
        a_set[i].append(a)
        b_set[i].append(b)
        Delta_set[i].append(Delta)
         
        text.append(('Center(%f,'%x0 + '%f)\n'%y0 + 'Axes(%f,'%a + '%f)\n'%b + '%f[deg]'%(Delta/np.pi*180)))        
            
    if 1:
        
        plt.figure('Calc.Intensity_%i'%m)
        plt.gcf().set_size_inches(18,10)
        
        plt.subplot(4,4,9)
        plt.text(0, 0, parameter_note, fontsize=10, style='oblique')
        plt.xticks([])
        plt.yticks([])
        plt.axis('off')
        
        for j in range(len(dx)*len(dy)):
            plt.subplot(4,4,j+1)
            plt.plot(D/Lamda, I_origin[j], color=fig_color[j], linestyle='solid', linewidth=1, marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    #         plt.xlabel('D/Lamda')
    #         plt.axis('scaled')
            plt.title('I_%i'%j)
    #     plt.title('Phase Delay = %f [deg]'%(k * np.pi/4 /(2*np.pi) * 360))
    
        for j in range(len(dx)*len(dy)-1):
            plt.subplot(4,4,j+10)
            plt.plot(I_origin[0], I_origin[j+1], color=fig_color[j+1], linestyle='solid', linewidth=1, marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
#             plt.contour(fit_xaxis, fit_yaxis, fit_ellipse[j], [1])     #x**2 + y**2 = 9 ��Բ��
            plt.axis('scaled')
            plt.title('I_%i/I_0'%(j+1))
        
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.2)
        plt.get_current_fig_manager().window.setGeometry(20, 50, 1600, 1000)
#         plt.savefig('C:/Users/yu03/Desktop/Gabor/Interferometer Sim/Simulation_Result/1. V_x_2 Modulation/Calc.Intensity_%i.jpg'%m, dpi=300)
        plt.savefig(export_path + 'Calc.Intensity_%i.jpg'%m, dpi=300)
        plt.close()
        
        
        plt.figure('Fit_Result_%i'%m)
        plt.gcf().set_size_inches(18,10)
        
        plt.subplot(2,4,1)
        plt.text(0, 0, parameter_note, fontsize=10, style='oblique')
        plt.xticks([])
        plt.yticks([])
        plt.axis('off')

        for j in range(len(dx)*len(dy)-1):
            plt.subplot(2,4,j+2)
    #         plt.plot(I_origin[0], I_origin[j+1], color=fig_color[j+1], linestyle='solid', linewidth=1, marker=' ', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
            plt.contour(fit_xaxis, fit_yaxis, fit_ellipse[j], [0], colors=fig_color[j+1])     #x**2 + y**2 = 9 ��Բ��
            plt.axis('scaled')
            plt.text(0.5, 1, text[j], fontsize=10, style='oblique')
            plt.title('I_%i/I_0'%(j+1))

        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
        plt.get_current_fig_manager().window.setGeometry(20, 50, 1800, 800)
#         plt.savefig('C:/Users/yu03/Desktop/Gabor/Interferometer Sim/Simulation_Result/1. V_x_2 Modulation/Fit_Result_%i.jpg'%m, dpi=300)
        plt.savefig(export_path + 'Fit_Result_%i.jpg'%m, dpi=300)

        plt.close()
        
    #     plt.show()


# plt.figure('Semi-axis length')
# plt.gcf().set_size_inches(18,10)
# 
# plt.subplot(2,4,1)
# plt.text(0, 0, (parameter_note + 'V_x_2 from 0 to %f'%V_x_2_set[-1]), fontsize=10, style='oblique')
# plt.text(0, 0, (parameter_note + 'V_y_2 from 0 to %f'%V_y_2_set[-1]), fontsize=10, style='oblique')
# plt.xticks([])
# plt.yticks([])
# plt.axis('off')
# 
# for i in range(len(dx)*len(dy)-1):
#     plt.subplot(2,4,i+2)
#     slope_a = (a_set[i][-1] - a_set[i][0]) / (V_x_2_set[-1] - V_x_2_set[0])
#     slope_b = (b_set[i][-1] - b_set[i][0]) / (V_x_2_set[-1] - V_x_2_set[0])
#     plt.plot(V_x_2_set, a_set[i], label='semi-major %f'%slope_a, color='red', linestyle='solid', linewidth=1, marker='o', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
#     plt.plot(V_x_2_set, b_set[i], label='semi-minior %f'%slope_b, color='blue', linestyle='solid', linewidth=1, marker='o', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
#     plt.xlabel('V_x_2')
# #     slope_a = (a_set[i][-1] - a_set[i][0]) / (V_y_2_set[-1] - V_y_2_set[0])
# #     slope_b = (b_set[i][-1] - b_set[i][0]) / (V_y_2_set[-1] - V_y_2_set[0])
# #     plt.plot(V_y_2_set, a_set[i], label='semi-major %f'%slope_a, color='red', linestyle='solid', linewidth=1, marker='o', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
# #     plt.plot(V_y_2_set, b_set[i], label='semi-minior %f'%slope_b, color='blue', linestyle='solid', linewidth=1, marker='o', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
# #     plt.xlabel('V_y_2')   
#     plt.legend(loc='best')
# plt.savefig(export_path + '1.jpg', dpi=300)
# plt.close()
# 
plt.figure('Delta')
plt.gcf().set_size_inches(18,10)
for i in range(len(dx)*len(dy)-1):
    plt.subplot(2,4,i+2)
    plt.plot(Delta_set[i], color='blue', linestyle='solid', linewidth=1, marker='o', fillstyle='full', markeredgecolor='red', markeredgewidth=0.0)
    print((Delta_set[i][0]-Delta_set[i][-1])/9)
#     plt.plot(np.arctan2((np.array(b_set[0])), np.array(a_set[0])), 'o')
plt.show()

    