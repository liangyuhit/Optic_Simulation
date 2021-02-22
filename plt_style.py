# -*- coding: utf-8 -*-
'''
Created on 13.03.2020

@author: yu03
'''
import matplotlib.pyplot as plt
import matplotlib
# from matplotlib.font_manager import _rebuild
# _rebuild() #reload一下
  
color_set_normal = {'black'    : (50, 50, 50),             
            'blue'     : (31, 119, 180),
            'orange'   : (255, 127, 14),
            'green'    : (44, 160, 44),
            'red'      : (214, 39, 40),
            'purple'   : (148, 103, 189),
            'brown'    : (140, 86, 75),
            'pink'     : (227, 119, 194),
             'grey'     : (127, 127, 127),
            'cyan'     : (23, 190, 207),
             }
color_set_medium = {'black'  : (80, 80, 80),
            'blue'     : (114, 158, 206),
            'orange'   : (255, 158, 74),
            'green'    : (103, 191, 92),
            'red'      : (237, 102, 93),
            'purple'   : (173, 139, 201),
            'brown'    : (168, 120, 110),
            'pink'     : (237, 151, 202),
            'grey'     : (162, 162, 162),
            'cyan'     : (109, 204, 218),
            'white'    : (250, 250, 250),
             }
color_set_light = {'black'  : (100, 100, 100),  
            'blue'     : (174, 199, 232),
            'orange'   : (255, 187, 120),
            'green'    : (152, 223, 138),
            'red'      : (255, 152, 150),
            'purple'   : (197, 176, 213),
            'brown'    : (196, 156, 148),
            'pink'     : (247, 182, 210),
             'grey'     : (199, 199, 199),
            'cyan'     : (158, 218, 229),
             }
for color_sets in [color_set_normal, color_set_medium, color_set_light]:
    for i in color_sets:
        r, g, b = color_sets['%s'%i]  
        color_sets['%s'%i] = (r / 255., g / 255., b / 255.)  

font = {'family' : 'serif', 
        'serif' : ['Simsun'], #带衬线 
#         'serif' : ['Times New Roman'],#不带衬线
        'style': 'normal',
        'weight' : 'normal',
        'size'   : 10.5,
        } ###18

matplotlib.rc('font', **font)
plt.rc('font', size=10.5, )          # controls default text sizes
plt.rc('axes', titlesize=10.5,labelsize=10.5, facecolor=color_set_medium['white'])     # fontsize of the axes title edgecolor=color_set_medium['black'],
plt.rc('xtick', labelsize=9, direction='in')    # fontsize of the tick labels
plt.rc('ytick', labelsize=9, direction='in')    # fontsize of the tick labels
plt.rc('legend', fontsize=10.5, framealpha=0.618)    # legend fontsize
plt.rc('figure', titlesize=10.5)  # fontsize of the figure title
plt.rcParams['axes.unicode_minus']=False #负号显示问题
plt.rcParams['mathtext.fontset']='cm'
# plt.rcParams['xtick.color']=color_set_normal['black']
# plt.rcParams['ytick.color']=color_set_normal['black']
plt.rcParams['text.color']='k' #color_set_normal['black']
plt.rcParams['axes.labelcolor']='k' #color_set_normal['black']
plt.rcParams['axes.spines.right']=False
plt.rcParams['axes.spines.top']=False
plt.rcParams['axes.edgecolor']= color_set_medium['grey']

bbox_abcd = {'boxstyle'     : 'round',
            'edgecolor'    : 'none',
#              'facecolor'    : 'none',
            'facecolor'    : color_set_medium['white'],
             'alpha'        : 1, 
             'pad'          : 0.4
             }
bbox_text = {'boxstyle':'round','edgecolor':color_set_medium['grey'],'facecolor': color_set_medium['white'],'alpha':0.8, 'pad': 0.4}