#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 22 17:12:08 2018

@author: asya
"""
import numpy as np
import sys
import os
import bloscpack as bp
sys.path.append('/home/asya/Documents/Yuncong_code/utilities')
from utilities2015 import *
from metadata import *
from data_manager import *
from annotation_utilities import *
from registration_utilities import *
from volume_display_utilities import *


from vis3d_utilities_stacy import *
import cPickle as pickle
import vtk
import scipy.ndimage as snd
import copy
from vis3d_utilities import *

from utils_for_plots import *

from shapely.geometry import Polygon

import pandas as pd
import seaborn as sns
import matplotlib.lines as mlines

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

include_stacks = ['RV13','RV4','RV19','RV9','RV10','RV14','RV15','RV16']#,'RV2_L','RV2_R']#,'RV2_L','RV2_R']#,'RV12','RV7','RV8']

time_by_stack = {'RV2':72,'RV4':67,'RV7':48,'RV8':48.1,'RV9':53,'RV10':51,'RV12':49,\
                 'RV13':64,'RV14':64.5,'RV15':77,'RV16':77.1,'RV19':61,'RV2_L':72,'RV2_R':72}


cc = pd.DataFrame()

for stack in include_stacks:
    fp = '/home/asya/Documents/Yuncong_code/data/cell_counts/Count_' + stack + '.pckl'
    new_frame = pd.read_pickle(fp)
    cc = cc.append(new_frame)
                    


for s in cc.struct.unique():
    cc = cc.append(pd.DataFrame([time_by_stack['RV7'],s,0],index=['timepoint','struct','count']).T, ignore_index=True)
    cc = cc.append(pd.DataFrame([time_by_stack['RV8'],s,0],index=['timepoint','struct','count']).T, ignore_index=True)
    cc = cc.append(pd.DataFrame([time_by_stack['RV12'],s,0],index=['timepoint','struct','count']).T, ignore_index=True)
    
## Manual prebotC
cc = cc.drop(1245)
cc = cc.drop(1430)

cc = cc.drop(1060)

cc = cc.drop(1235)
cc = cc.drop(1420)

cc = cc.append(pd.DataFrame([time_by_stack['RV16'],'PreBotC_L',52],index=['timepoint','struct','count']).T, ignore_index=True)
cc = cc.append(pd.DataFrame([time_by_stack['RV15'],'PreBotC_L',53],index=['timepoint','struct','count']).T, ignore_index=True)
cc = cc.append(pd.DataFrame([time_by_stack['RV16'],'PreBotC_R',20],index=['timepoint','struct','count']).T, ignore_index=True)
cc = cc.append(pd.DataFrame([time_by_stack['RV15'],'PreBotC_R',35],index=['timepoint','struct','count']).T, ignore_index=True)
cc = cc.append(pd.DataFrame([time_by_stack['RV14'],'PreBotC_L',4],index=['timepoint','struct','count']).T, ignore_index=True)


cc[cc==0] = 0.1
cc = cc.set_index('struct',drop = False)



#%%#Plotting
'''
fig = plt.figure(figsize = (20,10))
ax = fig.add_subplot(111)

Brainstem = ['Pr5_L',    'DMSp5_L',  'KF_L',              'SpVI_L', 
             'SpVC_L',   'SpVO_L',   'intertrigeminal_L', 'Amb_L',  
             'Vestibular_L','Sol_L',   'Mx_L',
             'PnC_L',   'PnO_L', 'SubC_L',
             'PCRt_L', 'PCRtA_L', 'IRt_L', 'Gi_L', 
             'LPGi_L', 'GiV_L',   'MdV_L', 'MdD_L'] 

sns.barplot(x = 'struct',y = 'count',hue = 'timepoint',data = cc.loc[Brainstem], log=True)#,palette=sns.cubehelix_palette(12,reverse=False))
ax.set_ylim([0.5,1e4])
#sns.set_style("whitegrid")
sns.despine()
filename = 'Bar1_brainstem_structs_L.pdf'
fig_title = '/home/asya/Documents/Yuncong_code/Figures_code/Atlas_figs/'+filename+'_L'
plt.savefig(fig_title, format='pdf',dpi=fig.dpi)


#%%#Plotting
fig = plt.figure(figsize = (20,10))
ax = fig.add_subplot(111)

Brainstem = ['Pr5_R',    'DMSp5_R',  'KF_R',              'SpVI_R', 
             'SpVC_R',   'SpVO_R',   'intertrigeminal_R', 'Amb_R',  
             'Vestibular_R','Sol_R',   'Mx_R',
             'PnC_R',   'PnO_R', 'SubC_R',
             'PCRt_R', 'PCRtA_R', 'IRt_R', 'Gi_R', 
             'LPGi_R', 'GiV_R',   'MdV_R', 'MdD_R'] # 

sns.barplot(x = 'struct',y = 'count',hue = 'timepoint',data = cc.loc[Brainstem], log=True)#,palette=sns.cubehelix_palette(12,reverse=False))
ax.set_ylim([0.5,1e4])
#sns.set_style("whitegrid")
sns.despine()
sns.despine(offset=10, trim=True);
filename = 'Bar1_brainstem_structs_R.pdf'
fig_title = '/home/asya/Documents/Yuncong_code/Figures_code/Atlas_figs/'+filename
plt.savefig(fig_title, format='pdf',dpi=fig.dpi)

'''



#sns.barplot(x = 'struct',y = 'count',hue = 'timepoint',data = cc.loc['7N_L','PreBotC_L'], log=True)
#%%#Plotting

fig = plt.figure(figsize = (20,10))
ax = fig.add_subplot(111)

Other1 =['3N_L', '4N_L', 'AP_L', 'Dk_L', 
        'DR_L', 'DpCe_L', 'Forel_L', 'GP_L',
        'HDB_L', 'LHab_L', 'IC_L', 
        'IMLF_L', 'LDT_L', 'LH_L',
        'LPB_L', 'MCPO_L', 'MPB_L', 'MiTg_L',
        'mRt_L',  'Nlatolfactorytract_L', 'NucleiPCom_L', 'Oculomotor_L']
        
Other2 =['PAG_L', 'PH_L', 'PL_L', 'PMnR_L', 
        'PPTg_L', 'PR_L',  'PaR_L','PreCommisural_L', 
        'RI_L',  'RMg_L', 'RN_L', 
        'RtTg_L', 'SCInG_L', 'SNC_L', 
        'SNL_L', 'SNR_L', 'SPTg_L', 'Su5_L', 
        'Tubercle_L', 'VP_L', 'VTA_L', 'ZI_L']


Other1R =['3N_R', '4N_R', 'AP_R', 'Dk_R', 
        'DR_R', 'DpCe_R', 'Forel_R', 'GP_R',
        'HDB_R', 'LHab_R', 'IC_R', 
        'IMLF_R', 'LDT_R', 'LH_R',
        'LPB_R', 'MCPO_R', 'MPB_R', 'MiTg_R',
        'mRt_R',  'Nlatolfactorytract_R', 'NucleiPCom_R', 'Oculomotor_R']
        
Other2R =['PAG_R', 'PH_R', 'PL_R', 'PMnR_R', 
        'PPTg_R', 'PR_R',  'PaR_R','PreCommisural_R', 
        'RI_R',  'RMg_R', 'RN_R', 
        'RtTg_R', 'SCInG_R', 'SNC_R', 
        'SNL_R', 'SNR_R', 'SPTg_R', 'Su5_R', 
        'Tubercle_R', 'VP_R', 'VTA_R', 'ZI_R']

unsided =['GiA','Cortex','RIP', 'ROB', 'RPa','PreBotC_L','PreBotC_R']
 # 

sns.barplot(x = 'struct',y = 'count',hue = 'timepoint',data = cc.loc[unsided], log=True)
plt.ylim(-1, 10e3)
#filename = 'Table1A_unsided.pdf'
#fig_title = '/home/asya/Documents/Yuncong_code/Figures_code/Atlas_figs/'+filename
#plt.savefig(fig_title, format='pdf',dpi=fig.dpi)


#%%% David's Exponent idea
'''
fig = plt.figure(figsize = (10,5))
ax = fig.add_subplot(111)

color_list = sns.husl_palette(11)
cc[cc==0.1] = 0

df = cc[cc['timepoint'] == 77.1]
df = df[df['count']>0]
plt.plot(log10(np.arange(len(df['count']))),log10(sorted(df['count'],reverse=True)),color = color_list[10])

df = cc[cc['timepoint'] == 77]
df = df[df['count']>0]
plt.plot(log10(sorted(df['count'],reverse=True)),color = color_list[9])

#df = cc[cc['timepoint'] == 72]
#df = df[df['count']>0]
#plt.plot(sorted(df['count'],reverse=True),'c')

df = cc[cc['timepoint'] == 67]
df = df[df['count']>0]
plt.plot(log10(sorted(df['count'],reverse=True)),color = color_list[8])

df = cc[cc['timepoint'] == 64.5]
df = df[df['count']>0]
plt.plot(log10(sorted(df['count'],reverse=True)),color = color_list[7])

df = cc[cc['timepoint'] == 64]
df = df[df['count']>0]
plt.plot(log10(sorted(df['count'],reverse=True)),color = color_list[6])

df = cc[cc['timepoint'] == 61]
df = df[df['count']>0]
plt.plot(log10(sorted(df['count'],reverse=True)),color = color_list[5])

df = cc[cc['timepoint'] == 53]
df = df[df['count']>0]
plt.plot(log10(sorted(df['count'],reverse=True)),color = color_list[4])

df = cc[cc['timepoint'] == 51]
df = df[df['count']>0]
plt.plot(log10(sorted(df['count'],reverse=True)),color = color_list[3])

#fig_title = '/home/asya/Documents/Yuncong_code/Figures_code/Atlas_figs/decay_labelling_log.pdf'
#plt.savefig(fig_title, format='pdf',dpi=fig.dpi)


'''