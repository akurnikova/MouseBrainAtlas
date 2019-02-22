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

import pandas as pd
import seaborn as sns
import matplotlib.lines as mlines

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

structs =['3N_L', '3N_R', '4N_L', '4N_R', 'AP_L', 'AP_R', 'Dk_L', 'Dk_R',
 'DR_L', 'DR_R', 'Forel_L', 'Forel_R', 'GP_L', #'GP_R',
 'HDB_L', #'HDB_R',
 'IC_L', 'IC_R', 'IMLF_L', 'IMLF_R', 
 'LDT_L', 'LDT_R', 'LH_L', 'LH_R',
 'LPB_L', 'LPB_R', 'MCPO_L','MPB_L', 'MPB_R', 'MiTg_L', 'MiTg_R', 
# 'Nlatolfactorytract_L', 'Nlatolfactorytract_R',
 'NucleiPCom_L', 'NucleiPCom_R', 'Oculomotor_L', 'Oculomotor_R',
 'PAG_L', 'PAG_R', 'PB_L', 'PB_R', 'PH_L', #'PH_R', 
 'PL_L', 'PL_R',
 'PMnR_L', 'PMnR_R', 'PPTg_L', 'PPTg_R',
 'PR_L', 'PR_R', #'PaR_L', 'PaR_R', #'PreCommisural_L', 'PreCommisural_R',
# 'RI_L', 'RI_R', 
 'RMg_L', 'RMg_R', 'RN_L', 'RN_R', 'RtTg_L', 'RtTg_R',
 'SCInG_L', 'SCInG_R',# 'SCSuG_L', 'SCSuG_R', 'SNC_L', 'SNC_R', 
 'SNL_L', 'SNL_R',
 'SNR_L', 'SNR_R', 'SPTg_L', 'SPTg_R', 'Su5_L', 'Su5_R', 'Tubercle_L', 'VLL_L', 'VLL_R',
 'VP_L', 'VTA_L', 'VTA_R', 'ZI_L', 'ZI_R',
  'Cortex']

structs = ['Cortex','SCInG_R']

structs
include_stacks = ['RV16','RV15','RV14','RV13','RV4','RV19','RV9','RV10','RV2_L','RV2_R']#,'RV12','RV7','RV8']
spacefill = False

atlas_name = 'Rat_brainstem_atlas'
time_by_stack = {'RV2':72,'RV4':67,'RV7':48,'RV8':48.1,'RV9':53,'RV10':51,'RV12':49,\
                 'RV13':64,'RV14':64.5,'RV15':77,'RV16':77.1,'RV19':61,'RV2_L':72,'RV2_R':72}
cell_count = {}

cc = pd.DataFrame()

actors = []

for stack in include_stacks:
    
    ## First load markers:
    
    tf_parameter_dict = load_alignment_parameters_v2(stack_f=atlas_name, stack_m=stack, warp_setting=24, 
                                             vol_type_f='annotationAsScore', vol_type_m='annotationAsScore', downscale=15)
    cf = tf_parameter_dict['centroid_f']
    cm = tf_parameter_dict['centroid_m']
    of = tf_parameter_dict['crop_origin_f']
    om = tf_parameter_dict['crop_origin_m']
    params = tf_parameter_dict['params']
    Rt = np.reshape(params, (3,4))
    R = Rt[:3,:3]
    t = Rt[:3,3]

    moving_brain_markers_raw = bp.unpack_ndarray_file(get_stacy_markers_filepath(stack=stack, structure='All'))
    brain_markers_aligned2atlas = np.dot(R, (moving_brain_markers_raw - om - cm).T).T + t + of + cf
    
    points = vtk.vtkPoints()
    for pt in brain_markers_aligned2atlas:
      points.InsertNextPoint(pt)
    
    pointsPolydata = vtk.vtkPolyData()
    pointsPolydata.SetPoints(points)
     
    selectEnclosedPoints = vtk.vtkSelectEnclosedPoints()
    selectEnclosedPoints.SetInputData(pointsPolydata)

    cell_count_single = {}
    
    for s in structs:
        fp_av = DataManager.get_mesh_filepath(stack_m=atlas_name, structure=s)    
        fp_spacefill = DataManager.get_mesh_filepath_spacefill(stack_m=atlas_name, structure=s)
        
        if not os.path.isfile(fp_av):
            print s
            continue
        if os.path.isfile(fp_spacefill) and spacefill:
            polydata = load_mesh_stl(fp_spacefill,return_polydata_only=True)
        else:
            polydata = load_mesh_stl(fp_av,return_polydata_only=True)
                
        ### TEST DISPLAY
        if stack == 'RV14':
#            moving_brain_markers_aligned2atlas_actors = [actor_sphere((x,y,z), radius = 2, color = (1.,0,0)) 
#                                                        for x,y,z in brain_markers_aligned2atlas]
#            actors.extend(a for a in moving_brain_markers_aligned2atlas_actors)
            A = actor_mesh(polydata,(0,0,0), wireframe=False,opacity = 0.2) 
            actors.append(A)
    
        selectEnclosedPoints.SetSurfaceData(polydata)
        selectEnclosedPoints.Update()
        
        numpoints_inside = 0
        for i in range(len(brain_markers_aligned2atlas)):
            if stack == 'RV14' and selectEnclosedPoints.IsInside(i):
                if s == 'Cortex':
                    pt_Actor = [actor_sphere((brain_markers_aligned2atlas[i]), radius = 6, color = (1.,0,0))]
                else:
                    pt_Actor = [actor_sphere((brain_markers_aligned2atlas[i]), radius = 6, color = (0.,0,1.))]
                actors.extend(pt_Actor)
                print brain_markers_aligned2atlas[i]
            numpoints_inside = numpoints_inside + selectEnclosedPoints.IsInside(i)
        cell_count_single[s] = numpoints_inside
        
        cc = cc.append(pd.DataFrame([time_by_stack[stack],s,numpoints_inside],index=['Timepoint','Struct','count']).T, ignore_index=True)
            
    cell_count[stack] = cell_count_single

## Append the stacks with no labelling
for s in structs:
    cc = cc.append(pd.DataFrame([time_by_stack['RV7'],s,0],index=['Timepoint','Struct','count']).T, ignore_index=True)
    cc = cc.append(pd.DataFrame([time_by_stack['RV8'],s,0],index=['Timepoint','Struct','count']).T, ignore_index=True)
    cc = cc.append(pd.DataFrame([time_by_stack['RV12'],s,0],index=['Timepoint','Struct','count']).T, ignore_index=True)
        
#%% Plotting
fig = plt.figure(figsize = (10,5))
ax = fig.add_subplot(111)

cc[cc==0] = 0.1
sns.barplot(x = 'Struct',y = 'count',hue = 'Timepoint',data = cc, log=True)

#fig_title = '/home/asya/Documents/Yuncong_code/Figures_code/Atlas_figs/Table1D_notbrainstem.pdf'
#plt.savefig(fig_title, format='pdf',dpi=fig.dpi)