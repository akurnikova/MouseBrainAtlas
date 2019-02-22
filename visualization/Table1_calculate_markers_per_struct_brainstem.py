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

structs = ['Pr5_L', 'Pr5_R', 'DMSp5_L',   'DMSp5_R',  'KF_L', 'KF_R', 'PR5VL_L', 'PR5VL_R', 'SPVmu_L', 'SPVmu_R', 'SpVI_L', 'SpVI_R', 'SpVC_L', 'SpVC_R', 'SpVO_L', 'SpVO_R', 'intertrigeminal_L', 'intertrigeminal_R']
filename = 'Table1C_brainstem_structs_trig.pdf'


'''
structs = ['Amb_L',  'Amb_R', 'IO_L', 'IO_R', 'LRT_L', 'LRT_R', 'Vestibular_L', 'Vestibular_R','RIP', 'ROB', 'RPa','Sol_L', 'Sol_R', 'Mx_L', 'Mx_R']
filename = 'Table1C_brainstem_structs.pdf'

structs = ['SubC_L','SubC_R','PnC_L', 'PnC_R', 'PnO_L', 'PnO_R',
           'PCRt_L','PCRt_R', 'PCRtA_L', 'PCRtA_R','IRt_L', 'IRt_R', 'Gi_L', 'Gi_R', 'LPGi_L','LPGi_R','GiV_L', 'GiV_R', 'GiA', 
           'MdV_L','MdV_R','MdD_L','MdD_R']

filename = 'Table1C_brainstem_structs_reticular.pdf'
'''


structs = ['Pr5_L', 'Pr5_R', 'DMSp5_L',   'DMSp5_R',  'KF_L', 'KF_R', 'PR5VL_L', 'PR5VL_R', 'SPVmu_L', 'SPVmu_R', 'SpVI_L', 'SpVI_R', 'SpVC_L', 'SpVC_R', 'SpVO_L', 'SpVO_R', 'intertrigeminal_L', 'intertrigeminal_R',
           'Amb_L',  'Amb_R', 'IO_L', 'IO_R', 'LRT_L', 'LRT_R', 'Vestibular_L', 'Vestibular_R','RIP', 'ROB', 'RPa','Sol_L', 'Sol_R', 'Mx_L', 'Mx_R',
           'SubC_L','SubC_R','PnC_L', 'PnC_R', 'PnO_L', 'PnO_R',
           'PCRt_L','PCRt_R', 'PCRtA_L', 'PCRtA_R','IRt_L', 'IRt_R', 'Gi_L', 'Gi_R', 'LPGi_L','LPGi_R','GiV_L', 'GiV_R', 'GiA', 
           'MdV_L','MdV_R','MdD_L','MdD_R']
filename = 'Table1C_brainstem_structs_ALL.pdf'

include_stacks = ['RV16','RV15','RV14','RV13','RV4','RV19','RV9','RV10','RV2_L','RV2_R']#,'RV12','RV7','RV8']
spacefill = True

atlas_name = 'Rat_brainstem_atlas'
time_by_stack = {'RV2':72,'RV4':67,'RV7':48,'RV8':48.1,'RV9':53,'RV10':51,'RV12':49,\
                 'RV13':64,'RV14':64.5,'RV15':77,'RV16':77.1,'RV19':61,'RV2_L':72,'RV2_R':72}
cell_count = {}

cc = pd.DataFrame()

    
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
            v,f = load_mesh_stl(fp_spacefill)
        else:
            v,f = load_mesh_stl(fp_av)
            
        polydata =mesh_to_polydata(v,f)
        
        selectEnclosedPoints.SetSurfaceData(polydata)
        selectEnclosedPoints.Update()
        
        numpoints_inside = 0
        for i in range(len(brain_markers_aligned2atlas)):
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
fig = plt.figure(figsize = (20,10))
ax = fig.add_subplot(111)

cc[cc==0] = 0.1
sns.barplot(x = 'Struct',y = 'count',hue = 'Timepoint',data = cc, log=True)


fig_title = '/home/asya/Documents/Yuncong_code/Figures_code/Atlas_figs/'+filename
plt.savefig(fig_title, format='pdf',dpi=fig.dpi)