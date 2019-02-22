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

include_stacks = ['RV15']#,RV4','RV19','RV16','RV14','RV13','RV4','RV19','RV9','RV10','RV2_L','RV2_R']#,'RV12','RV7','RV8']

spacefill = False

atlas_name = 'Rat_brainstem_atlas'
time_by_stack = {'RV2':72,'RV4':67,'RV7':48,'RV8':48.1,'RV9':53,'RV10':51,'RV12':49,\
                 'RV13':64,'RV14':64.5,'RV15':77,'RV16':77.1,'RV19':61,'RV2_L':72,'RV2_R':72}

cell_count = {}
cc = pd.DataFrame()
#%%
markers_count = dict()
for s in all_known_structures_sided:
    markers_count[s] = 0
markers_count['DpCe_L'] = 0
markers_count['DpCe_R'] = 0
markers_count['mRt_L'] = 0
markers_count['mRt_R'] = 0
markers_count['LHab_L'] = 0
markers_count['LHab_R'] = 0

scaling_factor = 15./1000

total_cells = 0

#%%## First load markers:
all_markers = dict()
for stack in include_stacks:
    tf_parameter_dict = load_alignment_parameters_v2(stack_f=atlas_name, stack_m=stack, warp_setting=24, 
                                             vol_type_f='annotationAsScore', vol_type_m='annotationAsScore', downscale=32)
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
    all_markers[stack] = brain_markers_aligned2atlas


#%%
thickness = 5
cut_plane_normal = (0.,0.,1.)

testplot = 0
markers_by_struct = pd.DataFrame()

for OZ in range(-440,440,10):
    cut_plane_origin = (0.,0.,OZ)
    
    ######
    ## Repeat for each slice
    ######
    
    plane=vtk.vtkPlane()
    
    ## Load in slice
    slice_coords_dict_atlas_top, slice_coords_dict_atlas_bot = load_atlas_slice_ranges(atlas_name = 'Rat_brainstem_atlas',cut_plane_origin = cut_plane_origin,cut_plane_normal = cut_plane_normal, spacefill = spacefill,thickness = thickness)
    
    
    if testplot == 1:
        fig = pl.figure(figsize=(10,10))
        ax = fig.add_subplot(111)
        ax.set_aspect(1.0)
        plot_atlas_slice(ax,slice_coords_dict_atlas_top,cut_plane_normal,usecolors = True)
        plot_atlas_slice(ax,slice_coords_dict_atlas_bot,cut_plane_normal,usecolors = True)
    ## Create Polygons
    all_polygons_T = dict()
    all_polygons_B = dict()
    
    X_brainstem = []
    Y_brainstem = []
    X_ml = []
    Y_ml = []
    X_ZI = []
    Y_ZI = []
    X_Tu = []
    Y_Tu = []
    X_fr = []
    Y_fr = []
    X_PR = []
    Y_PR = []
    
    
    for name_struct in slice_coords_dict_atlas_top.keys():
        
      #  if name_struct in ['IMLF_L','IMLF_R']:
      #      print 'IMLFFF'
        if name_struct in ['7N_L','7N_R']:
            continue
        if name_struct in ['5N_L','5N_R']:
            continue
        if name_struct in ['LRT_L','LRT_R']:
            continue
        if name_struct in ['IO_L','IO_R']:
            continue
        if name_struct in ['hipp_L','hipp_R']:
            continue
        if name_struct in ['PB_L','PB_R']:
            continue
        
        if name_struct == 'Brainstem':
            X_brainstem = np.hstack((slice_coords_dict_atlas_top[name_struct][:,0],slice_coords_dict_atlas_bot[name_struct][:,0]))
            Y_brainstem = -np.hstack((slice_coords_dict_atlas_top[name_struct][:,1],slice_coords_dict_atlas_bot[name_struct][:,1]))
        if name_struct in ['ml_L','ml_R']:
            X_ml = np.hstack((slice_coords_dict_atlas_top[name_struct][:,0],slice_coords_dict_atlas_bot[name_struct][:,0]))
            Y_ml = -np.hstack((slice_coords_dict_atlas_top[name_struct][:,1],slice_coords_dict_atlas_bot[name_struct][:,1]))
        if name_struct in ['ZI_L','ZI_R']:
            X_ZI = np.hstack((slice_coords_dict_atlas_top[name_struct][:,0],slice_coords_dict_atlas_bot[name_struct][:,0]))
            Y_ZI = -np.hstack((slice_coords_dict_atlas_top[name_struct][:,1],slice_coords_dict_atlas_bot[name_struct][:,1]))
        if name_struct in ['Tubercle_L','Tubercle_R']:
            X_Tu = np.hstack((slice_coords_dict_atlas_top[name_struct][:,0],slice_coords_dict_atlas_bot[name_struct][:,0]))
            Y_Tu = -np.hstack((slice_coords_dict_atlas_top[name_struct][:,1],slice_coords_dict_atlas_bot[name_struct][:,1]))    
        if name_struct in ['fr_L','fr_R']:
            X_fr = np.hstack((slice_coords_dict_atlas_top[name_struct][:,0],slice_coords_dict_atlas_bot[name_struct][:,0]))
            Y_fr = -np.hstack((slice_coords_dict_atlas_top[name_struct][:,1],slice_coords_dict_atlas_bot[name_struct][:,1]))
        if name_struct in ['PR_L','PR_R']:
            X_PR = np.hstack((slice_coords_dict_atlas_top[name_struct][:,0],slice_coords_dict_atlas_bot[name_struct][:,0]))
            Y_PR = np.hstack((-slice_coords_dict_atlas_top[name_struct][:,1],-slice_coords_dict_atlas_bot[name_struct][:,1]))   
        
        if name_struct in fiber_tracts_sided:
            continue
        
        X = slice_coords_dict_atlas_top[name_struct][:,0]
        Y = -slice_coords_dict_atlas_top[name_struct][:,1]
        poly = Polygon(zip(X,Y))
        all_polygons_T[name_struct] = poly
        
        X = slice_coords_dict_atlas_bot[name_struct][:,0]
        Y = -slice_coords_dict_atlas_bot[name_struct][:,1]
        poly = Polygon(zip(X,Y))
        all_polygons_B[name_struct] = poly
    
    
    ######
    ## Repeat for each stack
    ######
    
    ## Count cells
    markers_in_struct = list()
    for i, (x,y,z) in enumerate(all_markers[stack]):
        dist = np.abs(plane.DistanceToPlane((x,y,z),cut_plane_normal,cut_plane_origin))
        x,y,z = scaling_factor*x,-scaling_factor*y,scaling_factor*z ## rescale to mm
        p = Point(x,y)
        if dist < thickness:
            for name_struct in all_polygons_T.keys():
                T = p.within(all_polygons_T[name_struct])
                B = p.within(all_polygons_B[name_struct])
                Tdist = all_polygons_T[name_struct].boundary.distance(p)*(-2*int(T==True)+1)
                Bdist = all_polygons_B[name_struct].boundary.distance(p)*(-2*int(B==True)+1)
                markers_in_struct.append((i, name_struct,T,B,Tdist,Bdist,x,y,z))
    if len(markers_in_struct) == 0:
        continue
    mFr = pd.DataFrame.from_dict(markers_in_struct)
    mFr.columns = (['ncell','name_struct','T','B','Tdist','Bdist','x','y','z'])
    
    for c in mFr['ncell'].unique():
        A = mFr.loc[mFr['ncell'] == c]

        Cerebellum_cell = False
        Brainstem_cell = False
        LH_cell = False
        Habenula_cell = False
        if (A['name_struct']=='Brainstem').any():     
            ## Case 1, cell is part of brainstem:
            
            if A.loc[A['name_struct']=='Brainstem']['T'].values[0] or\
                A.loc[A['name_struct']=='Brainstem']['B'].values[0]:
                A = A[A.name_struct!='Brainstem']
                Brainstem_cell = True
            
                
            ## Case 2, cell is part of Cerebellar nuclei:
            
            if Brainstem_cell == False\
            and unique(A['x'])[0] < max(X_brainstem)-2\
            and unique(A['x'])[0] > min(X_brainstem)+2\
            and not A['T'].any():
                Cerebellum_cell = True
                if OZ < 0:
                    name = 'DpCe_R'
                else:
                    name = 'DpCe_L'
                
            if len(X_ZI) == 0:
                if not len(X_PR) == 0:
                    X_ZI = X_PR
                    Y_ZI = Y_PR
                else:
                    X_ZI = [-100,100]
                    Y_ZI = [-100,100]
            ## Case 3, cell is in the LH:
            
            if Brainstem_cell == False\
            and Cerebellum_cell == False\
            and unique(A['y'])[0] < min(Y_ZI)\
            and unique(A['x'])[0] > min(X_ZI)\
            and unique(A['x'])[0] < max(X_ZI)+2\
            and not A['T'].any():
                LH_cell = True
                if OZ < 0:
                    name = 'LH_R'
                else:
                    name = 'LH_L'
            
            ## Case 5, cell is part of Habenula:
            
            if not len(X_fr) == 0:
                if Brainstem_cell == False\
                and Cerebellum_cell == False\
                and LH_cell == False\
                and unique(A['x'])[0] > max(X_brainstem)\
                and unique(A['x'])[0] > (max(X_fr)-1)\
                and unique(A['x'])[0] < (max(X_fr)+2)\
                and unique(A['y'])[0] > (max(Y_fr))\
                and unique(A['y'])[0] < (max(Y_fr)+1)\
                and not A['T'].any():
                    Habenula_cell = True
                    if OZ < 0:
                        name = 'LHab_R'
                    else:
                        name = 'LHab_L'

        ## Case 4, cell is part of mRt:
        mRt_cell = False
        if len(X_fr) == 0:
            X_fr = [-100,100]
            Y_fr = [-100,100]
        if len(X_PR) == 0:
            X_PR = [-100,100]
            Y_PR = [-100,-50]
        if not len(X_ml) == 0:
            if Brainstem_cell == False\
            and Cerebellum_cell == False\
            and LH_cell == False\
            and Habenula_cell == False\
            and unique(A['x'])[0] > max(X_brainstem)-1\
            and unique(A['y'])[0] > min(Y_ml)\
            and unique(A['y'])[0] > max(Y_PR)\
            and unique(A['x'])[0] < 7\
            and not A['T'].any():
                mRt_cell = True
                if OZ < 0:
                    name = 'mRt_R'
                else:
                    name = 'mRt_L'
        

        ## Case 6, cell is in the VP:
        VP_cell = False
        if not len(X_Tu) == 0:
            if Brainstem_cell == False\
            and Cerebellum_cell == False\
            and LH_cell == False\
            and Habenula_cell == False\
            and mRt_cell == False\
            and abs(OZ) > 80\
            and unique(A['x'])[0] < max(X_Tu)\
            and unique(A['x'])[0] > min(X_Tu)\
            and unique(A['y'])[0] < max(Y_Tu)+2\
            and not A['T'].any():
                VP_cell = True
                if OZ < 0:
                    name = 'VP_R'
                else:
                    name = 'VP_L'
        Spinal_cell = False
        if OZ<40 and unique(A['x'])[0]<-4\
            and not A['T'].any():
                name = 'Spinal_cord'
                continue
            
        if Cerebellum_cell == False\
        and mRt_cell == False\
        and LH_cell == False\
        and Habenula_cell == False\
        and VP_cell == False:
            bot_min = min(A['Bdist'])
            top_min = min(A['Tdist'])
            if bot_min < top_min:
                name = A.loc[A.Bdist==bot_min]['name_struct'].values[0]
            else:
                name = A.loc[A.Tdist==top_min]['name_struct'].values[0]
        
        if testplot == 1:
            col = np.array(name_unsided_to_color[convert_to_original_name(name)])/255.
            x = unique(A['x'])[0]
            y = unique(A['y'])[0]
            ax.plot(x,y,'o',color = col,Markersize = 1)
        
        markers_by_struct = markers_by_struct.append(pd.DataFrame([unique(A['x'])[0],unique(A['y'])[0],unique(A['z'])[0],name]).T,ignore_index = True)
        markers_count[name] += 1
        total_cells+=1
        

cc_stack = pd.DataFrame.from_dict(markers_count,orient='index')
cc_stack['timepoint'] = time_by_stack[stack]
cc_stack = cc_stack.reset_index()
cc_stack.columns = ['struct', 'count','timepoint']

if 1:
    fp = '/home/asya/Documents/Yuncong_code/data/cell_counts/Count_' + stack + '.pckl'
    cc_stack.to_pickle(fp)
    
    
    markers_by_struct.columns = ['x', 'y','z','struct']
    fp = '/home/asya/Documents/Yuncong_code/data/cell_counts/Coords_' + stack + '.pckl'
    markers_by_struct.to_pickle(fp)

