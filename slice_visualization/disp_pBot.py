#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 15:06:19 2018

@author: asya
"""
import sys
import os
#import time

from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline

sys.path.append('/home/asya/Documents/Yuncong_code/utilities')
from utilities2015 import *
from metadata import *
# from data_manager import *
from annotation_utilities import *
# from registration_utilities import *
# from conversion import *
from vis3d_utilities import *

from volume_display_utilities import *
from byhand_alignment import *
from vis3d_utilities_stacy import *
from neurolucida_to_volume_utilities import *

from matplotlib import pyplot as plt
from shapely.geometry.polygon import Polygon
from descartes import PolygonPatch
import seaborn as sns
import matplotlib.patches as mpatches

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

    
fig_title = '/home/asya/Documents/Yuncong_code/Figures_code/Atlas_figs/PBot_1_S.pdf'

include_stacks = {'RV15','RV16'}
#include_stacks = {'RV4','RV14','RV13','RV19','RV9','RV10'}

dim1 = 0
dim2 = 1

stack_to_color = dict()
for i,stack in enumerate(include_stacks):
    stack_to_color[stack] = sns.color_palette()[i]

fig = figure()
ax = fig.add_subplot(111)


fp = DataManager.get_mesh_filepath(stack_m='Rat_brainstem_atlas', structure='PreBotC_R')
mesh = load_mesh_stl(fp,return_polydata_only=True)

ALL_STRUCT = dict()
ALL_PLANE_O = dict()
ALL_PLANE_N = dict()
i = 0

if dim1 == 0 and dim2 == 1: #Sagittal
    for orig in np.arange(-156,-130,2):
        contour_points,order, _ = polydata_cut(mesh, origin=(0,0,0),
                                           cut_plane_origin = (0,0,orig), 
                                           cut_plane_normal = (0,0,1))
        if len(order) ==0: continue
        print orig
        STRUCT = contour_points[order]
        STRUCT = vstack((STRUCT,STRUCT[0]))
        plt.plot(STRUCT[:,dim1],-STRUCT[:,dim2],'k')
        
        X = STRUCT[:,dim1]
        Y = -STRUCT[:,dim2]
        poly = Polygon(zip(X,Y))
        ALL_STRUCT[i] = poly
        ALL_PLANE_O[i] = (0,0,orig)
        ALL_PLANE_N[i] = (0,0,1)
        i += 1
        
        
if dim1 == 2 and dim2 == 1: #Coronal        
    for orig in np.arange(200,500,10):
        contour_points,order, _ = polydata_cut(mesh, origin=(0,0,0),
                                           cut_plane_origin = (orig,0,0), 
                                           cut_plane_normal = (1,0,0))
        if len(order) ==0: continue
        print orig
        STRUCT = contour_points[order]
        STRUCT = vstack((STRUCT,STRUCT[0]))
        plt.plot(STRUCT[:,dim1],-STRUCT[:,dim2],'k')
        
        X = STRUCT[:,dim1]
        Y = -STRUCT[:,dim2]
        poly = Polygon(zip(X,Y))
        ALL_STRUCT[i] = poly
        ALL_PLANE_O[i] = (orig,0,0)
        ALL_PLANE_N[i] = (1,0,0)
        i += 1
        

plane=vtk.vtkPlane()
thickness = 1

for stack in include_stacks:
    #C,M,D = get_sided_contours(stack)
    count_cells = 0
    tf_parameter_dict = load_alignment_parameters_v2(stack_f='Rat_brainstem_atlas', stack_m=stack, warp_setting=24, 
                                                     vol_type_f='annotationAsScore', vol_type_m='annotationAsScore',
                                                     downscale=15)
    cf = tf_parameter_dict['centroid_f']
    cm = tf_parameter_dict['centroid_m']
    of = tf_parameter_dict['crop_origin_f']
    om = tf_parameter_dict['crop_origin_m']
    params = tf_parameter_dict['params']
    Rt = np.reshape(params, (3,4))
    R = Rt[:3,:3]
    t = Rt[:3,3]
            
    
    markers = bp.unpack_ndarray_file(get_stacy_markers_filepath(stack=stack, structure='All'))
    markers_aligned2fixed = np.dot(R, (markers - om - cm).T).T + t + of + cf
    
    list_to_plot = list()
    for m in markers_aligned2fixed: 
       for i in ALL_STRUCT.keys():
           dist = np.dot(np.asarray(m),ALL_PLANE_N[i])-np.dot(ALL_PLANE_O[i],ALL_PLANE_N[i])
           if abs(dist) < thickness:    
                p = Point(m[dim1],-m[dim2])
                if p.within(ALL_STRUCT[i]):
                    if not ([m[dim1],-m[dim2]]) in list_to_plot:
                        list_to_plot.append([m[dim1],-m[dim2]])
                        count_cells+=1
    for m in list_to_plot:    
        plt.plot(m[0],m[1],color = stack_to_color[stack],marker = '.',linestyle = '')
    print stack
    print count_cells
        

ax.set_aspect('equal')

        
## create legend for the densities
patch_handles = []
for stack in include_stacks:
     patch_handles.append(mpatches.Patch(color=stack_to_color[stack], label='Density %s'%stack))
plt.legend(handles=patch_handles)

plt.savefig(fig_title, format='pdf',dpi=fig.dpi)