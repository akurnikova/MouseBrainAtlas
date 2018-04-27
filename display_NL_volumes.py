#%load_ext autoreload
#%autoreload 2

import sys
import os
#import time

from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline

sys.path.append('/home/asya/Documents/Yuncong_code/utilities')
from utilities2015 import *
# from metadata import *
# from data_manager import *
from annotation_utilities import *
# from registration_utilities import *
# from conversion import *
from vis3d_utilities import *

from volume_display_utilities import *
from byhand_alignment import *
from vis3d_utilities_stacy import *
from neurolucida_to_volume_utilities import *


load_volumes = 1
save_volumes = 1
display_volumes = 0

stack = 'BIRTles07' #,'RV14_65hrs','RV13_64hrs','RV19_61hrs'}:
cut_orientation = 'Sag'


stacks_to_colors = {'RV4_67hrs': (1,0,0), 'RV14_65hrs': (0,1,0), 'RV19_61hrs': (0,0,1), 'RV13_64hrs': (1,1,0),
                     'DIO6_R_align': (1,1,0), 'DIO6_L_align': (1,1,0),
                     'RFles02': (1,1,0), 'RFles03': (1,1,0)}
structs_to_colors = {'7n_R':(1,1,0), '7n_L':(1,1,0),'7N_R':(0.5,0,0),'7N_L':(0.5,0,0),'5N_R':(0,0,0.5),'5N_L':(0,0,0.5),'LRt_R':(0,0,0.8),'LRt_L':(0,0,0.8),'Amb_L':(1,.5,.5),'Amb_R':(1,.5,.5),\
                     'IO_L':(0,1,1), 'IO_R':(0,1,1), 'Lesion':(0.1,0.1,0.1)}

if load_volumes:    
    all_names_in_stack,all_mapped_names = get_list_of_contours(stack)
    
    ##load in contours
    contours, markers = get_stacy_contours(stack, cut_orientation = cut_orientation)
    
    ## 2) Contours to volumes
    if cut_orientation == 'Horizontal' or cut_orientation == 'H':
        vol_bbox_dict = contours_to_volume(contours,stack,interpolation_direction = 'y')
        if contours.has_key('Brainstem'):
            brainstem_contour_to_volume(contours,stack,interpolation_direction = 'y')
    else:
        vol_bbox_dict = contours_to_volume(contours,stack,interpolation_direction = 'z')
        if contours.has_key('Brainstem'):
            brainstem_contour_to_volume(contours,stack,interpolation_direction = 'z')

    
            
if save_volumes:
    ## 3) Save volumes
    save_vol_bboxes(vol_bbox_dict,stack)
    save_markers(markers,stack)

vol_bbox_dict_fixed = copy.deepcopy(vol_bbox_dict)
if display_volumes:
    ## Optional from byhand alignment:aligns to origin
    #ox,oy,oz = get_new_origin(vol_bbox_dict)
    ox,oy,oz = (0,0,0) #get_new_origin(vol_bbox_dict)
    
    ## Make polydata actors
    polydata_actor_list = []
    slice_actor_list = []
    for name_s in ['7n_R', '7n_L','7N_R','7N_L','LRt_R','LRt_L','Amb_R','Amb_L','Lesion']: #all_mapped_names:#['7n_R', '7n_L','7N_R','7N_L','LRt_R','LRt_L']: #'5N_R','5N_L'
        
        ### Fixed brain if loaded
        
        if name_s == 'Brainstem':
            continue
        
        if vol_bbox_dict_fixed.has_key(name_s):
            polydata = volume_to_polydata(vol_bbox_dict_fixed[name_s][0], num_simplify_iter=3, smooth=True,)
        
            xmin, _, ymin, _, zmin, _ = vol_bbox_dict_fixed[name_s][1]
            if not structs_to_colors.has_key(name_s):
                structs_to_colors[name_s] = (0.,0.,0.)
            polydata_actor = actor_mesh(polydata, color=structs_to_colors[name_s],origin=(xmin-ox,ymin-oy,zmin-oz),opacity = 0.4)
            polydata_actor_list.append(polydata_actor)
        
        ### newly loaded brain
        
        if not vol_bbox_dict.has_key(name_s):
            print '%s not in volume' %name_s
            continue
        if name_s == 'Brainstem':
            continue
        
        polydata = volume_to_polydata(vol_bbox_dict[name_s][0], num_simplify_iter=3, smooth=True,)
    
        xmin, _, ymin, _, zmin, _ = vol_bbox_dict[name_s][1]
        if not structs_to_colors.has_key(name_s):
            structs_to_colors[name_s] = (0.,0.,0.)
        polydata_actor = actor_mesh(polydata, color=structs_to_colors[name_s],origin=(xmin-ox,ymin-oy,zmin-oz),opacity = 0.4)
        polydata_actor_list.append(polydata_actor)
           
    if not len(markers)==0:
        for m_i in range(0,len(markers['All'])):
           nextpt = markers['All'][m_i,:]-(ox,oy,oz)
           marker_actor = actor_sphere(position= nextpt,radius=2,color=stacks_to_colors[stack])
           polydata_actor_list.append(marker_actor)
        

    launch_vtk(polydata_actor_list)
        
