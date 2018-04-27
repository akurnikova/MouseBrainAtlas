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


stacks_to_colors = {'RV4_67hrs': (1,0,0), 'RV14_65hrs': (0,1,0), 'RV19_61hrs': (0,0,1), 'RV13_64hrs': (1,1,0),
                     'DIO6_R_align': (1,1,0), 'DIO6_L_align': (1,1,0),
                     'RFles02': (1,1,0), 'RFles03': (1,1,0)}
structs_to_colors = {'7n_R':(1,1,0), '7n_L':(1,1,0),'7N_R':(0.5,0,0),'7N_L':(0.5,0,0),'5N_R':(0,0,0.5),'5N_L':(0,0,0.5),'LRt_R':(0,0,0.8),'LRt_L':(0,0,0.8),'Amb_L':(1,.5,.5),'Amb_R':(1,.5,.5),\
                     'IO_L':(0,1,1), 'IO_R':(0,1,1), 'Lesion':(0.1,0.1,0.1)}

stack = 'RFles03' #,'RV14_65hrs','RV13_64hrs','RV19_61hrs'}:
    
all_names_in_stack,all_mapped_names = get_list_of_contours(stack)

##load in contours
contours, markers = get_stacy_contours(stack)
vol_bbox_dict = contours_to_volume(contours,stack)

contours, markers = get_stacy_contours(stack)

## 2) Contours to volumes
if contours.has_key('Brainstem'):
    brainstem_contour_to_volume(contours,stack)
vol_bbox_dict = contours_to_volume(contours,stack)

## 3) Save volumes
save_vol_bboxes(vol_bbox_dict,stack)
save_markers(markers,stack)


# Now plot all the stuff     
## Optional from byhand alignment:aligns to origin
#ox,oy,oz = get_new_origin(vol_bbox_dict)
ox,oy,oz = (0,0,0) #get_new_origin(vol_bbox_dict)

polydata_actor_list = []
slice_actor_list = []

## Make polydata actors
for name_s in all_mapped_names:#['7n_R', '7n_L','7N_R','7N_L','LRt_R','LRt_L']: #'5N_R','5N_L'
    if name_s == 'Brainstem':
        continue
    polydata = volume_to_polydata(vol_bbox_dict[name_s][0], num_simplify_iter=3, smooth=True,)

    xmin, _, ymin, _, zmin, _ = vol_bbox_dict[name_s][1]
    polydata_actor = actor_mesh(polydata, color=structs_to_colors[name_s],origin=(xmin-ox,ymin-oy,zmin-oz),opacity = 0.4)
    polydata_actor_list.append(polydata_actor)
       
if not len(markers)==0:
    for m_i in range(0,len(markers['All'])):
       nextpt = markers['All'][m_i,:]-(ox,oy,oz)
       marker_actor = actor_sphere(position= nextpt,radius=2,color=stacks_to_colors[stack])
       polydata_actor_list.append(marker_actor)
        

launch_vtk(polydata_actor_list)
        
