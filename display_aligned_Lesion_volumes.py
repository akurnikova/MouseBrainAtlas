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
from metadata import *
from data_manager import *
from annotation_utilities import *
from registration_utilities import *
#from conversion import *
from vis3d_utilities import *

from volume_display_utilities import *
from vis3d_utilities_stacy import *
import vtk

import xml.etree.ElementTree as ET

generate_meshes = 0
save_meshes = 0
meshes_from_file = 0
generate_actors = 0
do_slice = 1

#Saggital through lesions
#cut_plane_origin = (0,0,140)
#cut_plane_normal = (0,0,1)

#Saggital through lesions
cut_plane_origin = (-8130,0,0)
cut_plane_normal = (1,0,0)
thickness = 1


include_stacks = {'RFles02','RFles03','RFles04','RFles05','RFles06','RFles07','RFles08'} #{'RV4_67hrs','RV13_64hrs','RV14_65hrs'}
stack_fixed = 'RV14_65hrs'

stacks_to_colors = {'RFles02': (0./255.,200./255.,10./255.), # green
                    'RFles03': (152./255.,78./255.,163./255.), # purple
                    'RFles04': (0./255.,200./255.,10./255.), # green
                    'RFles05':  (152./255.,78./255.,163./255.), # purple
                    'RFles06':  (152./255.,78./255.,163./255.), # purple
                    'RFles07': (0./255.,200./255.,10./255.), # green
                    'RFles08': (0./255.,0./255.,10./255.), # green
                    'BIRTles07': (55./255.,126./255.,184./255.), # blue
                    'BIRTles08': (55./255.,126./255.,184./255.), # blue
                    'BIRTles09': (55./255.,126./255.,184./255.), # blue
                    'BIRTles04': (55./255.,126./255.,184./255.), # blue
                    }

structs_to_colors = {'7n_R':(0,0,0), '7n_L':(0,0,0),'7N_R':(0,0,0),'7N_L':(0,0,0),'5N_R':(0,0,0),\
                     '5N_L':(0,0,0),'LRt_R':(0,0,0),'LRt_L':(0,0,0),'Amb_L':(0,0,0),'Amb_R':(0,0,0),\
                     'IO_L':(0,0,0), 'IO_R':(0,0,0),'Brainstem':(0,0,0)}

if generate_meshes:
    all_meshes, all_origin = generate_RV_meshes(include_stacks,stack_fixed = stack_fixed,warp_setting = 25)
    if save_meshes:
        save_RV_meshes(all_meshes, all_origin, folder = '/home/asya/Documents/Yuncong_code/data/meshes_save/')
elif meshes_from_file:
    all_meshes, all_origin = load_file_RV_meshes(folder = '/home/asya/Documents/Yuncong_code/data/meshes_save/')

if generate_actors:
#    brain_actor_list_fixed = generate_whole_actors(all_meshes,all_origin,structs_to_colors,include_stacks = {stack_fixed})
    brain_actor_list_all_stacks = generate_whole_actors_lesions(all_meshes,all_origin,stacks_to_colors,include_stacks = include_stacks)
    
if do_slice:
    slice_actor_list_fixed = generate_slice_actors(all_meshes,all_origin,structs_to_colors,include_stacks = {stack_fixed},\
                          doBrainstem = 0,wireframe=False,opacity= 0.3,\
                          cut_plane_origin = cut_plane_origin,cut_plane_normal = cut_plane_normal,thickness = thickness)
                          
    slice_actor_list_all_stacks = generate_slice_actors_lesions(all_meshes,all_origin,stacks_to_colors,include_stacks = include_stacks,\
                                  wireframe=False,opacity= 0.3,\
                                  cut_plane_origin = cut_plane_origin,cut_plane_normal = cut_plane_normal,thickness = thickness)


#%%
#vtk_3d_input_list = []
#vtk_3d_input_list.extend([a for actor_list in brain_actor_list_all_stacks.values() for a in actor_list.values()])
#vtk_3d_input_list.extend([a for actor_list in brain_actor_list_fixed.values() for a in actor_list.values()])

vtk_slice_input_list = []
vtk_slice_input_list.extend([a[0] for a in slice_actor_list_all_stacks.values()])

#del slice_actor_list_fixed[stack_fixed]['Brainstem']
vtk_slice_input_list.extend([a[0] for a in slice_actor_list_fixed[stack_fixed].values()])


#vtk_slice_input_list.extend([a for actor_list in brain_actor_list_fixed.values() for a in actor_list.values()])
launch_vtk(vtk_slice_input_list)


#launch_vtk(vtk_3d_input_list)
