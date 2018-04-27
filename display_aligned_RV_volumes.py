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

generate_meshes = 1
save_meshes = 1
meshes_from_file = 0
generate_actors = 1
domarkers = 0

include_stacks = {'RFles03','RFles08'} #{'RV4_67hrs','RV13_64hrs','RV14_65hrs'}
stack_fixed = 'RV14_65hrs'

stacks_to_colors = {'RV4_67hrs': (55./255.,126./255.,184./255.), # blue
                    'RV14_65hrs': (255./255.,255./255.,51./255.), # yellow
                    'RV19_61hrs': (152./255.,200./255.,163./255.), # green
                    'RV13_64hrs': (228./255.,26./255.,28./255.), # red
                    'RV9_53hrs':(152./255.,78./255.,163./255.),  # purple
                    'DIO6_R_align': (55./255.,126./255.,184./255.), # blue
                    'DIO6_L_align': (228./255.,26./255.,28./255.), # red
                    'RFles03': (55./255.,126./255.,184./255.), # blue
                    'RFles02': (228./255.,26./255.,28./255.), # red
                    }

structs_to_colors = {'7n_R':(1,1,0), '7n_L':(1,1,0),'7N_R':(0.5,0,0),'7N_L':(0.5,0,0),'5N_R':(0,0,0.5),'5N_L':(0,0,0.5),'LRt_R':(0,0,0.8),'LRt_L':(0,0,0.8),'Amb_L':(1,.5,.5),'Amb_R':(1,.5,.5),\
                     'IO_L':(0,1,1), 'IO_R':(0,1,1), 'Lesion':(0.1,0.1,0.1)}

if generate_meshes:
    all_meshes, all_origin = generate_RV_meshes(include_stacks,stack_fixed = stack_fixed,warp_setting = 25)
    if save_meshes:
        save_RV_meshes(all_meshes, all_origin, folder = '/home/asya/Documents/Yuncong_code/data/meshes_save_2/')
elif meshes_from_file:
    all_meshes, all_origin = load_file_RV_meshes(folder = '/home/asya/Documents/Yuncong_code/data/meshes_save/')

if generate_actors:
    brain_actor_list_all_stacks = generate_whole_actors(all_meshes,all_origin,structs_to_colors,include_stacks = include_stacks)

if domarkers == 1:
    fixed_brain_marker_actors, moving_brain_marker_actors = generate_RV_marker_actors(stacks_to_colors, include_stacks = include_stacks, stack_fixed = 'DIO6_L_align',warp_setting = 25)
#%%
vtk_3d_input_list = []

if domarkers == 1:
    vtk_3d_input_list.extend(fixed_brain_marker_actors)
    vtk_3d_input_list.extend([a for actor_list in moving_brain_marker_actors.values() for a in actor_list])

vtk_3d_input_list.extend([a for actor_list in brain_actor_list_all_stacks.values() for a in actor_list.values()])


launch_vtk(vtk_3d_input_list)
