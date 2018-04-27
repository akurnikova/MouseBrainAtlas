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
save_meshes = 0
meshes_from_file = 0
generate_actors = 1

slice_thickness = 15


cut_plane_origin = {}
cut_plane_normal = {}

def normalize_normal(a):
    N = (a[0]**2+a[1]**2+a[2]**2)**(0.5)
    return (a[0]/N,a[1]/N,a[2]/N)

cut_plane_origin[0] = (1400,0,90)
cut_plane_normal[0] =  normalize_normal((-0.2,0.,1.))

#cut_plane_origin_2 = (1400,650,0)
#cut_plane_normal_2 = (0.3,1.,0)

cut_plane_origin[1] = (1400,0.,0.)
cut_plane_normal[1] = normalize_normal((1.,0.,0.))

cut_plane_origin[2] = (1480.,0.,0.)
cut_plane_normal[2] =  normalize_normal((1.,0.,0.))
    
include_stacks = {'RV4_67hrs','RV13_64hrs','RV14_65hrs'}

stacks_to_colors = {'RV4_67hrs': (55./255.,126./255.,184./255.), # blue
                    'RV14_65hrs': (255./255.,255./255.,51./255.), # yellow
                    'RV19_61hrs': (152./255.,78./255.,163./255.), # green
                    'RV13_64hrs': (228./255.,26./255.,28./255.), # red
                    'RV9_53hrs':(152./255.,78./255.,163./255.)} # purple

structs_to_colors = {'7n_R':(1,1,0), '7n_L':(1,1,0),'7N_R':(0.5,0,0),'7N_L':(0.5,0,0),'5N_R':(0,0,0.5),'5N_L':(0,0,0.5),'LRt_R':(0,0,0.8),'LRt_L':(0,0,0.8),'Amb_L':(1,.5,.5),'Amb_R':(1,.5,.5),'Brainstem':(0.,0.,0.)}
structs_to_colors_slice = {'7n_R':(0.7,0.7,0), '7n_L':(0.7,0.7,0),'7N_R':(0.5,0,0),'7N_L':(0.5,0,0),'5N_R':(0,0,0.5),'5N_L':(0,0,0.5),'LRt_R':(0,0,0.8),'LRt_L':(0,0,0.8),'Amb_L':(1,.5,.5),'Amb_R':(1,.5,.5),'Brainstem':(0.7,0.7,0.7)}

# set normal length to 1

if generate_meshes:
    all_meshes, all_origin = generate_RV_meshes(include_stacks)
    if save_meshes:
        save_RV_meshes(all_meshes, all_origin, folder = '/home/asya/Documents/Yuncong_code/data/meshes_save_2/')
elif meshes_from_file:
    all_meshes, all_origin = load_file_RV_meshes(folder = '/home/asya/Documents/Yuncong_code/data/meshes_save/')


if generate_actors:
    brain_actor_list_all_stacks = generate_whole_actors(all_meshes,all_origin,structs_to_colors,include_stacks = include_stacks)

## Make slices
brain_slice_actor_list_all_stacks = dict({0:{},1:{},2:{}})
for name_stack, brain_meshes in all_meshes.iteritems():
    if name_stack in include_stacks:
        for i in range(0,3):
            brain_actor_list_slice = {name_s: actor_mesh_cut(mesh, color=structs_to_colors_slice[name_s],
                                                               origin=all_origin[name_stack][name_s], 
                                                               cut_plane_origin = cut_plane_origin[i], 
                                                               cut_plane_normal = cut_plane_normal[i]) 
                                  for name_s, mesh in brain_meshes.iteritems()}
        
            brain_slice_actor_list_all_stacks[i][name_stack] = brain_actor_list_slice

fixed_brain_marker_actors, moving_brain_marker_actors = generate_RV_marker_actors(stacks_to_colors, slice_thickness, include_stacks = include_stacks)

all_brain_marker_actors_slice = dict({0:{},1:{},2:{}})
for i in range(0,3):
    all_brain_marker_actors_slice[i] = generate_RV_marker_actors_slice(cut_plane_normal[i], cut_plane_origin[i], stacks_to_colors,slice_thickness,  include_stacks = include_stacks)

#%%
vtk_3d_input_list = []
vtk_slice_input_list = dict({0:[],1:[],2:[]})

vtk_3d_input_list.extend(fixed_brain_marker_actors)
vtk_3d_input_list.extend([a for actor_list in moving_brain_marker_actors.values() for a in actor_list])
vtk_3d_input_list.extend([a for actor_list in brain_actor_list_all_stacks.values() for a in actor_list.values()])

for i in range(0,3):
    for stack, val in brain_slice_actor_list_all_stacks[i].iteritems():
            vtk_3d_input_list.extend([slice_item[0] for slice_item in val.values()])
            vtk_slice_input_list[i].extend([slice_item[1] for slice_item in val.values()])
    vtk_slice_input_list[i].extend(all_brain_marker_actors_slice[i])


## Add in scale bars on each view
vtk_slice_input_list[0].extend([make_scaleBar_actor(cut_plane_origin[0],cut_plane_normal[0])])
vtk_slice_input_list[1].extend([make_scaleBar_actor(cut_plane_origin[1],cut_plane_normal[1])])
vtk_slice_input_list[2].extend([make_scaleBar_actor(cut_plane_origin[2],cut_plane_normal[2])])

launch_vtk_3view(vtk_3d_input_list,vtk_slice_input_list,cut_plane_normal)
