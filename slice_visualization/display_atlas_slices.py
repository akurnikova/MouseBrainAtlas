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

do_slice = 1
do_volume_display = 1


##Saggital
cut_plane_origin = (0,0,190)
cut_plane_normal = (0,0,1)

##Coronal
#cut_plane_origin = (-1000,0,0)
#cut_plane_normal = (1,0,0)

thickness = 5

include_stacks = {}#'RV14'} #
atlas_name = 'Rat_brainstem_atlas'
spacefill = False

#structs_to_plot = []


if do_slice:
    slice_actor_dict_atlas = generate_slice_actors_atlas(cut_plane_origin = cut_plane_origin,cut_plane_normal = cut_plane_normal, thickness = thickness,\
                                doBrainstem = 0,doCortex = 0,opacity= 0.3, structs_to_plot = all_known_structures_sided, constcolor = [],spacefill = spacefill,atlas_name = atlas_name)
    
    slice_marker_actors_all_stacks =  generate_RV_marker_actors_slice_for_atlas(cut_plane_normal, cut_plane_origin, slice_thickness = 10, 
                                              include_stacks=include_stacks, radius=2, stacks_to_colors = [])
    

if do_volume_display:
    volume_marker_actors_all_stacks =  generate_RV_marker_actors_for_atlas(include_stacks=include_stacks, radius=2, stacks_to_colors = [])
#%%

if do_slice:
    vtk_slice_input_list = []
    vtk_slice_input_list.extend([a for a in slice_marker_actors_all_stacks])   
    vtk_slice_input_list.extend([[] if a == [] else a[0] for a in slice_actor_dict_atlas.values()])
    
    vtk_slice_input_list = filter(None, vtk_slice_input_list)    
    vtk_slice_input_list.append(make_scaleBar_actor([0,-550,80],cut_plane_normal))
    
#%%
#####
    
if do_volume_display:
    structs_to_plot = ['7N','IO','LRT','Amb','5N','SpVI','SpVO','SpVC','PreBotC']#,'SCInG','SCSuG','IC','ZI','RN','fr','scp','SNR']

    structs_to_plot_sided = list()
    for n in structs_to_plot:
        if n in singular_structures:
            structs_to_plot_sided.append(n)
        else:
            structs_to_plot_sided.append(convert_to_left_name(n))
            structs_to_plot_sided.append(convert_to_right_name(n))
    structure_mesh_actors_rel2canon = list()
    
    for name_s in structs_to_plot_sided:
        fp = DataManager.get_mesh_filepath(stack_m=atlas_name, structure=name_s)
        fp_spacefill = DataManager.get_mesh_filepath_spacefill(stack_m=atlas_name, structure=name_s)
        if os.path.isfile(fp_spacefill) and spacefill:
            mesh_rel2canon = load_mesh_stl(fp_spacefill,return_polydata_only=True)      
        else:
            mesh_rel2canon = load_mesh_stl(fp,return_polydata_only=True)
        A = actor_mesh(mesh_rel2canon,np.array(name_unsided_to_color[convert_to_original_name(name_s)])/255., wireframe=False,opacity = 0.1) 
        structure_mesh_actors_rel2canon.append(A)
        structure_mesh_actors_rel2canon.append(slice_actor_dict_atlas[name_s][1])
        structure_mesh_actors_rel2canon.append(slice_actor_dict_atlas[name_s][2])
        structure_mesh_actors_rel2canon.append(slice_actor_dict_atlas[name_s][3])
    
  # structure_mesh_actors_rel2canon.extend([a for a in volume_marker_actors_all_stacks])

if do_slice and not do_volume_display:
    launch_vtk_slice(vtk_slice_input_list,cut_plane_normal)
if do_volume_display and not do_slice:
    launch_vtk([] \
            + structure_mesh_actors_rel2canon,
                   init_angle='sagittal',
            background_color=(1,1,1))
if do_slice and do_volume_display:
    launch_vtk_multi(structure_mesh_actors_rel2canon, vtk_slice_input_list, cut_plane_normal=cut_plane_normal,background_color=(1.,1.,1.))
