#%load_ext autoreload
#%autoreload 2

import sys
import os
#import time

from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline
import matplotlib.patches as mpatches


sys.path.append('/home/asya/Documents/Yuncong_code/utilities')
sys.path.append('/home/asya/Documents/data/utilities')

from utilities2015 import *
from metadata import *
from data_manager import *
from annotation_utilities import *
from registration_utilities import *
#from conversion import *
from vis3d_utilities import *

from volume_display_utilities import *
from vis3d_utilities_stacy import *
from utils_for_plots import *
import vtk

import xml.etree.ElementTree as ET

from skimage import measure

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

import seaborn as sns


do_volume_display = 1

## Saggital through RF
cut_plane_origin = (0,0,-115)
cut_plane_normal = (0,0,1)

## Coronal through RF
cut_plane_origin = (-75,0,0)
cut_plane_normal = (1,0,0)

## Coronal through nIRT
cut_plane_origin = (-138,0,0)
cut_plane_origin = (-150,0,0)
cut_plane_normal = (1,0,0)

## Saggital through nIRT
#cut_plane_origin = (0,0,-100)
#cut_plane_normal = (0,0,1)

#thickness = 1

include_stacks = {'RV4','RV14','RV13','RV19','RV9','RV10'}

stack_to_color = dict()
for i,stack in enumerate(include_stacks):
    stack_to_color[stack] = sns.color_palette()[i]

atlas_name = 'Rat_brainstem_atlas'
spacefill = True

contour_densities = {0.1} #{0.25,0.5,0.75}

if 0:
    P = dict()
    slice_coords_dict_atlas, temp = load_atlas_slice(atlas_name = 'Rat_brainstem_atlas',cut_plane_origin = cut_plane_origin,cut_plane_normal = cut_plane_normal, spacefill = True)
    for stack in include_stacks:
        fp = DataManager.get_density_pvol_filepath(stack,0.15)
        P[stack] = pickle.load(open(fp,'rb'))

## create legend for the densities
patch_handles = []
for stack in include_stacks:
     patch_handles.append(mpatches.Patch(color=stack_to_color[stack], label='Density %s'%stack))
        
### Plot densities as contours
fig = pl.figure(figsize=(10,10))
ax = fig.add_subplot(111)
plot_atlas_slice(ax,slice_coords_dict_atlas,cut_plane_normal,usecolors = False)

for stack in include_stacks:
    contours = load_density_contours_v2(P[stack],cut_plane_origin,cut_plane_normal, contour_densities = contour_densities)
    plot_density_contours(ax, contours = contours, contour_densities = contour_densities,basecolor = stack_to_color[stack])

ax.set_yticks([])
ax.set_xticks([-175,-150,-125,-100,-75,-50,-25,0,50,100])
ax.set_aspect(1.0)

plt.legend(handles=patch_handles)

### Plot densities as points
fig = pl.figure(figsize=(10,10))
ax = fig.add_subplot(111)
plot_atlas_slice(ax,slice_coords_dict_atlas,cut_plane_normal,usecolors = False)

for stack in include_stacks:
    plot_markers_slice(ax,stack,cut_plane_normal,cut_plane_origin,col = stack_to_color[stack],thickness = thickness)

ax.set_yticks([])
ax.set_xticks([-175,-150,-125,-100,-75,-50,-25,0,50,100])
ax.set_aspect(1.0)

plt.legend(handles=patch_handles)


if do_volume_display:
    slice_actor_dict_atlas = generate_slice_actors_atlas(cut_plane_origin = cut_plane_origin,cut_plane_normal = cut_plane_normal, thickness = thickness,\
                                doBrainstem = 1,doCortex = 0,opacity= 0.2, structs_to_plot = all_known_structures_sided, constcolor = [],spacefill = spacefill,atlas_name = atlas_name)

    structs_to_plot = ['Brainstem','7N','IO','LRT','Amb','5N','SpVI','SpVO','SpVC','PreBotC']

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
    
    launch_vtk([] \
            + structure_mesh_actors_rel2canon,
                   init_angle='sagittal',
            background_color=(1,1,1))
    