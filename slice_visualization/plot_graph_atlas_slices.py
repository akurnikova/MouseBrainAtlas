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
from alpha_hull_filter import *
import vtk

import xml.etree.ElementTree as ET

do_slice = 1
do_volume_display = 1


##Saggital
cut_plane_origin = (0,0,200)
cut_plane_normal = (0,0,1)

##Coronal
cut_plane_origin = (-250,0,0)
cut_plane_normal = (1,0,0)

thickness = 5

include_stacks = {'RV16'}#'RV14'} #
atlas_name = 'Rat_atlas'
spacefill = False


#%%
        
#def generate_slice_coords_atlas(cut_plane_origin = (20,0,0),cut_plane_normal = (0,0,1),thickness = 10,\
#                                doBrainstem = 0,doCortex = 0,opacity= 0.3, structs_to_plot = all_known_structures_sided, constcolor = [],spacefill = False,atlas_name = 'Rat_atlas'):
    
struct_pos = load_pickle(DataManager.get_structure_mean_positions_filepath(atlas_name=atlas_name))
slice_coords_dict_atlas = {}
slice_mesh_dict_atlas = {}

for name_s in all_known_structures_sided:
    ## load in file
    fp = DataManager.get_mesh_filepath(stack_m=atlas_name, structure=name_s)
    fp_spacefill = DataManager.get_mesh_filepath_spacefill(stack_m=atlas_name, structure=name_s)

    if os.path.isfile(fp_spacefill) and spacefill:
        mesh = load_mesh_stl(fp_spacefill,return_polydata_only=True)      
    elif os.path.isfile(fp):
        mesh = load_mesh_stl(fp,return_polydata_only=True)
    else:
        continue

    contour_points,order, P = polydata_cut(mesh, origin=(0,0,0),
                               cut_plane_origin = cut_plane_origin, 
                               cut_plane_normal = cut_plane_normal)
    
        
    slice_coords_dict_atlas[name_s]= contour_points[order]
    slice_mesh_dict_atlas[name_s]  = P

#%%


#%%
fig = pl.figure(figsize=(10,10))
ax = fig.add_subplot(111)

for name_s in slice_coords_dict_atlas.keys():
    name_u = convert_to_original_name(name_s)
    if len(slice_coords_dict_atlas[name_s])==0:
        continue
    print name_u
    if cut_plane_normal == (0,0,1):
        X = slice_coords_dict_atlas[name_s][:,0]
        Y = slice_coords_dict_atlas[name_s][:,1]
    if cut_plane_normal == (1,0,0):
        X = slice_coords_dict_atlas[name_s][:,2]
        Y = -slice_coords_dict_atlas[name_s][:,1]

        
    color_struct = np.array(name_unsided_to_color[name_u])/255.
    if name_u in dotted_structures:
        linestyle = '--'
    else:
        linestyle = '-'
    if name_u in fiber_tracts:
        color_struct = [0.9,0.9,0]
    
    ring = LinearRing(zip(X,Y))
    x, y = ring.xy

    ax.plot(x, y, color=color_struct, alpha=0.7,
        linewidth=1,linestyle = linestyle, solid_capstyle='round', zorder=2)

ax.set_yticks([])
ax.set_xticks([-100,0,100])
ax.set_aspect(1.0)