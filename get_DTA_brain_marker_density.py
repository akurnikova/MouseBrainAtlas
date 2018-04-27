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
from registration_utilities_2 import *
#from conversion import *
from vis3d_utilities import *

from vis3d_utilities_stacy import *
from volume_display_utilities import *
import cPickle as pickle
import vtk
import scipy.ndimage as snd
import scipy.stats as st


include_stacks = {'DIO6_R_align','DIO6_L_align'} #{'RV4_67hrs','RV13_64hrs','RV14_65hrs'}
stack_fixed = 'DIO6_L_align'

structs_to_colors_slice = {'7n_R':(0.7,0.7,0), '7n_L':(0.7,0.7,0),'7N_R':(0.5,0,0),'7N_L':(0.5,0,0),'5N_R':(0,0,0.5),'5N_L':(0,0,0.5),'LRt_R':(0,0,0.8),'LRt_L':(0,0,0.8),'Amb_L':(1,.5,.5),'Amb_R':(1,.5,.5),'Brainstem':(0.7,0.7,0.7),'Lesion':(1.,0.,0.)}


all_markers = load_aligned_marker_actors(include_stacks=include_stacks, stack_fixed = stack_fixed,warp_setting = 25, markertype = str('All'))


#%% Get density frrom marker data
calc_density = 1
if calc_density:
    density = dict()
    for name_stack, data in all_markers.iteritems():
        kde = st.gaussian_kde(data.T)
        # Create a regular 3D grid
        xmin, ymin, zmin = np.min(all_markers[stack_fixed],axis=0)
        xmax, ymax, zmax = np.max(all_markers[stack_fixed],axis=0)
        xstep,ystep,zstep = 5.,5.,5.
        xNstep = (xmax-xmin)/xstep
        yNstep = (ymax-ymin)/ystep
        zNstep = (zmax-zmin)/zstep
        xi, yi, zi = np.mgrid[xmin:xmax:xNstep*1j, ymin:ymax:yNstep*1j, zmin:zmax:zNstep*1j]
        # Evaluate the KDE on a regular grid...
        coords = np.vstack([item.ravel() for item in [xi, yi, zi]])
        density[name_stack] = kde(coords).reshape(xi.shape)


#    launch_vtk([actor_volume_respace(np.float32(10e4*density['DIO6_L_align']),'probability',spacing = [xstep,ystep,zstep] ,origin =  [xmin, ymin, zmin])])

#%%
frame = 1
plt.figure
plt.subplot(3,2,1)
plt.imshow(density['DIO6_L_align'][:,:,frame],vmin = 0, vmax = 10e-6 )
plt.title('DIO6_L_align')
plt.subplot(3,2,2)
plt.imshow(density['DIO6_R_align'][:,:,frame],vmin = 0, vmax = 10e-6 )
plt.title('DIO6_R_align')
plt.subplot(3,2,3)
plt.imshow(density['DIO6_R_align'][:,:,0]-density['DIO6_L_align'][:,:,0],vmin = -10e-6, vmax = 10e-6 )
plt.title('difference')
plt.subplot(3,2,4)
plt.imshow(density['DIO6_R_align'][:,:,1]-density['DIO6_L_align'][:,:,1],vmin = -10e-6, vmax = 10e-6 )
plt.title('difference')
plt.subplot(3,2,5)
plt.imshow(density['DIO6_R_align'][:,:,2]-density['DIO6_L_align'][:,:,2],vmin = -10e-6, vmax = 10e-6 )
plt.title('difference')
plt.subplot(3,2,6)
plt.imshow(density['DIO6_R_align'][:,:,3]-density['DIO6_L_align'][:,:,3],vmin = -10e-6, vmax = 10e-6 )
plt.title('difference')

#%%
## Make slices
def normalize_normal(a):
    N = (a[0]**2+a[1]**2+a[2]**2)**(0.5)
    return (a[0]/N,a[1]/N,a[2]/N)
cut_plane_origin = (0,0,50)
cut_plane_normal = normalize_normal((0,0.,1.))

generate_meshes = 1
save_meshes = 0
meshes_from_file = 0
generate_actors =1

slice_thickness = 15
if generate_meshes:
    all_meshes, all_origin = generate_RV_meshes(include_stacks,stack_fixed = stack_fixed, warp_setting = 25)
    if save_meshes:
        save_RV_meshes(all_meshes, all_origin, folder = '/home/asya/Documents/Yuncong_code/data/meshes_save/')
elif meshes_from_file:
    all_meshes, all_origin = load_file_RV_meshes(folder = '/home/asya/Documents/Yuncong_code/data/meshes_save/')


if generate_actors:
    brain_slice_actor_list_all_stacks = dict()
    for name_stack, brain_meshes in all_meshes.iteritems():
        if name_stack in include_stacks:
            brain_actor_list_slice = {name_s: actor_mesh_cut(mesh, color=structs_to_colors_slice[name_s],
                                                                   origin=all_origin[name_stack][name_s], 
                                                                   cut_plane_origin = cut_plane_origin, 
                                                                   cut_plane_normal = cut_plane_normal) 
                                      for name_s, mesh in brain_meshes.iteritems()}
            
            brain_slice_actor_list_all_stacks[name_stack] = brain_actor_list_slice

#%%
vtk_slice_input_list = []
for stack, val in brain_slice_actor_list_all_stacks.iteritems():
    vtk_slice_input_list.extend([slice_item[1] for slice_item in val.values()])

D = density['DIO6_R_align']-density['DIO6_L_align']
P =actor_mesh_cut_probability(np.float32(5*10e4*D),spacing = [xstep,ystep,zstep] ,origin =  [xmin, ymin, zmin], color=(1.,1.,1.),\
                   cut_plane_origin = cut_plane_origin, cut_plane_normal = cut_plane_normal)
vtk_slice_input_list.extend([P])
launch_vtk(vtk_slice_input_list)
