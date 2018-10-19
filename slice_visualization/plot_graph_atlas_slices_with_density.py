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

do_slice = 1
do_volume_display = 1


##Saggital
cut_plane_origin = (0,0,-90)
cut_plane_normal = (0,0,1)

##Coronal
#cut_plane_origin = (-200,0,0)
#cut_plane_normal = (1,0,0)

thickness = 20

stack = 'RV10' 
include_stacks = {'RV13','RV9','RV10'}#'RV14'} #
atlas_name = 'Rat_atlas'
spacefill = True

contour_densities = {0.25,0.5,0.75}

#%% 
## Load the stuff
slice_coords_dict_atlas, temp = load_atlas_slice(atlas_name = 'Rat_brainstem_atlas',cut_plane_origin = cut_plane_origin,cut_plane_normal = cut_plane_normal, spacefill = True)

if 1:
    fp = DataManager.get_density_pvol_filepath(stack)
    P = pickle.load(open(fp,'rb'))

#%%
fig = pl.figure(figsize=(10,10))
ax = fig.add_subplot(111)

plot_atlas_slice(ax,slice_coords_dict_atlas,cut_plane_normal,usecolors = False)
load_plot_density(ax,P,cut_plane_origin,cut_plane_normal)

contours = load_density_contours(P,cut_plane_origin,cut_plane_normal, contour_densities = contour_densities)

plot_density_contours(ax, contours = contours, contour_densities = contour_densities,basecolor = [1,0,1])
#plot_density_contours_from_Volume(ax,stack,cut_plane_origin,cut_plane_normal,contour_densities = contour_densities,step = 5)

plot_markers_slice(ax,stack,cut_plane_normal,cut_plane_origin,col = (0.,0.1,0.),thickness = thickness)


ax.set_yticks([])
ax.set_xticks([-100,0,100])
ax.set_aspect(1.0)

## create legend for the densities
basecolor = [1,0,1]
patch_handles = []
for d in contour_densities:
    colorscale = (1.2-d/max(contour_densities))
    color = np.asarray(basecolor)*colorscale
    patch_handles.append(mpatches.Patch(color=color, label='Density %s'%d))
plt.legend(handles=patch_handles)

