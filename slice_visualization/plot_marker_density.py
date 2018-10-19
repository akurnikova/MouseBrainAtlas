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

stack = 'RV14'
all_markers = {}
all_markers[stack] = load_atlas_marker_actors(stack = 'RV14', warp_setting = 24,markertype=str('All'))


#%% Get density frrom marker data

density = dict()
density_vol = {}
for name_stack, data in all_markers.iteritems():
    kde = st.gaussian_kde(data.T)
    # Create a regular 3D grid
    xmin, ymin, zmin = np.min(all_markers[stack],axis=0)
    xmax, ymax, zmax = np.max(all_markers[stack],axis=0)
    xstep,ystep,zstep = 5.,5.,5.
    xNstep = (xmax-xmin)/xstep
    yNstep = (ymax-ymin)/ystep
    zNstep = (zmax-zmin)/zstep
    xi, yi, zi = np.mgrid[xmin:xmax:xNstep*1j, ymin:ymax:yNstep*1j, zmin:zmax:zNstep*1j]
    # Evaluate the KDE on a regular grid...
    coords = np.vstack([item.ravel() for item in [xi, yi, zi]])
    density[name_stack] = kde(coords).reshape(xi.shape)
    
    ## Save density outlines:
    for d in np.arange(0.2,0.9,0.1):
        density_vol[name_stack] = (density[name_stack]/np.max(density[name_stack]))>d
        imagedata = volume_to_imagedata_respace(np.float32(density_vol[name_stack]), origin=[xmin, ymin, zmin], spacing = [xstep,ystep,zstep])
        surface = vtk.vtkMarchingCubes()
        surface.SetInputData(imagedata)
        surface.ComputeNormalsOn();
        surface.SetValue(0, 0.5);
        surface.Update()
        
        P = surface.GetOutput()
        
        fp = DataManager.get_density_mesh_filepath(stack,d)
        create_parent_dir_if_not_exists(fp)
        save_mesh_stl(P, fp)