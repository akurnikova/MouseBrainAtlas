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

include_stacks = {'RV4','RV10','RV13','RV19','RV9','RV14'}#,'RV15','RV16'}
#include_stacks = {}
kernel_bandwidth = 0.5
for stack in include_stacks:
    all_markers = {}
    all_markers[stack] = load_atlas_marker_actors(stack = stack, warp_setting = 24,markertype=str('All'))
    
    density = dict()
    density_vol = {}
    for name_stack, data in all_markers.iteritems():
        kde = st.gaussian_kde(data.T)
        kde.set_bandwidth(kernel_bandwidth)
        # Create a regular 3D grid
        xmin, ymin, zmin = np.min(all_markers[stack],axis=0)
        xmax, ymax, zmax = np.max(all_markers[stack],axis=0)
        xstep,ystep,zstep = 30.,30.,30.
        xNstep = (xmax-xmin)/xstep
        yNstep = (ymax-ymin)/ystep
        zNstep = (zmax-zmin)/zstep
        xi, yi, zi = np.mgrid[xmin:xmax:xNstep*1j, ymin:ymax:yNstep*1j, zmin:zmax:zNstep*1j]
        # Evaluate the KDE on a regular grid...
        coords = np.vstack([item.ravel() for item in [xi, yi, zi]])
        density['vol'] = kde(coords).reshape(xi.shape)
        
            
        fp = DataManager.get_density_pvol_filepath(stack,kernel_bandwidth)#int(round(kernel_bandwidth*1500)))
        create_parent_dir_if_not_exists(fp)
        density['xstep'] = xstep
        density['ystep'] = ystep
        density['zstep'] = zstep
        density['origin'] = (xmin,ymin,zmin)
        density['max'] = (xmax, ymax, zmax)
        pickle.dump(density, open(fp, "wb" ))
    