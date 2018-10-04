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

include_stacks = {'GR23_L',
              'GR24_L',
              'GR25_L',
              'GR26_L',
              'GR27_L',
           #   'GR28_L',
              'GR29_L',
              'GR30_L',
              'GR31_L',
             # 'GR32_L',
              'GR33_L',
            #  'GR34_L',
              'GR35_L',
            #  'GR36_L',
              'GR37_L',
              'GR38_L',
            #  'GR39_L',
            #  'GR40_L',
              'GR41_L',
              'GR42_L',
                       }

include_stacks = {'GR35_L'}
kernel_bandwidth = 0.266

all_markers = load_GR_markers(include_stacks=include_stacks, stack_fixed = 'GR35_L',warp_setting = 24, markertype = str('All'))


for stack, data in all_markers.iteritems():
    density = dict()
    kde = st.gaussian_kde(data.T)
    kde.set_bandwidth(kernel_bandwidth)
    # Create a regular 3D grid
    xmin, ymin, zmin = np.min(all_markers[stack],axis=0)
    xmax, ymax, zmax = np.max(all_markers[stack],axis=0)
    xstep,ystep,zstep = 1.,1.,1.
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
    