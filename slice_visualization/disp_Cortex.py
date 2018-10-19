#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 15:06:19 2018

@author: asya
"""
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
# from data_manager import *
from annotation_utilities import *
# from registration_utilities import *
# from conversion import *
from vis3d_utilities import *

from volume_display_utilities import *
from byhand_alignment import *
from vis3d_utilities_stacy import *
from neurolucida_to_volume_utilities import *

from matplotlib import pyplot as plt
from shapely.geometry.polygon import Polygon
from descartes import PolygonPatch
import seaborn as sns
import matplotlib.patches as mpatches

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

    
fig_title = '/home/asya/Documents/Yuncong_code/Figures_code/Atlas_figs/Cortex_1.pdf'


dim1 = 0
dim2 = 1


fp = DataManager.get_mesh_filepath(stack_m='Rat_brainstem_atlas', structure='Cortex')
mesh = load_mesh_stl(fp,return_polydata_only=True)

fp = '/home/asya/Documents/Yuncong_code/data/cell_counts/Coords_' + 'RV16' + '.pckl'
new_frame = pd.read_pickle(fp)
rv16 = new_frame[new_frame['struct']=='Cortex']        

fp = '/home/asya/Documents/Yuncong_code/data/cell_counts/Coords_' + 'RV15' + '.pckl'
new_frame = pd.read_pickle(fp)
rv15 = new_frame[new_frame['struct']=='Cortex']

#%%
actors_list = list()
actors_list.append(actor_mesh(mesh,opacity = 0.1))

scale = 1000./15
for row in rv15.itertuples(index=True, name='Pandas'):
    if getattr(row, "x") < 4: continue
    X = scale*(getattr(row, "x")-1)
    Y = -1*scale*(getattr(row, "y")+0.5)
    Z = scale*(1.2*getattr(row, "z")+0.5)
    actors_list.append(actor_sphere((X,Y,Z),radius = 4,color = (0,1,1)))
    
for row in rv16.itertuples(index=True, name='Pandas'):
    if getattr(row, "x") < 4: continue
    actors_list.append(actor_sphere((scale*(getattr(row, "x")),-scale*(getattr(row, "y")),scale*(getattr(row, "z"))),radius = 4,color = (0,0,0.5)))


A = make_scaleBar_actor((0,0,0),(0,1,0))
actors_list.append(A)
A = make_scaleBar_actor((0,0,0),(0,0,1))
actors_list.append(A)
A = make_scaleBar_actor((0,0,0),(1,0,0))
actors_list.append(A)

launch_vtk(actors_list)