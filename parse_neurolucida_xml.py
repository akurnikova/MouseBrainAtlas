#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 11:36:58 2017

@author: yuncong
rearranged functions by stacy
"""

import numpy as np
import sys
import os
import bloscpack as bp
sys.path.append('/home/asya/Documents/Yuncong_code/utilities')
from utilities2015 import *
from metadata import *
from data_manager import *
from annotation_utilities import *
from registration_utilities import *
from vis3d_utilities import *

from vis3d_utilities_stacy import *
from neurolucida_to_volume_utilities import *




stack = 'RFles08'
# cut_direction = 'Horizontal'

## 1) Load in contours
contours, markers = get_stacy_contours(stack)

## 2) Contours to volumes
if contours.has_key('Brainstem'):
    brainstem_contour_to_volume(contours,stack)
vol_bbox_dict = contours_to_volume(contours,stack)

## 3) Save volumes
save_vol_bboxes(vol_bbox_dict,stack)
save_markers(markers,stack)
