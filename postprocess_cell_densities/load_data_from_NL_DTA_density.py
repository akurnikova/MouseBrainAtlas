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



DTA_R_stacks = {'DTA03_R',  'DTA06_R',  'DTA05_R', 'DTA10_R','DTA04_R'}
DTA_L_stacks = {'DTA01_L', 'DTA05_L', 'DTA10_L', 'DTA02_L', 'DTA04_L', 'DTA07_L',  'DTA03_L',  'DTA06_L',  'DTA09_L',  'DTA11_L'}

def get_DTA_density(stack, downscale = 15, cut_orientation = 'Saggital'):

    xmlfile = '/home/asya/Documents/Yuncong_code/data/neurolucida_xml/DTA/%s.xml' % stack
    tree = ET.parse(xmlfile)
    root = tree.getroot()
    prefix = '{http://www.mbfbioscience.com/2007/neurolucida}'
    neurolucida_name_map = \
    {'FACIAL MOTOR': '7N_L',
    '7N': '7N_L',
    }
    #########################
    markers = defaultdict(list)

    index = 0
    for item in root.findall(prefix+'marker'):
        curr_markers = []
        for p in item.findall(prefix+'point'):
            curr_markers.append((float(p.attrib['x']), float(p.attrib['y']), float(p.attrib['z'])))
        if index == 0:
            name = 'Blue'
        else:
            name = 'Red'
        index += 1
        markers[name].append(np.array(curr_markers))
        
    markers = {name: np.concatenate(mkr_lists) for name, mkr_lists in markers.iteritems()}
    
    markers_orientationCorrected = {name_u: mkrs3d*[-1,-1,1]
                               for name_u, mkrs3d in markers.iteritems()}

    markers_atlasResol = {name: mkrs3d / (downscale)
                           for name, mkrs3d in markers_orientationCorrected.iteritems()}
    
    ###########
    z_plane = markers_orientationCorrected['Blue'][0][2]
    
    contours = defaultdict(list)

    for item in root.findall(prefix+'contour'):
        name = item.attrib['name']
        if name not in neurolucida_name_map:
            sys.stderr.write('Name %s in stack %s not recognized. Ignored.\n' % (name, stack))
            continue
        name = neurolucida_name_map[name]
        curr_contour = []
        for p in item.findall(prefix+'point'):
            z_curr = float(p.attrib['z'])
            curr_contour.append((float(p.attrib['x']), float(p.attrib['y']), float(p.attrib['z'])))
        if z_curr == z_plane:
            contours[name].append(np.array(curr_contour))

    contours.default_factory = None

    #################

    A, structure_subset = get_list_of_mapped_contours(stack)
    print structure_subset

    ####################

        
    contours_orientationCorrected = {name_u: [cnt*[-1,-1,-1]
                                         for cnt in cnts3d if len(cnt) > 2] 
                               for name_u, cnts3d in contours.iteritems()}

    contours_atlasResol = {name: [cnt / (downscale)
                                    for cnt in cnts3d] 
                           for name, cnts3d in contours_orientationCorrected.iteritems()
                                                if name in structure_subset}

    
    #####################
    
    return contours_atlasResol, markers_atlasResol


for stack in ['DTA03_R']:
    cut_orientation = 's'
        ## 1) Load in contours
    contours, markers = get_DTA_density(stack, downscale = 15, cut_orientation = cut_orientation)