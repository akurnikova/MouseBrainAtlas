#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 18 14:25:18 2018

@author: asya
"""

import sys
import os
import time
import matplotlib.pyplot as plt
import numpy as np
sys.path.append('/home/asya/Documents/Yuncong_code/utilities')
from utilities2015 import *
from registration_utilities import *
from annotation_utilities import *
from metadata import *
from data_manager import *
from pandas import read_csv
from volume_display_utilities import *



def generate_aligner_parameters_rigid(warp_setting, alignment_name_dict):
        ## LINE 1
    warp_properties = registration_settings.loc[warp_setting]
    print warp_properties
    
    upstream_warp_setting = warp_properties['upstream_warp_id']
    if upstream_warp_setting == 'None':
        upstream_warp_setting = None
    else:
        upstream_warp_setting = int(upstream_warp_setting)
        
    transform_type = warp_properties['transform_type']
    terminate_thresh = warp_properties['terminate_thresh']
    grad_computation_sample_number = warp_properties['grad_computation_sample_number']
    grid_search_sample_number = warp_properties['grid_search_sample_number']
    std_tx_um = warp_properties['std_tx_um']
    std_ty_um = warp_properties['std_ty_um']
    std_tz_um = warp_properties['std_tz_um']
    std_tx = std_tx_um/(XY_PIXEL_DISTANCE_LOSSLESS*32)
    std_ty = std_ty_um/(XY_PIXEL_DISTANCE_LOSSLESS*32)
    std_tz = std_tz_um/(XY_PIXEL_DISTANCE_LOSSLESS*32)
    std_theta_xy = np.deg2rad(warp_properties['std_theta_xy_degree'])
    
    try:
        surround_weight = float(warp_properties['surround_weight'])
        include_surround = surround_weight != 0 and not np.isnan(surround_weight)
    except:
        surround_weight = str(warp_properties['surround_weight'])
        include_surround = True
    
    reg_weight = warp_properties['regularization_weight']
    if np.isnan(reg_weight):
        reg_weights = np.zeros((3,))
    else:
        reg_weights = np.ones((3,))*reg_weight
    
    MAX_ITER_NUM = 10000
    HISTORY_LEN = 200
    lr1 = 10
    lr2 = 0.1
        
def generate_aligner_parameters(warp_setting, alignment_name_dict):
    """
    Args:
        warp_setting (int):
        alignment_name_dict (dict):
        
    Returns:
        - 'volume_moving': dict {ind_m: 3d array},
        - 'volume_fixed': dict {ind_m: 3d array},
        - 'structure_to_label_moving': dict {str: int},
        - 'label_to_structure_moving': dict {int: str},
        - 'structure_to_label_fixed': dict {str: int}, 
        - 'label_to_structure_fixed': dict {int: str},
        - 'label_weights_m': dict {int: float},
        - 'label_mapping_m2f': dict {int: int},
    """
    
    stack_m = alignment_name_dict['stack_m']
    stack_f = alignment_name_dict['stack_f']
    warp_setting = alignment_name_dict['warp_setting']
    vol_type_f = alignment_name_dict['vol_type_f']
    vol_type_m = alignment_name_dict['vol_type_m']
    
    ###Commented stacy
#    registration_settings = read_csv(REGISTRATION_SETTINGS_CSV, header=0, index_col=0)
#    warp_properties = registration_settings.to_dict()[str(warp_setting)]
    registration_settings = read_csv(REGISTRATION_SETTINGS_CSV, header=0, index_col=0)
    warp_properties = registration_settings.loc[warp_setting]
    print warp_properties


    ################################################################
    
    upstream_warp_setting = warp_properties['upstream_warp_id']
    if upstream_warp_setting == 'None':
        upstream_warp_setting = None
    else:
        upstream_warp_setting = int(upstream_warp_setting)
    
    transform_type = warp_properties['transform_type']
    terminate_thresh = warp_properties['terminate_thresh']
    grad_computation_sample_number = int(warp_properties['grad_computation_sample_number'])
    grid_search_sample_number = warp_properties['grid_search_sample_number']
    std_tx_um = warp_properties['std_tx_um']
    std_ty_um = warp_properties['std_ty_um']
    std_tz_um = warp_properties['std_tz_um']
    std_tx = std_tx_um/15
    std_ty = std_ty_um/15
    std_tz = std_tz_um/15
    std_theta_xy = np.deg2rad(warp_properties['std_theta_xy_degree'])


    surround_weight = warp_properties['surround_weight']
    if isinstance(surround_weight, float) or isinstance(surround_weight, int):
        surround_weight = float(surround_weight)
        include_surround = surround_weight != 0 and not np.isnan(surround_weight)
    elif isinstance(surround_weight, str):
        surround_weight = str(surround_weight)
        # Setting surround_weight as inverse is very important. Using -1 often gives false peaks.
        include_surround = True
        
    ## Added STACY 5/10/18
    #include_surround = False

    print surround_weight, include_surround

    reg_weight = warp_properties['regularization_weight']
    if np.isnan(reg_weight):
        reg_weights = np.zeros((3,))
    else:
        reg_weights = np.ones((3,))*reg_weight

    
    positive_weight = 'size'
    # positive_weight = 'inverse'
    
    structure_subset = ['RN_R']#alignment_structures_left
    #include_surround = False
 
    if include_surround:
        structure_subset = structure_subset + [convert_to_surround_name(s, margin=100) for s in structure_subset]
        
    ############################################################################
    
    volume_moving, volume_moving_bbox, structure_to_label_moving, label_to_structure_moving = \
    DataManager.load_original_volume_all_known_structures_v2(stack=stack_m, sided=True, 
                                                          volume_type=vol_type_m, 
                                                          include_surround=include_surround,
                                                          structures=structure_subset,
                                                            return_label_mappings=True, 
                                                             name_or_index_as_key='index',
                                                             common_shape=True)
    if len(volume_moving) == 0:
        sys.stderr.write("No moving volumes.\n")
    else:
        sys.stderr.write("Loaded moving volumes: %s.\n" % sorted(structure_to_label_moving.keys()))
    
    #############################################################################
    
    volume_fixed, volume_fixed_bbox, structure_to_label_fixed, label_to_structure_fixed = \
    DataManager.load_original_volume_all_known_structures_v2(stack=stack_f, sided=True, 
                                                          volume_type=vol_type_f,
                                                         structures=structure_subset,
                                                         return_label_mappings=True, 
                                                             name_or_index_as_key='index',
                                                         common_shape=True)
    if len(volume_fixed) == 0:
        sys.stderr.write("No fixed volumes.\n")
    else:
        sys.stderr.write("Loaded fixed volumes: %s.\n" % sorted(structure_to_label_fixed.keys()))
            
    ############################################################################
    
    
    label_mapping_m2f = {label_m: structure_to_label_fixed[name_m] 
                     for label_m, name_m in label_to_structure_moving.iteritems()
                     if name_m in structure_subset and name_m in structure_to_label_fixed}
    
    t = time.time()
    cutoff = .5 # Structure size is defined as the number of voxels whose value is above this cutoff probability.
#     pool = Pool(NUM_CORES)
#     volume_moving_structure_sizes = dict(zip(volume_moving.keys(), 
#                                              pool.map(lambda l: np.count_nonzero(volume_moving[l] > cutoff), 
#                                                       volume_moving.keys())))
#     pool.close()
#     pool.join()
    volume_moving_structure_sizes = dict(zip(volume_moving.keys(), 
                                             map(lambda l: np.count_nonzero(volume_moving[l] > cutoff), 
                                                      volume_moving.keys())))
    sys.stderr.write("Computing structure sizes: %.2f s\n" % (time.time() - t))
    
    label_weights_m = {}

    for label_m in label_mapping_m2f.iterkeys():
        name_m = label_to_structure_moving[label_m]
        if not is_surround_label(name_m):
            if positive_weight == 'size':
                label_weights_m[label_m] = 1.
            elif positive_weight == 'inverse':
                p = np.percentile(volume_moving_structure_sizes.values(), 50)
                label_weights_m[label_m] =  np.minimum(p / volume_moving_structure_sizes[label_m], 1.)
            else:
                sys.stderr.write("positive_weight %s is not recognized. Using the default.\n" % positive_weight)

    for label_m in label_mapping_m2f.iterkeys():
        name_m = label_to_structure_moving[label_m]
        if is_surround_label(name_m):
            label_ns = structure_to_label_moving[convert_to_nonsurround_name(name_m)]
            if surround_weight == 'inverse':
    #             label_weights_m[label_m] = - label_weights_m[label_ns] * volume_moving_structure_sizes[label_ns]/float(volume_moving_structure_sizes[label_m])
                label_weights_m[label_m] = label_weights_m[label_ns] * volume_moving_structure_sizes[label_ns]/float(volume_moving_structure_sizes[label_m])
            elif isinstance(surround_weight, int) or isinstance(surround_weight, float):
                label_weights_m[label_m] = surround_weight
            else:
                sys.stderr.write("surround_weight %s is not recognized. Using the default.\n" % surround_weight)
                
    ######################################################
    
    alinger_parameters = \
    {'label_weights_m': label_weights_m,
     'label_mapping_m2f': label_mapping_m2f,
     'volume_moving': volume_moving,
     'volume_fixed': volume_fixed,
     'structure_to_label_moving': structure_to_label_moving,
     'label_to_structure_moving': label_to_structure_moving,
     'structure_to_label_fixed': structure_to_label_fixed, 
     'label_to_structure_fixed': label_to_structure_fixed,
     'transform_type': transform_type,
     'grad_computation_sample_number': grad_computation_sample_number,
     'volume_moving_bbox': volume_moving_bbox,
     'volume_fixed_bbox': volume_fixed_bbox
    }
                
    return alinger_parameters

def compute_gradient(volumes, smooth_first=False):
    """
    Args:
        volumes (dict {int: 3d-array})
        smooth_first (bool): If true, smooth each volume before computing gradients. 
        This is useful if volume is binary and gradients are only nonzero at structure borders.
        
    Note:
        # 3.3 second - re-computing is much faster than loading
        # .astype(np.float32) is important;
        # Otherwise the score volume is type np.float16, np.gradient requires np.float32 and will have to convert which is very slow
        # 2s (float32) vs. 20s (float16)
    """
    gradients = {}
    for ind, v in volumes.iteritems():
        if smooth_first:
            gy_gx_gz = np.gradient(gaussian(v, 3).astype(np.float32), 3, 3, 3)
        else:
            gy_gx_gz = np.gradient(v.astype(np.float32), 3, 3, 3)
        gradients[ind] = np.array([gy_gx_gz[1], gy_gx_gz[0], gy_gx_gz[2]])
    return gradients


def load_aligned_marker_actors(include_stacks={'RV4_67hrs','RV19_61hrs','RV14_65hrs','RV13_64hrs'}, stack_fixed = 'RV4_67hrs',warp_setting = 24,markertype='All'):
    # Load fixed brains markers.
    
    fixed_brain_markers = bp.unpack_ndarray_file(get_stacy_markers_filepath(stack=stack_fixed, structure=markertype))
    
    
        
    if stack_fixed not in include_stacks:
        all_markers = dict()
        raw_markers = dict()
        lesion_markers = dict()
        aligned_lesion_markers = dict()
    else:
        all_markers = {stack_fixed:fixed_brain_markers}
        raw_markers = {stack_fixed:(fixed_brain_markers - np.min(fixed_brain_markers,axis = 0))}
        fn = get_stacy_markers_filepath(stack=stack_fixed, structure='1')
        if os.path.isfile(fn):
            lesion_brain_markers = bp.unpack_ndarray_file(fn)
            lesion_markers = {stack_fixed: (lesion_brain_markers - np.min(fixed_brain_markers,axis = 0))}
            aligned_lesion_markers = {stack_fixed: lesion_brain_markers}
        else:
            lesion_markers = dict()
            aligned_lesion_markers = dict()
        
    
    
    for stack_moving in include_stacks:
        if stack_moving == stack_fixed: continue
        tf_parameter_dict = load_alignment_parameters_v2(stack_f=stack_fixed, stack_m=stack_moving, warp_setting=warp_setting, 
                                                 vol_type_f='annotationAsScore', vol_type_m='annotationAsScore',
                                                 downscale=32)
        cf = tf_parameter_dict['centroid_f']
        cm = tf_parameter_dict['centroid_m']
        of = tf_parameter_dict['crop_origin_f']
        om = tf_parameter_dict['crop_origin_m']
        params = tf_parameter_dict['params']
        Rt = np.reshape(params, (3,4))
        R = Rt[:3,:3]
        t = Rt[:3,3]
    
        moving_brain_markers = bp.unpack_ndarray_file(get_stacy_markers_filepath(stack=stack_moving, structure='All'))
        
        moving_brain_markers_aligned2fixed = moving_brain_markers - om - cm + t + of + cf
        all_markers[stack_moving] = moving_brain_markers_aligned2fixed
        raw_markers[stack_moving] = moving_brain_markers - np.min(moving_brain_markers,axis = 0)
        
        fn = get_stacy_markers_filepath(stack=stack_moving, structure='1')
        if os.path.isfile(fn):
            lesion_brain_markers = bp.unpack_ndarray_file(fn)
            lesion_markers[stack_moving] = lesion_brain_markers - np.min(moving_brain_markers,axis = 0)
            lesion_markers_aligned2fixed = (lesion_brain_markers - om - cm) + t + of + cf
            aligned_lesion_markers[stack_moving] = lesion_markers_aligned2fixed
        else:
            lesion_markers[stack_moving]  = []
            aligned_lesion_markers[stack_moving]  = []
    
    return all_markers,raw_markers, aligned_lesion_markers, lesion_markers


def load_GR_markers(include_stacks={'GR35_L'}, stack_fixed = 'DTA04_L',warp_setting = 24,markertype='All'):
    # Load fixed brains markers.
    
    all_markers = dict()

    for stack_moving in include_stacks:
        if stack_moving == stack_fixed: 
            moving_brain_markers_aligned2fixed = bp.unpack_ndarray_file(get_stacy_markers_filepath(stack=stack_moving, structure='All'))
        else:
            tf_parameter_dict = load_alignment_parameters_v2(stack_f=stack_fixed, stack_m=stack_moving, warp_setting=warp_setting, 
                                                     vol_type_f='annotationAsScore', vol_type_m='annotationAsScore',
                                                     downscale=15)
            cf = tf_parameter_dict['centroid_f']
            cm = tf_parameter_dict['centroid_m']
            of = tf_parameter_dict['crop_origin_f']
            om = tf_parameter_dict['crop_origin_m']
            params = tf_parameter_dict['params']
            Rt = np.reshape(params, (3,4))
            R = Rt[:3,:3]
            t = Rt[:3,3]
        
            moving_brain_markers = bp.unpack_ndarray_file(get_stacy_markers_filepath(stack=stack_moving, structure='All'))
            
            moving_brain_markers_aligned2fixed = np.dot(R, (moving_brain_markers - om - cm).T).T + t + of + cf
        all_markers[stack_moving] = moving_brain_markers_aligned2fixed        
    
    return all_markers


def load_atlas_marker_actors(stack = 'RV14', warp_setting = 24,markertype='All'):
    # Load fixed brains markers.
    
    tf_parameter_dict = load_alignment_parameters_v2(stack_f='Rat_brainstem_atlas', stack_m=stack, warp_setting=warp_setting, 
                                                 vol_type_f='annotationAsScore', vol_type_m='annotationAsScore',
                                                 downscale=15)
    cf = tf_parameter_dict['centroid_f']
    cm = tf_parameter_dict['centroid_m']
    of = tf_parameter_dict['crop_origin_f']
    om = tf_parameter_dict['crop_origin_m']
    params = tf_parameter_dict['params']
    Rt = np.reshape(params, (3,4))
    R = Rt[:3,:3]
    t = Rt[:3,3]
    
    moving_brain_markers = bp.unpack_ndarray_file(get_stacy_markers_filepath(stack=stack, structure='All'))
    moving_brain_markers_aligned2fixed = np.dot(R, (moving_brain_markers - om - cm).T).T + t + of + cf
    
    return moving_brain_markers_aligned2fixed