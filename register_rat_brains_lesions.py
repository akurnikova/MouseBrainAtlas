#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 18 14:22:28 2018

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
from aligner import *
from registration_utilities_v0 import save_alignment_results


from registration_utilities_2 import *


## Set the stack moving, and fixed
stack_fixed = 'RV14_65hrs'
# stack_moving = 'RV14_65hrs'
#stack_moving = 'RV13_64hrs'
stack_moving = 'RFles08'

warp_setting = 25# 24 is affine transformation
# warp_setting = 16
# Adding surr is essential.

def run_aligner(stack_moving,stack_fixed,warp_setting, showplots = 'off', save = 'on'):
    alignment_name_dict=dict(stack_m=stack_moving, 
                                  stack_f=stack_fixed,
                                  warp_setting=warp_setting,
                              vol_type_f='annotationAsScore',
                              vol_type_m='annotationAsScore')
    
    
    aligner_parameters = generate_aligner_parameters(warp_setting=warp_setting, 
                                                     alignment_name_dict=alignment_name_dict)
    
    volume_fixed = aligner_parameters['volume_fixed']
    volume_moving = aligner_parameters['volume_moving']
    
    
    aligner = Aligner4(volume_fixed, volume_moving, labelIndexMap_m2f=aligner_parameters['label_mapping_m2f'])
    
    aligner.set_centroid(centroid_m='structure_centroid', centroid_f='structure_centroid', 
                         indices_m=[aligner_parameters['structure_to_label_moving']['7N_L']])
    
    aligner.set_label_weights(label_weights=aligner_parameters['label_weights_m'])

    
    gradients_f = compute_gradient(volume_fixed, smooth_first=True)
    aligner.load_gradient(gradients=gradients_f) # 120s-170 = 2 mins; all 28, 220s
    
    trial_num = 1
    
    T_all_trials = []
    scores_all_trials = []
    traj_all_trials = []
    
    for _ in range(trial_num):
    
        try:
            T, scores = aligner.optimize(tf_type=aligner_parameters['transform_type'], 
                                         max_iter_num=10,
                                         history_len=20, 
                                         terminate_thresh_rot=.002,
                                         terminate_thresh_trans=.2,
                                         grad_computation_sample_number=aligner_parameters['grad_computation_sample_number'],
                                         lr1=1, lr2=.1,
    #                                     init_T=grid_search_T, 
                                          affine_scaling_limits=(.9, 1.1)
                                        )
            T_all_trials.append(T)
            scores_all_trials.append(scores)
            traj_all_trials.append(aligner.Ts)
            
        except Exception as e:
            sys.stderr.write('%s\n' % e)
    
    ## select the best trial here
    Ts = np.array(aligner.Ts)
    best_trial = np.argsort([np.max(scores) for scores in scores_all_trials])[-1]
    
    # best_trial = 1
    T = T_all_trials[best_trial]
    scores = scores_all_trials[best_trial]
    print 'Best trial:', best_trial
    print max(scores), scores[-1]
    
    crop_origin_m = aligner_parameters['volume_moving_bbox'][[0,2,4]]
    print crop_origin_m
    
    crop_origin_f = aligner_parameters['volume_fixed_bbox'][[0,2,4]]
    print crop_origin_f
    
    if showplots == 'on':

        plt.plot(Ts[:, [0,1,2,4,5,6,8,9,10]]);
        plt.title('rotational params');
        plt.xlabel('Iteration');
        plt.show();
        
        plt.plot(Ts[:, [3,7,11]]);
        plt.title('translation params');
        plt.xlabel('Iteration');
        plt.show();
    
        print T.reshape((3,4))
        plt.figure();
        plt.plot(scores);
        plt.show();
    
    ## save the best trial here
    if save == 'on':
        save_alignment_results(T_all_trials[best_trial], 
                       aligner.centroid_m, aligner.centroid_f, 
                       crop_origin_m, crop_origin_f,
                       score_traj=scores_all_trials[best_trial],
                       parameter_traj=traj_all_trials[best_trial],
                      alignment_name_dict=alignment_name_dict)
    
 #   return aligner,T_all_trials, scores_all_trials,traj_all_trials

run_aligner(stack_moving,stack_fixed,warp_setting, showplots = 'on', save = 'on')

##################

'''
#%%
#Convert parameters
#%%

# If parameters are already saved.
tf_parameter_dict = DataManager.load_alignment_parameters_v2(stack_f=stack_fixed, stack_m=stack_moving, warp_setting=24, 
                                         vol_type_f='annotationAsScore', vol_type_m='annotationAsScore',
                                         downscale=32)
cf = np.array(tf_parameter_dict['centroid_f'])
cm = np.array(tf_parameter_dict['centroid_m'])
of = np.array(tf_parameter_dict['crop_origin_f'])
om = np.array(tf_parameter_dict['crop_origin_m'])
params = np.array(tf_parameter_dict['params'])
Rt = np.reshape(params, (3,4))
R = Rt[:3,:3]
t = Rt[:3,3]

G_movingvol2fixedvol = consolidate(params=params, 
                                   centroid_m=cm, 
                                   centroid_f=cf)

for ind_m, v in volume_moving.iteritems():    
    structure = aligner_parameters['label_to_structure_moving'][ind_m]
    transform_and_save_volume(v, structure, G_movingvol2fixedvol[:3], 
                              crop_origin_f=of,
                              alignment_name_dict=alignment_name_dict)
    
    
#%% DISPLAY

structure = '7n_L'
vol_m = volume_moving[aligner_parameters['structure_to_label_moving'][structure]]

volume_m2fg_in_bboxrel2fixedvol, volume_m2fg_bbox_rel2fixedvol = \
transform_volume_v2(vol=vol_m, tf_params=G_movingvol2fixedvol[:3].flatten())

ydim_f, xdim_f, zdim_f = volume_fixed[aligner_parameters['structure_to_label_fixed'][structure]].shape

v_m2fg = \
crop_and_pad_volume(volume_m2fg_in_bboxrel2fixedvol, in_bbox=volume_m2fg_bbox_rel2fixedvol,
                    out_bbox=(0, xdim_f-1, 0, ydim_f-1, 0, zdim_f-1))

display_volume_sections(v_m2fg, start_level=0)


# Warping all structures.

volume_m2fg = {}
for label_m, vol_m in volume_moving.iteritems():
    
    volume_m2fg_in_bboxrel2fixedvol, volume_m2fg_bbox_rel2fixedvol = \
        transform_volume_v2(vol=vol_m, tf_params=G_movingvol2fixedvol[:3].flatten())

    volume_m2fg[label_m] = crop_and_pad_volume(volume_m2fg_in_bboxrel2fixedvol, in_bbox=volume_m2fg_bbox_rel2fixedvol,
                    out_bbox=(0, xdim_f-1, 0, ydim_f-1, 0, zdim_f-1))
    
    
structures_to_draw = [l for l in volume_m2fg.keys() if not is_surround_label(aligner_parameters['label_to_structure_moving'][l])]

structures_to_draw = [3,4,5,6,9,10,11,12]
colors = {l: name_unsided_to_color[convert_to_original_name(aligner_parameters['label_to_structure_fixed'][l])]
for l in structures_to_draw}

draw_alignment(warped_atlas=volume_m2fg, fixed_volumes=volume_fixed, 
               zs=np.arange(0,volume_fixed.values()[0].shape[2],5), ncols=5,
              structures=structures_to_draw,
              colors=colors,
#               markers=np.array(marker_xyzs_rel2fixedvol.values())
              )
'''