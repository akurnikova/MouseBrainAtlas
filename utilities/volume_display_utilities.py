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

from vis3d_utilities_stacy import *
import cPickle as pickle
import vtk
import scipy.ndimage as snd
import copy
from vis3d_utilities import *

import xml.etree.ElementTree as ET
#folder = '/home/asya/Documents/Yuncong_code/data/meshes_save/'

#%%
def ConnectedVertices(mesh, id_in,id2):
    cellIdList =vtk.vtkIdList()
    mesh.GetPointCells(id_in, cellIdList)
    pointIdList = vtk.vtkIdList()
    mesh.GetCellPoints(cellIdList.GetId(id2), pointIdList)
    return pointIdList.GetId(0), pointIdList.GetId(1)


def polydata_cut(polydata, origin=(0,0,0),\
                   cut_plane_origin = (-60,0,0),cut_plane_normal = (1,0,0)):
    """
    Returns slice of a volume. Option to use different display parameters for planeCutActor2
    Args:
        color (float array): rgb between 0 and 1.
        origin: the initial shift for the mesh.
    """

    if polydata.GetNumberOfPoints() == 0:
        return None

    if origin[0] == 0 and origin[1] == 0 and origin[2] == 0:
        polydata_shifted = polydata
    else:
        polydata_shifted = move_polydata(polydata, origin)
        # Note that move_polydata() discards scalar data stored in polydata.

    ## PlaneCutActor1 will give contours directly through the central plane
    plane=vtk.vtkPlane()
    plane.SetOrigin(cut_plane_origin)
    plane.SetNormal(cut_plane_normal)
    
    cutter=vtk.vtkCutter()
    cutter.SetCutFunction(plane)
    cutter.SetInputData(polydata_shifted)
    cutter.Update()
    
    P = cutter.GetOutput()
    
    contour_points = []
    for i in range(P.GetNumberOfPoints()):
        contour_points.append(P.GetPoint(i))
    contour_points = np.asarray(contour_points)
    
    if P.GetNumberOfPoints() == 0:
        return contour_points,[], P
    order = []
    for k in [0,1]:
        for i in range(P.GetNumberOfPoints()):
            conV = ConnectedVertices(P, i,k)
            if not conV[0] == i and not conV[0] in order:
                order.append([i,conV[0]])
            else:
                order.append([i,conV[1]])
    
    ## Convert to order
    sorted_order = [0]
    r = 0
    start_array = np.asarray(order)[:,0]
    end_array = np.asarray(order)[:,1]
    for i in range(len(order)):
        idx_in_start = np.where(start_array == r)[0]
        idx_in_end = end_array[idx_in_start]
        if idx_in_end[0] not in sorted_order:
            sorted_order.append(idx_in_end[0])
            r = idx_in_end[0]
            continue
        elif idx_in_end[1] not in sorted_order:
            sorted_order.append(idx_in_end[1])
            r = idx_in_end[1]
            continue
        else:
            break
    return contour_points,sorted_order, P
    

def save_RV_meshes(all_meshes, all_origin, folder = '/home/asya/Documents/Yuncong_code/data/meshes_save/'):
    for name_stack, stack in all_meshes.iteritems():
        for name_s, item in stack.iteritems():
            if item == []:
                continue
            writer = vtk.vtkPolyDataWriter()
            writer.SetInputData(item)
            writer.SetFileName(folder+name_stack+'_'+name_s)
            writer.Write()
    pickle.dump(all_origin, open( folder+'origins_save.p', "wb" ) )
    
def load_file_RV_meshes(brains, folder = '/home/asya/Documents/Yuncong_code/data/meshes_save/'):
    A = os.listdir(folder)
    all_meshes = dict.fromkeys(brains)
    for name_stack in brains:
        newdict = dict()
        for name_s in alignment_structures_sided:
                if name_stack+'_'+name_s in A:
                    print name_stack+'_'+name_s
                    reader_poly = vtk.vtkPolyDataReader()
                    reader_poly.SetFileName(folder+name_stack+'_'+name_s)     
                    reader_poly.Update()
                    newdict[name_s] = reader_poly.GetOutput()
        all_meshes[name_stack] = newdict
    all_origin = pickle.load(open(folder+'origins_save.p','rb'))
    return all_meshes, all_origin
    
# Create actors - adjust colors and other rendering options.
def generate_whole_actors(all_meshes,all_origin,structs_to_colors,include_stacks = {'RV16'} ,doBrainstem = 0,wireframe=False,opacity= 0.05):
    brain_actor_list_all_stacks = {}
    for name_stack, brain_meshes in all_meshes.iteritems():
        if name_stack in include_stacks:
            brain_actor_list = {name_s: 
                                       [] if name_s == 'Brainstem' and not doBrainstem
                                       else actor_mesh(mesh, origin=all_origin[name_stack][name_s], 
                                                 wireframe=wireframe,
                                                 color=np.array(name_unsided_to_color[convert_to_original_name(name_s)])/255.,
                                                 #color=(0,0,0),
                                                 opacity= opacity) 
                                  for name_s, mesh in brain_meshes.iteritems()}
            if not doBrainstem and brain_actor_list.has_key('Brainstem'):
                del brain_actor_list['Brainstem']
            brain_actor_list_all_stacks[name_stack] = brain_actor_list
    return brain_actor_list_all_stacks

def generate_slice_actors(all_meshes,all_origin,structs_to_colors,include_stacks = {'RV4_67hrs','RV19_61hrs','RV14_65hrs','RV13_64hrs'},\
                          doBrainstem = 0,wireframe=False,opacity= 0.3,\
                          cut_plane_origin = (20,0,0),cut_plane_normal = (0,0,1),thickness = 10
                          ):
    brain_slice_actor_list_all_stacks = {}
    for name_stack, brain_meshes in all_meshes.iteritems():
        if name_stack in include_stacks:
            brain_actor_list_slice = {name_s: [] if mesh == [] else actor_mesh_cut(mesh, color=(0,0,0),
                                                                   origin=all_origin[name_stack][name_s], 
                                                                   cut_plane_origin = cut_plane_origin, 
                                                                   cut_plane_normal = cut_plane_normal,
                                                                   thickness = thickness,\
                                                                   opacity = opacity,opacityE = 0.5)
                                      for name_s, mesh in brain_meshes.iteritems()}
            
            brain_slice_actor_list_all_stacks[name_stack] = brain_actor_list_slice
    return brain_slice_actor_list_all_stacks

def generate_density_slice_actors(cut_plane_origin = (0,0,130),cut_plane_normal = (1,0,0), thickness = 2,\
                                opacity= 0.3, color = (1,0,0),linewidth = 5, include_stacks = {'RV14'},kernel_bandwidth = 0.133,include_d = {'0.5','0.7'}):
    max_d = float(max(include_d))
    density_slice_actor_list_all_stacks = {}
    for stack in include_stacks:
        density_actor_list_slice = {}
        for d in include_d:
            fp = DataManager.get_density_pvol_filepath(stack,kernel_bandwidth)
            P = pickle.load(open(fp,'rb'))
            
            density_vol = (P['vol']/np.max(P['vol'])) > d
            imagedata = volume_to_imagedata_respace(np.float32(density_vol), origin=P['origin'], spacing = [P['xstep'],P['ystep'],P['zstep']])
            surface = vtk.vtkMarchingCubes()
            surface.SetInputData(imagedata)
            surface.ComputeNormalsOn();
            surface.SetValue(0, 0.5);
            surface.Update()
            
            mesh = surface.GetOutput()
            #fp = DataManager.get_density_mesh_filepath(stack, d)
            #save_mesh_stl(mesh,fp)
            '''
            maxcol = np.max(density['vol'])
            V = density['vol']/np.sum(density['vol'])
            sorted_data = np.sort(np.ravel(V))[::-1]
            percentile_d = sorted_data[np.where(np.cumsum(sorted_data)>float(d))[0][0]]
            density_vol = V>percentile_d
            imagedata = volume_to_imagedata_respace(np.float32(density_vol), origin=[xmin, ymin, zmin], spacing = [xstep,ystep,zstep])
            surface = vtk.vtkMarchingCubes()
            surface.SetInputData(imagedata)
            surface.ComputeNormalsOn();
            surface.SetValue(0, 0.5);
            surface.Update()
        
            mesh = surface.GetOutput()
            #fp = DataManager.get_density_mesh_filepath_v2(stack, d)
            #save_mesh_stl(mesh,fp)
            '''
            
            color_scale = float(d)/max_d
            density_actor_list_slice[d] = actor_mesh_cut(mesh, color=np.asarray(color)*color_scale,
                                                                   origin=(0,0,0), 
                                                                   cut_plane_origin = cut_plane_origin, 
                                                                   cut_plane_normal = cut_plane_normal,
                                                                   thickness = thickness,\
                                                                   opacity = opacity, opacityE = 0.5,\
                                                                   linewidth = linewidth,dotted_line = 0)
            
        
        density_slice_actor_list_all_stacks[stack] = density_actor_list_slice
    
    return density_slice_actor_list_all_stacks

def get_density_slice(cut_plane_origin = (0,0,130),cut_plane_normal = (1,0,0), include_stacks = {'RV14'},include_d = {'0.5','0.7'}):
    max_d = float(max(include_d))
    density_slice_list_all_stacks = {}
    for name_stack in include_stacks:
        density_list_slice = {}
        for d in include_d:
            fp = DataManager.get_density_mesh_filepath(name_stack,d)
            v,f = load_mesh_stl(fp)
            mesh = mesh_to_polydata(v,f)
            color_scale = float(d)/max_d
            density_list_slice[d] = polydata_cut(mesh, origin=(0,0,0),\
                              cut_plane_origin = cut_plane_origin,cut_plane_normal = cut_plane_normal)    
            
        
        density_slice_list_all_stacks[name_stack] = density_list_slice
    
    return density_slice_list_all_stacks

def generate_slice_actors_atlas(cut_plane_origin = (20,0,0),cut_plane_normal = (0,0,1),thickness = 10,\
                                doBrainstem = 0,doCortex = 0,opacity= 0.3, structs_to_plot = all_known_structures_sided, constcolor = [],spacefill = False,atlas_name = 'Rat_atlas'):
    
    struct_pos = load_pickle(DataManager.get_structure_mean_positions_filepath(atlas_name=atlas_name))
    
    slice_actor_dict_atlas = {}
   
    for name_s in all_known_structures_sided:
        if not name_s in struct_pos: continue
        if name_s == 'Brainstem' and doBrainstem == 0: continue
        if name_s == 'Cortex' and doCortex == 0: continue
        if name_s == 'hipp_L' and doCortex == 0: continue
        if name_s == 'hipp_R' and doCortex == 0: continue
        
        print name_s
        
        fp = DataManager.get_mesh_filepath(stack_m=atlas_name, structure=name_s)
        fp_spacefill = DataManager.get_mesh_filepath_spacefill(stack_m=atlas_name, structure=name_s)

        if os.path.isfile(fp_spacefill) and spacefill:
            mesh_rel2canon = load_mesh_stl(fp_spacefill,return_polydata_only=True)      
        else:
            mesh_rel2canon = load_mesh_stl(fp,return_polydata_only=True)
        
        if constcolor == []:
            color_plot =  np.array(name_unsided_to_color[convert_to_original_name(name_s)])/255.
        else:
            color_plot = constcolor
            
        if name_s in fiber_tracts_sided:
            slice_actor_dict_atlas[name_s] = actor_mesh_cut(mesh_rel2canon,
                               color=color_plot,
                               origin=(0,0,0), 
                               cut_plane_origin = cut_plane_origin, 
                               cut_plane_normal = cut_plane_normal,
                               thickness = thickness,\
                               opacity = 0.3,opacityE = 0.5,
                               linewidth = 2,dotted_line = 0)
        if name_s in solid_structures_sided:
            slice_actor_dict_atlas[name_s] = actor_mesh_cut(mesh_rel2canon,
                               color=color_plot,
                               origin=(0,0,0), 
                               cut_plane_origin = cut_plane_origin, 
                               cut_plane_normal = cut_plane_normal,
                               thickness = thickness,\
                               opacity = 0.3,opacityE = 0.1,
                               linewidth = 2,dotted_line = 0)
        if name_s in dotted_structures_sided:
            slice_actor_dict_atlas[name_s] = actor_mesh_cut(mesh_rel2canon,
                               color=color_plot,
                               origin=(0,0,0), 
                               cut_plane_origin = cut_plane_origin, 
                               cut_plane_normal = cut_plane_normal,
                               thickness = thickness,\
                               opacity = 0.3,opacityE = 0.5,
                               linewidth = 2,dotted_line = 1)
        
    return slice_actor_dict_atlas



def generate_slice_actors_atlas_brainstem(cut_plane_origin = (20,0,0),cut_plane_normal = (0,0,1),thickness = 10,\
                                opacity= 0.3, constcolor = [],spacefill = False,atlas_name = 'Rat_brainstem_atlas'):
    
    struct_pos = load_pickle(DataManager.get_structure_mean_positions_filepath(atlas_name=atlas_name))
    
    slice_actor_dict_atlas = {}
   
    for name_s in ['5N_L', '5N_R', '7N_L','7N_R','7n_L',
                   '7n_R',  'Amb_L', 'Amb_R', 'IO_L', 'IO_R',
                   'LRT_L', 'LRT_R', 'Pr5_L', 'Pr5_R', 'RtTg_L',
                   'RtTg_R', 'SPVmu_L', 'SPVmu_R', 'SpVI_L', 'SpVI_R',
                   'SpVC_L', 'SpVC_R', 'SpVO_L', 'SpVO_R',
                   'Vestibular_L', 'Vestibular_R',
                   'Brainstem']:
        
        
        fp = DataManager.get_mesh_filepath(stack_m=atlas_name, structure=name_s)
        fp_spacefill = DataManager.get_mesh_filepath_spacefill(stack_m=atlas_name, structure=name_s)

        if os.path.isfile(fp_spacefill) and spacefill:
            mesh_rel2canon = load_mesh_stl(fp_spacefill,return_polydata_only=True)      
        else:
            mesh_rel2canon = load_mesh_stl(fp,return_polydata_only=True)
        
        if constcolor == []:
            color_plot =  np.array(name_unsided_to_color[convert_to_original_name(name_s)])/255.
        else:
            color_plot = constcolor
            
        slice_actor_dict_atlas[name_s] = actor_mesh_cut(mesh_rel2canon,
                               color=color_plot,
                               origin=(0,0,0), 
                               cut_plane_origin = cut_plane_origin, 
                               cut_plane_normal = cut_plane_normal,
                               thickness = thickness,\
                               opacity = 0.3,opacityE = 0.1,
                               linewidth = 2,dotted_line = 0,
                               lesion_actor = True)
        
    return slice_actor_dict_atlas


def generate_slice_actors_atlas_ms(all_meshes,all_origin, cut_plane_origin = (20,0,0),cut_plane_normal = (0,0,1),thickness = 10,\
                                opacity= 0.3, constcolor = [],spacefill = False):
        
    slice_actor_dict_atlas = {}
   
    for name_s in ['7N_L','7n_L', 'Amb_L',  'IO_L', 'LRT_L', 'Brainstem']:
        
        mesh = all_meshes[name_s]
        
        if constcolor == []:
            color_plot =  np.array(name_unsided_to_color[convert_to_original_name(name_s)])/255.
        else:
            color_plot = constcolor
            
        slice_actor_dict_atlas[name_s] = actor_mesh_cut(mesh,
                               color=color_plot,
                               origin=all_origin[name_s], 
                               cut_plane_origin = cut_plane_origin, 
                               cut_plane_normal = cut_plane_normal,
                               thickness = thickness,\
                               opacity = 0.3,opacityE = 0.1,
                               linewidth = 2,dotted_line = 0,
                               lesion_actor = True)
        
    return slice_actor_dict_atlas


def generate_slice_actors_atlas_ms_no_origin(all_meshes, cut_plane_origin = (20,0,0),cut_plane_normal = (0,0,1),thickness = 10,\
                                opacity= 0.3, constcolor = [],spacefill = False):
        
    slice_actor_dict_atlas = {}
   
    for name_s in ['7N_L','7n_L', 'Amb_L',  'IO_L', 'LRT_L']:
        
        mesh = all_meshes[name_s]
        
        if constcolor == []:
            color_plot =  np.array(name_unsided_to_color[convert_to_original_name(name_s)])/255.
        else:
            color_plot = constcolor
            
        slice_actor_dict_atlas[name_s] = actor_mesh_cut(mesh,
                               color=color_plot,
                               origin=(0,0,0),
                               cut_plane_origin = cut_plane_origin, 
                               cut_plane_normal = cut_plane_normal,
                               thickness = thickness,\
                               opacity = 0.3,opacityE = 0.1,
                               linewidth = 2,dotted_line = 0,
                               lesion_actor = True)
        
    return slice_actor_dict_atlas

def generate_whole_actors_lesions(all_meshes,all_origin,stacks_to_colors,include_stacks = {'RFles02','RFles03','RFles04','RFles05'}, wireframe=False,opacity= 0.1):
    brain_actor_list_all_stacks = {}
    for name_stack, brain_meshes in all_meshes.iteritems():
        if name_stack in include_stacks:
            brain_actor_list = {'Lesion': actor_mesh(brain_meshes, 
                                                     origin=all_origin[name_stack], 
                                                     wireframe=wireframe,
                                                     color=stacks_to_colors[name_stack],
                                                     opacity= opacity)}
            brain_actor_list_all_stacks[name_stack] = brain_actor_list
    return brain_actor_list_all_stacks

def generate_slice_actors_lesions(all_meshes,all_origin,stacks_to_colors,include_stacks = {'RFles02','RFles03','RFles04','RFles05'}, wireframe=False,opacity= 0.3,\
                                  cut_plane_origin = (0,0,10),cut_plane_normal = (0,0,1),thickness = 10):
    slice_actor_list_all_stacks = {}
    for name_stack, brain_meshes in all_meshes.iteritems():
        if name_stack in include_stacks:
            A,B,C,D = actor_mesh_cut(all_meshes[name_stack], color=stacks_to_colors[name_stack],\
                                         wireframe=False,\
                                         origin=all_origin[name_stack],
                                         cut_plane_origin = cut_plane_origin,cut_plane_normal = cut_plane_normal,thickness = thickness,\
                                     opacity = opacity,opacityE = opacity,lesion_actor =True)
            
            slice_actor_list_all_stacks[name_stack] = [A,B,C,D]
    return slice_actor_list_all_stacks

def generate_RV_marker_actors(stacks_to_colors, include_stacks={'RV4','RV19'}, stack_fixed = 'Rat_atlas',warp_setting=24):
    # Load fixed brains markers.
    fixed_brain_markers = bp.unpack_ndarray_file(get_stacy_markers_filepath(stack=stack_fixed, structure='All'))
    
    
    fixed_brain_marker_actors = [actor_sphere((x,y,z), radius=2, color=stacks_to_colors[stack_fixed]) 
                                 for x,y,z in fixed_brain_markers]
    
    
    if stack_fixed not in include_stacks:
        fixed_brain_marker_actors = []
    
    moving_brain_marker_aligned2fixed_actors_all_stacks = {}
    for stack_moving in include_stacks:
        if stack_moving == stack_fixed: continue
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
        moving_brain_marker_aligned2fixed_actors = [actor_sphere((x,y,z), radius=2, color=stacks_to_colors[stack_moving]) 
                                                    for x,y,z in moving_brain_markers_aligned2fixed]
        
        moving_brain_marker_aligned2fixed_actors_all_stacks[stack_moving] = moving_brain_marker_aligned2fixed_actors
    
    return fixed_brain_marker_actors, moving_brain_marker_aligned2fixed_actors_all_stacks

def generate_RV_marker_actors_for_atlas(include_stacks={'RV4','RV19'}, radius=2, stacks_to_colors = []):
    # Load fixed brains markers.
    
    if stacks_to_colors == []:
        col = (1,0,0)
    else:
        col = stacks_to_colors[stack_moving]
        
    marker_actors_all_stacks = list()
    for stack_moving in include_stacks:
        tf_parameter_dict = load_alignment_parameters_v2(stack_f='Rat_brainstem_atlas', stack_m=stack_moving, warp_setting=24, 
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
        moving_brain_markers_aligned2atlas = np.dot(R, (moving_brain_markers - om - cm).T).T + t + of + cf
        moving_brain_markers_aligned2atlas_actors = [actor_sphere((x,y,z), radius = radius, color = col) 
                                                    for x,y,z in moving_brain_markers_aligned2atlas]
        
        marker_actors_all_stacks.extend(a for a in moving_brain_markers_aligned2atlas_actors)
    
    return marker_actors_all_stacks

def generate_RV_marker_actors_slice(cut_plane_normal, cut_plane_origin, stacks_to_colors,slice_thickness = 10, include_stacks={'RV4_67hrs','RV19_61hrs','RV14_65hrs','RV13_64hrs'}, stack_fixed = 'RV4_67hrs',warp_setting=24):
    # Load fixed brains markers.
    plane=vtk.vtkPlane()
    fixed_brain_markers = bp.unpack_ndarray_file(get_stacy_markers_filepath(stack=stack_fixed, structure='All'))
    
    all_brain_marker_actors_slice = []
    
    if stack_fixed in include_stacks:
        for x,y,z in fixed_brain_markers:
            dist = np.abs(plane.DistanceToPlane((x,y,z),cut_plane_normal,cut_plane_origin))
            if dist < slice_thickness:
                all_brain_marker_actors_slice.append(actor_sphere((x,y,z), radius=2, color=stacks_to_colors[stack_fixed]))
    else:
        fixed_brain_marker_actors = []
    
    
    for stack_moving in include_stacks:
        if stack_moving == stack_fixed: continue
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
    
        for x,y,z in moving_brain_markers_aligned2fixed:
            dist = np.abs(plane.DistanceToPlane((x,y,z),cut_plane_normal,cut_plane_origin))
            if dist < slice_thickness:
                all_brain_marker_actors_slice.append(actor_sphere((x,y,z), radius=2, color=stacks_to_colors[stack_moving]))
    return all_brain_marker_actors_slice


def generate_RV_marker_actors_slice_for_atlas(cut_plane_normal, cut_plane_origin, slice_thickness = 10, 
                                              include_stacks={'RV4','RV19'}, radius=2, stacks_to_colors = []):
    # Load fixed brains markers.
    plane=vtk.vtkPlane()
    if stacks_to_colors == []:
        col = (1,0,0)
    else:
        col = stacks_to_colors[stack_moving]
    all_brain_marker_actors_slice = list()
    for stack_moving in include_stacks:
        tf_parameter_dict = load_alignment_parameters_v2(stack_f='Rat_atlas', stack_m=stack_moving, warp_setting=24, 
                                                 vol_type_f='annotationAsScore', vol_type_m='annotationAsScore', downscale=15)
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
    
        for x,y,z in moving_brain_markers_aligned2fixed:
            dist = np.abs(plane.DistanceToPlane((x,y,z),cut_plane_normal,cut_plane_origin))
            if dist < slice_thickness:
                all_brain_marker_actors_slice.append(actor_sphere((x,y,z), radius=radius, color=col))
    return all_brain_marker_actors_slice


def generate_RV_meshes(stack_moving_list, stack_fixed = 'RV16',warp_setting=24):
# Convert volumes to meshes, for ALL BRAINS
    new_stack_list = copy.deepcopy(stack_moving_list)
    all_meshes = dict()
    all_origin = dict()

    new_stack_list.discard(stack_fixed)
    brain_fixed = load_original_volume_all_known_structures_v2(stack=stack_fixed, 
                                                                sided=True, 
                                                               include_surround=False, 
                                                               common_shape=False)
    
    #fixed_brain_meshes = {name_s: volume_to_polydata_brainstem(brain_fixed[name_s][0], 
    #                                                 num_simplify_iter= 1, smooth=True) 
    #                            if name_s == 'Brainstem' else
    #                              volume_to_polydata(brain_fixed[name_s][0], 
    #                                                 num_simplify_iter= 3, smooth=True) 
    #                      for name_s in brain_fixed.iterkeys()}
    fixed_brain_meshes = dict()
    for name_s in brain_fixed.iterkeys():
        if name_s in alignment_structures_sided:
                fixed_brain_meshes[name_s] = volume_to_polydata(brain_fixed[name_s][0], 
                                                     num_simplify_iter= 3, smooth=True) 
                                
    
    all_meshes[stack_fixed] = fixed_brain_meshes
    all_origin[stack_fixed] =  {name_s: brain_fixed[name_s][1][[0,2,4]] for name_s in brain_fixed.iterkeys()}
    
    for stack_moving in new_stack_list:
        brain_moving = \
        load_original_volume_all_known_structures_v2(stack=stack_moving, sided=True, 
                                                              include_surround=False,
                                                                 common_shape=False)
        
        tf_parameter_dict = load_alignment_parameters_v2(stack_f=stack_fixed, stack_m=stack_moving, warp_setting=warp_setting, 
                                             vol_type_f='annotationAsScore', vol_type_m='annotationAsScore',
                                             downscale=15)
        cf = np.array(tf_parameter_dict['centroid_f'])
        cm = np.array(tf_parameter_dict['centroid_m'])
        of = np.array(tf_parameter_dict['crop_origin_f'])
        om = np.array(tf_parameter_dict['crop_origin_m'])
        params = np.array(tf_parameter_dict['params'])
    
        G_movingvol2fixedvol = consolidate(params=params, 
                                           centroid_m=cm, 
                                           centroid_f=cf)
        
        brain_movingAlignedToFixed = {}
        for name_s, (v, b) in brain_moving.iteritems():
                
            volume_m_warped_inbbox, volume_m_warped_bbox_rel2movingvol = \
            transform_volume_v3(vol=v, bbox=b-om[[0,0,1,1,2,2]], tf_params=G_movingvol2fixedvol[:3].flatten())
    
            volume_m_warped_bbox_rel2fixedvol = volume_m_warped_bbox_rel2movingvol
            volume_m_warped_bbox = volume_m_warped_bbox_rel2fixedvol + of[[0,0,1,1,2,2]]
        
            ## Fix holes due to rotation, added by stacy
            volume_m_warped_inbbox_fix = volume_m_warped_inbbox
            for i in range(0,np.shape(volume_m_warped_inbbox)[2]):
                volume_m_warped_inbbox_fix[:,:,i] = snd.filters.gaussian_filter(volume_m_warped_inbbox[:,:,i], 1.5, order=0)>0
        
            brain_movingAlignedToFixed[name_s] = (volume_m_warped_inbbox_fix, volume_m_warped_bbox)
            
            #moving_brain_meshes = {name_s:  volume_to_polydata_brainstem(brain_movingAlignedToFixed[name_s][0], 
            #                                         num_simplify_iter= 1, smooth=True) 
            #                                if name_s == 'Brainstem' else
            #                                volume_to_polydata(brain_movingAlignedToFixed[name_s][0], 
            #                                         num_simplify_iter= 3, smooth=True)
            #                   for name_s in brain_movingAlignedToFixed.iterkeys()}
            moving_brain_meshes = {name_s:  [] 
                                    if name_s not in alignment_structures_sided else
                                            volume_to_polydata(brain_movingAlignedToFixed[name_s][0], 
                                                     num_simplify_iter= 3, smooth=True)
                               for name_s in brain_movingAlignedToFixed.iterkeys()}
        bboxes = {name_s: brain_movingAlignedToFixed[name_s][1][[0,2,4]] for name_s in brain_movingAlignedToFixed.iterkeys()}
        all_meshes[stack_moving] = moving_brain_meshes
        all_origin[stack_moving] = bboxes
    return all_meshes, all_origin


def generate_and_save_RV_meshes(stack_moving_list, stack_fixed = 'Rat_brainstem_atlas',warp_setting=24, folder = '/home/asya/Documents/Yuncong_code/data/meshes_save/'):
    new_stack_list = copy.deepcopy(stack_moving_list)
    all_meshes = dict()
    all_origin = dict()
    
    structures_to_display =alignment_structures_sided
    
 #   for struct in structures_to_display:
 #       if struct[-1] == 'R':
 #           structures_to_display.remove(struct)
    
    structures_to_display.append('PreBotC_L')
    #structures_to_display.append('Lesion')
    
    new_stack_list.discard(stack_fixed)
    brain_fixed = load_original_volume_all_known_structures_v2(stack=stack_fixed, 
                                                                sided=True, 
                                                               include_surround=False, 
                                                               common_shape=False)
    
    if os.path.exists(folder+'origins_save.p'):
        all_origin = pickle.load(open(folder+'origins_save.p','rb'))
        if not all_origin.has_key(stack_fixed):
            all_origin[stack_fixed] = dict()
    else:
        all_origin = dict()
        all_origin[stack_fixed] = dict()

    
    for name_s in brain_fixed.iterkeys():
        
        all_origin[stack_fixed][name_s] = brain_fixed[name_s][1][[0,2,4]]

        if name_s not in structures_to_display:
            continue
        
        ## first check if file already exists
        filename = folder+stack_fixed+'_'+name_s
  #      if os.path.exists(filename):
  #          continue
            
        item = volume_to_polydata(brain_fixed[name_s][0],
                                                     num_simplify_iter= 3, smooth=True) 

        writer = vtk.vtkPolyDataWriter()
        writer.SetInputData(item)
        writer.SetFileName(filename)
        writer.Write()
                
    
    for stack_moving in new_stack_list:
        brain_moving = \
        load_original_volume_all_known_structures_v2(stack=stack_moving, sided=True, 
                                                              include_surround=False,
                                                                 common_shape=False)
        vol_f_type = 'annotationAsScore'
        tf_parameter_dict = load_alignment_parameters_v2(stack_f=stack_fixed, stack_m=stack_moving, warp_setting=warp_setting, 
                                             vol_type_f=vol_f_type, vol_type_m='annotationAsScore',
                                             downscale=15)
        cf = np.array(tf_parameter_dict['centroid_f'])
        cm = np.array(tf_parameter_dict['centroid_m'])
        of = np.array(tf_parameter_dict['crop_origin_f'])
        om = np.array(tf_parameter_dict['crop_origin_m'])
        params = np.array(tf_parameter_dict['params'])
    
#        params[0]  = 
     #   params[5]  = 0.7
     #   params[10] = 0.6
 #       cf = np.array([158.5, 113.1, 59.8])
 #       cf = np.array([158.5, 120, 60])
     #   of= [-380,   45,  -84]
        
        G_movingvol2fixedvol = consolidate(params=params, 
                                           centroid_m=cm, 
                                           centroid_f=cf)
        
        bboxes = dict()
        for name_s, (v, b) in brain_moving.iteritems():
                
            if name_s not in structures_to_display:
                continue
        
            volume_m_warped_inbbox, volume_m_warped_bbox_rel2movingvol = \
            transform_volume_v3(vol=v, bbox=b-om[[0,0,1,1,2,2]], tf_params=G_movingvol2fixedvol[:3].flatten())
    
            volume_m_warped_bbox_rel2fixedvol = volume_m_warped_bbox_rel2movingvol
            volume_m_warped_bbox = volume_m_warped_bbox_rel2fixedvol + of[[0,0,1,1,2,2]]
        
            ## Fix holes due to rotation, added by stacy
            volume_m_warped_inbbox_fix = volume_m_warped_inbbox
            for i in range(0,np.shape(volume_m_warped_inbbox)[2]):
                volume_m_warped_inbbox_fix[:,:,i] = snd.filters.gaussian_filter(volume_m_warped_inbbox[:,:,i], 1.5, order=0)>0
                               
            bboxes[name_s] = volume_m_warped_bbox[[0,2,4]]
            
             ## check if file already exists
            filename = folder+stack_moving+'_'+name_s
       #     if os.path.exists(filename):
       #         continue
                
            item = volume_to_polydata(volume_m_warped_inbbox_fix, 
                                                     num_simplify_iter= 3, smooth=True)
        
            writer = vtk.vtkPolyDataWriter()
            writer.SetInputData(item)
            writer.SetFileName(filename)
            writer.Write()
                        
        all_origin[stack_moving] = bboxes

    pickle.dump(all_origin, open(folder+'origins_save.p', "wb" ) )
    

def generate_lesion_meshes(stack_moving_list, stack_fixed = 'Rat_brainstem_atlas',warp_setting=24):
    all_meshes = dict()
    all_origin = dict()
    all_centroid = dict()
    all_volumes= dict()
              
    for stack_moving in stack_moving_list:
        brain_moving = \
        load_original_volume_all_known_structures_v2(stack=stack_moving, sided=True, 
                                                              include_surround=False,
                                                                 common_shape=False)
        
        if stack_moving == stack_fixed:
            name_s = 'Lesion'
            v, b = brain_moving['Lesion']
            polydata_1 = volume_to_polydata(v, num_simplify_iter= 3, smooth=True)
        
            instance_centroid = np.mean(np.nonzero(v), axis=1)[[1,0,2]] + b[[0,2,4]]
         
                            
            all_origin[stack_moving] = b[[0,2,4]]
            all_meshes[stack_moving] = polydata_1
            all_centroid[stack_moving]  = instance_centroid
            all_volumes[stack_moving] = v
            continue
        
        vol_f_type = 'annotationAsScore'
        tf_parameter_dict = load_alignment_parameters_v2(stack_f=stack_fixed, stack_m=stack_moving, warp_setting=warp_setting, 
                                             vol_type_f=vol_f_type, vol_type_m='annotationAsScore',
                                             downscale=15)
        cf = np.array(tf_parameter_dict['centroid_f'])
        cm = np.array(tf_parameter_dict['centroid_m'])
        of = np.array(tf_parameter_dict['crop_origin_f'])
        om = np.array(tf_parameter_dict['crop_origin_m'])
        params = np.array(tf_parameter_dict['params'])
    
        G_movingvol2fixedvol = consolidate(params=params, 
                                           centroid_m=cm, 
                                           centroid_f=cf)
        
        
        name_s = 'Lesion'
        print stack_moving
        v, b = brain_moving['Lesion']
        
        volume_m_warped_inbbox, volume_m_warped_bbox_rel2movingvol = \
            transform_volume_v3(vol=v, bbox=b-om[[0,0,1,1,2,2]], tf_params=G_movingvol2fixedvol[:3].flatten())
    
        volume_m_warped_bbox_rel2fixedvol = volume_m_warped_bbox_rel2movingvol
        volume_m_warped_bbox = volume_m_warped_bbox_rel2fixedvol + of[[0,0,1,1,2,2]]
        
        ## Fix holes due to rotation, added by stacy
        volume_m_warped_inbbox_fix = volume_m_warped_inbbox
        for i in range(0,np.shape(volume_m_warped_inbbox)[2]):
            volume_m_warped_inbbox_fix[:,:,i] = snd.filters.gaussian_filter(volume_m_warped_inbbox[:,:,i], 1.5, order=0)>0
                               
        bbox_lesion = volume_m_warped_bbox[[0,2,4]]
            
        
        polydata_1 = volume_to_polydata(volume_m_warped_inbbox_fix, 
                                                     num_simplify_iter= 3, smooth=True)
        
        instance_centroid = np.mean(np.nonzero(volume_m_warped_inbbox_fix), axis=1)[[1,0,2]] + bbox_lesion[[0,1,2]]
     
                        
        all_origin[stack_moving] = bbox_lesion
        all_meshes[stack_moving] = polydata_1
        all_centroid[stack_moving]  = instance_centroid
        all_volumes[stack_moving] = volume_m_warped_inbbox_fix

    return all_meshes,all_origin,all_centroid,all_volumes


def generate_ms_meshes():
# Convert volumes to meshes, for 'DTA04_L'
    structs_to_plot_sided = ['7N_L', '7n_L', 'Amb_L', 'LRT_L', 'IO_L']

    all_meshes = dict()
    all_origin = dict()

    brain_fixed = load_original_volume_all_known_structures_v2(stack='DTA04_L', 
                                                                sided=True, 
                                                               include_surround=False, 
                                                               common_shape=False)
    fixed_brain_meshes = dict()
    for name_s in brain_fixed.iterkeys():
        if name_s in structs_to_plot_sided:
                all_meshes[name_s] = volume_to_polydata(brain_fixed[name_s][0], brain_fixed[name_s][1][[0,2,4]],
                                                     num_simplify_iter= 3, smooth=True)
                all_origin[name_s] = brain_fixed[name_s][1][[0,2,4]]
                                
    
    ## GET BRAINSTEM FROM DTA06
    '''
    brain_moving = load_original_volume_all_known_structures_v2(stack='DTA06_L', sided=True, 
                                                              include_surround=False,
                                                                 common_shape=False)
        
    tf_parameter_dict = load_alignment_parameters_v2(stack_f='DTA04_L', stack_m='DTA06_L', warp_setting=24, 
                                             vol_type_f='annotationAsScore', vol_type_m='annotationAsScore',
                                             downscale=15)
    cf = np.array(tf_parameter_dict['centroid_f'])
    cm = np.array(tf_parameter_dict['centroid_m'])
    of = np.array(tf_parameter_dict['crop_origin_f'])
    om = np.array(tf_parameter_dict['crop_origin_m'])
    params = np.array(tf_parameter_dict['params'])
    
    G_movingvol2fixedvol = consolidate(params=params, 
                                           centroid_m=cm, 
                                           centroid_f=cf)
        
    v,b = brain_moving['Brainstem']
                
    volume_m_warped_inbbox, volume_m_warped_bbox_rel2movingvol = \
            transform_volume_v3(vol=v, bbox=b-om[[0,0,1,1,2,2]], tf_params=G_movingvol2fixedvol[:3].flatten())
    
    volume_m_warped_bbox_rel2fixedvol = volume_m_warped_bbox_rel2movingvol
    volume_m_warped_bbox = volume_m_warped_bbox_rel2fixedvol + of[[0,0,1,1,2,2]]
        
    ## Fix holes due to rotation, added by stacy
    volume_m_warped_inbbox_fix = volume_m_warped_inbbox
    for i in range(0,np.shape(volume_m_warped_inbbox)[2]):
        volume_m_warped_inbbox_fix[:,:,i] = snd.filters.gaussian_filter(volume_m_warped_inbbox[:,:,i], 1.5, order=0)>0
        
    v_new,b_new = (volume_m_warped_inbbox, volume_m_warped_bbox)
            
    all_meshes['Brainstem'] = volume_to_polydata(v_new,#b_new[[0,2,4]], 
                                                     num_simplify_iter= 3, smooth=True) 
    
    
    all_origin['Brainstem'] = b_new[[0,2,4]]
    '''
    
    return all_meshes, all_origin


def generate_ms_meshes_0origin(stack_fixed):
# Convert volumes to meshes, for 'DTA04_L'
    structs_to_plot_sided = ['7N_L', '7n_L', 'Amb_L', 'LRT_L', 'IO_L']

    all_meshes = dict()

    brain_fixed = load_original_volume_all_known_structures_v2(stack=stack_fixed, 
                                                                sided=True, 
                                                               include_surround=False, 
                                                               common_shape=False)
    fixed_brain_meshes = dict()
    for name_s in brain_fixed.iterkeys():
        if name_s in structs_to_plot_sided:
                all_meshes[name_s] = volume_to_polydata(brain_fixed[name_s][0], brain_fixed[name_s][1][[0,2,4]],
                                                     num_simplify_iter= 3, smooth=True)
                                             
    '''
    ## GET BRAINSTEM FROM DTA06
    brain_moving = load_original_volume_all_known_structures_v2(stack='DTA06_L', sided=True, 
                                                              include_surround=False,
                                                                 common_shape=False)
        
    tf_parameter_dict = load_alignment_parameters_v2(stack_f='GR35_L', stack_m='DTA06_L', warp_setting=24, 
                                             vol_type_f='annotationAsScore', vol_type_m='annotationAsScore',
                                             downscale=15)
    cf = np.array(tf_parameter_dict['centroid_f'])
    cm = np.array(tf_parameter_dict['centroid_m'])
    of = np.array(tf_parameter_dict['crop_origin_f'])
    om = np.array(tf_parameter_dict['crop_origin_m'])
    params = np.array(tf_parameter_dict['params'])
    
    G_movingvol2fixedvol = consolidate(params=params, 
                                           centroid_m=cm, 
                                           centroid_f=cf)
        
    v,b = brain_moving['Brainstem']
                
    volume_m_warped_inbbox, volume_m_warped_bbox_rel2movingvol = \
            transform_volume_v3(vol=v, bbox=b-om[[0,0,1,1,2,2]], tf_params=G_movingvol2fixedvol[:3].flatten())
    
    volume_m_warped_bbox_rel2fixedvol = volume_m_warped_bbox_rel2movingvol
    volume_m_warped_bbox = volume_m_warped_bbox_rel2fixedvol + of[[0,0,1,1,2,2]]
        
    ## Fix holes due to rotation, added by stacy
    volume_m_warped_inbbox_fix = volume_m_warped_inbbox
    for i in range(0,np.shape(volume_m_warped_inbbox)[2]):
        volume_m_warped_inbbox_fix[:,:,i] = snd.filters.gaussian_filter(volume_m_warped_inbbox[:,:,i], 1.5, order=0)>0
        
    v_new,b_new = (volume_m_warped_inbbox, volume_m_warped_bbox)
            
    all_meshes['Brainstem'] = volume_to_polydata(v_new,b_new[[0,2,4]], 
                                                     num_simplify_iter= 3, smooth=True) 
    '''
    return all_meshes

def volume_type_to_str(t):
    if t == 'annotationAsScore':
        return 'annotationAsScoreVolume'
    if t == 'score':
        return 'scoreVolume'
    else:
        raise Exception('Volume type %s is not recognized.' % t)
        
def get_original_volume_basename(stack, structure=None, volume_type='score', **kwargs):

    downscale = 15
    
    components = []
    if downscale is not None:
        components.append('down%(downscale)d' % {'downscale':downscale})
    tmp_str = '_'.join(components)

    basename = '%(stack)s_%(tmp_str)s_%(volstr)s' % \
        {'stack':stack, 'tmp_str':tmp_str, 'volstr':volume_type_to_str(volume_type)}
    if structure is not None:
        basename += '_' + structure
    return basename

def get_warped_volume_basename(stack_m,
                               stack_f=None,
                               warp_setting=None,
                               structure_m=None,
                               structure_f=None,
                               vol_type_m='score',
                               vol_type_f='score',
                               **kwargs):

    basename_m = get_original_volume_basename(stack=stack_m, volume_type=vol_type_m, structure=structure_m)

    if stack_f is None:
        assert warp_setting is None
        vol_name = basename_m
    else:
        basename_f = get_original_volume_basename(stack=stack_f, volume_type=vol_type_f, structure=structure_f)
        vol_name = basename_m + '_warp%(warp)d_' % {'warp':warp_setting} + basename_f

    return vol_name

def load_original_volume(stack, structure, volume_type=None):
    basename = get_original_volume_basename(stack=stack, volume_type=volume_type)
    vol_fp = os.path.join(VOLUME_ROOTDIR, '%(stack)s',
                          '%(basename)s',
                          'score_volumes',
                         '%(basename)s_%(struct)s.bp') % \
    {'stack':stack, 'basename':basename, 'struct':structure}
    
    volume = bp.unpack_ndarray_file(vol_fp).astype(np.float32)
    return volume

def load_original_volume_bbox(stack, structure, volume_type=None):
    basename = get_original_volume_basename(stack=stack, volume_type=volume_type)
    bbox_fp = os.path.join(VOLUME_ROOTDIR, '%(stack)s',
                          '%(basename)s',
                          'score_volumes',
                         '%(basename)s_%(struct)s_bbox.txt') % \
    {'stack':stack, 'basename':basename, 'struct':structure}
    
    volume_bbox = np.loadtxt(bbox_fp).astype(np.int)
    return volume_bbox    
    
def load_original_volume_all_known_structures_v2(stack,
                                                structures=None,
                                                sided=True,
                                                include_surround=False,
                                                 return_label_mappings=False,
                                                 name_or_index_as_key='name',
                                                 common_shape=True,
                                                 volume_type='annotationAsScore'):
    """
    Load original (un-transformed) volumes for all structures and optionally normalize them into a common shape.

    Args:
        common_shape (bool): If true, volumes are normalized to the same shape.

    Returns:
        If `common_shape` is True:
            if return_label_mappings is True, returns (volumes, common_bbox, structure_to_label, label_to_structure), volumes is dict.
            else, returns (volumes, common_bbox).
        If `common_shape` is False:
            if return_label_mappings is True, returns (dict of volume_bbox_tuples, structure_to_label, label_to_structure).
            else, returns volume_bbox_tuples.
    """

    if structures is None:
        if sided:
            if include_surround:
                structures = all_known_structures_sided_with_surround
            else:
                structures = all_known_structures_sided
        else:
            structures = all_known_structures

    loaded = False
    sys.stderr.write('Prior structure/index map not found. Generating a new one.\n')

    volumes = {}
    if not loaded:
        structure_to_label = {}
        label_to_structure = {}
        index = 1
    for structure in structures:
        try:

            if loaded:
                index = structure_to_label[structure]

            v = load_original_volume(stack=stack, structure=structure, volume_type=volume_type)
            b = load_original_volume_bbox(stack=stack, structure=structure, volume_type=volume_type)

            if name_or_index_as_key == 'name':
                volumes[structure] = (v,b)
            else:
                volumes[index] = (v,b)

            if not loaded:
                structure_to_label[structure] = index
                label_to_structure[index] = structure
                index += 1

        except Exception as e:
            sys.stderr.write('%s\n' % e)
            sys.stderr.write('Score volume for %s does not exist.\n' % structure)

    if common_shape:
        volumes_normalized, common_bbox = convert_vol_bbox_dict_to_overall_vol(vol_bbox_dict=volumes)

        if return_label_mappings:
            return volumes_normalized, common_bbox, structure_to_label, label_to_structure
        else:
            return volumes_normalized, common_bbox
    else:
        if return_label_mappings:
            return volumes, structure_to_label, label_to_structure
        else:
            return volumes
        
def load_transformed_volume(stack_m, stack_f,
                            warp_setting,
                            structure_f=None,
                            structure_m=None,
                            vol_type_m=None,
                            vol_type_f=None,
                            structure=None):
    basename = get_warped_volume_basename(**locals())
    if structure is not None:
        fn = basename + '_%s' % structure
    else:
        fn = basename
    fp = os.path.join(VOLUME_ROOTDIR, stack_m, basename, 'score_volumes', fn + '.bp')
    return bp.unpack_ndarray_file(fp)

def load_transformed_volume_bbox(stack_m, stack_f,
                                    warp_setting,
                                    structure_m=None,
                                    structure_f=None,
                                    vol_type_m='annotationAsScore',
                                    vol_type_f='annotationAsScore',
                                    structure=None,):
    basename = get_warped_volume_basename(**locals())
    if structure is not None:
        fn = basename + '_%s' % structure
    else:
        fn = basename
    fp = os.path.join(VOLUME_ROOTDIR, stack_m, basename, 'score_volumes', fn + '_bbox.txt')
    return np.loadtxt(fp)

def load_transformed_volume_all_known_structures_v2(stack_m,
                                                 stack_f,
                                                warp_setting,
                                                structures=None,
                                                sided=True,
                                                include_surround=False,
                                                 return_label_mappings=False,
                                                 name_or_index_as_key='name',
                                                 common_shape=True
):
    """
    Load transformed volumes for all structures and normalize them into a common shape.

    Args:
        common_shape (bool): If true, volumes are normalized to the same shape.

    Returns:
        If `common_shape` is True:
            if return_label_mappings is True, returns (volumes, common_bbox, structure_to_label, label_to_structure), volumes is dict.
            else, returns (volumes, common_bbox).
        If `common_shape` is False:
            if return_label_mappings is True, returns (dict of volume_bbox_tuples, structure_to_label, label_to_structure).
            else, returns volume_bbox_tuples.
    """

    if structures is None:
        if sided:
            if include_surround:
                structures = all_known_structures_sided_with_surround
            else:
                structures = all_known_structures_sided
        else:
            structures = all_known_structures

    loaded = False
    sys.stderr.write('Prior structure/index map not found. Generating a new one.\n')

    volumes = {}
    if not loaded:
        structure_to_label = {}
        label_to_structure = {}
        index = 1
    for structure in structures:
        try:

            if loaded:
                index = structure_to_label[structure]

            v = load_transformed_volume(stack_m=stack_m, vol_type_m='annotationAsScore',
                                                    stack_f=stack_f, vol_type_f='annotationAsScore',
                                                    warp_setting=warp_setting,
                                                    structure=structure)
            
            b = load_transformed_volume_bbox(stack_m=stack_m, vol_type_m='annotationAsScore',
                                                    stack_f=stack_f, vol_type_f='annotationAsScore',
                                                    warp_setting=warp_setting,
                                                    structure=structure)

            if name_or_index_as_key == 'name':
                volumes[structure] = (v,b)
            else:
                volumes[index] = (v,b)

            if not loaded:
                structure_to_label[structure] = index
                label_to_structure[index] = structure
                index += 1

        except Exception as e:
            sys.stderr.write('%s\n' % e)
            sys.stderr.write('Score volume for %s does not exist.\n' % structure)

    if common_shape:
        volumes_normalized, common_bbox = convert_vol_bbox_dict_to_overall_vol(vol_bbox_dict=volumes)

        if return_label_mappings:
            return volumes_normalized, common_bbox, structure_to_label, label_to_structure
        else:
            return volumes_normalized, common_bbox
    else:
        if return_label_mappings:
            return volumes, structure_to_label, label_to_structure
        else:
            return volumes

def get_alignment_result_filepath_v2(stack_f, stack_m, warp_setting, what, ext=None,
                                  detector_id_m=None, detector_id_f=None,
                                  prep_id_m=None, prep_id_f=None,
                                  structure_f=None, structure_m=None,
                                  vol_type_f='score', vol_type_m='score',
                                  downscale=15, trial_idx=None):
    reg_basename = get_warped_volume_basename(**locals())
    if what == 'parameters':
        ext = 'json'
    elif what == 'scoreHistory' or what == 'trajectory':
        ext = 'bp'
    elif what == 'scoreEvolution':
        ext = 'png'
    elif what == 'parametersWeightedAverage':
        ext = 'pkl'
    fp = os.path.join(REGISTRATION_PARAMETERS_ROOTDIR, stack_m, reg_basename, reg_basename + '_' + what + '.' + ext)
    return fp

def get_alignment_result_filepath_v3(stack_f, stack_m, warp_setting, what, ext=None,
                                  detector_id_m=None, detector_id_f=None,
                                  prep_id_m=None, prep_id_f=None,
                                  structure_f=None, structure_m=None,
                                  vol_type_f='score', vol_type_m='score',
                                  downscale=15, trial_idx=None):
    reg_basename = get_warped_volume_basename(**locals())
    if what == 'parameters':
        ext = 'json'
    elif what == 'scoreHistory' or what == 'trajectory':
        ext = 'bp'
    elif what == 'scoreEvolution':
        ext = 'png'
    elif what == 'parametersWeightedAverage':
        ext = 'pkl'
    fp = os.path.join(REGISTRATION_PARAMETERS_ROOTDIR, stack_m, reg_basename, reg_basename + '_' + what + '_RN' + '.' + ext)
    return fp

def load_alignment_parameters_v2(stack_f, stack_m, warp_setting,
                                  prep_id_m=None, prep_id_f=None,
                                  detector_id_m=None, detector_id_f=None,
                                  structure_f=None, structure_m=None,
                                  vol_type_f='score', vol_type_m='score',
                                  downscale=15, trial_idx=None):
    what = 'parameters'
    tf_param_fp = get_alignment_result_filepath_v2(**locals())
    return load_json(tf_param_fp)

def load_alignment_parameters_v3(stack_f, stack_m, warp_setting,
                                  prep_id_m=None, prep_id_f=None,
                                  detector_id_m=None, detector_id_f=None,
                                  structure_f=None, structure_m=None,
                                  vol_type_f='score', vol_type_m='score',
                                  downscale=15, trial_idx=None):
    what = 'parameters'
    tf_param_fp = get_alignment_result_filepath_v3(**locals())
    return load_json(tf_param_fp)
    
def get_stacy_markers_filepath(stack, structure):
    return os.path.join(DATA_ROOTDIR, 'markers', stack, stack + '_markers_%s.bp' % structure)