#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 20:01:58 2017

@author: asya
"""
import numpy as np
import sys
import time
from vis3d_utilities import *


try:
    import vtk
    from vtk.util import numpy_support
    import mcubes # https://github.com/pmneila/PyMCubes
except:
    sys.stderr.write('No vtk\n')

from skimage.measure import marching_cubes, correct_mesh_orientation, mesh_surface_area

import os
sys.path.append('/home/asya/Documents/Yuncong_code/utilities')


def actor_mesh_cut(polydata, color=(1.,1.,1.), wireframe=False, origin=(0,0,0),\
                   cut_plane_origin = (-60,0,0),cut_plane_normal = (1,0,0),thickness = 10,\
                   opacity = 0.3,opacityE = 0.5):
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
    cutterMapper=vtk.vtkPolyDataMapper()
    cutterMapper.SetInputConnection(cutter.GetOutputPort())
    
    planeCutActor=vtk.vtkActor()
    planeCutActor.GetProperty().SetColor(color[0],color[1],color[2])
    planeCutActor.GetProperty().SetLineWidth(1)
    planeCutActor.SetMapper(cutterMapper)
    
    
    ## clipActor will give contours sliced at the plane +- thickness
    plP = np.asarray(cut_plane_origin)+np.asarray(cut_plane_normal)*thickness
    plN = np.asarray(cut_plane_origin)-np.asarray(cut_plane_normal)*thickness
    planeP=vtk.vtkPlane() ## plane + thickness (upper bound)
    planeN=vtk.vtkPlane() ## plane - thickness (lower bound)
    planeP.SetOrigin(plP.tolist())
    planeN.SetOrigin(plN.tolist())
    planeP.SetNormal(cut_plane_normal)
    planeN.SetNormal(cut_plane_normal)
    
    clipper = vtk.vtkClipPolyData()
    clipper.SetInputData(polydata_shifted)
    clipper.SetClipFunction(planeP)
    clipper.GenerateClipScalarsOn()
    clipper.GenerateClippedOutputOn()
    clipper.SetValue(0.5)
    clipper.Update()
    
    clipperB = vtk.vtkClipPolyData()
    clipperB.SetInputConnection(clipper.GetClippedOutputPort())
    clipperB.SetClipFunction(planeN)
    clipperB.GenerateClipScalarsOn()
    clipperB.GenerateClippedOutputOn()
    clipperB.SetValue(0.5)
    
    clipMapper = vtk.vtkPolyDataMapper()
    clipMapper.SetInputConnection(clipperB.GetOutputPort())
    clipMapper.ScalarVisibilityOff()
    clipActor = vtk.vtkActor()
    clipActor.SetMapper(clipMapper)
    clipActor.GetProperty().SetOpacity(opacity)
    clipActor.GetProperty().SetColor(color[0],color[1],color[2])

    ## first surface
    cutEdges = vtk.vtkCutter()
    cutEdges.SetInputData(polydata_shifted)
    cutEdges.SetCutFunction(planeP)
    cutEdges.GenerateCutScalarsOn()
    cutEdges.SetValue(0, 0.5)
    cutStrips = vtk.vtkStripper()
    cutStrips.SetInputConnection(cutEdges.GetOutputPort())
    cutStrips.Update()
    cutPoly = vtk.vtkPolyData()
    cutPoly.SetPoints(cutStrips.GetOutput().GetPoints())
    cutPoly.SetPolys(cutStrips.GetOutput().GetLines())
    
    cutTriangles = vtk.vtkTriangleFilter()
    cutTriangles.SetInputData(cutPoly)
    cutMapper = vtk.vtkPolyDataMapper()
    cutMapper.SetInputData(cutPoly)
    cutMapper.SetInputConnection(cutTriangles.GetOutputPort())
    cutActorP = vtk.vtkActor()
    cutActorP.SetMapper(cutMapper)
    cutActorP.GetProperty().SetColor(0,color[1],color[2])
    cutActorP.GetProperty().SetOpacity(opacityE)
    
    ## second surface
    cutEdges = vtk.vtkCutter()
    cutEdges.SetInputData(polydata_shifted)
    cutEdges.SetCutFunction(planeN)
    cutEdges.GenerateCutScalarsOn()
    cutEdges.SetValue(0, 0.5)
    cutStrips = vtk.vtkStripper()
    cutStrips.SetInputConnection(cutEdges.GetOutputPort())
    cutStrips.Update()
    cutPoly = vtk.vtkPolyData()
    cutPoly.SetPoints(cutStrips.GetOutput().GetPoints())
    cutPoly.SetPolys(cutStrips.GetOutput().GetLines())
    
    cutTriangles = vtk.vtkTriangleFilter()
    cutTriangles.SetInputData(cutPoly)
    cutMapper = vtk.vtkPolyDataMapper()
    cutMapper.SetInputData(cutPoly)
    cutMapper.SetInputConnection(cutTriangles.GetOutputPort())
    cutActorN = vtk.vtkActor()
    cutActorN.SetMapper(cutMapper)
    cutActorN.GetProperty().SetColor(0,color[1],color[2])
    cutActorN.GetProperty().SetOpacity(opacityE)
    

    return (planeCutActor,clipActor,cutActorP,cutActorN)



def actor_mesh_cut_brainstem(polydata, color=(1.,1.,1.), wireframe=False, origin=(0,0,0),\
                   cut_plane_origin = (-60,0,0),cut_plane_normal = (1,0,0),cut_plane_width = 2):
    """
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
    
    # Set up cut plane
    plane=vtk.vtkPlane()
    plane.SetOrigin(cut_plane_origin)
    plane.SetNormal(cut_plane_normal)

    # Set up cutter for volume display
    cutter=vtk.vtkCutter()
    cutter.SetCutFunction(plane)
    cutter.SetInputData(polydata_shifted)
    cutter.Update()
    cutterMapper=vtk.vtkPolyDataMapper()
    cutterMapper.SetInputConnection(cutter.GetOutputPort())
    
    # Make actor for plane display in volume
    planeCutActor=vtk.vtkActor()
    planeCutActor.GetProperty().SetColor(color)
    planeCutActor.GetProperty().SetLineWidth(cut_plane_width)
    planeCutActor.SetMapper(cutterMapper)

    # Set up cutter for slice display
    cutEdges = vtk.vtkCutter()
    cutEdges.SetInputData(polydata_shifted)
    cutEdges.SetCutFunction(plane)
    cutStrips = vtk.vtkStripper()
    cutStrips.SetInputConnection(cutEdges.GetOutputPort())
    cutStrips.Update()
    cutPoly = vtk.vtkPolyData()
    cutPoly.SetPoints(cutStrips.GetOutput().GetPoints())
    cutPoly.SetPolys(cutStrips.GetOutput().GetLines())
    
    cutMapper = vtk.vtkPolyDataMapper()
    cutMapper.SetInputData(cutPoly)
    
    sliceCutActor = vtk.vtkActor()
    sliceCutActor.GetProperty().SetEdgeColor(color)
    sliceCutActor.GetProperty().SetLineWidth(2)
    sliceCutActor.GetProperty().EdgeVisibilityOn()
    sliceCutActor.GetProperty().SetSpecular(0.0)
    sliceCutActor.GetProperty().SetDiffuse(0.0)
    sliceCutActor.SetMapper(cutMapper)

    return (planeCutActor, sliceCutActor)


def actor_mesh_cut_probability(volume_in, spacing = (1.,1.,1.), origin=[0.,0.,0.], color=(1.,1.,1.),\
                   cut_plane_origin = (-60,0,0),cut_plane_normal = (1,0,0)):
    """
    Args:
        color (float array): rgb between 0 and 1.
        origin: the initial shift for the mesh.
    """

    plane=vtk.vtkPlane()
    plane.SetOrigin(cut_plane_origin)
    plane.SetNormal(cut_plane_normal)
    
    imagedata = volume_to_imagedata_respace(volume_in, origin=origin, spacing = spacing, auxdata=None)

    cutter=vtk.vtkCutter()
    cutter.SetCutFunction(plane)
    cutter.SetInputData(imagedata)
    cutter.Update()
    cutterMapper=vtk.vtkPolyDataMapper()
    cutterMapper.SetInputConnection(cutter.GetOutputPort())
    
    planeCutActor1=vtk.vtkActor()
    planeCutActor1.GetProperty().SetColor(color[0],color[1],color[2])
    planeCutActor1.GetProperty().SetLineWidth(1)
    planeCutActor1.SetMapper(cutterMapper)
    
    return planeCutActor1



def launch_vtk_multi(actors1, actors2,cut_plane_normal=(0,-1,0),background_color=(1.,1.,1.)):

    rw = vtk.vtkRenderWindow()
    rw.SetSize(1500,1080)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(rw)

    ren1 = vtk.vtkRenderer()
    ren1.SetBackground(background_color)
    rw.AddRenderer(ren1)
    ren1.SetViewport(0,0,0.5,1)
    
    for actor in actors1:
        if actor is not None:
            ren1.AddActor(actor)
            
    camera = vtk.vtkCamera()
    camera.SetViewUp(0, -1, 0)
    camera.SetPosition(-20, -20, -20)
    camera.SetFocalPoint(1, 1, 1)
    ren1.SetActiveCamera(camera)
    ren1.ResetCamera()
            
    ren1 = vtk.vtkRenderer()
    ren1.SetBackground(background_color)
    rw.AddRenderer(ren1)
    ren1.SetViewport(0.5,0,1,1)
    
    for actor in actors2:
        if actor is not None:
            ren1.AddActor(actor)
            
    camera = vtk.vtkCamera()
    camera.SetViewUp(0,-1,0)
    camera.SetPosition(cut_plane_normal)
    camera.SetFocalPoint(0.,0.,0.)
    
    ren1.SetActiveCamera(camera)
    ren1.ResetCamera()
    
    rw.Render()
    rw.SetWindowName('RW: Multiple ViewPorts')
    iren.Start()
    
    ren1 = vtk.vtkRenderer()
    
def launch_vtk_3view(vol_actors, slice_actors, cut_plane_normal,background_color=(1.,1.,1.)):

    rw = vtk.vtkRenderWindow()
    rw.SetSize(1500,1080)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(rw)

    ##   Volume
    ren1 = vtk.vtkRenderer()
    ren1.SetBackground(background_color)
    rw.AddRenderer(ren1)
    ren1.SetViewport(0,0.5,0.5,1)
    
    for actor in vol_actors:
        if actor is not None:
            ren1.AddActor(actor)
            
    camera = vtk.vtkCamera()
    camera.SetViewUp(0, -1, 0)
    camera.SetPosition(-20, -20, -20)
    camera.SetFocalPoint(1, 1, 1)
    ren1.SetActiveCamera(camera)
    ren1.ResetCamera()
            
    ##   Slice 1
    ren1 = vtk.vtkRenderer()
    ren1.SetBackground(background_color)
    rw.AddRenderer(ren1)
    ren1.SetViewport(0.5,0.5,1,1)
    
    for actor in slice_actors[0]:
        if actor is not None:
            ren1.AddActor(actor)
            
    camera = vtk.vtkCamera()
    camera.SetViewUp(0,-1,0)
    camera.SetPosition(cut_plane_normal[0])
    camera.SetFocalPoint(0.,0.,0.)
    ren1.SetActiveCamera(camera)
    ren1.ResetCamera()
    
     ##   Slice 2
    ren1 = vtk.vtkRenderer()
    ren1.SetBackground(background_color)
    rw.AddRenderer(ren1)
    ren1.SetViewport(0,0,0.5,0.5)
    
    for actor in slice_actors[1]:
        if actor is not None:
            ren1.AddActor(actor)
            
    camera = vtk.vtkCamera()
    camera.SetViewUp(0,-1,0)
    camera.SetPosition(cut_plane_normal[1])
    camera.SetFocalPoint(0.,0.,0.)
    ren1.SetActiveCamera(camera)
    ren1.ResetCamera()
    
    
     #   Slice 3
    ren1 = vtk.vtkRenderer()
    ren1.SetBackground(background_color)
    rw.AddRenderer(ren1)
    ren1.SetViewport(0.5,0,1,0.5)
    
    for actor in slice_actors[2]:
        if actor is not None:
            ren1.AddActor(actor)
            
    camera = vtk.vtkCamera()
    camera.SetViewUp(0,-1,0)
    camera.SetPosition(cut_plane_normal[2])
    camera.SetFocalPoint(0.,0.,0.)
    ren1.SetActiveCamera(camera)
    ren1.ResetCamera()
    
    rw.Render()
    rw.SetWindowName('RW: Multiple ViewPorts')
    iren.Start()
    
    ren1 = vtk.vtkRenderer()

def volume_to_polydata_brainstem(volume, origin=(0,0,0), num_simplify_iter=0, smooth=False, level=0., min_vertices=200):
    """
    Convert a volume to a mesh.
    
    Args:
        min_vertices (int): minimum number of vertices. Simplification will stop if the number of vertices drops below this value.
    """

    vol = volume > level

    # vol_padded = np.zeros(vol.shape+(10,10,10), np.bool)
    # vol_padded[5:-5, 5:-5, 5:-5] = vol
    vol_padded = np.pad(vol, ((5,5),(5,5),(5,5)), 'constant') # need this otherwise the sides of volume will not close and expose the hollow inside of structures

    t = time.time()
    vs, fs = mcubes.marching_cubes(vol_padded, 0) # more than 5 times faster than skimage.marching_cube + correct_orientation
    sys.stderr.write('marching cube: %.2f seconds\n' % (time.time() - t))

    # t = time.time()
    # vs, faces = marching_cubes(vol_padded, 0) # y,x,z
    # sys.stderr.write('marching cube: %.2f seconds\n' % (time.time() - t))

    # t = time.time()
    # fs = correct_mesh_orientation(vol_padded, vs, faces)
    # sys.stderr.write('correct orientation: %.2f seconds\n' % (time.time() - t))

    vs = vs[:, [1,0,2]] + origin - (5,5,5)
    # vs = vs[:, [1,0,2]] + origin

    # t = time.time()
    # area = mesh_surface_area(vs, fs)
    # # print 'area: %.2f' % area
    # sys.stderr.write('compute surface area: %.2f seconds\n' % (time.time() - t)) #

    t = time.time()

    polydata = mesh_to_polydata(vs, fs)
    
    for simplify_iter in range(num_simplify_iter):

        t = time.time()

        deci = vtk.vtkQuadricDecimation()
        deci.SetInputData(polydata)

        deci.SetTargetReduction(0.01)
        # 0.8 means each iteration causes the point number to drop to 20% the original

        deci.Update()

        polydata = deci.GetOutput()
        tsmooth = time.time()
        if smooth:

            smoother = vtk.vtkWindowedSincPolyDataFilter()
    #         smoother.NormalizeCoordinatesOn()
            smoother.SetPassBand(.1)
            smoother.SetNumberOfIterations(20)
            smoother.SetInputData(polydata)
            smoother.Update()

            polydata = smoother.GetOutput()

        sys.stderr.write('SMOOTH !!!!! %.2f seconds\n' % (time.time() - tsmooth)) #

        n_pts = polydata.GetNumberOfPoints()
        sys.stderr.write('simplify %d @ %d: %.2f seconds\n' % (simplify_iter, n_pts, time.time() - t)) #
        
        if polydata.GetNumberOfPoints() < min_vertices:
            break

    return polydata


def actor_volume_respace(volume_in, what,spacing = [1,1,1], auxdata=None, origin=(0,0,0), c=(1,1,1), tb_colors=None, tb_opacity=.05):
    """
    Add in a different spacing option, otherwise same as yunncong's function
    Args:
        what (str): tb, score or probability. A caveat when what="probability" is that zero-valued voxels are not transparent, so later actors will block previous actors.
        c (3-tuple): color
        tb_colors (dict {int: 3-tuple}): color transfer function that maps intensity value to color tuple.
        
    """

    imagedata = volume_to_imagedata_respace(volume_in, origin=origin, spacing = spacing, auxdata=auxdata)

    volumeMapper = vtk.vtkSmartVolumeMapper()
    volumeMapper.SetBlendModeToComposite()
    

    volumeMapper.SetInputData(imagedata)
    volumeProperty = vtk.vtkVolumeProperty()
    

    if what == 'tb':

        compositeOpacity = vtk.vtkPiecewiseFunction()
        compositeOpacity.AddPoint(0.0, 0.)

        if tb_colors is not None:
            for v, c in sorted(tb_colors.items()):
                vl = v - .5
                vr = v + .5
                cp1 = vl-.25
                cp2 = vr-.25
                compositeOpacity.AddPoint(cp1, .5*cp1/200., .5, 1.)
                compositeOpacity.AddPoint(v, 1., .5, 1.)
                compositeOpacity.AddPoint(cp2, .5*cp2/200., .5, 1.)
            compositeOpacity.AddPoint(vr, .5*vr/200.)
        compositeOpacity.AddPoint(240., tb_opacity)
        compositeOpacity.AddPoint(255.0, tb_opacity)
        # compositeOpacity.AddPoint(240., .5)
        # compositeOpacity.AddPoint(255.0, .5)

        color = vtk.vtkColorTransferFunction()
        color.AddRGBPoint(0.0, 0,0,0)

        if tb_colors is not None:
            for v, c in sorted(tb_colors.items()):
                vl = v - .5
                vr = v + .5
                cp1 = vl-.25
                cp2 = vr-.25
                color.AddRGBPoint(cp1, .5*cp1/200., .5*cp1/200., .5*cp1/200., .5, 1.)
                color.AddRGBPoint(v, c[0], c[1], c[2], .5, 1.)
                color.AddRGBPoint(cp2, .5*cp2/200., .5*cp2/200., .5*cp2/200., .5, 1.)
            color.AddRGBPoint(vr, .5*vr/200., .5*vr/200., .5*vr/200.)
        color.AddRGBPoint(200.0, .5,.5,.5)
        color.AddRGBPoint(255.0, 1,1,1)

    # volumeGradientOpacity = vtk.vtkPiecewiseFunction()
    # volumeGradientOpacity.AddPoint(0,   0.0)
    # volumeGradientOpacity.AddPoint(1,  0.5)
    # volumeGradientOpacity.AddPoint(2, 1.0)

    elif what == 'score':
        
        volumeProperty.IndependentComponentsOff()

        compositeOpacity = vtk.vtkPiecewiseFunction()
        compositeOpacity.AddPoint(0.0, 0.0)
        # compositeOpacity.AddPoint(0.95, 0.0)
        compositeOpacity.AddPoint(1.0, 1.0)

        color = vtk.vtkColorTransferFunction()
        
        color.AddRGBPoint(0.0, 0,0,0)

        if tb_colors is not None:
            for v, c in sorted(tb_colors.items()):
                vl = v - .5
                vr = v + .5
                cp1 = vl-.25
                cp2 = vr-.25
                color.AddRGBPoint(cp1, .5*cp1/200., .5*cp1/200., .5*cp1/200., .5, 1.)
                color.AddRGBPoint(v, c[0], c[1], c[2], .5, 1.)
                color.AddRGBPoint(cp2, .5*cp2/200., .5*cp2/200., .5*cp2/200., .5, 1.)
            color.AddRGBPoint(vr, .5*vr/200., .5*vr/200., .5*vr/200.)
        color.AddRGBPoint(200.0, .5,.5,.5)
        color.AddRGBPoint(255.0, 1,1,1)
        
#         color.AddRGBPoint(0.0, c[0], c[1], c[2])
#         color.AddRGBPoint(1.0, c[0], c[1], c[2])
        
        # volumeGradientOpacity = vtk.vtkPiecewiseFunction()
        # volumeGradientOpacity.AddPoint(0,  0.0)
        # volumeGradientOpacity.AddPoint(10,  1.0)
        # volumeGradientOpacity.AddPoint(2, 1.0)

    elif what == 'probability':
        compositeOpacity = vtk.vtkPiecewiseFunction()
        compositeOpacity.AddPoint(0.0, 0.)
        #compositeOpacity.AddPoint(.01, .9)
        compositeOpacity.AddPoint(1.0, 1.0)

        r,g,b = c

        color = vtk.vtkColorTransferFunction()
        color.AddRGBPoint(0.0, 1.,1.,1.)
        color.AddRGBPoint(.05, .5,.5,.5)
        color.AddRGBPoint(1., 0.,0.,0.)

        # lookupTable = vtkLookupTable()
        # lookupTable.SetNumberOfTableValues(2);
        # lookupTable.SetRange(0.0,1.0);
        # lookupTable.SetTableValue( 0, 0.0, 0.0, 0.0, 0.0 ); #label 0 is transparent
        # lookupTable.SetTableValue( 1, 0.0, 1.0, 0.0, 1.0 ); #label 1 is opaque and green
        # lookupTable.Build()
        #
        # mapTransparency = vtkImageMapToColors()
        # mapTransparency.SetLookupTable(lookupTable)
        # mapTransparency.PassAlphaToOutputOn()
        # mapTransparency.SetInputData(maskImage)
        #
        # maskActor = vtkImageActor()
        # maskActor.GetMapper().SetInputConnection(mapTransparency.GetOutputPort())
        #
        #
        # volumeProperty.SetScalarOpacity(compositeOpacity)
        # volumeProperty.SetColor(color)
        #
        # volume = vtk.vtkVolume()
        # volume.SetMapper(volumeMapper)
        # volume.SetProperty(volumeProperty)
        #
        # return volume

    else:
        sys.stderr.write('Color/opacity profile not recognized.\n')

    # volumeProperty.ShadeOff()
    volumeProperty.SetColor(color)
    volumeProperty.SetScalarOpacity(compositeOpacity)

    volume = vtk.vtkVolume()
    volume.SetMapper(volumeMapper)
    volume.SetProperty(volumeProperty)

    return volume

def volume_to_imagedata_respace(arr, origin=(0,0,0), spacing = [1.,1.,1.], auxdata = None):
    """
    Args:
        arr (3d-array of uint8 or float32):
        origin (3-tuple): the origin coordinate of the given volume
        
    Returns:
        (vtkImageData): Each point (in vtk parlance) gets ONE scalar value which is the value of an input volume voxel. Respects the (x,y,z) dimension ordering.
    """

    imagedata = vtk.vtkImageData()
    imagedata.SetDimensions([arr.shape[0], arr.shape[1], arr.shape[2]])
    imagedata.SetSpacing(spacing[0],spacing[1],spacing[2]) ###SPACING
    imagedata.SetOrigin(origin[0], origin[1], origin[2])

    v3 = np.transpose(arr,[2,1,0]) #[2,0,1])
    v3 = v3.flatten()
    if auxdata is not None:
        auxdata = np.transpose(auxdata, [2,0,1])
        v3 = np.column_stack([v3, auxdata.flatten()])
        
    if arr.dtype == np.uint8:
        t = vtk.VTK_UNSIGNED_CHAR
    elif arr.dtype == np.float32:
        t = vtk.VTK_FLOAT
    else:
        raise Exception('Data type must be uint8 or float32.')
    
    imagedata.GetPointData().SetScalars(numpy_support.numpy_to_vtk(v3, deep=True, array_type=t)) # deep copy must be true
    return imagedata


def make_scaleBar_actor(O,N):
    x0 = float(O[0])
    x1 = float(O[1]+750)
    x2 = float(O[2]+50)
    #(x0-y0)*N[0]+(x1-y1)*N[1] + (x2-y2)*N[2] = 0
    #(x0-y0)^2+(x1-y1)^2 + (x2-y2)^2 = 100^2
    #y1  = (x2-y2)*N[2]/N[1]+x1
    #y1  = x1 + sqrt(100^2 - (x2-y2)^2)
        
    if max(N)==N[0]:
        y0 = x0
        if N[1]==0. and N[2] == 0.:
            y2 = x2 + 100
            y1 = x1
        elif N[1] == 0:
            y1 = x1 + 100/(1+(N[1]/N[2]))
            y2 = (x1-y1)*(N[1]/N[2])+x1
        else:
            y2 = x2 + 100/(1+(N[2]/N[1]))
            y1 = (x2-y2)*(N[2]/N[1])+x1
            
    elif max(N)==N[1]:
        y1 = x1
        if N[0]==0. and N[2] == 0.:
            y2 = x2
            y0 = x0 + 100
        elif N[2] == 0:
            y2 = x2 + 100/(1+(N[2]/N[0]))
            y0 = (x2-y2)*(N[2]/N[0])+x0
        else:
            y2 = x2 + 100/(1+(N[2]/N[0]))
            y0 = (x2-y2)*(N[2]/N[0])+x0
    else:
        y2 = x2
        if N[1]==0. and N[0] == 0.:
            y1 = x1 + 100
            y0 = x0
        elif N[1] == 0:
            y1 = x1 + 100/(1+(N[1]/N[0]))
            y0 = (x1-y1)*(N[1]/N[0])+x0
        else:
            y1 = x1 + 100/(1+(N[1]/N[0]))
            y0 = (x1-y1)*(N[1]/N[0])+x0
    sbSource = vtk.vtkLineSource()
    sbSource.SetPoint1(x0,x1,x2)
    sbSource.SetPoint2(y0,y1,y2)
    sbMapper = vtk.vtkPolyDataMapper()
    sbMapper.SetInputConnection(sbSource.GetOutputPort())

    #create an actor
    sbActor = vtk.vtkActor()
    sbActor.SetMapper(sbMapper)
    sbActor.GetProperty().SetColor(0,0,0)
    sbActor.GetProperty().SetOpacity(1)
    sbActor.GetProperty().SetLineWidth(2)
    return sbActor
            
def make_scaleBar_coronal_actor(origin):
    sbSource = vtk.vtkLineSource()
    sbSource.SetPoint1(origin[0],750.,50.)
    sbSource.SetPoint2(origin[0],750.,150.)
    sbMapper = vtk.vtkPolyDataMapper()
    sbMapper.SetInputConnection(sbSource.GetOutputPort())

    #create an actor
    sbActor = vtk.vtkActor()
    sbActor.SetMapper(sbMapper)
    sbActor.GetProperty().SetColor(0,0,0)
    sbActor.GetProperty().SetOpacity(1)
    sbActor.GetProperty().SetLineWidth(2)
    return sbActor