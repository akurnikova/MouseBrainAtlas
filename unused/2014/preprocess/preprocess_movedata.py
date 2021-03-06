#!/usr/bin/python

import os
import sys
import argparse
import json
from operator import itemgetter
from subprocess import call

import numpy as np
from matplotlib import pyplot as plt
from skimage.color import rgb2gray, rgb2hsv
from skimage.io import imread, imsave
from skimage.morphology import remove_small_objects, binary_dilation, square
from skimage.measure import regionprops, label
from skimage.restoration import denoise_bilateral
from skimage.transform import rescale
from PIL import Image

from joblib import Parallel, delayed

from generate_mask import foreground_mask_morphsnakes
from preprocess_utility import create_if_not_exists, execute_command

# class NumpyAwareJSONEncoder(json.JSONEncoder):
# 	def default(self, obj):
# 	    if isinstance(obj, np.ndarray) and obj.ndim == 1:
# 	        return obj.tolist()
# 	    return json.JSONEncoder.default(self, obj)

def refine_save_highres_maskedimgs_parallel(stack, slide_ind, resol, 
											bbox_x03125, section_ind, slide_img):

	masked_sectionimg_path = os.path.join(os.environ['GORDON_SECTIONDATA_DIR'], stack, 
		'autogen_maskedimg_x0.3125', "%s_x0.3125_%02d_%d.tif" % (stack, slide_ind, section_ind))

	alpha = imread(masked_sectionimg_path)[..., -1]

	l = ['x0.3125', 'x1.25', 'x5'].index(resol)

	scaling = 4**l
	minr, minc, maxr, maxc = bbox_x03125 * scaling
	cur_resol_sectionimg = slide_img[minr:maxr, minc:maxc]
	cur_resol_alpha = rescale(alpha, scaling, order=0) > 0
	
	# if slide_ind == 10 and section_ind == 1:	
	# 	imsave('/tmp/cur_resol_sectionimg.tif', cur_resol_sectionimg)
	# 	imsave('/tmp/cur_resol_alpha.tif', cur_resol_alpha)

	sectionimg_hsv = rgb2hsv(cur_resol_sectionimg)

	# sys.stderr.write('start masking %d %d\n' % (slide_ind, section_ind))

	# cur_resol_alpha, new_bbox = foreground_mask_morphsnakes(sectionimg_hsv[...,1], 
	# 		levelset=cur_resol_alpha, max_iters=500, scaling=scaling)

	# sys.stderr.write('finish masking %d %d\n' % (slide_ind, section_ind))

	cur_resol_alpha = 255 * cur_resol_alpha.astype(np.uint8)

	cur_resol_masked_sectionimg = np.dstack([cur_resol_sectionimg, cur_resol_alpha])

	cur_resol_masked_sectionimg_path = os.path.join(os.environ['GORDON_SECTIONDATA_DIR'], stack, 
						resol, '%s_%s_%02d_%d.tif' % (stack, resol, slide_ind, section_ind))
	imsave(cur_resol_masked_sectionimg_path, cur_resol_masked_sectionimg)

	del cur_resol_sectionimg
	del cur_resol_alpha
	del cur_resol_masked_sectionimg
	# del sectionimg_hsv


def refine_save_highres_maskedimgs(stack, slide_ind):
	
	autogen_bbox_dir = os.path.join(os.environ['GORDON_SECTIONDATA_DIR'], stack, 'autogen_bbox_x0.3125')

	os.chdir(autogen_bbox_dir)

	# _, stack, _, slide_str = bbox_file[:-4].split('_')
	# slide_ind = int(slide_str)
# 
	print 'slide', slide_ind

	bbox_file = '_'.join(['bbox', stack, 'x0.3125', '%02d'%slide_ind]) + '.txt'

	bboxes_x03125 = np.loadtxt(bbox_file).astype(np.int)

	curr_resol_slide_image = {}

	# for resol in ['x0.3125', 'x1.25', 'x5', 'x20']:
	for resol in ['x0.3125', 'x1.25', 'x5']:
	# for resol in ['x20']:
		curr_resol_slide_path = os.path.join(os.environ['GORDON_SLIDEDATA_DIR'], stack, resol, "%s_%02d_%s_z0.tif" % (stack, slide_ind, resol))
		curr_resol_slide_image[resol] = np.asarray(Image.open(curr_resol_slide_path).convert('RGB'))
		# print curr_resol_slide_image[resol].shape

	# section_counter = first_section_ind - 1

	for section_ind, bbox in enumerate(bboxes_x03125):

		# section_counter += 1

		if (bbox < 0).all(): continue

		masked_sectionimg_path = os.path.join(os.environ['GORDON_SECTIONDATA_DIR'], stack, 'autogen_maskedimg_x0.3125', "%s_x0.3125_%02d_%d.tif" % (stack, slide_ind, section_ind))

		alpha = imread(masked_sectionimg_path)[..., -1]

		# for l, resol in enumerate(['x0.3125', 'x1.25', 'x5', 'x20']):
		for l, resol in enumerate(['x0.3125', 'x1.25', 'x5']):
			
			if l != 2: continue

			scaling = 4**l
			minr, minc, maxr, maxc = bbox * scaling
			cur_resol_sectionimg = curr_resol_slide_image[resol][minr:maxr+scaling, minc:maxc+scaling]
			cur_resol_alpha = rescale(alpha, scaling, order=0) > 0
			
			sectionimg_hsv = rgb2hsv(cur_resol_sectionimg)
			# if resol == 'x5':
				# cur_resol_alpha, new_bbox = foreground_mask_morphsnakes(sectionimg_hsv[...,1], 
				# levelset=cur_resol_alpha, max_iters=500, min_iters=10, diff_thresh=200)					
			cur_resol_alpha, new_bbox = foreground_mask_morphsnakes(sectionimg_hsv[...,1], 
					levelset=cur_resol_alpha, max_iters=500)
			# else:
				# cur_resol_alpha, new_bbox = foreground_mask_morphsnakes(sectionimg_hsv[...,1], 
				# levelset=cur_resol_alpha, max_iters=500, min_iters=10, diff_thresh=40)

			print 'bbox', new_bbox

			cur_resol_alpha = 255 * cur_resol_alpha.astype(np.uint8)

			cur_resol_masked_sectionimg = np.dstack([cur_resol_sectionimg, cur_resol_alpha])

			cur_resol_masked_sectionimg_path = os.path.join(os.environ['GORDON_SECTIONDATA_DIR'], stack, resol, '%s_%s_%02d_%d.tif' % (stack, resol, slide_ind, section_ind))
			imsave(cur_resol_masked_sectionimg_path, cur_resol_masked_sectionimg)

			del cur_resol_sectionimg
			del cur_resol_alpha
			del cur_resol_masked_sectionimg
			del sectionimg_hsv

	for x in curr_resol_slide_image:
		del x


# def refine_save_highres_maskedimgs(bbox_dir, mask_dir):
	# os.chdir(bbox_dir)

	# section_counter = -1
	
	# # summary_dict = {}

	# all_bboxfiles = sorted([bbox_file for bbox_file in os.listdir('.') if bbox_file.endswith('txt')])

	# for bbox_file in all_bboxfiles:

	# 	if section_counter > 3: break
		
	# 	_, stack, _, slide_str = bbox_file[:-4].split('_')
	# 	slide_ind = int(slide_str)

	# 	print 'slide', slide_ind

	# 	bboxes_x03125 = np.loadtxt(bbox_file).astype(np.int)

	# 	# slide_imagepath = os.path.join(os.environ['GORDON_SLIDEDATA_DIR'], stack, '', "CC35_%02d_x0.3125_z0.tif" % slide_ind)		# DavidData2015slides/CC35/x0.3125/CC35_45_x0.3125_z0.tif
	# 	# slide_image = np.array(Image.open(slide_imagepath).convert('RGB'))
	# 	# slide_width, slide_height = slide_image.shape[:2]

	# 	# section_imagepath = os.path.join(os.environ['GORDON_SECTIONDATA_DIR'], stack, resol, 'autogen_maskedimg_x0.3125', '%s_%s_%04d.tif' % (stack, resol, section_counter))		# DavidData2015slides/CC35/x0.3125/CC35_45_x0.3125_z0.tif

	# 	curr_resol_slide_image = {}

	# 	# for resol in ['x0.3125', 'x1.25', 'x5', 'x20']:
	# 	for resol in ['x0.3125', 'x1.25', 'x5']:
	# 	# for resol in ['x20']:
	# 		curr_resol_slide_path = os.path.join(os.environ['GORDON_SLIDEDATA_DIR'], stack, resol, "%s_%02d_%s_z0.tif" % (stack, slide_ind, resol))
	# 		curr_resol_slide_image[resol] = np.asarray(Image.open(curr_resol_slide_path).convert('RGB'))
	# 		# print curr_resol_slide_image[resol].shape

	# 	for section_ind, bbox in enumerate(bboxes_x03125):

	# 		section_counter += 1

	# 		if (bbox < 0).all(): continue

	# 		masked_sectionimg_path = os.path.join(os.environ['GORDON_SECTIONDATA_DIR'], stack, 'autogen_maskedimg_x0.3125', "%s_x0.3125_%02d_%d.tif" % (stack, slide_ind, section_ind))
	# 		# alpha = 255 * imread(masked_sectionimg_path)[..., -1].astype(np.uint8)

	# 		alpha = imread(masked_sectionimg_path)[..., -1]

	# 		# minr, minc, maxr, maxc = bbox
	# 		# section_dict = {}

	# 		# for l, resol in enumerate(['x0.3125', 'x1.25', 'x5', 'x20']):
	# 		for l, resol in enumerate(['x0.3125', 'x1.25', 'x5']):
				
	# 			if l != 2: continue

	# 			print resol

	# 			scaling = 4**l

	# 			curr_resol_bbox = bbox * scaling

	# 			# section_dict[resol] = {'slide_path': curr_resol_slide_path,
	# 			# 				# 'bbox': bboxes_x03125.astype(np.float) / np.r_[slide_image.shape[:2], slide_image.shape[:2]] * np.r_[curr_resol_slide_shape, curr_resol_slide_shape],
	# 			# 				'bbox': curr_resol_bbox,
	# 			# 				# 'mask_path': os.path.join(os.environ['GORDON_SECTIONDATA_DIR'], stack, resol, 'masks', '%s_%s_%04d_mask.png' % (stack, resol, section_counter)),
	# 			# 				'sectionimg_path': sectionimg_path
	# 			# 				}

	# 			minr, minc, maxr, maxc = curr_resol_bbox

	# 			# print bbox_file, slide_ind, minr, minc, maxr, maxc

	# 			# print bbox, curr_resol_bbox, curr_resol_slide_image[resol].shape, alpha.shape

	# 			cur_resol_sectionimg = curr_resol_slide_image[resol][minr:maxr+scaling, minc:maxc+scaling]
	# 			cur_resol_alpha = rescale(alpha, scaling, order=0) > 0
				
	# 			# if l >= 1:
	# 				# for dilation_i in range(3):
	# 				# 	cur_resol_alpha = binary_dilation(cur_resol_alpha, selem=square(10))
			
	# 			# if resol == 'x20':
	# 			# 	pass
	# 			# else:
	# 			sectionimg_hsv = rgb2hsv(cur_resol_sectionimg)
	# 			if resol == 'x5':
	# 				cur_resol_alpha, new_bbox = foreground_mask_morphsnakes(sectionimg_hsv[...,1], 
	# 				levelset=cur_resol_alpha, max_iters=500, min_iters=10, diff_thresh=200)					
	# 			else:
	# 				cur_resol_alpha, new_bbox = foreground_mask_morphsnakes(sectionimg_hsv[...,1], 
	# 				levelset=cur_resol_alpha, max_iters=500, min_iters=10, diff_thresh=40)

	# 			print 'bbox', new_bbox

	# 			cur_resol_alpha = 255 * cur_resol_alpha.astype(np.uint8)

	# 			cur_resol_masked_sectionimg = np.dstack([cur_resol_sectionimg, cur_resol_alpha])

	# 			# sectionimg_path = os.path.join(os.environ['GORDON_SECTIONDATA_DIR'], stack, resol, 'images', '%s_%s_%04d.tif' % (stack, resol, section_counter))
	# 			cur_resol_masked_sectionimg_path = os.path.join(os.environ['GORDON_SECTIONDATA_DIR'], stack, resol, '%s_%s_%04d.tif' % (stack, resol, section_counter))
	# 			imsave(cur_resol_masked_sectionimg_path, cur_resol_masked_sectionimg)

	# 			# execute_command('convert -resize 400\% %s %s' % (mask_path, section_dict[resol]['mask_path']))

	# 		# summary_dict[section_counter] = section_dict

	# 			del cur_resol_sectionimg
	# 			del cur_resol_alpha
	# 			del cur_resol_masked_sectionimg
	# 			del sectionimg_hsv

	# 	for x in curr_resol_slide_image:
	# 		del x

	# json.dump(summary_dict, open('/tmp/summary.json', 'w'), cls=NumpyAwareJSONEncoder)


# if __name__ == '__main__':

# 	parser = argparse.ArgumentParser()
# 	parser.add_argument("stack", type=str, help="choose what stack of images to crop and resolution, ex: RS141")
# 	args = parser.parse_args()
# 	stack = args.stack
	
# 	data_dir_stack = create_if_not_exists(os.path.join(os.environ['GORDON_SECTIONDATA_DIR'], stack))	# DavidData2015sections/CC35

# 	autogen_masked_img_dir = create_if_not_exists(os.path.join(data_dir_stack, 'autogen_maskedimg_x0.3125'))		# DavidData2015sections/CC35/autogen_masked_x0.3125
# 	autogen_bbox_dir = create_if_not_exists(os.path.join(data_dir_stack, 'autogen_bbox_x0.3125'))		# DavidData2015sections/CC35/autogen_bbox

# 	for l, resol in enumerate(['x0.3125', 'x1.25', 'x5', 'x20']):
# 		data_dir_resol = create_if_not_exists(os.path.join(data_dir_stack, resol))		# DavidData2015sections/CC35/x0.3125

# 	all_bboxfiles = sorted([bbox_file for bbox_file in os.listdir(autogen_bbox_dir) 
# 									if bbox_file.endswith('txt')])

# 	bf_sc_tuples = []

# 	os.chdir(autogen_bbox_dir)
# 	section_counter = 0
# 	for slide_ind, f in enumerate(all_bboxfiles):
# 		bboxes = np.loadtxt(f).astype(np.int)
# 		# existing_indices = np.where((bboxes == -1).all(axis=1).flat)[0]
# 		bf_sc_tuples.append((f, section_counter))
# 		section_counter += bboxes.shape[0]

# 	print bf_sc_tuples

# 	Parallel(n_jobs=16)(delayed(refine_save_highres_maskedimgs)(b,i, autogen_bbox_dir) for b, i in bf_sc_tuples)