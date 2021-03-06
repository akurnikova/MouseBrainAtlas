import sys
import os

import numpy as np
from skimage.io import imread, imsave
from matplotlib import pyplot as plt
from skimage.color import rgb2gray, rgb2hsv
from skimage.restoration import denoise_bilateral
from skimage.util import img_as_float, img_as_ubyte
from skimage.measure import label, regionprops, find_contours
from skimage.morphology import disk, remove_small_objects, binary_dilation, square
from skimage.filter.rank import median
from skimage.filter import gaussian_filter

from PIL import Image

def gen_mask(slide_ind=None, in_dir=None, stack=None, stain='nissl', 
		out_dir=None, num_section_per_slide=5, mirror=None, rotate=None, show=False, resol='x0.3125'):

	scaling = 4**['x0.3125', 'x1.25', 'x5'].index(resol)

	# scaling = 1

	try:
		img_rgb = np.array(Image.open(os.path.join(in_dir, '%s_%02d_%s_z0.tif'%(stack, slide_ind, resol))).convert('RGB'))
	except IOError as e:
		print e
		return

	img_gray = rgb2gray(img_rgb)

	if stain == 'nissl':
		imghsv = rgb2hsv(img_rgb)
		img = imghsv[..., 0]
	# imghue = imghsv[...,0]
	elif stain == 'fluorescent':
		img = img_gray

	# plt.imshow(imghsv[...,2], cmap=plt.cm.gray)
	# plt.colorbar()
	# plt.show()

	# plt.imshow(img_gray, cmap=plt.cm.gray)
	# plt.colorbar()
	# plt.show()

	# detect scanned zones on whole-slide images

	if stain == 'nissl':
		scanzone_mask = img_gray < 0.98
	elif stain == 'fluorescent':
		scanzone_mask = img_gray > 0.0001

	scanzone_mask = median(scanzone_mask.astype(np.float), disk(3*scaling))
	scanzone_mask = remove_small_objects(scanzone_mask.astype(np.bool), 
									min_size=50*scaling, connectivity=2, in_place=False)
									# min_size=2000*scaling, connectivity=2, in_place=False)

	if show:
		plt.imshow(scanzone_mask, cmap=plt.cm.gray)
		plt.title('detect scanned zones')
		plt.show()

	column_labels, n_columns = label(scanzone_mask, neighbors=8, return_num=True, background=0)
	column_props = regionprops(column_labels + 1)

	print len(column_props), 'scanned zones detected'

	mask = np.zeros_like(img_gray, dtype=np.bool)
	
	section_masks_all = []
	section_bboxes_all = []
	centers_slide_all = []

	for i, cp in enumerate(column_props):

		minr, minc, maxr, maxc = cp.bbox

		sys.stderr.write('start slide %d column %d\n' % (slide_ind, i))

		section_masks, section_bboxes = foreground_mask_simple(img[minr:maxr, minc:maxc],
																show=show, scaling=scaling)

		sys.stderr.write('finish slide %d column %d\n' % (slide_ind, i))
		
		if len(section_masks) > 0:
			section_masks_all += section_masks
			bboxes_slide = section_bboxes + np.r_[minr, minc, minr, minc][np.newaxis, :]

			section_bboxes_all += list(bboxes_slide)

			centers_slide = np.atleast_2d(np.column_stack([bboxes_slide[:,[1,3]].mean(axis=1),
												bboxes_slide[:,[0,2]].mean(axis=1)]))
			
			centers_slide_all += list(centers_slide)

	numbering, n_slots = renumber_blobs(np.array(centers_slide_all), img.shape)

	# if len(section_bboxes_all) != len(numbering):
	# 	raise Exception('Error: number of section boxes does not match length of numbering: slide %d' % (slide_ind))

	if len(section_bboxes_all) != len(numbering):
		print 'Error: number of section boxes does not match length of numbering: slide %d' % (slide_ind)
		return []

	for (minr, minc, maxr, maxc), mask, section_id in zip(section_bboxes_all, section_masks_all, numbering):

		img = img_rgb[minr:maxr, minc:maxc]

		masked_img = np.dstack([img, img_as_ubyte(mask)])

		if show:
			plt.imshow(masked_img)
			plt.title(section_id)
			plt.show()

		if out_dir is not None:
			image_filepath = os.path.join(out_dir, '%s_%s_%02d_%d.tif'%(stack, resol, slide_ind, section_id))
			imsave(image_filepath, masked_img)

	filled_bboxes = -1 * np.ones((n_slots, 4), dtype=np.int)
	filled_bboxes[numbering, :] = section_bboxes_all

	return filled_bboxes


def renumber_blobs(centers, img_shape, num_section_per_slide=5):
	"""
	centers n x (x,y)
	"""

	n = centers.shape[0]
	sorted_indices = np.argsort(centers[:,0])

	sorted_centers = centers[sorted_indices]
	height, width = img_shape

	if n > num_section_per_slide:
		indices_close = np.where(np.diff(sorted_centers[:,0]) < width * .1)[0]

		ups = []
		downs = []

		if len(indices_close) > 0:

			for i in range(n):
				if i - 1 in indices_close:
					continue
				elif i in indices_close:
					if sorted_centers[i,1] > height * .5:
						ups.append(sorted_indices[i + 1])
						downs.append(sorted_indices[i])
					else:
						ups.append(sorted_indices[i])
						downs.append(sorted_indices[i + 1])
				else:
					if sorted_centers[i,1] > height * .5:
						ups.append(99)
						downs.append(sorted_indices[i])
					else:
						ups.append(sorted_indices[i])
						downs.append(99)

		x = np.r_[ups, downs]
		numbering = np.argsort(x)[:n]

		return numbering, len(x)

	elif n < num_section_per_slide:

		# snap_to_columns = (np.round((sorted_centers[:,0] / width + 0.1) * num_section_per_slide) - 1).astype(np.int)
		# arrangement = -1 * np.ones((num_section_per_slide,), dtype=np.int)
		# arrangement[snap_to_columns] = sorted_indices
		# return arrangement[arrangement > -1], n
		return np.argsort(sorted_indices), n

	else:
		return np.argsort(sorted_indices), num_section_per_slide

# def foreground_mask_simple(img, show=False, hue=None, scaling=1):
def foreground_mask_simple(img, show=False, scaling=1):

	if show:
		plt.imshow(img, cmap=plt.cm.gray)
		plt.title('original')
		plt.show()

	img = denoise_bilateral(img, win_size=5*scaling)
	# img = gaussian_filter(img, sigma=.5)

	if img.shape[0] < 30*scaling or img.shape[1] < 30*scaling:
		return [], []

	# print img.shape

	bg_samples = [img[5,5],img[-5,-5],img[5,-5],img[-5,5],
					img[10,10],img[-10,-10],img[10,-10],img[-10,10]]

	bg = np.median(bg_samples)
	std = np.std(bg_samples)

	# if hue is not None:
	# 	bg_hue_samples = [hue[5,5],hue[-5,-5],hue[5,-5],hue[-5,5],
	# 					hue[10,10],hue[-10,-10],hue[10,-10],hue[-10,10]]
	# 	bg_hue = np.median(bg_hue_samples)
	# 	hue_std = np.std(bg_hue_samples)

	if show:
		plt.imshow(img, cmap=plt.cm.gray)
		plt.title('after denoising')
		plt.show()

	# if hue is not None:
	# 	img = (img < bg + 2 * 0.01) | (np.abs(hue - bg_hue) < 0.1)
	# else:
	img = img < bg + 2 * 0.01  # std computed from sample maybe too small, that is a risk

	if show:
		plt.imshow(img, cmap=plt.cm.gray)
		plt.title('thresh')
		plt.show()

	from skimage.morphology import binary_opening, disk, binary_closing

	img = binary_closing(img, selem=disk(1*scaling)) > 0

	if show:
		plt.imshow(img, cmap=plt.cm.gray)
		plt.title('closing')
		plt.show()

	# img = remove_small_objects(img, min_size=img.size * .001, connectivity=8, in_place=False)	

	# if show:
	# 	plt.imshow(img, cmap=plt.cm.gray)
	# 	plt.show()

	bg_region_labels, n_bg_regions = label(img, neighbors=8, return_num=True, background=0)

	bg_regions = regionprops(bg_region_labels + 1)
	
	border_indices = [i for i, r in enumerate(bg_regions) if (np.array(r.bbox, dtype=np.int) == np.r_[0, 0, img.shape[0], img.shape[1]]).any()]

	fg_mask = np.ones_like(img, dtype=np.bool)
	for i in border_indices:
		fg_mask[bg_region_labels==i] = 0

	if show:
		plt.imshow(fg_mask, cmap=plt.cm.gray)
		plt.title('fg_mask')
		plt.show()

	fg_mask = binary_closing(fg_mask, selem=disk(3*scaling)) > 0

	if show:
		plt.imshow(fg_mask, cmap=plt.cm.gray)
		plt.title('opening')
		plt.show()

	fg_mask_save = fg_mask.copy()

	fg_mask = remove_small_objects(fg_mask, min_size=int(img.size*.01), connectivity=8, in_place=True)

	if show:
		plt.imshow(fg_mask, cmap=plt.cm.gray)
		plt.title('remove small')
		plt.show()

	fg_region_labels, n_fg_regions = label(fg_mask, neighbors=8, return_num=True, background=0)
	fg_regions = regionprops(fg_region_labels + 1)

	if len(fg_regions) == 0:
		return [], []

	max_area = np.max([r.area for r in fg_regions])

	fg_mask = remove_small_objects(fg_mask_save, min_size=max_area*.5, connectivity=8, in_place=True)

	if show:
		plt.imshow(fg_mask, cmap=plt.cm.gray)
		plt.title('remove small given largest piece')
		plt.show()

	fg_region_labels, n_fg_regions = label(fg_mask, neighbors=8, return_num=True, background=0)
	fg_regions = regionprops(fg_region_labels + 1)

	masks = [r.image for r in fg_regions]

	bboxes = np.atleast_2d([r.bbox for r in fg_regions])

	return masks, bboxes


if __name__ == '__main__':

	# for i in range(23,48):
	# 	print i
	filled_bboxes = gen_mask(slide_ind=3, in_dir='/home/yuncong/DavidData2015slides/CC34/x0.3125', 
		# stack='CC28', stain='nissl', 
		stack='CC34', stain='fluorescent', 
		out_dir=None, num_section_per_slide=5, show=True)

	print filled_bboxes