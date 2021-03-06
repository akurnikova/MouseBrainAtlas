import sys
import os
import time

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from shapely.geometry import Polygon
from pandas import read_hdf, DataFrame, read_csv
from itertools import groupby
from collections import defaultdict
try:
    import mxnet as mx
except:
    sys.stderr.write("Cannot import mxnet.\n")

sys.path.append(os.environ['REPO_DIR'] + '/utilities')
from utilities2015 import *
from metadata import *
from data_manager import *
from visualization_utilities import *
from annotation_utilities import *

def plot_roc_curve(fp_allthresh, tp_allthresh, optimal_th, title=''):
    """
    Plot ROC curve.
    """

    plt.plot([fp_allthresh[th] for th in np.arange(0, 1, 0.01)],
             [tp_allthresh[th] for th in np.arange(0, 1, 0.01)]);

    plt.scatter(fp_allthresh[optimal_th], tp_allthresh[optimal_th], 
                marker='o', facecolors='none', edgecolors='k')
    
    plt.plot(np.arange(0, 1, 0.01), np.arange(0, 1, 0.01), c='k', linestyle='--');
    
    plt.legend();
    plt.axis('equal');
    plt.ylabel('True positive rate');
    plt.xlabel('False positive rate');
    plt.title(title);
    plt.show();
    
def plot_pr_curve(precision_allthresh, recall_allthresh, optimal_th, title=''):
    """
    Plot precision-recall curve. X axis is recall and y axis is precision.
    """

    plt.plot([recall_allthresh[th] for th in np.arange(0, 1, 0.01)],
             [precision_allthresh[th] for th in np.arange(0, 1, 0.01)]);

    plt.scatter(recall_allthresh[optimal_th], precision_allthresh[optimal_th], 
                marker='o', facecolors='none', edgecolors='k')
        
    plt.legend();
    plt.axis('equal');
    plt.ylabel('Precision');
    plt.xlabel('Recall');
    plt.title(title);
    plt.show();

def load_mxnet_model(model_dir_name, model_name, num_gpus=8, batch_size = 256):
    download_from_s3(os.path.join(MXNET_MODEL_ROOTDIR, model_dir_name), is_dir=True)
    model_iteration = 0
    output_symbol_name = 'flatten_output'
    output_dim = 1024
    mean_img = np.load(os.path.join(MXNET_MODEL_ROOTDIR, model_dir_name, 'mean_224.npy'))

    # Reference on how to predict with mxnet model:
    # https://github.com/dmlc/mxnet-notebooks/blob/master/python/how_to/predict.ipynb
    model0, arg_params, aux_params = mx.model.load_checkpoint(os.path.join(MXNET_MODEL_ROOTDIR, model_dir_name, model_name), 0)
    flatten_output = model0.get_internals()[output_symbol_name]
    # if HOST_ID == 'workstation':
    # model = mx.mod.Module(context=[mx.gpu(i) for i in range(1)], symbol=flatten_output)
    # else:
    model = mx.mod.Module(context=[mx.gpu(i) for i in range(num_gpus)], symbol=flatten_output)

    # Increase batch_size to 500 does not save any time.
    model.bind(data_shapes=[('data', (batch_size,1,224,224))], for_training=False)
    model.set_params(arg_params=arg_params, aux_params=aux_params, allow_missing=True)
    return model, mean_img

def convert_image_patches_to_features(patches, model, mean_img, batch_size):
    """
    Args:
        patches (list of (224,224)-ndarrays of uint8)
        model (mxnet.Module)
        mean_img ((224,224)-ndarray of uint8) : the mean image
        batch_size (int): batch size, must match the batch size specified for the model.
    """
    patches_mean_subtracted = patches - mean_img
    patches_mean_subtracted_input = patches_mean_subtracted[:, None, :, :] # n x 1 x 224 x 224

    if len(patches) < batch_size:
        n_pad = batch_size - len(patches)
        patches_mean_subtracted_input = np.concatenate([patches_mean_subtracted_input, np.zeros((n_pad,1) + mean_img.shape, mean_img.dtype)])
    
    print patches_mean_subtracted_input.shape
    data_iter = mx.io.NDArrayIter(
                    patches_mean_subtracted_input, 
                    batch_size=batch_size,
                    shuffle=False)
    outputs = model.predict(data_iter, always_output_list=True)
    features = outputs[0].asnumpy()
    
    if len(patches) < batch_size:
        features = features[:len(patches)]

    return features

def rotate_all_patches_variant(patches, variant):
    """
    Args:
        Variant (int): one of 0 to 7.
    """

    if variant == 0:
        patches_variant = patches
    elif variant == 1:
        patches_variant = [p.T[::-1] for p in patches]
    elif variant == 2:
        patches_variant = [p[::-1] for p in patches]
    elif variant == 3:
        patches_variant = [p[:, ::-1] for p in patches]
    elif variant == 4:
        patches_variant = [p[::-1, ::-1] for p in patches]
    elif variant == 5:
        patches_variant = [p.T for p in patches]
    elif variant == 6:
        patches_variant = [p.T[::-1, ::-1] for p in patches]
    elif variant == 7:
        patches_variant = [p.T[:, ::-1] for p in patches]
    return patches_variant

def rotate_all_patches(patches_enlarged, r, output_size=224):
    """
    Args:
        patches_enlarged: 
        r (int): rotation angle in degrees
        output_size (int): size of output patches
    """
    
    half_size = output_size/2
    patches_rotated = img_as_ubyte(np.array([rotate(p, angle=r)[p.shape[1]/2-half_size:p.shape[1]/2+half_size, 
                                                                p.shape[0]/2-half_size:p.shape[0]/2+half_size] 
                                             for p in patches_enlarged]))
    return patches_rotated


def load_dataset_addresses(dataset_ids, labels_to_sample=None, clf_rootdir=CLF_ROOTDIR):
    merged_addresses = {}

    for dataset_id in dataset_ids:
        addresses_fp = DataManager.get_dataset_addresses_filepath(dataset_id=dataset_id)
        download_from_s3(addresses_fp)
        addresses_curr_dataset = load_pickle(addresses_fp)
        
        if labels_to_sample is None:
            labels_to_sample = addresses_curr_dataset.keys()

        for label in labels_to_sample:
            try:
                if label not in merged_addresses:
                    merged_addresses[label] = addresses_curr_dataset[label]
                else:
                    merged_addresses[label] += addresses_curr_dataset[label]

            except Exception as e:
                continue

    return merged_addresses

def load_dataset_images(dataset_ids, labels_to_sample=None, clf_rootdir=CLF_ROOTDIR):
    
    merged_patches = {}
    merged_addresses = {}

    for dataset_id in dataset_ids:

        # load training addresses
        
#         addresses_fp = DataManager.get_dataset_addresses_filepath(dataset_id=dataset_id)
#         download_from_s3(addresses_fp)
#         addresses_curr_dataset = load_pickle(addresses_fp)
        
        # Load training features
        
        import re
        if labels_to_sample is None:
            labels_to_sample = []
            for dataset_id in dataset_ids:
                dataset_dir = DataManager.get_dataset_dir(dataset_id=dataset_id)
                download_from_s3(dataset_dir, is_dir=True)
                for fn in os.listdir(dataset_dir):
                    g = re.match('patch_images_(.*).hdf', fn).groups()
                    if len(g) > 0:
                        labels_to_sample.append(g[0])
        
        for label in labels_to_sample:
            try:
                # patches_curr_dataset_label_fp = os.path.join(CLF_ROOTDIR, 'datasets', 'dataset_%d' % dataset_id, 'patch_images_%s.hdf' % label)
                patches_curr_dataset_label_fp = DataManager.get_dataset_patches_filepath(dataset_id=dataset_id, structure=label)
                download_from_s3(patches_curr_dataset_label_fp)
                patches = bp.unpack_ndarray_file(patches_curr_dataset_label_fp)
                
                if label not in merged_patches:
                    merged_patches[label] = patches
                else:
                    merged_patches[label] = np.vstack([merged_patches[label], patches])
                
                addresses_fp = DataManager.get_dataset_addresses_filepath(dataset_id=dataset_id, structure=label)
                download_from_s3(addresses_fp)
                addresses = load_pickle(addresses_fp)
                    
                if label not in merged_addresses:
                    merged_addresses[label] = addresses
                else:
                    merged_addresses[label] += addresses

            except Exception as e:
                # sys.stderr.write("Cannot load dataset %d images for label %s: %s\n" % (dataset_id, label, str(e)))
                continue
                                
    return merged_patches, merged_addresses


def load_datasets_bp(dataset_ids, labels_to_sample=None, clf_rootdir=CLF_ROOTDIR):
    
    merged_features = {}
    merged_addresses = {}

    for dataset_id in dataset_ids:

        if labels_to_sample is None:
            import re
            labels_to_sample = []
            for dataset_id in dataset_ids:
                dataset_dir = DataManager.get_dataset_dir(dataset_id=dataset_id)
                #download_from_s3(dataset_dir, is_dir=True)
                for fn in os.listdir(dataset_dir):
                    g = re.match('patch_features_(.*).bp', fn).groups()
                    if len(g) > 0:
                        labels_to_sample.append(g[0])
        
        for label in labels_to_sample:
            try:
                # Load training features

                features_fp = DataManager.get_dataset_features_filepath(dataset_id=dataset_id, structure=label)
                #download_from_s3(features_fp)
                features = bp.unpack_ndarray_file(features_fp)
            
                # load training addresses

                addresses_fp = DataManager.get_dataset_addresses_filepath(dataset_id=dataset_id, structure=label)
                #download_from_s3(addresses_fp)
                addresses = load_pickle(addresses_fp)

                if label not in merged_features:
                    merged_features[label] = features
                else:
                    merged_features[label] = np.concatenate([merged_features[label], features])

                if label not in merged_addresses:
                    merged_addresses[label] = addresses
                else:
                    merged_addresses[label] += addresses

            except Exception as e:
                continue
                                
    return merged_features, merged_addresses

def load_datasets(dataset_ids, labels_to_sample=None, clf_rootdir=CLF_ROOTDIR, ext='hdf'):
    
    merged_features = {}
    merged_addresses = {}

    for dataset_id in dataset_ids:

        # load training addresses

        addresses_fp = DataManager.get_dataset_addresses_filepath(dataset_id=dataset_id)
        download_from_s3(addresses_fp)
        addresses_curr_dataset = load_pickle(addresses_fp)
        
        # Load training features

        features_fp = DataManager.get_dataset_features_filepath(dataset_id=dataset_id, ext=ext)
        download_from_s3(features_fp)
        features_curr_dataset = load_hdf_v2(features_fp).to_dict()
        
        if labels_to_sample is None:
            labels_to_sample = features_curr_dataset.keys()

        for label in labels_to_sample:
            try:
                if label not in merged_features:
                    merged_features[label] = features_curr_dataset[label]
                else:
                    merged_features[label] = np.vstack([merged_features[label], features_curr_dataset[label]])

                if label not in merged_addresses:
                    merged_addresses[label] = addresses_curr_dataset[label]
                else:
                    merged_addresses[label] += addresses_curr_dataset[label]

            except Exception as e:
                continue
                                
    return merged_features, merged_addresses


################
## Evaluation ##
################

def compute_accuracy(predictions, true_labels, exclude_abstained=True, abstain_label=-1):

    n = len(predictions)
    if exclude_abstained:
        predicted_indices = predictions != abstain_label
        acc = np.count_nonzero(predictions[predicted_indices] == true_labels[predicted_indices]) / float(len(predicted_indices))
    else:
        acc = np.count_nonzero(predictions == true_labels) / float(n)

    return acc

def compute_predictions(H, abstain_label=-1):
    predictions = np.argmax(H, axis=1)
    no_decision_indices = np.where(np.all(H == 0, axis=1))[0]

    if abstain_label == 'random':
        predictions[no_decision_indices] = np.random.randint(0, H.shape[1], len(no_decision_indices)).astype(np.int)
    else:
        predictions[no_decision_indices] = abstain_label

    return predictions, no_decision_indices

def compute_confusion_matrix(probs, labels, soft=False, normalize=True, abstain_label=-1, binary=True, decision_thresh=.5):
    """
    Args:
        probs ((n_example, n_class)-ndarray of float): prediction probabilities
        label ((n_example,)-ndarray of int): true example labels
    """

    n_labels = len(np.unique(labels))

    pred_is_hard = np.array(probs).ndim == 1

    if pred_is_hard:
        soft = False

    M = np.zeros((n_labels, n_labels))
    for probs0, tl in zip(probs, labels):
        if soft:
            M[tl] += probs0
        else:
            if pred_is_hard:
                if probs0 != abstain_label:
                    M[tl, probs0] += 1
            else:
                hard = np.zeros((n_labels, ))
                if binary:
                    if probs0[0] > decision_thresh:
                        hard[0] = 1.
                    else:
                        hard[1] = 1.
                else:
                    hard[np.argmax(probs0)] = 1.
                M[tl] += hard

    if normalize:
        M_normalized = M.astype(np.float)/M.sum(axis=1)[:, np.newaxis]
        return M_normalized
    else:
        return M

def plot_confusion_matrix(cm, labels, title='Confusion matrix', cmap=plt.cm.Blues, figsize=(4,4), text=True, axis=None,
xlabel='Predicted label', ylabel='True label', **kwargs):

    if axis is None:
        fig = plt.figure(figsize=figsize)
        axis = fig.add_subplot(1,1,1)
        return_fig = True
    else:
        return_fig = False

    axis.imshow(cm, interpolation='nearest', cmap=cmap, vmin=0, vmax=1, **kwargs)
    axis.set_title(title)
#     plt.colorbar()
    axis.set_xticks(np.arange(len(labels)))
    axis.set_yticks(np.arange(len(labels)))
    axis.set_xticklabels(labels)
    axis.set_yticklabels(labels)

    axis.set_ylabel(ylabel)
    axis.set_xlabel(xlabel)

    if cm.dtype.type is np.int_:
        fmt = '%d'
    else:
        fmt = '%.2f'

    if text:
        for x in xrange(len(labels)):
            for y in xrange(len(labels)):
                if not np.isnan(cm[y,x]):
                    axis.text(x,y, fmt % cm[y,x],
                             horizontalalignment='center',
                             verticalalignment='center');

    if return_fig:
        return fig

def locate_patches_given_addresses_v2(addresses):
    """
    addresses is a list of addresses.
    address: stack, section, framewise_index
    """

    from collections import defaultdict
    addresses_grouped = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for addressList_index, (stack, sec, i) in enumerate(addresses):
        addresses_grouped[stack][sec].append((i, addressList_index))

    patch_locations_all = []
    addressList_indices_all = []
    for stack, indices_allSections in addresses_grouped.iteritems():
        grid_spec = get_default_gridspec(stack)
        sample_locations = grid_parameters_to_sample_locations(grid_spec)
        for sec, fwInd_addrLstInd_tuples in indices_allSections.iteritems():
            frameWise_indices_selected, addressList_indices = map(list, zip(*fwInd_addrLstInd_tuples))
            patch_locations_all += [(stack, sec, loc, name) for loc in sample_locations[frameWise_indices_selected]]
            addressList_indices_all += addressList_indices

    patch_locations_all_inOriginalOrder = [patch_locations_all[i] for i in np.argsort(addressList_indices_all)]
    return patch_locations_all_inOriginalOrder

def extract_patches_given_locations(patch_size, locs, 
                                    img=None, 
                                    stack=None, sec=None, version=None, prep_id=2):
    """
    Extract patches from one image at given locations.
    The image can be given, or the user can provide stack,sec,version,prep_id.

    Args:
        img: the image
        locs ((n,2)-array): patch centers
        patch_size (int): size of a patch, assume square

    Returns:
        list of (patch_size, patch_size)-arrays.
    """

    if img is None:
        t = time.time()
        img = DataManager.load_image_v2(stack=stack, section=sec, prep_id=prep_id, version=version)
        sys.stderr.write('Load image: %.2f seconds.\n' % (time.time() - t)) 

    half_size = patch_size/2
    patches = [img[y-half_size:y+half_size, x-half_size:x+half_size].copy() for x, y in locs]
    return patches

def extract_patches_given_locations_multiple_sections(addresses, 
                                                      version=None, prep_id=2, 
                                                      win_id=None, patch_size=None, 
                                                      location_or_grid_index='location'):
    """
    Extract patches from multiple images.
    
    Args:
        addresses (list of tuples):
            if location_or_grid_index is 'location', then this is a list of (stack, section, location), patch_size is required.
            if location_or_grid_index is 'grid_index', then this is a list of (stack, section, framewise_index), win_id is required.
            
    Returns:
        list of (patch_size, patch_size)-arrays.
    """

    from collections import defaultdict

    locations_grouped = {}
    for stack_sec, list_index_and_address_grouper in groupby(sorted(enumerate(addresses), key=lambda (i, x): (x[0],x[1])),
        key=lambda (i,x): (x[0], x[1])):
        locations_grouped[stack_sec] = [(address[2], list_index) for list_index, address in list_index_and_address_grouper]

    patches_all = []
    list_indices_all = []
    for stack_sec, locations_allSections in locations_grouped.iteritems():
        stack, sec = stack_sec
        if location_or_grid_index == 'location':
            locs_thisSec, list_indices = map(list, zip(*locations_allSections))
        else:
            ginds_thisSec, list_indices = map(list, zip(*locations_allSections))
            assert win_id is not None
            patch_size = windowing_settings[win_id]['patch_size']
            locs_thisSec = grid_parameters_to_sample_locations(patch_size=patch_size, 
                                                       stride=windowing_settings[win_id]['spacing'],
                                                      w=metadata_cache['image_shape'][stack][0],
                                                       h=metadata_cache['image_shape'][stack][1])[ginds_thisSec]
        extracted_patches = extract_patches_given_locations(stack=stack, sec=sec, locs=locs_thisSec, version=version,
                                                               patch_size=patch_size)
        patches_all += extracted_patches
        list_indices_all += list_indices

    patch_all_inOriginalOrder = [patches_all[i] for i in np.argsort(list_indices_all)]
    return patch_all_inOriginalOrder

def get_names_given_locations_multiple_sections(addresses, location_or_grid_index='location', username=None):

    # from collections import defaultdict
    # locations_grouped = {stack_sec: (x[2], list_index) \
    #                     for stack_sec, (list_index, x) in group_by(enumerate(addresses), lambda i, x: (x[0],x[1]))}

    locations_grouped = defaultdict(lambda: defaultdict(list))
    for list_index, (stack, sec, loc) in enumerate(addresses):
        locations_grouped[stack][sec].append((loc, list_index))

    # if indices is not None:
    #     assert locations is None, 'Cannot specify both indices and locs.'
    #     if grid_locations is None:
    #         grid_locations = grid_parameters_to_sample_locations(grid_spec)
    #     locations = grid_locations[indices]

    names_all = []
    list_indices_all = []
    for stack, locations_allSections in locations_grouped.iteritems():

        structure_grid_indices = locate_annotated_patches(stack=stack, username=username, force=True,
                                                        annotation_rootdir=annotation_midbrainIncluded_v2_rootdir)

        for sec, loc_listInd_tuples in locations_allSections.iteritems():

            label_dict = structure_grid_indices[sec].dropna().to_dict()

            locs_thisSec, list_indices = map(list, zip(*loc_listInd_tuples))

            index_to_name_mapping = {}

            for name, full_indices in label_dict.iteritems():
                for i in full_indices:
                    index_to_name_mapping[i] = name

            if location_or_grid_index == 'location':
                raise Exception('Not implemented.')
            else:
                names = [index_to_name_mapping[i] if i in index_to_name_mapping else 'BackG' for i in locs_thisSec]

            names_all += names
            list_indices_all += list_indices

    names_all_inOriginalOrder = [names_all[i] for i in np.argsort(list_indices_all)]
    return names_all_inOriginalOrder


def get_default_gridspec(stack, patch_size=224, stride=56):
    # image_width, image_height = DataManager.get_image_dimension(stack)
    image_width, image_height = metadata_cache['image_shape'][stack]
    return (patch_size, stride, image_width, image_height)


def label_regions_multisections(stack, region_contours, surround_margins=None):
    """
    Args:
        region_contours: dict {sec: list of contours (nx2 array)}
    """

    # contours_df = read_hdf(ANNOTATION_ROOTDIR + '/%(stack)s/%(stack)s_annotation_v3.h5' % dict(stack=stack), 'contours')
    contours_df = DataManager.load_annotation_v4(stack=stack)
    contours = contours_df[(contours_df['orientation'] == 'sagittal') & (contours_df['downsample'] == 1)]
    contours = contours.drop_duplicates(subset=['section', 'name', 'side', 'filename', 'downsample', 'creator'])
    labeled_contours = convert_annotation_v3_original_to_aligned_cropped(contours, stack=stack)

    sections_to_filenames = metadata_cache['sections_to_filenames'][stack]

    labeled_region_indices_all_sections = {}

    for sec, region_contours_curr_sec in region_contours.iteritems():

        sys.stderr.write('Analyzing section %d..\n' % section)

        matching_contours = labeled_contours[labeled_contours['section'] == section]
        if len(matching_contours) == 0:
            continue
        polygons_this_sec = [(contour['name'], contour['vertices']) for contour_id, contour in matching_contours.iterrows()]
        mask_tb = DataManager.load_thumbnail_mask_v2(stack, section)
        labeled_region_indices = identify_regions_inside(region_contours=region_contours, mask_tb=mask_tb, polygons=polygons_this_sec, \
                                            surround_margins=surround_margins)

        labeled_region_indices_all_sections[sec] = labeled_region_indices

    return labeled_region_indices_all_sections


def label_regions(stack, section, region_contours, surround_margins=None, labeled_contours=None):
    """
    Identify regions (could be non-square) that belong to each class.
    
    Returns:
        dict of {label: region indices in input region_contours}
    """

    if is_invalid(sec=section, stack=stack):
        raise Exception('Section %s, %d invalid.' % (stack, section))

    if labeled_contours is None:
        # contours_df = read_hdf(ANNOTATION_ROOTDIR + '/%(stack)s/%(stack)s_annotation_v3.h5' % dict(stack=stack), 'contours')
        contours_df = DataManager.load_annotation_v4(stack=stack)
        contours = contours_df[(contours_df['orientation'] == 'sagittal') & (contours_df['downsample'] == 1)]
        contours = contours.drop_duplicates(subset=['section', 'name', 'side', 'filename', 'downsample', 'creator'])
        contours = convert_annotation_v3_original_to_aligned_cropped(contours, stack=stack)
        labeled_contours = contours[contours['section'] == section]

    sys.stderr.write('Analyzing section %d..\n' % section)

    if len(labeled_contours) == 0:
        return {}

    polygons_this_sec = [(contour['name'], contour['vertices']) for contour_id, contour in labeled_contours.iterrows()
                        if contour['name'] in all_known_structures]

    mask_tb = DataManager.load_thumbnail_mask_v2(stack, section)
    patch_indices = identify_regions_inside(region_contours=region_contours, mask_tb=mask_tb, polygons=polygons_this_sec, \
                                        surround_margins=surround_margins)
    return patch_indices

def identify_regions_inside(region_contours, stack=None, image_shape=None, mask_tb=None, polygons=None, surround_margins=None):
    """
    Return addresses of patches that are in polygons or on mask.
    - If mask is given, the valid patches are those whose centers are True. bbox and polygons are ANDed with mask.
    - If polygons is given, the valid patches are those whose bounding boxes.
        - polygons can be a dict, keys are structure names, values are x-y vertices (nx2 array).
        - polygons can also be a list of (name, vertices) tuples.
    if shrinked to 30%% are completely within the polygons.

    Args:
        surround_margins (list of int): list of surround margins in which patches are extracted, in unit of microns.
    """
    
    # rename for clarity
    surround_margins_um = surround_margins

    n_regions = len(region_contours)

    if isinstance(polygons, dict):
        polygon_list = [(name, cnt) for name, cnts in polygons.iteritems() for cnt in cnts] # This is to deal with when one name has multiple contours
    elif isinstance(polygons, list):
        assert isinstance(polygons[0], tuple)
        polygon_list = polygons
    else:
        raise Exception('Polygon must be either dict or list.')

    if surround_margins_um is None:
        surround_margins_um = [100,200,300,400,500,600,700,800,900,1000]

    # Identify patches in the foreground.
    region_centroids_tb = np.array([np.mean(cnt, axis=0)/32. for cnt in region_contours], np.int)
    indices_fg = np.where(mask_tb[region_centroids_tb[:,1], region_centroids_tb[:,0]])[0]
    # indices_fg = np.where(mask_tb[sample_locations[:,1]/32, sample_locations[:,0]/32])[0] # patches in the foreground
    # sys.stderr.write('%d patches in fg\n' % len(indices_fg))

    # Identify patches in the background.
    indices_bg = np.setdiff1d(range(n_regions), indices_fg)

    indices_inside = {}
    indices_allLandmarks = {}

    for label, poly in polygon_list:
        path = Path(poly)
        # indices_inside[label] = np.where([np.all(path.contains_points(cnt)) for cnt in region_contours])[0]

        # Being inside is defined as 70% of vertices are within the polygon
        indices_inside[label] = np.where([np.count_nonzero(path.contains_points(cnt)) >= len(cnt)*.7 for cnt in region_contours])[0]
        indices_allLandmarks[label] = indices_inside[label]
        # sys.stderr.write('%d patches in %s\n' % (len(indices_allLandmarks[label]), label))

    indices_allInside = np.concatenate(indices_inside.values())

    for label, poly in polygon_list:

        for margin_um in surround_margins_um:

            margin = margin_um / XY_PIXEL_DISTANCE_LOSSLESS
        
            surround = Polygon(poly).buffer(margin, resolution=2)

            path = Path(list(surround.exterior.coords))
            indices_sur = np.where([np.count_nonzero(path.contains_points(cnt)) >=  len(cnt)*.7  for cnt in region_contours])[0]

            # surround classes do not include patches of any no-surround class
            indices_allLandmarks[convert_to_surround_name(label, margin=margin_um, suffix='noclass')] = np.setdiff1d(indices_sur, np.r_[indices_bg, indices_allInside])

            for l, inds in indices_inside.iteritems():
                if l == label: continue
                indices = np.intersect1d(indices_sur, inds)
                if len(indices) > 0:
                    indices_allLandmarks[convert_to_surround_name(label, margin=margin_um, suffix=l)] = indices

        # Identify all foreground patches except the particular label's inside patches
        indices_allLandmarks[label+'_negative'] = np.setdiff1d(range(n_regions), np.r_[indices_bg, indices_inside[label]])
        # sys.stderr.write('%d patches in %s\n' % (len(indices_allLandmarks[label+'_negative']), label+'_negative'))

    indices_allLandmarks['bg'] = indices_bg
    indices_allLandmarks['noclass'] = np.setdiff1d(range(n_regions), np.r_[indices_bg, indices_allInside])

    return indices_allLandmarks


def locate_annotated_patches_v2(stack, grid_spec=None, sections=None, surround_margins=None):
    """
    Read in a structure annotation file generated by the labeling GUI. 
    Compute that each class label contains which grid indices.
    
    Returns:
        a DataFrame. Columns are structure names. Rows are section numbers. Cell is a 1D array of grid indices.
    """

    if grid_spec is None:
        grid_spec = get_default_gridspec(stack)

    contours_df = DataManager.load_annotation_v4(stack, annotation_rootdir=ANNOTATION_ROOTDIR)
    contours = contours_df[(contours_df['orientation'] == 'sagittal') & (contours_df['downsample'] == 1)]
    contours = contours.drop_duplicates(subset=['section', 'name', 'side', 'filename', 'downsample', 'creator'])
    contours_df = convert_annotation_v3_original_to_aligned_cropped(contours, stack=stack)
    
    contours_grouped = contours_df.groupby('section')

    patch_indices_allSections_allStructures = {}
    for sec, cnt_group in contours_grouped:
        if sections is not None and sec not in sections:
            continue
        sys.stderr.write('Analyzing section %d..\n' % sec)
        if is_invalid(sec=sec, stack=stack):
            continue
        polygons_this_sec = [(contour['name'], contour['vertices']) for contour_id, contour in cnt_group.iterrows()]
        mask_tb = DataManager.load_thumbnail_mask_v2(stack, sec)
        patch_indices_allSections_allStructures[sec] = \
        locate_patches_v2(grid_spec=grid_spec, mask_tb=mask_tb, polygons=polygons_this_sec, \
                                            surround_margins=surround_margins)
        
    return DataFrame(patch_indices_allSections_allStructures)


def generate_annotation_to_grid_indices_lookup(stack, by_human, win_id,
                                                   stack_m=None, 
                                                    classifier_setting_m=None, 
                                                    classifier_setting_f=None, 
                                                    warp_setting=None, trial_idx=None,
                                              surround_margins=None, suffix=None, timestamp=None, prep_id=2,
                                              ):
    """
    Load the structure annotation.
    Use the default grid spec.
    Find grid indices for each class label.
    
    Args:
        win_id (int): the spatial sample scheme
    
    Returns:
        DataFrame: Columns are class labels and rows are section indices.
    """
    
    windowing_properties = windowing_settings[win_id]
    patch_size = windowing_properties['patch_size']
    spacing = windowing_properties['spacing']
    w, h = metadata_cache['image_shape'][stack]
    half_size = patch_size/2
    grid_spec = (patch_size, spacing, w, h)
    
    contours_df = DataManager.load_annotation_v4(stack=stack, by_human=by_human, 
                                                      stack_m=stack_m, 
                                                           classifier_setting_m=classifier_setting_m,
                                                          classifier_setting_f=classifier_setting_f,
                                                          warp_setting=warp_setting,
                                                          trial_idx=trial_idx, suffix=suffix, timestamp=timestamp)
    contours = contours_df[(contours_df['orientation'] == 'sagittal') & (contours_df['downsample'] == 1)]
    contours = contours.drop_duplicates(subset=['section', 'name', 'side', 'filename', 'downsample', 'creator'])
    
    if by_human:    
        contours_df = convert_annotation_v3_original_to_aligned_cropped(contours, stack=stack)
    
    download_from_s3(DataManager.get_thumbnail_mask_dir_v3(stack=stack, prep_id=prep_id), is_dir=True)
    
    contours_grouped = contours_df.groupby('section')

    patch_indices_allSections_allStructures = {}
    for sec, cnt_group in contours_grouped:
        sys.stderr.write('Computing grid indices lookup for section %d...\n' % sec)
        if is_invalid(sec=sec, stack=stack):
            continue
        polygons_this_sec = [(contour['name'], contour['vertices']) for contour_id, contour in cnt_group.iterrows()]
        mask_tb = DataManager.load_thumbnail_mask_v3(stack=stack, section=sec, prep_id=prep_id)
        patch_indices_allSections_allStructures[sec] = \
        locate_patches_v2(grid_spec=grid_spec, mask_tb=mask_tb, polygons=polygons_this_sec, \
                          surround_margins=surround_margins)
    
    return DataFrame(patch_indices_allSections_allStructures).T
    
#     def locate_patches_worker(sec):
#         sys.stderr.write('Computing grid indices lookup for section %d...\n' % sec)
#         mask_tb = DataManager.load_thumbnail_mask_v3(stack=stack, section=sec)
#         contours = contours_df.loc[sec]
#         return locate_patches_v2(stack=stack, mask_tb=mask_tb, polygons=contours.to_dict(), 
#                                       surround_margins=surround_margins)
#     pool = Pool(NUM_CORES)
#     patch_indices_allSections_allStructures = dict(zip(contours_df.index, pool.map(locate_patches_worker, contours_df.index)))
#     pool.close()
#     pool.join()
    
    # patch_indices_allSections_allStructures = {}
    # for sec, contours in contours_df.iterrows():
    #     sys.stderr.write('Computing grid indices lookup for section %d...\n' % sec)
    #     mask_tb = DataManager.load_thumbnail_mask_v3(stack=stack, section=sec)
    #     patch_indices_allSections_allStructures[sec] = \
    #     locate_patches_v2(stack=stack, mask_tb=mask_tb, polygons=contours.dropna().to_dict(), 
    #                       surround_margins=surround_margins)
    
    # return DataFrame(patch_indices_allSections_allStructures).T
    
    
def sample_locations(grid_indices_lookup, labels, num_samples_per_polygon=None, num_samples_per_landmark=None):
    """
    !!! This should be sped up! It is now taking 100 seconds for one stack !!!
    
    Args:
        grid_indices_lookup: Rows are sections. Columns are class labels.
    
    Returns:
        address_list: list of (section, grid_idx) tuples.
    """
    
    location_list = defaultdict(list)
    for label in set(labels) & set(grid_indices_lookup.columns):
        for sec, grid_indices in grid_indices_lookup[label].dropna().iteritems():
            n = len(grid_indices)
            if n > 0:
                if num_samples_per_polygon is None:
                    location_list[label] += [(sec, i) for i in grid_indices]
                else:
                    random_sampled_indices = grid_indices[np.random.choice(range(n), min(n, num_samples_per_polygon), replace=False)]
                    location_list[label] += [(sec, i) for i in random_sampled_indices]

    if num_samples_per_landmark is not None:
        sampled_location_list = {}
        for name_s, addresses in location_list.iteritems():
            n = len(addresses)
            random_sampled_indices = np.random.choice(range(n), min(n, num_samples_per_landmark), replace=False)
            sampled_location_list[name_s] = [addresses[i] for i in random_sampled_indices]
        return sampled_location_list
    else:
        location_list.default_factory = None
        return location_list

def locate_patches_v2(grid_spec=None, stack=None, patch_size=None, stride=None, image_shape=None, 
                      mask_tb=None, polygons=None, bbox=None, bbox_lossless=None, surround_margins=None):
    """
    Return addresses of patches that are either in polygons or on mask.
    - If mask is given, the valid patches are those whose centers are True. bbox and polygons are ANDed with mask.
    - If bbox is given, valid patches are those entirely inside bbox. bbox = (x,y,w,h) in thumbnail resol.
    - If polygons is given, the valid patches are those whose bounding boxes.
        - polygons can be a dict, keys are structure names, values are x-y vertices (nx2 array).
        - polygons can also be a list of (name, vertices) tuples.
    if shrinked to 30%% are completely within the polygons.
            
    Args:
        grid_spec: If none, use the default grid spec.
        surround_margins: list of surround margin for which patches are extracted, in unit of microns. Default: [100,200,...1000]
    Returns:
        If polygons are given, returns dict {label: list of grid indices}.
        Otherwise, return a list of grid indices.
    """
    
    # rename for clarity
    surround_margins_um = surround_margins

    if grid_spec is None:
        grid_spec = get_default_gridspec(stack)

    if image_shape is not None:
        image_width, image_height = image_shape
    else:
        image_width, image_height = grid_spec[2:]
        # If use the following, stack might not be specified.
        # image_width, image_height = metadata_cache['image_shape'][stack]

    if patch_size is None:
        patch_size = grid_spec[0]

    if stride is None:
        stride = grid_spec[1]

    if surround_margins_um is None:
        surround_margins_um = [100,200,300,400,500,600,700,800,900,1000]

    sample_locations = grid_parameters_to_sample_locations(patch_size=patch_size, stride=stride, w=image_width, h=image_height)
    half_size = patch_size/2

    indices_fg = np.where(mask_tb[sample_locations[:,1]/32, sample_locations[:,0]/32])[0] # patches in the foreground
    indices_bg = np.setdiff1d(range(sample_locations.shape[0]), indices_fg) # patches in the background

    if polygons is None and bbox is None and bbox_lossless is None:
        return indices_fg

    if polygons is not None:
        if isinstance(polygons, dict):
            polygon_list = []
            for name, cnts in polygons.iteritems() :
                if np.asarray(cnts).ndim == 3: 
                    # If argument polygons is dict {label: list of nx2 arrays}. 
                    # This is the case that each label has multiple contours.
                    for cnt in cnts:
                        polygon_list.append((name, cnt))
                elif np.asarray(cnts).ndim == 2:
                    # If argument polygon is dict {label: nx2 array}
                    # This is the case that each label has exactly one contour.
                    polygon_list.append((name, cnts))
        elif isinstance(polygons, list):
            assert isinstance(polygons[0], tuple)
            polygon_list = polygons
        else:
            raise Exception('Polygon must be either dict or list.')

    if bbox is not None or bbox_lossless is not None:
        # Return grid indices in a bounding box
        assert polygons is None, 'Can only supply one of bbox or polygons.'

        if bbox is not None:

            box_x, box_y, box_w, box_h = bbox

            xmin = max(half_size, box_x*32)
            xmax = min(image_width-half_size-1, (box_x+box_w)*32)
            ymin = max(half_size, box_y*32)
            ymax = min(image_height-half_size-1, (box_y+box_h)*32)

        else:
            box_x, box_y, box_w, box_h = bbox_lossless
            # print box_x, box_y, box_w, box_h

            xmin = max(half_size, box_x)
            xmax = min(image_width-half_size-1, box_x+box_w-1)
            ymin = max(half_size, box_y)
            ymax = min(image_height-half_size-1, box_y+box_h-1)

        indices_roi = np.where(np.all(np.c_[sample_locations[:,0] > xmin, sample_locations[:,1] > ymin,
                                            sample_locations[:,0] < xmax, sample_locations[:,1] < ymax], axis=1))[0]
        indices_roi = np.setdiff1d(indices_roi, indices_bg)
        
        return indices_roi

    else:
        # Return grid indices in each polygon of the input dict
        assert polygons is not None, 'Can only supply one of bbox or polygons.'

        # This means we require a patch to have 30% of its radius to be within the landmark boundary to be considered inside the landmark
        tolerance_margin = int(.3*half_size)

        sample_locations_ul = sample_locations - (tolerance_margin, tolerance_margin)
        sample_locations_ur = sample_locations - (-tolerance_margin, tolerance_margin)
        sample_locations_ll = sample_locations - (tolerance_margin, -tolerance_margin)
        sample_locations_lr = sample_locations - (-tolerance_margin, -tolerance_margin)

        indices_inside = {}
        indices_allLandmarks = {}

        for label, poly in polygon_list:
            path = Path(poly)
            indices_inside[label] = np.where(path.contains_points(sample_locations_ll) &\
                                              path.contains_points(sample_locations_lr) &\
                                              path.contains_points(sample_locations_ul) &\
                                              path.contains_points(sample_locations_ur))[0]

            indices_allLandmarks[label] = indices_inside[label]

        indices_allInside = np.concatenate(indices_inside.values())

        for label, poly in polygon_list:

            for margin_um in surround_margins_um:
                margin = margin_um / XY_PIXEL_DISTANCE_LOSSLESS
                surround = Polygon(poly).buffer(margin, resolution=2)

                path = Path(list(surround.exterior.coords))
                indices_sur =  np.where(path.contains_points(sample_locations_ll) &\
                                        path.contains_points(sample_locations_lr) &\
                                        path.contains_points(sample_locations_ul) &\
                                        path.contains_points(sample_locations_ur))[0]

                # surround classes do not include patches of any no-surround class
                indices_allLandmarks[convert_to_surround_name(label, margin=margin_um, suffix='noclass')] = np.setdiff1d(indices_sur, np.r_[indices_bg, indices_allInside])
                
                for l, inds in indices_inside.iteritems():
                    if l == label: continue
                    indices = np.intersect1d(indices_sur, inds)
                    if len(indices) > 0:
                        indices_allLandmarks[convert_to_surround_name(label, margin=margin_um, suffix=l)] = indices

            # All foreground patches except the particular label's inside patches
            indices_allLandmarks[label+'_negative'] = np.setdiff1d(range(sample_locations.shape[0]), np.r_[indices_bg, indices_inside[label]])

        indices_allLandmarks['bg'] = indices_bg
        indices_allLandmarks['noclass'] = np.setdiff1d(range(sample_locations.shape[0]), np.r_[indices_bg, indices_allInside])

        return indices_allLandmarks

    
    
def generate_dataset_addresses(num_samples_per_label, stacks, labels_to_sample, grid_indices_lookup_fps):
    """
    Generate patch addresses grouped by label.
    
    Args:
        stacks (list of str):
        grid_indices_lookup_fps (dict of str):
    
    Returns:
        addresses
    """
    
    addresses = defaultdict(list)
    
    t = time.time()
    
    for stack in stacks:
        
        t1 = time.time()
        
        download_from_s3(grid_indices_lookup_fps[stack])
        if not os.path.exists(grid_indices_lookup_fps[stack]):
            raise Exception('Cannot find grid indices lookup. Use function generate_annotation_to_grid_indices_lookup to generate.')
        try:
            grid_indices_per_label = load_hdf_v2(grid_indices_lookup_fps[stack])
        except:
            grid_indices_per_label = read_hdf(grid_indices_lookup_fps[stack], 'grid_indices').T
                
        def convert_data_from_sided_to_unsided(data):
            """
            Args:
                data: rows are sections, columns are SIDED labels.
            """
            new_data = defaultdict(dict)
            for name_s in data.columns:
                name_u = convert_to_unsided_label(name_s)
                for sec, grid_indices in data[name_s].iteritems():
                    new_data[name_u][sec] = grid_indices
            return DataFrame(new_data)
        
        grid_indices_per_label = convert_data_from_sided_to_unsided(grid_indices_per_label)
        
        sys.stderr.write('Read: %.2f seconds\n' % (time.time() - t1))

        # Guarantee column is class name, row is section index.
        assert isinstance(grid_indices_per_label.columns[0], str) and isinstance(grid_indices_per_label.index[0], int), \
        "Must guarantee column is class name, row is section index."
        
        labels_this_stack = set(grid_indices_per_label.columns) & set(labels_to_sample)

        t1 = time.time()
        addresses_sec_idx = sample_locations(grid_indices_per_label, labels_this_stack, 
                                            num_samples_per_landmark=num_samples_per_label/len(stacks))
        sys.stderr.write('Sample addresses (stack %s): %.2s seconds.\n' % (stack, time.time() - t1))

        for label, addrs in addresses_sec_idx.iteritems():
            addresses[label] += [(stack, ) + addr for addr in addrs]

    addresses.default_factory = None
    
    sys.stderr.write('Sample addresses: %.2f seconds\n' % (time.time() - t))
        
    return addresses


def generate_dataset(num_samples_per_label, stacks, labels_to_sample, model_name, grid_indices_lookup_fps):
    """
    Generate a dataset using the following steps:
    - Load structure polygons
    - Sample addresses for each class (e.g. 7N, AP_surround_500_noclass, etc.)
    - Map addresses to conv net features
    - Remove None features
    
    Args:
        grid_indices_lookup_fps
    
    Returns:
        features, addresses
    """
    
    # Extract addresses
    
    addresses = defaultdict(list)

    t = time.time()
    
    for stack in stacks:
        
        t1 = time.time()
        download_from_s3(grid_indices_lookup_fps[stack])
        if not os.path.exists(grid_indices_lookup_fps[stack]):
            raise Exception('Cannot find grid indices lookup. Use function generate_annotation_to_grid_indices_lookup to generate.')
        grid_indices_per_label = load_hdf_v2(grid_indices_lookup_fps[stack])
        sys.stderr.write('Read grid indices lookup: %.2f seconds\n' % (time.time() - t1))
        
        def convert_data_from_sided_to_unsided(data):
            """
            Args:
                data: rows are sections, columns are SIDED labels.
            """
            new_data = defaultdict(dict)
            for name_s in data.columns:
                name_u = convert_to_unsided_label(name_s)
                for sec, grid_indices in data[name_s].iteritems():
                    new_data[name_u][sec] = grid_indices
            return DataFrame(new_data)
        
        grid_indices_per_label = convert_data_from_sided_to_unsided(grid_indices_per_label)
        labels_this_stack = set(grid_indices_per_label.columns) & set(labels_to_sample)
        
        t1 = time.time()
        addresses_sec_idx = sample_locations(grid_indices_per_label, labels_this_stack,
                                             num_samples_per_landmark=num_samples_per_label/len(stacks))
        sys.stderr.write('Sample addresses (stack %s): %.2f seconds.\n' % (stack, time.time() - t1))
        for label, addrs in addresses_sec_idx.iteritems():
            addresses[label] += [(stack, ) + addr for addr in addrs]

    addresses.default_factory = None
        
    sys.stderr.write('Sample addresses: %.2f seconds\n' % (time.time() - t))
    
    # Map addresses to features

    t = time.time()
    # test_features = apply_function_to_dict(lambda x: addresses_to_features(x, model_name=model_name), test_addresses)
    features = apply_function_to_dict(lambda x: addresses_to_features_parallel(x, model_name=model_name, n_processes=4), addresses)
    sys.stderr.write('Map addresses to features: %.2f seconds\n' % (time.time() - t))
    
    # Remove features that are None

    for name in features.keys():
        valid = [(ftr, addr) for ftr, addr in zip(features[name], addresses[name])
                    if ftr is not None]
        res = zip(*valid)
        features[name] = np.array(res[0])
        addresses[name] = res[1]
        
    return features, addresses
    

def grid_parameters_to_sample_locations(grid_spec=None, patch_size=None, stride=None, w=None, h=None):
    # patch_size, stride, w, h = grid_parameters.tolist()

    if grid_spec is not None:
        patch_size, stride, w, h = grid_spec

    half_size = patch_size/2
    ys, xs = np.meshgrid(np.arange(half_size, h-half_size, stride), np.arange(half_size, w-half_size, stride),
                     indexing='xy')
    sample_locations = np.c_[xs.flat, ys.flat]
    return sample_locations


def addresses_to_structure_distances(addresses, structure_centers_all_stacks_all_secs_all_names):
    """
    Return a list of dict {structure_name: distance}.
    """
    augmented_addresses = [addr + (i,) for i, addr in enumerate(addresses)]

    all_stacks = set([ad[0] for ad in addresses])
    grid_locations = {st: grid_parameters_to_sample_locations(get_default_gridspec(st)) for st in all_stacks}

    sorted_addresses = sorted(augmented_addresses, key=lambda (st,sec,gp_ind,list_ind): (st,sec))

    distances_to_structures = defaultdict(list)

    for (st, sec), address_group_ in groupby(sorted_addresses, lambda (st,sec,gp_ind,list_ind): (st,sec)):
        address_group = list(address_group_) # otherwise it is an iteraror, which can only be used once
        locations_this_group = np.array([grid_locations[st][gp_ind] for st,sec,gp_ind,list_ind in address_group])
        list_indices_this_group = [list_ind for st,sec,gp_ind,list_ind in address_group]

        for structure_name in structure_centers_all_stacks_all_secs_all_names[st][sec].keys():
            distances = np.sqrt(np.sum((locations_this_group - structure_centers_all_stacks_all_secs_all_names[st][sec][structure_name])**2, axis=1))
            distances_to_structures[structure_name] += zip(list_indices_this_group, distances)

    distances_to_structures.default_factory = None

    d = {structure_name: dict(lst) for structure_name, lst in distances_to_structures.iteritems()}

    import pandas
    df = pandas.DataFrame(d)
    return df.T.to_dict().values()

def addresses_to_locations(addresses):
    """
    Take a list of (stack, section, gridpoint_index),
    return x,y coordinates (lossless).
    """

    augmented_addresses = [addr + (i,) for i, addr in enumerate(addresses)]
    sorted_addresses = sorted(augmented_addresses, key=lambda (st,sec,gp_ind,list_ind): st)

    locations = []
    for st, address_group in groupby(sorted_addresses, lambda (st,sec,gp_ind,list_ind): st):
        grid_locations = grid_parameters_to_sample_locations(get_default_gridspec(st))
        locations_this_group = [(list_inf, grid_locations[gp_ind]) for st,sec,gp_ind,list_ind in address_group]
        locations.append(locations_this_group)

    locations = list(chain(*locations))
    return [v for i, v in sorted(locations)]

def addresses_to_features_parallel(addresses, model_name, n_processes=16):
    """
    If certain input address is outside the mask, the corresponding feature returned is None.
    """

    groups = [(st_se, list(group)) for st_se, group in groupby(sorted(enumerate(addresses), key=lambda (i,(st,se,idx)): (st, se)), key=lambda (i,(st,se,idx)): (st, se))]
    
    def f(x):
        st_se, group = x

        list_indices, addrs = zip(*group)
        sampled_grid_indices = [idx for st, se, idx in addrs]
        stack, sec = st_se
        fn = metadata_cache['sections_to_filenames'][stack][sec]

        if is_invalid(fn):
            sys.stderr.write('Image file is %s.\n' % (fn))
            features_ret = [None for _ in sampled_grid_indices]
        else:
            # Load mapping grid index -> location
            all_grid_indices, _ = DataManager.load_dnn_feature_locations(stack=stack, model_name=model_name, fn=fn)
            all_grid_indices = all_grid_indices.tolist()

            sampled_list_indices = []
            for gi, lst_idx in zip(sampled_grid_indices, list_indices):
                if gi in all_grid_indices:
                    sampled_list_indices.append(all_grid_indices.index(gi))
                else:
                    sys.stderr.write('Patch in annotation but not in mask: %s %d %s @%d\n' % (stack, sec, fn, gi))
                    sampled_list_indices.append(None)

            features = DataManager.load_dnn_features(stack=stack, model_name=model_name, fn=fn)
            features_ret = [features[i].copy() if i is not None else None for i in sampled_list_indices]
            del features

        return features_ret, list_indices
        
    from multiprocess import Pool
    pool = Pool(n_processes)
    res = pool.map(f, groups)    
    pool.close()
    pool.join()
    
    feature_list, list_indices = zip(*res)
    feature_list = list(chain(*feature_list))
    list_indices_all_stack_section = list(chain(*list_indices))

    return [feature_list[i] for i in np.argsort(list_indices_all_stack_section)]


def addresses_to_features(addresses, model_name='Sat16ClassFinetuned'):
    """
    If certain input address is outside the mask, the corresponding feature returned is None.
    """

    feature_list = []
    list_indices_all_stack_section = []

    for st_se, group in groupby(sorted(enumerate(addresses), key=lambda (i,(st,se,idx)): (st, se)), key=lambda (i,(st,se,idx)): (st, se)):

        print st_se

        list_indices, addrs = zip(*group)
        sampled_grid_indices = [idx for st, se, idx in addrs]

        stack, sec = st_se

        fn = metadata_cache['sections_to_filenames'][stack][sec]

        if is_invalid(fn):
            sys.stderr.write('Image file is %s.\n' % (fn))
            feature_list += [None for _ in sampled_grid_indices]
            list_indices_all_stack_section += list_indices
        else:
            # Load mapping grid index -> location
            all_grid_indices, _ = DataManager.load_dnn_feature_locations(stack=stack, model_name=model_name, fn=fn)
            all_grid_indices = all_grid_indices.tolist()

            sampled_list_indices = []
            for gi, lst_idx in zip(sampled_grid_indices, list_indices):
                if gi in all_grid_indices:
                    sampled_list_indices.append(all_grid_indices.index(gi))
                else:
                    sys.stderr.write('Patch in annotation but not in mask: %s %d %s @%d\n' % (stack, sec, fn, gi))
                    sampled_list_indices.append(None)

            features = DataManager.load_dnn_features(stack=stack, model_name=model_name, fn=fn)
            feature_list += [features[i].copy() if i is not None else None for i in sampled_list_indices]
            del features

            list_indices_all_stack_section += list_indices

    return [feature_list[i] for i in np.argsort(list_indices_all_stack_section)]


def visualize_filters(model, name, input_channel=0, title=''):

    filters = model.arg_params[name].asnumpy()

    n = len(filters)

    ncol = 16
    nrow = n/ncol

    fig, axes = plt.subplots(nrow, ncol, figsize=(ncol*.5, nrow*.5), sharex=True, sharey=True)

    fig.suptitle(title)

    axes = axes.flatten()
    for i in range(n):
        axes[i].matshow(filters[i][input_channel], cmap=plt.cm.gray)
        axes[i].tick_params(
            axis='both',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            left='off',
            right='off',
            labelbottom='off',
            labeltop='off',
            labelright='off',
            labelleft='off') # labels along the bottom edge are off
        axes[i].axis('equal')
    plt.show()
