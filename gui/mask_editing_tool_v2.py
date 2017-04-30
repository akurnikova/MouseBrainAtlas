#! /usr/bin/env python

import sys, os
import argparse

import matplotlib.pyplot as plt
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from multiprocess import Pool

sys.path.append(os.environ['REPO_DIR'] + '/utilities')
from utilities2015 import *
from metadata import *
from data_manager import *
from registration_utilities import find_contour_points
from gui_utilities import *
from qt_utilities import *
from preprocess_utilities import *
sys.path.append(os.path.join(os.environ['REPO_DIR'], 'web_services'))
from web_service import web_services_request

from ui.ui_MaskEditingGui5 import Ui_MaskEditingGui
from widgets.DrawableZoomableBrowsableGraphicsScene import DrawableZoomableBrowsableGraphicsScene
from widgets.DrawableZoomableBrowsableGraphicsScene_ForMasking import DrawableZoomableBrowsableGraphicsScene_ForMasking
from widgets.DrawableZoomableBrowsableGraphicsScene_ForSnake import DrawableZoomableBrowsableGraphicsScene_ForSnake
from widgets.ZoomableBrowsableGraphicsScene import ZoomableBrowsableGraphicsScene
from DataFeeder import ImageDataFeeder

#####################################################################

STR_USING_AUTO = 'Using AUTO (Click to switch to USER)'
STR_USING_USER = 'Using USER (Click to switch to AUTO)'

class MaskEditingGUI(QMainWindow):
    def __init__(self, parent=None, stack=None):
        QMainWindow.__init__(self, parent)

        self.stack = stack

        self.ui = Ui_MaskEditingGui()
        self.dialog = QDialog(self)
        self.ui.setupUi(self.dialog)

        self.ui.button_snake.clicked.connect(self.do_snake_current_section)
        self.ui.button_update_merged_mask.clicked.connect(self.update_merged_mask_button_clicked)
        self.ui.button_toggle_accept_auto.clicked.connect(self.toggle_accept_auto)
        self.ui.button_toggle_accept_auto.setText(STR_USING_AUTO)
        self.ui.button_autoSnake.clicked.connect(self.snake_all)
        self.ui.button_loadAnchorContours.clicked.connect(self.load_anchor_contours)
        self.ui.button_saveAnchorContours.clicked.connect(self.save_anchor_contours)
        self.ui.button_loadAllInitContours.clicked.connect(self.load_all_init_snake_contours)
        self.ui.button_saveAllInitContours.clicked.connect(self.save_all_init_snake_contours)
        self.ui.button_saveFinalMasks.clicked.connect(self.save_final_masks)

        self.ui.slider_snakeShrink.setSingleStep(1)
        self.ui.slider_snakeShrink.setMinimum(0)
        self.ui.slider_snakeShrink.setMaximum(20)
        self.ui.slider_snakeShrink.setValue(MORPHSNAKE_LAMBDA1)
        self.ui.slider_snakeShrink.valueChanged.connect(self.snake_shrinkParam_changed)

        self.sections_to_filenames = DataManager.load_sorted_filenames(stack)[1]
        # self.sections_to_filenames = {sec: fn for sec, fn in self.sections_to_filenames.iteritems() if sec > 263 and sec < 267}
        self.valid_sections_to_filenames = {sec: fn for sec, fn in self.sections_to_filenames.iteritems() if not is_invalid(fn)}
        self.valid_filenames_to_sections = {fn: sec for sec, fn in self.valid_sections_to_filenames.iteritems()}
        q = sorted(self.valid_sections_to_filenames.items())
        self.valid_sections = [sec for sec, fn in q]
        self.valid_filenames = [fn for sec, fn in q]

        ########################################################

        self.original_images = {}
        self.selected_channels = {}
        self.thresholded_images = {}
        self.contrast_stretched_images = {}
        self.selected_snake_lambda1 = {}
        self.init_submasks = {}
        self.init_submasks_vizs = {}
        # self.auto_submasks = {}
        self.user_submasks = {}
        self.accepted_final_masks = {}
        self.accept_which = {sec: 0 for sec in self.valid_sections}
        self.merged_masks = {}
        self.merged_mask_vizs = {}

        # self.auto_submask_decisions = {}
        self.user_submask_decisions = {}

        # anchor_fn = metadata_cache['anchor_fn'][stack]

        # Load decisions from final decision file.

        from pandas import read_csv

        auto_submask_rootdir = DataManager.get_auto_submask_rootdir_filepath(stack)
        for fn in self.valid_filenames:
            self.user_submasks[fn] = {}
            decision_fp = DataManager.get_auto_submask_filepath(stack=stack, what='decisions', fn=fn)
            decisions = read_csv(decision_fp).to_dict()
            self.user_submask_decisions[fn] = decisions
            for submask_ind, decision in decisions.iteritems():
                submask_fp = DataManager.get_auto_submask_filepath(stack=stack, what='submask', fn=fn, submask_ind=submask_ind)
                self.user_submasks[fn][submask_ind] = imread(submask_fp).astype(np.bool)

        # # Load decisions from final decision file.
        # try:
        #     accept_which, submask_decisions = load_final_decisions(self.stack)
        #     for fn, which in accept_which.iteritems():
        #         if fn not in self.valid_filenames:
        #             continue
        #         sec = self.valid_filenames_to_sections[fn]
        #         if which == 0:
        #             self.auto_submask_decisions[sec] = submask_decisions[fn]
        #         elif which == 1:
        #             self.user_submask_decisions[sec] = submask_decisions[fn]
        #         else:
        #             raise Exception('The type indicator %s must be either 0 or 1, but is %d' % (fn, which))
        #         self.accept_which[sec] = which
        #     sys.stderr.write('Decisions loaded from the final decision file.\n')
        # except Exception as e:
        #     sys.stderr.write('Error loading final decision file.\n')

        ######################################
        ## Generate submask review results. ##
        ######################################

        self.auto_submasks_gscene = DrawableZoomableBrowsableGraphicsScene_ForSnake(id='autoFinalMasks', gview=self.ui.init_snake_contour_gview)
        # self.auto_masks_feeder = ImageDataFeeder(name='autoFinalMasks', stack=self.stack, \
        #                             sections=self.valid_sections, use_data_manager=False, downscale=32,
        #                             labeled_filenames={sec: os.path.join(RAW_DATA_DIR, self.stack, fn + ".png")
        #                                 for sec, fn in self.valid_sections_to_filenames.iteritems()})
        self.auto_masks_feeder = ImageDataFeeder(name='autoFinalMasks', stack=self.stack, \
                                    sections=self.valid_sections, use_data_manager=True,
                                    downscale=32,
                                    version='aligned')
                                    # labeled_filenames={sec: os.path.join(RAW_DATA_DIR, self.stack, fn + ".png")
                                        # for sec, fn in self.valid_sections_to_filenames.iteritems()})
        self.auto_submasks_gscene.set_data_feeder(self.auto_masks_feeder)
        self.auto_submasks_gscene.active_image_updated.connect(self.auto_submasks_gscene_section_changed)
        self.auto_submasks_gscene.submask_decision_updated.connect(self.auto_submask_decision_updated)

        # submasks_rootdir = os.path.join(THUMBNAIL_DATA_DIR, self.stack, self.stack + '_submasks')

        def filter_by_keys(d, allowed_key_list):
            return {fn: v for fn, v in d.iteritems() if fn in allowed_key_list}

        def convert_keys_fn_to_sec(d):
            return {self.valid_filenames_to_sections[fn]: v for fn, v in d.iteritems()}

        # try:
        #     auto_submasks = load_submasks(submasks_rootdir=submasks_rootdir)
        #     self.auto_submasks = convert_keys_fn_to_sec(filter_by_keys(auto_submasks, self.valid_filenames))
        # except:
        #     self.auto_submasks = {}

        # # If user decisions exist
        # auto_submask_decisions = generate_submask_review_results(submasks_rootdir=submasks_rootdir, filenames=self.valid_filenames, which='user')
        # auto_submask_decisions = convert_keys_fn_to_sec(filter_by_keys(auto_submask_decisions, self.valid_filenames))
        # for sec, decisions in auto_submask_decisions.iteritems():
        #     if sec not in self.auto_submask_decisions or len(self.auto_submask_decisions[sec]) == 0:
        #         self.auto_submask_decisions[sec] = decisions

        # If no user decisions, load auto decisions
        # auto_submask_decisions = generate_submask_review_results(submasks_rootdir=submasks_rootdir, filenames=self.valid_filenames, which='auto')
        # auto_submask_decisions = convert_keys_fn_to_sec(filter_by_keys(auto_submask_decisions, self.valid_filenames))
        # for sec, decisions in auto_submask_decisions.iteritems():
        #     if sec not in self.auto_submask_decisions or len(self.auto_submask_decisions[sec]) == 0:
        #         self.auto_submask_decisions[sec] = decisions
        # self.auto_submasks_gscene.set_submasks_and_decisions(self.auto_submasks, self.auto_submask_decisions)

        self.anchor_fn = DataManager.load_anchor_filename(stack=self.stack)
        filenames_to_sections, _ = DataManager.load_sorted_filenames(stack=self.stack)
        self.auto_submasks_gscene.set_active_section(filenames_to_sections[self.anchor_fn], emit_changed_signal=False)

        ##########################
        ## User Submasks Gscene ##
        ##########################

        self.user_submasks_gscene = DrawableZoomableBrowsableGraphicsScene_ForMasking(id='userFinalMask', gview=self.ui.gview_final_masks_user)
        # self.user_submasks_feeder = ImageDataFeeder(name='autoFinalMasks', stack=self.stack, \
        #                             sections=self.valid_sections, use_data_manager=False, downscale=32,
        #                             labeled_filenames={sec: os.path.join(RAW_DATA_DIR, self.stack, fn + ".png")
        #                                 for sec, fn in self.valid_sections_to_filenames.iteritems()})
        self.user_submasks_feeder = ImageDataFeeder(name='autoFinalMasks', stack=self.stack, \
                                    sections=self.valid_sections, use_data_manager=True,
                                    downscale=32,
                                    version='aligned')
        self.user_submasks_gscene.set_data_feeder(self.user_submasks_feeder)
        self.user_submasks_gscene.submask_decision_updated.connect(self.user_submask_decision_updated)

        # # Load modified submasks and submask decisions.
        # user_modified_submasks_rootdir = create_if_not_exists(os.path.join(THUMBNAIL_DATA_DIR, self.stack, self.stack + '_submasks_user_modified'))
        #
        # user_submasks = load_submasks(submasks_rootdir=user_modified_submasks_rootdir)
        # self.user_submasks = convert_keys_fn_to_sec(filter_by_keys(user_submasks, self.valid_filenames))
        #
        # user_submask_decisions = generate_submask_review_results(submasks_rootdir=user_modified_submasks_rootdir, filenames=self.valid_filenames, which='user')
        # user_submask_decisions = convert_keys_fn_to_sec(filter_by_keys(user_submask_decisions, self.valid_filenames))
        # for sec, decisions in user_submask_decisions.iteritems():
        #     if sec not in self.user_submask_decisions or len(self.user_submask_decisions[sec]) == 0:
        #         self.user_submask_decisions[sec] = decisions

        # selected_snake_lambda1, selected_dissim_thresholds, selected_channels = load_masking_parameters(submasks_rootdir=user_modified_submasks_rootdir)
        # selected_snake_lambda1 = convert_keys_fn_to_sec(selected_snake_lambda1)
        # selected_dissim_thresholds = convert_keys_fn_to_sec(selected_dissim_thresholds)
        # selected_channels = convert_keys_fn_to_sec(selected_channels)
        #
        # for sec, v in selected_snake_lambda1.iteritems():
        #     self.selected_snake_lambda1[sec] = v
        # for sec, th in selected_dissim_thresholds.iteritems():
        #     self.selected_dissim_thresholds[sec] = th
        # for sec, ch in selected_channels.iteritems():
        #     self.selected_channels[sec] = ch

        self.user_submasks_gscene.set_submasks_and_decisions(self.user_submasks, self.user_submask_decisions)

        #########################################################

        self.ui.comboBox_channel.activated.connect(self.channel_changed)
        self.ui.comboBox_channel.addItems(['Red', 'Green', 'Blue'])

        #########################################################

        self.gscene_thresholded = ZoomableBrowsableGraphicsScene(id='thresholded', gview=self.ui.gview_thresholded)
        self.thresholded_image_feeder = ImageDataFeeder(name='thresholded', stack=self.stack, \
                                                        sections=self.valid_sections, use_data_manager=False,
                                                        downscale=32)
        self.gscene_thresholded.set_data_feeder(self.thresholded_image_feeder)

        #########################################################

        # self.gscene_slic = ZoomableBrowsableGraphicsScene(id='slic', gview=self.ui.gview_slic)
        # self.slic_image_feeder = ImageDataFeeder(name='slic', stack=self.stack, \
        #                                         sections=self.valid_sections, use_data_manager=False,
        #                                         downscale=32)
        # self.gscene_slic.set_data_feeder(self.slic_image_feeder)

        #########################################################

        # self.gscene_dissimmap = ZoomableBrowsableGraphicsScene(id='dissimmap', gview=self.ui.gview_dissimmap)
        # self.dissim_image_feeder = ImageDataFeeder(name='dissimmap', stack=self.stack, \
        #                                         sections=self.valid_sections, use_data_manager=False,
        #                                         downscale=32)
        # self.gscene_dissimmap.set_data_feeder(self.dissim_image_feeder)

        #########################################################

        # self.init_user_submasks_gscene = ZoomableBrowsableGraphicsScene(id='submasks', gview=self.ui.gview_submasks)
        # self.submask_image_feeder = ImageDataFeeder(name='submasks', stack=self.stack, \
        #                                         sections=self.valid_sections, use_data_manager=False,
        #                                         downscale=32)
        # self.init_user_submasks_gscene.set_data_feeder(self.submask_image_feeder)

        # self.init_user_submasks_gscene = DrawableZoomableBrowsableGraphicsScene_ForMasking(id='submasks', gview=self.ui.gview_submasks)
        # self.init_user_submasks_image_feeder = ImageDataFeeder(name='submasks', stack=self.stack, \
        #                                         sections=self.valid_sections, use_data_manager=False,
        #                                         downscale=32)
        # self.init_user_submasks_gscene.set_data_feeder(self.init_user_submasks_image_feeder)

        #########################################################

        self.gscene_merged_mask = ZoomableBrowsableGraphicsScene(id='mergedMask', gview=self.ui.gview_merged_mask)
        self.merged_masks_feeder = ImageDataFeeder(name='mergedMask', stack=self.stack, \
                                                sections=self.valid_sections, use_data_manager=False,
                                                downscale=32)
        self.gscene_merged_mask.set_data_feeder(self.merged_masks_feeder)

        #########################################################

        for sec in self.valid_sections:
            self.update_merged_mask(sec=sec)

        #########################################################

        self.dialog.showMaximized()


    # def upload_masks(self):
    #     transfer_data_synced(fp_relative=os.path.join(self.stack, self.stack + '_masks'),
    #                         from_hostname='localhost', to_hostname='dm')

    def load_anchor_contours(self):
        contours_on_anchor_sections = load_pickle(os.path.join(THUMBNAIL_DATA_DIR, self.stack, self.stack + '_alignedTo_' + self.anchor_fn + '_anchor_init_snake_contours.pkl'))
        for sec, vertices in contours_on_anchor_sections.iteritems():
            self.auto_submasks_gscene.set_init_snake_contour(section=sec, vertices=vertices)
            self.auto_submasks_gscene.set_section_as_anchor(section=sec)

    def save_anchor_contours(self):
        contours_on_anchor_sections = \
            {sec: vertices_from_polygon(self.auto_submasks_gscene.init_snake_contour_polygons[sec])
            for sec in self.auto_submasks_gscene.anchor_sections}
        save_pickle(contours_on_anchor_sections, os.path.join(THUMBNAIL_DATA_DIR, self.stack, self.stack + '_alignedTo_' + self.anchor_fn + '_anchor_init_snake_contours.pkl'))

    def load_all_init_snake_contours(self):
        init_snake_contours_on_all_sections = load_pickle(os.path.join(THUMBNAIL_DATA_DIR, self.stack, self.stack + '_alignedTo_' + self.anchor_fn + '_init_snake_contours.pkl'))
        for fn, vertices in init_snake_contours_on_all_sections.iteritems():
            self.auto_submasks_gscene.set_init_snake_contour(section=self.valid_filenames_to_sections[fn], vertices=vertices)

    def save_all_init_snake_contours(self):
        """Save initial snake contours for all sections."""
        init_snake_contours_on_all_sections = {}
        for sec, fn in self.valid_sections_to_filenames.iteritems():
            if sec in self.auto_submasks_gscene.init_snake_contour_polygons:
                init_snake_contours_on_all_sections[fn] = vertices_from_polygon(self.auto_submasks_gscene.init_snake_contour_polygons[sec])
            else:
                sys.stderr.write("Image %s (section %d) does not have any initial snake contour.\n" % (fn, sec))
        save_pickle(init_snake_contours_on_all_sections, os.path.join(THUMBNAIL_DATA_DIR, self.stack, self.stack + '_alignedTo_' + self.anchor_fn + '_init_snake_contours.pkl'))

    def save_final_masks(self):
        submasks_dir = create_if_not_exists(os.path.join(THUMBNAIL_DATA_DIR, self.stack, self.stack + '_alignedTo_' + self.anchor_fn + '_submasks_user_modified'))
        for sec in self.valid_sections:
            self.save_submasks_and_decisions(submasks_dir=submasks_dir, sec=sec)
        self.export_final_masks()

    def export_final_masks(self):
        final_masks_dir = create_if_not_exists(os.path.join(THUMBNAIL_DATA_DIR, self.stack, self.stack + '_alignedTo_' + self.anchor_fn + '_masks'))
        for sec, mask_viz in self.merged_mask_vizs.iteritems():
            fn = self.valid_sections_to_filenames[sec]
            imsave(os.path.join(final_masks_dir, fn + '_alignedTo_' + self.anchor_fn + '_mask.png'), mask_viz)

    def save_submasks_and_decisions(self, submasks_dir, sec):
        submasks = self.user_submasks
        submask_decisions = self.user_submask_decisions

        if sec not in submasks or sec not in submask_decisions:
            return

        fn = self.valid_sections_to_filenames[sec]

        submask_fn_dir = os.path.join(submasks_dir, fn)
        execute_command('rm -rf %(d)s; mkdir -p %(d)s' % {'d': submask_fn_dir})

        # Save submasks
        for submask_ind, m in submasks[sec].iteritems():
            submask_fp = os.path.join(submask_fn_dir, fn + '_alignedTo_' + self.anchor_fn + '_submask_%d.png' % submask_ind)
            imsave(submask_fp, np.uint8(m)*255)

        # Save submask contour vertices
        submask_contour_vertices_fp = os.path.join(submask_fn_dir, fn + '_alignedTo_' + self.anchor_fn + '_submask_contour_vertices.pkl')
        submask_contour_vertices_dict = {}
        for submask_ind, m in submasks[sec].iteritems():
            cnts = find_contour_points(m)[1]
            assert len(cnts) == 1, "Must have exactly one contour per submask."
            submask_contour_vertices_dict[submask_ind] = cnts[0]
        save_pickle(submask_contour_vertices_dict, submask_contour_vertices_fp)

        # Save submask decisions
        decisions_fp = os.path.join(submask_fn_dir, fn +'_alignedTo_' + self.anchor_fn +  '_submasksUserReview.txt')
        from pandas import Series
        Series(submask_decisions[sec]).to_csv(decisions_fp)

        if sec in self.selected_channels or sec in self.selected_snake_lambda1:
            parameters_fp = os.path.join(submask_fn_dir, fn + '_alignedTo_' + self.anchor_fn + '_maskingParameters.txt')
            with open(parameters_fp, 'w') as f:
                if sec in self.selected_snake_lambda1:
                    f.write('snake_lambda1 %d\n' % self.selected_snake_lambda1[sec])
                if sec in self.selected_channels:
                    f.write('channel %d\n' % self.selected_channels[sec])
        else:
            sys.stderr.write('Parameters for %s(%d) is not saved (no modification made ?)\n' % (fn, sec))

    def save_final_decisions(self):

        accept_which_fp = os.path.join(THUMBNAIL_DATA_DIR, self.stack, self.stack + '_submasks_finalDecisions.txt')
        with open(accept_which_fp, 'w') as f:
            for sec, accept_which in sorted(self.accept_which.items()):
                fn = self.valid_sections_to_filenames[sec]
                if accept_which == 0:
                    if sec not in self.auto_submask_decisions:
                        sys.stderr.write('No auto submask decisions for %s (%d)\n' % (fn, sec))
                        continue
                    decisions = self.auto_submask_decisions[sec]
                    if len(decisions) == 0:
                        sys.stderr.write('Auto submask decisions are empty for %s (%d)\n' % (fn, sec))
                elif accept_which == 1:
                    if sec not in self.user_submask_decisions:
                        sys.stderr.write('No user submask decisions for %s (%d)\n' % (fn, sec))
                        continue
                    decisions = self.user_submask_decisions[sec]
                    if len(decisions) == 0:
                        sys.stderr.write('User submask decisions are empty for %s (%d)\n' % (fn, sec))
                else:
                    raise

                f.write('%d %s %d %s\n' % (sec, fn, accept_which, ' '.join(map(lambda x: str(int(x)), decisions))))


    # def save_submasks_and_decisions(self, submasks_dir, which, fn=None, sec=None):
    #     """
    #     If which is auto, save submasks, submask decisions and parameters to both modified and final folder.
    #     If which is user, save submasks and submask decisions to both modified and final folder.
    #     """
    #
    #     if which == 'auto':
    #         submasks = self.auto_submasks
    #         submask_decisions = self.auto_submask_decisions
    #     elif which == 'user':
    #         submasks = self.user_submasks
    #         submask_decisions = self.user_submask_decisions
    #     else:
    #         raise
    #
    #     if sec is None:
    #         sec = self.auto_submasks_gscene.active_section
    #     if sec not in submasks or sec not in submask_decisions:
    #         return
    #     fn = self.valid_sections_to_filenames[sec]
    #
    #     submask_fn_dir = os.path.join(submasks_dir, fn)
    #     execute_command('rm -rf %(d)s; mkdir -p %(d)s' % {'d': submask_fn_dir})
    #
    #     # Save submasks
    #     for submask_ind, m in enumerate(submasks[sec]):
    #         submask_fp = os.path.join(submask_fn_dir, fn + '_submask_%d.png' % submask_ind)
    #         imsave(submask_fp, np.uint8(m)*255)
    #
    #     # Save submask decisions
    #     decisions_fp = os.path.join(submask_fn_dir, fn + '_submasksUserReview.txt')
    #     np.savetxt(decisions_fp, submask_decisions[sec], fmt='%d')
    #
    #     # Save masking parameters
    #     if which == 'user':
    #         if sec in self.selected_dissim_thresholds or \
    #             sec in self.selected_channels or \
    #             sec in self.selected_snake_lambda1:
    #
    #             parameters_fp = os.path.join(submask_fn_dir, fn + '_maskingParameters.txt')
    #             with open(parameters_fp, 'w') as f:
    #                 if sec in self.selected_snake_lambda1:
    #                     f.write('snake_lambda1 %d\n' % self.selected_snake_lambda1[sec])
    #                 if sec in self.selected_dissim_thresholds:
    #                     f.write('dissim_threshold %.2f\n' % self.selected_dissim_thresholds[sec])
    #                 if sec in self.selected_channels:
    #                     f.write('channel %d\n' % self.selected_channels[sec])
    #         else:
    #             sys.stderr.write('Parameters for %s(%d) is not saved (no modification made ?)\n' % (fn, sec))

    # def save(self, sec=None, fn=None):
    #
    #     if sec is None:
    #         sec = self.auto_submasks_gscene.active_section
    #     accept_which = ['auto', 'user'][self.accept_which[sec]]
    #
    #     if accept_which == 'user':
    #         submasks_dir = create_if_not_exists(os.path.join(THUMBNAIL_DATA_DIR, self.stack, self.stack + '_submasks_user_modified'))
    #         self.save_submasks_and_decisions(submasks_dir=submasks_dir, sec=sec, which=accept_which)

    def set_accept_auto_to_true(self, section):
        assert section in self.auto_submasks
        # Clear later stage images.
        self.accepted_final_masks[section] = self.auto_submasks[section]
        self.accept_which[section] = 0 # change to accept auto
        # self.ui.button_toggle_accept_auto.setText(STR_USING_AUTO)

    def set_accept_auto_to_false(self, section):
        assert section in self.user_submasks
        self.accepted_final_masks[section] = self.user_submasks[section]
        self.accept_which[section] = 1 # change to accept modified
        # self.ui.button_toggle_accept_auto.setText(STR_USING_USER)

    def toggle_accept_auto(self):

        sec = self.auto_submasks_gscene.active_section
        if self.accept_which[sec] == 0: # currently accepting auto
            self.set_accept_auto_to_false(section=sec)
        elif self.accept_which[sec] == 1: # currently accepting modified
            self.set_accept_auto_to_true(section=sec)
        self.update_merged_mask()

    @pyqtSlot(int)
    def user_submask_decision_updated(self, submask_ind):
        self.update_merged_mask()
        self.update_mask_gui_window_title()

    @pyqtSlot(int)
    def auto_submask_decision_updated(self, submask_ind):
        self.update_merged_mask()
        self.update_mask_gui_window_title()

    @pyqtSlot()
    def update_merged_mask_button_clicked(self):
        self.update_merged_mask()

    def update_merged_mask(self, sec=None):
        """
        Update merged mask. Change the image shown in "Merged Mask" panel.
        """

        if sec is None:
            sec = self.auto_submasks_gscene.active_section

        if sec not in self.user_submask_decisions:
            return

        fn = self.valid_sections_to_filenames[sec]
        accept_which = self.accept_which[sec]
        # try:
            # if accept_which == 0:
            #     if sec not in self.auto_submask_decisions or len(self.auto_submask_decisions[sec]) == 0:
            #         # sys.stderr.write('Error: section %d, %s, accept auto but auto decisions is empty.\n' % (sec, fn))
            #         raise Exception('Error: section %d, %s, accept auto but auto decisions is empty.' % (sec, fn))
            #     merged_mask = np.any([self.auto_submasks[sec][si] for si, dec in enumerate(self.auto_submask_decisions[sec]) if dec], axis=0)
            # elif accept_which == 1:
                # if sec not in self.user_submask_decisions or len(self.user_submask_decisions[sec]) == 0:
                #     raise Exception('Error: section %d, %s, accept user but user decisions is empty.' % (sec, fn))

        print self.user_submask_decisions[sec]
        merged_mask = np.any([self.user_submasks[sec][sm_i] for sm_i, dec in self.user_submask_decisions[sec].iteritems() if dec], axis=0)
        # else:
        #     raise Exception('accept_which is neither 0 or 1.')
        self.merged_masks[sec] = merged_mask
        self.merged_mask_vizs[sec] = img_as_ubyte(self.merged_masks[sec])
        self.merged_masks_feeder.set_image(sec=sec, numpy_image=self.merged_mask_vizs[sec])
        self.gscene_merged_mask.update_image(sec=sec)
        # except Exception as e:
        #     # sys.stderr.write('%s\n' % e)
        #     raise e

    # def update_slic(self):
    #     sec = self.auto_submasks_gscene.active_section
    #
    #     t = time.time()
    #     self.slic_labelmaps[sec] = slic(self.contrast_stretched_images[sec].astype(np.float),
    #                                 sigma=SLIC_SIGMA, compactness=SLIC_COMPACTNESS,
    #                                 n_segments=SLIC_N_SEGMENTS, multichannel=False, max_iter=SLIC_MAXITER)
    #     sys.stderr.write('SLIC: %.2f seconds.\n' % (time.time() - t)) # 10 seconds, iter=100, nseg=1000;
    #
    #     self.slic_boundary_images[sec] = img_as_ubyte(mark_boundaries(self.contrast_stretched_images[sec],
    #                                         label_img=self.slic_labelmaps[sec],
    #                                         background_label=-1, color=(1,0,0)))
    #
    #     self.slic_image_feeder.set_image(sec=sec, numpy_image=self.slic_boundary_images[sec])
    #     self.gscene_slic.update_image(sec=sec)
    #
    #     ####
    #
    #     # self.ncut_labelmaps[sec] = normalized_cut_superpixels(self.contrast_stretched_images[sec], self.slic_labelmaps[sec])
    #
    #     self.ncut_labelmaps[sec] = self.slic_labelmaps[sec]
    #     self.sp_dissim_maps[sec] = compute_sp_dissims_to_border(self.contrast_stretched_images[sec], self.ncut_labelmaps[sec])
    #     # self.sp_dissim_maps[sec] = compute_sp_dissims_to_border(self.thresholded_images[sec], self.ncut_labelmaps[sec])
    #     # self.border_dissim_images[sec] = generate_dissim_viz(self.sp_dissim_maps[sec], self.ncut_labelmaps[sec])
    #     # self.dissim_image_feeder.set_image(sec=sec, numpy_image=self.border_dissim_images[sec])
    #     # self.gscene_dissimmap.update_image(sec=sec)
    #
    #     self.selected_dissim_thresholds[sec] = determine_dissim_threshold(self.sp_dissim_maps[sec], self.ncut_labelmaps[sec])
    #     self.ui.slider_dissimThresh.setValue(int(self.selected_dissim_thresholds[sec]/0.01))
    #
    #     ######################################################
    #
    #     self.update_init_submasks_image()

    # def update_init_submasks_image(self):
    #     """
    #     Update the initial submasks for snake.
    #     """
    #
    #     sec = self.auto_submasks_gscene.active_section
    #
    #     self.selected_dissim_thresholds[sec] = self.ui.slider_dissimThresh.value() * 0.01
    #     # self.init_submasks[sec] = get_submasks(self.thresholded_images[sec])
    #     # self.init_submasks[sec] = get_submasks(ncut_labels=self.ncut_labelmaps[sec], sp_dissims=self.sp_dissim_maps[sec], dissim_thresh=self.selected_dissim_thresholds[sec])
    #     self.init_submasks[sec] = merge_overlapping_masks(get_submasks(ncut_labels=self.ncut_labelmaps[sec], sp_dissims=self.sp_dissim_maps[sec], dissim_thresh=self.selected_dissim_thresholds[sec]))
    #     self.init_submasks_vizs[sec] = generate_submasks_viz(self.contrast_stretched_images[sec], self.init_submasks[sec], color=(255,0,0))
    #     # self.init_submasks_vizs[sec] = generate_submasks_viz(self.thresholded_images[sec], self.init_submasks[sec], color=(255,0,0))
    #     self.init_user_submasks_image_feeder.set_image(sec=sec, numpy_image=self.init_submasks_vizs[sec])
    #     self.init_user_submasks_gscene.update_image(sec=sec)
    #
    #     if sec in self.selected_snake_lambda1:
    #         self.ui.slider_snakeShrink.setValue(self.selected_snake_lambda1[sec])
    #     else:
    #         self.ui.slider_snakeShrink.setValue(MORPHSNAKE_LAMBDA1)

    def snake_all(self):

        self.prepare_contrast_stretched_images()

        init_snake_contour_vertices = {sec: vertices_from_polygon(self.auto_submasks_gscene.init_snake_contour_polygons[sec])
                                        for sec in self.valid_sections}

        t = time.time()
        pool = Pool(NUM_CORES/2)
        # This will cause TypeError: can't pickle PyCapsule objects; Must not retain any "self" in the arguments.
        # submasks_all_sections = pool.map(lambda sec: snake(img=self.contrast_stretched_images[sec], init_contours=[init_snake_contour_vertices[sec]], lambda1=1.),
        #     self.valid_sections)
        images = {sec: self.contrast_stretched_images[sec] for sec in self.valid_sections}
        submasks_all_sections = pool.map(lambda sec: snake(img=images[sec], init_contours=[init_snake_contour_vertices[sec]], lambda1=1.),
            self.valid_sections)
        ### This works too. ###
        # init_snake_contour_vertices = [vertices_from_polygon(self.auto_submasks_gscene.init_snake_contour_polygons[sec])
        #                                 for sec in self.valid_sections]
        # images = [self.contrast_stretched_images[sec] for sec in self.valid_sections]
        # submasks_all_sections = pool.map(lambda (img, init_cnt): snake(img=img, init_contours=[init_cnt], lambda1=1.),
        #     zip(images, init_snake_contour_vertices))
        pool.close()
        pool.join()
        sys.stderr.write("Snake all: %.2f seconds.\n" % (time.time() - t))

        for sec, submasks in zip(self.valid_sections, submasks_all_sections):
            if len(submasks) == 0:
                sys.stderr.write('No submasks found for section %d.\n' % sec)
                return
            # Cast to dict - is this necessary ?
            self.user_submasks[sec] = {i: m for i, m in enumerate(submasks)}
            # self.user_submask_decisions[sec] = {i: d for i, d in enumerate(auto_judge_submasks(submasks))}
            self.user_submask_decisions[sec] = {i: True for i, d in enumerate(submasks)}

            self.user_submasks_gscene.update_image_from_submasks_and_decisions(sec=sec)
            self.set_accept_auto_to_false(section=sec)
            self.update_merged_mask(sec=sec)

    def prepare_contrast_stretched_images(self):
        for sec in self.valid_sections:
            if sec not in self.contrast_stretched_images:
                if sec not in self.original_images:
                    # img = imread(DataManager.get_image_filepath(stack=self.stack, section=sec, resol='thumbnail', version='original_png'))
                    img = imread(DataManager.get_image_filepath(stack=self.stack, section=sec, resol='thumbnail', version='aligned'))
                    self.original_images[sec] = brightfieldize_image(img)

                if sec in self.selected_channels:
                    channel = self.selected_channels[sec]
                else:
                    channel = 0
                self.contrast_stretched_images[sec] = contrast_stretch_image(self.original_images[sec][..., channel])
                self.update_thresholded_image(sec=sec)

    def do_snake(self, sec):
        # self.user_submasks_gscene.remove_submask_and_decisions_for_one_section(sec=sec)
        self.selected_snake_lambda1[sec] = self.ui.slider_snakeShrink.value()
        # self.user_submasks[sec] = snake(img=self.original_images[sec], submasks=self.init_submasks[sec])

        init_snake_contour_vertices = vertices_from_polygon(self.auto_submasks_gscene.init_snake_contour_polygons[sec])

        if sec not in self.contrast_stretched_images:
            if sec not in self.original_images:
                # img = imread(DataManager.get_image_filepath(stack=self.stack, section=sec, resol='thumbnail', version='original_png'))
                img = imread(DataManager.get_image_filepath(stack=self.stack, section=sec, resol='thumbnail', version='aligned'))
                border = np.median(np.concatenate([img[:10, :].flatten(), img[-10:, :].flatten(), img[:, :10].flatten(), img[:, -10:].flatten()]))
                if border < 123:
                    # dark background, fluorescent
                    img = img.max() - img # invert, make tissue dark on bright background
                self.original_images[sec] = img
            if sec in self.selected_channels:
                channel = self.selected_channels[sec]
            else:
                channel = 0
            self.contrast_stretched_images[sec] = contrast_stretch_image(self.original_images[sec][..., channel])

        submasks = snake(img=self.contrast_stretched_images[sec],
                        init_contours=[init_snake_contour_vertices],
                        lambda1=self.selected_snake_lambda1[sec])
        if len(submasks) == 0:
            sys.stderr.write('No submasks found for section %d.\n' % sec)
            return

        self.user_submasks[sec] = {i: m for i, m in enumerate(submasks)}
        self.user_submask_decisions[sec] = {i: d for i, d in enumerate(auto_judge_submasks(submasks))}

        self.user_submasks_gscene.update_image_from_submasks_and_decisions(sec=sec)
        # self.user_submasks_gscene.update_image(sec=sec)

        self.set_accept_auto_to_false(section=sec)
        self.update_merged_mask(sec=sec)

    def do_snake_current_section(self):
        self.do_snake(sec=self.auto_submasks_gscene.active_section)

    # def update_user_submasks_image(self, sec):
    #
    #     # sec = self.auto_submasks_gscene.active_section
    #     self.user_submasks_gscene.remove_submask_and_decisions_for_one_section(sec=sec)
    #     self.selected_snake_lambda1[sec] = self.ui.slider_snakeShrink.value()
    #     # self.user_submasks[sec] = snake(img=self.original_images[sec], submasks=self.init_submasks[sec])
    #     submasks = snake(img=self.contrast_stretched_images[sec], submasks=self.init_submasks[sec],
    #                                     lambda1=self.selected_snake_lambda1[sec])
    #     if len(submasks) == 0:
    #         return
    #     else:
    #         self.user_submasks[sec] = submasks
    #     # self.user_submasks[sec] = snake(img=self.thresholded_images[sec], submasks=self.init_submasks[sec])
    #     # self.final_submasks_vizs[sec] = generate_submasks_viz(self.original_images[sec], self.user_submasks[sec], color=(255,0,0))
    #     self.user_submasks_gscene.update_image(sec=sec)
    #     self.user_submask_decisions[sec] = auto_judge_submasks(self.user_submasks[sec])
    #     self.user_submasks_gscene.add_submask_and_decision_for_one_section(submasks=self.user_submasks[sec],
    #     submask_decisions=self.user_submask_decisions[sec], sec=sec)
    #
    #     self.set_accept_auto_to_false(section=sec)
    #     self.update_merged_mask()

    def change_channel(self, channel):
        """
        Compute contrast_stretch_image based on selected_channels.
        """

        print 'Changed to', channel
        sec = self.auto_submasks_gscene.active_section
        self.contrast_stretched_images[sec] = contrast_stretch_image(self.original_images[sec][..., channel])
        self.update_thresholded_image()

    def channel_changed(self, index):
        if index == self.selected_channels[self.auto_submasks_gscene.active_section]:
            return

        self.selected_channels[self.auto_submasks_gscene.active_section] = index

        channel_text = str(self.sender().currentText())
        if channel_text == 'Red':
            self.change_channel(0)
        elif channel_text == 'Green':
            self.change_channel(1)
        elif channel_text == 'Blue':
            self.change_channel(2)

    # def dissim_threshold_changed(self, value):
    #     self.ui.label_dissimThresh.setText(str(value * 0.01))

    def snake_shrinkParam_changed(self, value):
        self.ui.label_snakeShrink.setText(str(value))

    def update_thresholded_image(self, sec=None):
        """
        Update the image in the thresholded image gscene, based on contrast_stretched_images.
        """

        print "update_thresholded_image"
        if sec is None:
            sec = self.auto_submasks_gscene.active_section
        self.thresholded_image_feeder.set_image(sec=sec, qimage=numpy_to_qimage(self.contrast_stretched_images[sec]))
        self.gscene_thresholded.update_image(sec=sec)

    def auto_submasks_gscene_section_changed(self):
        """
        What happens when the image in "Automatic Masks" panel is changed.
        """

        self.update_mask_gui_window_title()

        sec = self.auto_submasks_gscene.active_section

        # if sec not in self.contrast_stretched_images:
        #     if sec not in self.original_images:
        #         # img = imread(DataManager.get_image_filepath(stack=self.stack, section=sec, resol='thumbnail', version='original_png'))
        #         img = imread(DataManager.get_image_filepath(stack=self.stack, section=sec, resol='thumbnail', version='aligned'))
        #         border = np.median(np.concatenate([img[:10, :].flatten(), img[-10:, :].flatten(), img[:, :10].flatten(), img[:, -10:].flatten()]))
        #         if border < 123:
        #             # dark background, fluorescent
        #             img = img.max() - img # invert, make tissue dark on bright background
        #         self.original_images[sec] = img

        # if self.accept_which[sec] == 1:
        #     self.ui.button_toggle_accept_auto.setText(STR_USING_USER)
        # elif self.accept_which[sec] == 0:
        #     self.ui.button_toggle_accept_auto.setText(STR_USING_AUTO)

        # Set parameters if those for the current section have been modified before.

        if sec not in self.selected_channels:
            self.selected_channels[sec] = 0

        self.ui.comboBox_channel.setCurrentIndex(self.selected_channels[sec])
        self.change_channel(self.selected_channels[sec])

        try:
            self.gscene_thresholded.set_active_section(sec)
        except: # The first time this will complain "Image not loaded" yet. But will not once update_thresholded_image() loads the image.
            pass

        # self.update_thresholded_image()

        # try:
        #     self.gscene_slic.set_active_section(sec)
        # except:
        #     pass
        #
        # try:
        #     self.init_user_submasks_gscene.set_active_section(sec)
        # except:
        #     pass
        #
        # if sec in self.selected_dissim_thresholds:
        #     self.ui.slider_dissimThresh.setValue(int(self.selected_dissim_thresholds[sec]/0.01))
        # else:
        #     self.ui.slider_dissimThresh.setValue(0)

        if sec in self.selected_snake_lambda1:
            self.ui.slider_snakeShrink.setValue(self.selected_snake_lambda1[sec])
        else:
            self.ui.slider_snakeShrink.setValue(MORPHSNAKE_LAMBDA1)

        try:
            self.user_submasks_gscene.set_active_section(sec)
        except:
            pass

        try:
            self.gscene_merged_mask.set_active_section(sec)
        except:
            pass

    def update_mask_gui_window_title(self):
        curr_sec = self.auto_submasks_gscene.active_section
        curr_fn = self.valid_sections_to_filenames[curr_sec]
        try:
            title = '%s (%d) - Active: %s - Alg:%s User:%s' % (curr_fn, curr_sec, ['Alg', 'User'][self.accept_which[curr_sec]], self.auto_submask_decisions[curr_sec], self.user_submask_decisions[curr_sec])
            self.dialog.setWindowTitle(title)
            print title
        except:
            pass

    # def generate_masks(self):
    #     web_services_request('generate_masks', stack=self.stack, filenames=self.valid_filenames, tb_fmt='png')
    #     transfer_data_synced(fp_relative=os.path.join(self.stack, self.stack + '_masks'))

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Mask Editing GUI')
    parser.add_argument("stack", type=str, help="stack name")
    args = parser.parse_args()
    stack = args.stack

    app = QApplication(sys.argv)

    m = MaskEditingGUI(stack=stack)
    # m.showMaximized()
    m.show()
    sys.exit(app.exec_())