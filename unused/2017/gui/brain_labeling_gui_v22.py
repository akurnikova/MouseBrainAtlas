#! /usr/bin/env python

import sip
sip.setapi('QVariant', 2) # http://stackoverflow.com/questions/21217399/pyqt4-qtcore-qvariant-object-instead-of-a-string

import sys
import os
import datetime
from random import random
import subprocess
import time
import json
from pprint import pprint
import cPickle as pickle
from itertools import groupby
from operator import itemgetter

import numpy as np

from matplotlib.backends import qt4_compat
use_pyside = qt4_compat.QT_API == qt4_compat.QT_API_PYSIDE
if use_pyside:
    #print 'Using PySide'
    from PySide.QtCore import *
    from PySide.QtGui import *
else:
    #print 'Using PyQt4'
    from PyQt4.QtCore import *
    from PyQt4.QtGui import *

from ui_BrainLabelingGui_v15 import Ui_BrainLabelingGui
# from ui_RectificationTool import Ui_RectificationGUI
# from rectification_tool import *

from shapely.geometry import Polygon as ShapelyPolygon
from shapely.geometry import Point as ShapelyPoint
from shapely.geometry import LineString as ShapelyLineString
from shapely.geometry import LinearRing as ShapelyLineRing

from skimage.color import label2rgb

from gui_utilities import *

sys.path.append(os.environ['REPO_DIR'] + '/utilities')
from utilities2015 import *
from data_manager import DataManager
from metadata import *

from collections import defaultdict, OrderedDict, deque

from operator import attrgetter

import requests

from joblib import Parallel, delayed

# from LabelingPainter import LabelingPainter
from custom_widgets import *
from SignalEmittingItems import *
from DataFeeder import ImageDataFeeder, VolumeResectionDataFeeder

from drawable_gscene import *

#######################################################################

class ReadImagesThread(QThread):
    def __init__(self, stack, sections):
        QThread.__init__(self)
        self.stack = stack
        self.sections = sections

    def __del__(self):
        self.wait()

    def run(self):
        for sec in self.sections:
            image = QImage(DataManager.get_image_filepath(stack=self.stack, section=sec, version='rgb-jpg', data_dir=data_dir))
            self.emit(SIGNAL('image_loaded(QImage, int)'), image, sec)

class BrainLabelingGUI(QMainWindow, Ui_BrainLabelingGui):
# class BrainLabelingGUI(QMainWindow, Ui_RectificationGUI):

    def __init__(self, parent=None, stack=None, first_sec=None, last_sec=None, downsample=None):
        """
        Initialization of BrainLabelingGUI.
        """

        # t0 = time.time()

        # self.app = QApplication(sys.argv)
        QMainWindow.__init__(self, parent)

        self.stack = stack
        self.sagittal_downsample = downsample

        self.setupUi(self)

        self.button_save.clicked.connect(self.save)
        self.button_load.clicked.connect(self.load)
        self.lineEdit_username.returnPressed.connect(self.username_changed)

        from collections import defaultdict
        self.structure_volumes = {}

        # self.volume_cache = {32: bp.unpack_ndarray_file(volume_dir + '/%(stack)s/%(stack)s_down%(downsample)dVolume.bp' % {'stack':self.stack, 'downsample':32}),
        #                     8: bp.unpack_ndarray_file(volume_dir + '/%(stack)s/%(stack)s_down%(downsample)dVolume.bp' % {'stack':self.stack, 'downsample':8})}

        self.volume_cache = {32: bp.unpack_ndarray_file(volume_dir + '/%(stack)s/%(stack)s_down%(downsample)dVolume.bp' % {'stack':self.stack, 'downsample':32})}

        # self.volume = self.volume_cache[self.downsample_factor]
        # self.y_dim, self.x_dim, self.z_dim = self.volume.shape

        self.coronal_gscene = DrawableGraphicsScene(id='coronal', gui=self, gview=self.coronal_gview)
        self.coronal_gview.setScene(self.coronal_gscene)

        self.horizontal_gscene = DrawableGraphicsScene(id='horizontal', gui=self, gview=self.horizontal_gview)
        self.horizontal_gview.setScene(self.horizontal_gscene)

        self.sagittal_gscene = DrawableGraphicsScene(id='sagittal', gui=self, gview=self.sagittal_gview)
        self.sagittal_gview.setScene(self.sagittal_gscene)

        self.gscenes = {'coronal': self.coronal_gscene, 'sagittal': self.sagittal_gscene, 'horizontal': self.horizontal_gscene}

        for gscene in self.gscenes.itervalues():
            gscene.drawings_updated.connect(self.drawings_updated)
            gscene.crossline_updated.connect(self.crossline_updated)
            gscene.active_image_updated.connect(self.active_image_updated)
            gscene.update_structure_volume_requested.connect(self.update_structure_volume_requested)
            gscene.set_structure_volumes(self.structure_volumes)
            # gscene.set_drawings(self.drawings)

        from functools import partial
        self.gscenes['sagittal'].set_conversion_func_section_to_z(partial(DataManager.convert_section_to_z, stack=self.stack))
        self.gscenes['sagittal'].set_conversion_func_z_to_section(partial(DataManager.convert_z_to_section, stack=self.stack))

        ##################
        # self.slider_downsample.valueChanged.connect(self.downsample_factor_changed)

        ###################
        self.contextMenu_set = True

        self.recent_labels = []

        self.structure_names = load_structure_names(os.environ['REPO_DIR']+'/gui/structure_names.txt')
        self.new_labelnames = load_structure_names(os.environ['REPO_DIR']+'/gui/newStructureNames.txt')
        self.structure_names = OrderedDict(sorted(self.new_labelnames.items()) + sorted(self.structure_names.items()))

        self.installEventFilter(self)

        # self.keyPressEvent = self.key_pressed
        # self.keyReleaseEvent = self.key_released

        # self.sections = range(127, 327)
        # self.sections = range(150, 304)
        # self.sections = range(150, 160)

        if first_sec is None and last_sec is None:
            first_sec, last_sec = section_range_lookup[self.stack]

        self.sections = range(first_sec, last_sec + 1)

        image_feeder = ImageDataFeeder('image feeder', stack=self.stack, sections=self.sections)
        image_feeder.set_orientation('sagittal')
        # image_feeder.set_downsample_factor(1)
        image_feeder.set_downsample_factor(self.sagittal_downsample)
        self.gscenes['sagittal'].set_data_feeder(image_feeder)

        volume_resection_feeder = VolumeResectionDataFeeder('volume resection feeder', self.stack)

        coronal_volume_resection_feeder = VolumeResectionDataFeeder('coronal resection feeder', self.stack)
        coronal_volume_resection_feeder.set_volume_cache(self.volume_cache)
        coronal_volume_resection_feeder.set_orientation('coronal')
        coronal_volume_resection_feeder.set_downsample_factor(32)
        self.gscenes['coronal'].set_data_feeder(coronal_volume_resection_feeder)

        horizontal_volume_resection_feeder = VolumeResectionDataFeeder('horizontal resection feeder', self.stack)
        horizontal_volume_resection_feeder.set_volume_cache(self.volume_cache)
        horizontal_volume_resection_feeder.set_orientation('horizontal')
        horizontal_volume_resection_feeder.set_downsample_factor(32)
        self.gscenes['horizontal'].set_data_feeder(horizontal_volume_resection_feeder)

        # self.gscenes['coronal'].set_downsample_factor(self.downsample_factor)
        # self.gscenes['sagittal'].set_downsample_factor(1)
        # self.gscenes['horizontal'].set_downsample_factor(self.downsample_factor)

        if self.gscenes['sagittal'].data_feeder.downsample == 1:
            self.read_images_thread = ReadImagesThread(self.stack, self.sections)
            self.connect(self.read_images_thread, SIGNAL("image_loaded(QImage, int)"), self.image_loaded)
            self.read_images_thread.start()
            self.button_stop.clicked.connect(self.read_images_thread.terminate)
        else:
            self.gscenes['sagittal'].data_feeder.load_images()
            self.gscenes['sagittal'].set_vertex_radius(3)
            self.gscenes['sagittal'].set_line_width(3)

        self.gscenes['coronal'].set_active_i(50)
        self.gscenes['sagittal'].set_active_section(self.sections[0])
        self.gscenes['horizontal'].set_active_i(150)

        # print time.time() - t0

    @pyqtSlot()
    def image_loaded(self, qimage, sec):
        self.gscenes['sagittal'].data_feeder.set_image(qimage, sec)
        print 'Image', sec, 'received.'
        if self.gscenes['sagittal'].active_section == sec:
            self.gscenes['sagittal'].load_histology()
        self.statusBar().showMessage('Image %d loaded.\n' % sec)

    @pyqtSlot()
    def username_changed(self):
        self.username = str(self.sender().text())
        print 'username changed to', self.username

    @pyqtSlot()
    def save(self):

        if not hasattr(self, 'username') or self.username is None:
            username, okay = QInputDialog.getText(self, "Username", "Please enter your username:", QLineEdit.Normal, 'anon')
            if not okay: return
            self.username = str(username)
            self.lineEdit_username.setText(self.username)

        # labelings_dir = create_if_not_exists('/home/yuncong/CSHL_labelings_new/%(stack)s/' % dict(stack=self.stack))
        labelings_dir = create_if_not_exists(os.path.join(annotation_midbrainIncluded_v2_rootdir, stack))

        timestamp = datetime.datetime.now().strftime("%m%d%Y%H%M%S")

        for gscene_id, gscene in self.gscenes.iteritems():
            # gscene.save_drawings(fn_template=os.path.join(labelings_dir, '%(stack)s_%(orientation)s_%(downsample)d_%(username)s_%(timstamp)s.pkl' % dict(username=self.username)))
            # gscene.save_drawings(fn_template=os.path.join(labelings_dir, '%(stack)s_%(orientation)s_downsample%(downsample)d_'+self.username+'_'+timestamp+'.pkl'))
            gscene.save_drawings(fn_template=os.path.join(labelings_dir, '%(stack)s_%(orientation)s_downsample%(downsample)d_%(username)s_%(timestamp)s.pkl'), timestamp=timestamp, username=self.username)

        self.statusBar().showMessage('Labelings saved to %s.' % labelings_dir)

        # pickle.dump(self.structure_volumes, open(os.path.join(labelings_dir, '%(stack)s_structure_volumes.pkl' % dict(stack=stack))))

        # # if sec is not None:
        # #
        # #     accepted_proposal_props = []
        # #     for polygon, props in self.accepted_proposals_allSections[sec].iteritems():
        # #
        # #         props_saved = props.copy()
        # #
        # #         # props_saved['vertices'] = [(v.scenePos().x(), v.scenePos().y()) for v in props['vertexCircles']]
        # #
        # #         path = polygon.path()
        # #
        # #         if path.elementCount() > 1 and polygon_is_closed(path=path):
        # #             props_saved['subtype'] = PolygonType.CLOSED
        # #             props_saved['vertices'] = [(int(path.elementAt(i).x), int(path.elementAt(i).y)) for i in range(path.elementCount()-1)]
        # #         else:
        # #             props_saved['subtype'] = PolygonType.OPEN
        # #             props_saved['vertices'] = [(int(path.elementAt(i).x), int(path.elementAt(i).y)) for i in range(path.elementCount())]
        # #
        # #         label_pos = props['labelTextArtist'].scenePos()
        # #         props_saved['labelPos'] = (label_pos.x(), label_pos.y())
        # #
        # #         props_saved.pop('vertexCircles')
        # #         props_saved.pop('labelTextArtist')
        # #
        # #         accepted_proposal_props.append(props_saved)
        #
        #     # print '#############'
        #     # print accepted_proposal_props
        #
        #     # labeling_path = self.dms[sec].save_annotation(accepted_proposal_props, self.username, timestamp)
        #     labeling_path = DataManager.save_annotation(accepted_proposal_props, self.stack, sec, self.username, timestamp,
        #     annotation_rootdir=annotation_midbrainIncluded_rootdir)
        #
        #     # print self.new_labelnames
        #     self.dms[sec].add_labelnames(self.new_labelnames, os.environ['REPO_DIR']+'/gui/newStructureNames.txt')
        #
        #     self.statusBar().showMessage('Labelings saved to %s' % (self.username+'_'+timestamp))
        #
        #     if sec in self.gscenes:
        #         pix = QPixmap(self.dms[sec].image_width/8, self.dms[sec].image_height/8)
        #         painter = QPainter(pix)
        #
        #         self.gscenes[sec].render(painter, QRectF(0,0,self.dms[sec].image_width/8, self.dms[sec].image_height/8),
        #                                 QRectF(0,0,self.dms[sec].image_width, self.dms[sec].image_height))
        #         pix.save(labeling_path[:-4] + '.jpg', "JPG")
        #         print 'Preview image saved to', labeling_path[:-4] + '.jpg'
        #         del painter
        #         del pix


    def load(self):
        # self.gscenes['sagittal'].load_drawings(username='Lauren', timestamp='latest', annotation_rootdir=annotation_midbrainIncluded_v2_rootdir)
        self.gscenes['sagittal'].load_drawings(username='yuncong', timestamp='latest', annotation_rootdir=annotation_midbrainIncluded_v2_rootdir)

    @pyqtSlot()
    def active_image_updated(self):
        self.setWindowTitle('BrainLabelingGUI, stack %(stack)s, section %(sec)d, z=%(z).2f, x=%(x).2f, y=%(y).2f' % \
        dict(stack=self.stack, sec=self.gscenes['sagittal'].active_section, z=self.gscenes['sagittal'].active_i, x=self.gscenes['coronal'].active_i, y=self.gscenes['horizontal'].active_i))

    @pyqtSlot(int, int, int, str)
    def crossline_updated(self, cross_x_lossless, cross_y_lossless, cross_z_lossless, source_gscene_id):
        print 'GUI: update all crosses to', cross_x_lossless, cross_y_lossless, cross_z_lossless, 'from', source_gscene_id

        for gscene_id, gscene in self.gscenes.iteritems():
            # if gscene_id != source_gscene_id:
            #     gscene.update_cross(cross_x_lossless, cross_y_lossless, cross_z_lossless)
            gscene.update_cross(cross_x_lossless, cross_y_lossless, cross_z_lossless)

    @pyqtSlot(object)
    def drawings_updated(self, polygon):
        print 'Drawings updated.'
        self.save()

    @pyqtSlot(object)
    def update_structure_volume_requested(self, polygon):

        name_u = polygon.label
        downsample = polygon.gscene.data_feeder.downsample

        matched_polygons = [p for i, polygons in polygon.gscene.drawings.iteritems() for p in polygons if p.label == name_u]

        if len(matched_polygons) < 2:
            return

        # NOTICE THE reconstructed VOLUME IS DOWNSAMPLED BY this number !!!!
        self.volume_downsample_factor = max(8, np.min([gscene.data_feeder.downsample for gscene in self.gscenes.values()]))
        contour_points_grouped_by_pos = {p.position*downsample/self.volume_downsample_factor: \
                                        [(c.scenePos().x()*downsample/self.volume_downsample_factor,
                                        c.scenePos().y()*downsample/self.volume_downsample_factor)
                                        for c in p.vertex_circles] for p in matched_polygons if p.type != 'interpolated'}

        print contour_points_grouped_by_pos.keys()

        if polygon.gscene.data_feeder.orientation == 'sagittal':
            volume, bbox = interpolate_contours_to_volume(contour_points_grouped_by_pos, 'z')
        elif polygon.gscene.data_feeder.orientation == 'coronal':
            volume, bbox = interpolate_contours_to_volume(contour_points_grouped_by_pos, 'x')
        elif polygon.gscene.data_feeder.orientation == 'horizontal':
            volume, bbox = interpolate_contours_to_volume(contour_points_grouped_by_pos, 'y')

        self.structure_volumes[name_u] = (volume, bbox)

        self.gscenes['coronal'].update_drawings_from_structure_volume(name_u)
        self.gscenes['horizontal'].update_drawings_from_structure_volume(name_u)
        self.gscenes['sagittal'].update_drawings_from_structure_volume(name_u)

        print '3D structure updated.'
        self.statusBar().showMessage('3D structure updated.')


    def eventFilter(self, obj, event):
        # print obj.metaObject().className(), event.type()

        if event.type() == QEvent.KeyPress:
            key = event.key()
            if key == Qt.Key_1:
                self.gscenes['sagittal'].show_previous()
            elif key == Qt.Key_2:
                self.gscenes['sagittal'].show_next()
            elif key == Qt.Key_3:
                self.gscenes['coronal'].show_previous()
            elif key == Qt.Key_4:
                self.gscenes['coronal'].show_next()
            elif key == Qt.Key_5:
                self.gscenes['horizontal'].show_previous()
            elif key == Qt.Key_6:
                self.gscenes['horizontal'].show_next()

            elif key == Qt.Key_Space:
                if not event.isAutoRepeat():
                    for gscene in self.gscenes.itervalues():
                        gscene.set_mode('crossline')

        elif event.type() == QEvent.KeyRelease:
            key = event.key()
            if key == Qt.Key_Space:
                if not event.isAutoRepeat():
                    for gscene in self.gscenes.itervalues():
                        gscene.set_mode('idle')

        return False


def load_structure_names(fn):
    names = {}
    with open(fn, 'r') as f:
        for ln in f.readlines():
            abbr, fullname = ln.split('\t')
            names[abbr] = fullname.strip()
    return names


if __name__ == "__main__":

    import argparse
    import sys
    import time

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Launch brain labeling GUI.')

    parser.add_argument("stack_name", type=str, help="stack name")
    parser.add_argument("-f", "--first_sec", type=int, help="first section")
    parser.add_argument("-l", "--last_sec", type=int, help="last section")
    parser.add_argument("-d", "--downsample", type=int, help="downsample", default=1)
    args = parser.parse_args()

    from sys import argv, exit
    appl = QApplication(argv)

    stack = args.stack_name
    downsample = args.downsample

    first_sec = section_range_lookup[stack][0] if args.first_sec is None else args.first_sec
    last_sec = section_range_lookup[stack][1] if args.last_sec is None else args.last_sec

    m = BrainLabelingGUI(stack=stack, first_sec=first_sec, last_sec=last_sec, downsample=downsample)

    m.showMaximized()
    m.raise_()
    exit(appl.exec_())
