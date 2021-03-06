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

from matplotlib.backend_bases import key_press_handler, MouseEvent, KeyEvent
from matplotlib.backends.backend_qt4agg import (
	FigureCanvasQTAgg as FigureCanvas,
	NavigationToolbar2QT as NavigationToolbar)
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

from ui_BrainLabelingGui_v10 import Ui_BrainLabelingGui

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, PathPatch
from matplotlib.colors import ListedColormap, NoNorm, ColorConverter
from matplotlib.path import Path
from matplotlib.text import Text

from shapely.geometry import Polygon as ShapelyPolygon
from shapely.geometry import Point as ShapelyPoint
from shapely.geometry import LineString as ShapelyLineString
from shapely.geometry import LinearRing as ShapelyLineRing

from skimage.color import label2rgb

from visualization_utilities import *

sys.path.append(os.environ['LOCAL_REPO_DIR'] + '/utilities')
from utilities2015 import *

from collections import defaultdict, OrderedDict, deque
from scipy.spatial.distance import cdist

from operator import attrgetter

import requests

from enum import Enum
class Mode(Enum):
	# PLACING_VERTICES = 'adding vertices'
	REVIEW_PROPOSAL = 'review proposal'
	SELECTING_ROI = 'selecting roi'
	# SELECTING_CONNECTION_TARGET = 'selecting the second vertex to connect'
	MOVING_POLYGON = 'moving polygon'
	MOVING_VERTEX = 'moving vertex'
	CREATING_NEW_POLYGON = 'create new polygon'
	ADDING_VERTICES_CONSECUTIVELY = 'adding vertices consecutively'
	ADDING_VERTICES_RANDOMLY = 'adding vertices randomly'
	KEEP_SELECTION = 'keep selection'

class ProposalType(Enum):
	GLOBAL = 'global'
	LOCAL = 'local'
	FREEFORM = 'freeform'
	ALGORITHM = 'algorithm'

class PolygonType(Enum):
	CLOSED = 'closed'
	OPEN = 'open'
	TEXTURE = 'textured'
	TEXTURE_WITH_CONTOUR = 'texture with contour'
	DIRECTION = 'directionality'

SELECTED_POLYGON_LINEWIDTH = 5
UNSELECTED_POLYGON_LINEWIDTH = 3
SELECTED_CIRCLE_SIZE = 30
UNSELECTED_CIRCLE_SIZE = 5
CIRCLE_PICK_THRESH = 1000.
PAN_THRESHOLD = 10

HISTORY_LEN = 5

AUTO_EXTEND_VIEW_TOLERANCE = 200

class ListSelection(QDialog):
	def __init__(self, item_ls, parent=None):
		super(ListSelection, self).__init__(parent)

		self.setWindowTitle('Detect which landmarks ?')

		self.selected = set([])

		self.listWidget = QListWidget()
		for item in item_ls:    
			w_item = QListWidgetItem(item)
			self.listWidget.addItem(w_item)

			w_item.setFlags(w_item.flags() | Qt.ItemIsUserCheckable)
			w_item.setCheckState(False)

		self.listWidget.itemChanged.connect(self.OnSingleClick)

		layout = QGridLayout()
		layout.addWidget(self.listWidget,0,0,1,3)

		self.but_ok = QPushButton("OK")
		layout.addWidget(self.but_ok ,1,1)
		self.but_ok.clicked.connect(self.OnOk)

		self.but_cancel = QPushButton("Cancel")
		layout.addWidget(self.but_cancel ,1,2)
		self.but_cancel.clicked.connect(self.OnCancel)

		self.setLayout(layout)
		self.setGeometry(300, 200, 460, 350)

	def OnSingleClick(self, item):
		if not item.checkState():
		#   item.setCheckState(False)
			self.selected = self.selected - {str(item.text())}
		#   print self.selected
		else:
		#   item.setCheckState(True)
			self.selected.add(str(item.text()))

		print self.selected

	def OnOk(self):
		self.close()

	def OnCancel(self):
		self.selected = set([])
		self.close()


class CustomQCompleter(QCompleter):
	# adapted from http://stackoverflow.com/a/26440173
	def __init__(self, *args):#parent=None):
		super(CustomQCompleter, self).__init__(*args)
		self.local_completion_prefix = ""
		self.source_model = None
		self.filterProxyModel = QSortFilterProxyModel(self)
		self.usingOriginalModel = False

	def setModel(self, model):
		self.source_model = model
		self.filterProxyModel = QSortFilterProxyModel(self)
		self.filterProxyModel.setSourceModel(self.source_model)
		super(CustomQCompleter, self).setModel(self.filterProxyModel)
		self.usingOriginalModel = True

	def updateModel(self):
		if not self.usingOriginalModel:
			self.filterProxyModel.setSourceModel(self.source_model)

		pattern = QRegExp(self.local_completion_prefix,
								Qt.CaseInsensitive,
								QRegExp.FixedString)

		self.filterProxyModel.setFilterRegExp(pattern)

	def splitPath(self, path):
		self.local_completion_prefix = path
		self.updateModel()
		if self.filterProxyModel.rowCount() == 0:
			self.usingOriginalModel = False
			self.filterProxyModel.setSourceModel(QStringListModel([path]))
			return [path]

		return []

class AutoCompleteComboBox(QComboBox):
	# adapted from http://stackoverflow.com/a/26440173
	def __init__(self, labels, *args, **kwargs):
		super(AutoCompleteComboBox, self).__init__(*args, **kwargs)

		self.setEditable(True)
		self.setInsertPolicy(self.NoInsert)

		self.comp = CustomQCompleter(self)
		self.comp.setCompletionMode(QCompleter.PopupCompletion)
		self.setCompleter(self.comp)#
		self.setModel(labels)

		self.clearEditText()

	def setModel(self, strList):
		self.clear()
		self.insertItems(0, strList)
		self.comp.setModel(self.model())

	def focusInEvent(self, event):
		# self.clearEditText()
		super(AutoCompleteComboBox, self).focusInEvent(event)

	def keyPressEvent(self, event):
		key = event.key()
		if key == 16777220:
			# Enter (if event.key() == QtCore.Qt.Key_Enter) does not work
			# for some reason

			# make sure that the completer does not set the
			# currentText of the combobox to "" when pressing enter
			text = self.currentText()
			self.setCompleter(None)
			self.setEditText(text)
			self.setCompleter(self.comp)

		return super(AutoCompleteComboBox, self).keyPressEvent(event)

class AutoCompleteInputDialog(QDialog):

	def __init__(self, labels, *args, **kwargs):
		super(AutoCompleteInputDialog, self).__init__(*args, **kwargs)
		self.comboBox = AutoCompleteComboBox(parent=self, labels=labels)
		va = QVBoxLayout(self)
		va.addWidget(self.comboBox)
		box = QWidget(self)
		ha = QHBoxLayout(self)
		va.addWidget(box)
		box.setLayout(ha)
		self.OK = QPushButton("OK", self)
		self.OK.setDefault(True)
		# cancel = QPushButton("Cancel", self)
		ha.addWidget(self.OK)
		# ha.addWidget(cancel)

	def set_test_callback(self, callback):
		self.OK.clicked.connect(callback)
		# OK.clicked.connect(self.accept)
		# cancel.clicked.connect(self.reject)

class BrainLabelingGUI(QMainWindow, Ui_BrainLabelingGui):
	def __init__(self, parent=None, stack=None):
		"""
		Initialization of BrainLabelingGUI.
		"""

		self.params_dir = '../params'

		# self.app = QApplication(sys.argv)
		QMainWindow.__init__(self, parent)

		# self.init_data(stack)
		self.stack = stack
		self.initialize_brain_labeling_gui()

		self.recent_labels = []

		self.history = deque(maxlen=HISTORY_LEN)

	def init_data(self, section):

		self.section = section
		self.section2 = section + 1

		self.dm = DataManager(
			data_dir=os.environ['LOCAL_DATA_DIR'], 
				 repo_dir=os.environ['LOCAL_REPO_DIR'], 
				 result_dir=os.environ['LOCAL_RESULT_DIR'], 
				 labeling_dir=os.environ['LOCAL_LABELING_DIR'],
			stack=self.stack, section=section, segm_params_id='tSLIC200')

		# self.dm2 = DataManager(
		#     data_dir=os.environ['LOCAL_DATA_DIR'], 
		#          repo_dir=os.environ['LOCAL_REPO_DIR'], 
		#          result_dir=os.environ['LOCAL_RESULT_DIR'], 
		#          labeling_dir=os.environ['LOCAL_LABELING_DIR'],
		#     stack=self.stack, section=section+1, segm_params_id='tSLIC200')
		# self.dm2._load_image(versions=['rgb-jpg'])

		# self.dms = dict([(i, DataManager(
		#     data_dir=os.environ['LOCAL_DATA_DIR'], 
		#          repo_dir=os.environ['LOCAL_REPO_DIR'], 
		#          result_dir=os.environ['LOCAL_RESULT_DIR'], 
		#          labeling_dir=os.environ['LOCAL_LABELING_DIR'],
		#     stack=stack, section=i, segm_params_id='tSLIC200')) for i in range(self.section-5, self.section) + range(self.section+1, self.section+6)])

		self.dms = dict([(section+1, DataManager(
			data_dir=os.environ['LOCAL_DATA_DIR'], 
				 repo_dir=os.environ['LOCAL_REPO_DIR'], 
				 result_dir=os.environ['LOCAL_RESULT_DIR'], 
				 labeling_dir=os.environ['LOCAL_LABELING_DIR'],
			stack=self.stack, section=section+1, segm_params_id='tSLIC200'))])

		for dm in self.dms.itervalues():
			dm._load_image(versions=['rgb-jpg'])

		self.dm2 = self.dms[self.section2]

		self.new_labelnames = {}
		if os.path.exists(self.dm.repo_dir+'/visualization/newStructureNames.txt'):
			with open(self.dm.repo_dir+'/visualization/newStructureNames.txt', 'r') as f:
				for ln in f.readlines():
					abbr, fullname = ln.split('\t')
					self.new_labelnames[abbr] = fullname.strip()
			self.new_labelnames = OrderedDict(sorted(self.new_labelnames.items()))

		self.structure_names = {}
		with open(self.dm.repo_dir+'/visualization/structure_names.txt', 'r') as f:
			for ln in f.readlines():
				abbr, fullname = ln.split('\t')
				self.structure_names[abbr] = fullname.strip()
		self.structure_names = OrderedDict(self.new_labelnames.items() + sorted(self.structure_names.items()))

		print self.dm.slice_ind

		t = time.time()
		self.dm._load_image(versions=['rgb-jpg'])
		print 1, time.time() - t

		required_results = [
		'segmentationTransparent', 
		'segmentation',
		'allSeedClusterScoreDedgeTuples',
		'proposals',
		'spCoveredByProposals',
		'edgeMidpoints',
		'edgeEndpoints',
		'spAreas',
		'spBbox',
		'spCentroids'
		]

		# t = time.time()
		# self.dm.download_results(required_results)
		# print 2, time.time() - t

		# t = time.time()
		# self.dm.load_multiple_results(required_results, download_if_not_exist=True)
		# print 3, time.time() - t

		self.segm_transparent = None
		self.under_img = None
		self.textonmap_vis = None
		self.dirmap_vis = None

		self.set_mode(Mode.REVIEW_PROPOSAL)

		self.boundary_colors = [(0,1,1), (1,0,0), (0,0,0),(0,0,1)] # unknown, accepted, rejected

		self.accepted_proposals = defaultdict(dict)

		self.selected_proposal_polygon = None
		self.alternative_global_proposal_ind = 0
		self.alternative_local_proposal_ind = 0

		# self.just_added_vertex = False

	def paramSettings_clicked(self):
		pass

	def openMenu(self, canvas_pos):

		self.deleteVertex_Action.setVisible(self.selected_circle is not None)

		self.movePolygon_Action.setVisible(self.selected_proposal_polygon is not None and \
			self.selected_proposal_polygon in self.accepted_proposals)

		# self.moveVertex_Action.setVisible(self.selected_proposal_polygon is not None and \
		#     self.selected_proposal_polygon in self.accepted_proposals)

		self.addVertex_Action.setVisible(self.selected_proposal_polygon is not None and \
			self.selected_proposal_polygon in self.accepted_proposals and self.mode != Mode.ADDING_VERTICES_RANDOMLY)

		self.extendPolygon_Action.setVisible(self.selected_circle is not None and \
			(self.selected_vertex_index == len(self.accepted_proposals[self.selected_proposal_polygon]['vertexPatches']) - 1 or \
			self.selected_vertex_index == 0))

		self.doneAddingVertex_Action.setVisible(self.mode == Mode.ADDING_VERTICES_CONSECUTIVELY or self.mode == Mode.ADDING_VERTICES_RANDOMLY)

		self.deleteVerticesROI_Action.setVisible(hasattr(self, 'roi_rectPatch') and self.roi_rectPatch is not None)
		self.deleteVerticesROIOpen_Action.setVisible(hasattr(self, 'roi_rectPatch') and self.roi_rectPatch is not None)
				
		self.breakEdge_Action.setVisible(self.selected_proposal_polygon is not None and \
			self.selected_proposal_polygon in self.accepted_proposals)

		self.newPolygon_Action.setVisible(self.mode == Mode.REVIEW_PROPOSAL)

		self.accProp_Action.setVisible(self.selected_proposal_polygon is not None and self.selected_proposal_polygon not in self.accepted_proposals)
		self.rejProp_Action.setVisible(self.selected_proposal_polygon is not None and self.selected_proposal_polygon in self.accepted_proposals)
		self.changeLabel_Action.setVisible(self.selected_proposal_polygon is not None and self.selected_proposal_polygon in self.accepted_proposals)

		print 'open menu'

		action = self.menu.exec_(self.cursor().pos())

		print 'action served'

		if action == self.doneAddingVertex_Action:

			self.set_mode(Mode.REVIEW_PROPOSAL)

			if 'label' not in self.accepted_proposals[self.selected_proposal_polygon]:
				self.acceptProposal_callback()

		elif action == self.deleteVertex_Action:
			self.remove_selected_vertex()

		# elif action == self.connectTo_Action:
		#     self.set_mode(Mode.SELECTING_CONNECTION_TARGET)

		elif action == self.selectROI_Action:
			self.set_mode(Mode.SELECTING_ROI)

		elif action == self.deleteVerticesROI_Action:
			self.remove_selected_vertices_in_region(link_endpoints=True)

		elif action == self.deleteVerticesROIOpen_Action:
			self.remove_selected_vertices_in_region(link_endpoints=False)
	
		# elif action == self.addVertex_Action:
		#     # self.add_vertex(canvas_pos)
		#     self.set_mode(Mode.ADDING_VERTICES_RANDOMLY)

		# elif action == self.keepSelection_Action:
		#     # self.add_vertex(canvas_pos)
		#     self.set_mode(Mode.KEEP_SELECTION)

		elif action == self.breakEdge_Action:
			self.break_edge(canvas_pos)

		elif action == self.accProp_Action:
			self.acceptProposal_callback()

		elif action == self.rejProp_Action:
			self.rejectProposal_callback()

		elif action == self.newPolygon_Action:
			
			self.cancel_current_selection()

			self.set_mode(Mode.CREATING_NEW_POLYGON)
			self.history.clear()

		elif action == self.changeLabel_Action:
			self.open_label_selection_dialog()

		elif action == self.movePolygon_Action:
			self.set_mode(Mode.MOVING_POLYGON)

		elif action == self.extendPolygon_Action:
			self.set_mode(Mode.ADDING_VERTICES_CONSECUTIVELY)
			
			if self.selected_vertex_index == 0:
				self.extend_head = True

			# self.previous_circle = self.selected_circle
			# self.previous_vertex_index = self.selected_vertex_index
			# self.previous_polygon = self.selected_proposal_polygon

			# print 'prev is set to', self.previous_polygon, self.selected_vertex_index

			# self.cancel_current_circle()
		
		elif action == self.addVertex_Action:
			self.set_mode(Mode.ADDING_VERTICES_RANDOMLY)

		else:
			pass


	def reload_brain_labeling_gui(self):

		self.picked_artists = []

		self.extend_head = False
		self.connecting_vertices = False

		self.seg_loaded = False
		self.superpixels_on = False
		self.labels_on = True
		self.contours_on = True
		self.object_picked = False
		self.shuffle_global_proposals = True # instead of local proposals

		self.base_scale = 1.2       # multiplication factor for zoom using scroll wheel
		self.pressed = False           # related to pan (press and drag) vs. select (click)
		
		self.selected_circle = None
		# self.selected_proposal_vertices = []
		# self.selected_proposal_vertexCircles = []

		self.setWindowTitle('BrainLabelingGUI, stack %s'%self.stack + ', section %d' %self.section)

		if hasattr(self, 'axis'):
			import copy
			for p in copy.copy(self.axis.patches):
				p.remove()
			for p in copy.copy(self.axis.artists):
				p.remove()
		else:
			self.axis = self.fig.add_subplot(111)
			self.axis.axis('off')

			self.axis2 = self.fig2.add_subplot(111)
			self.axis2.axis('off')

		t = time.time()

		if hasattr(self, 'orig_image_handle'):
			self.orig_image_handle.set_data(self.dm.image_rgb_jpg)
			self.orig_image_handle2.set_data(self.dm2.image_rgb_jpg)

			self.section_text.set_text(str(self.section))
			self.section_text2.set_text(str(self.section2))
		else:

			img = self.dm.image_rgb_jpg
			img[~self.dm.mask] = 0

			img2 = self.dm2.image_rgb_jpg
			img2[~self.dm2.mask] = 0

			self.orig_image_handle = self.axis.imshow(img, cmap=plt.cm.Greys_r, aspect='equal')
			self.orig_image_handle2 = self.axis2.imshow(img2, cmap=plt.cm.Greys_r, aspect='equal')

			# self.orig_image_handle = self.axis.imshow(self.dm.image_rgb_jpg, cmap=plt.cm.Greys_r, aspect='equal')
			# self.orig_image_handle2 = self.axis2.imshow(self.dm2.image_rgb_jpg, cmap=plt.cm.Greys_r, aspect='equal')

			self.section_text = self.axis.text(0.01, 0.99, str(self.section), transform=self.fig.transFigure, fontsize=16, fontweight='bold', va='top', color=(1,0,0))
			self.section_text2 = self.axis2.text(0.01, 0.99, str(self.section2), transform=self.fig2.transFigure, fontsize=16, fontweight='bold', va='top', color=(1,0,0))

			fig_right_on_display, fig_top_on_display = self.fig.transFigure.transform((1,1))
			fig_left_on_display, fig_bottom_on_display = self.fig.transFigure.transform((0,0))
			fig_width_on_display = fig_right_on_display - fig_left_on_display
			fig_height_on_display = fig_top_on_display - fig_bottom_on_display

			self.axis.set_position([0,0,1,1])
			self.axis.set_ylim([self.dm.image_height, 0])
			self.axis.set_xlim([0, fig_width_on_display * float(self.dm.image_height) / fig_height_on_display])

			self.axis2.set_position([0,0,1,1])
			self.axis2.set_ylim([self.dm.image_height, 0])
			self.axis2.set_xlim([0, fig_width_on_display * float(self.dm.image_height) / fig_height_on_display])

		print 4, time.time() - t

		self.setWindowTitle('BrainLabelingGUI, stack %s'%self.stack + ', section %d' %self.section + ', section %d'%self.section2)

		self.newxmin, self.newxmax = self.axis.get_xlim()
		self.newymin, self.newymax = self.axis.get_ylim()

		self.canvas.draw()
		self.canvas2.draw()

		self.show()

	def initialize_brain_labeling_gui(self):

		self.menu = QMenu()

		self.addVertex_Action = self.menu.addAction("Insert vertices")
		self.extendPolygon_Action = self.menu.addAction("Extend from this vertex")
		self.doneAddingVertex_Action = self.menu.addAction("Done adding vertex")

		self.movePolygon_Action = self.menu.addAction("Move")

		self.deleteVertex_Action = self.menu.addAction("Delete vertex")

		self.deleteVerticesROI_Action = self.menu.addAction("Delete vertices in ROI (close)")
		self.deleteVerticesROIOpen_Action = self.menu.addAction("Delete vertices in ROI (open)")

		self.selectROI_Action = self.menu.addAction("Select ROI")

		self.breakEdge_Action = self.menu.addAction("break edge")

		self.newPolygon_Action = self.menu.addAction("New polygon")

		self.accProp_Action = self.menu.addAction("Accept")
		self.rejProp_Action = self.menu.addAction("Reject")

		self.changeLabel_Action = self.menu.addAction('Change label')

		# A set of high-contrast colors proposed by Green-Armytage
		self.colors = np.loadtxt('100colors.txt', skiprows=1)
		self.label_cmap = ListedColormap(self.colors, name='label_cmap')

		self.setupUi(self)

		###########################

		self.fig = self.canvaswidget.fig
		self.canvas = self.canvaswidget.canvas

		self.canvas.setFocusPolicy( Qt.ClickFocus )
		self.canvas.setFocus()

		self.canvas.mpl_connect('scroll_event', self.on_zoom)
		self.bpe_id = self.canvas.mpl_connect('button_press_event', self.on_press)
		self.bre_id = self.canvas.mpl_connect('button_release_event', self.on_release)
		self.canvas.mpl_connect('motion_notify_event', self.on_motion)
		self.canvas.mpl_connect('key_press_event', self.on_key_press)
		self.canvas.mpl_connect('key_release_event', self.on_key_release)
		self.canvas.mpl_connect('pick_event', self.on_pick)

		###########################

		self.fig2 = self.canvaswidget2.fig
		self.canvas2 = self.canvaswidget2.canvas

		self.canvas2.setFocusPolicy( Qt.ClickFocus )
		self.canvas2.setFocus()

		self.canvas2.mpl_connect('scroll_event', self.on_zoom)
		self.bpe_id = self.canvas2.mpl_connect('button_press_event', self.on_press)
		self.bre_id = self.canvas2.mpl_connect('button_release_event', self.on_release)
		self.canvas2.mpl_connect('motion_notify_event', self.on_motion)
		self.canvas2.mpl_connect('key_press_event', self.on_key_press)
		self.canvas2.mpl_connect('pick_event', self.on_pick)

		###########################

		self.button_autoDetect.clicked.connect(self.autoDetect_callback)
		self.button_updateDB.clicked.connect(self.updateDB_callback)
		self.button_loadLabeling.clicked.connect(self.load_callback)
		self.button_saveLabeling.clicked.connect(self.save_callback)
		self.button_quit.clicked.connect(self.close)
		
		self.display_buttons = [self.img_radioButton, self.textonmap_radioButton, self.dirmap_radioButton]
		self.img_radioButton.setChecked(True)

		for b in self.display_buttons:
			b.toggled.connect(self.display_option_changed)

		self.radioButton_globalProposal.toggled.connect(self.mode_changed)
		self.radioButton_localProposal.toggled.connect(self.mode_changed)

		self.buttonSpOnOff.clicked.connect(self.display_option_changed)
		self.button_labelsOnOff.clicked.connect(self.toggle_labels)
		self.button_contoursOnOff.clicked.connect(self.toggle_contours)

		# self.thumbnail_list = QListWidget(parent=self)
		self.thumbnail_list.setIconSize(QSize(200,200))
		self.thumbnail_list.setResizeMode(QListWidget.Adjust)
		self.thumbnail_list.itemDoubleClicked.connect(self.section_changed)


		section_range_lookup = {'MD593': (41,176), 'MD594': (47,186), 'MD595': (35,164), 'MD592': (46,185), 'MD589':(49,186)}
		first_sec, last_sec = section_range_lookup[self.stack]
		for i in range(first_sec, last_sec):
			item = QListWidgetItem(QIcon("/home/yuncong/CSHL_data_processed/%(stack)s_thumbnail_aligned_cropped/%(stack)s_%(sec)04d_thumbnail_aligned_cropped.tif"%{'sec':i, 'stack': self.stack}), str(i))
			# item.setFont(QFont())
			self.thumbnail_list.addItem(item)

		self.thumbnail_list.resizeEvent = self.thumbnail_list_resized
		self.init_thumbnail_list_width = self.thumbnail_list.width()
		# print self.init_thumbnail_list_width

	def thumbnail_list_resized(self, event):
		new_size = 200 * event.size().width() / self.init_thumbnail_list_width
		self.thumbnail_list.setIconSize( QSize(new_size , new_size ) )

	def toggle_labels(self):

		self.labels_on = not self.labels_on

		if not self.labels_on:

			for patch, props in self.accepted_proposals.iteritems():
				props['labelTextArtist'].remove()

			self.button_labelsOnOff.setText('Turns Labels ON')

		else:
			for patch, props in self.accepted_proposals.iteritems():
				self.axis.add_artist(props['labelTextArtist'])

			self.button_labelsOnOff.setText('Turns Labels OFF')

		self.canvas.draw()

	def toggle_contours(self):

		self.contours_on = not self.contours_on

		if not self.contours_on:

			for patch, props in self.accepted_proposals.iteritems():
				patch.remove()
				for circ in props['vertexPatches']:
					circ.remove()

			self.button_contoursOnOff.setText('Turns Contours ON')

		else:
			for patch, props in self.accepted_proposals.iteritems():
				self.axis.add_patch(patch)
				for circ in props['vertexPatches']:
					self.axis.add_patch(circ)

			self.button_contoursOnOff.setText('Turns Contours OFF')

		self.canvas.draw()

	def updateDB_callback(self):
		cmd = 'rsync -az --include="*/" %(local_labeling_dir)s/%(stack)s yuncong@gcn-20-33.sdsc.edu:%(gordon_labeling_dir)s' % {'gordon_labeling_dir':os.environ['GORDON_LABELING_DIR'],
																			'local_labeling_dir':os.environ['LOCAL_LABELING_DIR'],
																			'stack': self.stack
																			}
		os.system(cmd)

		# cmd = 'rsync -az %(local_labeling_dir)s/labelnames.txt yuncong@gcn-20-33.sdsc.edu:%(gordon_labeling_dir)s' % {'gordon_labeling_dir':os.environ['GORDON_LABELING_DIR'],
		#                                                             'local_labeling_dir':os.environ['LOCAL_LABELING_DIR'],
		#                                                             }
		# os.system(cmd)
		self.statusBar().showMessage('labelings synced')

		# payload = {'section': self.dm.slice_ind}
		# r = requests.get('http://gcn-20-32.sdsc.edu:5000/update_db', params=payload)
		r = requests.get('http://gcn-20-32.sdsc.edu:5000/update_db')
		res = r.json()
		if res['result'] == 0:
			self.statusBar().showMessage('Landmark database updated')


	def detect_landmark(self, labels):

		payload = {'labels': labels, 'section': self.dm.slice_ind}
		r = requests.get('http://gcn-20-32.sdsc.edu:5000/top_down_detect', params=payload)
		print r.url
		return r.json()

	def autoDetect_callback(self):
		self.labelsToDetect = ListSelection([abbr + ' (' + fullname + ')' for abbr, fullname in self.structure_names.iteritems()], parent=self)
		self.labelsToDetect.exec_()

		if len(self.labelsToDetect.selected) > 0:
		
			returned_alg_proposal_dict = self.detect_landmark([x.split()[0] for x in list(self.labelsToDetect.selected)]) 
			# list of tuples (sps, dedges, sig)

			for label, (sps, dedges, sig) in returned_alg_proposal_dict.iteritems():

				props = {}

				props['vertices'] = self.dm.vertices_from_dedges(dedges)
				patch = Polygon(props['vertices'], closed=True, edgecolor=self.boundary_colors[0], fill=False, linewidth=UNSELECTED_POLYGON_LINEWIDTH)
				patch.set_picker(True)
				self.axis.add_patch(patch)

				props['vertexPatches'] = []
				for x,y in props['vertices']:
					vertex_circle = plt.Circle((x, y), radius=UNSELECTED_CIRCLE_SIZE, color=self.boundary_colors[1], alpha=.8)
					vertex_circle.set_picker(CIRCLE_PICK_THRESH)
					props['vertexPatches'].append(vertex_circle)
					self.axis.add_patch(vertex_circle)
					vertex_circle.set_picker(True)


				centroid = np.mean(props['vertices'], axis=0)
				props['labelTextArtist'] = Text(centroid[0], centroid[1], label, style='italic', bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
				self.axis.add_artist(props['labelTextArtist'])
				props['labelTextArtist'].set_picker(True)

				self.accepted_proposals[patch] = props

				props['sps'] = sps
				props['dedges'] = dedges
				props['sig'] = sig
				props['type'] = ProposalType.ALGORITHM
				props['label'] = label
		
		self.canvas.draw()

	def on_pick(self, event):

		self.object_picked = False

		if event.mouseevent.name == 'scroll_event':
			return

		print 'pick callback triggered'

		self.picked_artists.append(event.artist)


	def handle_picked_object(self, artist):

		matches = [patch for patch, props in self.accepted_proposals.iteritems() if 'labelTextArtist' in props and artist == props['labelTextArtist']]
		
		if len(matches) > 0:

			print 'clicked on a label'
			self.selected_proposal_polygon = matches[0]
			self.selected_proposal_polygon.set_linewidth(SELECTED_POLYGON_LINEWIDTH)

			if 'label' in self.accepted_proposals[self.selected_proposal_polygon]:
				self.statusBar().showMessage('picked %s proposal (%s, %s)' % (self.accepted_proposals[self.selected_proposal_polygon]['type'].value,
																	 self.accepted_proposals[self.selected_proposal_polygon]['label'],
																	 self.structure_names[self.accepted_proposals[self.selected_proposal_polygon]['label']]))

		else:

			patch_vertexInd_tuple = [(patch, props['vertexPatches'].index(artist)) for patch, props in self.accepted_proposals.iteritems() 
						if 'vertexPatches' in props and artist in props['vertexPatches']]
			
			if len(patch_vertexInd_tuple) == 1:

				# if self.mode == Mode.ADDING_VERTICES_CONSECUTIVELY or self.mode == Mode.ADDING_VERTICES_RANDOMLY:
				#     if patch_vertexInd_tuple[0][0] != self.selected_proposal_polygon:
				#         return

				self.cancel_current_selection()

				print 'clicked on a vertex circle'
				self.selected_proposal_polygon = patch_vertexInd_tuple[0][0]
				self.selected_vertex_index = patch_vertexInd_tuple[0][1]

				print 'vertex index', self.selected_vertex_index

				self.selected_circle = artist
				self.selected_circle.set_radius(SELECTED_CIRCLE_SIZE)
				
				self.selected_proposal_polygon.set_linewidth(SELECTED_POLYGON_LINEWIDTH)

				if 'label' in self.accepted_proposals[self.selected_proposal_polygon]:
					self.statusBar().showMessage('picked %s proposal (%s, %s), vertex %d' % (self.accepted_proposals[self.selected_proposal_polygon]['type'].value,
																		 self.accepted_proposals[self.selected_proposal_polygon]['label'],
																		 self.structure_names[self.accepted_proposals[self.selected_proposal_polygon]['label']],
																		 self.selected_vertex_index))

				self.selected_polygon_circle_centers_before_drag = [circ.center for circ in self.accepted_proposals[self.selected_proposal_polygon]['vertexPatches']]

			elif len(patch_vertexInd_tuple) == 0 and self.selected_circle is None:

				print 'clicked on a polygon'

				if self.mode == Mode.ADDING_VERTICES_CONSECUTIVELY or self.mode == Mode.ADDING_VERTICES_RANDOMLY:
					return
				else:
					self.cancel_current_selection()

					if artist in self.accepted_proposals:
						print 'this polygon has been accepted'

						self.selected_proposal_polygon = artist
						self.selected_proposal_polygon.set_linewidth(SELECTED_POLYGON_LINEWIDTH)

						# self.selected_polygon_xy_before_drag = self.selected_proposal_polygon.get_xy()
						self.selected_polygon_circle_centers_before_drag = [circ.center 
											for circ in self.accepted_proposals[self.selected_proposal_polygon]['vertexPatches']]

						self.selected_polygon_label_pos_before_drag = self.accepted_proposals[self.selected_proposal_polygon]['labelTextArtist'].get_position()

						if 'label' in self.accepted_proposals[self.selected_proposal_polygon]:
							self.statusBar().showMessage('picked %s proposal (%s, %s)' % (self.accepted_proposals[self.selected_proposal_polygon]['type'].value,
																				 self.accepted_proposals[self.selected_proposal_polygon]['label'],
																				 self.structure_names[self.accepted_proposals[self.selected_proposal_polygon]['label']]))

						self.selected_proposal_type = self.accepted_proposals[self.selected_proposal_polygon]['type']

		self.object_picked = True

		self.canvas.draw()

	def load_local_proposals(self):

		sys.stderr.write('loading local proposals ...\n')
		self.statusBar().showMessage('loading local proposals ...')
		
		cluster_tuples = self.dm.load_pipeline_result('allSeedClusterScoreDedgeTuples')
		self.local_proposal_tuples = [(cl, ed, sig) for seed, cl, sig, ed in cluster_tuples]
		self.local_proposal_clusters = [m[0] for m in self.local_proposal_tuples]
		self.local_proposal_dedges = [m[1] for m in self.local_proposal_tuples]
		self.local_proposal_sigs = [m[2] for m in self.local_proposal_tuples]

		self.n_local_proposals = len(self.local_proposal_tuples)
		
		if not hasattr(self, 'local_proposal_pathPatches'):
			self.local_proposal_pathPatches = [None] * self.n_local_proposals
			self.local_proposal_vertexCircles = [None] * self.n_local_proposals

		self.local_proposal_indices_from_sp = defaultdict(list)
		for i, (seed, _, _, _) in enumerate(cluster_tuples):
			self.local_proposal_indices_from_sp[seed].append(i)
		self.local_proposal_indices_from_sp.default_factory = None

		sys.stderr.write('%d local proposals loaded.\n' % self.n_local_proposals)
		self.statusBar().showMessage('Local proposals loaded.')

		self.local_proposal_labels = [None] * self.n_local_proposals


	def load_global_proposals(self):
		
		self.global_proposal_tuples =  self.dm.load_pipeline_result('proposals')
		self.global_proposal_clusters = [m[0] for m in self.global_proposal_tuples]
		self.global_proposal_dedges = [m[1] for m in self.global_proposal_tuples]
		self.global_proposal_sigs = [m[2] for m in self.global_proposal_tuples]

		self.n_global_proposals = len(self.global_proposal_tuples)

		if not hasattr(self, 'global_proposal_pathPatches'):
			self.global_proposal_pathPatches = [None] * self.n_global_proposals
			self.global_proposal_vertexCircles = [None] * self.n_global_proposals

		self.statusBar().showMessage('%d global proposals loaded' % self.n_global_proposals)

		self.sp_covered_by_proposals = self.dm.load_pipeline_result('spCoveredByProposals')
		self.sp_covered_by_proposals = dict([(s, list(props)) for s, props in self.sp_covered_by_proposals.iteritems()])

		self.global_proposal_labels = [None] * self.n_global_proposals

	def load_callback(self):

		fname = str(QFileDialog.getOpenFileName(self, 'Open file', self.dm.labelings_dir))
		stack, sec, username, timestamp, suffix = os.path.basename(fname[:-4]).split('_')

		self.accepted_proposals = {}

		_, _, _, accepted_proposal_props = self.dm.load_proposal_review_result(username, timestamp, suffix)

		for props in accepted_proposal_props:
			if props['subtype'] == PolygonType.CLOSED:
				patch = Polygon(props['vertices'], closed=True, fill=False, edgecolor=self.boundary_colors[1], linewidth=UNSELECTED_POLYGON_LINEWIDTH)
			elif props['subtype'] == PolygonType.OPEN:
				patch = Polygon(props['vertices'], closed=False, fill=False, edgecolor=self.boundary_colors[1], linewidth=UNSELECTED_POLYGON_LINEWIDTH)
			else:
				raise 'unknown polygon type'

			props['vertexPatches'] = []
			for x,y in props['vertices']:
				vertex_circle = plt.Circle((x, y), radius=UNSELECTED_CIRCLE_SIZE, color=self.boundary_colors[1], alpha=.8)
				vertex_circle.set_picker(CIRCLE_PICK_THRESH)
				props['vertexPatches'].append(vertex_circle)
				self.axis.add_patch(vertex_circle)
				vertex_circle.set_picker(True)

			self.axis.add_patch(patch)
			patch.set_picker(True)

			if 'labelPos' not in props:
				centroid = np.mean(props['vertices'], axis=0)
				props['labelTextArtist'] = Text(centroid[0], centroid[1], props['label'], style='italic', bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
			else:
				props['labelTextArtist'] = Text(props['labelPos'][0], props['labelPos'][1], props['label'], style='italic', bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})

			self.axis.add_artist(props['labelTextArtist'])
			props['labelTextArtist'].set_picker(True)

			props.pop('vertices')

			self.accepted_proposals[patch] = props

		self.canvas.draw()


	def open_label_selection_dialog(self):

		print 'open_label_selection_dialog'

		if hasattr(self, 'recent_labels') and self.recent_labels is not None and len(self.recent_labels) > 0:
			self.structure_names = OrderedDict([(abbr, fullname) for abbr, fullname in self.structure_names.iteritems() if abbr in self.recent_labels] + \
							[(abbr, fullname) for abbr, fullname in self.structure_names.iteritems() if abbr not in self.recent_labels])

		self.label_selection_dialog = AutoCompleteInputDialog(parent=self, labels=[abbr + ' (' + fullname + ')' for abbr, fullname in self.structure_names.iteritems()])
		# self.label_selection_dialog = QInputDialog(self)
		self.label_selection_dialog.setWindowTitle('Select landmark label')

		# if hasattr(self, 'invalid_labelname'):
		#     print 'invalid_labelname', self.invalid_labelname
		# else:
		#     print 'no labelname set'

		if 'label' in self.accepted_proposals[self.selected_proposal_polygon]:
			self.label_selection_dialog.comboBox.setEditText(self.accepted_proposals[self.selected_proposal_polygon]['label']+' ('+self.structure_names[self.accepted_proposals[self.selected_proposal_polygon]['label']]+')')
		else:
			self.accepted_proposals[self.selected_proposal_polygon]['label'] = ''

		self.label_selection_dialog.set_test_callback(self.label_dialog_text_changed)

		# self.label_selection_dialog.accepted.connect(self.label_dialog_text_changed)
		# self.label_selection_dialog.textValueSelected.connect(self.label_dialog_text_changed)

		self.label_selection_dialog.exec_()

	def label_dialog_text_changed(self):

		print 'label_dialog_text_changed'

		text = str(self.label_selection_dialog.comboBox.currentText())

		import re
		m = re.match('^(.+?)\s*\((.+)\)$', text)

		if m is None:
			QMessageBox.warning(self, 'oops', 'structure name must be of the form "abbreviation (full description)"')
			return

		else:
			abbr, fullname = m.groups()
			if not (abbr in self.structure_names.keys() and fullname in self.structure_names.values()):  # new label
				if abbr in self.structure_names:
					QMessageBox.warning(self, 'oops', 'structure with abbreviation %s already exists: %s' % (abbr, fullname))
					return
				else:
					self.structure_names[abbr] = fullname
					self.new_labelnames[abbr] = fullname

		self.accepted_proposals[self.selected_proposal_polygon]['label'] = abbr

		if 'labelTextArtist' in self.accepted_proposals[self.selected_proposal_polygon] and self.accepted_proposals[self.selected_proposal_polygon]['labelTextArtist'] is not None:
			self.accepted_proposals[self.selected_proposal_polygon]['labelTextArtist'].set_text(abbr)
		else:
			vertices = self.vertices_from_vertexPatches(self.selected_proposal_polygon)
			centroid = np.mean(vertices, axis=0)
			text_artist = Text(centroid[0], centroid[1], abbr, style='italic', bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
			self.accepted_proposals[self.selected_proposal_polygon]['labelTextArtist'] = text_artist
			self.axis.add_artist(text_artist)
			text_artist.set_picker(True)

		self.recent_labels.insert(0, abbr)
		# self.invalid_labelname = None

		self.label_selection_dialog.accept()


	def acceptProposal_callback(self):

		if self.selected_proposal_type == ProposalType.GLOBAL:
			self.accepted_proposals[self.selected_proposal_polygon] = {'sps': self.global_proposal_clusters[self.selected_proposal_id],
																	'dedges': self.global_proposal_dedges[self.selected_proposal_id],
																	'sig': self.global_proposal_sigs[self.selected_proposal_id],
																	'type': self.selected_proposal_type,
																	'id': self.selected_proposal_id,
																	# 'vertices': self.selected_proposal_polygon.get_xy(),
																	'vertexPatches': self.selected_proposal_vertexCircles}

		elif self.selected_proposal_type == ProposalType.LOCAL:
			self.accepted_proposals[self.selected_proposal_polygon] = {'sps': self.local_proposal_clusters[self.selected_proposal_id],
																	'dedges': self.local_proposal_dedges[self.selected_proposal_id],
																	'sig': self.local_proposal_sigs[self.selected_proposal_id],
																	'type': self.selected_proposal_type,
																	'id': self.selected_proposal_id,
																	# 'vertices': self.selected_proposal_polygon.get_xy(),
																	'vertexPatches': self.selected_proposal_vertexCircles}

		# self.selected_proposal_polygon.set_color(self.boundary_colors[1])
		# self.selected_proposal_polygon.set_picker(True)
		# for circ in self.selected_proposal_vertexCircles:
		# for circ in self.accepted_proposals[self.selected_proposal_polygon]['vertexPatches']:
			# circ.set_color(self.boundary_colors[1])
			# circ.set_picker(True)

		self.canvas.draw()

		self.open_label_selection_dialog()

		self.cancel_current_selection()

		self.save_callback()

		self.history.clear()

	def rejectProposal_callback(self):

		for circ in self.accepted_proposals[self.selected_proposal_polygon]['vertexPatches']:
			circ.remove()

		self.accepted_proposals[self.selected_proposal_polygon]['labelTextArtist'].remove()
		self.selected_proposal_polygon.remove()

		self.accepted_proposals.pop(self.selected_proposal_polygon)

		# self.selected_proposal_polygon.set_color(self.boundary_colors[0])
		# self.selected_proposal_polygon.set_picker(None)

		self.cancel_current_selection()
		
		self.canvas.draw()

	def shuffle_proposal_from_pool(self, sp_ind):

		if self.shuffle_global_proposals:   
			if not hasattr(self, 'sp_covered_by_proposals'):
				return
		else:
			if not hasattr(self, 'local_proposal_indices_from_sp'):
				return

		if self.shuffle_global_proposals:

			if sp_ind not in self.sp_covered_by_proposals or sp_ind == -1:
				self.statusBar().showMessage('No proposal covers superpixel %d' % sp_ind)
				return 
		else:
			if sp_ind == -1:
				return
		
		if self.object_picked:
			return

		self.cancel_current_selection()

		if self.shuffle_global_proposals:
			self.selected_proposal_type = ProposalType.GLOBAL

			self.alternative_global_proposal_ind = (self.alternative_global_proposal_ind + 1) % len(self.sp_covered_by_proposals[sp_ind])
			self.selected_proposal_id = self.sp_covered_by_proposals[sp_ind][self.alternative_global_proposal_ind]

			dedges = self.global_proposal_dedges[self.selected_proposal_id]
		else:

			self.selected_proposal_type = ProposalType.LOCAL

			self.alternative_local_proposal_ind = (self.alternative_local_proposal_ind + 1) % len(self.local_proposal_indices_from_sp[sp_ind])
			self.selected_proposal_id = self.local_proposal_indices_from_sp[sp_ind][self.alternative_local_proposal_ind]

			cl, dedges, sig = self.local_proposal_tuples[self.selected_proposal_id]


		if self.shuffle_global_proposals:
			proposal_pathPatches = self.global_proposal_pathPatches
			proposal_vertexCircles = self.global_proposal_vertexCircles
		else:
			proposal_pathPatches = self.local_proposal_pathPatches
			proposal_vertexCircles = self.local_proposal_vertexCircles

		if proposal_pathPatches[self.selected_proposal_id] is None:  
			vertices = self.dm.vertices_from_dedges(dedges)

			proposal_pathPatches[self.selected_proposal_id] = Polygon(vertices, closed=True, 
									edgecolor=self.boundary_colors[0], fill=False, linewidth=UNSELECTED_POLYGON_LINEWIDTH)
			proposal_vertexCircles[self.selected_proposal_id] = [plt.Circle(v, radius=UNSELECTED_CIRCLE_SIZE, color=self.boundary_colors[0], alpha=.8) for v in vertices]

		if self.shuffle_global_proposals:
			self.selected_proposal_polygon = self.global_proposal_pathPatches[self.selected_proposal_id]
			self.selected_proposal_vertexCircles = self.global_proposal_vertexCircles[self.selected_proposal_id]
		else:
			self.selected_proposal_polygon = self.local_proposal_pathPatches[self.selected_proposal_id]
			self.selected_proposal_vertexCircles = self.local_proposal_vertexCircles[self.selected_proposal_id]            

		if self.selected_proposal_polygon not in self.axis.patches:
			self.axis.add_patch(self.selected_proposal_polygon)

			for vertex_circ in self.selected_proposal_vertexCircles:
				self.axis.add_patch(vertex_circ)

		self.selected_proposal_polygon.set_picker(None)
		for vertex_circ in self.selected_proposal_vertexCircles:
			vertex_circ.set_picker(None)

		if self.selected_proposal_polygon in self.accepted_proposals:
			self.selected_proposal_polygon.set_linewidth(UNSELECTED_POLYGON_LINEWIDTH)
			label =  self.accepted_proposals[self.selected_proposal_polygon]['label']
		else:
			label = ''

		if self.shuffle_global_proposals:
			self.statusBar().showMessage('global proposal (%s) covering seed %d, score %.4f' % (label, sp_ind, self.global_proposal_sigs[self.selected_proposal_id]))
		else:
			self.statusBar().showMessage('local proposal (%s) from seed %d, score %.4f' % (label, sp_ind, sig))

		self.canvas.draw()


	def section_changed(self, item):

		self.statusBar().showMessage('Loading ....')

		if hasattr(self, 'global_proposal_tuples'):
			del self.global_proposal_tuples
		if hasattr(self, 'global_proposal_pathPatches'):
			for p in self.global_proposal_pathPatches:
				if p in self.axis.patches:
					p.remove()
			del self.global_proposal_pathPatches
		if hasattr(self, 'local_proposal_tuples'):
			del self.local_proposal_tuples
		if hasattr(self, 'local_proposal_pathPatches'):
			for p in self.local_proposal_pathPatches:
				if p in self.axis.patches:
					p.remove()
			del self.local_proposal_pathPatches

		sec = int(str(item.text()))
		self.init_data(section=sec)
		self.reload_brain_labeling_gui()

		self.mode_changed()
		# self.turn_superpixels_on()

		self.pixmap = QPixmap("/home/yuncong/CSHL_data_processed/%(stack)s_thumbnail_aligned_cropped/%(stack)s_%(sec)04d_thumbnail_aligned_cropped.tif"%{'sec':sec, 'stack':self.stack})
		self.pixmap_scaled = self.pixmap.scaledToHeight(self.bottom_panel.sizeHint().height())

		self.graphicsScene_navMap = QGraphicsScene(self.graphicsView_navMap)
		self.graphicsScene_navMap.addPixmap(self.pixmap_scaled)

		self.navRect = self.graphicsScene_navMap.addRect(10,10,200,200, QPen(QColor(255,0,0), 1))

		self.graphicsView_navMap.setScene(self.graphicsScene_navMap)
		self.graphicsView_navMap.show()

		self.navMap_scaling_x = self.pixmap_scaled.size().width()/float(self.dm.image_width)
		print self.navMap_scaling_x
		self.navMap_scaling_y = self.pixmap_scaled.size().height()/float(self.dm.image_height)

		self.update_navMap()

		self.graphicsScene_navMap.mousePressEvent = self.clicked_navMap

		self.statusBar().showMessage('Loading complete.')

	def clicked_navMap(self, event):
		# print event.scenePos()

		scene_x = event.scenePos().x()
		scene_y = event.scenePos().y()

		cur_data_left, cur_data_right = self.axis.get_xlim()
		cur_data_bottom, cur_data_top = self.axis.get_ylim()

		cur_data_center_x = .5 * cur_data_left + .5 * cur_data_right
		cur_data_center_y = .5 * cur_data_bottom + .5 * cur_data_top

		offset_x = scene_x / self.navMap_scaling_x - cur_data_center_x
		offset_y = scene_y / self.navMap_scaling_y - cur_data_center_y

		# print offset_x, offset_y

		new_data_left = cur_data_left + offset_x
		new_data_right = cur_data_right + offset_x
		new_data_bottom = cur_data_bottom + offset_y
		new_data_top = cur_data_top + offset_y

		# print new_data_left, new_data_right, new_data_bottom, new_data_top

		if new_data_right > self.dm.image_width:
			return
		if new_data_bottom > self.dm.image_height:
			return
		if new_data_left < 0:
			return
		if new_data_top < 0:
			return

		self.navRect.setRect(new_data_left * self.navMap_scaling_x, new_data_top * self.navMap_scaling_y, 
			self.navMap_scaling_x * (new_data_right - new_data_left), 
			self.navMap_scaling_y * (new_data_bottom - new_data_top))

		self.graphicsScene_navMap.update(0, 0, self.graphicsView_navMap.size().width(), self.graphicsView_navMap.size().height())
		self.graphicsView_navMap.setSceneRect(0, 0, self.dm.image_width*self.navMap_scaling_x, self.dm.image_height*self.navMap_scaling_y)

		self.axis.set_xlim([new_data_left, new_data_right])
		self.axis.set_ylim([new_data_bottom, new_data_top])

		self.canvas.draw()

	def save_callback(self):

		timestamp = datetime.datetime.now().strftime("%m%d%Y%H%M%S")

		if not hasattr(self, 'username') or self.username is None:
			username, okay = QInputDialog.getText(self, "Username", "Please enter your username:", QLineEdit.Normal, 'anon')
			if not okay: return
			self.username = str(username)

		accepted_proposal_props = []
		for patch, props in self.accepted_proposals.iteritems():

			props_saved = props.copy()

			props_saved['vertices'] = self.vertices_from_vertexPatches(patch)
			if patch.get_closed():
				props_saved['subtype'] = PolygonType.CLOSED
			else:
				props_saved['subtype'] = PolygonType.OPEN

			props_saved['labelPos'] = props['labelTextArtist'].get_position()

			props_saved.pop('vertexPatches')
			props_saved.pop('labelTextArtist')

			accepted_proposal_props.append(props_saved)

		self.dm.save_proposal_review_result(accepted_proposal_props, self.username, timestamp, suffix='consolidated')

		print self.new_labelnames
		self.dm.add_labelnames(self.new_labelnames, self.dm.repo_dir+'/visualization/newStructureNames.txt')

		self.statusBar().showMessage('Labelings saved to %s' % (self.username+'_'+timestamp))

		cur_xlim = self.axis.get_xlim()
		cur_ylim = self.axis.get_ylim()

		# self.axis.set_xlim([0, self.dm.image_width])
		# self.axis.set_ylim([self.dm.image_height, 0])

		self.fig.savefig('/tmp/preview.jpg', bbox_inches='tight')

		# self.axis.set_xlim(cur_xlim)
		# self.axis.set_ylim(cur_ylim)

		self.canvas.draw()

	def labelbutton_callback(self):
		pass

	############################################
	# matplotlib canvas CALLBACKs
	############################################


	def undo(self):

		if len(self.history) == 0:
			return

		history_item = self.history.pop()

		if history_item['type'] == 'drag_polygon':

			polygon = history_item['polygon']
			moved_x, moved_y = history_item['mouse_moved']

			for c in self.accepted_proposals[polygon]['vertexPatches']:
				c.center = (c.center[0] - moved_x, c.center[1] - moved_y)

			polygon.set_xy(self.vertices_from_vertexPatches(polygon))

			if 'labelTextArtist' in self.accepted_proposals[polygon]:
				curr_label_pos = self.accepted_proposals[polygon]['labelTextArtist'].get_position()
				self.accepted_proposals[polygon]['labelTextArtist'].set_position((curr_label_pos[0] - moved_x, 
									curr_label_pos[1] - moved_y))

			self.canvas.draw()

		elif history_item['type'] == 'drag_vertex':
			polygon = history_item['polygon']
			index = history_item['index']
			circ = history_item['circle']

			circ.center = (circ.center[0] - history_item['mouse_moved'][0], circ.center[1] - history_item['mouse_moved'][1])

			xys = self.vertices_from_vertexPatches(polygon)
			polygon.set_xy(xys)

			# if 'labelTextArtist' in self.accepted_proposals[polygon]:
			#     self.accepted_proposals[polygon]['labelTextArtist'].set_position(np.mean(xys, axis=0))

			self.canvas.draw()

		elif history_item['type'] == 'add_vertex':
			polygon = history_item['polygon']
			index = history_item['index']
			circ = history_item['circle']

			circ.remove()

			vertices = np.array(self.vertices_from_vertexPatches(polygon))
			new_vertices = np.delete(vertices, index, axis=0)
			polygon.set_xy(new_vertices)

			self.accepted_proposals[polygon]['vertexPatches'].remove(circ)

			if index == 0:
				polygon.remove()
				self.accepted_proposals.pop(polygon)
				# self.set_mode(Mode.REVIEW_PROPOSAL)

			self.canvas.draw()


	def set_mode(self, mode):
		if hasattr(self, 'mode'):
			if self.mode != mode:
				print self.mode, '=>', mode
		else:
			print mode
		self.mode = mode


	def on_key_release(self, event):
		pass

		# if event.key == 'control':
		#     if self.mode == Mode.ADDING_VERTICES_RANDOMLY:
		#         self.statusBar().showMessage('Done adding vertices.')
		#         self.set_mode(Mode.REVIEW_PROPOSAL)

		# if event.key == 'shift':
		#     if self.mode == Mode.ADDING_VERTICES_CONSECUTIVELY:
		#         self.statusBar().showMessage('Done adding vertices.')
		#         self.set_mode(Mode.REVIEW_PROPOSAL)

	def on_key_press(self, event):

		if event.key == ' ':
			print 'esc pressed'
			self.cancel_current_selection()

		if event.key == '[' or event.key == ']':
			if event.key == '[':
				if not hasattr(self, 'section2'):
					self.section2 = self.section - 1
				else:
					if self.section2 == self.section - 3:
						return
					else:
						self.section2 = self.section2 - 1
			else:
				if not hasattr(self, 'section2'):
					self.section2 = self.section + 1
				else:
					if self.section2 == self.section + 3:
						return
					else:
						self.section2 = self.section2 + 1

			self.section_text2.set_text(str(self.section2))

			print self.section, self.section2
			
			self.dm2 = self.dms[self.section2]

			t = time.time()
			self.dm2._load_image(versions=['rgb-jpg'])
			print 'load image', time.time() - t

			t = time.time()
			self.orig_image_handle2.set_data(self.dm2.image_rgb_jpg)
			print 'set data', time.time() - t

			self.setWindowTitle('BrainLabelingGUI, stack %s'%self.stack + ', section %d' %self.section + ', section %d'%self.section2)

		# if event.key == 'control':

		#     self.set_mode(Mode.ADDING_VERTICES_RANDOMLY)

		#     if self.selected_proposal_polygon is None:
		#         self.statusBar().showMessage('you must select a marking first')
		#     else:
		#         self.statusBar().showMessage('Random Mode: Left click to add vertices ...')
			

		# if event.key == 'shift':

		#     self.set_mode(Mode.ADDING_VERTICES_CONSECUTIVELY)

		#     if self.selected_circle is None:
		#         self.statusBar().showMessage('you must select a vertex first')
		#     else:
		#         self.statusBar().showMessage('Consecutive Mode: Left click to add vertices ...')
						

		if event.key == 'ctrl+z':
			self.undo()

		if event.key == '=' or event.key == '-':
			self.on_zoom(event)

		self.cur_xlim = self.axis.get_xlim()
		cur_xmin, cur_xmax = self.cur_xlim
		
		self.cur_ylim = self.axis.get_ylim()
		cur_ybottom, cur_ytop = self.cur_ylim

		if event.key == 'left' or event.key == 'right' or event.key == 'up' or event.key == 'down':

			if event.key == 'left':
				offset_x, offset_y = self.find_proper_offset(100, 0)
			elif event.key == 'right':
				offset_x, offset_y = self.find_proper_offset(-100, 0)
			elif event.key == 'down':
				offset_x, offset_y = self.find_proper_offset(0, -100)
			elif event.key == 'up':
				offset_x, offset_y = self.find_proper_offset(0, 100)

			print offset_x, offset_y


			t = time.time()
			self.axis.set_xlim([cur_xmin - offset_x, cur_xmax - offset_x])
			self.axis.set_ylim([cur_ybottom - offset_y, cur_ytop - offset_y])

			self.axis2.set_xlim([cur_xmin - offset_x, cur_xmax - offset_x])
			self.axis2.set_ylim([cur_ybottom - offset_y, cur_ytop - offset_y])
			sys.stderr.write('set lim in %f\n' % (time.time()-t))

		t = time.time()
		# if int(t) % 100 < 20:
		self.canvas.draw()
		self.canvas2.draw()
		sys.stderr.write('draw in %f\n' % (time.time()-t))

		t = time.time()
		self.update_navMap()
		sys.stderr.write('update %f\n' % (time.time()-t))

	def on_zoom(self, event):

		data_left, data_right = self.axis.get_xlim()
		data_bottom, data_top = self.axis.get_ylim()

		mouse_distance_to_left_data = event.xdata - data_left
		mouse_distance_to_right_data = data_right - event.xdata
		mouse_distance_to_bottom_data = data_bottom - event.ydata
		mouse_distance_to_top_data = event.ydata - data_top

		if event.xdata is None or event.ydata is None: # mouse position outside data region
			return

		if isinstance(event, MouseEvent):

			if event.button == 'up':
				# deal with zoom in
				scale_factor = 1/self.base_scale
			elif event.button == 'down':
				# deal with zoom out
				scale_factor = self.base_scale

		elif isinstance(event, KeyEvent):
			if event.key == '=':
				scale_factor = 1/self.base_scale
			elif event.key == '-':
				scale_factor = self.base_scale

		if event.xdata - mouse_distance_to_left_data * scale_factor < 0 or \
			event.xdata + mouse_distance_to_right_data * scale_factor > self.dm.image_width or \
			event.ydata - mouse_distance_to_top_data * scale_factor < 0 or \
			event.ydata + mouse_distance_to_bottom_data * scale_factor > self.dm.image_height:
			return

		new_left_data = max(0, event.xdata - mouse_distance_to_left_data * scale_factor)
		new_right_data = min(self.dm.image_width, event.xdata + mouse_distance_to_right_data * scale_factor)
		new_top_data = max(0, event.ydata - mouse_distance_to_top_data * scale_factor)
		new_bottom_data = min(self.dm.image_height, event.ydata + mouse_distance_to_bottom_data * scale_factor)

		# print new_left_data, new_right_data, new_bottom_data, new_top_data

		self.axis.set_xlim([new_left_data, new_right_data])
		self.axis.set_ylim([new_bottom_data, new_top_data])

		self.axis2.set_xlim([new_left_data, new_right_data])
		self.axis2.set_ylim([new_bottom_data, new_top_data])

		self.canvas.draw() # force re-draw
		self.canvas2.draw() # force re-draw

		self.update_navMap()
		
	def update_navMap(self):

		cur_xmin, cur_xmax = self.axis.get_xlim()
		cur_ybottom, cur_ytop = self.axis.get_ylim()
		self.navRect.setRect(cur_xmin * self.navMap_scaling_x, cur_ybottom * self.navMap_scaling_y, self.navMap_scaling_x * (cur_xmax - cur_xmin), self.navMap_scaling_y * (cur_ytop - cur_ybottom))
		self.graphicsScene_navMap.update(0, 0, self.graphicsView_navMap.size().width(), self.graphicsView_navMap.size().height())
		self.graphicsView_navMap.setSceneRect(0, 0, self.dm.image_width*self.navMap_scaling_x, self.dm.image_height*self.navMap_scaling_y)


	def on_press(self, event):

		if len(self.picked_artists) > 0:

			print self.picked_artists

			picked_artist = None

			for artist in self.picked_artists:

				for patch, props in self.accepted_proposals.iteritems():
					if 'labelTextArtist' in props and artist == props['labelTextArtist']:
						print 'clicked on a label'
						picked_artist = artist
						self.picked_object_type = 'label'
						break

				if picked_artist is not None: 
					break

				for patch, props in self.accepted_proposals.iteritems():
					if 'vertexPatches' in props and artist in props['vertexPatches']:
						print 'clicked on a vertex'
						picked_artist = artist
						self.picked_object_type = 'vertex'
						break

			if picked_artist is None: 

				matched_polygons = [artist for artist in self.picked_artists if artist in self.accepted_proposals]

				picked_artist = matched_polygons[np.argmin([(ShapelyLineString(self.vertices_from_vertexPatches(p) \
									if not p.get_closed() else ShapelyLineRing(self.vertices_from_vertexPatches(p)))).distance(\
										ShapelyPoint(event.xdata, event.ydata))
					for p in matched_polygons])]

				self.picked_object_type = 'polygon'

			print picked_artist

			if self.mode == Mode.ADDING_VERTICES_CONSECUTIVELY:
				if self.picked_object_type == 'vertex':
					self.previous_circle = self.selected_circle
					if hasattr(self, 'selected_vertex_index'):
						self.previous_vertex_index = self.selected_vertex_index
					else:
						self.previous_vertex_index = None

					self.previous_polygon = self.selected_proposal_polygon

					print 'prev is set to', self.previous_polygon, self.previous_vertex_index

					self.connecting_vertices = True
					print 'self.connecting_vertices set to', self.connecting_vertices

			self.handle_picked_object(picked_artist)


		self.picked_artists = []

		self.press_x = event.xdata
		self.press_y = event.ydata

		self.press_x_canvas = event.x
		self.press_y_canvas = event.y

		print 'press'
		self.pressed = True
		self.press_time = time.time()

		self.cur_xlim = self.axis.get_xlim()
		self.cur_ylim = self.axis.get_ylim()

		if event.button == 1:

			if self.mode == Mode.SELECTING_ROI:
				print event.xdata, event.ydata
				self.roi_xmin = event.xdata
				self.roi_ymin = event.ydata

			else:
				if hasattr(self, 'roi_rectPatch') and self.roi_rectPatch is not None and self.roi_rectPatch in self.axis.patches:
					self.roi_rectPatch.remove()
					self.roi_rectPatch = None
					self.roi_xmin = None
					self.roi_xmax = None
					self.roi_ymin = None
					self.roi_ymax = None
					self.canvas.draw()

	def on_motion(self, event):

		if not hasattr(self, 'mode'):
			return

		if self.mode == Mode.SELECTING_ROI and hasattr(self, 'roi_xmin') and self.roi_xmin is not None:

			self.roi_xmax = event.xdata
			self.roi_ymax = event.ydata

			if hasattr(self, 'roi_rectPatch') and self.roi_rectPatch is not None and self.roi_rectPatch in self.axis.patches:

				self.roi_rectPatch.set_width(self.roi_xmax-self.roi_xmin)
				self.roi_rectPatch.set_height(self.roi_ymax-self.roi_ymin)
			else:
				self.roi_rectPatch = Rectangle((self.roi_xmin, self.roi_ymin), self.roi_xmax-self.roi_xmin, self.roi_ymax-self.roi_ymin, edgecolor=(1,0,0), fill=False, linestyle='dashed')
				self.axis.add_patch(self.roi_rectPatch)

			self.canvas.draw()

			return

		if self.mode != Mode.ADDING_VERTICES_CONSECUTIVELY and self.mode != Mode.ADDING_VERTICES_RANDOMLY \
					and hasattr(self, 'selected_circle') and self.selected_circle is not None and self.pressed: # drag vertex
		# if self.mode == Mode.MOVING_VERTEX \


			self.set_mode(Mode.MOVING_VERTEX)

			print 'dragging vertex'

			self.selected_circle.center = event.xdata, event.ydata

			xys = self.selected_proposal_polygon.get_xy()
			xys[self.selected_vertex_index] = self.selected_circle.center

			if self.selected_proposal_polygon.get_closed():
				self.selected_proposal_polygon.set_xy(xys[:-1])
			else:
				self.selected_proposal_polygon.set_xy(xys)

			# self.accepted_proposals[self.selected_proposal_polygon]['vertices'] = xys

			# if 'labelTextArtist' in self.accepted_proposals[self.selected_proposal_polygon]:
			#     self.accepted_proposals[self.selected_proposal_polygon]['labelTextArtist'].set_position(xys.mean(axis=0))
			
			self.canvas.draw()

		elif self.mode == Mode.MOVING_POLYGON and hasattr(self, 'selected_proposal_polygon') and \
			self.pressed and self.selected_proposal_polygon in self.accepted_proposals:
			# drag polygon

			print 'dragging polygon'

			offset_x = event.xdata - self.press_x
			offset_y = event.ydata - self.press_y

			for c, center0 in zip(self.accepted_proposals[self.selected_proposal_polygon]['vertexPatches'], self.selected_polygon_circle_centers_before_drag):
				c.center = (center0[0] + offset_x, center0[1] + offset_y)

			xys = self.vertices_from_vertexPatches(self.selected_proposal_polygon)
			self.selected_proposal_polygon.set_xy(xys)

			if 'labelTextArtist' in self.accepted_proposals[self.selected_proposal_polygon]:
				self.accepted_proposals[self.selected_proposal_polygon]['labelTextArtist'].set_position((self.selected_polygon_label_pos_before_drag[0]+offset_x, 
																										self.selected_polygon_label_pos_before_drag[1]+offset_y))

			self.canvas.draw()

		elif hasattr(self, 'picked_object_type') and self.picked_object_type == 'label':

			print 'dragging label'
			self.accepted_proposals[self.selected_proposal_polygon]['labelTextArtist'].set_position((event.xdata, event.ydata))
			self.canvas.draw()


		elif hasattr(self, 'pressed') and self.pressed and (event.x - self.press_x_canvas) ** 2 + (event.y - self.press_y_canvas) ** 2 > PAN_THRESHOLD ** 2:
			# canvas move distance constraint is to prevent accidental panning caused by draggin stylus while clicking
			# and time.time() - self.press_time > .1:

			print 'panning canvas'
			
			if (event.xdata is None) | (event.ydata is None):
				return

			offset_x, offset_y = self.axis.transData.inverted().transform((event.x,event.y)) - self.axis.transData.inverted().transform((self.press_x_canvas,self.press_y_canvas))
			# offset_x > 0 move right
			# offset_y > 0 move up
			
			offset_x, offset_y = self.find_proper_offset(offset_x, offset_y)

			t = time.time()

			self.axis.set_ylim(bottom=self.cur_ylim[0] - offset_y, top=self.cur_ylim[1] - offset_y)
			self.axis.set_xlim(left=self.cur_xlim[0] - offset_x, right=self.cur_xlim[1] - offset_x)

			self.axis2.set_ylim(bottom=self.cur_ylim[0] - offset_y, top=self.cur_ylim[1] - offset_y)
			self.axis2.set_xlim(left=self.cur_xlim[0] - offset_x, right=self.cur_xlim[1] - offset_x)

			self.canvas.draw()
			self.canvas2.draw()

			sys.stderr.write('draw in %f\n' % (time.time()-t))

			self.update_navMap()

	def find_proper_offset(self, offset_x, offset_y):

		if self.cur_ylim[0] - offset_y > self.dm.image_height:
			offset_y = self.cur_ylim[0] - self.dm.image_height
		elif self.cur_ylim[1] - offset_y < 0:
			offset_y = self.cur_ylim[1]

		if self.cur_xlim[1] - offset_x > self.dm.image_width:
			offset_x = self.dm.image_width - self.cur_xlim[1]
		elif self.cur_xlim[0] - offset_x < 0:
			offset_x = self.cur_xlim[0]

		return offset_x, offset_y

	def on_release(self, event):
		print 'release'

		self.pressed = False

		self.release_x = event.xdata
		self.release_y = event.ydata

		self.release_x_canvas = event.x
		self.release_y_canvas = event.y

		if self.mode == Mode.MOVING_VERTEX:
			self.set_mode(Mode.REVIEW_PROPOSAL)
			self.history.append({'type': 'drag_vertex', 'polygon': self.selected_proposal_polygon, 'index': self.selected_vertex_index, 
								'circle': self.selected_circle, 'mouse_moved': (self.release_x - self.press_x, self.release_y - self.press_y)})

		elif self.mode == Mode.MOVING_POLYGON:
			self.set_mode(Mode.REVIEW_PROPOSAL) 
			self.history.append({'type': 'drag_polygon', 'polygon': self.selected_proposal_polygon, 'mouse_moved': (self.release_x - self.press_x, self.release_y - self.press_y)})

		self.release_time = time.time()

		if event.button == 1:

			if self.mode == Mode.SELECTING_ROI:
				self.set_mode(Mode.REVIEW_PROPOSAL)

		if abs(self.release_x_canvas - self.press_x_canvas) < 10 and abs(self.release_y_canvas - self.press_y_canvas) < 10:
			# short movement

			print 'short movement'

			if event.button == 1: # left click
							
				if self.mode == Mode.CREATING_NEW_POLYGON:
					self.add_vertex([event.xdata, event.ydata], create_if_no_selected=True, 
										consecutive=True)

				elif self.mode == Mode.ADDING_VERTICES_CONSECUTIVELY:

					# print 'prev', self.previous_polygon, self.previous_vertex_index, 'curr', self.selected_proposal_polygon, self.selected_vertex_index

					print 'self.connecting_vertices', self.connecting_vertices
					if self.connecting_vertices:

						if self.previous_polygon == self.selected_proposal_polygon:
							self.connect_two_vertices(self.selected_proposal_polygon)
							self.cancel_current_selection()
							self.set_mode(Mode.REVIEW_PROPOSAL)

						else:
							self.connect_two_vertices(self.selected_proposal_polygon, self.previous_polygon, 
														self.selected_vertex_index, self.previous_vertex_index)
							self.cancel_current_selection()
							self.set_mode(Mode.REVIEW_PROPOSAL)

						self.connecting_vertices = False
					else:
						self.add_vertex([event.xdata, event.ydata], consecutive=True)
				
				elif self.mode == Mode.ADDING_VERTICES_RANDOMLY:
					self.add_vertex([event.xdata, event.ydata], consecutive=False)

				else:

					self.cancel_current_selection()

					# if self.superpixels_on:
					print 'clicked a superpixel'
		
					if self.superpixels_on:
						self.clicked_sp = self.dm.segmentation[int(event.ydata), int(event.xdata)]
						sys.stderr.write('clicked sp %d\n'%self.clicked_sp)

						if self.mode == Mode.REVIEW_PROPOSAL:
							self.shuffle_proposal_from_pool(self.clicked_sp)

			elif event.button == 3: # right click
				canvas_pos = (event.xdata, event.ydata)
				self.openMenu(canvas_pos)

		else:
			print 'long movement'
			# long movement

			if self.mode != Mode.ADDING_VERTICES_CONSECUTIVELY and self.mode != Mode.ADDING_VERTICES_RANDOMLY:
				self.cancel_current_selection()        

		self.canvas.draw()

		print 'release finished'
		self.object_picked = False
		self.picked_object_type = None

		print self.selected_proposal_polygon, self.selected_circle, self.selected_vertex_index if hasattr(self, 'selected_vertex_index') else ''


	def find_vertex_insert_position(self, xys, pos, closed=True):

		n = len(xys)
		if n == 1:
			return 1

		xys_homo = np.column_stack([xys, np.ones(n,)])

		if closed:
			edges = np.array([np.cross(xys_homo[i], xys_homo[(i+1)%n]) for i in range(n)])
		else:
			edges = np.array([np.cross(xys_homo[i], xys_homo[i+1]) for i in range(n-1)])

		edges_normalized = edges/np.sqrt(np.sum(edges[:,:2]**2, axis=1))[:, np.newaxis]

		signed_dists = np.dot(edges_normalized, np.r_[pos,1])
		dists = np.abs(signed_dists)
		# sides = np.sign(signed_dists)

		projections = pos - signed_dists[:, np.newaxis] * edges_normalized[:,:2]

		endpoint = [None for _ in projections]
		for i, (px, py) in enumerate(projections):
			if (px > xys[i][0] and px > xys[(i+1)%n][0]) or (px < xys[i][0] and px < xys[(i+1)%n][0]):
				endpoint[i] = [i, (i+1)%n][np.argmin(np.squeeze(cdist([pos], [xys[i], xys[(i+1)%n]])))]
				dists[i] = np.min(np.squeeze(cdist([pos], [xys[i], xys[(i+1)%n]])))

		# print edges_normalized[:,:2]
		# print projections                
		# print dists
		# print endpoint
		nearest_edge_begins_at = np.argsort(dists)[0]

		if nearest_edge_begins_at == 0 and not closed and endpoint[0] == 0:
			new_vertex_ind = 0
		elif nearest_edge_begins_at == n-2 and not closed and endpoint[-1] == n-1:
			new_vertex_ind = n
		else:
			new_vertex_ind = nearest_edge_begins_at + 1  

		print 'nearest_edge_begins_at', nearest_edge_begins_at, 'new_vertex_ind', new_vertex_ind

		return new_vertex_ind


	def add_vertex(self, pos, create_if_no_selected=False, consecutive=False):
		print 'add vertex'

		# self.cancel_current_circle()
		self.auto_extend_view(pos[0], pos[1])

		if (not hasattr(self, 'selected_proposal_polygon') or self.selected_proposal_polygon is None) and not create_if_no_selected:
			print 'no proposal selected'
			return

		vertex_circle = plt.Circle(pos, radius=UNSELECTED_CIRCLE_SIZE, color=self.boundary_colors[1], alpha=.8)
		self.axis.add_patch(vertex_circle)
		vertex_circle.set_picker(CIRCLE_PICK_THRESH)

		if not hasattr(self, 'selected_proposal_polygon') or self.selected_proposal_polygon is None:
			
			self.selected_proposal_polygon = Polygon([pos], closed=False, edgecolor=self.boundary_colors[1], fill=False, linewidth=UNSELECTED_POLYGON_LINEWIDTH)

			self.selected_proposal_polygon.set_picker(True)
			self.axis.add_patch(self.selected_proposal_polygon)

			self.selected_proposal_type = ProposalType.FREEFORM
			self.accepted_proposals[self.selected_proposal_polygon] = {'type': ProposalType.FREEFORM,
																	# 'subtype': PolygonType.OPEN,
																	'vertexPatches': []
																	}

			new_vertex_ind = 0

			self.set_mode(Mode.ADDING_VERTICES_CONSECUTIVELY)
			# self.set_mode(Mode.REVIEW_PROPOSAL)

			self.accepted_proposals[self.selected_proposal_polygon]['vertexPatches'] = [vertex_circle]

		else:

			xys = self.vertices_from_vertexPatches(self.selected_proposal_polygon)

			n = len(xys)

			if consecutive:
				if hasattr(self, 'selected_vertex_index') and self.extend_head: # if a circle is selected, add after this vertex
					xys = np.vstack([xys[::-1], pos])
					self.selected_proposal_polygon.set_xy(xys)
					self.accepted_proposals[self.selected_proposal_polygon]['vertexPatches'] = self.accepted_proposals[self.selected_proposal_polygon]['vertexPatches'][::-1] + [vertex_circle]
				else:

					xys = np.vstack([xys, pos])
					self.selected_proposal_polygon.set_xy(xys)
					self.accepted_proposals[self.selected_proposal_polygon]['vertexPatches'].append(vertex_circle)
				
				new_vertex_ind = n

			else:

				new_vertex_ind = self.find_vertex_insert_position(xys, pos, self.selected_proposal_polygon.get_closed())
				xys = np.insert(xys, new_vertex_ind, pos, axis=0)
				self.selected_proposal_polygon.set_xy(xys)
				self.accepted_proposals[self.selected_proposal_polygon]['vertexPatches'].insert(new_vertex_ind, vertex_circle)

		self.history.append({'type': 'add_vertex', 'polygon': self.selected_proposal_polygon, 'index': new_vertex_ind, 
					'circle': vertex_circle})

		self.cancel_current_circle()

		self.extend_head = False

		self.canvas.draw()


	def connect_two_vertices(self, polygon1, polygon2=None, index1=None, index2=None):

		vertices1 = self.vertices_from_vertexPatches(polygon1)
		n1 = len(vertices1)

		if polygon2 is None: # connect two ends of a single polygon   
			
			print 'index1', index1, 'index2', index2
			assert not polygon1.get_closed()
			 # and ((index1 == 0 and index2 == n2-1) or (index1 == n1-1 and index2 == 0))
			polygon1.set_closed(True)
			if 'label' not in self.accepted_proposals[polygon1]:
				self.acceptProposal_callback()
		
		else:

			vertices2 = self.vertices_from_vertexPatches(polygon2)
			n2 = len(vertices2)
			print 'index1', index1, 'index2', index2
			assert not polygon1.get_closed() and index1 in [0, n1-1] and index2 in [0, n2-1]

			if index1 == 0 and index2 == 0:
				new_vertices = np.vstack([vertices1[::-1], vertices2])
			elif index1 != 0 and index2 != 0:
				new_vertices = np.vstack([vertices1, vertices2[::-1]])
			elif index1 != 0 and index2 == 0:
				new_vertices = np.vstack([vertices1, vertices2])
			elif index1 == 0 and index2 != 0:
				new_vertices = np.vstack([vertices1[::-1], vertices2[::-1]])

			props = {}

			patch = Polygon(new_vertices, closed=False, edgecolor=self.boundary_colors[1], fill=False, linewidth=UNSELECTED_POLYGON_LINEWIDTH)
			patch.set_picker(True)

			self.axis.add_patch(patch)

			props['vertexPatches'] = []
			for x,y in new_vertices:
				vertex_circle = plt.Circle((x, y), radius=UNSELECTED_CIRCLE_SIZE, color=self.boundary_colors[1], alpha=.8)
				vertex_circle.set_picker(CIRCLE_PICK_THRESH)
				props['vertexPatches'].append(vertex_circle)
				self.axis.add_patch(vertex_circle)
				vertex_circle.set_picker(True)

			if self.accepted_proposals[polygon1]['label'] == self.accepted_proposals[polygon2]['label']:
				props['label'] = self.accepted_proposals[polygon1]['label']
			else:
				props['label'] = self.accepted_proposals[polygon1]['label']
			# else:
				# self.acceptProposal_callback()

			props['type'] = self.accepted_proposals[polygon1]['type']

			centroid = np.mean(new_vertices, axis=0)
			props['labelTextArtist'] = Text(centroid[0], centroid[1], props['label'], style='italic', bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})

			self.axis.add_artist(props['labelTextArtist'])
			props['labelTextArtist'].set_picker(True)

			self.accepted_proposals[patch] = props
			
			for circ in self.accepted_proposals[polygon1]['vertexPatches']:
				circ.remove()
			
			for circ in self.accepted_proposals[polygon2]['vertexPatches']:
				circ.remove()

			self.accepted_proposals[polygon1]['labelTextArtist'].remove()
			self.accepted_proposals[polygon2]['labelTextArtist'].remove()
			
			polygon1.remove()
			polygon2.remove()

			self.accepted_proposals.pop(polygon1)
			self.accepted_proposals.pop(polygon2)

		self.cancel_current_selection()

		self.canvas.draw()


	def vertices_from_vertexPatches(self, polygon):
		if polygon not in self.accepted_proposals:
			assert 'polygon does not exist'
			return None
		else:
			xys = map(attrgetter('center'), self.accepted_proposals[polygon]['vertexPatches'])
			return xys


	def set_polygon_vertices(self, polygon, vertices):
		self.selected_proposal_polygon.set_xy(new_vertices)
		self.selected_proposal_polygon['vertexPatches'] = self.selected_proposal_polygon['vertexPatches'][v2:] + self.selected_proposal_polygon['vertexPatches'][:v2]


	def break_edge(self, pos):

		print 'break edge'

		xys = self.vertices_from_vertexPatches(self.selected_proposal_polygon)
		v2 = self.find_vertex_insert_position(xys, pos, self.selected_proposal_polygon.get_closed())

		print v2

		n = len(self.accepted_proposals[self.selected_proposal_polygon]['vertexPatches'])

		vertices = self.vertices_from_vertexPatches(self.selected_proposal_polygon)
		print 'vertices', vertices

		if self.selected_proposal_polygon.get_closed():

			self.selected_proposal_polygon.set_closed(False)

			if v2 == 0 or v2 == n:
				new_vertices = vertices
			else:
				new_vertices = np.vstack([vertices[v2:], vertices[:v2]])

			self.selected_proposal_polygon.set_xy(new_vertices)
			self.accepted_proposals[self.selected_proposal_polygon]['vertexPatches'] = self.accepted_proposals[self.selected_proposal_polygon]['vertexPatches'][v2:] + \
																self.accepted_proposals[self.selected_proposal_polygon]['vertexPatches'][:v2]

			self.accepted_proposals[self.selected_proposal_polygon]['labelTextArtist'].set_position(np.mean(new_vertices, axis=0))
			
		else:

			if v2 == 1 or v2 == n-1:

				if v2 == n-1:
					new_vertices = vertices[:-1]
					circ = self.accepted_proposals[self.selected_proposal_polygon]['vertexPatches'][v2]
				else:
					new_vertices = vertices[1:]
					circ = self.accepted_proposals[self.selected_proposal_polygon]['vertexPatches'][v2]

				self.selected_proposal_polygon.set_xy(new_vertices)
				
				# self.accepted_proposals[self.selected_proposal_polygon]['vertices'] = new_vertices
				
				circ.remove()
				self.accepted_proposals[self.selected_proposal_polygon]['vertexPatches'].remove(circ)

				centroid = np.mean(new_vertices, axis=0)
				self.accepted_proposals[self.selected_proposal_polygon]['labelTextArtist'].set_position(centroid)

				if self.selected_proposal_polygon.get_closed():
					self.selected_proposal_polygon.set_closed(False)

			else:
				new_vertices1 = vertices[:v2]
				new_vertices2 = vertices[v2:]

				for new_vertices in [new_vertices1, new_vertices2]:
					props = {}
					# props['vertices'] = new_vertices

					print v2
					print new_vertices

					patch = Polygon(new_vertices, closed=False, edgecolor=self.boundary_colors[1], 
							fill=False, linewidth=UNSELECTED_POLYGON_LINEWIDTH)

					patch.set_picker(True)
					self.axis.add_patch(patch)

					props['vertexPatches'] = []
					for x,y in new_vertices:
						vertex_circle = plt.Circle((x, y), radius=UNSELECTED_CIRCLE_SIZE, color=self.boundary_colors[1], alpha=.8)
						vertex_circle.set_picker(CIRCLE_PICK_THRESH)
						props['vertexPatches'].append(vertex_circle)
						self.axis.add_patch(vertex_circle)
						vertex_circle.set_picker(True)

					props['label'] = self.accepted_proposals[self.selected_proposal_polygon]['label']
					props['type'] = self.accepted_proposals[self.selected_proposal_polygon]['type']

					# props['subtype'] = PolygonType.OPEN

					centroid = np.mean(new_vertices, axis=0)
					props['labelTextArtist'] = Text(centroid[0], centroid[1], props['label'], style='italic', bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})

					self.axis.add_artist(props['labelTextArtist'])
					props['labelTextArtist'].set_picker(True)

					self.accepted_proposals[patch] = props

				for circ in self.accepted_proposals[self.selected_proposal_polygon]['vertexPatches']:
					circ.remove()

				self.accepted_proposals[self.selected_proposal_polygon]['labelTextArtist'].remove()

				self.accepted_proposals.pop(self.selected_proposal_polygon)
				self.selected_proposal_polygon.remove()
				self.selected_proposal_polygon = None


		self.canvas.draw()

	def remove_selected_vertex(self):

		# print self.selected_circle

		self.selected_circle.remove()
		self.accepted_proposals[self.selected_proposal_polygon]['vertexPatches'].remove(self.selected_circle)
		self.selected_circle = None

		p = self.selected_proposal_polygon

		xys = p.get_xy()
		xys = np.vstack([xys[:self.selected_vertex_index], xys[self.selected_vertex_index+1:]])

		vertices = xys[:-1] if p.get_closed() else xys

		self.selected_proposal_polygon.set_xy(vertices)

		self.accepted_proposals[self.selected_proposal_polygon]['vertices'] = vertices

		self.canvas.draw()


	def remove_selected_vertices_in_region(self, link_endpoints=True):

		done = False
		for patch, props in self.accepted_proposals.iteritems():
			if done: break
			for x,y in self.vertices_from_vertexPatches(patch):
				if x >= self.roi_xmin and x <= self.roi_xmax and y >= self.roi_ymin and y <= self.roi_ymax:
					self.selected_proposal_polygon = patch
					done = True
					break

		print self.selected_proposal_polygon

		n = len(self.accepted_proposals[self.selected_proposal_polygon]['vertexPatches'])

		in_roi_vertex_indices = []
		for vertex_ind, circ in enumerate(self.accepted_proposals[self.selected_proposal_polygon]['vertexPatches']):
			x, y = circ.center
			if x >= self.roi_xmin and x <= self.roi_xmax and y >= self.roi_ymin and y <= self.roi_ymax:
				# circ.remove()
				# self.accepted_proposals[self.selected_proposal_polygon]['vertexPatches'].remove(circ)
				in_roi_vertex_indices.append(vertex_ind)

		print in_roi_vertex_indices

		start_gappoints = [i for i in (np.array(in_roi_vertex_indices) + 1)%n if i not in in_roi_vertex_indices]
		end_gappoints = [i for i in (np.array(in_roi_vertex_indices) - 1)%n if i not in in_roi_vertex_indices]

		closed = self.selected_proposal_polygon.get_closed()

		if not closed:
			if 0 not in start_gappoints and n-1 not in end_gappoints:
				start_gappoints += [0]
				end_gappoints += [n-1]

		all_gappoints = np.sort(start_gappoints + end_gappoints)
		print start_gappoints, end_gappoints

		if all_gappoints[0] in start_gappoints:
			remaining_segments = zip(all_gappoints[::2], all_gappoints[1::2])
		else:
			remaining_segments = [(all_gappoints[-1], all_gappoints[0])] + zip(all_gappoints[1::2], all_gappoints[2::2])

		vertices = self.selected_proposal_polygon.get_xy()
		
		if closed:
			vertices = vertices[:-1]

		for start_gap, end_gap in remaining_segments:

			print start_gap, end_gap

			if end_gap < start_gap:
				new_vertices = np.vstack([vertices[start_gap:], vertices[:end_gap+1]])
			else:
				new_vertices = vertices[start_gap:end_gap+1]

			props = {}
			# props['vertices'] = new_vertices

			patch = Polygon(new_vertices, closed=link_endpoints, edgecolor=self.boundary_colors[1], fill=False, linewidth=UNSELECTED_POLYGON_LINEWIDTH)

			patch.set_picker(True)
			self.axis.add_patch(patch)

			props['vertexPatches'] = []
			for x,y in new_vertices:
				vertex_circle = plt.Circle((x, y), radius=UNSELECTED_CIRCLE_SIZE, color=self.boundary_colors[1], alpha=.8)
				vertex_circle.set_picker(CIRCLE_PICK_THRESH)
				props['vertexPatches'].append(vertex_circle)
				self.axis.add_patch(vertex_circle)
				vertex_circle.set_picker(True)

			props['label'] = self.accepted_proposals[self.selected_proposal_polygon]['label']
			props['type'] = self.accepted_proposals[self.selected_proposal_polygon]['type']

			# if link_endpoints:
			#     props['subtype'] = PolygonType.CLOSED
			# else:
			#     props['subtype'] = PolygonType.OPEN

			centroid = np.mean(new_vertices, axis=0)
			props['labelTextArtist'] = Text(centroid[0], centroid[1], props['label'], style='italic', bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})

			self.axis.add_artist(props['labelTextArtist'])
			props['labelTextArtist'].set_picker(True)

			self.accepted_proposals[patch] = props

		for circ in self.accepted_proposals[self.selected_proposal_polygon]['vertexPatches']:
			circ.remove()

		self.accepted_proposals[self.selected_proposal_polygon]['labelTextArtist'].remove()

		self.accepted_proposals.pop(self.selected_proposal_polygon)
		self.selected_proposal_polygon.remove()
		self.selected_proposal_polygon = None

		self.canvas.draw()


	def auto_extend_view(self, x, y):
		# always make just placed vertex at the center of the view
		cur_xmin, cur_xmax = self.axis.get_xlim()
		cur_ybottom, cur_ytop = self.axis.get_ylim()

		if abs(x - cur_xmin) < AUTO_EXTEND_VIEW_TOLERANCE or abs(x - cur_xmax) < AUTO_EXTEND_VIEW_TOLERANCE:
			cur_xcenter = cur_xmin * .6 + cur_xmax * .4 if abs(x - cur_xmin) < AUTO_EXTEND_VIEW_TOLERANCE else cur_xmin * .4 + cur_xmax * .6
			translation_x = cur_xcenter - x

			self.axis.set_xlim([cur_xmin-translation_x, cur_xmax-translation_x])
			self.axis2.set_xlim([cur_xmin-translation_x, cur_xmax-translation_x])

		if abs(y - cur_ybottom) < AUTO_EXTEND_VIEW_TOLERANCE or abs(y - cur_ytop) < AUTO_EXTEND_VIEW_TOLERANCE:
			cur_ycenter = cur_ybottom * .6 + cur_ytop * .4 if abs(y - cur_ybottom) < AUTO_EXTEND_VIEW_TOLERANCE else cur_ybottom * .4 + cur_ytop * .6
			translation_y = cur_ycenter - y

			self.axis.set_ylim([cur_ybottom-translation_y, cur_ytop-translation_y])
			self.axis2.set_ylim([cur_ybottom-translation_y, cur_ytop-translation_y])

		self.canvas.draw()
		self.canvas2.draw()

	# def place_vertex(self, x, y):
	#     self.selected_proposal_vertices.append([x, y])

	#     # curr_vertex_circle = plt.Circle((x, y), radius=10, color=self.colors[self.curr_label + 1], alpha=.8)
	#     curr_vertex_circle = plt.Circle((x, y), radius=UNSELECTED_CIRCLE_SIZE, color=self.boundary_colors[1], alpha=.8)
	#     self.axis.add_patch(curr_vertex_circle)
	#     self.selected_proposal_vertexCircles.append(curr_vertex_circle)

	#     curr_vertex_circle.set_picker(CIRCLE_PICK_THRESH)

	#     self.auto_extend_view(x, y)

	#     self.canvas.draw()
	#     self.canvas2.draw()

	#     self.history.append({'type': 'add_vertex', 'selected_proposal_vertexCircles': self.selected_proposal_vertexCircles,
	#         'selected_proposal_vertices': self.selected_proposal_vertices})

		# print self.history


	def load_segmentation(self):
		sys.stderr.write('loading segmentation...\n')
		self.statusBar().showMessage('loading segmentation...')

		self.dm.load_multiple_results(results=[
		  'segmentation', 
		  'edgeEndpoints', 'edgeMidpoints'])
		self.segmentation = self.dm.load_pipeline_result('segmentation')
		self.n_superpixels = self.dm.segmentation.max() + 1

		self.seg_loaded = True
		sys.stderr.write('segmentation loaded.\n')

		sys.stderr.write('loading sp props...\n')
		self.statusBar().showMessage('loading sp properties..')
		# self.sp_centroids = self.dm.load_pipeline_result('spCentroids')
		# self.sp_bboxes = self.dm.load_pipeline_result('spBbox')
		sys.stderr.write('sp properties loaded.\n')

		self.statusBar().showMessage('')

		# self.sp_rectlist = [None for _ in range(self.dm.n_superpixels)]


	def turn_superpixels_off(self):
		self.statusBar().showMessage('Supepixels OFF')

		self.buttonSpOnOff.setText('Turn Superpixels ON')

		self.segm_handle.remove()
		self.superpixels_on = False
		
		# self.axis.imshow(self.masked_img, cmap=plt.cm.Greys_r,aspect='equal')
		# self.orig_image_handle = self.axis.imshow(self.masked_img, aspect='equal')

	def turn_superpixels_on(self):
		self.statusBar().showMessage('Supepixels ON')

		self.buttonSpOnOff.setText('Turn Superpixels OFF')

		if self.segm_transparent is None:
			self.segm_transparent = self.dm.load_pipeline_result('segmentationTransparent')
			self.my_cmap = plt.cm.Reds
			self.my_cmap.set_under(color="white", alpha="0")

		if not self.seg_loaded:
			self.load_segmentation()

		self.superpixels_on = True
		
		if hasattr(self, 'segm_handle'):
			self.segm_handle.set_data(self.segm_transparent)
		else:
			self.segm_handle = self.axis.imshow(self.segm_transparent, aspect='equal', 
								cmap=self.my_cmap, alpha=1.)



	def cancel_current_circle(self):

		if self.selected_circle is not None:
			self.selected_circle.set_radius(UNSELECTED_CIRCLE_SIZE)
			self.selected_circle = None

			self.selected_vertex_index = None


	def cancel_current_selection(self):

		if self.selected_proposal_polygon is not None:

			if self.selected_proposal_polygon.get_linewidth() != UNSELECTED_POLYGON_LINEWIDTH:
				self.selected_proposal_polygon.set_linewidth(UNSELECTED_POLYGON_LINEWIDTH)

			if self.selected_proposal_polygon in self.axis.patches:
				if self.selected_proposal_polygon not in self.accepted_proposals:
					self.selected_proposal_polygon.remove()
					for vertex_circ in self.selected_proposal_vertexCircles:
						vertex_circ.remove()

		self.selected_proposal_polygon = None
		self.selected_proposal_vertexCircles = None

		if self.selected_circle is not None:
			self.selected_circle.set_radius(UNSELECTED_CIRCLE_SIZE)
			self.selected_circle = None

			self.selected_vertex_index = None


		self.canvas.draw()


	def mode_changed(self):

		self.cancel_current_selection()

		if self.radioButton_globalProposal.isChecked():

			self.shuffle_global_proposals = True

			if not self.superpixels_on:
				self.turn_superpixels_on()

			if not hasattr(self, 'global_proposal_tuples'):
				self.load_global_proposals()

		elif self.radioButton_localProposal.isChecked():

			if not self.superpixels_on:
				self.turn_superpixels_on()

			self.shuffle_global_proposals = False

			if not hasattr(self, 'local_proposal_tuples'):
				self.load_local_proposals()

		self.canvas.draw()

	def display_option_changed(self):
		if self.sender() == self.buttonSpOnOff:

			if not self.superpixels_on:
				self.turn_superpixels_on()
			else:
				self.turn_superpixels_off()
		else:
			print 'not implemented'
			return

			# if self.under_img is not None:
			#   self.under_img.remove()

			self.axis.clear()

			if self.sender() == self.img_radioButton:

				# self.axis.clear()
				# self.axis.axis('off')

				# self.under_img = self.axis.imshow(self.masked_img, aspect='equal', cmap=plt.cm.Greys_r)
				self.axis.imshow(self.dm.image_rgb_jpg, aspect='equal', cmap=plt.cm.Greys_r)
				# self.superpixels_on = False

			elif self.sender() == self.textonmap_radioButton:

				# self.axis.clear()
				# self.axis.axis('off')

				if self.textonmap_vis is None:
					self.textonmap_vis = self.dm.load_pipeline_result('texMapViz')

				# if self.under_img is not None:
				#   self.under_img.remove()

				# self.under_img = self.axis.imshow(self.textonmap_vis, cmap=plt.cm.Greys_r, aspect='equal')
				self.axis.imshow(self.textonmap_vis, cmap=plt.cm.Greys_r, aspect='equal')
				# self.superpixels_on = False

			elif self.sender() == self.dirmap_radioButton:

				# self.axis.clear()
				# self.axis.axis('off')

				if self.dirmap_vis is None:
					self.dirmap_vis = self.dm.load_pipeline_result('dirMap', 'jpg')
					self.dirmap_vis[~self.dm.mask] = 0


				# self.under_img = self.axis.imshow(self.dirmap_vis, aspect='equal')
				self.axis.imshow(self.dirmap_vis, aspect='equal')

				# if not self.seg_loaded:
				#   self.load_segmentation()

				# self.superpixels_on = False

			# elif self.sender() == self.labeling_radioButton:
			#   pass

		self.axis.axis('off')
		# self.axis.set_xlim([self.newxmin, self.newxmax])
		# self.axis.set_ylim([self.newymin, self.newymax])
		# self.fig.subplots_adjust(left=0, bottom=0, right=1, top=1)
		self.canvas.draw()

		self.axis2.axis('off')
		# self.axis2.set_xlim([self.newxmin, self.newxmax])
		# self.axis2.set_ylim([self.newymin, self.newymax])
		# self.fig2.subplots_adjust(left=0, bottom=0, right=1, top=1)
		self.canvas2.draw()

			   
if __name__ == "__main__":
	from sys import argv, exit
	appl = QApplication(argv)

	stack = sys.argv[1]
	m = BrainLabelingGUI(stack=stack)

	m.showMaximized()
	m.raise_()
	exit(appl.exec_())
