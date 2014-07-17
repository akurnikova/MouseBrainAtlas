"""This module implements a controller-viewer-model framework for
interactive labeling of brain slice images.

The initial input to the module is an image stack and for each image
in the stack: a segmentation map
"""

import numpy as np
import pylab
from matplotlib import pyplot as plt
import matplotlib as mpl
import Image
from time import time

class PBN_View:
    """ class for displaying the gui, coloring super-pixels, and capturing user input """
    def __init__(seg,seg_no,image):
        this.seg=seg
        this.image=image
        #labels is an array of pairs:
        # 0: the machine generated label (color number) of the super-pixel
        # 1: the human generate label

        this.labels=np.zeros([seg_no,2]) 

class PBN_State:
    """ Class for storing a labeling of the super-pixels and whether they are human labeled or 
    Machine generated
    """

class PBN_control:
    """Class for controlling the interactive labeling of super-pixels
    using information from the user, from the texton distributions, and
    from the directionality maps.
    """