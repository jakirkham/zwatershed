import numpy as np
import sys
import time
import os
import h5py
import os.path as op
import matplotlib.cm as cm
from matplotlib.widgets import Slider, Button, RadioButtons
import matplotlib.pyplot as plt
import array

sys.path.append('..')
from zwatershed import *
from edgelist_methods import *
from visualization.visualize_funcs import *


# -------------------------------- parameters ---------------------------------------
pred_file = '/groups/turaga/home/turagas/research/caffe_v2/processed/bock2/120000/cutout_3k.h5'
gt_file = '/groups/turaga/home/turagas/research/caffe_v2/processed/bock2/120000/sample_A_x1_y1_z1_xy1.h5'

print gt_file
# f = h5py.File(data_folder + 'im_uint8.h5', 'r')
# im = f[f.keys()[0]]