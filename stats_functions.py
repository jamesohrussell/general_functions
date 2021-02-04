#==================================================================
# Time functions
#==================================================================
# 
# Functions for statistics
# James Russell 2020
#
#==================================================================
# Import libraries
#==================================================================

import numpy as np

#==================================================================
# Calculate time since given units
#==================================================================

def rolling_window(arr, window):
  """
  Rolling window to run statistics on

  Input: 
   1) Array
   2) Integer indicating window size

  Output:
   

  Requires numpy
  """

  shape = arr.shape[:-1] + (arr.shape[-1] - window + 1, window)
  strides = arr.strides + (arr.strides[-1],)
  return np.lib.stride_tricks.as_strided(arr, shape=shape, 
   strides=strides)


