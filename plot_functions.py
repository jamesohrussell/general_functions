#==================================================================
# Create a truncated colormap  
#==================================================================
#
# 
# James Russell 2020
#
# * truncate_colormap
#   - Truncates an existing colormap
#
#==================================================================
# Import libraries
#==================================================================

from numba import jit
import matplotlib.colors as colors
import numpy as np

#==================================================================
# Truncates a colormap
#==================================================================

@jit(nopython=True)
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
  """
  Create a new colormap by truncating an existing colormap. Taken 
  from: 
  https://stackoverflow.com/questions/18926031/how-to-extract-a-subset-of-a-colormap-as-a-new-colormap-in-matplotlib.
  
  Inputs:
  1) Existing colormap e.g. matplotlib.pyplot.get_cmap('jet')
  2,3) Min and max of colorbar e.g. upper 50% of colormap 
   minval=0.5, maxval=1.0
  4) Number of levels for new colormap. Default = 100.
  
  Output is a new colormap truncated from the input colormap.

  Requires matplotlib.
  """

  # Make new colormap
  new_cmap = colors.LinearSegmentedColormap.from_list(
   'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, 
   b=maxval),cmap(np.linspace(minval, maxval, n)))

  return new_cmap



