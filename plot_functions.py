#==================================================================
# Plotting functions
#==================================================================
#
# James Russell, University of Utah, 2021
#
# * truncate_colormap
#   - Truncates an existing colormap
#
# * plot_pixel_axes_earth
#   - Plots pixel data on lat,lon grid and uses output from
#      shape_functions.fit_ellipse_svd_earth to plot fitted 
#      major and minor axes
#
# * plot_pixel
#   - Plots pixel level data
#
#==================================================================
# Import libraries
#==================================================================

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import shape_functions as sfns
import earth_functions as efns
import cmaps

#==================================================================
# Truncates a colormap
#==================================================================

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

#==================================================================
# Plot pixel data (optional addition of major and minor axes)
#==================================================================

def plot_pixel_axes_earth(x,y,z,center=None,axdir=None,axlen=None):
  """
  Plots pixel data along with the major and minor axes fitted using 
   fit_ellipse_svd_earth

  Input: 
   1,2) x (lon) and y (lat) coordinates of pixels as meshgrid
   3) values of pixels as grid (corresponding to x and y meshgrids)
   4) center of mass for data [lon, lat]
   5) major and minor axes orientation (as clockwise angle from north)
    as a list: [major axes orientation, minor axes orientation]
   6) major and minor axes lengths (km) as a list: 
    [major axes length, minor axes length]

  Shows figure with pixels and axes

  Requires numpy 1.16.3 (conda install -c anaconda numpy; 
   https://pypi.org/project/numpy/)
  """

  # Make plot with data
  cs = plt.pcolormesh(x, y, z,cmap=cmaps.WhiteBlueGreenYellowRed,
  vmin=0,vmax=17,shading="auto")
  cbar = plt.colorbar(cs, pad=.1, fraction=0.06)   
  plt.gca().set_aspect('equal', adjustable='box')
  plt.grid(True)

  # if inputs are given, plot axes too
  if center is not None:

    clon = np.radians(center[0])
    clat = np.radians(center[1])

    lon2mj,lon1mj,lat2mj,lat1mj = \
     efns.convert_cendirlen_latlon(
     [clon,clat],axdir[0],axlen[0])

    lon2mn,lon1mn,lat2mn,lat1mn = \
     efns.convert_cendirlen_latlon(
     [clon,clat],axdir[1],axlen[1])

    # Plot major axis
    plt.plot([lon1mj,lon2mj],[lat1mj,lat2mj],"k",linewidth=3)
    plt.plot([lon1mn,lon2mn],[lat1mn,lat2mn],"k",linewidth=3)

  # Show plot
  plt.show()

#==================================================================
# End functions
#==================================================================
