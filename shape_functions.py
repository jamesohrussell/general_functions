#==================================================================
# Shape functions
#==================================================================
#
# James Russell, University of Utah, 2021
#
# * label_wdiags
#   - Replicates the scipy.ndimage.label() function but 
#      connects features at diagonals
#
# * find_corners
#   - Finds all x,y coordinates of all corners given a set of 
#      pixel center coordinates
#
# * points_in_shape
#   - Finds all x,y coordinates within a shape as defined by its
#      vertices
#
# * fit_ellipse_svd_earth
#   - Uses singular value decomposition to define an ellipse given 
#      a set of points in 2d
#
# * periodic_cmass_earth
#   - Calculates a center of mass for a set of points on earth 
#      with the zonal direction periodic
#
# * convert_cendirlen_latlon_earth
#   - Takes a central point, direction/bearing, and length in km,
#      and converts it to a line with lat,lon coordinates for the 
#      two ends of the line
#
#==================================================================
# Import libraries
#==================================================================

import scipy.ndimage as ndimg
import numpy as np
from matplotlib.path import Path
import geophys_functions as gfns
import misc_functions as mfns
from geopy.distance import geodesic
import matplotlib.pyplot as plt
import cmaps

#==================================================================
# Identify contiguous areas in 2d field 
#==================================================================

def label_wdiags(array):
  """
  This is a small adjustment to the scipy.ndimage.label() function
   to connect pixels at diagonals as the same feature. Original 
   scipy.ndimage.label() documentation can be found here: 
   https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.label.html

  Input:
   1) A 2D list, numpy array, or pandas dataframe with the 
    data you wish to identify features in.

  Output is a 2D array corresponding to data but with features 
   labelled as integers.

  Requires scipy 1.2.1 (conda install -c anaconda scipy; 
   https://pypi.org/project/scipy/)
  """

  # Check array
  ndims = len(np.shape(array))

  # Setup a binary array to account for diagonals given 
  #  dimensions of array
  s = np.ones([3]*ndims)

  # Run scipy.ndimage.label() with binary structure
  labels, numL = ndimg.label(array, structure=s)

  # Returns array with labels and the number of features
  return(labels,numL)



#==================================================================
# Find unique corner points from a list of pixel centers
#==================================================================

def find_corners(coords,dx,dy,dp=2):
  """
  Defines a list of coordinates for corners given a list of pixel
   centers.

  Inputs:
   1) A list of coordinate pairs i.e. [[x1,y1],[x2,y2],...]
   2,3) Scalars defining the width and height of the pixels
   4) A scalar defining the decimal places to round the resulting 
    coordinates. Optional, default = 2.

  Output is a list of coordinate pairs for the corners.

  Requires numpy 1.16.3 (conda install -c anaconda numpy; 
   https://pypi.org/project/numpy/)

  Note, only returns unique corner coordinates (i.e. there are 
   no coordinates that are the same even if two pixels share a 
   corner).
  """

  # Define the distances to a corner from center
  dx2=dx/2; dy2=dy/2

  # Define all corners
  allcorners = [(round(x+dx2,dp),round(y+dy2,dp))
                              for x,y in coords] + \
               [(round(x-dx2,dp),round(y+dy2,dp)) 
                              for x,y in coords] + \
               [(round(x+dx2,dp),round(y-dy2,dp)) 
                              for x,y in coords] + \
               [(round(x-dx2,dp),round(y-dy2,dp)) 
                              for x,y in coords]

  # Get rid of all doubles and return list of corner coordinates
  return([list(rows) for rows in 
          np.unique(allcorners, axis=0)])



#==================================================================
# Find which points are inside a shape
#==================================================================

def points_in_shape(verts,points,widen=0):
  """
  Finds and returns all grid points within a list of vertices 
   defining a shape. 

  Input: 
   1) A list of tuples of the x,y coordinates of the vertices
   2) A list of tuples of the x,y coordinates of the points to 
       check whether they are in the shape
   3) The amount to widen the shape by. Most common use - due to
       computational issues, sometimes the shape doesn't include
       points on the edge. By widening the shape slightly, you 
       can ensure points are included. Positive values widen the 
       shape when the vertcies are ordered counterclockwise and
       negative widens the shape when the vertices are ordered
       clockwise. Default is 0 - no widening.

  Output 

  Requires maplotlib 3.0.3 (conda install -c conda-forge 
   matplotlib; https://pypi.org/project/matplotlib/)
  """

  # Define the shape
  p = Path(verts)

  # Get list of booleans for if inside array
  grid = p.contains_points(points,radius=widen)

  # Return list of coordinate tuples for points inside shape
  return([points[i] for i, x in enumerate(grid) if x])



#==================================================================
# Fit an ellipse to a set of 2D data points with SVD
#==================================================================

def fit_ellipse_svd_earth(x,y,center,fit=False,plot=False):
  """
  Fit an ellipse to a set of 2D data points with singular value 
   decomposition (SVD)

  Input: 
   1,2) Lists of x and y coordinates to fit the ellipse to
   3) Center lon and lat as pair in tuple
   4) Output the 1000 point ellipse fit
   5) Make a figure

  Output:
   1) A list with x and y location of ellipse center
   2) A list with major and minor axes angle from north
   3) A list with major and minor axes geographic distances

  Requires numpy 1.16.3 (conda install -c anaconda numpy; 
   https://pypi.org/project/numpy/)
  """

  # Get number of points
  N = len(x)

  # Get distance to center for all points in cartesian space
  xc = [geodesic((center[1],ln),(center[1],center[0])).m if ln>center[0] else
       -geodesic((center[1],ln),(center[1],center[0])).m if ln<center[0] else
       0. for ln in x]
  yc = [geodesic((lt,center[0]),(center[1],center[0])).m if lt>center[1] else
       -geodesic((lt,center[0]),(center[1],center[0])).m if lt<center[1] else
       0. for lt in y]

  # Do singular value decomposition
  U,S = np.linalg.svd(np.stack((xc, yc)))[0:2]

  # Calculate length in degrees of axes
  axlen = [2*np.sqrt(2/N)*l for l in S]

  # Calculate angle of axes
  tt = np.linspace(0, 2*np.pi, 5)
  circle = np.stack((np.cos(tt), np.sin(tt)))
  transform = np.sqrt(2/N) * U.dot(np.diag(S))
  fitxy = transform.dot(circle)
  oa   = [[fitxy[0][0],fitxy[1][0]],
          [fitxy[0][1],fitxy[1][1]]]
  axdir = [np.degrees(np.arctan(fitxy[0][0]/fitxy[1][0])),
           np.degrees(np.arctan(fitxy[0][1]/fitxy[1][1]))]
  axdir = [d-360 if d>180 else d for d in axdir]

  # Ensure major and minor axes in correct order
  if axlen[0]<axlen[1]:
    axlen = axlen[::-1]
    axdir = axdir[::-1]

  # Get 1000 point data fit for the ellipse
  if fit or plot:
    # Define a unit circle
    tt = np.linspace(0, 2*np.pi, 1000)
    circle = np.stack((np.cos(tt), np.sin(tt)))

    # Define transformation matrix
    transform = np.sqrt(2/N) * U.dot(np.diag(S))
    fitxy = transform.dot(circle)

  # Plot ellipse and data
  if plot:

    # Make plot with data
    print("Plotting ellipse")
    import matplotlib.pyplot as plt
    plt.plot(xc, yc, '.')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(True)

    # Plot ellipse
    plt.plot(fitxy[0, :], fitxy[1, :], 'r')

    # Plot center
    plt.plot(0,0,"r")

    # Plot major axis
    x1 = -((axlen[0]/2)*np.sin(np.deg2rad(axdir[0])))
    x2 = +((axlen[0]/2)*np.sin(np.deg2rad(axdir[0])))
    y1 = -((axlen[0]/2)*np.cos(np.deg2rad(axdir[0])))
    y2 = +((axlen[0]/2)*np.cos(np.deg2rad(axdir[0])))
    plt.plot([x1,x2],[y1,y2])

    # Plot minor axis
    x1 = -((axlen[1]/2)*np.sin(np.deg2rad(axdir[1])))
    x2 = +((axlen[1]/2)*np.sin(np.deg2rad(axdir[1])))
    y1 = -((axlen[1]/2)*np.cos(np.deg2rad(axdir[1])))
    y2 = +((axlen[1]/2)*np.cos(np.deg2rad(axdir[1])))
    plt.plot([x1,x2],[y1,y2])

    # Show plot
    plt.show()

  if fit:
    # Also return fit
    return(axdir,axlen,fitxy)
  else:
    # Return center, and direction and length of major and minor
    #  axes
    return(axdir,axlen)


#==================================================================
# Peridiodic center of mass on earth calculation
#==================================================================

def periodic_cmass_earth(lon):
  """
  Calculates the center of mass for a periodic domain. Input is x
   but this can be repeated for any coordinate if you have multiple
   periodic coordinates.

  Input: 
   1) Lists of 1D coordinates in -180->179.999... or 0->359.999...

  Output:
   1) Center of mass location in same coordinate system as input.

  Requires numpy 1.16.3 (conda install -c anaconda numpy; 
   https://pypi.org/project/numpy/)

  Method from:
  https://en.wikipedia.org/wiki/Center_of_mass
  """

  # Calculate variables
  xiibar   = np.cos(np.radians(lon))
  zetaibar = np.sin(np.radians(lon))
  
  # Check coordinate system and return center of mass in same
  #  coordinate system
  cmass = np.degrees(np.pi+
   np.arctan2(-np.mean(zetaibar),-np.mean(xiibar)))
  if min(lon)<0:
    if 180<cmass<360:
      return(cmass-360)
    elif cmass==180:
      return(cmass-360)
    else:
      return(cmass)
  else:
    if cmass==360:
      return(cmass-360)
    else:
      return(cmass)

#==================================================================
# Plot PF and ellipse
#==================================================================

def convert_cendirlen_latlon_earth(center,direction,length):
  """
  Converts a line with bearing, length, and a center to lat and 
   lon coordinates.

  Input: 
   1) center of mass for data [lon, lat] converted to radians
   2) direction/bearing (as clockwise angle from north)
   3) lengths (km)

  Outputs lat and lon points at ends of line

  Requires numpy 1.16.3 (conda install -c anaconda numpy; 
   https://pypi.org/project/numpy/)
  """

  # Constants
  R = 6378.1 #Radius of the Earth km

  # Convert direction to radians
  direction = np.radians(direction)

  # Get northern most latitude of line
  lat2 = np.degrees(
   np.arcsin(np.sin(center[1])*np.cos(length*.0005/R) + \
   np.cos(center[1])*np.sin(length*.0005/R)*\
   np.cos(direction)))

  # Get southern most latitude of line
  lat1 = np.degrees(
   np.arcsin(np.sin(center[1])*np.cos(length*.0005/R) + \
   np.cos(center[1])*np.sin(length*.0005/R)*\
   np.cos(direction-np.pi)))

  # Get eastern and western most longitude of line
  faclon = np.arctan2(np.sin(direction)*\
   np.sin(length*.0005/R)*np.cos(center[1]),\
   np.cos(length*.0005/R)-np.sin(center[1])*np.sin(center[1]))
  lon2 = np.degrees(center[0] + faclon)
  lon1 = np.degrees(center[0] - faclon)

  return(lon2,lon1,lat2,lat1)

#==================================================================
# End functions
#==================================================================
