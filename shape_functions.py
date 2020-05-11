#==================================================================
# Shape functions
#==================================================================
#
# Functions to calculate various shape parameters
# James Russell 2020
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
# * fit_ellipse_svd
#   - Uses singular value decomposition to define an ellipse given 
#      a set of points in 2d
#
# * plot_pf_ellipse
#   - Uses output from the above function to plot the PF and the 
#      fitted ellipse
#
#==================================================================
# Import libraries
#==================================================================

import scipy.ndimage as ndimg
import numpy as np
from matplotlib.path import Path
import geophys_functions as gfns

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

  # Setup a 2x2 binary array to account for diagonals
  s = ndimg.generate_binary_structure(2,2)

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

def fit_ellipse_svd(x,y,fit=False,plot=False):
  """
  Fit an ellipse to a set of 2D data points with singular value 
   decomposition (SVD)

  Input: 
   1,2) Lists of x and y coordinates to fit the ellipse to

  Output:
   1) A list with x and y location of ellipse center
   2) A list with major and minor axes angle from north
   3) A list with major and minor axes geographic distances

  Requires numpy 1.16.3 (conda install -c anaconda numpy; 
   https://pypi.org/project/numpy/)
  """

  # Get number of points
  N = len(x)

  # Calculate center of mass of shape
  center = [np.mean(x),np.mean(y)]

  # Get distance to center for all points
  xc = x - center[0]
  yc = y - center[1]

  # Do singular value decomposition
  U, S = np.linalg.svd(np.stack((xc, yc)))[0:2]

  # Calculate angle of axes
  axdir = [gfns.calc_direction(center[0],center[1],
           center[0]+a[0],center[1]+a[1]) for a in U]
  axdir = [d-180 if d>180 else d for d in axdir]

  # Calculate length in degrees of axes
  axlen = [2*np.sqrt(2/N)*l for l in S]

  # Calculate geopgraphic distance of axes
  axdist =  [gfns.calc_distance(center[0] - ((axlen[i]/2)* \
                       np.sin(np.deg2rad(axdir[i]))), \
                       center[1] - ((axlen[i]/2)* \
                       np.cos(np.deg2rad(axdir[i]))), \
                       center[0] + ((axlen[i]/2)* \
                       np.sin(np.deg2rad(axdir[i]))), \
                       center[1] + ((axlen[i]/2)* \
                       np.cos(np.deg2rad(axdir[i])))) \
             for i in range(len(axlen))]

  # Ensure major and minor axes in correct order
  if axdist[0]<axdist[1]:
    axdist = axdist[::-1]
    axdir  = axdir[::-1]

  # Get 1000 point data fit for the ellipse
  if fit or plot:
    # Define a unit circle
    tt = np.linspace(0, 2*np.pi, 1000)
    circle = np.stack((np.cos(tt), np.sin(tt)))
  
    # Define transformation matrix
    transform = np.sqrt(2/N) * U.dot(np.diag(S))
    fitxy = transform.dot(circle)+ \
     np.array([[center[0]],[center[1]]])

  # Plot ellipse and data
  if plot:

    # Make plot with data
    print("Plotting ellipse")
    import matplotlib.pyplot as plt
    plt.plot(x, y, '.')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(True)

    # Plot ellipse
    plt.plot(fitxy[0, :], fitxy[1, :], 'r')

    # Plot center
    plt.plot(center[0],center[1],".r")

    # Plot major axis
    x1 = center[0]-((axlen[0]/2)*np.sin(np.deg2rad(axdir[0])))
    x2 = center[0]+((axlen[0]/2)*np.sin(np.deg2rad(axdir[0])))
    y1 = center[1]-((axlen[0]/2)*np.cos(np.deg2rad(axdir[0])))
    y2 = center[1]+((axlen[0]/2)*np.cos(np.deg2rad(axdir[0])))
    plt.plot([x1,x2],[y1,y2])

    # Plot minor axis
    x1 = center[0]-((axlen[1]/2)*np.sin(np.deg2rad(axdir[1])))
    x2 = center[0]+((axlen[1]/2)*np.sin(np.deg2rad(axdir[1])))
    y1 = center[1]-((axlen[1]/2)*np.cos(np.deg2rad(axdir[1])))
    y2 = center[1]+((axlen[1]/2)*np.cos(np.deg2rad(axdir[1])))
    plt.plot([x1,x2],[y1,y2])

    # Show plot
    plt.show()

  if fit:
    # Also return fit
    return(center,axdir,axlen,fitxy)
  else:
    # Return center, and direction and length of major and minor
    #  axes
    return(center,axdir,axdist)


#==================================================================
# Plot PF and ellipse
#==================================================================

def plot_pf_ellipse(x,y,z,center,axdir,axlen,fit):
  """
  

  Input: 
   1,2) 

  Output:
   1) 
   2) 
   3) 

  Requires numpy 1.16.3 (conda install -c anaconda numpy; 
   https://pypi.org/project/numpy/)
  """

  # Make plot with data
  print("Plotting ellipse")
  import matplotlib.pyplot as plt
  import nclcmaps
  cmap = nclcmaps.cmap('WhiteBlueGreenYellowRed')
  cs = plt.pcolormesh(x, y, z,cmap=cmap,vmin=0,vmax=16)
  cbar = plt.colorbar(cs, pad=.1, fraction=0.06)   
  plt.gca().set_aspect('equal', adjustable='box')
  plt.grid(True)

  # Plot ellipse
  plt.plot(fit[0, :], fit[1, :], 'r')

  # Plot center
  plt.plot(center[0],center[1],".r")

  # Plot major axis
  x1 = center[0]-((axlen[0]/2)*np.sin(np.deg2rad(axdir[0])))
  x2 = center[0]+((axlen[0]/2)*np.sin(np.deg2rad(axdir[0])))
  y1 = center[1]-((axlen[0]/2)*np.cos(np.deg2rad(axdir[0])))
  y2 = center[1]+((axlen[0]/2)*np.cos(np.deg2rad(axdir[0])))
  plt.plot([x1,x2],[y1,y2])

  # Plot minor axis
  x1 = center[0]-((axlen[1]/2)*np.sin(np.deg2rad(axdir[1])))
  x2 = center[0]+((axlen[1]/2)*np.sin(np.deg2rad(axdir[1])))
  y1 = center[1]-((axlen[1]/2)*np.cos(np.deg2rad(axdir[1])))
  y2 = center[1]+((axlen[1]/2)*np.cos(np.deg2rad(axdir[1])))
  plt.plot([x1,x2],[y1,y2])

  # Show plot
  plt.show()

#==================================================================
# Plot PF and ellipse
#==================================================================

def plot_pf(x,y,z):
  """
  

  Input: 
   1,2) 

  Output:
   1) 
   2) 
   3) 

  Requires numpy 1.16.3 (conda install -c anaconda numpy; 
   https://pypi.org/project/numpy/)
  """

  # Make plot with data
  print("Plotting ellipse")
  import matplotlib.pyplot as plt
  import nclcmaps
  cmap = nclcmaps.cmap('WhiteBlueGreenYellowRed')
  cs = plt.pcolormesh(x, y, z,cmap=cmap,vmin=0,vmax=16)
  cbar = plt.colorbar(cs, pad=.1, fraction=0.06)   
  plt.gca().set_aspect('equal', adjustable='box')
  plt.grid(True)

  # Show plot
  plt.show()



