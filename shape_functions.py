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
#==================================================================
# Import libraries
#==================================================================

import scipy.ndimage as ndimg
import numpy as np
from matplotlib.path import Path
import matplotlib.pyplot as plt
from geopy.distance import geodesic
import earth_functions as efns
import misc_functions as mfns

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

def points_in_shape(verts,points,widen=0,output="coordsTrue"):
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
       shape when the vertices are ordered counterclockwise and
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
  if output=="coordsTrue":
    return([points[i] for i, x in enumerate(grid) if x])
  elif output=="TrueFalse":
    return(grid)


#==================================================================
# Fit an ellipse to a set of 2D data points with SVD
#==================================================================

def fit_ellipse_svd_earth(x,y,center,goodness=False,
  dx=None,dy=None,plot=False):
  """
  Fit an ellipse to a set of 2D data points with singular value 
   decomposition (SVD)

  Input: 
   1,2) x,y: Lists of x and y coordinates to fit the ellipse to
   3) center: Center lon and lat as pair in tuple
   Optional:
   goodness: Calculate goodness of the ellipse fit. Requires dx,dy.
   dx,dy: For the goodness of fit for the ellipse, we calculate 
    the number of points within the ellipse and divide this by
    the number of points possible within the ellipse. The 
    algorithm takes the max and min x,y points and then calculates
    a grid based on a grid spacing dx,dy for this.
   plot: Make a figure

  Output:
   1) A list with major and minor axes angle from north
   2) A list with major and minor axes geographic distances
   Optional:
   3) goodness (scalar): This is adapted from a typical goodness
      of fit for an ellipse to account for a filled ellipse. In a
      typical goodness of fit, it is simply the sum of the 
      distance to the ellipse edge. However, this treats points in
      the ellipse equally to those outside which is undesirable. 
      The goodness here is calculated as follows:
       fill factor = (n points in ellipse / n points possible
        within ellipse) * (n points in ellipse/ total n points)
       residual outside = (1 - (mean of distances between points 
        outside ellipse to nearest point on ellipse / mean radius
        of ellipse )) * (n points outside ellipse/ total n points)
       goodness = fill factor * residual outside

  Requires numpy 1.16.3 (conda install -c anaconda numpy; 
   https://pypi.org/project/numpy/)
  """

  # Get number of points
  N = len(x)

  # Get distance to center for all points in cartesian space
  xc = [geodesic((center[1],ln),(center[1],center[0])).m
       if ln>center[0] else
       -geodesic((center[1],ln),(center[1],center[0])).m 
       if ln<center[0] else 0. for ln in x]
  yc = [geodesic((lt,center[0]),(center[1],center[0])).m 
       if lt>center[1] else
       -geodesic((lt,center[0]),(center[1],center[0])).m 
       if lt<center[1] else 0. for lt in y]

  # Do singular value decomposition
  U,S = np.linalg.svd(np.stack((xc, yc)))[0:2]

  # Calculate length of axes in degrees
  axlen = [2*np.sqrt(2/N)*l for l in S]

  # Calculate angle of axes in -180-180
  tt = np.linspace(0, 2*np.pi, 5)
  circle = np.stack((np.cos(tt), np.sin(tt)))
  transform = np.sqrt(2/N) * U.dot(np.diag(S))
  fitxy = transform.dot(circle)
  oa   = [[fitxy[0][0],fitxy[1][0]],
          [fitxy[0][1],fitxy[1][1]]]
  axdir = [np.degrees(np.arctan(fitxy[0][0]/fitxy[1][0])),
           np.degrees(np.arctan(fitxy[0][1]/fitxy[1][1]))]
  axdir = [d-360 if d>180 else d for d in axdir]

  # Ensure major and minor axes are in correct order
  if axlen[0]<axlen[1]:
    axlen = axlen[::-1]
    axdir = axdir[::-1]

  # Get 50 point data fit for the ellipse
  if plot or goodness:

    # Define a unit circle
    tt = np.linspace(0, 2*np.pi, 50)
    circle = np.stack((np.cos(tt), np.sin(tt)))

    # Define transformation matrix
    transform = np.sqrt(2/N) * U.dot(np.diag(S))
    fitxy  = transform.dot(circle)

  # Get goodness of fit
  if goodness:

    # Get vertices of ellipse
    verts  = list(zip(fitxy[0],fitxy[1]))

    # Get all possible grid points in lat lon
    nx = [round(i,2) for i in 
     np.arange(min(x),max(x)+dx,dx)]
    ny = [round(i,2) for i in 
     np.arange(min(y),max(y)+dy,dy)]
    grid = np.meshgrid(nx,ny)

    # Get cartesian coordinates of all grid points
    xg = [geodesic((center[1],ln),
     (center[1],center[0])).m if ln>center[0] else
     -geodesic((center[1],ln),(center[1],center[0])).m 
     if ln<center[0] else 0. for ln in 
     list(np.array(grid[0]).flatten())]
    yg = [geodesic((lt,center[0]),
     (center[1],center[0])).m if lt>center[1] else
     -geodesic((lt,center[0]),(center[1],center[0])).m 
     if lt<center[1] else 0. for lt in
     list(np.array(grid[1]).flatten())]
    glls = list(zip(xg,yg))

    # Get first factor (fraction of points within ellipse 
    #  that are filled scaled by the fraction of all the 
    #  points they represent)
    points = list(zip(xc,yc))
    pointsin = points_in_shape(verts,points,output="TrueFalse")
    fac1 = sum(pointsin)**2/\
     (len(points_in_shape(verts,glls))*len(pointsin)) 

    # Get second factor ()
    pointsout = [points[i] for i, x in enumerate(pointsin) 
     if not x]     
    fac2 = (1-np.mean([min([mfns.cartesian_distance(
     p[0],v[0],p[1],v[1])
     for v in verts]) for p in pointsout])/(.5*np.mean(axlen)))\
     *(sum(pointsin==False)/len(pointsin))
 
    # Get goodness factor
    goodness = fac1+fac2

  # Plot ellipse and data
  if plot:

    # Make plot with data
    print("Plotting ellipse")
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

  if goodness:
    # Also return fit
    return(axdir,axlen,goodness)
  else:
    # Return center, and direction and length of major and minor
    #  axes
    return(axdir,axlen)


#==================================================================
# End functions
#==================================================================
