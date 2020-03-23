#==========================================================
# Check a PF for upscale growth
#==========================================================

def upscale_growth_check(area,volr):
  "Takes precipitation feature time series of area and volumetric rain rate as inputs and outputs a boolean if precipitation feature passes all tests indicating whether upscale growth occurred or not. True indicates upscale growth occurred, False indicates no upscale growth."

  # Check if system does not last long enough (must last 3 hours)
  if len(volr)<6:
    print("Ignoring object because it does not live for long enough")
    return(False)

  # Check if PF reaches a minimum threshold to be defined as an MCS
  if max(area)<1000:
    print("Ignoring object because it does not reach a threshold size for an MCS")
    return(False)

  # Find indices of maxima in times series
  idta = (np.abs(area-max(area))).argmin()
  idtv = (np.abs(volr-max(volr))).argmin()

  # Check if object maximum is at start
  if idta==0 or idtv==0:
    print("Ignoring object because it's already past it's prime")
    return(False)

  # Check if area change is large enough
  if max(area)/min(area[0:idta])<2:
    print("Ignoring object because it does not double in area")
    return(False)

  # Check if VRR change is large enough
  if max(volr)/min(volr[0:idtv])<2:
    print("Ignoring object because it does not double in VRR")
    return(False)

  # If all tests passed, output True
  return(True)

#==========================================================
# Fit Ellipse Function
#==========================================================

import numpy as np
from numpy.linalg import eig, inv

def fitEllipse(x,y):
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T,D)
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    E, V =  eig(np.dot(inv(S), C))
    n = np.argmax(np.abs(E))
    a = V[:,n]
    return a

def ellipse_center(a):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    num = b*b-a*c
    x0=(c*d-b*f)/num
    y0=(a*f-b*d)/num
    return np.array([x0,y0])


def ellipse_angle_of_rotation( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    return 0.5*np.arctan(2*b/(a-c))


def ellipse_axis_length( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    res1=np.sqrt(up/down1)
    res2=np.sqrt(up/down2)
    return np.array([res1, res2])

def ellipse_angle_of_rotation2( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    if b == 0:
        if a > c:
            return 0
        else:
            return np.pi/2
    else: 
        if a > c:
            return np.arctan(2*b/(a-c))/2
        else:
            return np.pi/2 + np.arctan(2*b/(a-c))/2

#==========================================================
# Least-Squares Ellipse Fitting Function
#==========================================================

import numpy
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

"""Demonstration of least-squares fitting of ellipses
    __author__ = "Ben Hammel, Nick Sullivan-Molina"
    __credits__ = ["Ben Hammel", "Nick Sullivan-Molina"]
    __maintainer__ = "Ben Hammel"
    __email__ = "bdhammel@gmail.com"
    __status__ = "Development"
    Requirements 
    ------------
    Python 2.X or 3.X
    numpy
    matplotlib
    References
    ----------
    (*) Halir, R., Flusser, J.: 'Numerically Stable Direct Least Squares 
        Fitting of Ellipses'
    (**) http://mathworld.wolfram.com/Ellipse.html
    (***) White, A. McHale, B. 'Faraday rotation data analysis with least-squares 
        elliptical fitting'
"""

class LSqEllipse:

    def fit(self, data):
        """Least Squares fitting algorithm 
        Theory taken from (*)
        Solving equation Sa=lCa. with a = |a b c d f g> and a1 = |a b c> 
            a2 = |d f g>
        Args
        ----
        data (list:list:float): list of two lists containing the x and y data of the
            ellipse. of the form [[x1, x2, ..., xi],[y1, y2, ..., yi]]
        Returns
        ------
        coef (list): list of the coefficients describing an ellipse
           [a,b,c,d,f,g] corresponding to ax**2+2bxy+cy**2+2dx+2fy+g
        """
        x, y = numpy.asarray(data, dtype=float)

        #Quadratic part of design matrix [eqn. 15] from (*)
        D1 = numpy.mat(numpy.vstack([x**2, x*y, y**2])).T
        #Linear part of design matrix [eqn. 16] from (*)
        D2 = numpy.mat(numpy.vstack([x, y, numpy.ones(len(x))])).T
        
        #forming scatter matrix [eqn. 17] from (*)
        S1 = D1.T*D1
        S2 = D1.T*D2
        S3 = D2.T*D2  
        
        #Constraint matrix [eqn. 18]
        C1 = numpy.mat('0. 0. 2.; 0. -1. 0.; 2. 0. 0.')

        #Reduced scatter matrix [eqn. 29]
        M=C1.I*(S1-S2*S3.I*S2.T)

        #M*|a b c >=l|a b c >. Find eigenvalues and eigenvectors from this equation [eqn. 28]
        eval, evec = numpy.linalg.eig(M) 

        # eigenvector must meet constraint 4ac - b^2 to be valid.
        cond = 4*numpy.multiply(evec[0, :], evec[2, :]) - numpy.power(evec[1, :], 2)
        a1 = evec[:, numpy.nonzero(cond.A > 0)[1]]
        
        #|d f g> = -S3^(-1)*S2^(T)*|a b c> [eqn. 24]
        a2 = -S3.I*S2.T*a1
        
        # eigenvectors |a b c d f g> 
        self.coef = numpy.vstack([a1, a2])
        self._save_parameters()
            
    def _save_parameters(self):
        """finds the important parameters of the fitted ellipse
        
        Theory taken form http://mathworld.wolfram
        Args
        -----
        coef (list): list of the coefficients describing an ellipse
           [a,b,c,d,f,g] corresponding to ax**2+2bxy+cy**2+2dx+2fy+g
        Returns
        _______
        center (List): of the form [x0, y0]
        width (float): major axis 
        height (float): minor axis
        phi (float): rotation of major axis form the x-axis in radians 
        """

        #eigenvectors are the coefficients of an ellipse in general form
        #a*x^2 + 2*b*x*y + c*y^2 + 2*d*x + 2*f*y + g = 0 [eqn. 15) from (**) or (***)
        a = self.coef[0,0]
        b = self.coef[1,0]/2.
        c = self.coef[2,0]
        d = self.coef[3,0]/2.
        f = self.coef[4,0]/2.
        g = self.coef[5,0]
        
        #finding center of ellipse [eqn.19 and 20] from (**)
        x0 = (c*d-b*f)/(b**2.-a*c)
        y0 = (a*f-b*d)/(b**2.-a*c)
        
        #Find the semi-axes lengths [eqn. 21 and 22] from (**)
        numerator = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
        denominator1 = (b*b-a*c)*( (c-a)*numpy.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
        denominator2 = (b*b-a*c)*( (a-c)*numpy.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
        width = numpy.sqrt(numerator/denominator1)
        height = numpy.sqrt(numerator/denominator2)

        # angle of counterclockwise rotation of major-axis of ellipse to x-axis [eqn. 23] from (**)
        # or [eqn. 26] from (***).
        phi = .5*numpy.arctan((2.*b)/(a-c))

        self._center = [x0, y0]
        self._width = width
        self._height = height
        self._phi = phi

    @property
    def center(self):
        return self._center

    @property
    def width(self):
        return self._width

    @property
    def height(self):
        return self._height

    @property
    def phi(self):
        """angle of counterclockwise rotation of major-axis of ellipse to x-axis 
        [eqn. 23] from (**)
        """
        return self._phi

    def parameters(self):
        return self.center, self.width, self.height, self.phi

#==================================================================
# Calculate eigenvalues and eigenvectors
#==================================================================

def calc_mjrmnrax(lons,lats):
  """
  Finds the major and minor axes of a set of points on earth

  Input: 
   1,2) Lists of longitude and latitude coordinates to fit the 
    major and minor axes too

  Output 

  Requires numpy 1.16.3 (conda install -c anaconda numpy; 
   https://pypi.org/project/numpy/)
  """

  # Import libraries
  import numpy as np

  # Calculate center of mass of shape
  center = [np.mean(lons),np.mean(lats)]

  # Calculate eigen-value/vector pairs for largest piece
  eigvals, eigvecs = np.linalg.eig(np.cov((lons[::-1],lats)))

  # Calculate coordinates of axes
  locseig0 = (center+(eigvals[0]*2)*eigvecs[0,:], 
              center-(eigvals[0]*2)*eigvecs[0,:])
  locseig1 = (center+(eigvals[1]*2)*eigvecs[1], 
              center-(eigvals[1]*2)*eigvecs[1])
  lonseig0 = [x for x,y in locseig0]
  lonseig1 = [x for x,y in locseig1]
  latseig0 = [y for x,y in locseig0]
  latseig1 = [y for x,y in locseig1]

  # Calculate lengths and angles of the axes
  lengths = np.zeros(2); angles = np.zeros(2)
  lengths[0],angles[0] = calc_distandangle(
   lonseig0[1],latseig0[1],lonseig0[0],latseig0[0])
  lengths[0],angles[0] = calc_distandangle(
   lonseig0[0],latseig0[0],lonseig0[1],latseig0[1])
  lengths[1],angles[1] = calc_distandangle(
   lonseig1[0],latseig1[0],lonseig1[1],latseig1[1])

  # Adjust angle since direction doesn't matter
  for i in range(len(angles)):
    if angles[i]>=180:
      angles[i] = angles[i]-180

  # Calculate which is the major and minor axes
  mjrind1 = np.argmax(lengths)
  if hasattr(mjrind1, "__len__"): mjrind=mjrind1[0]
  else: mjrind=mjrind1
  if mjrind==0: mnrind=1
  if mjrind==1: mnrind=0

  # Assign axes
  mjrax_len = lengths[mjrind]
  mnrax_len = lengths[mnrind]
  mjrax_ang = angles[mjrind]
  mnrax_ang = angles[mnrind]

  # Return
  return(center,mjrax_len,mjrax_ang,mnrax_len,mnrax_ang)

#==================================================================
# Time variables
#==================================================================


#  if nl.addnormtime=="True":
#    # Calculate normalized time variable
#    np.seterr(divide='ignore', invalid='ignore')
#    normalizedtime = [i for i in range(len(datakeys))]/\
#     np.float32(len(datakeys)-1)


      # Calculate local solar hour
      localsolarhour[c] = round(float(localsolartime[c][8:10])+ \
       float(localsolartime[c][10:12])/60,2)

