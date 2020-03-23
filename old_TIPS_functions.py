#==================================================================
# Functions used to add variables to IMERG PFs
# Author: James Russell 2019
#
# Tested with:
# * numpy 1.16.3
# * scipy 1.2.1
# * matplotlib 3.0.3
# * pandas 0.24.2
# * cartopy 0.17.0
# * shapely 1.6.4
# * fiona 1.8.6
# * area 1.1.1
# * pyproj 1.9.6
# * geopy 1.20.0
#
# Description of functions:
#
# * calc_area
#   - Calculates area on earth from a list of pixels 
#      locations
#
# * calc_area_and_volrainrate
#   - Calculates these variables on earth from a list of 
#      pixels locations and their rain rates
#
# * calc_distance
#   - Calculates distance between two points on earth 
#
# * calc_distance_and_angle
#   - Calculates distance between two points on earth 
#      and the angle the vector those points make from 
#      north
#
# * calc_propagation
#   - Calculate propagation speed and direction between 
#      two points on earth
#
# * interp_TC
#   - Reads IBtracs database and interpolates position of 
#      all TCs within 3 hours of input time, to the input 
#      time. Also outputs a list of largest possible TC 
#      radii.
#
# * calc_if_TC
#   - Takes output from interp_TC and a latitude, longitude 
#      location and calculates whether a TC is close to 
#      that location. Outputs information on the TC that is 
#      closest if it is within largest possible TC radius.
#
# * load_land
#   - Loads a land shape file
# 
# * is_land
#   - Takes a latitude and longitude location, reads in a 
#      land area shape file, it checks location against 
#      shape file to return True if location is over land.
#
# * calc_local_solar_time
#   - Takes a date and time, and a longitude, and adds an 
#      offset factor to give a local solar time (i.e. for a
#      diurnal cycle). Not the actual local time.
#
# * create_2d_dataframe
#   - Generates a dataframe from lists of coordinates and 
#      the corresponding data.
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
# * calc_mjrmnrax
#   - Finds the major and minor axis of a set of points on earth
#
# * k_closest
#   - Finds the k closest values in a list to a given value
#
# * time_since
#   - From a string indicating the time and date and a units 
#     since reference time and date, calculate a time number
#
# * write_var
#   - Defines/opens a variable and writes the data to a 
#      netcdf file.
#
# * write_group
#   - Defines/opens a group and writes data to a data group
#      within that netcdf group as attribute and value 
#      pairs. This is a means for writing a python 
#      dictionary to a netcdf file.
#
#==================================================================

#==================================================================
# Calculate area 
#==================================================================

def calc_area(lons,lats,dx,dy):
  """
  Calculates the area on earth given a set of pixels.
  
  Inputs:
  1,2) Lists of latitude and longitude coordinates of pixel 
   centers
  3,4) Scalars of the pixel width (dx) and height (dy). 
  All inputs should be in units of degrees longitude or latitude.
  
  Output is an area in units of m**2.

  Requires area 1.1.1 (https://pypi.org/project/area/).
   Note: At time of writing this was not available 
   through conda and had to be installed manually.
  """

  # Import libraries
  from area import area
  import numpy

  # Loop over all pixels
  for l in range(0,len(lats)):
    # Create geojson format variable for each pixel. 
    geom = {'type': 'Polygon',
            'coordinates': [[[lons[l]-(dx/2.), lats[l]-(dy/2.)], 
                             [lons[l]+(dx/2.), lats[l]-(dy/2.)],
                             [lons[l]+(dx/2.), lats[l]+(dy/2.)],
                             [lons[l]-(dx/2.), lats[l]+(dy/2.)],
                             [lons[l]-(dx/2.), lats[l]-(dy/2.)]]]}

    # Calculate area for each individual pixel
    proj_area = area(geom)

    # Sum area of all pixels
    if l==0:
      area1 = proj_area
    else:
      area1 = area1 + proj_area

  # Return area
  return(area1)



#==================================================================
# Calculate area and volumetric rain rate
#==================================================================

def calc_area_and_volrainrate(lons,lats,rain,dx,dy):
  """
  Calculates the area on earth and the volumetric rain 
   rate of that area given a set of pixels.
  
  Inputs:
   1,2) Lists of latitude and longitude coordinates of pixel 
    centers.
   3) List of corresponding rainfalls for each pixel
   4,5) Scalars of the pixel width (dx) and height (dy). 
   All inputs should be in units of degrees longitude or 
    latitude, with rainfall in mm/hr.
  
  Output is an area in units of m**2 and a volumetric rain 
   rate in units of mm/hr m**2.

  Requires area 1.1.1 (https://pypi.org/project/area/).
   Note: At time of writing this was not available 
   through conda and had to be installed manually.
  """

  # Import libraries
  from area import area

  # Loop over all pixels
  for l in range(0,len(lats)):
    # Create geojson format variable for each pixel. 
    geom = {'type': 'Polygon',
            'coordinates': [[[lons[l]-(dx/2.), lats[l]-(dy/2.)], 
                             [lons[l]+(dx/2.), lats[l]-(dy/2.)],
                             [lons[l]+(dx/2.), lats[l]+(dy/2.)],
                             [lons[l]-(dx/2.), lats[l]+(dy/2.)],
                             [lons[l]-(dx/2.), lats[l]-(dy/2.)]]]}

    # Calculate area for each individual pixel
    proj_area = area(geom)

    # Calculate volumetric rain rate per pixel
    volrainratel = proj_area * rain[l]

    # Add area and volumetric rain rates of all pixels
    if l==0:
      area1 = proj_area
      volrainrate = volrainratel
    else:
      area1 = area1 + proj_area
      volrainrate = volrainrate + volrainratel

  # Return area and volumetric rain rate
  return(area1,volrainrate)



#==================================================================
# Calculate distance between two points on earth
#==================================================================

def calc_distance(lon1,lat1,lon2,lat2):
  """
  Calculates distance on earth between two sets of 
   latitude, longitude coordinates.

  Input:
   1,2) lon1, lat1
   3,4) lon2, lat2

  Output is the distance in units of m.

  Requires geopy 1.20.0 (conda install -c conda-forge geopy; 
   https://pypi.org/project/geopy/).
  """

  # Import libraries
  from geopy.distance import geodesic

  # Calculate distance on earth
  dist = geodesic((lat1,lon1),(lat2,lon2)).m

  # Return distance
  return(dist)



#==================================================================
# Calculate direction on earth
#==================================================================

def calc_direction(lon1,lat1,lon2,lat2):
  """
  Calculates direction on earth between two sets of 
   latitude, longitude coordinates.

  Inputs:
  1,2) lat1, lon1
  3,4) lat2, lon2

  Output is the direction in degrees from north.

  Requires geopy 1.20.0 (conda install -c conda-forge geopy; 
   https://pypi.org/project/geopy/) and numpy.
  """

  # Import libraries
  from geopy.distance import geodesic
  import numpy as np

  # Case where both locations are the same
  if abs(lat1-lat2)<0.0001 and abs(lon1-lon2)<0.0001:
    direction = float('NaN')
  # Case where vector is directly east or west
  elif abs(lat1-lat2)<0.0001:
    if lon1<lon2: direction=90
    if lon1>lon2: direction=270
  # Case where vector is directly north or south
  elif abs(lon1-lon2)<0.0001:
    if lat1<lat2: direction=0
    if lat1>lat2: direction=180
  # All other cases
  else:
    direction = np.rad2deg(np.arctan(
     calc_distance(lon1,lat1,lon2,lat2)/
     calc_distance(lon1,lat1,lon2,lat2)))

  # Adjust for quadrant
  if lat1>lat2 and lon1<lon2: direction = 180-direction
  elif lat1>lat2 and lon1>lon2: direction = 180+direction
  elif lat1<lat2 and lon1>lon2: direction = 360-direction

  # Return direction
  return(direction)



#==================================================================
# Calculate distance and angle of a vector on earth
#==================================================================

def calc_distdir(lon1,lat1,lon2,lat2):
  """
  Calculates distance on earth between two sets of latitude, 
   longitude coordinates, and direction that vector makes from 
   north.

  Inputs:
  1,2) lat1, lon1
  3,4) lat2, lon2

  Output is the distance in units of m and direction 
   the vector makes from north.

  Requires geopy 1.20.0 (conda install -c conda-forge geopy; 
   https://pypi.org/project/geopy/).
  """

  # Calculate distance on earth
  distance = calc_distance(lon1,lat1,lon2,lat2)

  # calculate direction on earth
  direction = calc_direction(lon1,lat1,lon2,lat2)

  # Return distance and direction
  return(distance,direction)



#==================================================================
# Calculate propagation
#==================================================================

def calc_propagation(date1,lon1,lat1,date2,lon2,lat2):
  """
  Calculates propagation speed and direction on earth 
   given two dates and times with corresponding latitudes 
   and longitudes.

  Inputs:
   1,2,3) date and time 1, lat1, lon1
   4,5,6) date and time 2, lat2, lon2. 
   Dates and times should be in string format YYYYMMDDhhmmss. 
   All others should be floats or integers.

  Output is the speed in units of m/s and angle the 
   propagation vector makes from north.

  Requires pyproj 1.9.6 (conda install -c conda-forge pyproj; 
   https://pypi.org/project/pyproj/) and geopy 1.20.0
   (conda install -c conda-forge geopy; 
   https://pypi.org/project/geopy/).
  """

  # Import libraries
  import datetime as dt

  # Calculate distance and propagation direction
  distance,propdir = calc_distdir(lon1,lat1,lon2,lat2)

  # Make time objects
  dtob1 = dt.datetime(int(date1[0:4]),int(date1[4:6]),
                      int(date1[6:8]),int(date1[8:10]),
                      int(date1[10:12]),int(date1[12:14]))
  dtob2 = dt.datetime(int(date2[0:4]),int(date2[4:6]),
                      int(date2[6:8]),int(date2[8:10]),
                      int(date2[10:12]),int(date2[12:14]))

  # Calculate propagation speed
  propspd = distance/((dtob2-dtob1).seconds)

  # Return speed and direction
  return(propspd,propdir)



#==================================================================
# Interpolate TC information to a time
#==================================================================

def interp_TC(dtim,fTC):
  """
  Interpolates information on a TC to a given a time. 

  Inputs:
  1) A date and time in YYYYMMDDhhmm string format to interpolate
  2) File handle (fTC) for IBTrACS data.

  Output is  a dictionary of information about the TC after
   interpolating that information to the same time as that 
   input. In that dictionary, it returns lists of TC radii,
   lat and lon position of the TC center, and the 
   corresponding indices of the data in the ibtracs data. If
   TC is not within 1.5 hours of the given time, returns 
   dictionary with one key and value pair TC:False. If within 
   1.5 hours this TC key is True and all other information is 
   provided."

  Requires geopy 1.20.0 (conda install -c conda-forge geopy; 
   https://pypi.org/project/geopy/).
  """

  # Import libraries 
  import datetime
  import time
  from geopy.distance import geodesic 
  import numpy as np

  # Convert current time to TC time units
  #  *Hard coded for IBTrACS default unit of days
  timeunits = fTC.variables["time"].units
  d1   = datetime.datetime.strptime(timeunits[11:21], 
                                    "%Y-%m-%d")
  d2   = datetime.datetime.strptime(str(dtim),
                                    "%Y%m%d%H%M")
  d1ts = time.mktime(d1.timetuple())
  d2ts = time.mktime(d2.timetuple())
  TCtimenow = (d2ts-d1ts)/(3600.*24)

  # Check for TC times within just over 1.5 hours
  #  i.e. since TCs observations are typically 3 hours apart
  TCtimelow = TCtimenow-0.08
  TCtimehig = TCtimenow+0.08
  indTC = np.asarray(np.where(
    (fTC.variables["time"][:]>TCtimelow) 
    & (fTC.variables["time"][:]<TCtimehig))).T.tolist()

  if not indTC:
    TCinfo = {"TC":False}
    return(TCinfo)

  TCinfo = {"TC":True}

  # Preallocate some output arrays
  TCradmax = [0.]*len(indTC)
  TClatnow = [0.]*len(indTC)
  TClonnow = [0.]*len(indTC)

  # Loop over all TCs with a close time
  for iTC in range(len(indTC)):
    # Find TC times and locations and interpolate them to
    #  the time of the PF

    # Instance where time of PF and TC ob are same
    if abs(TCtimenow-
     fTC.variables["time"][indTC[iTC][0],indTC[iTC][1]])<0.0001:
      TClatnow[iTC] = fTC.variables["lat"][indTC[iTC][0],
                                           indTC[iTC][1]]
      TClonnow[iTC] = fTC.variables["lon"][indTC[iTC][0],
                                           indTC[iTC][1]]

    # Instance where time of PF is earlier than first TC ob 
    elif ((TCtimenow<fTC.variables["time"][indTC[iTC][0],
              indTC[iTC][1]]) & (indTC[iTC][1]==0)):
      TClatnow[iTC] = fTC.variables["lat"][indTC[iTC][0],
                                           indTC[iTC][1]]
      TClonnow[iTC] = fTC.variables["lon"][indTC[iTC][0],
                                           indTC[iTC][1]]

    # Instance where time of PF is earlier than TC ob 
    elif TCtimenow<fTC.variables["time"][indTC[iTC][0],
                                         indTC[iTC][1]]:
      TClatnow[iTC] = np.interp(TCtimenow,
        fTC.variables["time"][indTC[iTC][0],
        indTC[iTC][1]-1:indTC[iTC][1]+1],
        fTC.variables["lat"][indTC[iTC][0],
        indTC[iTC][1]-1:indTC[iTC][1]+1])
      TClonnow[iTC] = np.interp(TCtimenow,
        fTC.variables["time"][indTC[iTC][0],
        indTC[iTC][1]-1:indTC[iTC][1]+1],
        fTC.variables["lon"][indTC[iTC][0],
        indTC[iTC][1]-1:indTC[iTC][1]+1])

    # Instance where time of PF is later than first TC ob 
    elif TCtimenow>fTC.variables["time"][indTC[iTC][0],
                                         indTC[iTC][1]]:
      TCtim = fTC.variables["time"][indTC[iTC][0],
        indTC[iTC][1]:indTC[iTC][1]+2]
      TCtim = TCtim[~TCtim.mask]
      if len(TCtim)==2:
        TClat = fTC.variables["lat"][indTC[iTC][0],
          indTC[iTC][1]:indTC[iTC][1]+2]
        TClon = fTC.variables["lon"][indTC[iTC][0],
          indTC[iTC][1]:indTC[iTC][1]+2]
        TClatnow[iTC] = np.interp(TCtimenow,TCtim,TClat)
        TClonnow[iTC] = np.interp(TCtimenow,TCtim,TClon)
      else:
        TClatnow[iTC] = fTC.variables["lat"][indTC[iTC][0],
                                             indTC[iTC][1]]
        TClonnow[iTC] = fTC.variables["lon"][indTC[iTC][0],
                                             indTC[iTC][1]]

    # Read in the widest TC radius possible
    TCrad = np.ma.empty(6)
    TCrad[0] = fTC.variables["td9635_roci"][indTC[iTC][0],
                                            indTC[iTC][1]]
    TCrad[1] = fTC.variables["bom_roci"][indTC[iTC][0],
                                         indTC[iTC][1]]
    TCrad[2] = max(fTC.variables["bom_r34"][indTC[iTC][0],
                                            indTC[iTC][1],:])
    TCrad[3] = max(fTC.variables["reunion_r34"][indTC[iTC][0],
                                                indTC[iTC][1],:])
    TCrad[4] = fTC.variables["tokyo_r30_long"][indTC[iTC][0],
                                               indTC[iTC][1]]
    TCrad[5] = max(fTC.variables["usa_r34"][indTC[iTC][0],
                                            indTC[iTC][1],:])

    if TCrad.count()==0:
    # If there is no information on a minimum radius select a 
    #  typical average radius of 0 m/s winds(source: Chavas et
    #  al. 2016: Observed Tropical Cyclone Size Revisited. JCLI)
      TCradmax[iTC] = 600.
    else:
    # If there is information, select the largest and 
    #  convert for nautical miles to km and double
    #  to ensure all of TC is accounted for (i.e. lighter 
    #  than 30 kt winds. Source above shows this is likely 
    #  underestimate)
      TCradmax[iTC] = (np.ma.MaskedArray.max(TCrad)*1.852)*2

  TCinfo["TCrad"] = TCradmax
  TCinfo["TClat"] = TClatnow
  TCinfo["TClon"] = TClonnow
  TCinfo["indTC"] = indTC

  return(TCinfo)



#==================================================================
# Calculate TC information
#==================================================================

def calc_if_TC(lon,lat,TCinfo,fTC):
  """
  Calculates proximity to TC. 

  Inputs:
   1,2) lat and lon locations for position to check
   3) A dictionary with:
    TClat and TClon - position of TC center
    TCrad - list of radiuses 
   4) Open file handle of the IBTrACS dataset

  Outputs are the distance between the lat, lon location 
   and the active TC centers at that time, a boolean 
   indicating if it is within a radius of the TC, the name 
   of the TC it is within radius of, and the radius used to
   estimate if location is within the TC.

  Requires geopy 1.20.0 (conda install -c conda-forge geopy; 
   https://pypi.org/project/geopy/).
  """

  # Import libraries
  from geopy.distance import geodesic 
  import numpy as np
  import netCDF4 as nc

  # Identify closest TC
  dist_loc_cTC = [0.]*len(TCinfo["TClat"])
  for i in range(len(TCinfo["TClat"])):
    # Check if location is within distance to TC center
    dist_loc_cTC[i] = geodesic((lat,lon),(TCinfo["TClat"][i],
                                          TCinfo["TClon"][i])).km
  ind_closest_TC = dist_loc_cTC.index(min(dist_loc_cTC))
  dist_cTC = dist_loc_cTC[ind_closest_TC]
  TCradius = TCinfo["TCrad"][ind_closest_TC]

  # Check if within radius
  in_TC = 0
  TCname = ""
  if dist_cTC<TCradius:
    in_TC    = 1
    TCname   = str(nc.chartostring(
                  fTC.variables["name"][TCinfo["indTC"][
                  ind_closest_TC][0]]))

  return(dist_cTC,in_TC,TCname,TCradius)



#==================================================================
# Function to load land shape file
#==================================================================

def load_land(res='50m'):
  """
  Load land shape file

  Inputs:
   1) Optional resolution arguement, string format. Default = 50m.

  Outputs land shapoe file.

  Requires fiona 1.8.6 (conda install -c conda-forge fiona;
   https://pypi.org/project/Fiona/), cartopy 0.17.0
   (conda install -c conda-forge cartopy; 
   https://pypi.org/project/Cartopy/), and shapely 1.6.4
   (conda install -c conda-forge shapely; 
   https://pypi.org/project/Shapely/)
   
  """

  # Import libraries
  import fiona
  import cartopy.io.shapereader as shpreader
  import shapely.geometry as sgeom
  from shapely.prepared import prep
  
  # Pull in land shape files
  geoms = fiona.open(shpreader.natural_earth(resolution=res,
                           category='physical', name='land'))
  land_geom = sgeom.MultiPolygon([sgeom.shape(geom['geometry'])
                                for geom in geoms])
  land = prep(land_geom)
  
  # Check if point is over land and return information
  return(land)



#==================================================================
# Function to check if a point is over land
#==================================================================

def is_land(lon, lat, res='50m'):
  """
  Check if a location is over land. 

  Inputs:
   1,2) Longitude (x) and latitude (y) of position to check.
   3) Optional resolution arguement, string format. Default = 50m.

  Outputs True if over land. Outputs False if not.

  Requires fiona 1.8.6 (conda install -c conda-forge fiona;
   https://pypi.org/project/Fiona/), cartopy 0.17.0
   (conda install -c conda-forge cartopy; 
   https://pypi.org/project/Cartopy/), and shapely 1.6.4
   (conda install -c conda-forge shapely; 
   https://pypi.org/project/Shapely/)
   
  """

  # Import libraries
  import fiona
  import cartopy.io.shapereader as shpreader
  import shapely.geometry as sgeom
  from shapely.prepared import prep
  
  # Pull in land shape files
  geoms = fiona.open(shpreader.natural_earth(resolution=res,
                           category='physical', name='land'))
  land_geom = sgeom.MultiPolygon([sgeom.shape(geom['geometry'])
                                for geom in geoms])
  land = prep(land_geom)
  
  # Check if point is over land and return information
  return(land.contains(sgeom.Point(lon, lat)))



#==================================================================
# Calculate local solar time based on longitude
#==================================================================

def calc_local_solar_time(yr,mo,dy,hr,mn,sc,lon):
  """
  Calculates the local solar time given a date and time and 
   location. Calculated as the UTC time plus an offset based 
   on the longitude. The offset is calculated by multiplying 
   the longitude by 24/360. 

  Inputs:
   1) The date and time (dt). dt is a string indicating UTC
    time in format YYYYMMDDhhmmss.
   2) The longitude position (lon). lon is the longitude 
    position as a float.

  Output is a string in YYYYMMDDhhmmss format indicating the
   local solar time.

  Note: This is not the actual local time. This should 
   typically only be used to calculate times for the 
   diurnal cycle.
  """

  # Import libraries
  import datetime

  # Make a datetime object for the current date and time
  timenow = datetime.datetime(int(yr),int(mo),int(dy),
                              int(hr),int(mn),int(sc))
 
  # Calculate time with offset added
  newtime = timenow + datetime.timedelta(hours=lon*(24./360.))
  
  # Return time in same format as input
  return(str(newtime.year).zfill(4)+\
         str(newtime.month).zfill(2)+\
         str(newtime.day).zfill(2)+\
         str(newtime.hour).zfill(2)+\
         str(newtime.minute).zfill(2)+\
         str(newtime.second).zfill(2))



#==================================================================
# Create 2D dataframe given a list of coordinates and data
#==================================================================

def create_2d_dataframe(x,y,dx,dy,data):
  """
  Creates a dataframe given a set of pixel locations, pixel sizes,
   and the data in each pixel. Fills any coordinates that are in 
   the new dataframe but were not in the original lists with zeros.

  Inputs:
   1,2) Lists of pixel locations (indexes) x and y
   3,4) Scalars of pixel width (dx) and pixel height (dy)
   5) A corresponding list of the data values.

  Output is a 2D dataframe with x and y the indexes, and data as 
   the values

  Requires numpy 1.16.3 and pandas 0.24.2 (conda install -c 
   conda-forge pandas; https://pypi.org/project/pandas/)
  """

  # Import libraries
  import numpy as np
  import pandas as pd

  # Create coordinates
  ny = [round(i,2) for i in np.arange(min(y)-2*dy,max(y)+2*dy,dy)]
  nx = [round(i,2) for i in np.arange(min(x)-2*dx,max(x)+2*dx,dx)]

  # Generate dataframe
  df = pd.DataFrame(np.zeros([len(ny),len(nx)],float),
            [str(i) for i in ny],[str(i) for i in nx])

  # Populate dataframe
  for i in range(len(y)): df.loc[str(y[i]),str(x[i])] = data[i]

  # Return dataframe
  return(df)



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

  # Import libraries
  import scipy.ndimage as ndimg

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

  # Import libraries
  import numpy as np

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

  # Import libraries
  from matplotlib.path import Path

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

  # Import libraries
  import numpy as np

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
  axdir = [calc_direction(center[0],center[1],
           center[0]+a[0],center[1]+a[1]) for a in U]
  axdir = [d-180 if d>180 else d for d in axdir]

  # Calculate length in degrees of axes
  axlen = [2*np.sqrt(2/N)*l for l in S]

  # Calculate geopgraphic distance of axes
  axdist =  [calc_distance(center[0] - ((axlen[i]/2)* \
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

  if fit:
    # Also return fit
    return(center,axdir,axlen,fitxy)
  else:
    # Return center, and direction and length of major and minor
    #  axes
    return(center,axdir,axdist)

#==================================================================
# Find indices of k closest values in list
#==================================================================

def k_closest(lst,value,k):
  """
  Find the indices of the k closest values in list to a given value

  Input: 
   1) A list of values
   2) The value you wish to find the closest value to in lst
   3) The number of closest values you wish to find.

  Output:
   A list of the indices in lst corresponding to the k closest
    values to value

  Requires numpy 1.16.3 (conda install -c anaconda numpy; 
   https://pypi.org/project/numpy/)
  """

  import numpy as np

  # Find the absolute differences between the value and all values
  #  in the list
  diff  = [abs(l-value) for l in lst]

  # Sort these differences
  diffs = sorted(diff)

  # Set the k smallest values in diff equal to zero
  for d in diffs[0:k]:
    diff = np.where(diff==d,0,diff)

  # Return indices where diff equals zero
  return([i for i, x in enumerate(diff) if x==0])



#==================================================================
# Calculate time since given units
#==================================================================

def time_since(dateandtime,reftimeandunits):
  """
  Calculate the date and time in units since a date

  Input: 
   1) The date and time in string format: "yyyy-mm-dd HH:MM:SS"
   2) The units and reference time to convert to in string format:
       <units> since yyyy-mm-dd HH:MM:SS
       <units> can be any of the below:
       * seconds or secs
       * minutes or mins
       * hours or hrs
       * days or dys
       * weeks or wks

  Output:
   A float indicating the current time in the units specified

  Requires datetime
  """

  # Import libraries
  import datetime as dt

  # Split units string and make date and time object
  splitrtu = reftimeandunits.split(" ")
  unit = splitrtu[0]
  refdtob = dt.datetime(int(splitrtu[2][0:4]),
                        int(splitrtu[2][5:7]),
                        int(splitrtu[2][8:10]),
                        int(splitrtu[3][0:2]),
                        int(splitrtu[3][3:5]),
                        int(splitrtu[3][6:8]))

  # Make current date and time object
  curdtob = dt.datetime(int(dateandtime[0:4]),
                        int(dateandtime[5:7]),
                        int(dateandtime[8:10]),
                        int(dateandtime[11:13]),
                        int(dateandtime[14:16]),
                        int(dateandtime[17:19]))

  # Return current time given units
  if unit=="seconds" or unit=="secs":
    return((curdtob-refdtob).total_seconds())
  elif unit=="minutes" or unit=="mins":
    return((curdtob-refdtob).total_seconds()/60)
  elif unit=="hours" or unit=="hrs":
    return((curdtob-refdtob).total_seconds()/3600)
  elif unit=="days" or unit=="dys":
    return((curdtob-refdtob).total_seconds()/86400)
  elif unit=="weeks" or unit=="wks":
    return((curdtob-refdtob).total_seconds()/604800)



#==================================================================
# Calculate time since given units
#==================================================================

def time_since_inv(datenumber,reftimeandunits):
  """
  Convert a date number (e.g. time units since a date) to a date 
   in string format

  Input: 
   1) The date number
   2) The units and reference time of the date number:
       <units> since yyyy-mm-dd HH:MM:SS
       <units> can be any of the below:
       * seconds or secs
       * minutes or mins
       * hours or hrs
       * days or dys
       * weeks or wks

  Output:
   A float indicating the current time in the units specified

  Requires datetime
  """

  # Import libraries
  import datetime as dt

  # Split units string and make date and time object
  splitrtu = reftimeandunits.split(" ")
  unit = splitrtu[0]
  refdtob = dt.datetime(int(splitrtu[2][0:4]),
                        int(splitrtu[2][5:7]),
                        int(splitrtu[2][8:10]),
                        int(splitrtu[3][0:2]),
                        int(splitrtu[3][3:5]),
                        int(splitrtu[3][6:8]))

  # Get date number in seconds
  if unit=="minutes" or unit=="mins":
    datenumber = datenumber*60
  elif unit=="hours" or unit=="hrs":
    datenumber = datenumber*3600
  elif unit=="days" or unit=="dys":
    datenumber = datenumber*86400
  elif unit=="weeks" or unit=="wks":
    datenumber = datenumber*604800


  # Using timedelta return a string of the current 
  return(str(refdtob+dt.timedelta(seconds=datenumber)))
  
  

#==================================================================
# Calculate local time
#==================================================================

def utc_to_local(utcdatetime,timezone):
 
  """
  Calculate the date and time in units since a date

  Input: 
   1) The UTC date and time in string format: "yyyy-mm-dd HH:MM:SS"
   2) A time zone for pytz such as America/Denver or US/Mountain. 
   A list of possible options can be found here:
   https://gist.github.com/heyalexej/8bf688fd67d7199be4a1682b3eec7568

  Output:
   The local date and time in string format: "yyyy-mm-dd HH:MM:SS"

  Requires datetime
  """

  # Import libraries
  from datetime import datetime
  from dateutil import tz

  # Get time zones 
  from_zone = tz.gettz('UTC')
  to_zone = tz.gettz(timezone)

  # Strip time and assign time zone
  utc = datetime.strptime(utcdatetime,
         '%Y-%m-%d %H:%M:%S').replace(tzinfo=from_zone)
  
  # Convert to local time
  lt = utc.astimezone(to_zone)

  # Return local time string
  return(
   str(lt.year).zfill(4)+"-"+str(lt.month).zfill(2)+"-"+\
   str(lt.day).zfill(2)+" "+str(lt.hour).zfill(2)+":"+\
   str(lt.minute).zfill(2)+":"+str(lt.second).zfill(2))



#==================================================================
# Write netcdf variable
#==================================================================

def write_var(varname,long_name,description,dimname,dtype,units,
              fileout,datain,f):
  """
  Generic script for defining/opening a variable in a 
   netcdf file and then writing the associated data to the 
   variable.
  """

  try:
    datafile = fileout.createVariable(varname, dtype, (dimname))
  except:
    #print(varname + " already defined")
    datafile = fileout.variables[varname]
  #print("Writing "+varname+" to "+f)
  datafile.long_name   = long_name
  datafile.description = description
  datafile.units = units
  datafile[:] = datain 



#==================================================================
# Write netcdf groups as attribute and value pairs 
#==================================================================

def write_group(groupname,long_name,description,units,
                format1,fileout,datain,f):
   """
   Generic script for defining/opening a group in a netcdf
    file and then writing a python dictionary to that group
    as attribute:value pairs. Note: value can be a string, 
    integer, float, list, or anything else that can be 
    taken as an attribute in a netcdf file.
   """

   try:
     datafile  = fileout.createGroup(groupname)
   except:
     #print(groupname+" already defined")
     datafile = fileout.group[groupname]
   try:
     datadatafile = datafile.createGroup('data')
   except:
     #print("data"+groupname+" already defined")
     datadatafile = fileout.group[groupname].groups["data"]

   #print("Writing "+groupname+" to "+f)
   datafile.long_name   = long_name
   datafile.description = description
   datafile.units = units
   datafile.format = format1
   for k,v in datain.items():
     setattr(datadatafile, k,  v)



#==================================================================
# End
#==================================================================
