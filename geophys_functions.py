#==================================================================
# Geophysical functions
#==================================================================
# 
# Functions to calculate various geophysical things
# James Russell 2020
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
# * calc_direction
#   - Calculates a direction given points defining two vectors
#
# * calc_distdir
#   - Combines the above two functions
#
# * calc_propagation
#   - Calculate propagation speed and direction between 
#      two points on earth
#
# * load_land
#   - Loads a land shape file
# 
# * is_land
#   - Takes a latitude and longitude location, reads in a 
#      land area shape file, it checks location against 
#      shape file to return True if location is over land.
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
