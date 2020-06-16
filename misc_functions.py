#==================================================================
# Miscellaneous functions
#==================================================================
#
# James Russell 2020
#
# * k_closest
#   - Finds the k closest values in a list to a given value
#
# * create_2d_dataframe
#   - Generates a dataframe from lists of coordinates and 
#      the corresponding data.
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
# Import libraries
#==================================================================

from numba import jit
import numpy as np
import pandas as pd

#==================================================================
# Find indices of k closest values in list
#==================================================================

def cartesian_distance(x1,x2,y1,y2):
  """
  Propagation 

  Input: 
   1) x locations
   2) y locations
   3) time 1 and 2

  Output:
   
  """

  # Calculate speed
  return(np.sqrt((x2-x1)**2+(y2-y1)**2))



#==================================================================
# Find indices of k closest values in list
#==================================================================

def cartesian_speed(x1,x2,y1,y2,t1,t2):
  """
  Propagation 

  Input: 
   1) x locations
   2) y locations
   3) time 1 and 2

  Output:
   
  """

  # Calculate speed
  return(cartesian_distance(x1,x2,y1,y2)/(t2-t1))



#==================================================================
# Calculate direction
#==================================================================

def cartesian_direction(x1,y1,x2,y2):
  """
  Calculates direction on earth between two sets of y, x 
   coordinates.

  Inputs:
  1,2) x1, y1
  3,4) x2, y2

  Output is the direction in degrees from north.

  Requires geopy 1.20.0 (conda install -c conda-forge geopy; 
   https://pypi.org/project/geopy/) and numpy.
  """

  # Case where both locations are the same
  if abs(y1-y2)<0.0001 and abs(x1-x2)<0.0001:
    direction = float('NaN')
  # Case where vector is directly east or west
  elif abs(y1-y2)<0.0001:
    if x1<x2: direction=90
    if x1>x2: direction=270
  # Case where vector is directly north or south
  elif abs(x1-x2)<0.0001:
    if y1<y2: direction=0
    if y1>y2: direction=180
  # All other cases
  else: direction = np.rad2deg(np.arctan((x1-x2)/(y1-y2)))
  if direction<0: direction = direction + 360

  # Return direction
  return(direction)



#==================================================================
# Find indices of k closest values in list
#==================================================================

def divzero(n,d):
  """
  Division where divide by zero equals 0

  Input: 
   1) Numerator
   2) Denominator

  Output:
   Division calculation or zero if d is zero
  """

  return(n/d if d else 0)



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
# Find indices of k closest values in list
#==================================================================

def k_closest_ma(lst,value,k):
  """
  Find the indices of the k closest values in list to a given value

  Input: 
   1) A masked list of values
   2) The value you wish to find the closest value to in lst
   3) The number of closest values you wish to find.

  Output:
   A list of the indices in lst corresponding to the k closest
    values to value

  Requires numpy 1.16.3 (conda install -c anaconda numpy; 
   https://pypi.org/project/numpy/)
  """

  # Find the absolute differences between the value and all values
  #  in the list
  valuea = np.array([value]*len(lst))
  diff = np.absolute(np.subtract(lst,valuea))

  # Sort these differences
  diffs = sorted(diff.filled(np.nan))
  diffs = [x for x in diffs if np.isnan(x)==False]

  # Set the k smallest values in diff equal to zero
  for d in diffs[0:k]:
    diff = np.ma.where(diff==d,0,diff)
 
  # Return indices where diff equals zero
  return([i for i, x in enumerate(diff) if x==0])



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
# Write netcdf variable
#==================================================================

def write_var(varname,long_name,description,dimname,dtype,units,
              fileout,datain,f,fv):
  """
  Generic script for defining/opening a variable in a 
   netcdf file and then writing the associated data to the 
   variable.
  """

  try:
    datafile = fileout.createVariable(varname, dtype, (dimname),
     fill_value=fv)
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



