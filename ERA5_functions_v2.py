def driver_addvars(fn):
  """
  Driver to calculate various variables for a single PF.
  fn corresponds to the line number in filenames_av.txt.
  """

#==================================================================
# Begin script
#==================================================================

  # Import libraries
  import TIPS_functions as fns
  import ERA5_functions_v2 as E5fns
  from netCDF4 import Dataset
  import numpy as np
  import datetime as dt

  # Read in namelist variables
  nl = Dataset("namelist_av.nc","r")

  # Read in filename for PF
  for i, row in enumerate(open("filenames_av.txt")):
    if i==fn:
      f = row[:-1]

  print("Working with file "+str(f))

  # Open file and assign data
  fd = Dataset(f)
  datalat  = fd.groups["lats"].groups["data"]
  datalon  = fd.groups["lons"].groups["data"]
  dataclat = fd.variables["centrallat"][:]
  dataclon = fd.variables["centrallon"][:]
  datadtim = fd.variables["datetime"][:]
