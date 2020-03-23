#==================================================================
# Time functions
#==================================================================
# 
# Functions to calculate various time related variables.
# James Russell 2020
#
# * time_since
#   - From a string indicating the time and date and a units 
#     since reference time and date, calculate a time number.
#
# * time_since_inv
#   - Converts a time number and its units into a date and time 
#      string.
#
# * calc_local_time
#   - Takes a date and time, and a location, and calculates
#      the local time.
#
# * calc_local_solar_time
#   - Takes a date and time, and a longitude, and adds an 
#      offset factor to give a local solar time (i.e. for a
#      diurnal cycle). Not the actual local time.
#
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
