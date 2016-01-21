#!/usr/bin/env python

from netCDF4 import Dataset, num2date
from datetime import datetime, timedelta
import numpy as np

def get_gefs_ensemble(variables, bounds=None, start=datetime.utcnow()-timedelta(hours=12),
                      end=datetime.utcnow()+timedelta(hours=48), writeout=True):

    """
    Retrieves the latest ("best") ensemble forecast valid at a single point from the Unidata THREDDS server using
    the Unidata siphon library.
    
    Requires:
    point -> A tuple of (lat, lon) of the point we are trying to retrieve
    variables -> A list of variables we want to retrieve.  Check this page for a full list:
            http://thredds.ucar.edu/thredds/metadata/grib/NCEP/GEFS/Global_1p0deg_Ensemble/members/Best?metadata=variableMap
    start -> A datetime object of the earliest time to look for an ensemble initialization,
            default is current time minus 12 hours
    end -> The last time for which we want ensemble forecast output.  Default is current time plus 48 hours.
    
    Returns:
    A dictionary with one item being the list of valid times in the data ('times') and the rest of the items
    being numpy arrays of nTimes x nEnsmems for each variable requested
        
    """
    # Import the Siphon utilities 
    from siphon.catalog import TDSCatalog
    from siphon.ncss import NCSS
    
    # In Siphon, we connect to a thredds catalog.  Here's the address for the GEFS
    catalog = 'http://thredds.ucar.edu/thredds/catalog/grib/NCEP/GEFS/Global_1p0deg_Ensemble/members/catalog.xml' 
    best_model = TDSCatalog(catalog)
    
    # We select a specific dataset in this catalog, in this case the "best" (most recent) ensemble run 
    best_ds = list(best_model.datasets.values())[2]
    ncss = NCSS(best_ds.access_urls['NetcdfSubset'])

    
    
    # Here we format our subsetting query.  We specify the exact point we want,
    # the time range, and the variables we are requesting.  We're also going
    # to retrieve the data in a netcdf-like format
    query = ncss.query()
    #print dir(query)
    #return
    if bounds is not None:
        # Assume bounds are West, East, South, North 
        query.lonlat_box(*bounds)

    #query.lonlat_point(point[1], point[0])
    query.time_range(start, end)
    query.variables(*variables)
    query.accept('netcdf')

    # Actually get the data
    data = ncss.get_data(query)
    print data


    # If we are writing this out...
    if writeout:
        ntime, nmem, dum, nlat, nlon = data.variables[variables[0]].shape
        # To get valid times, need to add the deltas to time_coverage_start
        #init_time = datetime.strptime(str(data.time_coverage_start), '%Y-%m-%dT%H:%M:%SZ')
        #valid_times = [init_time + timedelta(hours=int(x)) for x in list(data.variables['time2'][:])]
        time = data.variables['time2']
        valid_times = num2date(time[:], time.units)
        valid_timestamps = [(x - datetime(1970,1,1)).total_seconds() for x in list(valid_times)]

        with Dataset('gefs_forecast.nc','w',format="NETCDF3_CLASSIC") as dset:
          dset.createDimension('time', None)
          dset.createDimension('ens', nmem)
          dset.createDimension('lat', nlat)
          dset.createDimension('lon', nlon)
          dset.createVariable('time','i4', ('time',))
          dset.createVariable('ens', 'i4', ('ens',))
          dset.createVariable('lat', 'f8', ('lat',))
          dset.createVariable('lon', 'f8', ('lon',))
          dset.variables['time'][:] = np.array(valid_timestamps)
          dset.variables['ens'][:] = np.arange(1,nmem+1)
          dset.variables['lat'][:] = data.variables['lat'][:]
          dset.variables['lon'][:] = data.variables['lon'][:]
          for v in variables:
            dset.createVariable(v, 'f8', ('time','ens','lat','lon',))
            outvar = np.squeeze(data.variables[v][:])
            dset.variables[v][:] = outvar[:]   

    else:
        return data

if __name__ == '__main__':
    get_gefs_ensemble(['Temperature_height_above_ground_ens'], bounds=(-124.9,-66.8,24.3,49.4))
