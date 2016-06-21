#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import pandas as pd
import xarray
import xarray.ufuncs as xu
import netCDF4
import matplotlib
import matplotlib.pyplot as plt
from datetime import datetime, timedelta, time
import pytz
from copy import deepcopy


class EnsembleState(xarray.Dataset):
    """Define an ensemble state vector"""
    
    #def __init__(self):
    #    """ Initialize based on the input """
    #    # Quick check to remove duplicates
    #    #xarray.Dataset.__init__(self, vardict, coords=coorddict)
    #    # Get dimension lengths here for quick reference
    #    pass

    @classmethod
    def from_vardict(cls, vardict, coorddict):
        """
        Initializes a new EnsembleClass object when given
        a dictionary of variables (vardict) and a dictionary
        of coordinates (coorddict).  See documentation for
        xarray on how to create an xarray object for how
        these dictionaries are to be formatted
        """
        newstate = xarray.Dataset(vardict, coords=coorddict)
        newstate.__class__ = cls
        return newstate


    # Functions to get various useful attributes of the state
    def nmems(self):
        return len(self.coords['mem'])
    def ny(self):
        return len(self.coords['y'])
    def nx(self):
        return len(self.coords['x'])
    def ntimes(self):    
        return len(self.coords['validtime'])
    def vars(self):
        return [x for x in self.variables.keys() if x not in ['validtime','lat','lon','mem','x','y']]
    def nvars(self):
        return len(self.vars())
    def nstate(self):    
        return self.ntimes() * self.ny() * self.nx() * self.nvars()
    def shape(self):
        """ Returns the full shape of the DataArray """
        return self.to_array().shape


    def split_state(self,nChunks):
        """ Function to split the state xarray object into nChunks number of
        smaller xarray objects for multiprocessing.  Returns a dictionary of state
        chunks.  This dictionary may be re-fed into the master state using the
        function "reintegrate_state" which will overwrite the master state xarray
        object with the separate parts """
        
        state_chunks = {}
        print("In split state!")
        # Figure out how big each section must be
        bounds = self.chunk_bounds(nChunks)

        # Now separate along the location dimension according to the bounds
        for cnum, bnds in bounds.items():
            state_chunks[cnum] = \
                    xarray_Ensemble_State(state=self.state.isel(location=slice(bnds[0],bnds[1])))
        return state_chunks

    def reintegrate_state(self, state_chunks):
        """ Reintegrate the state vector from various chunks.  The opposite of
        split_state.  This will overwrite self.state """
        num_chunks = len(state_chunks.keys())

        # Get the bounds
        bounds = self.chunk_bounds(num_chunks)

        # Now reset the master state vector
        for cnum, bnds in bounds.items():
            self.state[dict(location=slice(bnds[0],bnds[1]))] = state_chunks[cnum].state



    def chunk_bounds(self, nChunks):
        """ Function to compute the bounding array locations when dividing up a
        state array.  Returns a dictionary."""
        chunk_bounds = {}
        num_locs = self.num_locs()

        # Divide along the locations dimension
        chunk_length = num_locs / (nChunks)
        print("Num locs:", num_locs, "Chunk_length:", chunk_length)

        # Now set up the chunks
        for x in xrange(nChunks):
            if x != (nChunks-1):
                chunk_bounds[x] = (x*chunk_length,(x+1)*chunk_length)
            else:
                chunk_bounds[x] = (x*chunk_length, None)
        return chunk_bounds

        
    def to_vect(self):
        """ Returns an array of the values in a shape of 
        Nstate x Nmems """
        # Make sure mems dimension is last
        return np.reshape(self.transpose('validtime','y','x','mem').to_array().values, (self.nstate(), self.nmems()))
    
    def from_vect(self,instate):
        """ Takes an Nstate x Nmems ndarray and updates the state accordingly"""
        instate = np.reshape(instate,self.shape())
        statearr = self.to_array()
        statearr.values = instate
        self.update(statearr.to_dataset(dim='variable'))

    def ensemble_mean(self):
        """Returns the ensemble mean of the state as a DataArray"""
        return self.mean(dim='mem')
    
    def ensemble_perts(self):
        """Removes the ensemble mean and returns an xarray DataArray        of the perturbations from the ensemble mean"""
        #emean = self.ensemble_mean()
        return self - self.ensemble_mean()
        #return self.state.values

    def ensemble_times(self):
        """ Return the values of the time dimension AS DATETIME OBJECTS """
        return self['validtime'].values

   
    def project_coordinates(self,m):
        """ Function to return projected coordinates given:
            m --> a Basemap instance defining the projection
            Returns:
                gy,gx --> The projected coordinate arrays
        """
        # Get the grid lat lons
        # Make negative lons because Basemap is like that
        lons = deepcopy(np.array(self['lon'].values[:]))
        #lons = self['lon'].values
        lons[lons > 180] = lons[lons > 180] - 360
        gx,gy = m(lons, self['lat'].values)
        return gx, gy
    
    def nearest_points(self, lat, lon, npt=1):
        """
        Use the lat-lon arrays to return a list of indices
        of the nearest npt points to the given lat-lon
        """
        # Use sin of lat lon to handle periodic
        # and not worry about if we are in negative
        # degrees
        dist = xu.hypot(xu.sin(xu.radians(self['lat'].values)) -
                 xu.sin(xu.radians(lat)),\
                 xu.cos(xu.radians(self['lon'].values)) - 
                 xu.cos(xu.radians(lon)))
        # Get indices of the flattened array
        nearest_raw = dist.argsort(axis=None)[:npt]
        # Convert back to 2-d coords
        nearest = np.unravel_index(nearest_raw, self['lat'].shape)
        return nearest

    def interpolate(self, var, time, lat, lon):
        """
        Given a variable, lat, lon and time,
        interpolate the state to that point
        """

        # Get the nearest four points in space
        # Check to see if lat/lons are 2d or 1d
        if len(self['lat'].shape) == 2:
            closey, closex = self.nearest_points(lat, lon, npt=4)
            # Distances in km
            distances = np.array([self.haversine(
                            (self['lat'][y,x].values, self['lon'][y,x].values),
                               (lat, lon)) for y,x in 
                               zip(list(closey), list(closex))])
        else:
            closen = self.nearest_points(lat, lon, npt=4)
            closey = closen
            closex = closen
            # Distances in km
            distances = np.array([self.haversine(
                            (self['lat'][n].values, self['lon'][n].values),
                               (lat, lon)) for n in list(closen)])
        # Check for exact match (within some tolerance)
        spaceweights = np.zeros(distances.shape)
        if (distances < 1.0).sum() > 0:
            spaceweights[:,distances.argmin()] = 1
        else:
        # Here, inverse distance weighting (for simplicity)
            spaceweights = 1.0 / distances
            spaceweights /= spaceweights.sum()
        # Get weights in time
        time64 = np.datetime64(time)
        valids = self['validtime'].values
        timeweights = np.zeros(valids.shape)
        # Check if we are outside the valid time range
        if (time64 < valids[0]) or (time64 > valids[-1]):
            print("Interpolation is outside of time range in state!")
            return None
        # Find where we are in this list
        lastdex = (valids >= time64).argmax()
        # If we match a particular time value, then
        # this is just an identity
        if valids[lastdex] == time64:
            # Just make a one at this time
            timeweights[lastdex] = 1
        else:
            # Linear interpolation
            diff  = (valids[lastdex] - valids[lastdex-1])
            totsec = np.abs(diff / np.timedelta64(1, 's'))
            thisdiff = time64 - valids[lastdex]
            thissec = np.abs(thisdiff / np.timedelta64(1,'s'))
            # Put in appropriate weights
            timeweights[lastdex] = float(thissec) / totsec
            timeweights[lastdex-1] = 1.0 - (float(thissec)/totsec)
        # Now that we have the weights, do the interpolation
        interp = self.variables[var].values[:,closey,closex,:]
        # Do a dot product with the time weights
        # And with the space weights
        if len(interp.shape) == 3:
            interp = (timeweights[:,None,None] * interp).sum(axis=0)
        else:
            interp = (timeweights[:,None,None,None] * interp).sum(axis=0)
            
        if len(interp.shape) == 3:
           interp = (spaceweights[:,:,None] * interp).sum(axis=1)
        else:
            interp = (spaceweights[:,None] * interp).sum(axis=0)
        # Return estimate from all ensemble members
        return interp

    def haversine(self,loc1,loc2):
        """ Use Haversine formula to compute the distance between two lat-lon
        coordinate pairs """
        R = 6371. # Radius of earth in kilometers
        lat1 = np.radians(loc1[0])
        lat2 = np.radians(loc2[0])
        dlat = lat2 - lat1
        dlon = np.radians(loc2[1] - loc1[1])

        a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
        c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
        return R * c

    def distance_to_point(self, lat, lon):
        """
        Use Haversine formula to estimate distances from all
        gridpoints to a given location (lat, lon)
        """
        R = 6371. # Radius of earth in km
        lat = np.radians(lat)
        lon = np.radians(lon)
        dlat = lat - xu.radians(self['lat'].values)
        dlon = lon - xu.radians(self['lon'].values)
        a = xu.sin(dlat/2)**2 + xu.cos(lat) * xu.cos(xu.radians(self['lat'].values)) * \
                xu.sin(dlon/2)**2
        c = 2 * xu.arctan2(xu.sqrt(a), xu.sqrt(1.0-a))
        return R*c

    def save_to_disk(self, filename='ens_state.nc'):
        """
        Dump this object to disk
        """
        self.to_netcdf(filename)
