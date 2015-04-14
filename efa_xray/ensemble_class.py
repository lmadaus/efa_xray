#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import pandas as pd
import xray
from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
from datetime import datetime, timedelta, time
import pytz



class Xray_Ensemble_State:
    """Define an ensemble state vector"""
    def __init__(self, state=None, meta=None, usevars=None, usedims=None):
        """ Initialize based on the input.  We are either given a
        netCDF4 "Dataset" object as "state" with a list of variables and dimensions
        to use OR we are given a numpy ndarray as a state with the variables and
        dimensions specified in the variable "meta" """
        
        if isinstance(state, Dataset):
            """ Read in the dataset here """
            pass
        
        else:
            meta_names = meta.keys()
            meta_names.sort()
            meta_titles = [m[1] for m in meta_names]
            # Be sure this matches the number of dimensions
            # In the array
            assert len(meta_titles) == len(state.shape)
            
            # Be sure that there is a dimension called
            # "mem" so that we know how to format the
            # state later for assimiation
            assert "mem" in meta_titles
            
            
            
            # Make sure that the 'mem' dimension is the
            # last one -- VERY IMPORTANT
            if meta_titles[-1] != 'mem':
                print("Dimension 'mem' is not the last dimension!")
                return None
            
            # We grab all the coordinate values from the
            # dictionary of metadata
            coords = [meta[m] for m in meta_names]
            # Make a DataArray (basically, a labeled
            # numpy ndarray) from the state data with
            # the dimensions and coordinate values specified
            # in "meta"
            self.state = xray.DataArray(state,
                                    dims=meta_titles,
                                    coords=coords)
            
            #Convert self.state to a Dataset instead?



        
    def state_to_array(self):
        """ Returns an array of the values in a shape of 
        Nstate x Nmems """
        # This assumes that the mems dimension is last
        return np.reshape(self.state.values,(self.num_state(), self.num_mems()))
    
    def update_state_from_array(self,instate):
        """ Takes an Nstate x Nmems ndarray and rewrites the state accordingly """
        instate = np.reshape(instate,self.shape())
        self.state.values = instate
    
    def shape(self):
        """ Returns the full shape of the DataArray """
        return self.state.shape
    
    def num_mems(self):
        """Returns number of ensemble members"""
        return self.state.coords['mem'].size
    def num_times(self):
        """ Returns number of times in the ensemble"""
        return self.state.coords['time'].size
    def num_vars(self):
        """ Returns number of variables in the ensemble """
        return self.state.coords['var'].size


    def num_state(self):
        """Returns length of state vector"""
        coord_lengths = [s.shape for v,s in self.state.coords.items()]
        return np.product(coord_lengths)/self.num_mems()
    
    def ensemble_mean(self):
        """Returns the ensemble mean of the state as a DataArray"""
        return self.state.mean(dim='mem')
    
    def ensemble_perts(self):
        """Removes the ensemble mean and returns an Xray DataArray        of the perturbations from the ensemble mean"""
        #emean = self.ensemble_mean()
        return self.state - self.ensemble_mean()
        #return self.state.values

    def ensemble_times(self):
        """ Return the values of the time dimension AS DATETIME OBJECTS """
        tstamps = [(x - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's') \
            for x in list(self.state['time'].values)]
        return [datetime.utcfromtimestamp(x) for x in tstamps]

    def locations(self):
        """ Return an array of latitude/longitude pairs for each location in the
        state """
        # Get the locations from the state
        locations = self.state['location'].values
        latlons = zip([float(x.split(',')[0]) for x in locations],\
                      [float(x.split(',')[1]) for x in locations])
        return latlons
   
    def project_coordinates(self,m,ny,nx):
        """ Function to return projected coordinates given:
            m --> a Basemap instance defining the projection
            ny --> Number of y points
            nx --> Number of x points
            Returns:
                gy,gx --> The projected coordinate arrays
        """
        # Get the grid lat lons
        gridlocs = self.locations()
        latgrid = np.reshape([x[0] for x in gridlocs], (ny,nx))
        longrid = np.reshape([x[1] for x in gridlocs], (ny,nx))
        gx,gy = m(longrid, latgrid)
        return gx, gy


