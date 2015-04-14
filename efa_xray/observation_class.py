#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import pandas as pd
import xray
import cPickle
from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
from datetime import datetime, timedelta, time
import pytz
import sys, os
sys.path.append('../data_parsers')


class Observation:
    def __init__(self,value=None,obtype=None,time=None,error=None,location=None,
                 prior_mean=None, post_mean=None, prior_var=None, post_var=None,
                 assimilate_this=False,description=None,localize_radius=None):
        self.value = value
        self.obtype = obtype
        self.time = time
        self.error = error
        self.location = location
        self.prior_mean = prior_mean
        self.post_mean = post_mean
        self.prior_var = prior_var
        self.post_var = post_var
        self.assimilate_this = assimilate_this
        self.assimilated = False
        self.description = description
        self.localize_radius = localize_radius
        


    def H(self,state):
        """ Given an ensemble state class, compute an H """
        # Make an empty array filled with zeros
        state_shape = state.shape()
        H_values = np.zeros(state_shape[:-1]) # Ignore the members dimension
        #print(H_values.shape)
        #print(state_shape)

        # Now make a new DataArray of this zeros
        # array with the same dimensions as the state
        state_dims = state.state.coords
        state_dim_names = list(state.state.dims[:])
        # Remove the member dimension here
        state_dim_names.remove('mem')
        new_dims = {}
        for dname,dvals in state_dims.items():
            if dname != 'mem':
                new_dims[dname] = dvals
                
        # Make the data array
        H_da = xray.DataArray(H_values, coords=new_dims, dims=state_dim_names)
       

        # Now figure out the weights
        # Time and variable are straightforward (must be identity)
        # Location is the only question
        # Do a simple distance-weighted average of nearest four points
        # Get distance to all points in state
        distances = self.distance_to_state(state)
        # Find the closest four points
        closest_indices = distances.argsort()[:4]
        closest_distances = distances[closest_indices]
        # Develop a weight using inverse distance squared
        # First, see if any points are within 1 meter of the ob
        if np.min(closest_distances) < 0.01:
            useloc = closest_indices[np.where(closest_distances < 0.01)[0]]
            H_da.loc[dict(location=state['location'].values[useloc], time=self.time, var=self.obtype)] = 1.0
        # Otherwise, do the inverse distance squared
        else:
            weights = np.divide(1.0,np.power(closest_distances,2))
            sum_weights = np.sum(weights)
            for index, w in zip(list(closest_indices),list(weights)):
                H_da.loc[dict(location=state.state['location'].values[index],
                              time=self.time, var=self.obtype)] = w / sum_weights
               

        # Return a flattened array of length Nstate containing H
        return np.ravel(H_da.values)

    def distance_to_state(self,state):
        """ Return the distance from this ob to all locations in the state
        vector (in km) """
        # Get all locations in the state
        latlons = state.locations()
        # Compute distances to all points in the state using Haversine formula
        ourloc = [np.float(x) for x in self.location.split(',')]
        distances = np.array([haversine(ourloc, s) for s in latlons])
        # Return the array of distances
        return distances


    def localize(self,state,type='GC',full_state=False):
        """ Given a state vector object, assume location is in lat/lon and compute a
        Gaspari-Cohn weighting function with the specified halfwidth (in km) """
        # Get the localization halfwidth from 
        halfwidth = self.localize_radius

        # Get distance to all points in the state
        distances = self.distance_to_state(state)
        
        # If halfwidth is None, return an array of ones
        if halfwidth is None:
            localization = np.ones(distances.shape)
        else:
        
            r = np.divide(distances, abs(halfwidth))

            # For Gaspari-Cohn
            if type == 'GC':
                # Do the Gaspari Cohn weighting
                localization = np.zeros(r.shape)
                # Less than halfwidth
                localization[r <= 1.0] = ((((-0.25*r+0.5)*r+0.625)*r-5.0/3.0) * r**2 + 1.0)[r <= 1.0]
                # Between halfwidth and fullwidth
                localization[(r > 1.0) & (r < 2.0)] = (((((r/12.0 - 0.5)*r + 0.625) *r+\
                                                5.0/3.0)*r-5.0)*r + 4.0 -\
                                                2.0 / (3.0 * r))[(r > 1.0) & (r < 2.0)]

        if full_state:
            # Return the localization applied to the whole state vector
            ntotal = state.num_times() * state.num_vars()
            return np.tile(localization, ntotal)

        else:
            # Return localization just over xy 
            return localization


    def map_localization(self, state, m, ny, nx, type='GC'):
        """ Function to map localization radius 
        Requires:
            state --> The state vector we are using
            m     --> A basemap instance for projecting the map

            """
        # Get the localization weights and reshape to ny x nx
        localization = self.localize(state, type=type)
        local_shaped = np.reshape(localization, (ny,nx))
        # Get map projected coordinates from the state
        gx, gy = state.project_coordinates(m,ny,nx)
        # Make the plot
        F = plt.figure()
        plt.pcolormesh(gx,gy,local_shaped,vmin=0.0,vmax=1.0)
        m.drawcoastlines()
        m.drawcountries()
        m.drawstates()
        plt.title('Localization Weights for {:s} ({:s})'.format(self.description, self.location))
        plt.show()


    
    def H_Xb(self, state):
        """ Return the ensemble (state) estimate of the ob """
        return np.dot(self.H(state),state.state_to_array())


def haversine(loc1,loc2):
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

