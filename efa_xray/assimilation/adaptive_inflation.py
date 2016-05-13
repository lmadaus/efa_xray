#!/usr/bin/env python
from __future__ import print_function
from efa_xray.state.ensemble import EnsembleState
import xarray
import numpy as np
from pandas import to_datetime

class AdaptiveInflation():
    """
    Class to handle adaptive spatial (and temporal) inflation
    in accordance with:
        Anderson (2009): Spatially and temporally varying adaptive covariance
        inflation for ensemble filters.  Tellus. 61A, 72-83.

    """
    def __init__(self, priorstate, priorinf):
        """
        Tries to load in prior inflation file
        """
        assert isinstance(priorstate, EnsembleState)
        # A couple of options for priorinf
        inftype, infile, initvals = priorinf
        # Try to open the prior file
        # if this fails, build a new inflation
        try:
            self.inflation = xarray.open_dataset(infile)
        except:
            self.inflation = self.build_initial_inflation(priorstate, initvals)



    def build_initial_inflation(self, priorstate, initvals):
        """
        Takes the shape of the priorstate and builds
        a corresponding xarray DataSet object with
        two "members": the mean and variance of 
        the initial inflation given as a tuple (priorval)
        """
        # Get the leadtimes
        valids = priorstate['validtime'].values
        firsttime = to_datetime(valids[0])
        leads = [(to_datetime(x) - firsttime).total_seconds()/3600. for x in list(valids)]

        # Get the variables
        infvars = {}
        for var in priorstate.vars():
            nt,ny,nx,nm = priorstate.variables[var].shape
            outfield = np.ones((nt,ny,nx,2))
            # Set the initial mean
            outfield[:,:,:,0] *= initvals[0]
            # And standard deviation
            outfield[:,:,:,1] *= initvals[1]
            infvars[var] = (['validtime','y','x','moment'], outfield)
        # Build the xarray dataset
        return xarray.Dataset(infvars, {'validtime': leads, 'lat': (['y','x'], priorstate['lat'].values),
                                        'lon' : (['y','x'], priorstate['lon'].values), 'moment': ['mean','std']})
        

    def inflate_state(self, priorstate):
        """
        Uses the currently defined inflation mean to multiply the 
        perturbations from priorstate's mean
        """
        # Need to convert leadtime to validtime in accordance with the priorstate
        self.inflation.coords['validtime'] = priorstate.coords['validtime']

        # Get the mean inflation
        this_inflate = self.inflation[dict(moment=0)]
        # Inflate the perturbations
        this_pert = priorstate - priorstate.mean(dim='mem')
        this_mean = priorstate.mean(dim='mem')
        priorstate = (this_inflate * this_pert) + this_mean
        # Return the inflated state
        return priorstate

    def save_to_disk(self, filename='prior_inflation.nc'):
        """
        Writes the current inflation state out to disk as a netcdf file
        """
        self.inflation.to_netcdf(filename)





