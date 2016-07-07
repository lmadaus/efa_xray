#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from copy import deepcopy
import multiprocessing as mp
from xarray import DataArray, open_dataset



class Assimilation():
    """
    This class contains methods for computing obs priors and
    splitting up state
    """
    def __init__(self, state, obs,nproc=1, inflation=None, verbose=False):
        """
        Set prior, post and obs
        
        inflate option can either be None (no inflation), a float (all state variables
        will be inflated by that amount), a dictionary, or a string file name pointing
        to an npz file of inflation factors.  Dictionary keys may be
        variable names with values as the inflation factor for that variable, or
        dimension names as keys with the values being a numpy array the same length
        as that dimension, with each element being the inflation factor for a given time.

        """
        self.prior = state
        self.obs = obs
        self.post = deepcopy(state)
        self.verbose = verbose
        self.nproc = nproc
        self.inflation = inflation
        self.is_inflated = False


    def compute_ob_priors(self):
        """
        Loop through the observations and compute
        prior means and perts
        """
        nobs = len(self.obs)
        nmems = self.prior.nmems()
        means = np.zeros(nobs)
        perts = np.zeros((nobs, nmems))
        for obnum, ob in enumerate(self.obs):
            ye = ob.estimate(self.prior)
            means[obnum] = ye.mean()
            perts[obnum,:] = ye - ye.mean()
        return means, perts


    def inflate_state(self):
        """
        Inflate the state as specified
        """
        # Check if we're already inflated
        if self.is_inflated:
            print("State already inflated.  Skipping additional inflation.")
            return

        varnames = self.prior.vars()
        if isinstance(self.inflation, float):
            if self.verbose: print("Inflating all variables by factor: {:3.2f}".format(self.inflation))
            # Inflate all variables by this number
            for v in varnames:
                if self.verbose: print(v, "BEFORE stdev:", np.mean(np.std(self.prior.variables[v], axis=-1), axis=None).values)
                self.prior.variables[v][:] = self.prior.ensemble_perts().variables[v] * self.inflation + \
                        self.prior.ensemble_mean().variables[v]
                if self.verbose: print(v, "AFTER stdev:", np.mean(np.std(self.prior.variables[v], axis=-1), axis=None).values)

        elif isinstance(self.inflation, str):
            if self.verbose: print("Trying to load inflation from file: {:s}".format(self.inflation))
            # Load the npy state vector inflation
            with open_dataset(self.inflation) as infile:
                # Update with inflated values.  The beauty of xarray is that it will handle
                # dimension broadcasting with this.  Should error out if the inflation
                # file dimensions don't match the prior state's dimensions
                self.prior = self.prior.ensemble_perts() * infile + self.prior.ensemble_mean()
            if self.verbose: print("Succeeded inflation from file: {:s}".format(self.inflation))
        else:
            # Assume this is a dictionary, so try it as such
            for k, v in self.inflation.items():
                if k in ['validtime','lat','lon','x','y']:
                    if self.verbose: print("Inflating all variables along {:s} dimension".format(k))
                    # This is a dimension. Check to be sure that the inflation size matches the
                    # dimension length
                    Ninflate = v.shape[0]
                    assert Ninflate == len(self.prior.coords[k])
                    # Make the numpy array into a data array
                    thisInflate = DataArray(v, [(k, self.prior.coords[k].values)])
                    # Get the first variable
                    if self.verbose: 
                        priorvars = np.mean(np.mean(np.std(self.prior.variables[varnames[0]],axis=-1), axis=-1),axis=1)
                        print(varnames[0], "BEFORE stdev:", priorvars)
                    # Now multiply these together---xarray should handle broadcasting based on array dimensions
                    self.prior = self.prior.ensemble_perts() * thisInflate + self.prior.ensemble_mean()
                    if self.verbose: 
                        aftervars = np.mean(np.mean(np.std(self.prior.variables[varnames[0]],axis=-1), axis=-1),axis=-1)
                        print(varnames[0], "AFTER stdev:", aftervars)
                        print(aftervars / priorvars)


                else:
                    # This must be a variable name.  Inflate that variable by the
                    # specified value
                    assert isinstance(v, float)
                    if k not in self.prior.variables.keys():
                        print("Unable to find variable {:s} to inflate.  Skipping...")
                        continue
                    if self.verbose: print("Inflating variable {:s} by factor: {:3.2f}".format(k, v))
                    if self.verbose: print(k, "BEFORE stdev:", np.mean(np.std(self.prior.variables[k], axis=-1), axis=None).values)
                    # Update the perturbations for this variable accordingly
                    self.prior.variables[k][:] = self.prior.ensemble_perts().variables[k] * v + self.prior.ensemble_mean().variables[k]
                    if self.verbose: print(k, "AFTER stdev:", np.mean(np.std(self.prior.variables[k], axis=-1), axis=None).values)
                

        # Set the flag to inflated
        self.is_inflated = True

    def format_prior_state(self):
        """
        Builds state into vector
        """
        Nens = self.prior.nmems()   # Number of ensemble members
        Nstate = self.prior.nstate()   # Number of variables in state vector
        Nobs = len(self.obs)     # Total number of observations
        # Get the full state shape
        statesize = self.prior.shape()
        
        
        # Option here to inflate state
        if self.inflation is not None:
            if self.verbose: print("Inflating Prior State")
            self.inflate_state()
                
        # Get ob priors
        if self.verbose: print("Computing observation priors")
        obmeans, obperts = self.compute_ob_priors()

        # Convert state to vector
        if self.verbose: print("Converting state to vector")
        prior = self.prior.to_vect()

        # Divide based on how many processors we need
        #if self.nproc <= 1:
        xbm = prior.mean(axis=1)
        Xbp = prior - xbm[:,None]
        # Compute observation priors and append to state
        xbm = np.hstack((xbm, obmeans))
        Xbp = np.vstack((Xbp, obperts))
        # Need other case here for multiprocessing
        
        
        return xbm, Xbp


    def format_posterior_state(self, xam, Xap):
        """
        Creates posterior state object from analysis mean
        and perturbations
        """
        # Additional work here for repopulating in multiprocessing
        if self.verbose: print("Formatting posterior")
        # Reset the state
        post_state = deepcopy(self.prior)
        # Rebuild the full values
        Nstate = self.prior.nstate()   # Number of variables in state vector
        post = (xam[:,None] + Xap)[:Nstate]
        post_state.from_vect(post)
        # Return the assimilated observations
        return post_state, self.obs



        
def update(prior_state,obs,inflate=None,loc=False,nproc=1, verbose=False):

    # If there is 1 (or none) processors given, then 
    # Don't worry about splitting up the state
    if nproc <= 1:
        posterior_state, posterior_obs = enkf_update(prior_state, obs, inflate=inflate,
                                                     loc=loc, verbose=verbose)
        return posterior_state, posterior_obs


    # If we are splitting up the state, first need to calculate observation
    # priors
    print("Computing observation priors...")
    ob_priors = np.array([ob.H_Xb(prior_state) for ob in obs])
    
    # Splitting function here
    print("Splitting up state among processors...")
    old_state = prior_state.split_state(nproc)

    # Set up the processes
    print("Beginning update")
    processes = []
    output = mp.Queue()
    for p in xrange(nproc):
        process = mp.Process(target=worker, args=(p, old_state[p], obs,\
                                                         inflate, loc,\
                                                    ob_priors, output))
        processes.append(process)
        process.start()

    new_state = []
    stopcount = 0
    for newchunk in iter(output.get, 'STOP'):
        new_state.append(newchunk)
    new_state = dict(new_state)

    # Wait for all to finish
    for p in processes:
        p.join()
    print("Done with update")
    print("Re-assembling state")
    # Get the updated state from queue

    posterior_state = deepcopy(prior_state)
    posterior_state.reintegrate_state(new_state)
    posterior_obs = prior_obs
    return posterior_state, posterior_obs
    

def worker(num, statechunk, allobs, inflate, loc, obs_in_state, output):
    updated_chunk, obout = enkf_update(statechunk,allobs,inflate,loc,obs_in_state=obs_in_state)
    output.put((num, updated_chunk))
    output.put('STOP')
    print("Worker Done!")
    return

