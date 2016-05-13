#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from copy import deepcopy
import multiprocessing as mp



class Assimilation():
    """
    This class contains methods for computing obs priors and
    splitting up state
    """
    def __init__(self, state, obs, nproc=1, verbose=False):
        # Set prior, post and obs
        self.prior = state
        self.obs = obs
        self.post = deepcopy(state)
        self.verbose = verbose
        self.nproc = nproc

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


    def format_prior_state(self):
        """
        Builds state into vector
        Splits the state into chunks for potential multiprocessing
        """
        # Get ob priors first
        if self.verbose: print("Computing observation priors")
        obmeans, obperts = self.compute_ob_priors()
        Nens = self.prior.nmems()   # Number of ensemble members
        Nstate = self.prior.nstate()   # Number of variables in state vector
        Nobs = len(self.obs)     # Total number of observations
        
        # Convert state to vector
        if self.verbose: print("Converting state to vector")
        prior = self.prior.to_vect()
        
        # Divide based on how many processors we need
        if self.nproc <= 1:
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

