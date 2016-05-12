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










class EnSRFUpdate(Assimilation):
    """
    Class to do the EnSRF update
    Originator: G. J. Hakim

     Input variables are
     Xbp => Array of ensemble estimates of state perturbations from mean (num_state_vars x num_ens_mems)
     xbm => ensemble mean (vector of length num_state_vars)
     Y => Observations (num_observations)
     H => Translation matrix from state vector to obs (num_observations x num_state_vars)
     R => Observation error covariance matrix (num_observations x num_observations)
     loc => Localization
    #print "Y shape:", np.shape(Y)
    #print "H shape:", np.shape(H)
    #print "R shape:", np.shape(R)
    #print "loc shape:", np.shape(loc)

    # Modified to work with XRAY enkf --> L. Madaus 2/20/2015
    """
    
    def __init__(self, state, obs, nproc=1, verbose=False, inflate=None, loc=False):
        # Initialize the Assimilation inheritance
        Assimilation.__init__(self, state, obs, nproc, verbose)
        self.inflate = inflate
        self.loc = loc

    def update(self):
        if self.verbose: print("Beginning update sequence")
        # Make a dummy localization here to allocate this
        # array only once
        state_shape = self.prior.shape()[:-1] # Don't include mems
        dum_localize = np.ones(state_shape)
        Nstate = self.prior.nstate()
        Nens = self.prior.nmems()

        # Do pre-processing to estimate obs and format
        # as state vector
        xam, Xap = self.format_prior_state()

        # Now loop over all observations
        if self.verbose: print("Beginning observation loop")
        for obnum,ob in enumerate(self.obs):
            if (obnum % 100==0) and self.verbose: print("    On ob:", obnum)
            # Reset the mean and perturbations
            xbm = xam
            Xbp = Xap
            # Make a vector of all of the ensemble members
            # estimates of state translated into observation
            # space (by H)
            #print "Mean", np.tile(xbm,(Nens,1))
            #print np.transpose(np.tile(xbm,(Nens,1))) + Xbp
            #print H
            H = np.zeros(xam.shape)
            H[Nstate+obnum] = 1.0
            mye = np.dot(H, xbm)
            ye = np.dot(H, Xbp)

            ob.prior_mean = mye
            #print "ye", ye
            # Find the variance among the ensemble members
            varye = np.var(ye)
            ob.prior_var = varye

            # IMPORTANT --- here we check to see if we should actually
            # assimilate this ob
            if not ob.assimilate_this:
                ob.assimilated = False
                continue


            # And find the observation error variance from the R matrix
            # (Assumes ob errors are uncorrelated)
            obs_err = ob.error

            # Find the innovation --the difference between the ob value
            # and the ensemble mean estimate of the ob
            # This is y-HXb
            innov = ob.value - mye

            # Now find the innovation variance -- the sum of the variance of the ob
            # and the varaiance of the ensemble estimate
            # This goes into the denominator of the Kalman gain
            kdenom = (varye + obs_err)

            # The numerator of the Kalman gain is the covariance between
            # the ensemble members and the obs-transformed ensemble members
            kcov = np.dot(Xbp,np.transpose(ye)) / (Nens-1)

            # Option to inflate the covariances by a certain factor
            if self.inflate is not None:
                kcov = self.inflate * kcov

            # Option to localize the gain
            if self.loc not in [None, False]:
                # Project the localization
                state_localize = ob.localize(self.prior, type=self.loc)
                # This needs to be projected into the full state vector
                # Currently is (ny,nx)
                # Project this across all vars, all times
                state_localize = (state_localize[None,None,:,:] * dum_localize).flatten()

                # Now need to localize for obs
                obs_localize = ob.localize(self.obs, type=self.loc)
                state_localize = np.hstack((state_localize, obs_localize))
                kcov = np.multiply(state_localize,kcov)
                #kcov = np.dot(kcov,np.transpose(loc[ob,:]))
   
            # Compute the Kalman gain
            kmat = np.divide(kcov, kdenom)
            #kmat = np.divide(kcov,kdenom)
            #print "kmat", kmat.shape
            #print "innov", innov.shape

            # Now do the updates
            # First update the mean
            #xam = xbm + np.dot(np.dot(H,kmat),innov)
            #print "kmat", np.shape(kmat)
            #print "innov", np.shape(innov)
            #xam = xbm + np.dot(kmat,innov)
            xam = xbm + np.multiply(kmat,innov)

            # And each ensemble member perturbation
            # This is the "Square Root" 
            # step in the Kalman filter equations
            beta = 1./(1. + np.sqrt(obs_err/(varye+obs_err)))
            kmat = np.multiply(beta,kmat)

            ye = np.array(ye)[np.newaxis]
            kmat = np.array(kmat)[np.newaxis]

            Xap = Xbp - np.dot(kmat.T, ye)

            # For reference, grab the post mean and variance
            post_ye = np.dot(H,xam)
            post_var = np.var(np.dot(H,Xap))
            ob.post_mean = post_ye
            ob.post_var = post_var
            # Record that this ob was assimilated
            ob.assimilated = True
        # After having assimilated everything, rebuild the state
        return self.format_posterior_state(xam, Xap)
        
        
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

