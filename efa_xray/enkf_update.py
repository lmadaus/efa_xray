#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from copy import deepcopy
import multiprocessing as mp


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






def enkf_update(prior_state,obs,inflate=None,loc=False,verbose=False):
    """
    Function to do the EnSRF update
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
   


    #Nens = np.shape(Xbp)[1]   # Number of ensemble members
    #Ndim = np.shape(Xbp)[0]   # Number of variables in state vector
    #Nobs = np.shape(Y)[0]     # Total number of observations
    if verbose: print("Converting state to vector")
    Nens = prior_state.nmems()   # Number of ensemble members
    Nstate = prior_state.nstate()   # Number of variables in state vector
    Nobs = len(obs)     # Total number of observations
    

    #Xa = state.state_to_array()
    # Convert state to vector
    prior = prior_state.to_vect()
    xam = prior.mean(axis=1)
    xbm = prior.mean(axis=1)
    Xap = prior - xam[:,None]
    Xbp = prior - xbm[:,None]

    # Make a dummy localization here to allocate this
    # array only once
    state_shape = prior_state.shape()[:-1] # Don't include mems
    dum_localize = np.ones(state_shape)
    
    # Check to see if we are appending the obs to the state
    """
    if obs_in_state is not None:
        ob_xam = np.mean(obs_in_state, axis=1)
        ob_Xap = np.subtract(obs_in_state,ob_xam[:,None])
        xam = np.hstack((xam, ob_xam))
        xbm = np.hstack((xbm, ob_xam))
        #print Xap.shape, ob_Xap.shape
        Xap = np.vstack((Xap, ob_Xap))
        Xbp = np.vstack((Xbp, ob_Xap))
    """
    if verbose: print("Computing observation priors")
    # Compute observation priors and append to state
    obmeans = np.zeros(Nobs)
    obperts = np.zeros((Nobs, Nens))
    for obnum, ob in enumerate(obs):
        ye = ob.estimate(prior_state)
        obmeans[obnum] = ye.mean()
        obperts[obnum,:] = ye - ye.mean()
    xam = np.hstack((xam, obmeans))
    xbm = np.hstack((xbm, obmeans))
    Xap = np.vstack((Xap, obperts))
    Xbp = np.vstack((Xbp, obperts))


    if verbose: print("Beginning observation loop")
    # Now loop over all observations
    for obnum,ob in enumerate(obs):
        if (obnum % 100==0) and verbose: print("    On ob:", obnum)
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
        if inflate is not None:
            kcov = inflate * kcov

        # Option to localize the gain
        if loc not in [None, False]:
            # Project the localization
            state_localize = ob.localize(prior_state, type=loc)
            # This needs to be projected into the full state vector
            # Currently is (ny,nx)
            # Project this across all vars, all times
            state_localize = (state_localize[None,None,:,:] * dum_localize).flatten()

            # Now need to localize for obs
            obs_localize = ob.localize(obs, type=loc)
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
    if verbose: print("Ob loop done.  Formatting posterior.")
    # Reset the state
    post_state = deepcopy(prior_state)
    # Rebuild the full values
    post = (xbm[:,None] + Xbp)[:Nstate]
    post_state.from_vect(post)
    # Return the assimilated observations
    return post_state, obs
    # Return the analysis mean and perturbations
    #return xam, Xap

        
        
