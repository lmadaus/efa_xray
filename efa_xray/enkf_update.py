#!/usr/bin/env python

import numpy as np
from copy import deepcopy
import multiprocessing as mp


def update(prior_state,obs,inflate=None,loc=False,nproc=1):

    # If there is 1 (or none) processors given, then 
    # Don't worry about splitting up the state
    if nproc <= 1:
        posterior_state, posterior_obs = enkf_update(prior_state, obs, inflate=inflate, loc=loc)
        return posterior_state, posterior_obs


    # If we are splitting up the state, first need to calculate observation
    # priors
    ob_priors = np.array([ob.H_Xb(prior_state) for ob in obs])
    
    # Splitting function here
    old_state = prior_state.split_state(nproc)

    


    posterior_state = deepcopy(prior_state)
    posterior_state.reintegrate_state(old_state)
    return posterior_state, posterior_obs
    


def worker(num, statechunk, allobs, nflate, loc):
    updated_chunk, obout = enkf_update(statechunk,allobs,inflate,loc)
    return num, updated_chunk





def enkf_update(prior_state,obs,inflate=None,loc=False,obs_in_state=None):
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
    Nens = prior_state.num_mems()   # Number of ensemble members
    Nstate = prior_state.num_state()   # Number of variables in state vector
    Nobs = len(obs)     # Total number of observations
    

    #Xa = state.state_to_array()

    # Copy the xray state for the output ensemble
    post_state = deepcopy(prior_state)
    xam = np.reshape(post_state.ensemble_mean().values,(Nstate))
    xbm = np.reshape(prior_state.ensemble_mean().values,(Nstate))
    print xam.shape, Nstate
    Xap = np.reshape(post_state.ensemble_perts().values,(Nstate,Nens))
    Xbp = np.reshape(prior_state.ensemble_perts().values,(Nstate,Nens))
    print Xap.shape

    # Check to see if we are appending the obs to the state
    if obs_in_state is not None:
        ob_xam = np.mean(obs_in_state, axis=1)
        ob_Xap = np.subtract(obs_in_state,ob_xam[:,None])
        xam = np.concatenate(xam, ob_xam)
        xbm = np.concatenate(xbm, ob_xam)
        Xap = np.concatenate(Xap, ob_Xap)
        Xbp = np.concatenate(Xbp, ob_Xap)


    # Now loop over all observations
    for obnum,ob in enumerate(obs):
        # Reset the mean and perturbations
        xbm = xam
        Xbp = Xap
        # Make a vector of all of the ensemble members
        # estimates of state translated into observation
        # space (by H)
        #print "Mean", np.tile(xbm,(Nens,1))
        #print np.transpose(np.tile(xbm,(Nens,1))) + Xbp
        #print H
        if obs_in_state is not None:
            H = np.zeros(xam.shape)
            H[Nstate+obnum] = 1.0
        else:
            H = ob.H(prior_state)
        Ye = np.dot(H,np.add(xbm[:,None], Xbp))

        #Ye = ob.H_Xb(post_state)
        # Ye now has shape (1 x num_members)
        #For now, obs are just the first value
        #Ye = np.tile(xbm[0],Nens) + Xbp[0,:]
        #print "After:", Ye
        #raw_input()

        # The ensemble mean of the model estimate (in obs space)
        #mye = np.mean(Ye,axis=1)
        mye = np.mean(Ye)
        # Remove the mean from the model estimate of ob
        ye = np.subtract(Ye, mye)
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
            state_localize = ob.localize(prior_state, type=loc,
                                             full_state=True)
            if obs_in_state is not None:
                # Now need to localize for obs
                obs_localize = ob.localize(obs, type=loc, localizing_obs = True)
                state_localize = np.concatenate(state_localize, obs_localize)
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

    # Reset the state
    post_state.update_state_from_array(np.add(xam[:Nstate,None],Xap[:Nstate,:]))
    # Return the assimilated observations
    return post_state, obs
    # Return the analysis mean and perturbations
    #return xam, Xap

        
        
