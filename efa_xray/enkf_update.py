#!/usr/bin/env python

import numpy as np
from copy import deepcopy

#def enkf_update(xbm, Xbp, Y, H, R, loc=None, inflate=None):
def enkf_update(prior_state,obs,inflate=None,loc=False):
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
            kcov = np.multiply(ob.localize(prior_state, type=loc, full_state=True),kcov)
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
    post_state.update_state_from_array(np.add(xam[:,None],Xap))
    # Return the assimilated observations
    return post_state, obs
    # Return the analysis mean and perturbations
    #return xam, Xap

        
        
