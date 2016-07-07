#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from copy import deepcopy
from efa_xray.assimilation.assimilation import Assimilation


class EnSRF(Assimilation):
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
    
    def __init__(self, state, obs, nproc=1, inflation=None, verbose=True, loc=False):
        # Initialize the Assimilation inheritance
        Assimilation.__init__(self, state, obs, nproc, inflation, verbose)
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


            # Option to localize the gain
            if self.loc not in [None, False]:
                # Project the localization
                state_localize = ob.localize(self.prior, type=self.loc)
                #print(state_localize.shape)
                #print(dum_localize.shape)
                # LEM---CHECK TO BE SURE THIS LOGIC WORKS FOR 1-D LATLON!!!
                # This needs to be projected into the full state vector
                # Check for (ny,nx) or just 1-d
                # Project this across all vars, all times
                if len(state_localize.shape) == 2:
                    state_localize = (state_localize[None,None,:,:] * dum_localize).flatten()
                else:
                    state_localize = (state_localize[None,None,None,:] * dum_localize).flatten()
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
        
