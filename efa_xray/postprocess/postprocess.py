#!/usr/bin/env python

import pandas as pd
from efa_xray.state.ensemble import EnsembleState
#from efa_xray.observation.observation import Observation
import numpy as np

def obs_assimilation_statistics(prior, post, obs):
    """
    Loops through all observations and builds a pandas dataframe
    with statistical info about the observations
    """
    # Make sure these states are the right kind of object
    assert isinstance(prior, EnsembleState)
    assert isinstance(post, EnsembleState)

    # Build a list of dictionaries
    oblist = []
    for ob in obs:
        obd = {}
        obd['validtime'] = ob.time
        obd['flead'] = (ob.time - pd.to_datetime(prior['validtime'].values[0])).total_seconds()/3600
        obd['lat'] = ob.lat
        obd['lon'] = ob.lon
        obd['obtype'] = ob.obtype
        obd['description'] = ob.description
        obd['ob error'] = ob.error
        obd['value'] = ob.value
        obd['assimilated'] = ob.assimilated
        prior_ye = ob.estimate(prior)
        post_ye = ob.estimate(post)
        obd['prior mean'] = prior_ye.mean()
        obd['post mean'] = post_ye.mean()
        obd['prior variance'] = prior_ye.var()
        obd['post variance'] = post_ye.var()
        oblist.append(obd)
    # Build a dataframe from this list of objects
    df = pd.DataFrame(oblist)
    return df

