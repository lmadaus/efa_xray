#!/usr/bin/env python

from datetime import datetime, timedelta
from MesoPy import Meso
from itertools import compress
import json
import urllib2
from efa_xray.observation_class import Observation


token= '####################'



def get_mesowest_obs(time, bounds=None, vars=['t2']):

    # Create mesowest object
    m = Meso(token=token)
    # Check to be sure we have strings since that's what
    # the api wants
    if isinstance(time, datetime):
        timestr = time.strftime('%Y%m%d%H%M')
    else:
        timestr = time

    # Find stations within these bounds
    stations = m.metadata(bbox=bounds)
    #print help(m.metadata()
    # Give us ASOS
    #print stations
    #nets = m.networks()
    #shorts = [n['SHORTNAME'] for n in nets['MNET']]
    #ids = [n['ID'] for n in nets['MNET']]
    #for s,i in zip(shorts, ids):
    #    print s, i
    netname = 'NWS/FAA'
    netid = 1
    #asos = [s['STID'] for s in  stations['STATION'] if s['STID'].startswith('K') and not any(char.isdigit() for char in s['STID'])]
    asos = [s['STID'] for s in stations['STATION'] if int(s['MNET_ID']) == netid and not any(char.isdigit() for char
                                                                                                     in s['STID'])]
    # Thin
    asos = asos[::3]
    #print stations['STATION'][0]
    #print asos

    # Set the variable list
    translated_vars = [ob_translator[v] for v in vars]


    # Make the call
    obs = m.attime(stid=asos, within='15', vars=translated_vars, attime=timestr,\
                       units='temp|K')
    #print obs.keys()
    # Extract the information
    observations = []
    for x in obs['STATION']:
        obd = {}
        stid = x['STID']
        lat = float(x['LATITUDE'])
        lon = float(x['LONGITUDE'])
        elev = float(x['ELEVATION'])
        for v in vars:
            try:
                obval = x['OBSERVATIONS'][ob_translator[v]+'_value_1']['value']
                #time = x['OBSERVATIONS'][vars[0]+'_value_1']['date_time']
                #obtime = datetime.strptime(time, '%Y-%m-%dT%H:%M:%SZ')
                obsject = Observation(value=float(obval),obtype=v,time=time,error=1.0,\
                                      location='{:3.7f},{:3.7f}'.format(lat,lon),
                                        description=stid,localize_radius=1000.0)
            
            except KeyError:
                continue
            observations.append(obsject)
    return observations


ob_translator = {'t2' : 'air_temp',
                 'psfc' : 'pressure',
                 'slp' : 'sea_level_pressure'}


if __name__ == '__main__':
    get_mesowest_obs(datetime(2016,3,8,12), bounds=(-124.9,24.3,-66.8,49.4))
