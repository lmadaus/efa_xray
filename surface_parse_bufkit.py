#!/usr/bin/env python
import itertools
from math import sqrt


def split_sref(file,psu_format = False):
    """ Simple routine to split a single sref file into one for each member """
    if psu_format:
        memids = {'MEAN MEAN MEAN MEAN' : 'srefmean',
              'CTL 0 ARW KF' : 'em_ctl',
              '-PER 1 ARW KF' : 'em_n1',
              '+PER 1 ARW KF' : 'em_p1',
              '-PER 2 ARW KF' : 'em_n2',
              '+PER 2 ARW KF' : 'em_p2',
              '-PER 3 ARW KF' : 'em_n3',
              '+PER 3 ARW KF' : 'em_p3',
              'CTL 0 NMB BMJ' : 'nmb_ctl',
              '-PER 1 NMB BMJ' : 'nmb_n1',
              '+PER 1 NMB BMJ' : 'nmb_p1',
              '-PER 2 NMB BMJ' : 'nmb_n2',
              '+PER 2 NMB BMJ' : 'nmb_p2',
              '-PER 3 NMB BMJ' : 'nmb_n3',
              '+PER 3 NMB BMJ' : 'nmb_p3',
              'CTL 0 NMM BMJ' : 'nmm_ctl',
              '-PER 1 NMM BMJ' : 'nmm_n1',
              '+PER 1 NMM BMJ' : 'nmm_p1',
              '-PER 2 NMM BMJ' : 'nmm_n2',
              '+PER 2 NMM BMJ' : 'nmm_p2',
              '-PER 3 NMM BMJ' : 'nmm_n3',
              '+PER 3 NMM BMJ' : 'nmm_p3'}
    else:

        memids = {'MEAN MEAN MEAN MEAN' : 'srefmean',
              'CTL 0 ARW KF # RAP' : 'em_ctl',
              '-PER 1 ARW KF # RAP' : 'em_n1',
              '+PER 1 ARW KF # RAP' : 'em_p1',
              '-PER 2 ARW BMJ # RAP' : 'em_n2',
              '+PER 2 ARW BMJ # RAP' : 'em_p2',
              '-PER 3 ARW BMJ # RAP' : 'em_n3',
              '+PER 3 ARW BMJ # RAP' : 'em_p3',
              'CTL 0 NMMB BMJ # NDAS' : 'nmb_ctl',
              '-PER 1 NMMB BMJ # NDAS' : 'nmb_n1',
              '+PER 1 NMMB BMJ # NDAS' : 'nmb_p1',
              '-PER 2 NMMB SAS # NDAS' : 'nmb_n2',
              '+PER 2 NMMB SAS # NDAS' : 'nmb_p2',
              '-PER 3 NMMB WSM6 # NDAS' : 'nmb_n3',
              '+PER 3 NMMB WSM6 # NDAS' : 'nmb_p3',
              'CTL 0 NMM BMJ # GFS' : 'nmm_ctl',
              '-PER 1 NMM BMJ # GFS' : 'nmm_n1',
              '+PER 1 NMM BMJ # GFS' : 'nmm_p1',
              '-PER 2 NMM SAS # GFS' : 'nmm_n2',
              '+PER 2 NMM SAS # GFS' : 'nmm_p2',
              '-PER 3 NMM KF # GFS' : 'nmm_n3',
              '+PER 3 NMM KF # GFS' : 'nmm_p3'}

    siteid = file[5:9]
    # Open the main file
    f = open(file, 'r')
    outfile = None
    for line in f:
        if line.strip() in memids.keys():
            if outfile != None:
                outfile.close()
                del outfile
            outfile = open('%s_%s.buf' % (siteid.upper(), memids[line.strip()]), 'w')
            print >>outfile, line,
        else:
            if outfile != None:
                print >>outfile, line,
    outfile.close()
    f.close()






def lamp_parser(siteid):
    """ Goes to web and grabs the latest LAMP GFS text output for the
    site and then parses it, returning a dictionary of forecast dates
    and values """
    
    import urllib2, re   # urllib2 opens web data, re is regular expressions
    from datetime import datetime, timedelta # Our dates will be formatted as datetimes
    from numpy import nan,sin,cos,pi  # Some times don't have a value
    from profile_class import Profile


    # Start by manually specifying the web adress
    lamp_addr = 'http://www.nws.noaa.gov/cgi-bin/lamp/getlav.pl?sta=%s' % siteid.upper()

    # Open the latest LAMP forecast and read it 
    website = urllib2.urlopen(lamp_addr)
    html = website.readlines()

    data_lines = {}
    in_block = False
    for line in html:
        # We want to find the line that starts with the siteid
        # and append all subsequent lines to the data_lines dictionary
        # Until we hit the </PRE> tag
        if not in_block:
            if line.strip().startswith(siteid.upper()):
                in_block = True
                #linesplit = line.split()
                # Get the model initialization time
                # Update --10/8/2012-- uses re now
                #startdate = ' '.join(['0'+linesplit[5],linesplit[6]])
                #startdt = datetime.strptime(startdate,'%m/%d/%Y %H00')
                smonth,sday,syear,shour,smin = re.search('(\d{1,2})/(\d{1,2})/(\d{4})  (\d{1,2})(\d{2})', line).groups()
                startdt = datetime(int(syear),int(smonth),int(sday),int(shour),int(smin))
        else:
            # The next blank line we encounter is the end of
            # the useful data
            if len(line.strip().split()) == 0:
                in_block = False
            else:
                # Split the data line and store the list of
                # values keyed by the first element of the list
                linesplit = line.strip().split()
                data_lines[linesplit[0]] = linesplit[1:]
    # Now we have the data, we can re-arrrange the dictionary to make it nicer
    outdict = {}
    # Each subsequent time is actually 1 hour interval, so build on the initialization
    # time to get an hour list
    datelist = []
    p06_list = []
    p06_ind = 0
    for fhour,hourind in zip(data_lines['UTC'],range(len(data_lines['UTC']))):
        # Add on the date
        curdate = startdt+timedelta(hours=1*(hourind+1))
        datelist.append(curdate)

        outdict[curdate] = Profile()
        # If the current hour is in 0,6,12,18, then add the P06 value
        # Else, put np.nan
        if int(fhour) in (0,6,12,18):
            try:
                p06_list.append(int(data_lines['P06'][p06_ind]))
                outdict[curdate].p01m = int(data_lines['P06'][p06_ind])
                p06_ind = p06_ind + 1
            except:
                p06_list.append(nan)
                outdict[curdate].p01m = nan
        else:
            p06_list.append(nan)
            outdict[curdate].p01m = nan

    # Now loop through each variable list and change things to ints if they
    # should be ints
    for var in data_lines.keys():
        if var == 'TMP':
            for time,value in zip(datelist,data_lines[var]):
                outdict[time].tmpc = (float(value) - 32.) * 5. / 9.
        elif var == 'DPT':
            for time,value in zip(datelist,data_lines[var]):
                outdict[time].dwpc = (float(value) - 32.) * 5. / 9.
        elif var == 'WDR':
            for time,value in zip(datelist,data_lines[var]):
                outdict[time].wdir = float(value)
        elif var == 'WSP':
            for time,value in zip(datelist,data_lines[var]):
                outdict[time].sknt = float(value)
        elif var == 'WGS':
            for time,value in zip(datelist,data_lines[var]):
                if value == 'NG':
                    pass
                else:
                    outdict[time].sknt = float(value)
        else:
            pass

    
        """
        if var in ('TMP','DPT','WDR','WSP','PPO','TP2','POZ','POS','CIG','VIS','CVS'):
            newvars = [int(v) for v in data_lines[var]]
            outdict[var] = newvars
        elif var == 'P06':
            outdict['P06'] = p06_list
        elif var == 'UTC':
            outdict['FDATE'] = datelist
        elif var == 'WGS':
            # Wind gusts should change 'NG' to nan
            newvars = [nan if v == 'NG' else int(v) for v in data_lines[var]]
            outdict['WGS'] = newvars
        else:
            # For all others, just return the text
            outdict[var] = data_lines[var]
        """

    # Quickly compute the uv winds
    for d in datelist:
        outdict[d].uwnd = -1 * outdict[d].sknt * (-1 * sin(outdict[d].wdir * \
                                                      pi/180.))
        outdict[d].vwnd = -1 * outdict[d].sknt * (-1 * cos(outdict[d].wdir * \
                                                      pi/180.))
    # Now return the outdict
    return outdict


def climo_parser(siteid,wfoid):
    # Get the verification values for the 24 hour period ending
    # on the most recent day.  For the forecast competition,
    # only the wind comes from this report, but this will
    # grab winds and temperature anyway

    import urllib2, re   # urllib2 to open the web data efficiently, re is regular expressions
    from datetime import datetime, timedelta, date   # datetime tools for efficient handling of dates
    from profile_class import Profile   # the profile class written to store the ata

    # First, need a month dictionary to get
    # the correct date number from a literal month
    monthd = {'JAN' : 1,
              'FEB' : 2,
              'MAR' : 3,
              'APR' : 4,
              'MAY' : 5,
              'JUN' : 6,
              'JUL' : 7,
              'AUG' : 8,
              'SEP' : 9,
              'OCT' : 10,
              'NOV' : 11,
              'DEC' : 12}


    # Have to grab climo data from a website
    baseurl = 'http://www.nws.noaa.gov/view/validProds.php?prod=CLI&node=%s' % wfoid.upper()
    # Experimental website that requires form submission (not working...)
    #baseurl = 'http://www.nws.noaa.gov/climate/index_nonjs.php?wfo=cle'
    


    # Open the latest batch of climo reports
    website = urllib2.urlopen(baseurl)
    html = website.read()
    
    # Each climo report begins with the phrase 'CDUS', so split
    # the output buffer into blocks based on that
    blocks = html.split('CDUS')
    #print blocks
    #raw_input()

    # Set up a regular expression string to find the date line
    dateline = '(\d{3,4}) (AM|PM) \D{3} \D{3} (\D{3}) (\d{1,2}) (\d{4})'

    # Start with an empty profile and profile dictionary
    # Profile is a user-defined class imported from another file (see above)
    cur_profile = Profile()
    profile_dir = {}
   
    for block in blocks:
        # Loop through each block
        if re.search('CLI%s' % (siteid.upper()[1:]),block):
            #print block
            #raw_input()
            # We have the correct site, now only grab it if
            # it doesn't contain "VALID TODAY"--those are not the final reports
            # CHANGE THIS
            if not re.search('VALID AS OF',block) and not re.search('VALID TODAY AS OF',block):
                # Find the date line as specified above
                dateg = re.search(dateline,block).groups()
                # Make the date and subtract a day -- these are the results from the previous day
                validdate = date(int(dateg[4]),monthd[dateg[2]],int(dateg[3]))
                validdate = validdate - timedelta(days=1)

                # Make a new, clean profile to store the data
                del cur_profile
                cur_profile = Profile()
  
                # Now find the maximum wind in the block
                try:
                    maxwl = re.search('HIGHEST WIND SPEED    ( \d{1}|\d{2})',block).groups()
                    maxwind = int(round(float(maxwl[0]) * 0.868976))
                except:
                    maxwind = ''
                cur_profile.sknt = maxwind
              
                # And the maximum temp
                print siteid.upper(), validdate
                maxtl = re.search('MAXIMUM        ( -?\d{2}|  -?\d{1}|-?\d{3})', block).groups()
                maxtemp = int(maxtl[0])
                cur_profile.maxT = maxtemp

                # And the minimum temp
                
                maxtl = re.search('MINIMUM        ( \d{2}|\d{3}| -\d{1}|  \d{1}|-\d{2})', block).groups()
                mintemp = int(maxtl[0])
                cur_profile.minT = mintemp

                # Print the found values just to be sure this worked
                #print validdate.strftime('%Y%m%d'), 'MAX:', maxtemp, 'MIN:', mintemp, 'WIND:', maxwind
                # Store this set of values in the profile dictionary keyed by the date
                profile_dir[validdate] = cur_profile
    # Return the profile dictionary
    return profile_dir


def mos_parser(siteid,model):
    # Parse out the latest mos forecasts for a site
    from profile_class import Profile
    import re,os
    from datetime import datetime,timedelta

    # Different mos to load -- these are the commands to use
    mosd = {'gfs' : '/usr/local/wx/bin/gmos',
            'nam' : '/usr/local/wx/bin/emos'}

    # Load the data by opening a process on the system using
    # the correct command from the dictionary above and reading the output
    p = os.popen('%s %s' % (mosd[model.lower()],siteid))
    output = p.readlines()
    #print output

    # Go through each line and split the line into individual elements.
    # All elements are separated by spaces except for the
    # Date line which is separated by slashes
    
    lined = {}
    for line in output:
        if line.startswith('DT'):
            spt = line.split('/')       
        else:
            spt = line.split()
        if len(spt) > 1:
            # We don't need to keep the row identifier
            # in each row, but store the rest in a dictionary
            # keyed by that line identifier
            if spt[0] in ['TMP','DPT']:
                # Check for misaligned negatives
                if len(spt[1:]) < len(lined['HR']):
                    newlist = []
                    for r in spt[1:]:
                        if r.startswith('-') and len(r) > 3:
                            vals = r.split('-')
                            for v in vals:
                                if v != '':
                                    newlist.append('-' + v)
                        else:
                            newlist.append(r)
                else:
                    newlist = spt[1:]
            else:
                newlist = spt[1:]

            lined[spt[0]] = newlist

    # Uncomment to print out the data from each line in list form
    #for key in lined.keys():
    #    print key, lined[key]

    # Strip out the date of initialization
    longdate = lined[siteid.upper()][3] + lined[siteid.upper()][4]
    init_date = datetime.strptime(longdate,'%m/%d/%Y%H%M')
    #print init_date

    # We know output is every three hours, so go through
    # every time in the HR list
    profile_dir = {}
    validtime = init_date
    cur_profile = Profile()
    for hr in range(len(lined['HR'])):
        del cur_profile
        cur_profile = Profile()
        # Check to see if it's a new day
        if lined['HR'][hr] == '00':
            validtime = validtime + timedelta(hours=24)
        # Update the date
        validtime = validtime.replace(hour=int(lined['HR'][hr]))
        #print validtime
        #raw_input()
        # Set the cur_profile values accordingly for this date and time
        tmpf = float(lined['TMP'][hr])
        dwpf = float(lined['DPT'][hr])
        cur_profile.tmpc = (tmpf-32.) * 5. / 9.
        cur_profile.dwpc = (dwpf-32.) * 5. / 9.
        cur_profile.sknt = float(lined['WSP'][hr])
        # Store this profile in the profile_dir keyed by the date
        # at which the forecast is valid
        profile_dir[validtime] = cur_profile
        
    # Sort out the highs and lows and precip
    # Precip has an odd coded structure with maxes and mins
    # As of now, this script doesn't actually handle that
    precip_vals = {0 : (0.0,0.0),
                   1 : (0.01,0.09),
                   2 : (0.1,0.24),
                   3 : (0.25,0.49),
                   4 : (0.50,0.99),
                   5 : (1.00,1.99),
                   6 : (2.00,2.00)}
 
    precip_num = 0
    # Numpy for access to the concept of nan
    import numpy as np

    # Blank lists will be populated later
    precip_hours = []
    hilo_hours = []
    # Get a sorted list of all the forecast times in the profile_dir
    sorted = profile_dir.keys()
    sorted.sort()
    for date in sorted:
        # Loop through all the times
        # and find which hours are associated with precip
        # forecasts and hi/lo temp forecasts
        # Precip every 6 hours
        # hi/lo every 12 hours on 00 and 12
        if date.hour in [0,6,12,18]:
            precip_hours.append(date)
            if date.hour in [0,12]:
                hilo_hours.append(date)


    # Set up a max/min dictionary
    maxmin_dir = {}
    # Same format as in the obs, but including
    # the lead time in days as the third value
    # an entry would be maxmin_dir[date] = [maxT,minT,lead]
    # May change this later to include precip values and wind values
    # instead

    # Make a blank list of three elements for each date in the forecast
    # and then put these in the maxmin dictionary, keyed by valid date
    for n in range(len(hilo_hours)):
        date = hilo_hours[n].date()
        maxmin_dir[date] = [np.nan,np.nan,np.nan]

    #print hilo_hours
    # Now remove the first of both lists -- mos doesn't give
    # the value for the first time (no idea why...)
    precip_hours = precip_hours[1:]
    hilo_hours = hilo_hours[1:] 

    init_date = init_date.date()
    hilo_num = 0

    # We have to figure out if the key is 'N/X' or 'X/N'
    if 'N/X' in lined.keys():
        xn = 'N/X'
    else:
        xn = 'X/N'
    # Now loop through and start setting things
    for n in range(len(hilo_hours)):
        date = hilo_hours[n]
        if date.hour == 0:
            # The maximum temperature is stored in the 0 hour
            maxmin_dir[date.date()-timedelta(days=1)][0] = float(lined[xn][hilo_num])
        elif date.hour == 12:
            # The minimum temperature is in the 12 hour
            maxmin_dir[date.date()][1] = float(lined[xn][hilo_num])
        else:
            print "ERROR!"
            exit(1)
        hilo_num = hilo_num + 1
     
    # Now loop through and write the correct lead time
    for date in maxmin_dir.keys():
        delta = date - init_date
        maxmin_dir[date][2] = delta.days

    # Addition here to check to see that there is no higher or lower temperature
    # on each day
    for day in maxmin_dir.keys():
        start_time = datetime(day.year,day.month,day.day,6)
        end_time = start_time + timedelta(hours=24)
        # Get all temperatures from times valid for that day
        #good_dates = [d for d in hourly.keys() if d >= start_time and d <= end_time]
        #good_dates.sort()
        #temps = [(hourly[d].tmpc * 9./5.) + 32. for d in good_dates]
        #for d,t in zip(good_dates,temps):
        #    print d, t
        temps = [(data.tmpc * 9./5.) + 32. for d, data in profile_dir.iteritems() if d >= start_time and d <= end_time]

        if temps == []:
            # Hourly data doesn't include this day
            continue
        # Find the max and min
        mintemp = min(temps)
        maxtemp = max(temps)
        given_max = maxmin_dir[day][0]
        given_min = maxmin_dir[day][1]
        # Adjust these values if needed
        if maxtemp > given_max:
            print "   WARNING: Adjusting high for %s: %3.0f to %3.0f" %\
                (day.strftime('%m/%d/%Y'),given_max, maxtemp)
            maxmin_dir[day][0] = maxtemp
        if mintemp < given_min:
            print "   WARNING: Adjusting low for %s: %3.0f to %3.0f" %\
                (day.strftime('%m/%d/%Y'),given_min, mintemp)
            maxmin_dir[day][1] = mintemp





    # Option to print the max/min results
    #print lined['N/X']
    #raw_input()
    # Return both the full 3-hourly forecasts and the max/mins
    return profile_dir,maxmin_dir

        



def bufkit_parser(file, include_severe=False):
    # Parses out the surface data from bufkit profiles
    # of various models.
    from profile_class import Profile
    import re
    from datetime import datetime
    
    # Load the file
    infile = open(file,'r')

    profile_dir = {}
    cur_profile = Profile()
    validtime = '' 

    # Find the block that contains the description of
    # what everything is (header information)
    if include_severe:
        severe_values = {}
        severe_values['CAPE'] = []
        severe_values['CINH'] = []
        severe_values['LCLP'] = []
    block_lines = []
    inblock = False
    for line in infile:
        if re.search('SELV',line):
            elev = re.search('SELV = (\d{1,4})',line).groups()[0]
            elev = float(elev)
        if re.search('CAPE',line) and 'STNPRM' not in line and include_severe:
            cape = re.search('CAPE = (\d{1,4}.\d{1,2})',line).groups()[0]
            severe_values['CAPE'].append(float(cape))
            lclp = re.search('LCLP = (\d{1,4}.\d{1,2})',line).groups()[0]
            severe_values['LCLP'].append(float(lclp))
        if re.search('CINS',line) and 'STNPRM' not in line and include_severe:
            cinh = re.search('CINS = (-?\d{1,4}.\d{1,2})',line).groups()[0]
            severe_values['CINH'].append(float(cinh))

        if line.startswith('STN YY'):
            # We've found the line that starts the header info
            inblock = True
            block_lines.append(line)
        elif inblock:
            # Keep appending lines until we start hitting numbers
            if re.match('^\d{6}',line):
                inblock = False
            else:
                block_lines.append(line)

    #print block_lines
    # Get the station elevation


    # Build an re search pattern based on this
    # We know the first two parts of the section are station id num and date
    re_string = "(\d{6}) (\d{6})/(\d{4})"

    # Now compute the remaining number of variables
    dum_num = len(block_lines[0].split()) - 2
    for n in range(dum_num):
        re_string = re_string + " (-?\d{1,4}.\d{2})"
    re_string = re_string + '\r\n'    
    for line in block_lines[1:]:
        dum_num = len(line.split())
        for n in range(dum_num):
            re_string = re_string + '(-?\d{1,4}.\d{2}) '
        re_string = re_string[:-1]  # Get rid of the trailing space
        re_string = re_string + '\r\n'

    # If you want to see what the above code managed to put together
    # as a regular expression search pattern, uncomment this:
    #print re_string
    #raw_input()

    # Compile this re_string for more efficient re searches
    block_expr = re.compile(re_string)

    # Now get corresponding indices of the
    # variables we need
    full_line = ''
    for r in block_lines:
        full_line = full_line + r[:-2] + ' '
    # Now split it
    varlist = re.split('[ /]',full_line)
    
    # To see the variable list, uncomment
    #print varlist
    #raw_input()

    with open(file) as infile:
       # Now loop through all blocks that match the
       # search pattern we definied above
       blocknum = -1 
       for block_match in block_expr.finditer(infile.read()):
            blocknum += 1
            #print "Match found"
            # For each match, make a fresh profile
            del cur_profile
            cur_profile = Profile()
            cur_profile.elev = elev 
            # Split out the match into each component number
            vals = block_match.groups()
            # Set the time
            dt = '20' + vals[varlist.index('YYMMDD')] + vals[varlist.index('HHMM')]
            validtime = datetime.strptime(dt,'%Y%m%d%H%M')
            # Have to manually compute the wind
            # What's nice is because we made a list of variable names from the header
            # information that exactly matches the number of values we get after splitting
            # each matched block into its component numbers, the index of each variable name
            # in the varlist list corresponds with the index of the corresponding value in the
            # list of components in each block.  This makes the script flexible.  Also explains the
            # vals[varlist.index()]] notation--get the value from the index that matches the index of
            # the variable name we want
            uwind = float(vals[varlist.index('UWND')])
            vwind = float(vals[varlist.index('VWND')])
            wspd = sqrt(uwind**2 + vwind**2)
            cur_profile.sknt = wspd
            cur_profile.uwnd = uwind
            cur_profile.vwnd = vwind
            cur_profile.press = float(vals[varlist.index('PRES')])
            cur_profile.tmpc = float(vals[varlist.index('T2MS')])
            cur_profile.dwpc = float(vals[varlist.index('TD2M')])
            hcld = float(vals[varlist.index('HCLD')])
            mcld = float(vals[varlist.index('MCLD')])
            lcld = float(vals[varlist.index('LCLD')])
            cur_profile.cfrl = int((hcld + mcld + lcld) / 3.0)
            # Could be 3 hour or 1 hour precip -- store them both as p01m
            try:
                cur_profile.p01m = float(vals[varlist.index('P01M')])
            except:
                cur_profile.p01m = float(vals[varlist.index('P03M')])
            # If including severe parameters, grab them
            if include_severe:
                cur_profile.CAPE = severe_values['CAPE'][blocknum]
                cur_profile.CINH = severe_values['CINH'][blocknum]
                cur_profile.LCLP = severe_values['LCLP'][blocknum]
            # Store the profile in the profile_dir, keyed by the valid time            
            profile_dir[validtime] = cur_profile
            
    #raw_input()
    return profile_dir
    #import pickle
    #pickle.dump(profile_dir,open('profiletest.pickle','w'))



def bufkit_parser_time_height(file):
    import re
    import numpy as np
    from scipy import interpolate
    from datetime import datetime

    # Open the file  
    infile = open(file,'r')

    profile_dir = {}

    # Find the block that contains the description of
    # what everything is (header information)
    block_lines = []
    inblock = False
    block_found = False
    for line in infile:
        if line.startswith('PRES TMPC') and not block_found:
            # We've found the line that starts the header info
            inblock = True
            block_lines.append(line)
        elif inblock:
            # Keep appending lines until we start hitting numbers
            if re.match('^\d{3}|^\d{4}',line):
                inblock = False
                block_found = True
            else:
                block_lines.append(line)

    #print block_lines
    # Now compute the remaining number of variables
    re_string = ''
    for line in block_lines:
        dum_num = len(line.split())
        for n in range(dum_num):
            re_string = re_string + '(-?\d{1,5}.\d{2}) '
        re_string = re_string[:-1]  # Get rid of the trailing space
        re_string = re_string + '\r\n'

    # If you want to see what the above code managed to put together
    # as a regular expression search pattern, uncomment this:
    #print re_string
    #raw_input()

    # Compile this re_string for more efficient re searches
    block_expr = re.compile(re_string)

    # Now get corresponding indices of the
    # variables we need
    full_line = ''
    for r in block_lines:
        full_line = full_line + r[:-2] + ' '
    # Now split it
    varlist = re.split('[ /]',full_line)
    # Get rid of trailing space
    varlist = varlist[:-1]
    # To see the variable list, uncomment
    #print varlist
    #raw_input()

    # Make a dictionary for each variable
    #variables = {}
    #variables[var] = [] for var in varlist
    # Pressure levels to interpolate to
    interp_res = 5 
    plevs = range(200,1050,interp_res)
    
    def min_gt(seq,val):
         # Returns the smallest value greater than some other value
         seq.sort()
         for v in seq:
             if v > val:
                 return v
         return None
    def max_lt(seq,val):
         seq.sort()
         seq.reverse()
         for v in seq:
             if v < val:
                 return v
         return None


    # We now need to break everything up into a chunk for each
    # forecast date and time
    with open(file) as infile:
        blocks = infile.read().split('STID')
        for block in blocks:
            cur_plevs = []
            interp_plevs = []
            header = block
            #print 'HEADER:', header
            
            if header.split()[0] != '=':
                continue
            fcst_time = re.search('TIME = (\d{6}/\d{4})', header).groups()[0]
            fcst_dt = datetime.strptime(fcst_time,'%y%m%d/%H%M')
            #print "FCST:", fcst_dt
            #raw_input()
            temp_vars  = {}
            for var in varlist:
                temp_vars[var] = []
            for block_match in block_expr.finditer(block):
                 vals = block_match.groups()
                 #print "vals:", vals
                 #raw_input()
                 for val, name in zip(vals,varlist):
                      temp_vars[name].append(float(val))
                      #if name == 'PSFC':
                      #    print temp_vars[name]

            #print temp_vars['PRES']
            # Unfortunately, bufkit values aren't always uniformly distributed.
            # Have to interpolate
            cur_plevs = temp_vars['PRES']
            cur_plevs.reverse()
            for var in varlist[1:]:
                values = temp_vars[var]
                values.reverse()
                #print cur_plevs
                #print values
                # Need to re-orient which values we're interpolating to
                if min(plevs) > min(cur_plevs):
                    min_plev = min(plevs)
                else:
                    min_plev = min_gt(plevs,min(cur_plevs))
              
                if max(plevs) < max(cur_plevs):
                    max_plev = max(plevs)
                else:
                    max_plev = max_lt(plevs,max(cur_plevs))

           
                interp_plevs = range(min_plev,max_plev+interp_res, interp_res)

                f = interpolate.interp1d(cur_plevs,values)
                #print "cur_plevs:", cur_plevs
                #print "interp_plevs:", interp_plevs
   
                interp_vals = f(interp_plevs)

                interp_vals = list(interp_vals)
                interp_plevs = list(interp_plevs)
                interp_vals.reverse()
                interp_plevs.reverse()
                temp_vars[var] = interp_vals
            temp_vars['PRES'] = interp_plevs
            profile_dir[fcst_dt] = temp_vars

    return profile_dir




def obs_parser(siteid, starttime = 0, endtime = 0):
    # Script to parse out an observation list from our
    # tdd/td program
    # Start time and end time default to 0.  If not 0, then
    # only parse observations from the time period specified
    # and compute max, min and precip from that.  The start and
    # end times should be specified as datetime objects

    from profile_class import Profile
    import re,os
    from datetime import datetime,timedelta
    from math import sin, cos,atan
    from numpy import around
    # Be able to find elevation
    lat_pat = re.compile('Lat: (-?\d{1,3}.\d{1,4})')   
    lon_pat = re.compile('Lon: (-?\d{1,3}.\d{1,4})')
    elev_pat = re.compile('Elev: (\d{1,4}) m')

    # Know column headers
    cols = ['DY','STA','TP','HHMM','N','CIG','VSBY','WEA','SLP','T','TD','DIR','SP','GS','PK','ALT','TEND','1H','3H6H','24H','SNW','Tmx','Tmn']

    # We know we're looking for lines that match a certain pattern
    re_pat = "(\d{2}) (\D{4})  (\D{2}) (\d{4}) (\d| |\D)(.{4}| \D{3}) (   \d|  \d{1,2}|\d{1}.\d{2})    (\D{1,3}) ( \d{3}.\d{1}|\d{4}.\d{1}) (-\d{2}| \d{2}| -\d{1}|  \d{1})"
    re_pat = re_pat + " (-\d{2}| \d{2}| -\d{1}|  \d{1}) ( \d{2}|\d{3}|  \d{1}|\D{3}) ( \d|\d{2}) (\d{2}|  ) (\d{2}|  ) (\d{2}.\d{2}) (\d{1} \d{2}|    ) (\d{2}| \d| \D)"
    re_pat = re_pat + "  (\d{2}| \D| \d)   (\d{2}| \D| \d)  (\d{2}| \D| \d) (-\d{1,2}| \d{1,2}|  ) (-\d{1,2}| \d{1,2}|  ) (AO2 .+)?"
    block_expr = re.compile(re_pat)    


    # Another compiled regex to find max/min temps
    xn_pat = "Mx=(-?\d{1,3}) Mn=(-?\d{1,3})"
    xn_expr = re.compile(xn_pat)

    # Load the data using tdd or td
    # Depending on if a time was specified
    if starttime:
        p = os.popen('echo "%s\n%s\n\n\n" | /usr/local/wx/bin/td -no10 -q %s' % (starttime.strftime('%H %d %m %Y'),endtime.strftime('%H %d %m %Y'),siteid.lower()))
    else:
        p = os.popen('/usr/local/wx/bin/tdd -no10 %s' % siteid)
    #indata = p.read()
    profile_dir = {}
    maxmin_dir = {}
    cur_profile = Profile()
    #print indata
    #raw_input()
    with p as bufr:
        br = bufr.read()
        for elev_match in elev_pat.finditer(br):
            #print elev_match.groups()
            elevation = float(elev_match.groups()[0])
        for lon_match in lon_pat.finditer(br):
            longitude = float(lon_match.groups()[0])
        for lat_match in lat_pat.finditer(br):
            latitude = float(lat_match.groups()[0])
        curmonth = True
    
    if starttime:
        # If we are looking for longer than 6 days, use the long-term archive
        duration = endtime - starttime
        if duration.days >= 6:
            command = 'echo "%s\n%s\n\n\n" | /usr/local/wx/bin/td -d %s -q %s' % (starttime.strftime('%H %d %m %Y'),endtime.strftime('%H %d %m %Y'),'/home/disk/archive/saous/decoded,YYYYMM/YYYYMMDDHH.saous', siteid.lower())
        else:
            command = 'echo "%s\n%s\n\n\n" | /usr/local/wx/bin/td -q %s' % (starttime.strftime('%H %d %m %Y'),endtime.strftime('%H %d %m %Y'), siteid.lower())
            
        p = os.popen(command)
    else:
        p = os.popen('/usr/local/wx/bin/tdd %s' % siteid)
    with p as bufr:
        for block_match in block_expr.finditer(bufr.read()):
            #print block_match.groups()
            #raw_input()
            obline = block_match.groups()
            del cur_profile
            cur_profile = Profile()

            # Time gets tricky
            nowtime = datetime.utcnow()
            day = int(obline[cols.index('DY')])

            # Go back in time until the day matches since the
            # output of td or tdd doesn't give you the month
            curday = nowtime
            while curday.day != day:
                curday = curday - timedelta(hours=24)
            year = curday.year
            month = curday.month
            hm = obline[cols.index('HHMM')]

            # These are official times, round to nearest hour
            # But let datetime handle the dirty work
            hour = int(hm[0:2])
            minute = int(hm[2:])
            minute = 0
            validtime = datetime(year,month,day,hour,minute)
            validtime = validtime + timedelta(hours=1)
            #print obline            
            # Now assign the appropriate values
            tmpc = (float(obline[cols.index('T')])-32)* 5. / 9.
            dwpc = (float(obline[cols.index('TD')])-32)* 5. / 9.
            cur_profile.tmpc = tmpc
            cur_profile.dwpc = dwpc
            cur_profile.sknt = float(obline[cols.index('SP')])
            try:
                cur_profile.cfrl = int(float(obline[cols.index('N')])/8.0*100)
            except:
                cur_profile.cfrl = 0
            try:
                cur_profile.drct = float(obline[cols.index('DIR')])
            except:
                cur_profile.drct = 0.0
            # Get u and v winds
            rad = 4.0 * atan(1.0) / 180.
            cur_profile.uwnd = -1 * cur_profile.sknt * sin(rad * cur_profile.drct)
            cur_profile.vwnd = -1 * cur_profile.sknt * cos(rad * cur_profile.drct)
            cur_profile.press = float(obline[cols.index('SLP')])
            cur_profile.alt = float(obline[cols.index('ALT')]) * 33.86  # In mb
            cur_profile.elev = elevation
            cur_profile.lat = latitude
            cur_profile.lon = longitude
            maxT = obline[cols.index('Tmx')]
            minT = obline[cols.index('Tmn')]
            if maxT not in ['', '  ',' ']:
                # If the 6-hour max min temp fields are
                # not blank, store them
                cur_profile.maxT = float(maxT)
                cur_profile.minT = float(minT)
                #print "MAX", maxT
                #print "MIN", minT


            # Trace of precip = 0
            if obline[cols.index('1H')] in ['T',' T',' ','  ','']:
                cur_profile.p01m = 0.0 
            else:
                #print obline[cols.index('1H')]
                # Only record precip if the "weather" value is something non-empty (LEM 9/27/2012)
                # Amended to look for the "RAB" in the comment section
                if obline[cols.index('WEA')].strip() not in ['','H','F']: 
                    cur_profile.p01m = float(obline[cols.index('1H')])
                elif obline[-1] != None:
                    if re.search('RAB',obline[-1]) or re.search('RAE',obline[-1]):
                        cur_profile.p01m = float(obline[cols.index('1H')])
                    else:
                        cur_profile.p01m = 0.0

                else:
                    # Assume this was a bad bucket tip
                    cur_profile.p01m = 0.0
            profile_dir[validtime] = cur_profile    
 
            # Find the extreme values for temp only if start and end times have not been given
            if not starttime:
                if obline[-1] is not None:
                    if xn_expr.search(obline[-1]):
                        maxT,minT = xn_expr.search(obline[-1]).groups()
                        #print "Max Min found!", maxT, minT
                        # The date will be for the previous day
                        xntime = validtime - timedelta(hours=24)
                        xntime = xntime.replace(hour=0,minute=0)
                        maxmin_dir[xntime.date()] = [int(maxT),int(minT)]
             
            #print validtime.strftime('%Y%m%d%H'), cur_profile.tmpc, cur_profile.sknt, cur_profile.p01m
    # Before returning the dictionaries, if a starttime and endtime where specified then
    # compute the max and min temperatures and total precip over the period to make the
    # max/min dictionary
    if starttime and len(maxmin_dir.keys()) == 0:
        print "Computing max and min from time series values"
        # Get a sorted list of keys
        sorted = profile_dir.keys()
        sorted.sort()
        minT = 150.
        maxT = -150.
        maxwind = 0.0
        precip = 0.0
        first_record = True
        # Loop through each date 
        for record in sorted:
            if first_record:
                # Skip the first record as it includes data
                # from the previous period
                first_record = False
                # But could contain wind
                if profile_dir[record].sknt > maxwind:
                    maxwind = profile_dir[record].sknt 
            else:
                # Start accumulating precip
                precip = precip + profile_dir[record].p01m
                # Check wind
                if profile_dir[record].sknt > maxwind:
                    maxwind = profile_dir[record].sknt 
                # Check temps
                if profile_dir[record].maxT != '':
                    if profile_dir[record].maxT > maxT:
                        maxT = profile_dir[record].maxT
                    if profile_dir[record].minT < minT:
                        minT = profile_dir[record].minT
                # Now check to see if any intermediate temperature values
                # were higher or lower (if the time period given is less than a day, then
                # perhaps we are between 6-hour writeouts)
                if around(profile_dir[record].tmpc * 9./5. + 32.) > maxT:
                    maxT = around(profile_dir[record].tmpc * 9./5. + 32)
                if around(profile_dir[record].tmpc * 9./5. + 32.) < minT:
                    minT = around(profile_dir[record].tmpc * 9./5. + 32)


        # Now make the maxmin_dir, complete with the 
        # date of starttime.  We prefer to reserve the
        # fourth position for climo report winds if possible,
        # so make the metar max winds the fifth field
        #print "HERE"
        #print maxT, minT, maxwind, precip/100.
        maxmin_dir[starttime.date()] = [maxT,minT,'',maxwind,precip/100.]    

    return profile_dir,maxmin_dir


def mesonet_obs_parser(siteid, starttime = 0, endtime = 0):
    # Script to parse out an observation list
    # Oklahoma Mesonet mts files
    # Start time and end time default to 0.  If not 0, then
    # only parse observations from the time period specified
    # and compute max, min and precip from that.  The start and
    # end times should be specified as datetime objects

    from profile_class import Profile
    import re,os
    from datetime import datetime,timedelta
    from math import sin, cos,atan
    import urllib2
    from numpy import exp,log,around


    def RH_to_Td(T,RH):
        """ Convert T (in Kelvin) and RH to Td (in Kelvin) """
        if T > 273:
            es = 611.0 * exp(5423. * (1/273. - 1/T))
            e = RH * es
            Td = ((1/273.) - (461./2.5E6) * log(e/611.)) ** (-1)
        else:
            es = 611.0 * exp(6139. * (1/273. - 1/T))
            e = RH * es
            Td = ((1/273.) - (461./2.83E6) * log(e/611.)) ** (-1)
        return Td


    if starttime != 0 and endtime != 0:
        specified_time = True
    else:
        specified_time = False
    # Load the mesonet data files
    if endtime == 0:
        endtimefull = datetime.utcnow()
    else:
        endtimefull = endtime
    endmin = endtimefull.minute
    while (endmin % 5) != 0:
        endmin -= 1
    endtime = datetime(endtimefull.year,endtimefull.month,endtimefull.day,endtimefull.hour,endmin)

    if starttime == 0:
        starttime = endtime - timedelta(hours=48)

    # Find the files to load
    indata = ''
    nowdate = datetime(starttime.year,starttime.month,starttime.day)
    nowdates = []
    while nowdate <= datetime(endtime.year,endtime.month,endtime.day):
        nowdates.append(nowdate)
        nowdate += timedelta(hours=24)
    profile_dir = {}
    maxmin_dir = {}
    cur_profile = Profile()
    for nowdate in nowdates:
        webdir = 'http://www.mesonet.org/data/public/mesonet/mts/%d/%02d/%02d/%s%s.mts' %\
            (nowdate.year,nowdate.month,nowdate.day,nowdate.strftime('%Y%m%d'),siteid.lower())
        print nowdate, webdir
        # Convert to UTC
        #nowdate = nowdate + timedelta(hours=6)
        req = urllib2.Request(webdir)
        response = urllib2.urlopen(req)
        indata = response.readlines()
        for line in indata:
            if line.startswith(' STID'):
                cols = line.split()
                continue
            elif not line.startswith(' %s' % siteid.upper()):
                continue
            del cur_profile
            cur_profile = Profile()
            obline = line.split()
            # Get the ob's time in minutes
            min_since = int(obline[cols.index('TIME')])
            validtime = nowdate + timedelta(minutes=min_since)
            if validtime >= (endtimefull-timedelta(minutes=10)):
                continue
            if validtime <= starttime:
                continue
            #print obline            
            # Now assign the appropriate values
            #tmpc = (float(obline[cols.index('TAIR')])-32)* 5. / 9.
            tmpc = float(obline[cols.index('TAIR')])
    
            #dwpc = (float(obline[cols.index('TD')])-32)* 5. / 9.
            rh = (float(obline[cols.index('RELH')]) * 0.01)
            dwpK = RH_to_Td(tmpc+273,rh)
            dwpc = dwpK - 273
            cur_profile.tmpc = tmpc
            cur_profile.dwpc = dwpc
            cur_profile.sknt = float(obline[cols.index('WSPD')]) *1.94384 
            try:
                cur_profile.drct = float(obline[cols.index('WDIR')])
            except:
                cur_profile.drct = 0.0
            # Get u and v winds
            rad = 4.0 * atan(1.0) / 180.
            cur_profile.uwnd = -1 * cur_profile.sknt * sin(rad * cur_profile.drct)
            cur_profile.vwnd = -1 * cur_profile.sknt * cos(rad * cur_profile.drct)

            #cur_profile.press = float(obline[cols.index('PRES')])
            #if maxT not in ['', '  ',' ']:
            #    # If the 6-hour max min temp fields are
            #    # not blank, store them
            #    cur_profile.maxT = float(maxT)
            #    cur_profile.minT = float(minT)
            #    #print "MAX", maxT
            #    #print "MIN", minT


            # Trace of precip = 0
            if obline[cols.index('RAIN')] in ['T',' T',' ','  ','']:
                cur_profile.p01m = 0.0 
            else:
                #print obline[cols.index('1H')]
                # Only record precip if the "weather" value is something non-empty (LEM 9/27/2012)
                # Amended to look for the "RAB" in the comment section
                cur_profile.p01m = float(obline[cols.index('RAIN')]) * 0.0393701
            profile_dir[validtime] = cur_profile    
 
    # Now that we have everything, sort out the maxes and mins
    # Find the max and min of each period
    pkeys = profile_dir.keys()
    pkeys.sort()
    #adjusted_keys = [p+timedelta(hours=1) for p in pkeys]
    #adjusted_keys = pkeys[13:]
    #pkeys = pkeys[13:]
    if specified_time:
        # Find the max and min of time period specified
            # Find the max and min for this day
            temps = [profile_dir[t].tmpc for t in pkeys]
            # Remove bad values
            newtemps = [t for t in temps if t != -999.0]
            temps = newtemps
            # Windspeeds
            wspds = [profile_dir[t].sknt for t in pkeys]
            # Precip values 
            precips = [profile_dir[t].p01m for t in pkeys]
            # First precip is actually the end of the day before
            precips = precips[1:]
            maxT = max(temps) * 9./5. + 32
            minT = min(temps) * 9./5. + 32
            #print "MAX: C:", max(temps), "  F:",maxT,"  Fint:", around(maxT) 
            #print "MIN: C:", min(temps), "  F:",minT,"  Fint:", around(minT) 
            minT = around(minT)
            maxT = around(maxT)
            maxwind = max(wspds)
            # Precip is a bit tricky
            prevmax = 0.
            for pdex in range(1,len(precips)):
                if precips[pdex] < precips[pdex-1]:
                    # Must have hit the end of the day.  Add previous value
                    prevmax = precips[pdex-1]
            precip = precips[-1] + prevmax
            precip = precip - precips[0]
            maxmin_dir[nowdates[-1].date()] = [int(maxT), int(minT),'',maxwind,precip]
    
    else:
        # Find the max and min of each period
        pkeys = profile_dir.keys()
        pkeys.sort()
        #adjusted_keys = [p-timedelta(hours=5) for p in pkeys]
        adjusted_keys = pkeys
        for nowdate in nowdates:
            # Find the max and min for this day
            temps = [profile_dir[t].tmpc for t in pkeys if adjusted_keys[pkeys.index(t)].date() == nowdate.date()]
            # Windspeeds
            wspds = [profile_dir[t].sknt for t in pkeys if adjusted_keys[pkeys.index(t)].date() == nowdate.date()]
            # Precip values 
            precips = [profile_dir[t].p01m for t in pkeys if adjusted_keys[pkeys.index(t)].date() == nowdate.date()]
            precips = precips[1:]
            maxT = max(temps) * 9./5. + 32
            minT = min(temps) * 9./5. + 32
            minT = around(minT)
            maxT = around(maxT)
            maxwind = max(wspds)
            # Precip is a bit tricky
            prevmax = 0.
            for pdex in range(1,len(precips)):
                if precips[pdex] < precips[pdex-1]:
                    # Must have hit the end of the day.  Add previous value
                    prevmax = precips[pdex-1]
            precip = precips[-1] + prevmax
            #print maxT, minT, maxwind, precip
            maxmin_dir[nowdate.date()] = [int(maxT), int(minT),'',maxwind,precip]


    return profile_dir,maxmin_dir


def get_wxchall_verification(siteid,wfoid,dateval):
    # Uses the above parsers to return a one-member dictionary of all the values
    # needed to verify the wxchallenge from 6Z to 6Z for the date specified

    # Convert to the appropriate time
    from datetime import date,datetime,time,timedelta
    tempval = datetime.combine(dateval,time())
    starttime = tempval.replace(hour=6)
    endtime = starttime + timedelta(hours=24)

    verify_dict = {}
    # First, get the max Temp, min Temp, and total precip from td
    profdir, maxmindir = obs_parser(siteid, starttime, endtime)
    for key in maxmindir.keys():
        # Should only be one key, so not a problem
        value_list = maxmindir[key]
    #print "HERE"
    #print value_list
    # Now for the winds from climo parser (if it exists)
    winddict = climo_parser(siteid,wfoid) 
    # Now we have to find the right date
    if starttime.date() in winddict.keys():
        wspd = winddict[starttime.date()].sknt
        value_list[2] = wspd
    else:
        print "No wind value found in climo file!  Date not in record."
        value_list[2] = ''

    # Return the dictionary
    verify_dict[starttime.date()] = value_list
    return verify_dict

def get_todays_verification(siteid,wfoid):
    # Uses the above parsers to return a one-member dictionary of all the values
    # needed to verify the wxchallenge from 6Z this morning through now
    # Won't need to pull the climo file as it hasn't been written yet 

    # Convert to the appropriate time
    from datetime import date,datetime,time,timedelta
    tempval = datetime.utcnow()
    if tempval.hour < 6:
        # We're actually on the next day, so subtract 1 day
        starttime = datetime.strptime(((tempval-timedelta(hours=24)).replace(hour=6)).strftime('%Y%m%d%H'),'%Y%m%d%H')
    else:
        starttime = datetime(tempval.year,tempval.month,tempval.day,6)
    
    #endtime = starttime + timedelta(hours=24)
    #endtime = tempval - timedelta(hours=1)
    endtime = datetime.strptime(tempval.strftime('%Y%m%d%H'),'%Y%m%d%H')

    #print starttime
    #print endtime
    #raw_input()
    verify_dict = {}
    # First, get the max Temp, min Temp, and total precip from td
    if siteid.lower() == 'koun':
        # Endtime can be as recent as 10 minutes ago
        
        while (tempval.minute % 5) != 0:
            tempval -= timedelta(minutes=1)
        endtime = datetime(tempval.year,tempval.month,tempval.day,tempval.hour,tempval.minute)
        #print endtime
        starttime = starttime + timedelta(hours=1)
        profdir, maxmindir = mesonet_obs_parser('nrmn', starttime, endtime)
    else:
        profdir, maxmindir = obs_parser(siteid, starttime, endtime)
    #profdirkeys = profdir.keys()
    #profdirkeys.sort()
    #print "Start:", profdirkeys[0], "End:", profdirkeys[-1]
    for key in maxmindir.keys():
        # Should only be one key, so not a problem
        value_list = maxmindir[key]
    
    #print "HERE"
    #print value_list
    # Now for the winds from climo parser (if it exists)
    winddict = climo_parser(siteid,wfoid) 
    # Now we have to find the right date
    if starttime.date() in winddict.keys():
        wspd = winddict[starttime.date()].sknt
    else:
        print "No wind value found in climo file!  Date not in record."
        value_list[2] = ''
    if value_list[3] != '': 
        # Unfortunately, need to re-convert wind speed back to 
        # mph in this case
        #value_list[3] = value_list[3] * 1.15077945
        pass
    # Return the dictionary
    verify_dict[starttime.date()] = value_list
    return verify_dict


