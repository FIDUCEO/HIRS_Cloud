#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
import time
from datetime import date,timedelta

months = ['01','02','03','04','05','06','07','08','09','10','11','12']
days = ['01','02','03','04','05','06','07','08','09','10','11','12', \
            '13','14','15','16','17','18','19','20','21','22','23','24', \
            '25','26','27','28','29','30','31']

def get_new_date(in_date,delta):

    """ Take the input date and time delta (in days) and find the requested
    new date """

    t = time.strptime(in_date,'%Y%m%d')
    new_date = date(t.tm_year,t.tm_mon,t.tm_mday)+timedelta(delta)

    return new_date

def get_dates(in_date):

    """ Function determines two new dates - the day previous and the day after
    the input date.  This calls 'get_new_date' and then converts the output
    to a string."""

    p1 = get_new_date(in_date,1)
    p2 = get_new_date(in_date,-1)

    #Convert to string.
    p1 = str(p1.year)+months[p1.month-1]+days[p1.day-1]
    p2 = str(p2.year)+months[p2.month-1]+days[p2.day-1]

    return p1,p2

def time_mask(time_hirs,time_gac):

    """Mask the HIRS and AVHRR GAC arrays to cover a consistent period 
    in time. """

    #There is only one time per scanline so keep the arrays two dimensional.
    hirs_t = time_hirs[:,0]
    gac_t = time_gac[:,0]
    gd_gac = gac_t > np.min(hirs_t)
    gd_gac &= gac_t < np.max(hirs_t)
    minimum = np.argwhere(gac_t == gac_t[gd_gac][0])
    maximum = np.argwhere(gac_t == gac_t[gd_gac][-1])    
    gac_min = int(minimum[0])
    gac_max = int(maximum[0])

    gd_hirs = hirs_t > np.min(gac_t[gd_gac])
    gd_hirs &= hirs_t < np.max(gac_t[gd_gac])
    if len(hirs_t[gd_hirs]) == 0:
        #Add loop to identify out places where times overlap but no data is available
        print "Your loop worked James"
        hirs_min = 0.0
        hirs_max = 0.0
        data_check = False
    else:
        minimum = np.argwhere(hirs_t == hirs_t[gd_hirs][0])
        maximum = np.argwhere(hirs_t == hirs_t[gd_hirs][-1])
        hirs_min = int(minimum[0])
        hirs_max = int(maximum[0])
        data_check = True

    return gac_min,gac_max,hirs_min,hirs_max,data_check
