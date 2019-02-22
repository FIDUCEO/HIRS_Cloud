#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
import os
import fnmatch

def read_time(filename):

    """ Read time from a given file. Assumes time is saved as an offset
    and dtime array. """

    data = Dataset(filename)
    time = data.variables['time'][:]
    dtime = data.variables['sst_dtime'][0,:,:]
    data.close()
    time = dtime+time

    return time

def read_geo(filename):
    
    """ Read the geolocation data from a given file. """

    data = Dataset(filename)
    lat = data.variables['lat'][:]
    lon = data.variables['lon'][:]
    data.close()

    return lat,lon

def read_pclear(filename,mini,maxi):

    """ Read the probability of clear-sky from a given file.
    Note that for HIRS this is the clear-sky likelihood (probability density)
    from comparing the observations with the RTTOV simulations. """

    data = Dataset(filename)
    pclr = data.variables['probability_clear'][0]
    pclr = pclr[mini:maxi,:]
    data.close()

    return pclr

def concatenate_2d(stack,array,count):

    if count == 0:
        stack = np.copy(array)
        stack = np.reshape(stack,[1,stack.shape[0],stack.shape[1]])
    else:
        array = np.reshape(array,[1,array.shape[0],array.shape[1]])
        stack = np.concatenate([stack,array])

    return stack

def read_ffm(filename,mini,maxi):

    """Read in the dY arrays (Observation - Background) for a selection
    of channels."""

    dy = []

    data = Dataset(filename)
    dy = concatenate_2d(dy,data.variables['ffm.brightness_temperature_ch1'][0,mini:maxi,:],0)
    dy = concatenate_2d(dy,data.variables['ffm.brightness_temperature_ch2'][0,mini:maxi,:],1)
    dy = concatenate_2d(dy,data.variables['ffm.brightness_temperature_ch3'][0,mini:maxi,:],2)
    dy = concatenate_2d(dy,data.variables['ffm.brightness_temperature_ch4'][0,mini:maxi,:],3)
    dy = concatenate_2d(dy,data.variables['ffm.brightness_temperature_ch5'][0,mini:maxi,:],4)
    dy = concatenate_2d(dy,data.variables['ffm.brightness_temperature_ch6'][0,mini:maxi,:],5)
    dy = concatenate_2d(dy,data.variables['ffm.brightness_temperature_ch7'][0,mini:maxi,:],6)
    dy = concatenate_2d(dy,data.variables['ffm.brightness_temperature_ch8'][0,mini:maxi,:],7)
    dy = concatenate_2d(dy,data.variables['ffm.brightness_temperature_ch9'][0,mini:maxi,:],8)
    dy = concatenate_2d(dy,data.variables['ffm.brightness_temperature_ch10'][0,mini:maxi,:],9)
    dy = concatenate_2d(dy,data.variables['ffm.brightness_temperature_ch11'][0,mini:maxi,:],10)
    dy = concatenate_2d(dy,data.variables['ffm.brightness_temperature_ch12'][0,mini:maxi,:],11)
    dy = concatenate_2d(dy,data.variables['ffm.brightness_temperature_ch13'][0,mini:maxi,:],12)
    dy = concatenate_2d(dy,data.variables['ffm.brightness_temperature_ch14'][0,mini:maxi,:],13)
    dy = concatenate_2d(dy,data.variables['ffm.brightness_temperature_ch15'][0,mini:maxi,:],14)
    dy = concatenate_2d(dy,data.variables['ffm.brightness_temperature_ch16'][0,mini:maxi,:],15)
    dy = concatenate_2d(dy,data.variables['ffm.brightness_temperature_ch17'][0,mini:maxi,:],16)
    dy = concatenate_2d(dy,data.variables['ffm.brightness_temperature_ch18'][0,mini:maxi,:],17)
    dy = concatenate_2d(dy,data.variables['ffm.brightness_temperature_ch19'][0,mini:maxi,:],18)
    data.close()  
   
    return dy

def read_flags(filename):

    """Read in the L2 flags from the CCI processing to identify ocean only. """

    data = Dataset(filename)
    flags = data.variables['l2p_flags'][0]
    data.close()

    return flags

def read_hirs_obs(filename,mini,maxi):

    """ Read in the brightness temperatures from the original FIDUCEO files. """

    data = Dataset(filename)
    data.set_auto_mask(False)
    hirs_bt11 = data.variables['bt'][7,mini:maxi,:]
    hirs_noise = data.variables['u_independent'][mini:maxi,:,7]
    hirs_obs = data.variables['bt'][:,mini:maxi,:]
    hirs_ind = data.variables['u_independent'][mini:maxi,:,:]
    hirs_str = data.variables['u_structured'][mini:maxi,:,:]
    hirs_common = data.variables['u_common'][mini:maxi,:,:]
    data.close

    return hirs_bt11, hirs_noise, hirs_obs, hirs_ind, hirs_str, hirs_common

def read_avhrr_obs(filename,mini,maxi,sensor):

    print filename

    """ Read in the brightness temperatures from the GAC L1C file. """

    avhrr_v16 = -1
    avhrr_bt12 = -1
    
    avhrr_obs = []
    avhrr_noise = []

    data = Dataset(filename)
    avhrr_obs = concatenate_2d(avhrr_obs,data.variables['ch1'][0,mini:maxi,:],0)
    avhrr_obs = concatenate_2d(avhrr_obs,data.variables['ch2'][0,mini:maxi,:],1)
    avhrr_noise = concatenate_2d(avhrr_noise,data.variables['ch1_noise']\
                                     [0,mini:maxi,:],0)
    avhrr_noise = concatenate_2d(avhrr_noise,data.variables['ch2_noise']\
                                     [0,mini:maxi,:],1)

    if sensor.lower() in ['noaa15','noaa16','noaa17','noaa18','noaa19','metopa']:
        avhrr_obs = concatenate_2d(avhrr_obs,data.variables['ch3a'][0,mini:maxi,:],2)
        avhrr_noise = concatenate_2d(avhrr_noise,data.variables['ch3a_noise']\
                                         [0,mini:maxi,:],2)
        i = 3
    else:
        i = 2
    avhrr_obs = concatenate_2d(avhrr_obs,data.variables['ch3b'][0,mini:maxi,:],i)
    avhrr_obs = concatenate_2d(avhrr_obs,data.variables['ch4'][0,mini:maxi,:],i+1)
    avhrr_noise = concatenate_2d(avhrr_noise,data.variables['ch3b_nedt']\
                                     [0,mini:maxi,:],i)
    avhrr_noise = concatenate_2d(avhrr_noise,data.variables['ch4_nedt']\
                                     [0,mini:maxi,:],i+1)

    if sensor.lower() in ['noaa07','noaa09','noaa11','noaa12','noaa14','noaa15',\
                              'noaa16','noaa17','noaa18','noaa19','metopa']:
        avhrr_obs = concatenate_2d(avhrr_obs,data.variables['ch5'][0,mini:maxi,:],i+2)
        avhrr_noise = concatenate_2d(avhrr_noise,data.variables['ch5_nedt']\
                                         [0,mini:maxi,:],i+2)

        
    #Retain original functionality
    avhrr_bt11 = data.variables['ch4'][0,mini:maxi,:]

    data.close()

    return avhrr_bt11,avhrr_obs,avhrr_noise 

def construct_fiduceo(fiduceo_dir,string_date_hirs,infile):

    instrument = infile[28:35]
    if instrument == 'HIRSMTA': 
        instrument = 'METOPA'
    else:
        instrument = instrument[:-1]
    acq_time = infile[0:10]

    #print instrument,acq_time
    #print fiduceo_dir+string_date_hirs

    pattern = 'FIDUCEO_FCDR_L1C_HIRS*_'+instrument+'_'+acq_time+'*.nc'
    #print pattern
    #print 'list of files: ',os.listdir(fiduceo_dir+string_date_hirs)
    filename = fnmatch.filter(os.listdir(fiduceo_dir+string_date_hirs),pattern)[0]

    return filename

def construct_l1c(l1c_dir,string_date_gac,infile,gac_dir):

    data = Dataset(gac_dir+string_date_gac+infile)
    l1c_file = data.source_file
    data.close()

    return l1c_file


def refine_geo(filename,mini,maxi):

    """ Refine the geo, lat and lon arrays according to the overlap
    limits between the AVHRR and HIRS files. """

    lat,lon = read_geo(filename)
    flags = read_flags(filename)
    lat = lat[mini:maxi,:]
    lon = lon[mini:maxi,:]
    flags = flags[mini:maxi,:]

    return lat,lon,flags

def read_dy(filename,mini,maxi):

    """Read the 11 micron dy array for AVHRR data. """

    data = Dataset(filename)
    dy = data.variables['ffm.brightness_temperature_ch4'][0,mini:maxi,:]
    data.close()

    return dy
