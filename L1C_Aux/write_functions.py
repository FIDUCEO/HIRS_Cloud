#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
import os
import time
import matplotlib.pyplot as plt

out_root = '/gws/nopw/j04/fiduceo/Data/CDR/hirs_cloud/'

def get_output_filename(args,fiduceo_fn):

    outfile = fiduceo_fn[30:59]+'-FIDUCEO-L1C-AUX-'+args.hirs_sensor.upper()+'-fv01.0.nc'

    outdir = out_root+args.hirs_sensor.upper()+'/'+fiduceo_fn[30:34]+'/'+fiduceo_fn[34:36] \
        +'/'+fiduceo_fn[36:38]+'/'

    if not os.path.isdir(outdir):
        os.makedirs(outdir,0777)


    out = outdir+outfile

    return out

def write_metadata(obj):

    obj.description = 'L1C auxiliary data for HIRS cloud detection and stratospheric aerosol'\
        +' retrieval.'
    obj.licence = 'This dataset is released for use under CC-BY licence and was developed' \
        + ' in the EC FIDUCEO project "Fidelity and Uncertainty in Climate Data Records from'\
        ' Earth Observations". Grant Agreement 638822.'
    obj.contact = 'Claire Bulgin: c.e.bulgin@reading.ac.uk'
    obj.history = 'Created '+time.ctime(time.time())

    return

def var_metadata(var_rtm,var_lat,var_lon,var_flag,var_hlat,var_hlon,var_hflag,var_a_obs,\
                     var_a_noise,var_cloud_frac,var_cloud_height,var_cloud_std,\
                     var_obs_std,var_n,var_likelihood,var_rtm_land,var_likelihood_land,\
                     var_landmask):

    #HIRS RTM
    var_rtm.units = 'Kelvin'
    var_rtm.long_name = 'RTTOV clear-sky simulations for HIRS channels'
    var_rtm_land.units = 'Kelvin'
    var_rtm_land.long_name = 'RTTOV clear-sky simulations for HIRS channels over land'
    #GAC Lat
    var_lat.units = 'degrees North'
    var_lat.long_name = 'Closest GAC pixel to HIRS footprint latitude'
    #GAC Lon
    var_lon.units = 'degrees East'
    var_lon.long_name = 'Closest GAC pixel to HIRS footprint longitude'
    #GAC L2P Flags
    var_flag.long_name = 'L2P flags for GAC pixel closest to HIRS footprint'
    #HIRS Lat
    var_hlat.units = 'degrees North'
    var_hlat.long_name = 'Latitude'
    #HIRS Lon
    var_hlon.units = 'degrees East'
    var_hlon.long_name = 'Longitude'
    #HIRS L2P Flags
    var_hflag.long_name = 'L2P flags for HIRS footprint'
    #AVHRR Obs
    var_a_obs.long_name = 'AVHRR GAC Observations at HIRS resolution'
    var_a_obs.comment = 'Channel dimension assigned as follows: 4 channels' \
        + ' [0.6,0.8,3.7,11 micron], 5 channels includes 12 micron, 6 channels' \
        + ' also includes 1.6 micron.'
    #AVHRR Noise
    var_a_noise.long_name = 'AVHRR GAC Noise at HIRS resolution'
    var_a_noise.comment = 'Channel dimension assigned as follows: 4 channels' \
        + ' [0.6,0.8,3.7,11 micron], 5 channels includes 12 micron, 6 channels' \
        + ' also includes 1.6 micron.'
    #Cloud fraction
    var_cloud_frac.valid_min ='0'
    var_cloud_frac.valid_max = '1'
    var_cloud_frac.long_name = 'Definitely cloudy fraction calculated from AVHRR GAC data'\
        +' (pclr threshold = 0.1).'
    #Cloud top height
    var_cloud_height.units = 'km'
    var_cloud_height.long_name = 'Cloud top height calcuated using BT difference from'\
        +' surface simulations and the adiabatic lapse rate'
    #Cloud Std
    var_cloud_std.units = 'Kelvin'
    var_cloud_std.long_name = '11 micron standard deviation of AVHRR definitely cloudy pixels '\
        +'(or three coldest pixels) within HIRS footprint'
    #Obs Std
    var_obs_std.long_name = 'Standard deviation of AVHRR GAC observations within '\
        +'HIRS footprint'
    #N
    var_n.long_name = 'Number of valid GAC pixels'
    #HIRS Likelihood
    var_likelihood.long_name = 'Clear-sky likelihood'
    var_likelihood.comment = 'Spectral only clear-sky likelihood, not normalised.'
    var_likelihood_land.long_name = 'Clear-sky likelihood over land'
    var_likelihood_land.comment = 'Spectral only clear-sky likelihood, not normalised.'
    #Landmask
    var_landmask.long_name = 'AVHRR as HIRS landmask'
    var_landmask.comment = 'If any AVHRR pixel in those collocated to HIRS footprint is land'\
        + ' this mask as true.' 


    return var_rtm,var_lat,var_lon,var_flag,var_hlat,var_hlon,var_hflag,var_a_obs,\
        var_a_noise,var_cloud_frac,var_cloud_height,var_cloud_std,\
        var_obs_std,var_n,var_likelihood,var_rtm_land,var_likelihood_land,var_landmask


def create_output(obs,a_obs,time,outfile):

    output = Dataset(outfile,'w')
    write_metadata(output)

    #Create dimensions
    output.createDimension('hirs_channel',obs.shape[0])
    output.createDimension('avhrr_channel',a_obs.shape[0])
    #Create dimensions that are equal to the entire HIRS orbit size
    output.createDimension('y',time.shape[0])
    output.createDimension('x',time.shape[1])
    
    var_rtm = output.createVariable('HIRS_RTM','f',['hirs_channel','y','x'],\
                                        zlib=True,complevel=6,shuffle=True)
    var_rtm_land = output.createVariable('HIRS_RTM_Land','f',['hirs_channel','y','x'],\
                                        zlib=True,complevel=6,shuffle=True)
    var_lat = output.createVariable('AVHRR_Lat','f',['y','x'],\
                                        zlib=True,complevel=6,shuffle=True)
    var_lon = output.createVariable('AVHRR_Lon','f',['y','x'],\
                                        zlib=True,complevel=6,shuffle=True)
    var_flag = output.createVariable('AVHRR_Flag','f',['y','x'],\
                                        zlib=True,complevel=6,shuffle=True)
    var_hlat = output.createVariable('HIRS_Lat','f',['y','x'],\
                                        zlib=True,complevel=6,shuffle=True)
    var_hlon = output.createVariable('HIRS_Lon','f',['y','x'],\
                                        zlib=True,complevel=6,shuffle=True)
    var_hflag = output.createVariable('HIRS_Flag','f',['y','x'],\
                                        zlib=True,complevel=6,shuffle=True)
    var_a_obs = output.createVariable('AVHRR_Observations','f',\
                                          ['avhrr_channel','y','x'],\
                                        zlib=True,complevel=6,shuffle=True)
    var_a_noise = output.createVariable('AVHRR_Noise','f',['avhrr_channel','y','x'],\
                                        zlib=True,complevel=6,shuffle=True)
    var_cloud_frac = output.createVariable('Cloud_Fraction','f',['y','x'],\
                                        zlib=True,complevel=6,shuffle=True)
    var_cloud_height = output.createVariable('Cloud_Height','f',['y','x'],\
                                        zlib=True,complevel=6,shuffle=True)
    var_cloud_std = output.createVariable('Cloudy_STD','f',['y','x'],\
                                        zlib=True,complevel=6,shuffle=True)
    var_obs_std = output.createVariable('AVHRR_Obs_STD','f',['avhrr_channel','y','x'],\
                                        zlib=True,complevel=6,shuffle=True)
    var_n = output.createVariable('N_AVHRR','d',['y','x'],\
                                        zlib=True,complevel=6,shuffle=True)
    var_likelihood = output.createVariable('HIRS_likelihood','f',['y','x'],\
                                        zlib=True,complevel=6,shuffle=True)
    var_likelihood_land = output.createVariable('HIRS_likelihood_land','f',['y','x'],\
                                        zlib=True,complevel=6,shuffle=True)
    var_landmask = output.createVariable('AVHRR_Landmask','d',['y','x'],\
                                        zlib=True,complevel=6,shuffle=True)

    var_rtm,var_lat,var_lon,var_flag,var_hlat,var_hlon,var_hflag,var_a_obs,\
        var_a_noise,var_cloud_frac,var_cloud_height,var_cloud_std,var_obs_std,\
        var_n,var_likelihood,var_rtm_land,var_likelihood_land,var_landmask \
        = var_metadata(var_rtm,var_lat,var_lon,var_flag,var_hlat,var_hlon,\
                           var_hflag,var_a_obs,var_a_noise,var_cloud_frac,\
                           var_cloud_height,var_cloud_std,var_obs_std,var_n,\
                           var_likelihood,var_rtm_land,var_likelihood_land,\
                           var_landmask)


    return output,var_rtm,var_lat,var_lon,var_flag,var_hlat,var_hlon,\
        var_hflag,var_a_obs,var_a_noise,var_cloud_frac,var_cloud_height,\
        var_cloud_std,var_obs_std,var_n,var_likelihood,var_rtm_land,var_likelihood_land,\
        var_landmask

def write_data(obs,dy,lat,lon,flag,hirs_lat,hirs_lon,hirs_flag,a_obs,a_noise,\
                   output,hirs_min,hirs_max,cloud_frac,\
                   cloud_std,obs_std,n,hirs_prob,landmask,var_rtm,\
                   var_lat,var_lon,var_flag,var_hlat,var_hlon,var_hflag,var_a_obs,\
                   var_a_noise,var_cloud_frac,var_cloud_height,var_cloud_std,\
                   var_obs_std,var_n,var_likelihood,var_rtm_land,var_likelihood_land,\
                   var_landmask):


    rtm = np.zeros([obs.shape[0],obs.shape[1],obs.shape[2]])
    rtm[:,:,:] = -1
    mask = (obs > 0)*(obs < 400)
    mask &= (dy > -150)
    mask &= (dy < 20)
    rtm[mask] = obs[mask]-dy[mask]

    a_obs = np.swapaxes(a_obs,0,1)
    a_noise = np.swapaxes(a_noise,0,1)
    obs_std = np.swapaxes(obs_std,0,1)

    hirs_flag = np.float64(hirs_flag)

    for i in np.arange(0,19):
        rtm[i,:,:] = np.fliplr(rtm[i,:,:])

    var_rtm[:,hirs_min:hirs_max,:] = rtm
    var_lat[hirs_min:hirs_max,:] = lat
    var_lon[hirs_min:hirs_max,:] = lon
    var_flag[hirs_min:hirs_max,:] = flag
    var_hlat[hirs_min:hirs_max,:] = np.fliplr(hirs_lat)
    var_hlon[hirs_min:hirs_max,:] = np.fliplr(hirs_lon)
    var_hflag[hirs_min:hirs_max,:] = hirs_flag
    var_a_obs[:,hirs_min:hirs_max,:] = a_obs
    var_a_noise[:,hirs_min:hirs_max,:] = a_noise
    var_cloud_frac[hirs_min:hirs_max,:] = cloud_frac
    var_cloud_height[hirs_min:hirs_max,:] = 0.0
    var_cloud_std[hirs_min:hirs_max,:] = cloud_std
    var_obs_std[:,hirs_min:hirs_max,:] = obs_std
    var_n[hirs_min:hirs_max,:] = n
    var_likelihood[hirs_min:hirs_max,:] = np.fliplr(hirs_prob)
    #Add empty arrays over land
    var_rtm_land[:,:,:] = 0.0
    var_likelihood_land[:,:] = 0.0
    var_landmask[hirs_min:hirs_max,:] = landmask

    return output

