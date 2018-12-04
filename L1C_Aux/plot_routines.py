#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot_cloudmasks(gac_as_hirs_min,hirs_prob,mask_gac,mask_hirs,ch7,ch8,ch10,ch11,ch12,hirs_lim,hirs_five,gac_as_hirs_mean,cloud_frac):

    """ Visualise the cloudmask data and dY arrays. """

    plt.figure(figsize=[18,10])
    plt.subplot(251)
    plt.imshow(np.fliplr(gac_as_hirs_min),cmap='seismic')
    plt.title('AVHRR GAC PClr (min)')
    plt.colorbar()
    plt.subplot(252)
    plt.imshow(np.fliplr(gac_as_hirs_mean),cmap='seismic')
    plt.title('AVHRR GAC PClr (mean)')
    plt.colorbar()
    plt.subplot(253)
    plt.imshow(np.log10(hirs_lim),cmap='gist_ncar')
    plt.title('HIRS Spec PDF (8,10)')
    plt.colorbar()
    plt.subplot(254)
    plt.imshow(np.log10(hirs_five),cmap='gist_ncar')
    plt.title('HIRS Spec PDF (7,8,10,11,12)')
    plt.colorbar()
    plt.subplot(255)
    plt.imshow(np.log10(hirs_prob),cmap='gist_ncar')
    plt.title('HIRS Spec PDF (7,8,10,11)')
    plt.colorbar()

    plt.subplot(256)
    plt.imshow(ch7,cmap='bwr',vmin=-5,vmax=5)
    plt.title('CH7 (Obs-Bkg)')
    plt.colorbar()
    plt.subplot(257)
    plt.imshow(ch8,cmap='bwr',vmin=-5,vmax=5)
    plt.title('CH8 (Obs-Bkg)')
    plt.colorbar()
    plt.subplot(258)
    plt.imshow(ch10,cmap='bwr',vmin=-5,vmax=5)
    plt.title('CH10 (Obs-Bkg)')
    plt.colorbar()
    plt.subplot(259)
    plt.imshow(ch11,cmap='bwr',vmin=-5,vmax=5)
    plt.title('CH11 (Obs-Bkg)')
    plt.colorbar()
    plt.subplot(2,5,10)
#    plt.imshow(ch12,cmap='bwr',vmin=-5,vmax=5)
#    plt.title('CH12 (Obs-Bkg)')
    plt.imshow(np.fliplr(cloud_frac),vmin=0,cmap='jet')
    plt.title('Cloud Fraction')

    plt.colorbar()


    plt.tight_layout()
    plt.show()
    
    return


def extract_landmask(flag_centre,hirs_flags):

    mask_hirs = hirs_flags&2 != 0

    mask = np.where(flag_centre >= 128)
    flag_centre[mask] = flag_centre[mask]-128
    mask = (flag_centre >= 64)
    flag_centre[mask] = flag_centre[mask]-64
    mask = (flag_centre >= 32)
    flag_centre[mask] = flag_centre[mask]-32
    mask = (flag_centre >= 16)
    flag_centre[mask] = flag_centre[mask]-16
    mask = (flag_centre >= 8)
    flag_centre[mask] = flag_centre[mask]-8
    mask = (flag_centre >= 4)
    flag_centre[mask] = flag_centre[mask]-4
    mask_gac = (flag_centre == 2)

    return mask_hirs,mask_gac


def plot_geoloc(lat_centre,lon_centre,flag_centre,hirs_lat,hirs_lon,hirs_flags):

    """ Plot the geolocation information for the HIRS pixels and nearest
    AVHRR GAC pixel. """

    mask_hirs,mask_gac = extract_landmask(flag_centre,hirs_flags)
    
    plt.plot(lon_centre[~mask_gac],lat_centre[~mask_gac],'o',\
                 fillstyle=u'none',ms=0.1,mec='k',label='GAC')
    plt.plot(hirs_lon[~mask_hirs],hirs_lat[~mask_hirs],'o',\
                 fillstyle=u'none',ms=0.1,mec='r',label='HIRS')
    plt.legend()
    plt.show()

    return mask_gac,mask_hirs

def plot_cloud_char(cloud_frac,std,bt,dy):

    """ Plot cloud characteristics for input into SST retrievals. 
    """

    plt.figure(figsize=(12,8))

    plt.subplot(141)
    plt.imshow(cloud_frac,vmin=0,cmap='jet')
    plt.title('Cloud Fraction')
    plt.colorbar()
    
    plt.subplot(142)
    plt.imshow(std,cmap='jet',vmin=0)
    plt.title('11 $\mu$m uncertainty (K)')
    plt.colorbar()

    plt.subplot(143)
    plt.imshow(np.fliplr(dy))
    plt.title('11 $\mu$m dY')
    plt.colorbar()

    
    cloud_height = np.zeros([dy.shape[0],dy.shape[1]])
    mask = (dy <= 0.)
    cloud_height[mask] = abs(dy[mask])*(1./6.)
    cloud_height[dy.mask] = -np.nan

    plt.subplot(144)
    plt.imshow(np.fliplr(cloud_height))
    plt.title('Cloud Height (km)')
    plt.colorbar()
    plt.tight_layout()
    plt.show()

    return

def plot_sim_comp(gac_obs,gac_dy,hirs_obs,hirs_dy):

    obs_diff = np.fliplr(gac_obs)-hirs_obs
    gac_rtm = gac_obs-gac_dy
    hirs_rtm = hirs_obs-hirs_dy
    rtm_diff = np.fliplr(gac_rtm)-hirs_rtm

    plt.figure(figsize=(20,8))
    plt.subplot(191)
    plt.imshow(np.fliplr(gac_obs),vmin=210,vmax=300)
    plt.title('GAC 11 $\mu$m obs')
    plt.colorbar()

    plt.subplot(192)
    plt.imshow(hirs_obs,vmin=210,vmax=300)
    plt.title('HIRS 11 $\mu$m obs')
    plt.colorbar()

    plt.subplot(193)
    plt.imshow(np.fliplr(gac_dy),vmin=-60,vmax=20)
    plt.title('GAC 11 $\mu$m dY')
    plt.colorbar()

    plt.subplot(194)
    plt.imshow(hirs_dy,vmin=-60,vmax=20)
    plt.title('HIRS 11 $\mu$m dY')
    plt.colorbar()

    plt.subplot(195)
    plt.imshow(np.fliplr(gac_rtm),vmin=210,vmax=300)
    plt.title('GAC 11 $\mu$m RTM')
    plt.colorbar()

    plt.subplot(196)
    plt.imshow(hirs_rtm,vmin=210,vmax=300)
    plt.title('HIRS 11 $\mu$m RTM')
    plt.colorbar()

    plt.subplot(197)
    plt.imshow(obs_diff,vmin=-3,vmax=3,cmap='bwr')
    plt.title('GAC - HIRS Obs')
    plt.colorbar()

    plt.subplot(198)
    plt.imshow(rtm_diff,vmin=-3,vmax=3,cmap='bwr')
    plt.title('GAC-HIRS RTM')
    plt.colorbar()

    plt.subplot(199)
    plt.imshow(obs_diff-rtm_diff,vmin=-3,vmax=3,cmap='bwr')
    plt.title('Double Difference')
    plt.colorbar()

    plt.tight_layout()

    plt.show()
                 


    return
