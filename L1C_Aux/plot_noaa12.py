#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import os
import time
from datetime import date,timedelta
import argparse
import collocate_hirs_avhrr as cha

months = ['01','02','03','04','05','06','07','08','09','10','11','12']
days = ['01','02','03','04','05','06','07','08','09','10','11','12', \
            '13','14','15','16','17','18','19','20','21','22','23','24', \
            '25','26','27','28','29','30','31']

#Base directory paths
hirs_dir = '/group_workspaces/cems2/esacci_sst/output/v2.7.1-50-g2effc4b/l1c/NOAA12/'
gac_dir = '/group_workspaces/cems2/esacci_sst/output/CDR2.0/l2p/AVHRR12_G/'


def plot_cloudmasks(gac_as_hirs_min,hirs_prob,mask_gac,mask_hirs,ch7,ch8,ch10,ch11,ch12,gac_as_hirs_mean):


    plt.figure(figsize=[18,10])
    plt.subplot(251)
    plt.imshow(np.fliplr(gac_as_hirs_min),cmap='seismic')
    plt.title('AVHRR GAC PClr (min)')
    plt.colorbar()
    plt.subplot(252)
    plt.imshow(np.fliplr(gac_as_hirs_mean),cmap='seismic')
    plt.title('AVHRR GAC PClr (mean)')
    plt.colorbar()
    plt.subplot(254)
    plt.imshow(np.log10(hirs_prob),cmap='gist_ncar')
    plt.title('HIRS Spec PDF (7,8,11,12)')
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
    plt.imshow(ch12,cmap='bwr',vmin=-5,vmax=5)
    plt.title('CH12 (Obs-Bkg)')
    plt.colorbar()

    plt.tight_layout()
    plt.show()
    
    return


if __name__ in '__main__':

    #Pass in the date to process
    parser = argparse.ArgumentParser()
    parser.add_argument("year")
    parser.add_argument("month")
    parser.add_argument("day")
    args = parser.parse_args()

    in_date = args.year+args.month+args.day
    prev,next = cha.get_dates(in_date)

    string_date_hirs = in_date[0:4]+'/'+in_date[4:6]+'/'+in_date[6:8]+'/'

    for f in os.listdir(hirs_dir+string_date_hirs):
        
        time_hirs = cha.read_time(hirs_dir+string_date_hirs+f)
        
        plt.figure()
        plt.xlabel('Lon')
        plt.ylabel('Lat')

        for date in [prev,in_date,next]:

            string_date = date[0:4]+'/'+date[4:6]+'/'+date[6:8]+'/'

            for g in os.listdir(gac_dir+string_date):
            
                time_gac = cha.read_time(gac_dir+string_date+g)

                #Identify files where the HIRS/AVHRR times overlap
                if np.min(time_gac) > np.max(time_hirs):
                    continue
                elif np.max(time_gac) < np.min(time_hirs):
                    continue
                else:
                    print 'overlap'
                    gac_min,gac_max,hirs_min,hirs_max \
                        = cha.time_mask(time_hirs,time_gac)
                    gac_lat,gac_lon,gac_flags \
                        = cha.refine_geo(gac_dir+string_date+g,gac_min,gac_max)
                    hirs_lat,hirs_lon,hirs_flags \
                        = cha.refine_geo(hirs_dir+string_date_hirs+f,\
                                         hirs_min,hirs_max)
                    gac_prob = cha.read_pclear(gac_dir+string_date+g,gac_min,gac_max)
            
                    ch7,ch8,ch10,ch11,ch12 = cha.read_ffm(hirs_dir+string_date_hirs+f,\
                                                          hirs_min,hirs_max)

                    hirs_prob = cha.read_pclear(hirs_dir+string_date_hirs+f,hirs_min,\
                                                hirs_max)
        
                    lat_centre,lon_centre,flag_centre,gac_as_hirs_min,gac_as_hirs_mean,\
                        = cha.collocate_gac_hirs(time_gac,time_hirs,gac_min,gac_max,\
                                                 hirs_min,hirs_max,gac_lat,\
                                                 gac_lon,gac_flags,gac_prob)
                    
                    mask_gac,mask_hirs \
                        = cha.plot_geoloc(lat_centre,lon_centre,flag_centre,hirs_lat,\
                                          hirs_lon,hirs_flags)

                    plot_cloudmasks(gac_as_hirs_min,hirs_prob,mask_gac,mask_hirs,\
                                        ch7,ch8,ch10,ch11,ch12,gac_as_hirs_mean)

    
