#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
import plot_routines as pr
import hirs_footprint as hf
import math
import read_functions as rf
import matplotlib.pyplot as plt

#Identify pixels across track that best match HIRS pixel locations.
#Andy's locations
pix_mask = [21,28,35,41,48,55,61,68,75,81,88,95,101,108,115,121,128,135,141,148,154,161,168,174,181,188,194,201,208,214,221,228,234,241,248,254,261,268,274,281,288,294,301,308,314,321,328,334,341,348,354,361,367,374,381,387]

pix_mask_prev = list(np.asarray(pix_mask)-1)
pix_mask_next = list(np.asarray(pix_mask)+2)
#Modificantions to scanline number to track HIRS data.
#Andy's locations
scanline = [12,11,11,11,11,11,10,10,10,10,10,9,9,9,9,9,9,8,8,8,8,8,7,7,7,7,7,6,6,6,6,5,5,5,5,5,4,4,4,4,4,3,3,3,3,3,2,2,2,2,2,1,1,1,1,1,0]


#Allocate some arrays to contain the data
extract_lat = np.zeros([56])
extract_lon = np.zeros([56])
extract_flag = np.zeros([56])
cloud_mask_min = np.zeros([56])
cloud_mask_mean = np.zeros([56])
bt_std_cloud = np.zeros([56])
cloud_frac = np.zeros([56])
extract_dy = np.zeros([56])
extract_n = np.zeros([56])

def extract_3_3(input_array,j,baseline):

    """
    Will extract a 3x3 array from input data centred on the input coordinates.
    Input: input_array: array from which to extract data.
           j: hirs_spot
           baseline: array of closest GAC scanlines for each spot in HIRS scanline
    Output: output_array: 3x3 array extracted from input data, centred on closest GAC pixel
                          to HIRS spot
    """
    
    output_array = np.zeros([3,3])

    output_array[0,:] \
        = input_array[int(baseline[j])-1,pix_mask_prev[j]:pix_mask_next[j]]
    output_array[1,:] \
        = input_array[int(baseline[j]),pix_mask_prev[j]:pix_mask_next[j]]
    output_array[2,:] \
        = input_array[int(baseline[j])+1,pix_mask_prev[j]:pix_mask_next[j]]

    return output_array

def extract_5_3(input_array,j,baseline):

    """
    Will extract a 5x3 array from input data centred on the input coordinates.
    Input: input_array: array from which to extract data.
           j: hirs_spot
           baseline: array of closest GAC scanlines for each spot in HIRS scanline
    Output: output_array: 5x3 array extracted from input data, centred on closest GAC pixel
                          to HIRS spot
    """
 
    output_array = np.zeros([5,3])

    output_array[0,:] \
        = input_array[int(baseline[j])-2,pix_mask_prev[j]:pix_mask_next[j]]
    output_array[1,:] \
        = input_array[int(baseline[j])-1,pix_mask_prev[j]:pix_mask_next[j]]
    output_array[2,:] \
        = input_array[int(baseline[j]),pix_mask_prev[j]:pix_mask_next[j]]
    output_array[3,:] \
        = input_array[int(baseline[j])+1,pix_mask_prev[j]:pix_mask_next[j]]
    output_array[4,:] \
        = input_array[int(baseline[j])+2,pix_mask_prev[j]:pix_mask_next[j]]

    return output_array


def concatenate(stack,array,count):

    """
    Will concatenate a one dimensional array.  Will take empty list as input for count = 0
    and reassign.  For count > 0, will concatenate along first dimension.
    Input: stack: array to add data too
           array: data to add
           count: number of times concatenate has been called for this output
    Output: Concatenated array
    """

    if count == 0:
        stack = np.copy(array)
    elif count == 1:
        stack = np.concatenate([[stack],[array]])
    else:
        array = np.reshape(array,[1,array.shape[0]])
        stack = np.concatenate([stack,array])
    

    return stack

def eval_cloud_mask(cloud_array,bt_array):
        
    """
    Function will calculate the minimum and mean clear-sky probabilities, the cloud fraction
    and mean brightness temperature.
    Inputs: cloud_array: array of clear-sky probabilities
            bt_array: correpsonding array of 11 micron BTs
    Output: cloud_mask_min: minimum clear-sky probability
            cloud_mask_mean: mean clear-sky probability
            cloud_frac: cloud fraction in input array using clear-sky probability threshold
                        of 0.0001
            mean: mean 11 micron brightness temperature
            std: standard deviation of cloudy pixels or coldest three pixels where ncloud < 3.
            all_std: standard deviation of all GAC pixels in HIRS footprint
            n: number of valid GAC pixels in HIRS footprint
    """
    
    #Set invalid pclear pixels = nan
    mask = (cloud_array < 0)
    cloud_array[mask] = np.nan
    bt_array[:,mask] = np.nan
    bt_mask = bt_array < 0.0
    bt_array[bt_mask] = np.nan
    #Calculate pclear min and mean.
    cloud_mask_min = np.nanmin(cloud_array)
    cloud_mask_mean = np.nanmean(cloud_array)

    #Extract 11 micron temp:
    if bt_array.shape[0] == 4:
        bt_11 = bt_array[-1,:]
    else:
        bt_11 = bt_array[-2,:]

    all_std = np.zeros([bt_array.shape[0]])
    mean = np.zeros([bt_array.shape[0]])

    #If there are enough valid data calculate cloud fraction and cloud BT std.
    if len(cloud_array[~mask]) == 0:
        cloud_frac = -1
        std = -1
    else:
        mask2 = (cloud_array[~mask] < 0.0001)
        cloud_frac = float(len(cloud_array[~mask][mask2]))/float(len(cloud_array[~mask]))
        std = -1
        if len(cloud_array[~mask][mask2]) >= 3:
            std = np.std(bt_11[~mask][mask2])
        if std == -1:
            if len(bt_11[~mask]) >= 3:
                std = np.std(np.sort(bt_11[~mask])[0:4])
    

    #Calculate BT mean
    for i in np.arange(0,all_std.shape[0],1):
        all_std[i] = np.nanstd(bt_array[i,:])
        mean[i] = np.nanmean(bt_array[i,:])
#        if i == 3:
#            print all_std[i],bt_array[i,:]
    n = len(cloud_array[~mask])

    return cloud_mask_min,cloud_mask_mean,cloud_frac,mean,std,all_std,n


def get_landmask(hirs_flags,gac_flags):

    """
    Evaluate and return the HIRS and AVHRR landmasks.
    Input: hirs_flags: HIRS L2P flags from GBCS code
           gac_flags: GAC L2P flags from GBCS code
    Output: hirs_landmask (=1 for land pixels)
            avhrr_landmask (=1 for land pixels)
    """

    hirs_landmask = hirs_flags&2 != 0

    #Evaluate the L2P flags to extact land mask
    mask = np.where(gac_flags >= 128)
    gac_flags[mask] = gac_flags[mask]-128
    mask = (gac_flags >= 64)
    gac_flags[mask] = gac_flags[mask]-64
    mask = (gac_flags >= 32)
    gac_flags[mask] = gac_flags[mask]-32
    mask = (gac_flags >= 16)
    gac_flags[mask] = gac_flags[mask]-16
    mask = (gac_flags >= 8)
    gac_flags[mask] = gac_flags[mask]-8
    mask = (gac_flags >= 4)
    gac_flags[mask] = gac_flags[mask]-4
    avhrr_landmask = (gac_flags == 2)

    return hirs_landmask,avhrr_landmask


def collocate_gac_hirs(gac_t,hirs_t,gac_min,gac_max,hirs_min,hirs_max,\
                           gac_lat,gac_lon,gac_flags,hirs_flags,gac_prob,avhrr_obs,\
                           avhrr_noise,gac_dy=-1):

    """
    Collocate GAC and HIRS using time arrays and Andy's function for locating GAC 
    pixels in HIRS footprint.
    Inputs: gac_t: gac_scanline_time
            hirs_t: hirs_scanline_time
            gac_min: lower limit for gac overlap extraction
            gac_max: upper limit for gac overlap extraction
            hirs_min: lower limit for HIRS overlap extraction
            hirs_max: upper limit for HIRS overlap extraction
            gac_lat: GAC latitudes
            gac_lon: GAC longitudes
            gac_flags: GAC L2P flags from GBCS processing
            hirs_flags: HIRS L2P flags from GBCS processing
            gac_prob: GAC clear-sky probability
            gac_bt: GAC 11um brightness temperature
            gac_dy: GAC 11um obs-RTM
    Ouptuts: lat_centre: GAC latitude for centre of HIRS pixel
             lon_centre: GAC longitude for centre of HIRS pixel
             flag_centre: L2P flags for closes GAC pixels
             gac_as_hirs_min: GAC minimum pclear in HIRS footprint
             gac_as_hirs_mean: GAC mean pclear in HIRS footprint
             bt_as_hirs_mean: Mean GAC 11um BT in HIRS footprint
             bt_as_hirs_cloud_std: Cloud BT standard deviation or coldest 3 pixels
             bt_as_hirs_all_std: BT standard deviation across all valid GAC pixels
             gac_as_hirs_cloud_frac: GAC cloud fraction over HIRS pixel
             gac_as_hirs_dy: GAC mean obs-RTM for HIRS pixel
             gac_n: Number of valid GAC pixels within HIRS footprint
    """
    bt_mean = np.zeros([avhrr_obs.shape[0],56])
    bt_std_all = np.zeros([avhrr_obs.shape[0],56])
    bt_noise =  np.zeros([avhrr_obs.shape[0],56])

    #Create space for output data
    lat_centre = []
    lon_centre = []
    flag_centre = []
    gac_as_hirs_min = []
    gac_as_hirs_mean = []
    bt_as_hirs_mean = []
    bt_as_hirs_cloud_std = []
    bt_as_hirs_all_std = []
    gac_as_hirs_cloud_frac = []
    gac_as_hirs_n = []
    gac_as_hirs_noise = []
    #GAC_dy is an optional input argument
    if len(gac_dy) != 1:
         gac_as_hirs_dy = []

    #Extract the relevant parts of the HIRS and GAC orbits where they overlap.
    gac_t = gac_t[gac_min:gac_max,:]
    hirs_t = hirs_t[hirs_min:hirs_max,:]   

    #Extract the landmask
    hirs_landmask,avhrr_landmask = get_landmask(hirs_flags,gac_flags)

    #Only relevant for ocean processing...
    gac_prob[avhrr_landmask] = -1
    avhrr_obs[:,avhrr_landmask] = -1
    
    #HIRS and GAC data both have a single time per scanline.
    for i in np.arange(0,hirs_t.shape[0],1):

        #Identify GAC scanline that closest to HIRS scanline
        diff = gac_t[:,0] - hirs_t[i,0]
        #Slight shift where the gap between identified base scanline is > 12 
        #GAC scanlines.
        if i != 0:
            scanline_diff =  abs(np.argmin(abs(diff))-ref)
            if scanline_diff > 12:
                scanline_use = list(np.asarray(scanline)-1)
            else:
                scanline_use = scanline
        else:
            scanline_use = scanline
        #Extract the data.

        baseline = np.argmin(abs(diff))+scanline_use
        for j in np.arange(0,56,1):
            #Extract lat, lon, flags
            extract_lat[j] = gac_lat[baseline[j],pix_mask[j]]
            extract_lon[j] = gac_lon[baseline[j],pix_mask[j]]
            extract_flag[j] = gac_flags[int(baseline[j]),pix_mask[j]]

            hirs_spot = 55-j
            hirs_time = (hirs_spot*0.1) + gac_t[baseline[-1]][0]
#            dtime = gac_t[baseline[j]][0]-hirs_time
            dtime = hirs_time-gac_t[baseline[j]][0]
            
            use_arr,x0,a,b,nx,ny,ix0 = hf.getgacpix(j,dtime)
            width = int(np.floor(use_arr.shape[1]/2.))
            height = int(np.floor(use_arr.shape[0]/2.))

            #Extract the data that matches the use_arr
            extract_cloud_mask = gac_prob[int(baseline[j])-height:int(baseline[j])+height+1,\
                                              pix_mask[j]-width:pix_mask[j]+width+1]
            extract_bt_mask = avhrr_obs[:,int(baseline[j])-height:int(baseline[j])+height+1,\
                                              pix_mask[j]-width:pix_mask[j]+width+1]

            extract_noise = avhrr_noise[:,int(baseline[j])-height:int(baseline[j])+height+1,\
                                              pix_mask[j]-width:pix_mask[j]+width+1]

            noise_data = np.zeros([avhrr_noise.shape[0]])
            mask_noise = extract_noise < 0.
            extract_noise[mask_noise] = np.nan
            for k in np.arange(0,avhrr_noise.shape[0]):
                squares = extract_noise[k,:]**2
                bt_noise[k,j] = np.sqrt(np.nansum(squares))/len(squares)
                if k == 3:
                    if bt_noise[k,j] > 1:
                        print extract_noise[k,:],bt_noise[k,j]

            if len(gac_dy) != 1:
                ex = gac_dy[int(baseline[j])-height:int(baseline[j])+height+1,\
                                pix_mask[j]-width:pix_mask[j]+width+1]
            
            mask1 = use_arr == 1
            #Calculate 'gac_as_hirs' from extract
            if len(gac_dy) != 1:
                if len(ex[0]) > 0:
                    extract_dy[j] = np.mean(ex[mask1])

            if use_arr.shape == extract_cloud_mask.shape:
                if len(extract_cloud_mask[mask1]) > 0:
                    try: 
                        cloud_mask_min[j],cloud_mask_mean[j],\
                            cloud_frac[j] ,bt_mean[:,j],bt_std_cloud[j],bt_std_all[:,j], \
                            extract_n[j] = eval_cloud_mask(extract_cloud_mask[mask1],\
                                                               extract_bt_mask[:,mask1])
                    except:
                        pass

        #Concatenate output
        lon_centre = concatenate(lon_centre,extract_lon,i)
        lat_centre = concatenate(lat_centre,extract_lat,i)
        flag_centre = concatenate(flag_centre,extract_flag,i)
        gac_as_hirs_min = concatenate(gac_as_hirs_min,cloud_mask_min,i)
        gac_as_hirs_mean = concatenate(gac_as_hirs_mean,cloud_mask_mean,i)
        bt_as_hirs_mean = rf.concatenate_2d(bt_as_hirs_mean,bt_mean,i)
        gac_as_hirs_noise = rf.concatenate_2d(gac_as_hirs_noise,bt_noise,i)
        bt_as_hirs_cloud_std =concatenate(bt_as_hirs_cloud_std,bt_std_cloud,i)
        bt_as_hirs_all_std = rf.concatenate_2d(bt_as_hirs_all_std,bt_std_all,i)
        gac_as_hirs_cloud_frac = concatenate(gac_as_hirs_cloud_frac,cloud_frac,i)
        gac_as_hirs_n = concatenate(gac_as_hirs_n,extract_n,i)
        if len(gac_dy) != 1:
            gac_as_hirs_dy = concatenate(gac_as_hirs_dy,extract_dy,i)
        else:
            gac_as_hirs_dy = -1

        #Reference to calculate the scanline difference.
        ref = np.argmin(abs(diff))


    return lat_centre,lon_centre,flag_centre,gac_as_hirs_min,gac_as_hirs_mean,\
        bt_as_hirs_mean,bt_as_hirs_cloud_std,bt_as_hirs_all_std,gac_as_hirs_cloud_frac,\
        gac_as_hirs_dy,gac_as_hirs_n,gac_as_hirs_noise
