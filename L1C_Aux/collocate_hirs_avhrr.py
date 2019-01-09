#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
import os
import time
from datetime import date,timedelta
import argparse
import read_functions as rf
import plot_routines as pr
import collocation as col
import time_and_date_functions as td
import compare_11_micron as cm
import write_functions as wf
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

#Base directory paths
root = '/group_workspaces/cems2/esacci_sst/'
fiduceo_root = '/gws/nopw/j04/fiduceo/Data/FCDR/HIRS/v0.8pre2_no_harm/easy/'
avhrr_sim_dir = ''

#avhrr_sim_dir = '/group_workspaces/cems2/esacci_sst/output/v2.7.1-50-g2effc4b/l1c/AVHRRMTA_G/'
#avhrr_sim_dir = '/group_workspaces/cems2/esacci_sst/output/v2.7.1-50-g2effc4b/l1c/AVHRR10_G/'

def create_filepaths(args):

    """
    Define the filepaths required for processing from the input arguments.
    """

    hirs_dir = root+'/output/'+args.processing_version+'/l1c/'+args.hirs_sensor.upper()+'/'
    fiduceo_dir = fiduceo_root+args.hirs_sensor.lower()+'/'
    gac_dir = root+'output/CDR2.0/l2p/'+args.avhrr_sensor+'/'
    l1c_dir = root+'input/avhrr/l1c/'+args.avhrr_sensor+'/v01.6/'

    return hirs_dir,fiduceo_dir,gac_dir,l1c_dir

if __name__ in '__main__':

    #Pass in the date to process
    parser = argparse.ArgumentParser()
    parser.add_argument("hirs_sensor")
    parser.add_argument("processing_version")
    parser.add_argument("avhrr_sensor")
    parser.add_argument("year")
    parser.add_argument("month")
    parser.add_argument("day")
    args = parser.parse_args()

    #Determine file paths
    hirs_dir,fiduceo_dir,gac_dir,l1c_dir = create_filepaths(args)

    in_date = args.year+args.month+args.day
    prev,next = td.get_dates(in_date)

    string_date_hirs = in_date[0:4]+'/'+in_date[4:6]+'/'+in_date[6:8]+'/'

    for f in os.listdir(hirs_dir+string_date_hirs):

        fiduceo_filename = rf.construct_fiduceo(fiduceo_dir,string_date_hirs,f)
        time_hirs = rf.read_time(hirs_dir+string_date_hirs+f)
        
        plt.figure()
        plt.xlabel('Lon')
        plt.ylabel('Lat')

        fileobj = -1

        for date in [prev,in_date,next]:

            string_date = date[0:4]+'/'+date[4:6]+'/'+date[6:8]+'/'

            for g in os.listdir(gac_dir+string_date):
            
                time_gac = rf.read_time(gac_dir+string_date+g)

                #Identify files where the HIRS/AVHRR times overlap
                if np.min(time_gac) > np.max(time_hirs):
                    continue
                elif np.max(time_gac) < np.min(time_hirs):
                    continue
                else:
                    print 'overlap'
                    #Identify overlap times for GAC and HIRS
                    gac_min,gac_max,hirs_min,hirs_max \
                        = td.time_mask(time_hirs,time_gac)
                    #Read in for GAC and HIRS: lat,lon and l2 flags
                    gac_lat,gac_lon,gac_flags \
                        = rf.refine_geo(gac_dir+string_date+g,gac_min,gac_max)
                    hirs_lat,hirs_lon,hirs_flags \
                        = rf.refine_geo(hirs_dir+string_date_hirs+f,\
                                         hirs_min,hirs_max)

                    #Read in GAC and HIRS BT
                    hirs_bt,hirs_noise,hirs_obs,hirs_ind,hirs_str,hirs_com\
                        = rf.read_hirs_obs(fiduceo_dir+string_date_hirs \
                                               +fiduceo_filename,hirs_min,hirs_max)
                    l1c_filename = rf.construct_l1c(l1c_dir,string_date,g,gac_dir)
                    print l1c_filename
                    gac_bt,avhrr_obs,avhrr_noise \
                        = rf.read_avhrr_obs(l1c_dir+string_date+l1c_filename,\
                                                gac_min,gac_max,args.hirs_sensor)

                    if os.path.isfile(avhrr_sim_dir+string_date+g):
                        gac_dy = rf.read_dy(avhrr_sim_dir+string_date+g,gac_min,\
                                                gac_max)
                    else:
                        gac_dy=np.zeros([1])

                    #Read in the PClear arrays
                    gac_prob = rf.read_pclear(gac_dir+string_date+g,gac_min,gac_max)
            
                    hirs_prob = rf.read_pclear(hirs_dir+string_date_hirs+f,hirs_min,\
                                                hirs_max)

                    #### Extra probabilities from different runs - move from
                    #### core code??
#                    hirs_prob_lim = rf.read_pclear(comp_dir+string_date_hirs+f,hirs_min,\
#                                                hirs_max)
#                    hirs_prob_five =  rf.read_pclear(all_dir+string_date_hirs+f,hirs_min,\
#                                                hirs_max)

                    #Read in the dY arrays
                    dy = rf.read_ffm(hirs_dir+string_date_hirs+f,hirs_min,hirs_max)

                    #Co-locate the GAC and HIRS arrays
                    lat_centre,lon_centre,flag_centre,gac_as_hirs_min,\
                        gac_as_hirs_mean,gac_bt_mean,gac_bt_cloud_std,gac_bt_all_std,\
                        cloud_frac,gac_as_hirs_dy,gac_as_hirs_n,gac_as_hirs_noise, \
                        gac_as_hirs_landmask \
                        = col.collocate_gac_hirs(time_gac,time_hirs,gac_min,gac_max,\
                                                     hirs_min,hirs_max,gac_lat,gac_lon,\
                                                     gac_flags,hirs_flags,gac_prob,\
                                                     avhrr_obs,avhrr_noise,args.hirs_sensor,gac_dy)


                    #Calculate the cloud height
  #                  try:
  #                      cloud_height = np.zeros([lat_centre.shape[0],lat_centre.shape[1]])
  #                      mask = (dy[7,:,:] <= 0.)
  #                      cloud_height[mask] = abs(dy[7,:,:][mask])*(1./6.)
  #                      try:
  #                          cloud_height[dy[7,:,:].mask] = -np.nan
  #                      except:
  #                          pass
  #                  except:
  #                      pass

                    #Write out data
                    outfile = wf.get_output_filename(args,fiduceo_filename)

                    if fileobj == -1:
                        fileobj,var_rtm,var_lat,var_lon,var_flag,\
                            var_hlat,var_hlon,var_hflag,var_a_obs,var_a_noise, \
                            var_cloud_frac,var_cloud_height,var_cloud_std,\
                            var_obs_std,var_n,var_likelihood,var_rtm_land, \
                            var_likelihood_land,var_landmask \
                            = wf.create_output(hirs_obs,avhrr_obs,\
                                                   time_hirs,outfile)

 
                    fileobj = wf.write_data(hirs_obs,dy,lat_centre,lon_centre,\
                                                flag_centre,hirs_lat,hirs_lon,\
                                                hirs_flags,gac_bt_mean,gac_as_hirs_noise,\
                                                fileobj,hirs_min,hirs_max,cloud_frac,\
                                                gac_bt_cloud_std,\
                                                gac_bt_all_std,gac_as_hirs_n,hirs_prob,\
                                                gac_as_hirs_landmask,var_rtm,var_lat,var_lon,\
                                                var_flag,var_hlat,var_hlon,\
                                                var_hflag,var_a_obs,var_a_noise,\
                                                var_cloud_frac,var_cloud_height,\
                                                var_cloud_std,var_obs_std,var_n,var_likelihood,\
                                                var_rtm_land,var_likelihood_land,var_landmask)
                    

                    #Make some comparison plots

                    mask_gac,mask_hirs \
                        = pr.plot_geoloc(lat_centre,lon_centre,flag_centre,hirs_lat,\
                                          hirs_lon,hirs_flags)
#                    pr.plot_cloud_char(cloud_frac,gac_bt_std,gac_bt_mean,ch8)


#                    pr.plot_sim_comp(gac_bt_mean,gac_as_hirs_dy,hirs_bt,ch8)
#                    cm.explore_bias(hirs_bt,gac_bt_mean,hirs_noise,gac_bt_all_std,\
#                                        gac_as_hirs_n)

  #                  cm.calc_stats(hirs_bt,hirs_noise,gac_bt_mean,gac_as_hirs_mean,\
  #                                    hirs_prob,gac_bt_std,cloud_frac,gac_bt_mean_ext,\
  #                                    gac_as_hirs_mean_ext,gac_bt_std_ext,cloud_frac_ext)

  #                  pr.plot_cloudmasks(gac_as_hirs_min,hirs_prob,mask_gac,mask_hirs,\
  #                                      ch7,ch8,ch10,ch11,ch12,hirs_prob_lim,\
  #                                      hirs_prob_five,gac_as_hirs_mean,cloud_frac)

                    
        #Close output file
        fileobj.close()
