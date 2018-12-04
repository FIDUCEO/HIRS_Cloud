#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def plot_uncertainties(gac_prob,cloud_frac,met,hirs_prob,colour,string):

    intervals = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.98,1]
    locs = [0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.925,0.965,0.99]
    met_int = np.zeros(len(intervals)-1)
    met_frac = np.zeros(len(intervals)-1)
    hirs_frac = np.zeros(len(intervals)-1)
    percentile = np.zeros(len(intervals)-1)

    for x in np.arange(0,12):

        mask = (gac_prob > intervals[x])*(gac_prob <= intervals[x+1])
        met_int[x] = np.nanmean(met[mask])
        mask = (cloud_frac  > intervals[x])*(cloud_frac <= intervals[x+1])
        met_frac[x] = np.nanmean(met[mask])
        mask &= hirs_prob != -999.
        hirs_frac[x] = np.nanmean(hirs_prob[mask])
 

        mask = hirs_prob != -999.
        lower_lim =  np.percentile(hirs_prob[mask],intervals[x]*100)
        upper_lim =  np.percentile(hirs_prob[mask],intervals[x+1]*100)
        mask &= hirs_prob > lower_lim
        mask &= hirs_prob <= upper_lim
        percentile[x] = np.nanmean(met[mask])
        if lower_lim == upper_lim:
            mask = hirs_prob != -999.
            mask &= hirs_prob == lower_lim
            percentile[x] = np.nanmean(met[mask])

    plt.subplot(221)
    plt.plot(locs,met_int,color=colour,label=string)
    plt.xlabel('Mean GAC Clear-Sky Probability')
    plt.ylabel('|gac-hirs| / sqrt( ((hirs_noise^2+gac_std^2) / 3 ))')
    plt.subplot(222)
    plt.plot(locs,met_frac,color=colour,label=string)
    plt.xlabel('GAC Cloud Fraction')
    plt.ylabel('|gac-hirs| / sqrt( ((hirs_noise^2+gac_std^2) / 3 ))')
    plt.subplot(223)
    plt.plot(locs,hirs_frac,color=colour,label=string)
    plt.xlabel('GAC Cloud Fraction')
    plt.ylabel('Mean HIRS Clear-Sky Likelihood')
    plt.subplot(224)
    plt.plot(locs,percentile,color=colour,label=string)
    plt.xlabel('HIRS likelihood percentile')
    plt.ylabel('|gac-hirs| / sqrt( ((hirs_noise^2+gac_std^2) / 3 ))')
    plt.legend()
    plt.tight_layout()
    
    return


def calc_stats(hirs_bt,hirs_noise,gac_bt,gac_prob,hirs_prob,gac_std,cloud_frac,\
                   gac_bt_ext,gac_prob_ext,gac_std_ext,cloud_frac_ext):

    #Calculate metrics
    mod = abs(np.fliplr(gac_bt)-hirs_bt)
    unc = np.sqrt(((hirs_noise**2)+(np.fliplr(gac_std)**2/3.)))
    plt.figure()
    plt.subplot(141)
    plt.imshow(np.fliplr(gac_bt),vmin=210,vmax=290,cmap='jet')
    plt.colorbar()
    plt.subplot(142)
    plt.imshow(hirs_bt,vmin=210,vmax=290,cmap='jet')
    plt.colorbar()
    plt.subplot(143)
    plt.imshow(np.fliplr(gac_bt)-hirs_bt,vmax=3,vmin=-3,cmap='bwr')
    plt.colorbar()
    plt.subplot(144)
    plt.imshow(unc,vmax=10)
    plt.colorbar()
    plt.show()

    met = mod/unc

    mod = abs(gac_bt_ext-hirs_bt)
    unc = np.sqrt(((hirs_noise**2)+(gac_std_ext**2/3.)))
    met_ext = mod/unc

    plt.figure(figsize=(10,10))

    plot_uncertainties(gac_prob,cloud_frac,met,hirs_prob,'b','3x3 GAC')
    plot_uncertainties(gac_prob_ext,cloud_frac_ext,met_ext,hirs_prob,'r','5x3 GAC')

    plt.show()

    return

def explore_bias(hirs_bt,gac_bt,hirs_noise,gac_std,gac_n):

    samp = np.zeros([1000])
    diff = np.fliplr(hirs_bt)-gac_bt
    for i in np.arange(0,1000,1):
        samp[i] = np.nanmedian(np.random.choice(np.ravel(diff),size=1000))
    mean_median = np.nanmean(samp)
    print abs(np.nanmedian(diff)-mean_median)
    bias = (2*np.nanmedian(diff)) - mean_median
    diff = diff-bias

#    diff = diff-np.nanmedian(diff)
#    diff = diff-np.nanmean(diff)

    gac_unc = (gac_std/np.sqrt(gac_n-1))**2

    unc = np.sqrt((np.fliplr(hirs_noise)**2)+(gac_unc**2))

    metric = diff/unc
    
    plt.figure(figsize=[10,5])
    plt.subplot(131)
    plt.imshow(diff,cmap='bwr',vmin=-8,vmax=8)
    plt.title('T$_{HIRS}$-<T$_{GAC}$>')
    plt.colorbar()
    plt.subplot(132)
    plt.imshow(unc,vmax=5)
    plt.title(r'($\sigma_{HIRS}^{2}+( \frac{SD(GAC)}{(n-1)^{1/2}} )^2)^{1/2}$')
    plt.colorbar()
    plt.subplot(133)
    plt.imshow(metric,vmin=-10,vmax=10,cmap='bwr')
    plt.colorbar()
    plt.title('Difference/Uncertainty')

    gd = np.isfinite(metric)

    print np.nanmedian(metric[gd])

    plt.figure()
    plt.hist(metric[gd],bins=100,range=(-50,50))
    plt.xlabel('Diff/Unc Metric')
    plt.ylabel('Frequency')



    plt.show()

    return




