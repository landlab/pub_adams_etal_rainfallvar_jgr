 #/!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Figure 5 from Adams et al., "The competition between frequent and rare flood events:
 impacts on erosion rates and landscape form"

Written by Jordan M. Adams

Updated April 14, 2020
"""
from landlab.io import read_esri_ascii, write_esri_ascii
from landlab.components import OverlandFlow, PrecipitationDistribution, DepthSlopeProductErosion, LinearDiffuser
import numpy as np
from landlab.utils.depth_dependent_roughness import depth_dependent_mannings_n
from matplotlib import pyplot as plt
from landlab.plot import imshow_grid
from scipy.stats import norm
import time
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Helvetica"
plt.rcParams['font.size'] = 10

## Generate 1000 year precipitation time series
total_t = 1000.*365.25*24

lowRvar_PrecipDist = PrecipitationDistribution(mean_storm_duration = 11.75,
                               mean_interstorm_duration = 146.25,
                               mean_storm_depth = 4.775,
                               total_t=total_t)

highRvar_PrecipDist = PrecipitationDistribution(mean_storm_duration = 10.75,
                               mean_interstorm_duration = 433.58,
                               mean_storm_depth = 9.62,
                               total_t=total_t)


thresh = 0.5

# Get actual time series for the high Rvar case
highRvar_storm_arr = np.array(highRvar_PrecipDist.get_storm_time_series())
highRvar_intensity_threshold, = np.where(highRvar_storm_arr[:, 2] > thresh)
highRvar_durations = highRvar_storm_arr[highRvar_intensity_threshold][:,1]-highRvar_storm_arr[highRvar_intensity_threshold][:,0]
highRvar_durations_s = [x *3600. for x in highRvar_durations]
highRvar_intensities = highRvar_storm_arr[highRvar_intensity_threshold][:,2]
highRvar_depth = highRvar_durations * highRvar_intensities
highRvar_durations_raw  = highRvar_storm_arr[:,1]-highRvar_storm_arr[:,0]
highRvar_intensities_raw = highRvar_storm_arr[:,2]
highRvar_depth_raw = highRvar_durations_raw * highRvar_intensities_raw

# Get actual time series for the low Rvar case
lowRvar_storm_arr = np.array(lowRvar_PrecipDist.get_storm_time_series())
lowRvar_intensity_threshold, = np.where(lowRvar_storm_arr[:, 2] > thresh)
lowRvar_durations = lowRvar_storm_arr[lowRvar_intensity_threshold][:,1]-lowRvar_storm_arr[lowRvar_intensity_threshold][:,0]
lowRvar_durations_s = [x *3600. for x in lowRvar_durations]
lowRvar_intensities = lowRvar_storm_arr[lowRvar_intensity_threshold][:,2]
lowRvar_depth = lowRvar_durations * lowRvar_intensities
lowRvar_durations_raw  = lowRvar_storm_arr[:,1]-lowRvar_storm_arr[:,0]
lowRvar_intensities_raw = lowRvar_storm_arr[:,2]
lowRvar_depth_raw = lowRvar_durations_raw * lowRvar_intensities_raw

## Histogram weighting
highRvar_raw_weights = np.ones_like(highRvar_intensities_raw)/float(len(highRvar_intensities_raw))
highRvar_weights = np.ones_like(highRvar_intensities)/float(len(highRvar_intensities))
lowRvar_raw_weights = np.ones_like(lowRvar_intensities_raw)/float(len(lowRvar_intensities_raw))
lowRvar_weights = np.ones_like(lowRvar_intensities)/float(len(lowRvar_intensities))

dhighRvar_raw_weights = np.ones_like(highRvar_depth_raw)/float(len(highRvar_depth_raw))
dhighRvar_weights = np.ones_like(highRvar_depth)/float(len(highRvar_depth))
dlowRvar_raw_weights = np.ones_like(lowRvar_depth_raw)/float(len(lowRvar_depth_raw))
dlowRvar_weights = np.ones_like(lowRvar_depth)/float(len(lowRvar_depth))

## CDF creation
lowRvar_dur_cdf = norm.cdf(lowRvar_durations, np.average(lowRvar_durations), np.std(lowRvar_durations))
highRvar_dur_cdf = norm.cdf(highRvar_durations, np.average(highRvar_durations), np.std(highRvar_durations))

lowRvar_depth_cdf = norm.cdf(lowRvar_depth, np.average(lowRvar_depth), np.std(lowRvar_depth))
highRvar_depth_cdf = norm.cdf(highRvar_depth, np.average(highRvar_depth), np.std(highRvar_depth))

lowRvar_int_cdf = norm.cdf(lowRvar_intensities, np.average(lowRvar_intensities), np.std(lowRvar_intensities))
highRvar_int_cdf = norm.cdf(highRvar_intensities, np.average(highRvar_intensities), np.std(highRvar_intensities))

plt.figure("Panel A")
int_bins = np.arange(0, 5, 0.1)
ax = plt.gca()
plt.title('Comparison sampled high/low $R_{var}$ intensity probabilities')
plt.hist(highRvar_intensities, int_bins, color='saddlebrown', label='Sampled high $R_{var}$ Data', weights =highRvar_weights, histtype='step',linewidth=1.5)#alpha=0.4)#  histtype='step')
plt.hist(lowRvar_intensities, int_bins, color='teal', label='Sampled low $R_{var}$ Data', histtype='step',weights=lowRvar_weights,linewidth=1.5)
plt.legend()
ax.tick_params(direction='in',bottom='on', top='on', right='on', left='on')
plt.xlabel('Intensity (mm/hr)')
plt.ylabel('Probability')
plt.ylim(0, 0.3)
ax.tick_params(direction='in',bottom='on', top='on', right='on', left='on')
ax.set_xticks(np.arange(0, 5, 1))
plt.xlim(0, 4)
plt.show()

inds_i = highRvar_intensities.argsort()
inds_bi = lowRvar_intensities.argsort()
plt.figure("Panel B")
ax = plt.gca()

plt.loglog(np.sort(highRvar_intensities),norm.pdf(highRvar_intensities, scale=np.average(highRvar_intensities))[inds_i],'-' ,label='Sampled high R$_{var}$', color='saddlebrown')
plt.loglog(np.sort(lowRvar_intensities), norm.pdf(lowRvar_intensities, scale=np.average(lowRvar_intensities))[inds_bi], label='Sampled low R$_{var}$', color='teal')
plt.legend()
plt.ylim(10**-3, 10**-0)
plt.xlim(5*(10**-1), 10**1)
plt.xlabel('Rainfall intensity for an event (mm hr$^{-1}$)')
plt.ylabel('Relative frequency')
ax.tick_params(direction='in',bottom='on', top='on', right='on', left='on')
plt.show()

depth_bins = np.arange(0, 100, 2)
plt.figure("Panel C")
ax = plt.gca()
plt.title('Comparison sampled high/low $R_{var}$ depth probabilities')
plt.hist(highRvar_depth, depth_bins, color='saddlebrown', label='Sampled high $R_{var}$ Data',weights=dhighRvar_weights,  histtype='step', linewidth=1.5)# alpha=0.4)
plt.hist(lowRvar_depth, depth_bins, color='teal', label='Sampled low $R_{var}$ Data', weights=dlowRvar_weights, histtype='step', linewidth=1.5)
plt.legend()
plt.xlabel('Depth (mm)')
plt.ylabel('Probability')
plt.ylim(0, 0.3)
ax.tick_params(direction='in',bottom='on', top='on', right='on', left='on')
ax.set_xticks(np.arange(0, 90, 15))
plt.xlim(0, 60)
plt.show()


inds = highRvar_depth.argsort()
inds_b = lowRvar_depth.argsort()
plt.figure("Panel D")
ax = plt.gca()
plt.loglog(np.sort(highRvar_depth),norm.pdf(highRvar_depth, scale=np.average(highRvar_depth))[inds],'-' ,label='Sampled high R$_{var}$', color='saddlebrown')
plt.loglog(np.sort(lowRvar_depth), norm.pdf(lowRvar_depth, scale=np.average(lowRvar_depth))[inds_b], label='Sampled low R$_{var}$', color='teal')
plt.legend()
plt.ylim(10**-6, 10**-1)
plt.xlim(10**0, 10**2)
ax.tick_params(direction='in',bottom='on', top='on', right='on', left='on')
plt.xlabel('Rainfall depth for an event (mm)')
plt.ylabel('Relative frequency')
plt.show()
