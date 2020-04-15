#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Figure 9 from Adams et al., "The competition between frequent and rare flood events:
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
import time

(rmg, z) = read_esri_ascii('_data/topographic_elevation_4.asc', name='topographic__elevation')
(rmg, A) = read_esri_ascii('_data/drainage_area_4.asc', name='drainage_area')
log_bins = np.loadtxt('_data/log_bins.txt')

# Read in the data output
highRvar = 'highRvar/nonsteady'
lowRvar = 'lowRvar/nonsteady'

twenty_m_erosion_highRvar = np.loadtxt('_data/' + highRvar +
            '/20m_HalfMtrThresh_Diffusion_Tc0/total_eroded_depth.txt')
twenty_m_erosion_highRvar_tc1 = np.loadtxt('_data/' + highRvar +
            '/20m_HalfMtrThresh_Diffusion_Tc1/total_eroded_depth.txt')
twenty_m_erosion_highRvar_tc5 = np.loadtxt('_data/' + highRvar +
            '/20m_HalfMtrThresh_Diffusion_Tc5/total_eroded_depth.txt')
twenty_m_erosion_highRvar_tc10 = np.loadtxt('_data/' + highRvar +
            '/20m_HalfMtrThresh_Diffusion_Tc10/total_eroded_depth.txt')

twenty_m_erosion_lowRvar = np.loadtxt('_data/' + lowRvar +
            '/20m_HalfMtrThresh_Diffusion_Tc0/total_eroded_depth.txt')
twenty_m_erosion_lowRvar_tc1 = np.loadtxt('_data/' + lowRvar +
            '/20m_HalfMtrThresh_Diffusion_Tc1/total_eroded_depth.txt')
twenty_m_erosion_lowRvar_tc5 = np.loadtxt('_data/' + lowRvar +
            '/20m_HalfMtrThresh_Diffusion_Tc5/total_eroded_depth.txt')
twenty_m_erosion_lowRvar_tc10 = np.loadtxt('_data/' + lowRvar +
            '/20m_HalfMtrThresh_Diffusion_Tc10/total_eroded_depth.txt')

# Lists for the log binned average depths.
(bin_ids) = np.digitize(A, log_bins)
unique_bins = np.unique(bin_ids)

avg_A = []
avg_twenty_m_erosion_highRvar = []
avg_twenty_m_erosion_highRvar_tc1 = []
avg_twenty_m_erosion_highRvar_tc5 = []
avg_twenty_m_erosion_highRvar_tc10 = []

avg_twenty_m_erosion_lowRvar = []
avg_twenty_m_erosion_lowRvar_tc1 = []
avg_twenty_m_erosion_lowRvar_tc5 = []
avg_twenty_m_erosion_lowRvar_tc10 = []

# Binning each the erosion depths by drainage area
for each in unique_bins:
    locs = np.where(bin_ids == each)

    avg_A.append(np.median(A[locs]))

    avg_twenty_m_erosion_highRvar.append(np.median(twenty_m_erosion_highRvar[locs]))
    avg_twenty_m_erosion_lowRvar.append(np.median(twenty_m_erosion_lowRvar[locs]))

    avg_twenty_m_erosion_highRvar_tc1.append(np.median(twenty_m_erosion_highRvar_tc1[locs]))
    avg_twenty_m_erosion_highRvar_tc5.append(np.median(twenty_m_erosion_highRvar_tc5[locs]))
    avg_twenty_m_erosion_highRvar_tc10.append(np.median(twenty_m_erosion_highRvar_tc10[locs]))

    avg_twenty_m_erosion_lowRvar_tc1.append(np.median(twenty_m_erosion_lowRvar_tc1[locs]))
    avg_twenty_m_erosion_lowRvar_tc5.append(np.median(twenty_m_erosion_lowRvar_tc5[locs]))
    avg_twenty_m_erosion_lowRvar_tc10.append(np.median(twenty_m_erosion_lowRvar_tc10[locs]))

# Conversion from meters to millimeters
avg_twenty_m_erosion_highRvar = [x * 1000 for x in avg_twenty_m_erosion_highRvar]
avg_twenty_m_erosion_highRvar_tc1 = [x * 1000 for x in avg_twenty_m_erosion_highRvar_tc1]
avg_twenty_m_erosion_highRvar_tc5 = [x * 1000 for x in avg_twenty_m_erosion_highRvar_tc5]
avg_twenty_m_erosion_highRvar_tc10 = [x * 1000 for x in avg_twenty_m_erosion_highRvar_tc10]

avg_twenty_m_erosion_lowRvar = [x * 1000 for x in avg_twenty_m_erosion_lowRvar]
avg_twenty_m_erosion_lowRvar_tc1 = [x * 1000 for x in avg_twenty_m_erosion_lowRvar_tc1]
avg_twenty_m_erosion_lowRvar_tc5 = [x * 1000 for x in avg_twenty_m_erosion_lowRvar_tc5]
avg_twenty_m_erosion_lowRvar_tc10 = [x * 1000 for x in avg_twenty_m_erosion_lowRvar_tc10]

# Generating all of the plots for Figure 9.
plt.figure(1, figsize=(8,8))
ax = plt.gca()
ax.tick_params(direction='in',bottom='on', top='on', right='on', left='on')

plt.title('High R$_{var}$ - Changing Tc')
plt.loglog(avg_A, avg_twenty_m_erosion_highRvar, 'x',
    label='20 m rainfall depth, high R$_{var}$, Tc = 0 Pa', color='#dfc27d')
plt.loglog(avg_A, avg_twenty_m_erosion_highRvar_tc1,'x',
    label='20 m rainfall depth, high R$_{var}$, Tc = 1 Pa', color='#bf812d')
plt.loglog(avg_A, avg_twenty_m_erosion_highRvar_tc5,'x',
    label='20 m rainfall depth, high R$_{var}$, Tc = 5 Pa', color='#8c510a')
plt.loglog(avg_A, avg_twenty_m_erosion_highRvar_tc10, 'x',
    label='20 m rainfall depth, high R$_{var}$, Tc = 10 Pa', color='#543005')

plt.xlim(10**3, 10**9)
plt.ylim(10**-6, 10**0)
plt.legend()
plt.xlabel('Drainage area (m$^{2}$)')
plt.ylabel('Total eroded depth (mm)')

plt.figure(2,figsize=(8,8))
plt.title('Low R$_{var}$ - Changing Tc')
plt.loglog(avg_A, avg_twenty_m_erosion_lowRvar, 'o',
    label='20 m rainfall depth, low R$_{var}$, Tc = 0 Pa',
    color='lightskyblue', mfc='none')
plt.loglog(avg_A, avg_twenty_m_erosion_lowRvar_tc1,'o',
    label='20 m rainfall depth, low R$_{var}$, Tc = 1 Pa',
    color='cornflowerblue', mfc='none')
plt.loglog(avg_A, avg_twenty_m_erosion_lowRvar_tc5,'o',
    label='20 m rainfall depth, low R$_{var}$, Tc = 5 Pa',
    color='royalblue', mfc='none')
plt.loglog(avg_A, avg_twenty_m_erosion_lowRvar_tc10, 'o',
    label='20 m rainfall depth, low R$_{var}$, Tc = 10 Pa',
    color='navy', mfc='none')
plt.xlim(10**3, 10**9)
plt.ylim(10**-6, 10**0)
plt.xlabel('Drainage area (m$^{2}$)')
plt.ylabel('Total eroded depth (mm)')
plt.legend()

plt.figure(3,figsize=(4,4))
plt.title('High R$_{var}$ vs. Low R$_{var}$ Tc = 0 Pa')
plt.loglog(avg_A, avg_twenty_m_erosion_highRvar, 'x', label='20 m rainfall depth, high R$_{var}$', color='black')
plt.loglog(avg_A, avg_twenty_m_erosion_lowRvar,'o', label='20 m rainfall depth, low R$_{var}$', color='black', mfc='none')
plt.xlim(10**3, 10**9)
plt.ylim(10**-2, 10**0)
plt.legend()
plt.xlabel('Drainage area (m$^{2}$)')
plt.ylabel('Total eroded depth (mm)')

plt.figure(4,figsize=(4,4))
plt.title('High R$_{var}$ vs. Low R$_{var}$ Tc = 1 Pa')
plt.loglog(avg_A, avg_twenty_m_erosion_highRvar_tc1, 'x', label='20 m rainfall depth, high R$_{var}$', color='black')
plt.loglog(avg_A, avg_twenty_m_erosion_lowRvar_tc1,'o', label='20 m rainfall depth, low R$_{var}$', color='black', mfc='none')
plt.xlim(10**3, 10**9)
plt.ylim(10**-2, 10**0)
plt.legend()
plt.xlabel('Drainage area (m$^{2}$)')
plt.ylabel('Total eroded depth (mm)')

plt.figure(5,figsize=(4,4))
plt.title('High R$_{var}$ vs. Low R$_{var}$ Tc = 5 Pa')
plt.loglog(avg_A, avg_twenty_m_erosion_highRvar_tc5, 'x', label='20 m rainfall depth, high R$_{var}$', color='black')
plt.loglog(avg_A, avg_twenty_m_erosion_lowRvar_tc5,'o', label='20 m rainfall depth, low R$_{var}$', color='black', mfc='none')
plt.xlim(10**3, 10**9)
plt.ylim(10**-4, 10**0)
plt.legend()
plt.xlabel('Drainage area (m$^{2}$)')
plt.ylabel('Total eroded depth (mm)')

plt.figure(6,figsize=(4,4))
plt.title('High R$_{var}$ vs. Low R$_{var}$ Tc = 10 Pa')
plt.loglog(avg_A, avg_twenty_m_erosion_highRvar_tc10, 'x', label='20 m rainfall depth, high R$_{var}$', color='black')
plt.loglog(avg_A, avg_twenty_m_erosion_lowRvar_tc10,'o', label='20 m rainfall depth, low R$_{var}$', color='black', mfc='none')
plt.xlim(10**3, 10**9)
plt.ylim(10**-5, 10**0)
plt.legend()
plt.xlabel('Drainage area (m$^{2}$)')
plt.ylabel('Total eroded depth (mm)')
plt.show()
