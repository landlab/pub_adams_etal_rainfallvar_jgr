#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Figure 11 from Adams et al., "The competition between frequent and rare flood events:
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
start = time.time()

(rmg, z) = read_esri_ascii('_data/topographic_elevation_4.asc', name='topographic__elevation')
(rmg, A) = read_esri_ascii('_data/drainage_area_4.asc', name='drainage_area')
log_bins = np.loadtxt('_data/log_bins.txt')

# Read in the model output data for the steady hydrology cases
highRvar = 'highRvar'
lowRvar = 'lowRvar'

twenty_m_erosion_highRvar = np.loadtxt('_data/' + highRvar +
        '/steadystate/20m_HalfMtrThresh_Diffusion_Tc0/total_eroded_depth.txt')
twenty_m_erosion_highRvar_tc10 = np.loadtxt('_data/' + highRvar +
        '/steadystate/20m_HalfMtrThresh_Diffusion_Tc1/total_eroded_depth.txt')
twenty_m_erosion_highRvar_tc25 = np.loadtxt('_data/' + highRvar +
        '/steadystate/20m_HalfMtrThresh_Diffusion_Tc5/total_eroded_depth.txt')
twenty_m_erosion_highRvar_tc50 = np.loadtxt('_data/' + highRvar +
        '/steadystate/20m_HalfMtrThresh_Diffusion_Tc10/total_eroded_depth.txt')

twenty_m_erosion_lowRvar = np.loadtxt('_data/' + lowRvar +
        '/steadystate/20m_HalfMtrThresh_Diffusion_Tc0/total_eroded_depth.txt')
twenty_m_erosion_lowRvar_tc10 = np.loadtxt('_data/' + lowRvar +
        '/steadystate/20m_HalfMtrThresh_Diffusion_Tc1/total_eroded_depth.txt')
twenty_m_erosion_lowRvar_tc25 = np.loadtxt('_data/' + lowRvar +
        '/steadystate/20m_HalfMtrThresh_Diffusion_Tc5/total_eroded_depth.txt')
twenty_m_erosion_lowRvar_tc50 = np.loadtxt('_data/' + lowRvar +
        '/steadystate/20m_HalfMtrThresh_Diffusion_Tc10/total_eroded_depth.txt')

# Read in the model output data for the nonsteady hydrology cases
highRvar = 'highRvar/nonsteady'
lowRvar = 'lowRvar/nonsteady'

nstwenty_m_erosion_highRvar = np.loadtxt('_data/' + highRvar +
        '/20m_HalfMtrThresh_Diffusion_Tc0/total_eroded_depth.txt')
nstwenty_m_erosion_highRvar_tc1 = np.loadtxt('_data/' + highRvar +
        '/20m_HalfMtrThresh_Diffusion_Tc1/total_eroded_depth.txt')
nstwenty_m_erosion_highRvar_tc5 = np.loadtxt('_data/' + highRvar +
        '/20m_HalfMtrThresh_Diffusion_Tc5/total_eroded_depth.txt')
nstwenty_m_erosion_highRvar_tc10 = np.loadtxt('_data/' + highRvar +
        '/20m_HalfMtrThresh_Diffusion_Tc10/total_eroded_depth.txt')

nstwenty_m_erosion_lowRvar = np.loadtxt('_data/' + lowRvar +
        '/20m_HalfMtrThresh_Diffusion_Tc0/total_eroded_depth.txt')
nstwenty_m_erosion_lowRvar_tc1 = np.loadtxt('_data/' + lowRvar +
        '/20m_HalfMtrThresh_Diffusion_Tc1/total_eroded_depth.txt')
nstwenty_m_erosion_lowRvar_tc5 = np.loadtxt('_data/' + lowRvar +
        '/20m_HalfMtrThresh_Diffusion_Tc5/total_eroded_depth.txt')
nstwenty_m_erosion_lowRvar_tc10 = np.loadtxt('_data/' + lowRvar +
        '/20m_HalfMtrThresh_Diffusion_Tc10/total_eroded_depth.txt')

# Creating the log-binned lists for drainage area and erosion data
(bin_ids) = np.digitize(A, log_bins)
unique_bins = np.unique(bin_ids)

avg_A = []
avg_twenty_m_erosion_highRvar = []
avg_twenty_m_erosion_highRvar_tc10 = []
avg_twenty_m_erosion_highRvar_tc25 = []
avg_twenty_m_erosion_highRvar_tc50 = []

avg_twenty_m_erosion_lowRvar = []
avg_twenty_m_erosion_lowRvar_tc10 = []
avg_twenty_m_erosion_lowRvar_tc25 = []
avg_twenty_m_erosion_lowRvar_tc50 = []

nsavg_twenty_m_erosion_highRvar = []
nsavg_twenty_m_erosion_highRvar_tc1 = []
nsavg_twenty_m_erosion_highRvar_tc5 = []
nsavg_twenty_m_erosion_highRvar_tc10 = []

nsavg_twenty_m_erosion_lowRvar = []
nsavg_twenty_m_erosion_lowRvar_tc1 = []
nsavg_twenty_m_erosion_lowRvar_tc5 = []
nsavg_twenty_m_erosion_lowRvar_tc10 = []

# Log-bin average the erosion data by drainage area
for each in unique_bins:
    locs = np.where(bin_ids == each)

    avg_A.append(np.median(A[locs]))

    avg_twenty_m_erosion_highRvar.append(np.median(twenty_m_erosion_highRvar[locs]))
    avg_twenty_m_erosion_lowRvar.append(np.median(twenty_m_erosion_lowRvar[locs]))
    avg_twenty_m_erosion_highRvar_tc10.append(np.median(twenty_m_erosion_highRvar_tc10[locs]))
    avg_twenty_m_erosion_highRvar_tc25.append(np.median(twenty_m_erosion_highRvar_tc25[locs]))
    avg_twenty_m_erosion_highRvar_tc50.append(np.median(twenty_m_erosion_highRvar_tc50[locs]))
    avg_twenty_m_erosion_lowRvar_tc10.append(np.median(twenty_m_erosion_lowRvar_tc10[locs]))
    avg_twenty_m_erosion_lowRvar_tc25.append(np.median(twenty_m_erosion_lowRvar_tc25[locs]))
    avg_twenty_m_erosion_lowRvar_tc50.append(np.median(twenty_m_erosion_lowRvar_tc50[locs]))

    nsavg_twenty_m_erosion_highRvar.append(np.median(nstwenty_m_erosion_highRvar[locs]))
    nsavg_twenty_m_erosion_lowRvar.append(np.median(nstwenty_m_erosion_lowRvar[locs]))
    nsavg_twenty_m_erosion_highRvar_tc1.append(np.median(nstwenty_m_erosion_highRvar_tc1[locs]))
    nsavg_twenty_m_erosion_highRvar_tc5.append(np.median(nstwenty_m_erosion_highRvar_tc5[locs]))
    nsavg_twenty_m_erosion_highRvar_tc10.append(np.median(nstwenty_m_erosion_highRvar_tc10[locs]))
    nsavg_twenty_m_erosion_lowRvar_tc1.append(np.median(nstwenty_m_erosion_lowRvar_tc1[locs]))
    nsavg_twenty_m_erosion_lowRvar_tc5.append(np.median(nstwenty_m_erosion_lowRvar_tc5[locs]))
    nsavg_twenty_m_erosion_lowRvar_tc10.append(np.median(nstwenty_m_erosion_lowRvar_tc10[locs]))

# Convert meters to millimeters
avg_twenty_m_erosion_highRvar = [x * 1000 for x in avg_twenty_m_erosion_highRvar]
avg_twenty_m_erosion_highRvar_tc10 = [x * 1000 for x in avg_twenty_m_erosion_highRvar_tc10]
avg_twenty_m_erosion_highRvar_tc25 = [x * 1000 for x in avg_twenty_m_erosion_highRvar_tc25]
avg_twenty_m_erosion_highRvar_tc50 = [x * 1000 for x in avg_twenty_m_erosion_highRvar_tc50]

avg_twenty_m_erosion_lowRvar = [x * 1000 for x in avg_twenty_m_erosion_lowRvar]
avg_twenty_m_erosion_lowRvar_tc10 = [x * 1000 for x in avg_twenty_m_erosion_lowRvar_tc10]
avg_twenty_m_erosion_lowRvar_tc25 = [x * 1000 for x in avg_twenty_m_erosion_lowRvar_tc25]
avg_twenty_m_erosion_lowRvar_tc50 = [x * 1000 for x in avg_twenty_m_erosion_lowRvar_tc50]

nsavg_twenty_m_erosion_highRvar = [x * 1000 for x in nsavg_twenty_m_erosion_highRvar]
nsavg_twenty_m_erosion_highRvar_tc1 = [x * 1000 for x in nsavg_twenty_m_erosion_highRvar_tc1]
nsavg_twenty_m_erosion_highRvar_tc5= [x * 1000 for x in nsavg_twenty_m_erosion_highRvar_tc5]
nsavg_twenty_m_erosion_highRvar_tc10 = [x * 1000 for x in nsavg_twenty_m_erosion_highRvar_tc10]

nsavg_twenty_m_erosion_lowRvar = [x * 1000 for x in nsavg_twenty_m_erosion_lowRvar]
nsavg_twenty_m_erosion_lowRvar_tc1 = [x * 1000 for x in nsavg_twenty_m_erosion_lowRvar_tc1]
nsavg_twenty_m_erosion_lowRvar_tc5 = [x * 1000 for x in nsavg_twenty_m_erosion_lowRvar_tc5]
nsavg_twenty_m_erosion_lowRvar_tc10 = [x * 1000 for x in nsavg_twenty_m_erosion_lowRvar_tc10]

# Plotting erosion rate data by drainage area - we are normalizing by the outlet
# data, which is index number 37 in the list.
f, (ax1, ax2) = plt.subplots(2, sharex=True, figsize=(10, 5))
ax1.tick_params(axis='y', which='minor', direction='in',bottom='on', top='on', right='on', left='on')

ax1.loglog(avg_A, np.ones(len(avg_A)), 'k-')
ax1.loglog(avg_A, avg_twenty_m_erosion_highRvar/(avg_twenty_m_erosion_highRvar[37]), 'x', label='Steady-state, a = 1', color='grey')
ax1.loglog(avg_A, nsavg_twenty_m_erosion_highRvar/(nsavg_twenty_m_erosion_highRvar[37]), 'x', label='Nonsteady, a = 1', color='black')#, mfc='none')
ax1.loglog(avg_A, avg_twenty_m_erosion_lowRvar/(avg_twenty_m_erosion_lowRvar[37]), 'o', label='Steady-state, a = 1', color='grey', mfc='none')
ax1.loglog(avg_A, nsavg_twenty_m_erosion_lowRvar/(nsavg_twenty_m_erosion_lowRvar[37]), 'o', label='Nonsteady, a = 1', color='black', mfc='none')
ax1.set_ylim(10**-2, 2)

ax2.tick_params(direction='in',bottom='on', top='on', right='on', left='on')
ax2.tick_params(axis='y', which='minor', direction='in',bottom='on', top='on', right='on', left='on')
ax2.loglog(avg_A, np.ones(len(avg_A)), 'k-')
ax2.loglog(avg_A, avg_twenty_m_erosion_highRvar_tc10/(avg_twenty_m_erosion_highRvar_tc10[37]), 'x', label='Steady-state, a = 1', color='grey')
ax2.loglog(avg_A, nsavg_twenty_m_erosion_highRvar_tc10/(nsavg_twenty_m_erosion_highRvar_tc10[37]), 'x', label='Nonsteady, a = 1', color='black')#, mfc='none')
ax2.loglog(avg_A, avg_twenty_m_erosion_lowRvar_tc50/(avg_twenty_m_erosion_lowRvar_tc50[37]), 'o', label='Steady-state, a = 1', color='grey', mfc='none')
ax2.loglog(avg_A, nsavg_twenty_m_erosion_lowRvar_tc10/(nsavg_twenty_m_erosion_lowRvar_tc10[37]), 'o', label='Nonsteady, a = 1', color='black', mfc='none')
ax2.set_xlim(10**3, 10**9)
ax2.set_ylim(10**-7, 10)
ax2.set_xlabel('Drainage area (m$^{2}$)')
ax2.set_ylabel('Total eroded depth (mm)')

plt.show()
