#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Figure 12 from Adams et al., "The competition between frequent and rare flood events:
 impacts on erosion rates and landscape form"

Written by Jordan M. Adams
Updated April 14, 2020
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import norm

# Read in output peak discharge data for the nonsteady and steady cases
low_rvar_qp = np.loadtxt('_data/lowRvar/nonsteady/20mGeneral/each_storm_peak_discharge_20m.txt')
high_rvar_qp = np.loadtxt('_data/highRvar/nonsteady/20mGeneral/each_storm_peak_discharge_20m.txt')

high_rvar_ss = np.loadtxt('_data/highRvar/steadystate/20mGeneral/peak_q_20m.txt')
low_rvar_ss = np.loadtxt('_data/lowRvar/steadystate/20mGeneral/peak_q_20m.txt')

high_rvar_ss = high_rvar_ss[:,0]
low_rvar_ss = low_rvar_ss[:,0]

high_rvar_outlet_med = np.loadtxt('_data/highRvar/nonsteady/20mGeneral/outlet_median_peak_discharges.txt')
low_rvar_outlet_med = np.loadtxt('_data/lowRvar/nonsteady/20mGeneral/outlet_median_peak_discharges.txt')

high_rvar_upstream_med = np.loadtxt('_data/highRvar/nonsteady/20mGeneral/upstream_median_peak_discharges.txt')
low_rvar_upstream_med = np.loadtxt('_data/lowRvar/nonsteady/20mGeneral/upstream_median_peak_discharges.txt')

inds = high_rvar_outlet_med.argsort()
inds_b = low_rvar_outlet_med.argsort()
inds_ss = high_rvar_ss.argsort()
inds_b_ss = low_rvar_ss.argsort()

# Create Figure 12 Panel A
plt.figure("PDFs of Discharge", figsize=(5, 5))
plt.title('Comparison of pdfs')
ax = plt.gca()
plt.loglog(np.sort(low_rvar_outlet_med), norm.pdf(low_rvar_outlet_med, scale=np.nanmean(low_rvar_outlet_med))[inds_b], label='Low R$_{var}$', color='navy', linewidth=2)
plt.loglog(np.sort(high_rvar_outlet_med),norm.pdf(high_rvar_outlet_med, scale=np.nanmean(high_rvar_outlet_med))[inds],'-' ,label='High R$_{var}$', color='#D9B39B', linewidth=2)
plt.legend()
plt.ylim(10**-5, 10**-1)
plt.xlim(10**0, 200)
ax.tick_params(direction='in',bottom='on', top='on', right='on', left='on')
plt.xlabel('Median event discharge (m$^{3}$s$^{-1}$)')
plt.ylabel('Relative frequency')
ax = plt.gca()
plt.loglog(np.sort(high_rvar_ss), norm.pdf(high_rvar_ss, scale=np.nanmean(high_rvar_ss))[inds_ss], '--',color='#D9B39B',  alpha=0.25, linewidth=4)
plt.loglog(np.sort(low_rvar_ss), norm.pdf(low_rvar_ss, scale=np.nanmean(low_rvar_ss))[inds_b_ss], '--', color='navy',alpha=0.15, linewidth=4)
plt.legend()
plt.ylim(10**-5, 10**-1)
plt.xlim(10**0, 200)
ax.tick_params(direction='in',bottom='on', top='on', right='on', left='on')
plt.xlabel('Median event discharge (m$^{3}$s$^{-1}$)')
plt.ylabel('Relative frequency')

highRvar_folder = 'highRvar/nonsteady'
lowRvar_folder = 'lowRvar/nonsteady'

highRvarQP = np.loadtxt('_data/' + highRvar_folder + '/20mGeneral/each_storm_peak_discharge_20m.txt')
lowRvarQP = np.loadtxt('_data/' + lowRvar_folder + '/20mGeneral/each_storm_peak_discharge_20m.txt')

highRvar_ss_folder = 'highRvar/steadystate'
lowRvar_ss_folder = 'lowRvar/steadystate'
highRvarQP_ss = np.loadtxt('_data/' + highRvar_ss_folder + '/20mGeneral/peak_q_20m.txt')
lowRvarQP_ss = np.loadtxt('_data/' + lowRvar_ss_folder + '/20mGeneral/peak_q_20m.txt')

highRvar_outlet = highRvarQP[:,0]
highRvar_nearoutlet = highRvarQP[:,1]
highRvar_midstream = highRvarQP[:,2]
highRvar_upstream = highRvarQP[:,3]

lowRvar_outlet = lowRvarQP[:,0]
lowRvar_nearoutlet = lowRvarQP[:,1]
lowRvar_midstream = lowRvarQP[:,2]
lowRvar_upstream = lowRvarQP[:,3]

highRvar_outlet_ss = highRvarQP_ss[:,0]
highRvar_nearoutlet_ss = highRvarQP_ss[:,1]
highRvar_midstream_ss = highRvarQP_ss[:,2]
highRvar_upstream_ss = highRvarQP_ss[:,3]

lowRvar_outlet_ss = lowRvarQP_ss[:,0]
lowRvar_nearoutlet_ss = lowRvarQP_ss[:,1]
lowRvar_midstream_ss = lowRvarQP_ss[:,2]
lowRvar_upstream_ss = lowRvarQP_ss[:,3]

plt.figure("Outlet", figsize=(5, 5))
ax = plt.gca()
plt.plot(np.arange(0, 150), np.arange(0, 150), '-', color='k', linewidth=1)
plt.plot(high_rvar_outlet_med, highRvar_outlet_ss[:1001], 'x', mfc='none', color='#D9B39B', markersize=6, alpha=0.4)
plt.plot(low_rvar_outlet_med, lowRvar_outlet_ss[:1001], '.', mfc='none', color='navy', markersize=10, alpha=0.4)
plt.ylabel('Steady-state discharge (m$^{3}$ s$^{-1}$)')
plt.xlabel('Nonsteady median discharge (m$^{3}$ s$^{-1}$)')
ax.set_xticks(np.arange(0, 150, 25))
ax.set_yticks(np.arange(0, 150, 25))
plt.xlim(0, 125)
plt.ylim(0, 125)

plt.figure("Upstream", figsize=(5, 5))
ax = plt.gca()
plt.plot(np.arange(0, 100), np.arange(0, 100), '-', color='k', linewidth=1)
plt.plot(high_rvar_upstream_med, highRvar_upstream_ss[:1001], 'x', mfc='none', color='#D9B39B', markersize=6)
plt.plot(low_rvar_upstream_med, lowRvar_upstream_ss[:1001], '.', mfc='none', color='navy', markersize=10)
plt.ylabel('Steady-state discharge (m$^{3}$ s$^{-1}$)')
ax.set_xticks(np.arange(0, 30, 5))
ax.set_yticks(np.arange(0, 30, 5))
plt.xlim(0, 25)
plt.ylim(0, 25)
plt.show()
