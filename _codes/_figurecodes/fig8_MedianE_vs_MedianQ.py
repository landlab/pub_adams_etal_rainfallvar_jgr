#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Figure 8 from Adams et al., "The competition between frequent and rare flood events:
 impacts on erosion rates and landscape form"

Written by Jordan M. Adams
Updated April 14, 2020
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import norm

# Peak discharges for nonsteady hydrology cases
low_rvar_qp = np.loadtxt('_data/lowRvar/nonsteady/20mGeneral/each_storm_peak_discharge_20m.txt')
high_rvar_qp = np.loadtxt('_data/highRvar/nonsteady/20mGeneral/each_storm_peak_discharge_20m.txt')

# Median discharges for nonsteady hydrology cases
high_rvar_outlet_med = np.loadtxt('_data/highRvar/nonsteady/20mGeneral/outlet_median_peak_discharges.txt')
low_rvar_outlet_med = np.loadtxt('_data/lowRvar/nonsteady/20mGeneral/outlet_median_peak_discharges.txt')

high_rvar_upstream_med = np.loadtxt('_data/highRvar/nonsteady/20mGeneral/upstream_median_peak_discharges.txt')
low_rvar_upstream_med = np.loadtxt('_data/lowRvar/nonsteady/20mGeneral/upstream_median_peak_discharges.txt')

high_rvar_outlet_emedian = np.loadtxt('_data/highRvar/nonsteady/20mGeneral/outlet_median_erosion_rates_tc0.txt')
low_rvar_outlet_emedian = np.loadtxt('_data/lowRvar/nonsteady/20mGeneral/outlet_median_erosion_rates_tc0.txt')

# Median erosion rates for nonsteady hydrology cases, Tc = 0 Pa
high_rvar_upstream_emedian = np.loadtxt('_data/highRvar/nonsteady/20mGeneral/upstream_median_erosion_rates_tc0.txt')
low_rvar_upstream_emedian = np.loadtxt('_data/lowRvar/nonsteady/20mGeneral/upstream_median_erosion_rates_tc0.txt')

high_rvar_midstream_med = np.loadtxt('_data/highRvar/nonsteady/20mGeneral/downstream1_median_peak_discharges.txt')
low_rvar_midstream_med = np.loadtxt('_data/lowRvar/nonsteady/20mGeneral/downstream1_median_peak_discharges.txt')

high_rvar_midstream_emedian = np.loadtxt('_data/highRvar/nonsteady/20mGeneral/downstream1_median_erosion_rates_tc0.txt')
low_rvar_midstream_emedian = np.loadtxt('_data/lowRvar/nonsteady/20mGeneral/downstream1_median_erosion_rates_tc0.txt')

# Median erosion rates for nonsteady hydrology cases, Tc = 10 Pa
high_rvar_outlet_emedian_tc10 = np.loadtxt('_data/highRvar/nonsteady/20mGeneral/outlet_median_erosion_rates_tc10.txt')
low_rvar_outlet_emedian_tc10 = np.loadtxt('_data/lowRvar/nonsteady/20mGeneral/outlet_median_erosion_rates_tc10.txt')

high_rvar_upstream_emedian_tc10 = np.loadtxt('_data/highRvar/nonsteady/20mGeneral/upstream_median_erosion_rates_tc10.txt')
low_rvar_upstream_emedian_tc10 = np.loadtxt('_data/lowRvar/nonsteady/20mGeneral/upstream_median_erosion_rates_tc10.txt')

high_rvar_midstream_emedian_tc10 = np.loadtxt('_data/highRvar/nonsteady/20mGeneral/downstream1_median_erosion_rates_tc10.txt')
low_rvar_midstream_emedian_tc10 = np.loadtxt('_data/lowRvar/nonsteady/20mGeneral/downstream1_median_erosion_rates_tc10.txt')

# Sorting indicies for plotting
inds = high_rvar_outlet_med.argsort()
inds_b = low_rvar_outlet_med.argsort()

## Creation of the 6 subpots.
f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharey=False, sharex=False, figsize=(8, 6))
f.add_subplot(111, frameon=False)
plt.ylabel('Median erosion rate in nonsteady method (m s$^{-1}$)')

ax1.plot(high_rvar_upstream_med[:500], high_rvar_upstream_emedian[:500], 'x', color='#D9B39B')
ax1.plot(low_rvar_upstream_med[:500], low_rvar_upstream_emedian[:500], 'o', color='navy', mfc='none')
ax1.set_ylim(0, 1.5e-11)
ax1.set_xlim(0, 20)
ax1.text(1, 1.40e-11, 'upstream, Tc = 0 Pa', fontsize=8)

ax2.plot(high_rvar_midstream_med[:500], high_rvar_midstream_emedian[:500], 'x', color='#D9B39B')
ax2.plot(low_rvar_midstream_med[:500], low_rvar_midstream_emedian[:500], 'o', color='navy', mfc='none')
ax2.set_ylim(0, 2.0e-11)
ax2.set_xlim(0, 100)
ax2.text(5, 1.85e-11, 'midstream, Tc = 0 Pa', fontsize=8)

ax3.plot(high_rvar_outlet_med[:500], high_rvar_outlet_emedian[:500], 'x', color='#D9B39B')
ax3.plot(low_rvar_outlet_med[:500], low_rvar_outlet_emedian[:500], 'o', color='navy', mfc='none')
ax3.text(5, 2.75e-11, 'outlet, Tc = 0 Pa', fontsize=8)
ax3.set_xlim(0, 100)
ax3.set_ylim(0, 3.0e-11)

ax4.plot(high_rvar_upstream_med[:500], high_rvar_upstream_emedian_tc10[:500], 'x', color='#D9B39B', alpha=0.25)
ax4.plot(low_rvar_upstream_med[:500], low_rvar_upstream_emedian_tc10[:500], 'o', color='navy', mfc='none', alpha=0.25)
ax4.set_ylim(0, 1.5e-11)
ax4.set_xlim(0, 20)
ax4.text(1, 1.40e-11, 'upstream, Tc = 10 Pa', fontsize=8)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)

ax5.plot(high_rvar_midstream_med[:500], high_rvar_midstream_emedian_tc10[:500], 'x', color='#D9B39B', alpha=0.25)
ax5.plot(low_rvar_midstream_med[:500], low_rvar_midstream_emedian_tc10[:500], 'o', color='navy', mfc='none', alpha=0.25)
ax5.set_ylim(0, 2.0e-11)
ax5.set_xlim(0, 100)
ax5.text(5, 1.85e-11, 'midstream, Tc = 10 Pa', fontsize=8)
ax5.set_xlabel('Median discharge in nonsteady method (m$^{3}$ s$^{-1}$)')

ax6.plot(high_rvar_outlet_med[:500], high_rvar_outlet_emedian_tc10[:500], 'x', color='#D9B39B', alpha=0.25)
ax6.plot(low_rvar_outlet_med[:500], low_rvar_outlet_emedian_tc10[:500], 'o', color='navy', mfc='none', alpha=0.25)
ax6.text(5, 2.75e-11, 'outlet, Tc = 10 Pa', fontsize=8)
ax6.set_xlim(0, 100)
ax6.set_ylim(0, 3.0e-11)

plt.show()
