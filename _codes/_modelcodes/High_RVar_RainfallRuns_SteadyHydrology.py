#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
High_RVar_RainfallRuns_SteadyHydrology.py

This code can be used to simulate a rainfall time series, calculating discharge
and erosion rates until 20 m total rainfall depth is reached, using a steady
discharge method (Q = PA).

It uses parameters from a High RVar Site (see Istanbulluoglu and Bras [2006]).

Created by Jordan M. Adams
Updated in April 2020

"""
from landlab.io import read_esri_ascii, write_esri_ascii
from landlab.components import FlowRouter, PrecipitationDistribution, LinearDiffuser, StreamPowerEroder
import numpy as np
from matplotlib import pyplot as plt
from landlab.plot import imshow_grid
import time
import copy

# SET UP GRID
(rmg, z) = read_esri_ascii('_input/topographic_elevation_4.asc', name='topographic__elevation')
(rmg, A) = read_esri_ascii('_input/drainage_area_4.asc', name='drainage_area')
(rmg1, z1) = read_esri_ascii('_input/topographic_elevation_4.asc', name='topographic__elevation')
(rmg5, z5) = read_esri_ascii('_input/topographic_elevation_4.asc', name='topographic__elevation')
(rmg10, z10) = read_esri_ascii('_input/topographic_elevation_4.asc', name='topographic__elevation')

# SET BOUNDARY CONDITIONS
rmg.set_closed_boundaries_at_grid_edges(True, True, True, True)
rmg.status_at_node[1] = 1

rmg1.set_closed_boundaries_at_grid_edges(True, True, True, True)
rmg1.status_at_node[1] = 1

rmg5.set_closed_boundaries_at_grid_edges(True, True, True, True)
rmg5.status_at_node[1] = 1

rmg10.set_closed_boundaries_at_grid_edges(True, True, True, True)
rmg10.status_at_node[1] = 1

link_to_sample = 200
sample_das = [201, 24048, 14834, 14268, 12097, 8035, 5022] # Sampling sites
peak_q = []
ORIGINAL_TOPO = copy.deepcopy(z)

## GENERATE AND SET PRECIPITATION TIME SERIES
total_t = 1000.* 365.25 * 24               # 1,000 years of rainfall generated.
thresh = 0.5                               # Intensity threshold (lower limit)

PD = PrecipitationDistribution(mean_storm_duration = 10.75,
                               mean_interstorm_duration=433.58,
                               mean_storm_depth = 9.62,
                               total_t=total_t)

storm_arr = np.array(PD.get_storm_time_series())
intensity_threshold, = np.where(storm_arr[:, 2] > thresh)
interstorm_durs = (storm_arr[intensity_threshold][:,0][np.arange(1,
                        len(storm_arr[intensity_threshold]) - 1)] -
                        storm_arr[intensity_threshold][:, 1][np.arange(0,
                        len(storm_arr[intensity_threshold]) - 2)])

durations = (storm_arr[intensity_threshold][:,1] -
             storm_arr[intensity_threshold][:,0])

durations_s = [x * 3600. for x in durations]

intensities = storm_arr[intensity_threshold][:,2]

## SET INITIAL FIELDS, ARRAYS AND ITERATORS
rmg['node']['topographic__elevation'] = z
rmg1['node']['topographic__elevation'] = z1
rmg5['node']['topographic__elevation'] = z5
rmg10['node']['topographic__elevation'] = z10

total_incision_depth = np.zeros(rmg.number_of_nodes)
total_incision_depth1 = np.zeros(rmg.number_of_nodes)
total_incision_depth5 = np.zeros(rmg.number_of_nodes)
total_incision_depth10 = np.zeros(rmg.number_of_nodes)

depth = 0
i = 0
current_model_time = 0.

## SET LISTS FOR SAVING DATA
real_time_list=[]
discharge_at_sampling_sites = []
sampled_peak_discharges = []
sampled_erosion_rates = []
sampled_erosion_rates1 =[]
sampled_erosion_rates5 = []
sampled_erosion_rates10 = []

modeled_time =[]
storm_events = []

a = 1

if a == 1:
    m = 0.3
    n = 0.7
    k_tau = 8.357150818380437e-15

tc0 = (k_tau * (0.)**a)
tc1 = (k_tau * (10.)**a)
tc5 = (k_tau * (25.)**a)
tc10 = (k_tau * (50.)**a)

fr = FlowRouter(rmg, method='D4')
fr1 = FlowRouter(rmg1, method='D4')
fr5 = FlowRouter(rmg5, method='D4')
fr10 = FlowRouter(rmg10, method='D4')
dle = StreamPowerEroder(rmg, K_sp = (1e-10), threshold_sp=float(tc0), use_Q='water__discharge', m_sp = m, n_sp = n)#(3.4e-11)
dle1 = StreamPowerEroder(rmg1,K_sp = (1e-10), threshold_sp=float(tc1), use_Q='water__discharge', m_sp = m, n_sp = n)#(3.4e-11)
dle5 = StreamPowerEroder(rmg5, K_sp = (1e-10), threshold_sp=float(tc5), use_Q='water__discharge', m_sp = m, n_sp = n)#(3.4e-11)
dle10 = StreamPowerEroder(rmg10, K_sp = (1e-10), threshold_sp=float(tc10), use_Q='water__discharge', m_sp = m, n_sp = n)#(3.4e-11)

lin_diffuse = LinearDiffuser(rmg, linear_diffusivity=0.1, deposit=False)
lin_diffuse1 = LinearDiffuser(rmg1, linear_diffusivity=0.1, deposit=False)
lin_diffuse5 = LinearDiffuser(rmg5, linear_diffusivity=0.1, deposit=False)
lin_diffuse10 = LinearDiffuser(rmg10, linear_diffusivity=0.1, deposit=False)

uplift_rate = 0.001

dictionary_data = np.array(['a', 'm','n','mean storm depth',
                'mean storm dur ','mean interstorm dur'])

equals = np.array([" = "," = "," = "," = "," = "," = "])

intensities_m = [x / 1000. for x in intensities]
erosion_depths = np.zeros(rmg.number_of_nodes)
while i < 1558: #depth <= 20000.:

    peak_Q = np.zeros(rmg.number_of_nodes)
    peak_I = np.zeros(rmg.number_of_nodes)
    peak_Itau1 = np.zeros(rmg.number_of_nodes)
    peak_Itau5 = np.zeros(rmg.number_of_nodes)
    peak_Itau10 = np.zeros(rmg.number_of_nodes)

    incision_rate_at_sampling_sites = []
    incision_rate_at_sampling_sites1 = []
    incision_rate_at_sampling_sites5 = []
    incision_rate_at_sampling_sites10 = []

    storm_duration = durations_s[i]
    model_run_time = durations_s[i] + (interstorm_durs[i] * 3600.)

    storm_events.append([durations[i], intensities[i], durations[i] * intensities[i]])
    dt = float(durations_s[i])
    uplift_dt = float(durations_s[i] + (interstorm_durs[i] * 3600.)) * 3.17098e-8

    rainfall_ms = intensities[i] * 2.777777777778e-7

    fr = FlowRouter(rmg, method='D4', runoff_rate= float(rainfall_ms))
    fr1 = FlowRouter(rmg1, method='D4', runoff_rate= float(rainfall_ms))
    fr5 = FlowRouter(rmg5, method='D4', runoff_rate= float(rainfall_ms))
    fr10 = FlowRouter(rmg10, method='D4', runoff_rate= float(rainfall_ms))

    fr.run_one_step()
    fr1.run_one_step()
    fr5.run_one_step()
    fr10.run_one_step()

    last_z = copy.deepcopy(z)
    last_z1 = copy.deepcopy(z1)
    last_z5 = copy.deepcopy(z5)
    last_z10 = copy.deepcopy(z10)

    rmg, z, _ = dle.erode(rmg, dt=dt)
    rmg1,z1, e1 = dle1.erode(rmg1, dt=dt)
    rmg5,z5, e5 = dle5.erode(rmg5, dt=dt)
    rmg10,z10, e10 = dle10.erode(rmg10, dt=dt)

    total_incision_depth += last_z-z
    total_incision_depth1 += last_z1-z1
    total_incision_depth5 += last_z5-z5
    total_incision_depth10 += last_z10-z10

    rmg['node']['topographic__elevation'][rmg.core_nodes] += (uplift_rate * uplift_dt)
    rmg1['node']['topographic__elevation'][rmg1.core_nodes] += (uplift_rate * uplift_dt)
    rmg5['node']['topographic__elevation'][rmg5.core_nodes] += (uplift_rate * uplift_dt)
    rmg10['node']['topographic__elevation'][rmg10.core_nodes] += (uplift_rate * uplift_dt)

    lin_diffuse.run_one_step(dt = uplift_dt)
    lin_diffuse1.run_one_step(dt = uplift_dt)
    lin_diffuse5.run_one_step(dt = uplift_dt)
    lin_diffuse10.run_one_step(dt = uplift_dt)

    peak_q.append(rmg.at_node['surface_water__discharge'][sample_das])

    print ("Storm #", i, " Duration: ", round(durations[i], 2),
          " Intensity: ", round(intensities[i], 2))

    depth += (durations[i]*intensities[i])

    # if i == 70 or 198 or 414 or 609 or 802 or 1149:
    #         #these i values represent storms that hit 1 m, 2.5 m, 5 m, 7.5 m
    #         # 10 m, 12.5 m and 15 m.
    #     dictionary_data2 = (np.array([a, m, n, PD.mean_storm_depth,
    #                                 PD.mean_storm_duration,
    #                                 PD.mean_interstorm_duration]))
    #
    #     g = zip(dictionary_data, equals, dictionary_data2)

    i+=1
