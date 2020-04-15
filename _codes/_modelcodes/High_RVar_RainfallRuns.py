#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
High_RVar_RainfallRuns.py

This code can be used to simulate a rainfall time series, calculating discharge
and erosion rates until 20 m total rainfall depth is reached, using a
nonsteady hydrology method.

It uses parameters from a High RVar Site (see Istanbulluoglu and Bras [2006]).

Please note:
The nonsteady codes can take several days to run a total of 20 m rainfall.

Created by Jordan M. Adams
Updated in April 2020

"""
from landlab.io import read_esri_ascii, write_esri_ascii
from landlab.components import OverlandFlow, PrecipitationDistribution, DepthSlopeProductErosion, LinearDiffuser
import numpy as np
from landlab.utils.depth_dependent_roughness import depth_dependent_mannings_n
from matplotlib import pyplot as plt
from landlab.plot import imshow_grid
import time

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
sample_das = [201, 24048, 14834, 14268, 12097, 8035, 5022] # sample drainage areas

## GENERATE AND SET PRECIPITATION TIME SERIES
total_t = 1000.* 365.25 * 24               # 1,000 years of rainfall generated.
thresh = 0.5                               # Intensity threshold (lower limit, mm/hr)


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
tau_at_sampling_sites = []
#depth_arr = []
modeled_time =[]
storm_events = []

## SET UP COMPONENTS AND FIXED VALUES.
dle = DepthSlopeProductErosion(rmg, k_e = (1e-13), tau_crit=0., a_exp=(1.))
dle1 = DepthSlopeProductErosion(rmg1, k_e = (1e-13), tau_crit=1., a_exp=(1.))
dle5 = DepthSlopeProductErosion(rmg5, k_e = (1e-13), tau_crit=5., a_exp=(1.))
dle10 = DepthSlopeProductErosion(rmg10, k_e = (1e-13), tau_crit=10., a_exp=(1.))

lin_diffuse = LinearDiffuser(rmg, linear_diffusivity=0.1, deposit=False)
lin_diffuse1 = LinearDiffuser(rmg1, linear_diffusivity=0.1, deposit=False)
lin_diffuse5 = LinearDiffuser(rmg5, linear_diffusivity=0.1, deposit=False)
lin_diffuse10 = LinearDiffuser(rmg10, linear_diffusivity=0.1, deposit=False)

uplift_rate = 0.001  # m/yr

# If you want to save the dictionary of parameters for your run, uncomment
# these lines, and the lines before i+=1
# dictionary_data = np.array(['k_e', 'alpha','h_init','mannings_n',
#                 'mean storm dur ','mean interstorm dur', 'mean storm depth',
#                 'mean intensity', 'threshold','tau_crit','diffusion coeff',
#                 'total rain depth ', 'summed erosion depth',
#                 'real modeled time'])
#
# equals = np.array([" = "," = "," = "," = "," = "," = "," = "," = "," = "," = ",
#                    " = "," = "," = "])

total_start_time = time.time()
while depth <= 20000.:

    peak_Q = np.zeros(rmg.number_of_nodes)
    peak_I = np.zeros(rmg.number_of_nodes)
    peak_Itau1 = np.zeros(rmg.number_of_nodes)
    peak_Itau5 = np.zeros(rmg.number_of_nodes)
    peak_Itau10 = np.zeros(rmg.number_of_nodes)

    incision_rate_at_sampling_sites = []
    incision_rate_at_sampling_sites1 = []
    incision_rate_at_sampling_sites5 = []
    incision_rate_at_sampling_sites10 = []

    start = time.time()
    rmg['node']['surface_water__depth'] = np.zeros(rmg.number_of_nodes)
    rmg1['node']['surface_water__depth'] = rmg['node']['surface_water__depth']
    rmg5['node']['surface_water__depth'] = rmg['node']['surface_water__depth']
    rmg10['node']['surface_water__depth'] = rmg['node']['surface_water__depth']

    rmg['link']['mannings_n'] = np.ones(rmg.number_of_links)
    rmg['link']['mannings_n'] *= 0.03

    of = OverlandFlow(rmg, steep_slopes=True, alpha=0.20, mannings_n='mannings_n',
                  h_init=0.00001)
    rmg['node']['water_surface__slope'] = np.zeros(rmg.number_of_nodes)
    #
    rmg['node']['surface_water__discharge'] = np.zeros(rmg.number_of_nodes)

    elapsed_time = 0.
    storm_duration = durations_s[i]
    model_run_time = durations_s[i] + (interstorm_durs[i] * 3600.)

    storm_events.append([durations[i], intensities[i], durations[i] * intensities[i]])
    # these values above are in mm/hr

    rainfall_ms = intensities[i] * 2.777777777778e-7 #convert to m/s

    while (elapsed_time < model_run_time):

        of.dt = of.calc_time_step()

        ## The storm starts when the model starts. While the elapsed time is
        ## less than the storm duration, we add water to the system as rainfall.
        if elapsed_time < (storm_duration):

            of.rainfall_intensity =  rainfall_ms

        ## Then the elapsed time exceeds the storm duration, rainfall ceases.
        else:

            of.rainfall_intensity = 0.0

        # calculating depth depth_dependent_roughness
        depth_dependent_mannings_n(rmg, min_mannings_n=0.03,
                                index_flow_depth=0.05)

        rmg['link']['mannings_n'] = rmg.map_min_of_link_nodes_to_link(
                                                 'mannings_n')

        # calculating discharge
        of.overland_flow(dt=of.dt)

        # calculating Q at nodes using link discharge and node slopes
        node_slope = ((of.water_surface_slope[rmg.links_at_node] *
                     rmg.active_link_dirs_at_node))

        incision_Q = np.abs(of.q * rmg.dx)[rmg.links_at_node]

        rmg['node']['surface_water__discharge'] = (incision_Q[np.arange(len(
                 node_slope)), np.argmax(node_slope, axis=1)])

        peak_Q = np.maximum(peak_Q, np.abs(rmg['node']['surface_water__discharge']))

        node_slope = node_slope.max(axis=1)

        rmg['node']['water_surface__slope'] = node_slope
        rmg1['node']['water_surface__slope'] = rmg['node']['water_surface__slope']
        rmg5['node']['water_surface__slope'] = rmg['node']['water_surface__slope']
        rmg10['node']['water_surface__slope'] = rmg['node']['water_surface__slope']

        dle.erode(of.dt, slope = 'water_surface__slope')
        dle1.erode(of.dt, slope = 'water_surface__slope')
        dle5.erode(of.dt, slope = 'water_surface__slope')
        dle10.erode(of.dt, slope = 'water_surface__slope')

        total_incision_depth += (np.abs(dle.dz))
        total_incision_depth1 += (np.abs(dle1.dz))
        total_incision_depth5 += (np.abs(dle5.dz))
        total_incision_depth10 += (np.abs(dle10.dz))

        modeled_time.append(current_model_time + elapsed_time)

        discharge_at_sampling_sites.append(np.abs(
                rmg['node']['surface_water__discharge'][sample_das]))

        incision_rate_at_sampling_sites.append(np.abs(dle.E[sample_das]))
        incision_rate_at_sampling_sites1.append(np.abs(dle1.E[sample_das]))
        incision_rate_at_sampling_sites5.append(np.abs(dle5.E[sample_das]))
        incision_rate_at_sampling_sites10.append(np.abs(dle10.E[sample_das]))

        tau_at_sampling_sites.append(dle.tau[sample_das])

        elapsed_time += of.dt

        #Some small amount of water may keep the rainfall routing going
        #This can slow down your model big time.
        #Assumes that if Q is less than 0.5 m3/s, and 3x the storm duration
        #has elapsed, the storm is ended (very, very low water = negligible I)
        if (of.rainfall_intensity == 0.0) and (np.abs(of.q[link_to_sample] *
        rmg.dx) <= 0.5) and (elapsed_time > (3.0 * storm_duration)):
            elapsed_time = model_run_time

    end = time.time()

    # This updates model time and will be used to ensure uplift is happening
    # even if rainfall isn't occuring (Tb+Tr)
    old_time = current_model_time
    current_model_time = current_model_time + model_run_time
    new_dt = current_model_time - old_time
    new_dt = new_dt * 3.17098e-8 # seconds to years

    # This lets you know what event number you are at.
    print ("Storm #", i, " Duration: ", round(durations[i], 2),
    " Intensity: ", round(intensities[i], 2), "      Time in loop: ",
    round(end-start, 2), "seconds")

    # Updating total depth of the time series
    depth += (durations[i]*intensities[i])

    # uplifting over the whole storm and interstorm
    rmg['node']['topographic__elevation'][rmg.core_nodes] += uplift_rate * new_dt
    rmg1['node']['topographic__elevation'][rmg1.core_nodes] += uplift_rate * new_dt
    rmg5['node']['topographic__elevation'][rmg5.core_nodes] += uplift_rate * new_dt
    rmg10['node']['topographic__elevation'][rmg10.core_nodes] += uplift_rate * new_dt

    lin_diffuse.run_one_step(dt = new_dt)
    lin_diffuse1.run_one_step(dt = new_dt)
    lin_diffuse5.run_one_step(dt = new_dt)
    lin_diffuse10.run_one_step(dt = new_dt)

    real_time_list.append(round(end - start, 2))

    # Keeping track of erosion rates and discharges at the study locations
    # For each of the tau_c values
    sampled_peak_discharges.append(peak_Q[sample_das])
    sampled_erosion_rates.append(np.amax(np.abs(incision_rate_at_sampling_sites), axis=0))
    sampled_erosion_rates1.append(np.amax(np.abs(incision_rate_at_sampling_sites1), axis=0))
    sampled_erosion_rates5.append(np.amax(np.abs(incision_rate_at_sampling_sites5), axis=0))
    sampled_erosion_rates10.append(np.amax(np.abs(incision_rate_at_sampling_sites10), axis=0))

    # These "i" values represent the events where rainfall hits 1 m, 2.5 m,
    #   5 m, 7.5 m, 10 m and 15 m. The end of the loop is 20 m.
    # This can be used to save your values
    #if i == 70 or 198 or 414 or 609 or 802 or 1149:
    #
    #     dictionary_data2 = np.array([dle.k_e, of.alpha, of.h_init, min(rmg.at_node[
    #     'mannings_n']), PD.mean_storm_duration, PD.mean_interstorm_duration,
    #     PD.mean_storm_depth, PD.mean_intensity, thresh, 'variable tau crit', 0.1,
    #     depth, sum(total_incision_depth), sum(real_time_list)])
    #
    #     g = zip(dictionary_data, equals, dictionary_data2)
        # with open('/Users/Jordan/Desktop/ID_Results_a1/20mGeneral/param_summary.txt', 'w') as f:
        #     for tuple in g:
        #         f.write('%s %s %s\n' % tuple)

    i+=1
