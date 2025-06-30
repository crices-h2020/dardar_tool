#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# In[]:
# ########################################################################
# Plotting tools
# ########################################################################

"""
Created on Thu Oct 31 14:23:10 2024

@author: vaisanea
"""
# In[]:
# ########################################################################
# Pakages
# ########################################################################

import os
import sys

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import numpy as np

current_path = os.path.abspath(__file__)
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(current_path)))
working_dir = os.path.dirname(root_dir)

source_path = os.path.join(root_dir, "source")
sys.path.append(source_path)
import functions_dardar_validation as functions

# In[]:
# ########################################################################
# Functions
# ########################################################################


def de_mask(tim, lat, lon, profile_val, profile_height):
    if isinstance(tim, np.ma.MaskedArray):
        tim = tim.data
    if isinstance(lat, np.ma.MaskedArray):
        lat = lat.filled(np.nan)
    if isinstance(lon, np.ma.MaskedArray):
        lon = lon.filled(np.nan)
    if isinstance(profile_val, np.ma.MaskedArray):
        profile_val = profile_val.filled(np.nan)
    if isinstance(profile_height, np.ma.MaskedArray):
        profile_height = profile_height.filled(np.nan)
    return tim, lat, lon, profile_val, profile_height


def plot_profile(fig, ax, tim, lat, lon, profile_val, profile_height, colormap="viridis", long_ax_loc=-0.2, time_ax_loc=-0.4):
    # Pre process
    tim, lat, lon, profile_val, profile_height = de_mask(tim, lat, lon, profile_val, profile_height)
    nan_mask = np.logical_not(np.isnan(profile_height))
    profile_val = profile_val[:, nan_mask]
    profile_height = profile_height[nan_mask]
    if profile_height.shape[0] == profile_val.shape[1]:
        profile_val = np.swapaxes(profile_val, axis1=0, axis2=1)
    # Plotting
    plt.figure(fig)
    #fig, ax = plt.subplots()
    cax = ax.imshow(profile_val, cmap=colormap, aspect="auto")
    ax.set_xticks([])
    ax.set_yticks([])
    fig.colorbar(cax, ax=ax)
    # y-axis
    ax_height = ax.twinx()
    ax_height.set_ylim(profile_height[0], profile_height[-1])
    ax_height.set_ylabel("Height")
    ax_height.spines["left"].set_position(("axes", 0))
    ax_height.yaxis.set_ticks_position("left")
    ax_height.yaxis.set_label_position("left")
    # x-axis
    ax_lat = ax.twiny()
    ax_lat.set_xlim(lat[0], lat[-1])
    ax_lat.set_xlabel("Latitude")
    ax_lat.spines["bottom"].set_position(("axes", 0.0))
    ax_lat.xaxis.set_ticks_position("bottom")
    ax_lat.xaxis.set_label_position("bottom")
    ax_lon = ax.twiny()
    ax_lon.set_xlim(lon[0], lon[-1])
    ax_lon.set_xlabel("Longitude")
    ax_lon.spines["bottom"].set_position(("axes", long_ax_loc))
    ax_lon.xaxis.set_ticks_position("bottom")
    ax_lon.xaxis.set_label_position("bottom")
    ax_time = ax.twiny()
    ax_time.set_xlim(tim[0], tim[-1])
    ax_time.set_xlabel("Time")
    ax_time.spines["bottom"].set_position(("axes", time_ax_loc))
    ax_time.xaxis.set_ticks_position("bottom")
    ax_time.xaxis.set_label_position("bottom")
    ax_time.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M')) # <--- seconds can be added
    for label in ax_time.get_xticklabels():
        label.set_rotation(-80)
    return fig, ax, ax_height, ax_lat, ax_lon, ax_lon


def plot_flight(fig, ax, tim, lat, lon, north_pad, south_pad, west_max, east_pad):
    if not hasattr(ax, 'projection') or not isinstance(ax.projection, ccrs.PlateCarree):
        fig.delaxes(ax)
        ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
    lat_min, lat_max = lat.min(), lat.max()
    lon_min, lon_max = lon.min(), lon.max()
    extent = [lon_min - east_pad, lon_max + west_max, lat_min - south_pad, lat_max + north_pad]
    ax.set_extent(extent)
    for i in range(len(tim)):
        plt.plot(lon[i], lat[i], marker='o', color='red', markersize=1, transform=ccrs.PlateCarree())
        if (i == 0) or (i == (len(tim)-1)):
            time_formating = tim[i] # <----
            plt.text(lon[i], lat[i], time_formating.astype("datetime64[m]"), fontsize=18, transform=ccrs.PlateCarree())
    ax.coastlines()
    gl = ax.gridlines(draw_labels=True)
    gl.right_labels = False
    gl.top_labels = False
    return fig, ax


def plot_box_on_fig(ax, lat_1, lat_2, lon_1, lon_2):
    # Draw a box on the map by plotting the four corner points
    ax.plot([lon_1, lon_1, lon_2, lon_2, lon_1],
            [lat_1, lat_2, lat_2, lat_1, lat_1],
            color='blue', linewidth=2, transform=ccrs.PlateCarree())


def plot_flight_box(fig, ax, tim, lat, lon, zoom_padding=5):
    #lat_north_limit, lat_south_limit, lon_west_limit, lon_east_limit, 
    #north_limit, south_limit, west_limit, east_limit, 
    if isinstance(zoom_padding, int):
        north_pad = south_pad = west_pad = east_pad = zoom_padding
    elif isinstance(zoom_padding, (tuple, list)) and len(zoom_padding) == 4:
        north_pad = zoom_padding[0]
        south_pad = zoom_padding[1]
        west_pad = zoom_padding[2]
        east_pad = zoom_padding[3]
    else:
        raise ValueError("zoom_padding must be an integer or a tuple of four values")
    fig, ax = plot_flight(fig, ax, tim, lat, lon, north_pad, south_pad, west_pad, east_pad)
    #lon_min = max(np.max(lon), lon_east_limit)
    #lon_max = min(np.min(lon), lon_west_limit)
    #lat_min = min(np.min(lat), lat_south_limit)
    #lat_max = max(np.max(lat), lat_north_limit)
    lon_min = np.max(lon)
    lon_max = np.min(lon)
    lat_min = np.min(lat)
    lat_max = np.max(lat)
    extent = [lon_min + east_pad, lon_max - west_pad, lat_min - south_pad, lat_max + north_pad]
    ax.set_extent(extent)
    #plot_box_on_fig(ax, lat_north_limit, lat_south_limit, lon_west_limit, lon_east_limit)
    return fig, ax


def plot_flight_cumulative_sum(fig, ax, tim, lat, lon, cumulative_sum_values, north_pad, south_pad, west_max, east_pad, cmap_name="viridis"):
    if not hasattr(ax, 'projection') or not isinstance(ax.projection, ccrs.PlateCarree):
        fig.delaxes(ax)
        ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
    lat_min, lat_max = lat.min(), lat.max()
    lon_min, lon_max = lon.min(), lon.max()
    extent = [lon_min - east_pad, lon_max + west_max, lat_min - south_pad, lat_max + north_pad]
    ax.set_extent(extent)
    norm = mcolors.Normalize(vmin=np.min(cumulative_sum_values), vmax=np.max(cumulative_sum_values))
    cmap = cm.get_cmap(cmap_name)
    for i in range(len(tim)):
        color = cmap(norm(cumulative_sum_values[i]))
        plt.plot(lon[i], lat[i], marker='o', color=color, markersize=1, transform=ccrs.PlateCarree())
        if (i == 0) or (i == (len(tim)-1)):
            time_formating = tim[i] # <----
            plt.text(lon[i], lat[i], time_formating.astype("datetime64[m]"), fontsize=18, transform=ccrs.PlateCarree())

    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation='vertical', shrink=0.7, pad=0.05)
    cbar.set_label('Cumulative Sum')


    ax.coastlines()
    gl = ax.gridlines(draw_labels=True)
    gl.right_labels = False
    gl.top_labels = False
    return fig, ax


def plot_flight_box_cumulative_sum(fig, ax, tim, lat, lon, cum_sum_vals, zoom_padding=5, colormap="viridis"):
    #lat_north_limit, lat_south_limit, lon_west_limit, lon_east_limit, 
    #north_limit, south_limit, west_limit, east_limit, 
    if isinstance(zoom_padding, int):
        north_pad = south_pad = west_pad = east_pad = zoom_padding
    elif isinstance(zoom_padding, (tuple, list)) and len(zoom_padding) == 4:
        north_pad = zoom_padding[0]
        south_pad = zoom_padding[1]
        west_pad = zoom_padding[2]
        east_pad = zoom_padding[3]
    else:
        raise ValueError("zoom_padding must be an integer or a tuple of four values")
    fig, ax = plot_flight_cumulative_sum(fig, ax, tim, lat, lon, cum_sum_vals, north_pad, south_pad, west_pad, east_pad, colormap)
    #lon_min = max(np.max(lon), lon_east_limit)
    #lon_max = min(np.min(lon), lon_west_limit)
    #lat_min = min(np.min(lat), lat_south_limit)
    #lat_max = max(np.max(lat), lat_north_limit)
    lon_min = np.max(lon)
    lon_max = np.min(lon)
    lat_min = np.min(lat)
    lat_max = np.max(lat)
    extent = [lon_min + east_pad, lon_max - west_pad, lat_min - south_pad, lat_max + north_pad]
    ax.set_extent(extent)
    #plot_box_on_fig(ax, lat_north_limit, lat_south_limit, lon_west_limit, lon_east_limit)
    return fig, ax


def plot_combined_model_comparison(tim, lat, lon, m_profile_vals, profile_vals, vert_profile, zoom_padding=5, flight_name_str="", colormap="viridis"):
    profile_vals_diff = m_profile_vals - profile_vals
    cumsum_prof = np.sum(profile_vals, axis=1)
    cumsum_m_prof = np.sum(m_profile_vals, axis=1)
    cumsum_diff = cumsum_m_prof - cumsum_prof
    megafig = plt.figure(figsize=(18,9), constrained_layout=True)
    megafig.suptitle(f"Flight Path Cumulative Sum and Vertical Profile Comparison, Flight:{flight_name_str}", fontsize=30)
    subfigs = megafig.subfigures(nrows=1, ncols=2, width_ratios=(0.49, 0.51))
    fig1 = subfigs[0]
    ax1 = fig1.subplots(nrows=1, ncols=1)
    _ = plot_flight_box_cumulative_sum(fig1, ax1, tim, lat, lon, cumsum_diff, zoom_padding, colormap=colormap)
    fig2 = subfigs[1]
    ax2 = fig2.subplots()
    _ = plot_profile(fig2, ax2, tim, lat, lon, profile_vals_diff, vert_profile, colormap, long_ax_loc=-0.09, time_ax_loc=-0.18)
    plt.show()


def plot_combined(tim, lat, lon, profile_vals, vert_profile, zoom_padding=5, flight_name_str="", colormap="viridis"):
    megafig = plt.figure(figsize=(18,9), constrained_layout=True)
    megafig.suptitle(f"Flight Path and Vertical Profile, Flight:{flight_name_str}", fontsize=30)
    subfigs = megafig.subfigures(nrows=1, ncols=2, width_ratios=(0.49, 0.51))
    fig1 = subfigs[0]
    ax1 = fig1.subplots(nrows=1, ncols=1)
    _ = plot_flight_box(fig1, ax1, tim, lat, lon, zoom_padding)
    fig2 = subfigs[1]
    ax2 = fig2.subplots()
    _ = plot_profile(fig2, ax2, tim, lat, lon, profile_vals, vert_profile, colormap, long_ax_loc=-0.09, time_ax_loc=-0.18)
    plt.show()
    

def separate_flights(tim, lat, lon, vertical_profile, flight_separation_time=1, flight_separation_unit="m"):
    time_max_separation = np.timedelta64(flight_separation_time, flight_separation_unit)
    time_differences = np.diff(tim)
    time_over_mask = time_differences > time_max_separation
    ind_locs = np.where(time_over_mask)[0] + 1
    tim_groups = np.split(tim, ind_locs)
    lat_groups = np.split(lat, ind_locs)
    lon_groups = np.split(lon, ind_locs)
    profile_groups = np.split(vertical_profile, ind_locs, axis=0)
    return tim_groups, lat_groups, lon_groups, profile_groups


def plot_comparison(tim1, lat1, lon1, prof_vals1, vert_prof1, tim2, lat2, lon2, prof_vals2, vert_prof2, flight_name_str="", colormap="viridis"):
    megafig = plt.figure(figsize=(18,9), constrained_layout=True)
    megafig.suptitle(f"Vertical Profile Comparison, Flight:{flight_name_str}", fontsize=30)
    subfigs = megafig.subfigures(nrows=1, ncols=2, width_ratios=(0.49, 0.51))
    fig1 = subfigs[0]
    ax1 = fig1.subplots(nrows=1, ncols=1)
    ax1.set_title("User vertical profile")
    _ = plot_profile(fig1, ax1, tim1, lat1, lon1, prof_vals1, vert_prof1, colormap, long_ax_loc=-0.09, time_ax_loc=-0.18)
    fig2 = subfigs[1]
    ax2 = fig2.subplots()
    ax2.set_title("DARDAR vertical profile")
    _ = plot_profile(fig2, ax2, tim2, lat2, lon2, prof_vals2, vert_prof2, colormap, long_ax_loc=-0.09, time_ax_loc=-0.18)
    plt.show()


def plot_data_comparison(tim1, lat1, lon1, vert_alt1, vert_prof1, tim2, lat2, lon2, vert_alt2, vert_prof2, flight_separation_time=1, flight_separation_unit="m", colormap="viridis"):
    final_time1, final_lat1, final_lon1, final_profile1 = separate_flights(tim1, lat1, lon1, vert_prof1, flight_separation_time, flight_separation_unit)
    final_time2, final_lat2, final_lon2, final_profile2 = separate_flights(tim2, lat2, lon2, vert_prof2, flight_separation_time, flight_separation_unit)
    for ind, (i1_tim, i1_lat, i1_lon, i1_prof, i2_tim, i2_lat, i2_lon, i2_prof) in enumerate(zip(final_time1, final_lat1, final_lon1, final_profile1, final_time2, final_lat2, final_lon2, final_profile2)):
        plot_comparison(i1_tim, i1_lat, i1_lon, i1_prof, vert_alt1, i2_tim, i2_lat, i2_lon, i2_prof, vert_alt2, str(ind+1), colormap=colormap)


# In[]:
# ########################################################################
# Main Function
# ########################################################################

def plot_saved_data(filename, zoom_map_in_out=5, flight_separation_time=1, flight_separation_unit="m", colormap="viridis"):
    tim, lat, lon, profile, alt_profile = functions.read_from_results_file_time_conversion(filename)
    final_time, final_lat, final_lon, final_profile = separate_flights(tim, lat, lon, profile, flight_separation_time, flight_separation_unit)
    for ind, (i_time, i_lat, i_lon, i_profile) in enumerate(zip(final_time, final_lat, final_lon, final_profile)):
        plot_combined(i_time, i_lat, i_lon, i_profile, alt_profile, zoom_map_in_out, str(ind+1), colormap=colormap)


def plot_local_data(time, latitude, longitude, vertical_altitudes, vertical_profile, zoom_map_in_out=5, flight_separation_time=1, flight_separation_unit="m", colormap="viridis"):
    final_time, final_lat, final_lon, final_profile = separate_flights(time, latitude, longitude, vertical_profile, flight_separation_time, flight_separation_unit)
    for ind, (i_time, i_lat, i_lon, i_profile) in enumerate(zip(final_time, final_lat, final_lon, final_profile)):
        plot_combined(i_time, i_lat, i_lon, i_profile, vertical_altitudes, zoom_map_in_out, str(ind+1), colormap=colormap)


def plot_model_vs_dardar_data_vertical(m_tim, m_lat, m_lon, m_vertical_altitudes, m_vertical_profile, vertical_profile, flight_separation_time=1, flight_separation_unit="m", zoom_map_in_out=5, colormap="viridis"):
    final_time, final_lat, final_lon, m_final_profile = separate_flights(m_tim, m_lat, m_lon, m_vertical_profile, flight_separation_time, flight_separation_unit)
    _, _, _, final_profile = separate_flights(m_tim, m_lat, m_lon, vertical_profile, flight_separation_time, flight_separation_unit)    
    for ind, (i_time, i_lat, i_lon, m_i_profile, i_profile) in enumerate(zip(final_time, final_lat, final_lon, m_final_profile, final_profile)):
        plot_combined_model_comparison(i_time, i_lat, i_lon, m_i_profile, i_profile, m_vertical_altitudes, zoom_map_in_out, str(ind+1), colormap=colormap)


# def plot_model_vs_dardar_data_vertical(m_tim, m_lat, m_lon, m_vertical_altitudes, m_vertical_profile, vertical_profile, flight_separation_time=1, flight_separation_unit="m", colormap="viridis"):
#     prof_vals = m_vertical_profile - vertical_profile
#     final_time, final_lat, final_lon, final_profile = separate_flights(m_tim, m_lat, m_lon, prof_vals, flight_separation_time, flight_separation_unit)
#     for ind, (i_time, i_lat, i_lon, i_profile) in enumerate(zip(final_time, final_lat, final_lon, final_profile)):
#         megafig = plt.figure(figsize=(18,9), constrained_layout=True)
#         #megafig.suptitle(f"Difference = Model - DARDAR data Vertical Profile, Flight:{ind+1}", fontsize=30)
#         ax1 = megafig.subplots()
#         ax1.set_title(f"Difference = Model - DARDAR data Vertical Profile, Flight:{ind+1}")
#         _ = plot_profile(megafig, ax1, i_time, i_lat, i_lon, i_profile, m_vertical_altitudes, colormap, long_ax_loc=-0.09, time_ax_loc=-0.18)
#         plt.show()


# def plot_model_vs_dardar_data_cumulative_sum_trajectory(m_tim, m_lat, m_lon, m_vertical_altitudes, m_vertical_profile, vertical_profile, flight_separation_time=1, flight_separation_unit="m", zoom_map_in_out=5, colormap="viridis"):
#     prof_vals = m_vertical_profile - vertical_profile
#     final_time, final_lat, final_lon, final_profile = separate_flights(m_tim, m_lat, m_lon, prof_vals, flight_separation_time, flight_separation_unit)
#     for ind, (i_time, i_lat, i_lon, i_profile) in enumerate(zip(final_time, final_lat, final_lon, final_profile)):
#         cumulative_profile_vals = np.sum(i_profile, axis=1)
#         megafig = plt.figure(figsize=(18,9), constrained_layout=True)
#         #megafig.suptitle("Difference = Model - DARDAR data Vertical Profile", fontsize=30)
#         ax1 = megafig.subplots()
#         ax1.set_title(f"Difference = Model - DARDAR data Vertical Profile, Flight:{ind+1}")
#         plot_flight_box_cumulative_sum(megafig, ax1, i_time, i_lat, i_lon, cumulative_profile_vals, zoom_padding=zoom_map_in_out, colormap="viridis")
#         plt.show()


# In[]:
# ########################################################################
# Run Main
# ########################################################################

if __name__ == "__main__":
    print("Test")
    #plot_local_data(tim, lat, lon, vert_profile, profile_vals)



