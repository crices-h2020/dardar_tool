#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# In[]:
# ########################################################################
# 
# ########################################################################
"""
Created on Fri Jun  6 13:32:26 2025

@author: vaisanea
"""
# In[]:
# ########################################################################
# Pakages
# ########################################################################

import glob
import warnings
import sys
import os

import numpy as np


# In[]:
# ########################################################################
# Paths
# ########################################################################

current_path = os.path.abspath(__file__)
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(current_path)))
working_dir = os.path.dirname(root_dir)

default_dardar_path = os.path.join(working_dir, "validation_data")

# In[]:
# ########################################################################
# Functions
# ########################################################################


def check_mask(vert_profile, vert_vals, mask):
    valid_indices = np.where(~mask)[0]
    last_valid_index = valid_indices[-1]
    return vert_profile[:last_valid_index + 1], vert_vals[:, :last_valid_index + 1]


def round_to_resolution(array, resolution):
    return np.round(array / resolution) * resolution


def create_tag_array(values):
    changes = np.diff(values) != 0
    tag_array = np.cumsum(np.insert(changes, 0, True))
    return tag_array


def calc_groups(lat_tag, lon_tag, time_tag=None):
    if time_tag is None:
        tags = np.vstack([lat_tag, lon_tag])
    else:
        tags = np.vstack([lat_tag, lon_tag, time_tag])
    sorted_indices = np.lexsort(tags[::-1])
    sorted_tags = tags[:, sorted_indices]
    diffs = np.any(np.diff(sorted_tags, axis=1) != 0, axis=0)
    group_starts = np.insert(diffs, 0, True)
    group_ids = np.cumsum(group_starts)
    return group_ids


def round_datetime_array(array, resolution):
    orig_unit = array.dtype.name.split('[')[1].strip(']')
    arr_ns = array.astype('datetime64[ns]').astype('int64')
    res_ns = resolution.astype('timedelta64[ns]').astype('int64')
    rounded_ns = (np.round(arr_ns / res_ns).astype('int64') * res_ns)
    final = rounded_ns.astype('datetime64[ns]').astype(f'datetime64[{orig_unit}]')
    return final


def coord_scaling(tags, data):
    tags_zero_based = tags - tags.min()
    sum_per_group = np.bincount(tags_zero_based, weights=data)
    count_per_group = np.bincount(tags_zero_based)
    average_per_group = sum_per_group / count_per_group
    return average_per_group


def pick_times(tags, time_data):
    reference_time = time_data.min()
    time_seconds = (time_data - reference_time) / np.timedelta64(1, 's')
    tags_zero_based = tags - tags.min()
    sum_time_per_group = np.bincount(tags_zero_based, weights=time_seconds)
    count_per_group = np.bincount(tags_zero_based)
    average_time_per_group = sum_time_per_group / count_per_group
    average_times = reference_time + average_time_per_group.astype('timedelta64[s]')
    return average_times


def mean_profile_vals(tags, data):
    unique_tags, tags_zero_based = np.unique(tags, return_inverse=True)
    num_groups = len(unique_tags)
    sum_per_group = np.zeros((num_groups, data.shape[1]), dtype=data.dtype)
    count_per_group = np.zeros(num_groups, dtype=int)
    np.add.at(sum_per_group, tags_zero_based, data)
    np.add.at(count_per_group, tags_zero_based, 1)
    mean_per_group = sum_per_group / count_per_group[:, None]
    return mean_per_group


def mean_profile_height(tags, data):
    unique_tags, tags_zero_based = np.unique(tags, return_inverse=True)
    num_tags = len(unique_tags)
    x = data.shape[0]
    sum_per_tag = np.zeros((x, num_tags), dtype=data.dtype)
    count_per_tag = np.zeros((num_tags,), dtype=int)
    np.add.at(count_per_tag, tags_zero_based, 1)
    np.add.at(sum_per_tag, (np.arange(x)[:, None], tags_zero_based[None, :]), data)
    mean_data = sum_per_tag / count_per_tag[None, :]
    return mean_data


# In[]:
# ########################################################################
# Main
# ########################################################################


def calculate_cumulative_values_for_column(profile_values):
    cumulative_values = np.sum(profile_values, axis=1)
    return cumulative_values


def scale_data_resolution(time_array, latitude_array, longitude_array, profile_values, vertical_profile, latitude_resolution=None, longitude_resolution=None, height_resolution=None, time_resolution=None):
    vertical_height_profile, profile_values = check_mask(vertical_profile.data, profile_values, vertical_profile.mask)
    if height_resolution != None:
        new_reso_height = round_to_resolution(vertical_height_profile, height_resolution)
        height_tags = create_tag_array(new_reso_height)
        vertical_height_profile = coord_scaling(height_tags, new_reso_height)
        profile_values = mean_profile_height(height_tags, profile_values)
    if latitude_resolution == None:
        latitude_resolution = np.min(np.abs(np.diff(latitude_array)))
    if longitude_resolution == None:
        longitude_resolution = np.min(np.abs(np.diff(longitude_array)))
    new_reso_latitudes = round_to_resolution(latitude_array, latitude_resolution)
    new_reso_longitudes = round_to_resolution(longitude_array, longitude_resolution)
    lat_tag_array = create_tag_array(new_reso_latitudes)
    lon_tag_array = create_tag_array(new_reso_longitudes)
    if time_resolution != None:
        new_reso_time = round_datetime_array(time_array, time_resolution)
        tim_tag_array = create_tag_array(new_reso_time)
        id_arr = calc_groups(lat_tag_array, lon_tag_array, tim_tag_array)
        tim_picked = pick_times(id_arr, new_reso_time)
        time_array = tim_picked
    else:
        id_arr = calc_groups(lat_tag_array, lon_tag_array)
    lat_scaled = coord_scaling(id_arr, new_reso_latitudes)
    lon_scaled = coord_scaling(id_arr, new_reso_longitudes)
    profile_scaled = mean_profile_vals(id_arr, profile_values)
    return time_array, lat_scaled, lon_scaled, profile_scaled, vertical_height_profile


# In[]:
# ########################################################################
# Run Main
# ########################################################################

if __name__ == "__main__":
    
    try:
        latitude_array = lat
        longitude_array = lon
        time_array = tim
        profile_values = profile_vals
        vertical_height_profile = vert_profile
    except NameError:
        print("Lataa nää jostain! vaikke test.py!")

    latitude_resolution = 0.1
    longitude_resolution = 0.1
    time_resolution = np.timedelta64(1, "m")
    height_resolution = 1500

    t, la, lo, pr, h = scale_data_resolution(time_array, latitude_array, longitude_array, profile_values, vertical_height_profile, latitude_resolution, longitude_resolution, height_resolution, time_resolution)

    import matplotlib.pyplot as plt

    x_min = h[0]
    x_max = h[-1]
    y_min = latitude_array[0]
    y_max = latitude_array[-1]
    plt.imshow(profile_values, aspect="auto", extent=[x_min, x_max, y_min, y_max])
    #plt.imshow(profile_values, aspect="auto")
    plt.colorbar()
    plt.show()

    x_min = h[0]
    x_max = h[-1]
    y_min = latitude_array[0]
    y_max = latitude_array[-1]
    plt.imshow(pr, aspect="auto", extent=[x_min, x_max, y_min, y_max])
    #plt.imshow(pr, aspect="auto")
    plt.colorbar()
    plt.show()
