#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# In[]:
# ########################################################################
# Validation 
# ########################################################################
"""
Created on Tue May  2 14:10:06 2023

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

current_path = os.path.abspath(__file__)
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(current_path)))
working_dir = os.path.dirname(root_dir)

source_path = os.path.join(root_dir, "source")
sys.path.append(source_path)
import functions_dardar_validation as functions

download_path = os.path.join(root_dir, "main/download_function")
sys.path.append(download_path)
import download_data_based_on_input as load_functions

plotting_path = os.path.join(root_dir, "main/other_tools")
sys.path.append(plotting_path)
import plotting_functions as plt_tools

default_dardar_path = os.path.join(working_dir, "validation_data")

warnings.filterwarnings('ignore')

# In[]:
# ########################################################################
# Control paramns
# ########################################################################

#variable_name = "extinction"
#variable_name = "vis_optical_depth"

# Test periodes
# NP? 2016-05-02:04-11 -- 2016-05-15:12-00
# SP 2014-12-15 -- 2015-01-15

#start_date = "2016-05-02:04-11"
#end_date = "2016-05-15:12-00"

#start_d = "2006-07-01:04-11"
#end_d = "2006-07-04:15-00"

#raw_data_loc = ""


# Instruction for coordinate selection

# ##### lat_1 ##### #
# lon_1 ##### lon_2 #
# ##### lat_2 ##### #

# Latitude range [-90, 90]
#lat1 = 50
#lat2 = 20

# Longitude range [-180, 180]
#lon1 = 20
#lon2 = 70


#vert_prof_1 = np.arange(-150, 1500, 150)
#vert_prof_2 = np.arange(1500, 5000, 250)
#vert_prof_3 = np.arange(5000, 20000, 500)
#vert_prof_4 = np.arange(20000, 51000, 1000)

#vertical_profile = np.concatenate((vert_prof_1, vert_prof_2, vert_prof_3, vert_prof_4))


# In[]:
# ########################################################################
# Functions for file handling
# ########################################################################


def file_list_filter(all_files_list, relavant_files_list):
    final_list = list()
    for filenum in relavant_files_list:
        filenum = functions.filenum_check_parser(int(filenum))
        filenum = f"_{filenum}"
        for file_path in all_files_list:
            if filenum in file_path:
                final_list.append(file_path)
    return final_list


def check_for_missing_file(user_files_list, relevant_files_list):
    missing_files_counter = 0
    for filenum in relevant_files_list:
        filenumber_int = filenum[0]
        file_does_not_exists = True
        file_num = functions.filenum_check_parser(int(filenumber_int))
        file_num = f"_{file_num}"
        for file_path in user_files_list:
            if file_num in file_path:
                file_does_not_exists = False
                continue
        if file_does_not_exists:
            if missing_files_counter == 0:
                missing_files_list = np.expand_dims(filenum, axis=0)
            else:
                missing_files_list = np.concatenate((missing_files_list, np.expand_dims(filenum, axis=0)), axis=0)
            missing_files_counter = missing_files_counter + 1
    return missing_files_list

# In[]:
# ########################################################################
# Functions for filtering
# ########################################################################


def create_time_filter(start_date, end_date, time):
    start_date_obj = functions.parse_to_datetime_obj(start_date)
    end_date_obj = functions.parse_to_datetime_obj(end_date)
    time_filt = (start_date_obj <= time) & (time <= end_date_obj)
    return time_filt


def create_location_filter(data_lat, data_lon, lat_1, lat_2, lon_1, lon_2):
    if lat_2 <= lat_1:
        lat_filt = (lat_2 <= data_lat) & (data_lat <= lat_1)
    else:
        print("Something went wrong with latitude definition.\nLat1 >= Lat2")
    if (lon_1 < lon_2):
        lon_filt = (lon_1 <= data_lon) & (data_lon <= lon_2)
    elif (lon_1 > lon_2):
        lon_filt = (lon_1 <= data_lon) | (data_lon <= lon_2)
    else:
        print("Something went wrong with longitude definition!")
    final_location_filter = lat_filt & lon_filt
    return final_location_filter


def de_mask(filter_mask):
    if isinstance(filter_mask, np.ma.MaskedArray):
        filter_mask = filter_mask.filled(False)
    return filter_mask


# In[]:
# ########################################################################
# Functions vertical profile calculation
# ########################################################################


def calc_weighted_average_from_user_profile(user_f, dardar_h, user_half_d, user_half_u, user_p, var, top_on, bot_on):
    user_f = de_mask(user_f)
    inds = np.where(user_f)
    min_ind = np.min(inds)
    max_ind = np.max(inds)
    basic_weight_size = np.abs(dardar_h[max_ind] - dardar_h[min_ind])/(np.sum(user_f) - 1)
    main_weights = basic_weight_size*np.ones((np.sum(user_f)))
    if (min_ind > 0) and bot_on:
        bottom_ind = min_ind - 1
        bottom_weight = np.abs(np.abs(user_half_d) - np.abs(dardar_h[min_ind]))
        main_weights = np.insert(main_weights, 0, bottom_weight)
        user_f[bottom_ind] = True
    if (max_ind < (user_p.size - 2)) and top_on:
        top_ind = max_ind + 1
        top_weight = np.abs(np.abs(user_half_u) - np.abs(dardar_h[max_ind]))
        main_weights = np.append(main_weights, top_weight)
        user_f[top_ind] = True
    w_ave = np.sum(main_weights*var[user_f])/np.sum(main_weights)
    return w_ave


def find_closest(vertical_p, dardar_h, dar_var):
    closest_point_ind = np.argsort(np.abs(vertical_p - dardar_h))[0]
    w_ave = dar_var[closest_point_ind]
    return w_ave


def calculate_vertical_profile(variable, dardar_height, user_profile, exclude_mode_top_on=True, exclude_mode_bottom_on=True):
    dardar_height = np.sort(dardar_height)
    user_profile = np.sort(user_profile)
    weighted_average_new_profile = np.zeros((user_profile.shape))
    for ind, vertical_point in enumerate(user_profile):
        if ind == 0:
            # This part calculates the bottom part of profile
            if exclude_mode_bottom_on:
                user_half_up = vertical_point + np.abs(np.abs(user_profile[ind]) - np.abs(user_profile[ind+1]))/2
                user_half_dist_down = vertical_point - np.abs(np.abs(user_profile[ind]) - np.abs(user_profile[ind+1]))/2
                user_filter = (user_half_dist_down <= dardar_height) & (dardar_height <= user_half_up)
                if np.sum(user_filter) > 0:
                    # Calculates excluded bottom part if excess of dardar data points
                    weighted_ave = calc_weighted_average_from_user_profile(user_filter, dardar_height, user_half_dist_down, user_half_up, user_profile, variable, True, True)
                else:
                    # Calculates excluded bottom if sortage of dardar points
                    weighted_ave = find_closest(vertical_point, dardar_height, variable)
            else:
                user_half_up = vertical_point + np.abs(np.abs(user_profile[ind]) - np.abs(user_profile[ind+1]))/2
                user_filter = (dardar_height <= user_half_up)
                if np.sum(user_filter) > 0:
                    # Calculates included bottom part if excess of dardar data points
                    weighted_ave = calc_weighted_average_from_user_profile(user_filter, dardar_height, None, user_half_up, user_profile, variable, True, False)
                else:
                    # Calculates included bottom if sortage of dardar points
                    weighted_ave = find_closest(vertical_point, dardar_height, variable)
        elif ind == (user_profile.size - 1):
            # This part calculates the top part of profile
            if exclude_mode_top_on:
                user_half_dist_up = vertical_point + np.abs(np.abs(user_profile[ind]) - np.abs(user_profile[ind-1]))/2
                user_half_down = vertical_point - np.abs(np.abs(user_profile[ind]) - np.abs(user_profile[ind-1]))/2
                user_filter = (user_half_down <= dardar_height) & (dardar_height <= user_half_dist_up)
                if np.sum(user_filter) > 0:
                    # Calculates excluded top part if excess of dardar data points
                    weighted_ave = calc_weighted_average_from_user_profile(user_filter, dardar_height, user_half_down, user_half_dist_up, user_profile, variable, True, True)
                else:
                    # Calculates excluded top if sortage of dardar points
                    weighted_ave = find_closest(vertical_point, dardar_height, variable)
            else:
                user_half_down = vertical_point - np.abs(np.abs(user_profile[ind]) - np.abs(user_profile[ind-1]))/2
                user_filter = (user_half_down <= dardar_height)
                if np.sum(user_filter) > 0:
                    # Calculates included top part if excess of dardar data points
                    weighted_ave = calc_weighted_average_from_user_profile(user_filter, dardar_height, user_half_down, None, user_profile, variable, False, True)
                else:
                    # Calculates included top if sortage of dardar points
                    weighted_ave = find_closest(vertical_point, dardar_height, variable)
        else:
            user_half_up = vertical_point + np.abs(np.abs(user_profile[ind]) - np.abs(user_profile[ind+1]))/2
            user_half_down = vertical_point - np.abs(np.abs(user_profile[ind]) - np.abs(user_profile[ind-1]))/2
            user_filter = (user_half_down <= dardar_height) & (dardar_height <= user_half_up)
            
            if np.sum(user_filter) > 0:
                # This part calculates weighted average when there is excess of dardar data points compared to user profile
                weighted_ave = calc_weighted_average_from_user_profile(user_filter, dardar_height, user_half_down, user_half_up, user_profile, variable, True, True)
            else:
                # This part calculates weighted average when there is shortage of dardar data points compared to user profile
                vertical_point_diff = np.abs(vertical_point - dardar_height)
                closenes_index_order = np.argsort(vertical_point_diff)
                closest = closenes_index_order[0]
                if dardar_height[closest] < vertical_point:
                    for close_index in closenes_index_order:
                        if vertical_point < dardar_height[close_index]:
                            second_closest_other_direction = close_index
                            break
                else:
                    for close_index in closenes_index_order:
                        if vertical_point > dardar_height[close_index]:
                            second_closest_other_direction = close_index
                            break
                if "second_closest_other_direction" not in locals():
                    weighted_ave = variable[closest]
                else:
                    weight_closest = np.abs(dardar_height[closest] - vertical_point)
                    weight_second = np.abs(dardar_height[second_closest_other_direction] - vertical_point)
                    weighted_ave = (weight_closest*variable[closest] + weight_second*variable[second_closest_other_direction])/(weight_closest + weight_second)
        weighted_average_new_profile[ind] = weighted_ave
    return weighted_average_new_profile, user_profile


def categorical_variable_type_check(var_name):
    if var_name == "DARMASK_Simplified_Categorization":
        category_true = True
    else:
        category_true = False    
    return category_true


# In[]:
# ########################################################################
# Functions for error checking
# ########################################################################


def data_integrity_check(time,lat,lon,var):
    time_s = time.size
    lat_s = lat.size
    lon_s = lon.size
    var_s = var.shape    
    if (time_s != lat_s) | (lat_s != lon_s) | (lon_s != var_s[0]):
        #print("Filttering failed somethings wrong in data file, skipping file!")
        integrity = True
    elif (time_s == 0) | (lat_s == 0) | (lon_s == 0) | (var_s[0] == 0):
        #print("No data in file, skipping file!")
        integrity = True
    else:
        #print("All ok!")
        integrity = False
    return integrity


# In[]:
# ########################################################################
# Functions misc
# ########################################################################

def print_message(message1, message2, print_on):
    if print_on == 1:
        print(message1, end="")
    elif print_on == 2:
        print(message2, end="")


# In[]:
# ########################################################################
# Main function
# ########################################################################


def raw_validation_data_creator(start_date, end_date, lat_1, lat_2, lon_1, lon_2, variable_name, vertical_profile=None, exclude_top_on=True, exclude_bottom_on=True, verbose=1, auto_load_missing_files=True, version="V30", save_results=False, save_results_name="cloud_data", raw_data_loc=default_dardar_path):
    if version == "V2":
        file_list = sorted(glob.glob(f"{raw_data_loc}/*.hdf"))
    else:
        if version == "V31":
            file_list = sorted(glob.glob(f"{raw_data_loc}/*_V3-10.nc"))
        elif version == "V30":
            file_list = sorted(glob.glob(f"{raw_data_loc}/*_V3-00.nc"))
        else:
            print("\nProblems with chosen version. Check input: version, accepted inputs: V2, V30, V31\n")
    
    relevant_files_list = load_functions.select_files_to_download_load(start_date, end_date, lat_1, lat_2, lon_1, lon_2, version)
    file_list = file_list_filter(file_list, relevant_files_list[:, 0])

    if len(file_list) == 0:
        print(f"No files to process!\nCheck if file location path to dardar files is correct. Current path={raw_data_loc}. Only the main folder containing subfolders needs to be provided.\nCheck that time range matches your datas time range.")
        sys.exit(0)

    missing_files_count = relevant_files_list[:, 0].size - len(file_list)

    if (missing_files_count > 0) and auto_load_missing_files:
        print("Files on your data folder didn't include all files that are relevant for chosen spatiotemporal location. \nDownloading missing files.")
        missing_files_list = check_for_missing_file(file_list, relevant_files_list)
        load_functions.download_based_on_filenumber(missing_files_list, 0, raw_data_loc, version)

    elif (missing_files_count > 0) and (auto_load_missing_files == False):
        print("Files on your data folder didn't include all files that are relevant for chosen spatiotemporal location.")
        print(f"You are missing: {missing_files_count} files.")
        print("You can download with download function with current spatiotemporal location or set auto_load_missing_files=True.")
        print("Continuing preprocessing...")

    for ind, file_path in enumerate(file_list):
        print_message(f"Files done {ind+1}/{len(file_list)}\n", f"Files done {ind+1}/{len(file_list)}\n", verbose)
        
        lat, lon, tim, var, height = functions.open_cloud_file_including_variable(file_path, variable_name, version)
        var_dim = var.ndim
        
        if (ind == 0) or (ind == (len(file_list)-1)):
            time_filter = create_time_filter(start_date, end_date, tim)            
            lat = lat[time_filter]
            lon = lon[time_filter]
            tim = tim[time_filter]
            
            if var_dim == 1:
                var = var[time_filter]
            elif var_dim == 2:
                var = var[time_filter, :]
            else:
                print("Variable time filtering failed!")
            if data_integrity_check(tim, lat, lon, var):
                continue

        location_filter = create_location_filter(lat, lon, lat_1, lat_2, lon_1, lon_2)
        lat = lat[location_filter]
        lon = lon[location_filter]
        tim = tim[location_filter]
        
        if var_dim == 1:
            var = var[location_filter]
        elif var_dim == 2:
            var = var[location_filter, :]
        else:
            print("Variable location filtering failed!")
        if data_integrity_check(tim,lat,lon,var):
            continue

        if var_dim == 1:
            var_user_profile_array = np.expand_dims(var, axis=1)
            vert_profile = None
        elif var_dim == 2:
            if vertical_profile is None:
                print("\nDardar data variable contains vertical profiles. No user profile provided. Using dardar profile.\n")
                var_user_profile_array = var
                vert_profile = np.sort(height)
            else:
                vert_profile = vertical_profile
                for ii, dar_profile in enumerate(var):
                    if np.mod(ii , int(var.shape[0]*0.1)) == 0:
                        print_message("", f"File {ind+1}/{len(file_list)} Done: {ii+1}/{var.shape[0]}\n", verbose)
                    var_user_profile, _ = calculate_vertical_profile(dar_profile, height, vertical_profile, exclude_top_on, exclude_bottom_on)
                    
                    if categorical_variable_type_check(variable_name):
                        var_user_profile = np.round(var_user_profile, 0)
                    
                    if ii == 0:
                        var_user_profile_array = np.expand_dims(var_user_profile, axis=0)
                    else:
                        var_user_profile_array = np.concatenate((var_user_profile_array, np.expand_dims(var_user_profile, axis=0)), axis=0)

        if save_results:
            final_data = None
            if ind == 0 or ("final_time" not in locals()):
                final_time = tim
                functions.create_results_file(final_time, lat, lon, vert_profile, var_user_profile_array, start_date, end_date, lat_1, lat_2, lon_1, lon_2, save_results_name)
            else:
                functions.write_to_results_file(tim, lat, lon, var_user_profile_array, save_results_name)
        else:
            if ind == 0 or ("final_time" not in locals()):
                final_time = tim
                final_lat = lat
                final_lon = lon
                final_profile = var_user_profile_array
            else:
                final_time = np.concatenate((final_time, tim), axis=0)
                final_lat = np.concatenate((final_lat, lat), axis=0)
                final_lon = np.concatenate((final_lon, lon), axis=0)
                final_profile = np.concatenate((final_profile, var_user_profile_array), axis=0)


    if ("final_time" not in locals()):
        print("Files didn't contain any information for your specified spatio temporal location.")
        sys.exit(0)
    
    if save_results:
        print("Results saved!")
        return
    
    return final_time, final_lat, final_lon, final_profile, vert_profile


def extractor_plotting_tool(start_date, end_date, lat_1, lat_2, lon_1, lon_2, variable_name, vertical_profile, exclude_top_on=True, exclude_bottom_on=True, flight_separation_time=1, flight_separation_unit="m", colormap="viridis", verbose=1, auto_load_missing_files=True, version="V30", raw_data_loc=default_dardar_path):
    save_results=False
    save_results_name="" 
    profile_dardar = None
    print_message("Processing Custom Profile...\n", "Processing Custom Profile...\n", verbose)
    tim1, lat1, lon1, user_profile_values, user_profile_alts = raw_validation_data_creator(start_date, end_date, lat_1, lat_2, lon_1, lon_2, variable_name, vertical_profile, exclude_top_on, exclude_bottom_on, verbose, auto_load_missing_files, version, save_results, save_results_name, raw_data_loc)
    print_message("Processing dardar Profile...\n", "Processing dardar Profile...\n", verbose)
    tim2, lat2, lon2, dardar_profile_values, dardar_profile_alts = raw_validation_data_creator(start_date, end_date, lat_1, lat_2, lon_1, lon_2, variable_name, profile_dardar, exclude_top_on, exclude_bottom_on, verbose, auto_load_missing_files, version, save_results, save_results_name, raw_data_loc)
    plt_tools.plot_data_comparison(tim1, lat1, lon1, user_profile_alts, user_profile_values, tim2, lat2, lon2, dardar_profile_alts, dardar_profile_values, flight_separation_time=1, flight_separation_unit="m", colormap=colormap)



# In[]:
# ########################################################################
# Run Main program
# ########################################################################

if __name__ == "__main__":
    
    variable_name = "extinction"
    #variable_name = "vis_optical_depth"

    # Test periodes
    # NP? 2016-05-02:04-11 -- 2016-05-15:12-00
    # SP 2014-12-15 -- 2015-01-15

    #start_date = "2016-05-02:04-11"
    #end_date = "2016-05-15:12-00"

    start_d = "2017-08-10:05-00"
    end_d = "2017-08-11:01-00"

    raw_data_loc = ""


    # Instruction for coordinate selection

    # ##### lat_1 ##### #
    # lon_1 ##### lon_2 #
    # ##### lat_2 ##### #

    # Latitude range [-90, 90], Use ONE pole at time!
    lat1 = 50
    lat2 = 20

    # Longitude range [-180, 180]
    lon1 = 20
    lon2 = 70


    #warnings.filterwarnings('ignore')

    vert_prof_1 = np.arange(-150, 1500, 150)
    vert_prof_2 = np.arange(1500, 5000, 250)
    vert_prof_3 = np.arange(5000, 20000, 500)
    vert_prof_4 = np.arange(20000, 51000, 1000)

    vertical_profile = np.concatenate((vert_prof_1, vert_prof_2, vert_prof_3, vert_prof_4))

    
    #ddvd, vertical_profile_out = raw_validation_data_creator(start_d, end_d, lat1, lat2, lon1, lon2, variable_name, vertical_profile, save_results=True, save_results_name="asdasd", verbose=2)
    
    ddvd, vertical_profile_out = raw_validation_data_creator(start_d, end_d, lat1, lat2, lon1, lon2, variable_name, vertical_profile, save_results=False, raw_data_loc="/Users/vaisanea/Desktop/asd", verbose=2)
    
    ddvdd = ddvd[:, 3:].astype("float")
    
    x_points = vertical_profile_out
    y_scale = np.arange(np.shape(ddvd)[0])
    import matplotlib.pyplot as plt
    # Plot the 2D array
    X, Y = np.meshgrid(x_points, y_scale)

    # Plot the 2D array with custom axis scales
    plt.pcolormesh(X, Y, ddvdd, cmap='viridis')
    plt.colorbar(label='Value')
    

    #data_asd = functions.read_from_results_file_time_conversion("asdasd")
    
    
    
    