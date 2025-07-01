#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# In[]:
# ########################################################################
# Functions for DARDAR validation
# ########################################################################
"""
Created on Fri Mar 24 15:40:36 2023
This script contains functions which dardar data validation tool uses
@author: vaisanea
"""
# In[]:
# ########################################################################
# Pakages
# ########################################################################

import datetime
import os
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import getpass
import pickle

from netCDF4 import Dataset
from cftime import num2date, date2num
from cryptography.fernet import Fernet
from cartopy.feature import COASTLINE


# In[]:
# ########################################################################
# General Functions
# ########################################################################


def time_formulation(data_time, year_doy):
    file_epoch = datetime.datetime.strptime(year_doy, "%Y%j")
    file_epoch = np.datetime64(file_epoch)
    data_time = data_time*1e3
    data_time = data_time.astype("timedelta64[ms]")
    data_time = file_epoch + data_time
    return data_time


def convert_time_to_int(start_time):
    time_obj = start_time.astype(datetime.datetime)
    intgtime = date2num(time_obj,"milliseconds since 01-01-2006", calendar="gregorian")
    return intgtime


def convert_int_to_time(intgtime):
    date_obj = num2date(intgtime, "milliseconds since 01-01-2006", calendar="gregorian")
    return date_obj


def filenum_check_parser(filenum):
    filenum_str = str(filenum)
    if len(filenum_str) == 1:
        filenum_str = "0000" + filenum_str
    elif len(filenum_str) == 2:
        filenum_str = "000" + filenum_str
    elif len(filenum_str) == 3:
        filenum_str = "00" + filenum_str
    elif len(filenum_str) == 4:
        filenum_str = "0" + filenum_str
    elif len(filenum_str) > 5:
        print("Unexpected behavior with filenumbers!")
    return filenum_str



# In[]:
# ########################################################################
# Parse Strings and Dates Functions
# ########################################################################


def parse_to_datetime_obj(date_str):
    length_of_date = len(date_str)
    if length_of_date == 10:
        try:
            time_str = datetime.datetime.strptime(date_str, "%Y-%m-%d")
        except ValueError:
            print("INPUT ERROR!: Could not parse date given!")
    elif (length_of_date > 11) and (length_of_date < 14):
        try:
            time_str = datetime.datetime.strptime(date_str, "%Y-%m-%d:%H")
        except ValueError:
            print("INPUT ERROR!: Could not parse date given!")
    elif (length_of_date > 14) and (length_of_date < 17):
        try:
            time_str = datetime.datetime.strptime(date_str, "%Y-%m-%d:%H-%M")
        except ValueError:
            print("INPUT ERROR!: Could not parse date given!")
    else:
        print("INPUT ERROR!: Could not parse date given!")
    return time_str


# Permited input formates for: parse_date_string
# YYYY-MM-DD
# YYYY-MM-DD:hh
# YYYY-MM-DD:hh:mm
# Parses 3 strings needed for the path from string
def parse_date_string(date_str):
    length_of_date = len(date_str)
    if length_of_date == 10:
        try:
            time_str = datetime.datetime.strptime(date_str, "%Y-%m-%d")
            year_str = time_str.strftime("%Y")
            day_str = time_str.strftime("%Y_%m_%d")
            file_str = time_str.strftime("%Y%j")
        except ValueError:
            print("INPUT ERROR!: Could not parse date given!")
    elif (length_of_date > 11) and (length_of_date < 14):
        try:
            time_str = datetime.datetime.strptime(date_str, "%Y-%m-%d:%H")
            year_str = time_str.strftime("%Y")
            day_str = time_str.strftime("%Y_%m_%d")
            file_str = time_str.strftime("%Y%j%H")
        except ValueError:
            print("INPUT ERROR!: Could not parse date given!")
    elif (length_of_date > 14) and (length_of_date < 17):
        try:
            time_str = datetime.datetime.strptime(date_str, "%Y-%m-%d:%H-%M")
            year_str = time_str.strftime("%Y")
            day_str = time_str.strftime("%Y_%m_%d")
            file_str = time_str.strftime("%Y%j%H%M")
        except ValueError:
            print("INPUT ERROR!: Could not parse date given!")
    else:
        print("INPUT ERROR!: Could not parse date given!")
    return year_str, day_str, file_str


# Permited input formates for: parse_datetime_obj
# YYYY-MM-DD
# Parses 3 strings needed for the path from datetime object
def parse_datetime_obj(datetime_obj):
    year_str = datetime_obj.strftime("%Y")
    day_str = datetime_obj.strftime("%Y_%m_%d")
    file_str = datetime_obj.strftime("%Y%j")
    return year_str, day_str, file_str


# Permited inputs
# start_date: YYYY-MM-DD
# end_date: YYYY-MM-DD
# Parses datetime list based on give date strings
def parse_time_range_list(start_date, end_date):
    start_time_str = datetime.datetime.strptime(start_date, "%Y-%m-%d")
    end_time_str = datetime.datetime.strptime(end_date, "%Y-%m-%d")
    time_diff = end_time_str - start_time_str
    date_list = [start_time_str + datetime.timedelta(days=x) for x in range(time_diff.days+1)]
    return date_list


# Parses path needed to get specific DARDAR Mask files
def parse_mask_path_str(load_str):
    year_str, day_str, file_sub_str = parse_datetime_obj(load_str)
    file_str = f"DARDAR-MASK_v1.1.4_{file_sub_str}*.hdf"
    path_str = f"/SPACEBORNE/MULTI_SENSOR/DARDAR_MASK.v1.1.4/{year_str}/{day_str}/{file_str}"
    return path_str


# Parses path needed to get specific DARDAR cloud files
def parse_cloud_path_str(load_str):
    year_str, day_str, file_sub_str = parse_datetime_obj(load_str)
    file_str = f"DARDAR-CLOUD_v2.1.1_{file_sub_str}*.hdf"
    path_str = f"/SPACEBORNE/MULTI_SENSOR/DARDAR_CLOUD.v2.1.1/{year_str}/{day_str}/{file_str}"
    return path_str


# Parse year_doy to np.datetime object
def year_doy_to_datetime(year, doy, add_time_sec):
    year_date = np.datetime64(year, "Y") + np.timedelta64(doy-1, "D")
    time = np.array(add_time_sec, dtype='timedelta64[s]')
    year_date = year_date + time
    return year_date


# In[]:
# ########################################################################
# Download Functions
# ########################################################################


# Permited inputs
# start_date: YYYY-MM-DD
# end_date: YYYY-MM-DD
def download_list_of_dardar_mask_files(start_date, end_date, save_path):
    time_range_list = parse_time_range_list(start_date, end_date)    
    parsed_paths_list = list()
    for date_str in time_range_list:
        path = parse_mask_path_str(date_str)
        parsed_paths_list.append(path)
        os.system(f"wget -nv -r -nc -np -nH --cut-dirs=3 ftp://ftp.icare.univ-lille1.fr{path} -P {save_path} --accept=hdf -e robots=off --level=0 --wait=2s")
    return parsed_paths_list


# Permited inputs
# start_date: YYYY-MM-DD
# end_date: YYYY-MM-DD
def download_list_of_dardar_cloud_files(start_date, end_date, save_path):
    time_range_list = parse_time_range_list(start_date, end_date)    
    parsed_paths_list = list()
    for date_str in time_range_list:
        path = parse_cloud_path_str(date_str)
        parsed_paths_list.append(path)
        os.system(f"wget -nv -r -nc -np -nH --cut-dirs=3 ftp://ftp.icare.univ-lille1.fr{path} -P {save_path} --accept=hdf -e robots=off --level=0 --wait=2s")
    return parsed_paths_list


# In[]:
# ########################################################################
# Filters
# ########################################################################


def lat_lon_selector(lat, lon):
    lat_int = np.ceil(lat/30 + 3)
    lat_int[lat_int == 0] = 1
    lat_int = lat_int - 1
    lat_int = np.expand_dims(lat_int, axis=1)
    lon_int = np.ceil(lon/30 + 6)
    lon_int[lon_int == 0] = 1
    lon_int = lon_int - 1
    lon_int = np.expand_dims(lon_int, axis=1)
    coords = np.concatenate((lat_int, lon_int), axis=1)
    grid_cells = np.unique(coords, axis=0)
    if grid_cells.ndim > 2:
        grid_cells = np.squeeze(grid_cells)    
    grid = np.reshape(np.arange(0, 72), (6, 12))
    for ind, cell in enumerate(grid_cells):
        point = grid[int(cell[0]), int(cell[1])]
        if ind == 0:
            grid_points = point
        else:
            grid_points = np.append(grid_points, point)
    final_array = np.zeros((72))
    final_array[grid_points] = 1
    return final_array


def location_filter(lat1, lat2, lon1, lon2, data):
    lat_r = np.linspace(lat1, lat2, 181)
    if lon1 > lon2:
        lon_pos = np.linspace(lon1, 180, 361)
        lon_neg = np.linspace(-180, lon2, 361)
        lon_r = np.append(lon_pos, lon_neg)
    else:    
        lon_r = np.linspace(lon1, lon2, 361)
    lat_arr, lon_arr = np.meshgrid(lat_r, lon_r)
    lat_arr = lat_arr.ravel()
    lon_arr = lon_arr.ravel()
    locations_points = lat_lon_selector(lat_arr, lon_arr)
    locations_points = np.append(np.array([0, 0, 0]), locations_points)
    locations_points = locations_points.astype(bool)
    filter_cols = data[:, locations_points]
    location_filter_arr = np.sum(filter_cols, axis=1)
    location_filter_arr[location_filter_arr > 0] = 1
    location_filter_arr = location_filter_arr.astype(bool)
    return location_filter_arr


# In[]:
# ########################################################################
# Open file functions
# ########################################################################


def open_mask_file(filename):
    year_doy = filename.split('_')[-2][0:7]
    with Dataset(filename) as rootgrp:
        lat = rootgrp["CLOUDSAT_Latitude"][:]
        lon = rootgrp["CLOUDSAT_Longitude"][:]
        time = rootgrp["CLOUDSAT_UTC_Time"][:]
    time = time_formulation(time, year_doy)
    return lat, lon, time


def open_cloud_file(filename):
    year_doy = filename.split('_')[-2][0:7]
    with Dataset(filename) as rootgrp:
        lat = rootgrp["latitude"][:]
        lon = rootgrp["longitude"][:]
        time = rootgrp["time"][:]
    time = time_formulation(time, year_doy)
    return lat, lon, time


def open_cloud_file_including_variable(filename, varname, version):
    if version == "V2":
        year_doy = filename.split('_')[-2][0:7]
    else:
        year_doy = filename.split('_')[-3][0:7]
    with Dataset(filename) as rootgrp:
        lat = rootgrp["latitude"][:]
        lon = rootgrp["longitude"][:]
        time = rootgrp["time"][:]
        height = rootgrp["height"][:]
        try:
            var = rootgrp[varname][:]
        except:
            print("\nProblem with variable name! Check: variable_name\n")
    time = time_formulation(time, year_doy)
    return lat, lon, time, var, height


# In[]:
# ########################################################################
# LUT Functions
# ########################################################################
"""
LUT Tools
Primitive lut is the pre alpha version including only polar areas
General lut is the alpha version including hole earth
"""
# In[]:
# ##########################
# LUT primitive
# ##########################

def select_primitive_data(lat, time):
    if np.sum(lat <= -60) > 0:
        s_p = 1
    else:
        s_p = 0
    if np.sum(lat >= 60) > 0:
        n_p = 1
    else:
        n_p = 0
    start_time = time[0]
    end_time = time[-1]
    return start_time, end_time, n_p, s_p


def create_primitive_lut(data, filename='spatiotemporal_primitive_lut'):
    with Dataset(f"{filename}.nc", "w", format="NETCDF4") as nc:
        # Global attributes
        nc.description = "Primitive time location LUT"
        # Dimensions
        nc.createDimension("rows", None)
        nc.createDimension("cols", 5)
        # Variables
        nc_out = nc.createVariable("primitive_data","f8",("rows", "cols"))
        # Attributes
        nc_out.longname = "Lookup table for primitive locations"
        # Assign values
        nc_out[:, :] = data


def write_to_primitive_lut(data, filename):
    with Dataset(f"{filename}.nc", "a") as nc:
        nc["primitive_data"][nc["primitive_data"].shape[0], :] = data


def read_from_primitive_lut(filename):
    with Dataset(f"{filename}", "r") as nc:
        data = nc.variables["primitive_data"][:, :]
    return data


def read_from_primitive_lut_time_conversion(filename):
    with Dataset(f"{filename}", "r") as nc:
        data = nc.variables["primitive_data"][:, :]
    start_time_vec = convert_int_to_time(data[:, 1])
    end_time_vec = convert_int_to_time(data[:, 2])
    return np.swapaxes(np.array([data[:, 0], start_time_vec, end_time_vec, data[:, 3], data[:, 4]]), 0, 1)


# In[]:
# ##########################
# LUT general, alpha version 
# ##########################


def select_general_lut_data(time_data, lat, lon):
    start_time = time_data[0]
    end_time = time_data[-1]
    grid_locations = lat_lon_selector(lat, lon)
    return start_time, end_time, grid_locations
    
    
def create_general_lut(data, filename='spatiotemporal_general_lut'):
    with Dataset(f"{filename}.nc", "w", format="NETCDF4") as nc:
        # Global attributes
        nc.description = "General time location LUT"
        # Dimensions
        nc.createDimension("rows", None)
        nc.createDimension("cols", 75)
        # Variables
        nc_out = nc.createVariable("general_data","f8",("rows", "cols"))
        # Attributes
        nc_out.longname = "Lookup table for whole earth locations"
        # Assign values
        nc_out[:, :] = data


def write_to_general_lut(data, filename):
    with Dataset(f"{filename}.nc", "a") as nc:
        nc["general_data"][nc["general_data"].shape[0], :] = data


def read_from_general_lut(filename):
    with Dataset(f"{filename}", "r") as nc:
        data = nc.variables["general_data"][:, :]
    return data


def read_from_general_lut_time_conversion_to_gregorian(filename):
    with Dataset(f"{filename}", "r") as nc:
        data = nc.variables["general_data"][:, :]
    start_time_vec = convert_int_to_time(data[:, 1])
    end_time_vec = convert_int_to_time(data[:, 2])
    return np.concatenate((np.expand_dims(data[:, 0], axis=1), np.expand_dims(start_time_vec, axis=1), np.expand_dims(end_time_vec, axis=1), data[:, 3:]), axis=1)


# In[]:
# ########################################################################
# Saving data
# ########################################################################


def create_results_file(time, latitude, longitude, vertical_profile, vertical_column_var, tim1, tim2, lat1, lat2, lon1, lon2, filename="cloud_data"):
    time = convert_time_to_int(time)
    with Dataset(f"{filename}.nc", "w", format="NETCDF4") as nc:
        # Describtion
        nc.description = f"Dardar cloud data from time interval {tim1} - {tim2} from latitudes {lat1} to {lat2} and longitudes {lon1} to {lon2}. The first column is time, second latitude, third longitude and the rest are either a vertical profile or just one measured variable in that spatiotemporal location."
        # Dimensions
        nc.createDimension("column", None)
        nc.createDimension("time", None)
        # Variables
        nc_tim = nc.createVariable("time", "f8", ("time"))
        nc_lat = nc.createVariable("latitude", "f8", ("time"))
        nc_lon = nc.createVariable("longitude", "f8", ("time"))
        nc_prof = nc.createVariable("vertical_profile", "f8", ("column"))
        nc_var = nc.createVariable("profile_values", "f8", ("time", "column"))
        # Attributes
        nc_tim.longname = "Time"
        nc_lat.longname = "Latitudes"
        nc_lon.longname = "Longitudes"
        nc_prof.longname = "Vertical profile of selected spatiotemporal location."
        nc_var.longname = "Cloud data of selected spatiotemporal location."
        # Assigning values
        nc_tim[:] = time
        nc_lat[:] = latitude
        nc_lon[:] = longitude
        nc_prof[:] = vertical_profile
        nc_var[:][:] = vertical_column_var


def write_to_results_file(time, latitude, longitude, vertical_column_var, filename):
    time = convert_time_to_int(time)
    with Dataset(f"{filename}.nc", "a") as nc:
        begin_ind = nc["data"].shape[0]
        end_ind = begin_ind + latitude.shape[0]
        nc["time"][begin_ind:end_ind] = time
        nc["latitude"][begin_ind:end_ind] = latitude
        nc["longitude"][begin_ind:end_ind] = longitude
        nc["profile_values"][begin_ind:end_ind][:] = vertical_column_var


def read_from_results_file_time_conversion(filename):
    if filename.split(".")[-1] == "nc":
        filename = filename.rsplit(".", 1)[0]
    with Dataset(f"{filename}.nc", "r") as nc:
        final_time = nc.variables["time"][:]
        final_lat = nc.variables["latitude"][:]
        final_lon = nc.variables["longitude"][:]
        final_profile = nc.variables["profile_values"][:][:]
        vert_profile = nc.variables["vertical_profile"][:]
    final_time = convert_int_to_time(final_time)
    return final_time, final_lat, final_lon, final_profile, vert_profile

# In[]:
# ########################################################################
# User setup
# ########################################################################

def write_pickle(filename, data):
    with open(f"{filename}.pickle", 'wb') as f:
        pickle.dump(data, f)


def setup():
    print("This setup process encrypts and saves credentials localy to your computer. It also deletes localy saved old credentials. Do you want to proceed?")
    proceed = input("yes/[no]:")
    if proceed != "yes" and proceed != "Yes" and proceed != "YES":
        print("Exiting setup without modifications")
        return
    print("Proceeding with setup")
    username = input("Enter ICARE username: ")
    passwor = getpass.getpass("Enter ICARE Pasword: ")
    keyname = Fernet.generate_key()
    cipher_n = Fernet(keyname)
    encrypted_name = cipher_n.encrypt(username.encode('utf-8'))
    keypass = Fernet.generate_key()
    cipher_p = Fernet(keypass)
    encrypted_pass = cipher_p.encrypt(passwor.encode('utf-8'))
    setup_path = input("Please enter save path for keys:")
    write_pickle(f"{setup_path}/setup_file_00", keyname)    
    write_pickle(f"{setup_path}/setup_file_01", encrypted_name)
    write_pickle(f"{setup_path}/setup_file_10", keypass)
    write_pickle(f"{setup_path}/setup_file_11", encrypted_pass)
    print("Setup completed!")


# In[]:
# ########################################################################
# Random Misc
# ########################################################################

def print_flight_path(grid_locations):
    main_str = ""
    print_str = ""
    for ind, cell_val in enumerate(grid_locations):
        if cell_val == 1:
            print_str = print_str + " X "
        else:
            print_str = print_str + " o "
        if (ind == 11 or ind == 23 or ind == 35 or ind == 47 or ind == 59):
            if ind == 11:
                main_str = print_str
            else:
                main_str = print_str + "\n" + main_str
            print_str = ""
    main_str = print_str + "\n" + main_str
    print("")
    print(main_str)
    print("")
    
    
def plot_coord_boxes(ax, grid_arr):
    grid_arr = grid_arr.astype(bool)
    up_left_pick_inds_remove = np.arange(25, 91, 13) - 13
    up_right_pick_inds_remove = np.arange(13, 91, 13) - 13
    down_left_pick_inds_remove = np.arange(12, 78, 13)
    down_right_pick_inds_remove = np.arange(0, 78, 13)
    up_arr = np.arange(13, 91)
    down_arr = np.arange(0, 78)
    up_left = np.delete(up_arr, up_left_pick_inds_remove)
    up_right = np.delete(up_arr, up_right_pick_inds_remove)
    down_left = np.delete(down_arr, down_left_pick_inds_remove)
    down_right = np.delete(down_arr, down_right_pick_inds_remove)
    coord_ind_up_left = up_left[grid_arr]
    coord_ind_up_right = up_right[grid_arr]
    coord_ind_down_left = down_left[grid_arr]
    coord_ind_down_right = down_right[grid_arr]
    lats = np.linspace(-90, 90, 7)
    lons = np.linspace(-180, 180, 13)
    ylon, xlat = np.meshgrid(lons, lats)
    ylon = ylon.flatten()
    xlat = xlat.flatten()
    xlat_up_left = xlat[coord_ind_up_left]
    xlat_up_right = xlat[coord_ind_up_right]
    xlat_down_left = xlat[coord_ind_down_left]
    xlat_down_right = xlat[coord_ind_down_right]
    ylon_up_left = ylon[coord_ind_up_left]
    ylon_up_right = ylon[coord_ind_up_right]
    ylon_down_left = ylon[coord_ind_down_left]
    ylon_down_right = ylon[coord_ind_down_right]
    coords_lat = np.concatenate((xlat_up_left[:, None], xlat_up_right[:, None], xlat_down_right[:, None], xlat_down_left[:, None]), axis=1)
    coords_lon = np.concatenate((ylon_up_left[:, None], ylon_up_right[:, None], ylon_down_right[:, None], ylon_down_left[:, None]), axis=1)
    ax.add_feature(COASTLINE, linewidth=0.5)
    for i, (lats, lons) in enumerate(zip(coords_lat, coords_lon)):
        ax.fill(lons, lats, alpha=0.4, color='red', transform=ccrs.PlateCarree())
    return ax


def plot_box_on_fig(ax, lat_1, lat_2, lon_1, lon_2):
    # Draw a box on the map by plotting the four corner points
    ax.plot([lon_1, lon_1, lon_2, lon_2, lon_1],
            [lat_1, lat_2, lat_2, lat_1, lat_1],
            color='blue', linewidth=2, transform=ccrs.PlateCarree())        


def plot_flight_box(filenum, start_date, end_date, grid_data, lat_north_limit, lat_south_limit, lon_west_limit, lon_east_limit):
    fig, ax = plt.subplots(figsize=(16, 8), subplot_kw={'projection': ccrs.PlateCarree()})
    fig.suptitle(f"File:{int(filenum)}, Time interval: {start_date}-{end_date}", fontsize=20)
    ax = plot_coord_boxes(ax, grid_data)
    extent = [-180, 180, -90, 90]
    ax.set_extent(extent)
    plot_box_on_fig(ax, lat_north_limit, lat_south_limit, lon_west_limit, lon_east_limit)
    plt.show()
