#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# In[]:
# ########################################################################
# Download data Based on input
# ########################################################################
"""
Created on Tue Apr 11 10:15:43 2023
Downloads data based on location and time input. Uses LUT to download specific files.
@author: vaisanea
"""
# In[]:
# ########################################################################
# Notes
# ########################################################################
"""
Laita save pathi kondikseen kanssa silleen

"""
# In[]:
# ########################################################################
# Pakages
# ########################################################################

import os
import sys
import paramiko
import pickle
import inspect
import numpy as np

from cryptography.fernet import Fernet

current_path = os.path.abspath(__file__)
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(current_path)))
sys.path.append(f"{root_dir}/source/")

import functions_dardar_validation as functions

# In[]:
# ########################################################################
# Control params
# ########################################################################

# Accepted time string formats
# YYYY-MM-DD
# YYYY-MM-DD:hh
# YYYY-MM-DD:hh-mm
# ATTENTION: If just dates are used without hours or minuter code will interpret them just as first minute of day

#start_date = "2006-07-05"
#end_date = "2006-07-06"


# Instruction for coordinate selection

##### ##### lat_1 ##### #####
##### lon_1 ##### lon_2 #####
##### ##### lat_2 ##### #####

# Latitude range [-90, 90]
#lat_1 = 20
#lat_2 = 40

# Longitude range [-180, 180] 
# Can handel also boundary points eg. lon_1 = 160, lon_2 = -180
#lon_1 = 20
#lon_2 = 50

# Downloads cloud or mask file
# Cloud download = 0
# Mask download = 1
#cloud_0_mask_1 = 0

#save_path = ""


# In[]:
# ########################################################################
# Functions
# ########################################################################


def select_files_to_download_load(start_date, end_date, lat_1, lat_2, lon_1, lon_2, version):
    start_dt_obj = functions.parse_to_datetime_obj(start_date)
    end_dt_obj = functions.parse_to_datetime_obj(end_date)
    # Find main LUT
    current_path = os.path.abspath(__file__)
    root_dir = os.path.dirname(os.path.dirname(os.path.dirname(current_path)))
    if version == "V2":
        lut_name = "main_lut_V2.nc"
    else:
        lut_name = "main_lut_V3.nc"
    path = f"{root_dir}/source/{lut_name}"
    data = functions.read_from_general_lut_time_conversion_to_gregorian(path)
    start_time_range = (data[:, 1] >= start_dt_obj) & (data[:, 1] <= end_dt_obj)
    end_time_range = (data[:, 2] >= start_dt_obj) & (data[:, 2] <= end_dt_obj)
    time_filt = start_time_range | end_time_range
    data = data[time_filt, :]
    loc_filter = functions.location_filter(lat_1, lat_2, lon_1, lon_2, data)
    data = data[loc_filter, :]
    return data


def lut_duplicate_filter(lut):
    sorted_indices = np.argsort(lut[:, 0])
    sorted_arr = lut[sorted_indices]
    unique_values, unique_indices = np.unique(sorted_arr[:, 0], return_index=True)
    filtered_arr = sorted_arr[unique_indices]
    return filtered_arr


def read_pickle(filename):
    with open(f"{filename}.pickle", 'rb') as f:
        loaded_data = pickle.load(f)
    return loaded_data


def load_user_config(key_location):
    if key_location == "":
        key_location = "./"
    user_path = key_location.rstrip("/")
    uk = read_pickle(os.path.join(user_path, "setup_file_00"))
    un = read_pickle(os.path.join(user_path, "setup_file_01"))
    pk = read_pickle(os.path.join(user_path, "setup_file_10"))
    pn = read_pickle(os.path.join(user_path, "setup_file_11"))
    cipher_n = Fernet(uk)    
    username = cipher_n.decrypt(un)
    cipher_p = Fernet(pk)    
    password = cipher_p.decrypt(pn)
    return username, password


def download_t_v2(file_list, data_product, save_path, key_location):
    hostname = "sftp.icare.univ-lille.fr"
    port = 22
    username, password = load_user_config(key_location)
    with paramiko.Transport((hostname, port)) as transport:
        transport.connect(username=username, password=password)
        with paramiko.SFTPClient.from_transport(transport) as sftp:
            for file_info in file_list:
                filenum = int(file_info[0])
                filenum = functions.filenum_check_parser(filenum)
                year = file_info[1].strftime("%Y")
                date = file_info[1].strftime("%Y_%m_%d")
                remote_dir_path = f"/SPACEBORNE/MULTI_SENSOR/DARDAR_{data_product}/{year}/{date}/"
                directory_contents = sftp.listdir(remote_dir_path)
                file_str = f"_{filenum}.hdf"
                file_to_download = [file for file in directory_contents if file.endswith(file_str)]
                if len(file_to_download) == 1:
                    file_to_download = file_to_download[0]
                    print("Loading:", file_to_download)
                    remote_file_path = f"/SPACEBORNE/MULTI_SENSOR/DARDAR_{data_product}/{year}/{date}/{file_to_download}"
                    local_path = save_path + "/" + file_to_download
                    sftp.get(remote_file_path, local_path)


def download_t_v3(file_list, version_number, save_path, key_location):
    hostname = "sftp.icare.univ-lille.fr"
    port = 22
    username, password = load_user_config(key_location)
    with paramiko.Transport((hostname, port)) as transport:
        transport.connect(username=username, password=password)
        with paramiko.SFTPClient.from_transport(transport) as sftp:
            for file_info in file_list:
                filenum = int(file_info[0])
                filenum = functions.filenum_check_parser(filenum)
                year = file_info[1].strftime("%Y")
                date = file_info[1].strftime("%Y_%m_%d")
                remote_dir_path = f"/SPACEBORNE/CLOUDSAT/DARDAR-CLOUD.v3.{version_number}0/{year}/{date}/"
                directory_contents = sftp.listdir(remote_dir_path)
                file_str = f"_{filenum}_V3-{version_number}0.nc"
                file_to_download = [file for file in directory_contents if file.endswith(file_str)]
                if len(file_to_download) == 1:
                    file_to_download = file_to_download[0]
                    print("Loading:", file_to_download, end=" ")
                    remote_file_path = f"/SPACEBORNE/CLOUDSAT/DARDAR-CLOUD.v3.{version_number}0/{year}/{date}/{file_to_download}"
                    local_path = save_path + "/" + file_to_download
                    sftp.get(remote_file_path, local_path)
                print("Done.")


def download_based_on_filenumber(file_list, cloud_0_mask_1, save_path, version, key_location):
    if version == "V2":
        if cloud_0_mask_1 == 0:
            data_product = "CLOUD.v2.1.1"
        elif cloud_0_mask_1 == 1:
            data_product = "MASK.v1.1.4"
        else:
            print("Define data type cloud or mask!")
            data_product = None
        download_t_v2(file_list, data_product, save_path, key_location)
    else:
        if cloud_0_mask_1 == 1:
            print("Only cloud files available for V3!\nCheck function input: cloud_0_mask_1")
        else:
            version_num = version.split("V3")[-1]
            if version_num == "0" or version_num == "1":
                download_t_v3(file_list, version_num, save_path, key_location)
            else:
                print("Problems with version information!\nCheck function inputs!")


# In[]:
# ########################################################################
# Main program
# ########################################################################

def download_based_on_input(start_date, end_date, lat_1, lat_2, lon_1, lon_2, version="V30", cloud_0_mask_1=0, save_path=None, key_location=""):
    if save_path:
        save_path = os.path.normpath(save_path)
        save_path = os.path.basename(save_path)
        print(f"Using user provided directory {save_path}")
    else:
        #caller_frame = inspect.stack()[1]
        #caller_file = caller_frame.filename
        #save_path = os.path.dirname(os.path.abspath(caller_file))
        save_path = os.getcwd()
    file_list = select_files_to_download_load(start_date, end_date, lat_1, lat_2, lon_1, lon_2, version)
    file_list = lut_duplicate_filter(file_list)
    print(f"{len(file_list)} files found for this spatio temporal location.")
    print(f"Saving location: {save_path} \nProceeding with downloading...")
    download_based_on_filenumber(file_list, cloud_0_mask_1, save_path, version, key_location)


def check_overpass(start_date, end_date, lat_1, lat_2, lon_1, lon_2, version="V30", plotting_on=True):
    file_list = select_files_to_download_load(start_date, end_date, lat_1, lat_2, lon_1, lon_2, version)
    file_list = lut_duplicate_filter(file_list)
    if plotting_on:
        for file in file_list:
            functions.plot_flight_box(file[0], file[1], file[2], file[3:],lat_1, lat_2, lon_1, lon_2)
            print(f"File: {int(file[0])} containing information")
    print(f"\n{len(file_list)} files containing information of the specified location")

# In[]:
# ########################################################################
# Run Main program
# ########################################################################

if __name__ == "__main__":
    start_date = "2017-08-05"
    end_date = "2017-08-05:23-59"
    
    lat_1 = 90
    lat_2 = -90
    
    lon_1 = -180
    lon_2 = 180
    
    save_dir = "asd"
    
    #download_based_on_input(start_date, end_date, lat_1, lat_2, lon_1, lon_2, cloud_0_mask_1=0, save_path=save_dir)

