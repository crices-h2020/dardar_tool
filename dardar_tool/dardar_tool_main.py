#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# In[]:
# ########################################################################
# DarDar Main
# ########################################################################
"""
Created on Mon Feb 12 17:53:32 2024

@author: Arttu Väisänen
"""

# In[]:
# ########################################################################
# Info
# ########################################################################

"""
Modules used in this whole project are:

Standard librarys:
    - os
    - glob
    - warnings
    - sys
    - shutil
    - datetime
    - getpass
    - pickle
    - cryptography    

External libraries:
    - numpy
    - netCDF4
    - cftime


"""

# In[]:
# ########################################################################
# Function
# ########################################################################

def check_module(module_name):
    try:
        __import__(module_name)
        return True
    except ImportError:
        print(f"Module '{module_name}' is not available. Please install module.")
        return False

# In[]:
# ########################################################################
# Check Dependencies and import
# ########################################################################

if check_module("sys") & check_module("os"):
    import sys
    import os
    current_path = os.path.abspath(__file__)
    sys.path.append(f"{current_path}/dardar_tool/")
else:
    sys.exit(0)

if check_module("glob") & check_module("warnings") & check_module("shutil") & check_module("datetime") & check_module("numpy") & check_module("netCDF4") & check_module("cftime") & check_module("pickle") & check_module("getpass") & check_module("cryptography") & check_module("paramiko"):
    from dardar_tool.main.download_function.download_data_based_on_input import download_based_on_input as download
    from dardar_tool.main.download_function.download_data_based_on_input import check_overpass
    from dardar_tool.main.validation_data_creator.validation_data_creator_main import raw_validation_data_creator as extractor
    from dardar_tool.main.validation_data_creator.validation_data_creator_main import extractor_plotting_tool as extract_compare
    from dardar_tool.main.other_tools.plotting_functions import plot_saved_data
    from dardar_tool.main.other_tools.plotting_functions import plot_local_data
    from dardar_tool.main.other_tools.plotting_functions import plot_model_vs_dardar_data_vertical as compare_user_data
    from dardar_tool.source.functions_dardar_validation import setup
    from dardar_tool.source.functions_dardar_validation import read_from_results_file_time_conversion as read_data
    from dardar_tool.main.other_tools.data_manipulation_tools import calculate_cumulative_values_for_column as column_cumsum
    from dardar_tool.main.other_tools.data_manipulation_tools import scale_data_resolution as scale_resolution
else:
    sys.exit(0)







