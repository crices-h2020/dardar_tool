# dardar_tool Python Pakage
This pakage provides tools for model validation with DARDAR-product. For more information about the functionalities and other features we recomend the user manual.

## Installation
**1.** Create a icare account: https://www.icare.univ-lille.fr (Not needed if download functionalities are not used.)

**2.** Install python3.  **Important:** Pakage is locked in python version 3.10.18 so that it will always work with that installation.

**3.** Install pakage with the following line:
  **pip install git+https://github.com/ADD GITHUB REPO HERE.git**

**4.** Import dardar tool with following line:
```python
import dardar_tool as dt 
```
After import run function from package to create setup files.
```python
dt.setup()
```
**Note:** download functionalities of this package need to know the key location, ergo the
location of these setup files

**5.**
You can test that everything works by printing all function handles in following way.
```python
import dardar_tool as dt
print(dt.download)
print(dt.extractor)
print(dt.extract_compare)
print(dt.check_overpass)
print(dt.read_data)
print(dt.plot_saved_data)
print(dt.plot_local_data)
print(dt.compare_user_data)
print(dt.column_cumsum)
print(dt.scale_resolution)
```

## DARDAR info
DARDAR (raDAR/liDAR) is an algorithm that combines data from NASA's CloudSat radar and the CALIPSO lidar to create comprehensive cloud products by leveraging the strengths of both instruments. By applying the optimal estimation, DARDAR compares simulated radar and lidar backscatter with actual measurements to determine the best estimate of cloud properties, even when one instrument's signal is attenuated or missing.
The algorithm (Delanoe and Hogan, 2008, 2010; Cazenave et al., 2019) retrieves cloud properties, such as particle size and ice water content, and produces detailed masks for cloud and aerosol classification with high vertical and horizontal resolution. The combination of radar's ability to penetrate deep clouds and lidar's sensitivity to thin clouds and liquid water offers a more complete picture of cloud structure and microphysics than either instrument could provide alone.
DARDAR products are DARDAR\_CLOUD and DARDAR\_MASK. DARDAR\_CLOUD focuses on retrieving ice cloud properties and microphysics from the collocated CloudSat and CALIPSO measurements. DARDAR\_MASK provides a detailed cloud and aerosol classification mask. 

DARDAR V3.00 uses the usual Brown and Francis modified mass-size relationship (https://doi.org/10.25326/450), while DARDAR V3.10 uses the Heymfield’s composite mass-size relationship (https://doi.org/10.25326/449).

dardar_tool Python Package software provides functionality to download DARDAR\_CLOUD and DARDAR\_MASK product data files from ICARE portal (https://www.icare.univ-lille.fr/dardar/). Extraction, plotting, and data resolution modification functionalities are only provided for DARDAR\_CLOUD files. Detailed descriptions of this functions with examples are provided in section \ref{sec:function_descriptions}.

### References
Delanoë, J., and R. J. Hogan, 2008: 
A variational scheme for retrieving ice cloud properties from combined radar, lidar, and infrared radiometer,  J. Geophys. Res., 113, D07204, doi:10.1029/2007JD009000.

Delanoë, J., and R. J. Hogan, 2010: 
Combined CloudSat-CALIPSO-MODIS retrievals of the properties of ice clouds., J. Geophys. Res., 115, D00H29, doi:10.1029/2009JD012346.

Cazenave, Q., Ceccaldi, M., Delanoë, J., Pelon, J., Groß, S., and Heymsfield, A.: Evolution of DARDAR-CLOUD ice cloud retrievals: new parameters and impacts on the retrieved microphysical properties, Atmos. Meas. Tech., 12, 2819–2835, https://doi.org/10.5194/amt-12-2819-2019, 2019.

## How to cite
Arttu Väisänen, & Larisa Sogacheva. (2025). A tool to download and process DARDAR data for model evaluation. Finnish Meteorological Institute. https://doi.org/10.57707/FMI-B2SHARE.593C60A403C844EE857A21775C97838D

## Acknowledgements
The development of this software has been supported by European Union’s Horizon 2020 research and innovation programme under grant agreement No. 101003826 via project CRiceS. Additionally, we wish to thank Antti Kukkurainen for his work testing the software and ICARE for providing the data and download service (https://www.icare.univ-lille.fr/dardar/).
