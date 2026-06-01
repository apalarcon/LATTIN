```
============================================================================================
||            ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  ~ ~ ~ _                               ||
||            ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~                              ||
||            ++           ++++    +++++++++++ ++++++++++  [~~]  ++++++   ++++             ||
||            ++          ++  ++   +   ++    + +   ++   +  [~~]   ++ ++    ++              ||
||            ++         ++    ++      ++          ++      [~~]   ++  ++   ++              ||
||            ++        ++++++++++     ++          ++      [~~]   ++   ++  ++              ||
||            ++       ++        ++    ++          ++      [~~]   ++    ++ ++              ||
||            +++++++ ++          ++   ++          ++      [~~]  ++++   ++++++             ||
||       <<======================================================================>>        ||
============================================================================================
```
# Lagrangian Atmospheric moisTure and heaT trackINg 

LATTIN is a Python-based tool for Lagrangian atmospheric moisture and heat tracking

[![License: GPLv3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

 [![Current Version: ](https://img.shields.io/badge/Current_Version-1.1.11-blue)](https://anaconda.org/tramo-ephyslab/lattin)


This version includes:

- Functions for reading HYSPLIT (trajectory mode) data
- Update requirements to solve security issues



If you use LATTIN, please cite it as follows:

Pérez-Alarcón, A., Nieto, R., Gimeno, L. (2026). Version (1.1.10) — LATTIN: A Python-based tool for Lagrangian atmospheric moisture and heat tracking. Software Impacts, 28, 100834. https://doi.org/10.1016/j.simpa.2026.100834

Pérez-Alarcón, A.; Fernández-Alvarez, J.C.; Nieto, R.; Gimeno, L. (2024). LATTIN: A Python-based tool for Lagrangian atmospheric moisture and heat tracking. Software Impacts, 20, 100638. https://doi.org/10.1016/j.simpa.2024.100638

# What do I need to get and run LATTIN?

## To run LATTIN, you need

   * [![python](https://img.shields.io/badge/Python-3-3776AB.svg?style=flat&logo=python&logoColor=white)](https://www.python.org)
   * [![Git: ](https://img.shields.io/badge/Git--blue)](https://git-scm.com/)

and
 * [![Anaconda 3](https://img.shields.io/badge/Anaconda-3-green.svg)](https://www.anaconda.com/) (or similar to manage python packages)

or

  *  [![python](https://img.shields.io/badge/Python-3-3776AB.svg?style=flat&logo=python&logoColor=white)](https://www.python.org) and the required modules on a cluster

## The packages required to run LATTIN are:
  
```
- time
- struct
- datetime
- functools
- pathlib 
- gzip
- shutil
- math 
- fnmatch
- sys
- os
- importlib
- aiohttp==3.13.5
- dask==2023.9.2
- mpi4py==3.1.4
- netcdf4==1.7.2
- numpy==1.26.4
- pandas==2.2.3
- psutil==7.0.0
- requests==2.34.2
- scipy==1.12.0
- setuptools==82.0.1
- tenacity==9.1.2
- tqdm==4.67.1
- xarray==2023.6.0 
- zarr==3.0.5
- matplotlib==3.10.3
- FORTRAN 90 Compiler
```
# Installation

### First Method
```
conda install -c tramo-ephyslab lattin
```

If you install LATTIN using ```conda```, the first time you import lattin, it automatically will compile the FORTRAN subroutines. 
Therefore, we strongly recommend typing the following sentences in terminal after the LATTIN installation.

```
$ python
>> import lattin
```
If you have problems compiling FORTRAN subroutines when importing ```lattin``` the first time, you can manually compile these subroutines. Go to your Anaconda installation and run the build_lattin_so.sh script.

```
cd  patho_to_anaconda_installation/.../site-packages/lattin
sh build_lattin_so.sh
```
To knnow the exactly path of your Anaconda installation (```patho_to_anaconda_installation/.../site-packages/```), you can follow these instructions:

```
$python
>>  from distutils.sysconfig import get_python_lib
>> print(get_python_lib())
```


### Second Method 

1 - Clone LATTIN repository.

  ```
git clone https://github.com/apalarcon/LATTIN.git
  ```

2 - Verify you have installed all packages requiered for LATTIN (see LATTIN requirements section). If you use an Anaconda environment, please be sure you have activated the environment.

```
cd  src
```

3 - run install_lattin.sh.

### Third Method

1 - Clone LATTIN repository

  ```
git clone https://github.com/apalarcon/LATTIN.git
  ```
2 - Verify you have installed all packages requiered for LATTIN (see LATTIN requirements section). If you use an Anaconda environment, please be sure you have activated the environment

3 - Enter the src/lattin/ directory and run  buid_lattin_so.sh,  which will compile LATTIN subroutines in FORTRAN 90.

4 - Copy the lattin directory to your Anaconda instalation
```
cp -r src/lattin path_to_anaconda_installation/.../site-packages/
````

### NOTE
If you have a problem with the mpi4py library, try these steps:

* Remove the ```mpi4py``` library ```conda remove mpi4py```
* Install the ```openmpi``` library ```conda install conda-forge::openmpi ```
* Install again the ```mpi4py``` library ```conda install mpi4py ```

If the problem continue (the problem is frequently related with the ```libmpi.so.12 ``` or similar), you can also try
* Search the mising library on your system and link it to your Anaconda lib path.
  ```
  ln -s path_to_missing_library/libmpi.so.12 patho_to_anaconda_installation/lib/
  ```
or 

* Contact your system administrator

If you have a problem compiling the FORTRAN package in Python 3.12+, the error is caused by the shift in how Python and NumPy handle compilation in newer versions. 

The Problem:

- Python 3.12+ removed the distutils package.

- f2py (part of NumPy) now defaults to using Meson as the build backend instead of distutils.

Try these steps to solve it:

* Install ```meson``` and  ```ninja```   (```conda install meson ninja```)

* Install ```gcc_linux-64``` and ```gfortran_linux-64``` (```conda install -c conda-forge gcc_linux-64 gfortran_linux-64```)




# LATTIN namelist file configuration
```
===========================================================================================================
                                        - OPEN CONFIGURATION FILE -
===========================================================================================================
Parameter                        Value                      Description
verbose                       =  'True'/'False'             -> Print run logs
runID                         =  "Experimental"             -> Run name

==========================================================================================
      - PATHS -
==========================================================================================
raw_partposit_path            = "path"                      -> Directory containing input data
basedir_lara                  = 'path'                      -> Directory to read or save LARA reanalysis files
file_gz                       = 'True' / 'False'            -> Checking if model data is compresssed in gz format
output_path                   = "path"                      -> Directory to save LATTIN outputs

==========================================================================================
    - MODEL DETAILS -
==========================================================================================
model                           = 'FLEXPART'/'FLEXPART-WRF'  -> Lagrangian model
                                  'FLEXPART11'/'LARA' 
check_lara_files                = 'True'/'False'             -> To check if LARA files exists
lara_from_https                 = 'True'/'False'             -> To download LARA files if they are missing.
          
total_emited_mass               = value                      -> Total emited mass in model simulation
total_release_parcels           = value                      -> Total number of released parcels in model simulation

==========================================================================================
    - LATTIN RUN CONFIGURATION-
==========================================================================================
mode                            = "backward"                  -> Run mode  
year                            = int value  or list          -> Start year. E.g. year=2015 or year=[2015,2016,2017] 
month                           = int value or list           -> Start month. E.g. month=1 or month=[1,7,3]. 
day                             = int value or list           -> Start day. E.g. day=1 or day=[1,2,3]. 
hour                            = int value or list           -> Start hour. E.g. hour=0 or day=[0,6,12]. 
minutes                         = int value or list           -> Start minutes. E.g. minutes=0 or minutes=[0,10,20].
ndays                           = int value                   -> Number of continuos days to start the simulation.
time_step                       = int value                   -> Temporal resolution of input data [minutes]
tracking_time                   = int value                   -> Total simulation time for tracking [minutes]
calendar                        = '365d'/'366d'               -> Calendar type.
                                                                 Use calendar="365d" to discard Feb 29 in leap years

lon_left_lower_corner           = value                       -> Domain limits for regional partposit files.
lat_left_lower_corner           = value
lon_right_upper_corner          = value
lat_right_upper_corner          = value
save_full_parts_position        = True / False                -> To save raw parcel trajectories

==========================================================================================
    - MASK FILE DETAILS -
==========================================================================================
file_mask                       = 'path'                      -> Path to mask file (netcdf format)  
maskname                        = 'mask'                      -> Name of mask variable in the mask file 
maskvar_lat                     = 'lat'                       -> Latitude variable name in the mask file 
maskvar_lon                     = 'lon'                       -> Longitude variable name in the mask file
mask_value                      = value                       -> Mask value for filterirng parcels in the target region

==========================================================================================
    - OUTPUT DOMAIN RESOLUTION -
==========================================================================================
resolution                      = 0.5                         -> Output resolutiom
numPdX                          = 720                         -> Number of grid points in x-direction
numPdY                          = 360                         -> Number of grid points in y-direction
lon_lower_left                  = -180                        -> Longitude in lower left corner
lat_lower_left                  = -90                         -> Latitude in lower left corner
area_conserving_grid            = 'True'/'False'              -> To generate an area-conserving grid output

==========================================================================================
    - SPECIFIC FOR HEAT TRACKING -
==========================================================================================
tracking_heat                   = 'True'/'False'              -> Activate heat tracking

heat_tracking_method            = 'SCH19'                     -> Heat tracking method [SCH19, SCH20, JK22, CUSTOM].
                                                                 If you select one of this method [SCH19, SCH20, JK22],
                                                                 you do not to specify the next parameters.
                                                                 WARNING: The default values inly work for
                                                                 time_step=360 minutes

var_heat_track                  = 'potTemp'/'dse'             -> Variable for heat tracking
                                                                dse: Dry static energy
                                                                potTemp: Potential Temperature  

dvarheatthreshold               = value                       -> Minimun change in tracking var to be considered an uptake
                                                                 If tracking var is potential temperature,
                                                                 dvarheatthreshold is in Kelvin
                                                                 If tracking var is dry static energy,
                                                                 dvarheatthreshold is in kJ

filter_pbl_parcels              = 'True'/'False'              -> Filter parcels within the target region within the PBL
heat_custom_limits_highs        = [lower_limit, upper_limit]  -> Custom limits for filtering parcel within the target region [m]
                                                                  Set heat_custom_limits_highs = [0,0] to use PBL highs for filtering.
                                                                  Only it works if filter_pbl_parcels=True

pblcheck                        = int value                    -> checking PBL condition along the parcels trajectories
                                                                  0: no PBL check, use everything
                                                                  1: at least one location within the PBL
                                                                  2: both locations within the PBL

pbl_method                      = "maxval"                     -> PBL method for PBL check. [maxval, meanval, actualval] 
trk_rh_check                    = 'True'/'False'               -> Check relative humidity
rh_threshold                    = value                        -> Allowed relative humidity changes.
                                                                  Only needed if trk_rh_check=True

dqcheck                         = 'True'/'False'               -> Checking changes in specific humidty along the parcels trajectory.
dqthreshold                     = value                        -> Allowed changes in specific humidity.
                                                                  Only needed if dqcheck=True

heat_linear_adjustment          = 'True'/'False'               -> Apply linear adjusment to detected uptakes

save_heat_parts_position        = True / False                 -> To save processed parcels trajectories (heat)

==========================================================================================
    - SPECIFIC FOR MOISTURE TRACKING -
==========================================================================================
tracking_moisture               = 'True' /'False'             -> Activate moisture tracking
moisture_tracking_method        =  "SOD08"                    -> Misture tracking method [SOD08, SJ05, FAS19, JK22, APA22, APA25,CUSTOM]
                                                                 If you select one of this method[SOD08, SJ05, FAS19, JK22, APA22, APA25],
                                                                 you do not to specify the next parameters.
                                                                 WARNING: The default values only work for
                                                                 time_step=360 minutes


filter_dqdt_parcels             = 'True' /'False'              -> Only track precipitating parcesl within the target region
filter_pbl_dq_parcels           = 'True'/'False'               -> Filter parcels within the target region within the PBL
moist_custom_limits_highs       = [lower_limit, upper_limit]   -> Custom limits for filtering parcel within the target region [m]
                                                                  Set moist_custom_limits_highs = [0,0] to use PBL highs for filtering.
                                                                  Only it works if filter_pbl_dq_parcels=True

dqdt_threshold                  = value                        -> Change in specific humidity for considering that a precipitation
                                                                   event occurred within the target region.
                                                                   Only needed if filter_dqdt_parcels=True

precip_minrh                    = value [float]                -> Minumim relative humidity to account for precipitation  [%]             
                                                                  Set precip_minrh=0 to do not apply

dqpblcheck                      = value                        -> checking PBL condition along the parcels trajectories
                                                                  0: no PBL check, use everything
                                                                  1: at least one location within the PBL
                                                                  2: both locations within the PBL

dqpbl_method                    = 'maxval'                     -> PBL method for PBL check [maxval, meanval, actualval] 
trkdq_rh_check                  = 'True'/'False'               -> Check relative humidity 
dqrh_threshold                  = value                        -> Allowed relative humidity changes
                                                                  Only needed if trkdq_rh_check=True
mindq_gain                      = value [float]                -> Minimun change in specific humidity to be considered an uptake 
check_RH_route_precip           = 'True' /'False'              -> To evaluate RH during moisture diagnostic                                                           
precip_minrh_en_route           = value [float]                -> Minumim relative humidity to account for precipitation in route [%]  
mindq_loss                      = value [float]                -> Minimum change in specific humidity to account for moisture loss [kg/kg]


moisture_linear_adjustment      = 'True'/'False'               -> Apply linear adjusment to detected uptakes
save_moisture_parts_position    = True / False                 -> To save processed parcels trajectories (moisture)

------------------------ Bias correction approach------------------------------------
moist_bias_correction           = True / False                 -> To apply bias correction to Lagrangian precipitation estimate and moiture sources
precip_fname_prefix             = value [str]                  -> Prefix of the precipitation filename
precip_date_format              = value [str]                  -> Format of the date in the name of the precipitation file
                                                                  * "yyyymmddHM"  
                                                                  * 'yyyy-mm-dd H:M'
                                                                  * 'yyyymmddHM'
                                                                  * 'yyyy-mm-dd_H'
                                                                  * 'yyyymmdd_H'
                                                                  * "mmdd_H"
                                                                  * "mmddH"
                                                                  * "yyyymmdd_HM"

precip_path                     = value                         -> Path to the precipitation file (netcdf format).
                                                                   It should contains the accumulated preciptation over a time period equal to `time_step`
precip_lat                      = value                         -> Name of the latitude variable in the precipitation file
precip_lon                      = value                         -> Name of the longitude variable in the precipitation file
precip_var                      = value                         -> Name of the precipitation variable in the precipitation file

==========================================================================================
    - SPECIFIC FOR DRY INTRUSION ANALYSIS -
==========================================================================================
DI_analysis                     = 'True'/'False'                -> Activate dry intrusion (DI) analysis
DI_dt_thresh                    = value [float]                 -> Time period to check for DI [minutes]
DI_step                         = value [float]                 -> Time step to define start backward time to checking DI. [minutes]
DI_dp_change                    = value [float]                 -> Pressure change to account for DI  [Pa] 
DI_start_pressure_level         = value [float]                 -> Parcel should be above this level at the start point to be considered for a Dry Intrusion [Pa] 
DI_end_pressure_level           = value [float]                 -> Parcel should do downn this level to be considered for a Dry Intrusion (e.g., to PBL)  [Pa] 
DI_time_limits_checking         = list [backward minutes]       -> Minimum and maximum value of track length to check for dry intrusion. 
                                                                   Negative value for backward tracking
DI_moisture_tracking            = 'True' / 'False'              -> To account for moisture uptake of parcels with DI
DI_save_raw_parcels             = 'True' / 'False'              -> To save processed parcels trajectories (DI) 

==========================================================================================
    - SPECIFIC FOR TEMPERATURE ANOMALY TRACKING -
==========================================================================================
tracking_Tanom                    = 'True' / 'False'              -> Activating temperature anomaly tracking
Tanom_tracking_method             = value [str]                   -> Temperature anomaly tracking method: ['RP23','PR23']
path_clim_temperature             = Path [str]                    -> Path to netCDF files containing climatological temperature and other variables if needed
climT_fname_prefix                = value [str]                   -> Prefix of the name of the climatological netCDf files
climT_date_format                 = value [str]                   -> Format of the date in the name of the climatological temperature data file
                                                                    * "yyyymmddHM"  
                                                                    * 'yyyy-mm-dd H:M'
                                                                    * 'yyyymmddHM'
                                                                    * 'yyyy-mm-dd_H'
                                                                    * 'yyyymmdd_H'
                                                                    * "mmdd_H"
                                                                    * "mmddH"
                                                                    * "yyyymmdd_HM"

Tlat_var_name                     = value [str]                   -> Latitude variable name in the climatological temperature files
Tlon_var_name                     = value [str]                   -> Longitude variable name in the climatological temperature files
climTvar_name                     = value [str]                   -> Climatological air temperature variable name
dTdp_var_name                     = value [str]                   -> dTdp variable name
dTdt_var_name                     = value [str]                   -> dTdt variable name
Tplves_var_name                   = value [str]                   -> Pressure levels variable name in the climatological temperature files
Tanom_threshold                   = value list[float]             -> Temperature anomaly threshold for filtering air parcels trajectories [K]
analysis_levels                   = value list                    -> Atmopsheric levels to perform the analysis
                                                                     sfc: for surface
                                                                     pbl: for PBL
                                                                     Except for 'sfc' and 'pbl', levels should be in [hPa]
dp_sfc                            = value [float]                 -> To filter air parcels close to surface [hPa]. 
                                                                     It will retain all parcels betweenn psfc and psfc - dp_sfc

dp_upper                          = value [float]                 -> To filter air parcels at specific upper level. 
                                                                     It will retain all parcels between upper_level-dp_upper and upper_level+dp_upper
Tanom_linear_adjustment           = 'True' / 'False'              -> Apply linear adjusment to temperature changes along parcel trajectories
                                                                     It is only required for PR23 method.
save_Tanom_parts_position         = 'True' / 'False'              -> To save processed parcels trajectories (temperature anomaly tracking) 

interpolate_parcel_temperature    = 'True'/ 'False'               -> To interpotate air temperature to parcel trajectories
interpolate_sfc                   = 'True'/ 'False'               -> To interpotate surface pressure to parcel trajectories
path_to_meteodata                 = value [str]                   -> Path to meteorological data
meteodata_fname_prefix            = value [str]                   -> Prefix of the name of the meteorologocal netCDf files
meteodata_date_format             = value [str]                   -> Format of the date in the name of the meteorological data files
                                                                    * "yyyymmddHM"  
                                                                    * 'yyyy-mm-dd H:M'
                                                                    * 'yyyymmddHM'
                                                                    * 'yyyy-mm-dd_H'
                                                                    * 'yyyymmdd_H'
                                                                    * "mmdd_H"
                                                                    * "mmddH"
                                                                    * "yyyymmdd_HM"
meteolat_var_name                 = value [str]                   -> Latitude variable name in the meteorological data files
meteolon_var_name                 = value [str]                   -> Longitude variable name in the meteorological data files
psfc_var_name                     = value [str]                   -> Surface pressure variable name
Tvar_name                         = value [str]                   -> Air temperature variable name 
meteo_plves_var_name              = value [str]                   -> Pressure levels variable name in the meteorological files

===========================================================================================================
                                        - CLOSE CONFIGURATION FILE -
===========================================================================================================

```

### NOTE
- Precipitation and climatological temperature files should be one for each time step (`time_step`) of the tracking period

- We provide an example of LATTIN namelist input file in the ```testing_LATTIN folder```, but t is only focus on moisture and heat tracking. You can modifiy it or create a new one based on the description above

## Tracking methodologies

### Heat Tracking
   * SCH19: <a href="https://doi.org/10.1038/s41561-019-0431-6" target="blank"> Schumacher et al. (2019)</a>
   * SCH20: <a href="https://doi.org/10.1111/nyas.14357 " target="blank"> Schumacher et al. (2020) </a>
   * JK22:  <a href="https://doi.org/10.5194/gmd-15-1875-2022" target="blank">  Keune et al. (2022) </a>

### Moisture tracking
   * SJ05:  <a href="https://doi.org/10.1175/JHM470.1" target="blank"> Stohl and James (2004, 2005) </a>
   * SOD08: <a href="https://doi.org/10.1029/2007JD008503" target="blank"> Sodemann et al. (2008) </a>
   * FAS19: <a href="https://doi.org/10.5194/hess-23-2525-2019" target="blank"> Freme and Sodemann (2019) </a>
   * JK22:  <a href="https://doi.org/10.5194/gmd-15-1875-2022" target="blank">  Keune et al. (2022) </a>
   * APA22: <a href="https://doi.org/10.1175/JHM-D-21-0117.1" target="blank"> Pérez-Alarcón et al. (2022) </a>
   * APA25  <a href="https://doi.org/10.1016/j.atmosres.2024.107822" target="blank"> Pérez-Alarcón et al. (2025) </a>

### Temperature anomaly tracking
   * RP23:  <a href="https://doi.org/10.1038/s41561-023-01126-1" target="blank"> Röthlisberger and Papritz (2023) </a>
   * PR23:  <a href="https://doi.org/10.1029/2023GL105641" target="blank"> Papritz and Röthlisberger (2023) </a>

### Dry Intrusion analysis
   * Following <a href="https://doi.org/10.1175/JCLI-D-16-0782.1" target="blank"> Raveh-Rubin (2017) </a>

# Input data

LATTIN can read files
* from  FLEXPARTv9+ (Piso et al., 2019), FLEXPARTv11 (Bakels et al., 2024) and  FLEXPART-WRFv3.3.2 (Brioude et al., 2013) outputs in binary file format.
* from the global Lagrangian Reanalysis (LARA; <a href="https://doi.org/10.5194/essd-2025-26" target="blank"> Bakels et al., 2025 </a>)

Mask of target region for moisture and heat tracking in netCDF format

# LATTIN outputs
A netCDF file containg the spatial distribution of moisture and heat sources


# Running LATTIN

* On Linux PC
  
1 - By using run_lattin.py
```
python run_lattin.py input_file
```
* Using MPI
```
mpirun -n N_proc python run_lattin.py input_file
```

2 - You can import lattin package in your code
```
import lattin as lt

lt.lattin_main(input_file)

```

* On a HPC with Linux:

1 - Create a bash script (run_lattin.sh). This example is valid for FINESTARRAE III cluster at the Galician Supercomputing Center.
  
  ```
#!/bin/bash -l

#SBATCH --mem=512GB
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -t 7-00:00:00

module --purge
module load cesga/2020
module load miniconda3/4.9.2
conda activate envname


srun -n $SLURM_NTASKS  --mpi=pmi2 python run_lattin.py input_file
  ``` 
2 - Submit run script

```
sbatch run_lattin.sh
```

# Testing LATTIN
Please, follow the intructions at

[![Code Ocean: 10.24433/CO.4396135.v1](https://img.shields.io/badge/CodeOcean-10.24433/CO.4396135.v1-blue)](https://doi.org/10.24433/CO.4396135.v1)

[![Git: Testing_LATTIN ](https://img.shields.io/badge/GitHub-Testing_LATTIN-blue)](https://github.com/apalarcon/LATTIN/tree/main/testing_LATTIN)

# Contact and Support

* This code is not bug-free. Please report any bugs through 'Issues': https://github.com/apalarcon/LATTIN/issues

or

* Contact to Albenis Pérez Alarcón:
  
  apalarcon1991[a]gmail.com
  
  albenis.perez.alarcon[a]uvigo.es


# LICENSE

This software is published under the GPLv3 license. This means: 
1. Anyone can copy, modify and distribute this software. 
2. You have to include the license and copyright notice with each and every distribution.
3. You can use this software privately.
4. You can use this software for commercial purposes.
5. If you dare build your business solely from this code, you risk open-sourcing the whole code base.
6. If you modify it, you have to indicate changes made to the code.
7. Any modifications of this code base MUST be distributed with the same license, GPLv3.
8. This software is provided without warranty.
9. The software author or license can not be held liable for any damages inflicted by the software.

# References
* Bakels, L., Blaschek, M., Dütsch, M., Plach, A., Lechner, V., Brack, G., ... & Stohl, A. (2025). LARA: a Lagrangian Reanalysis based on ERA5 spanning from 1940 to 2023. Earth System Science Data, 17, 4569–4585, https://doi.org/10.5194/essd-17-4569-2025.
* Bakels, L., Tatsii, D., Tipka, A., Thompson, R., Dütsch, M., Blaschek, M., ... & Stohl, A. (2024). FLEXPART version 11: Improved accuracy, efficiency, and flexibility. Geoscientific Model Development, 17(21), 7595-7627. https://doi.org/10.5194/gmd-17-7595-2024.
* Brioude, J., Arnold, D., Stohl, A., Cassiani, M., Morton, D., Seibert, P., et al. 2013. The Lagrangian particle dispersion model FLEXPART-WRF version 3.1. Geosci. Model Dev., 6(6), 1889-1904. https://doi.org/10.5194/gmd-6-1889-2013.
* Keune, J., Schumacher, D.L., Miralles, D.G. 2022. A unified framework to estimate the origins of atmospheric moisture and heat using Lagrangian models. Geoscientific Model Development, 15(5), 1875-1898. Geosci. Model Dev., 15, 1875–1898. https://doi.org/10.5194/gmd-15-1875-2022.
* Papritz, L., & Röthlisberger, M. (2023). A novel temperature anomaly source diagnostic: Method and application to the 2021 heatwave in the Pacific Northwest. Geophysical Research Letters, 50(23), e2023GL105641. https://doi.org/10.1029/2023GL105641.
* Pérez-Alarcón, A., Sorí, R., Fernández-Alvarez, J.C., Nieto, R., Gimeno, L. 2022. Where does the moisture for North Atlantic tropical cyclones come from? J. Hydrometeorol., 23(3), 457–472. https://doi.org/10.1175/JHM-D-21-0117.1 
* Pisso, I., Sollum, E., Grythe, H., Kristiansen, N. I., Cassiani, M., Eckhardt, S., et al. 2019. The Lagrangian particle dispersion model FLEXPART version 10.4. Geosci. Model Dev., 12(12), 4955-4997. https://doi.org/10.5194/gmd-12-4955-2019.
* Raveh-Rubin, S. (2017). Dry intrusions: Lagrangian climatology and dynamical impact on the planetary boundary layer. Journal of Climate, 30(17), 6661-6682.  https://doi.org/10.1175/JCLI-D-16-0782.1.
* Röthlisberger, M., & Papritz, L. (2023). Quantifying the physical processes leading to atmospheric hot extremes at a global scale. Nature Geoscience, 16(3), 210-216. https://doi.org/10.1038/s41561-023-01126-1.
* Schumacher, D.L., Keune, J., Van Heerwaarden, C.C., Vilà-Guerau de Arellano, J., Teuling, A.J.,  Miralles, D.G. 2019. Amplification of mega-heatwaves through heat torrents fuelled by upwind drought. Nat. Geosci., 12, 712–717. https://doi.org/10.1038/s41561-019-0431-6.
* Schumacher, D. L., Keune, J.,  Miralles, D. G. 2020. Atmospheric heat and moisture transport to energy‐and water‐limited ecosystems. Ann. NY Acad. Sci., 1472, 123–138. https://doi.org/10.1111/nyas.14357
* Sodemann H, Schwierz C, Wernli H. 2008. Interannual variability of Greenland winter precipitation sources: Lagrangian moisture diagnostic and North Atlantic Oscillation influence. J. Geophys. Res.-Atmos., 113,D03107. https://doi.org/10.1029/2007JD008503. 
* Stohl, A., James, P. 2005. A Lagrangian analysis of the atmospheric branch of the global water cycle. Part II: Moisture transports between Earth’s ocean basins and river catchments. J. Hydrometeorol. 6(6), 961-984. https://doi.org/10.1175/JHM470.1

 