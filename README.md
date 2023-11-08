# LATTIN: LAgrangian moisTure and heaT trackINg
It is a Python-based tool for Lagrangian moisture and heat tracking

# What do I need to get and run LATTIN?

<table>
<thead>
<tr>
<th>python</th>
<th>status</th>
</tr>
</thead>
<tbody>
<tr>
<td>Python3</td>
<td> tested</td>
</tr>
</tbody>
</table>

## To run LATTIN, you need

   * python 3
   * git

and

  *  anaconda (or similar to manage python packages)

or

  *  python 3 and the required modules on a cluster

## The packages required to run LATTIN are:
  
```
- netCDF4
- numpy 
- scipy 
- mpi4py
- numpy 
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
- matplotlib
- imp
- FORTRAN 90 Compiler
```
# Installation

### First Method
  
1 - Clone LATTIN repository

  ```
https://github.com/apalarcon/LATTIN.git
  ```
2 - Verify you have installed all packages requiered for LATTIN (see LATTIN requirements section). If you use an Anaconda environment, please be sure you have activate the environment

3 - Enter the lattin/ directory and run  buid_lattin_so.sh,  which will compile LATTIN subroutines in FORTRAN 90.

4 - Copy the lattin directory to your Anaconda instalation
```
cp -r lattin path_to_anaconda_installation/lib/python3.x/site-packages/
````

### Second Method

1 - Clone LATTIN repository.

2 - Verify you have installed all packages requiered for LATTIN (see LATTIN requirements section). If you use an Anaconda environment, please be sure you have activate the environment.

3 - run install_lattin.sh.

## LATTIN namelist file configuration
```
===========================================================================================================
||                    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  ~ ~ ~ _                                      ||
||                    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~                                     ||
||                    ++           ++++    +++++++++++ ++++++++++  [~~]  ++++++   ++++                    ||
||                    ++          ++  ++   +   ++    + +   ++   +  [~~]   ++ ++    ++                     ||
||                    ++         ++    ++      ++          ++      [~~]   ++  ++   ++                     ||
||                    ++        ++++++++++     ++          ++      [~~]   ++   ++  ++                     ||
||                    ++       ++        ++    ++          ++      [~~]   ++    ++ ++                     ||
||                    +++++++ ++          ++   ++          ++      [~~]  ++++   ++++++                    ||
||               <<======================================================================>>               ||
===========================================================================================================
                                        - OPEN CONFIGURATION FILE -
===========================================================================================================
Parameter                        Value                          Description
verbose                    =  'True'/'False'                 -> Print run logs
runID                      =  "Experimental"                 -> Run name

==========================================================================================
      - PATHS -
==========================================================================================
raw_partposit_path         = "path"                          -> Directory containing input data
file_gz                    = 'True' / 'False'                -> Checking if model data is compresssed in gz format
output_path                = "path"                          -> Directory to save LATTIN outputs


==========================================================================================
    - MODEL DETAILS -
==========================================================================================
model                      = 'FLEXPART'/'FLEXPART-WRF'       -> Lagrangian model
total_emited_mass          = value                           -> Total emited mass in model simulation
total_release_parcels      = value                           -> Total number of released parcels in model simulation

==========================================================================================
    - LATTIN RUN CONFIGURATION-
==========================================================================================
mode                       = "backward"                      -> Run mode  
year                       = int value  or list              -> Start year. E.g. year=2015 or year=[2015,2016,2017] 
month                      = int value or list               -> Start month. E.g. month=1 or month=[1,7,3]. 
day                        = int value or list               -> Start day. E.g. day=1 or day=[1,2,3]. 
hour                       = int value or list               -> Start hour. E.g. hour=0 or day=[0,6,12]. 
minutes                    = int value or list               -> Start minutes. E.g. minutes=0 or minutes=[0,10,20].
ndays                      = int value                       -> Number of continuos days to start the simulation.
time_step                  = int value                       -> Temporal resolution of input data [minutes]
tracking_time              = int value                       -> Total simulation time for tracking [minutes]

lon_left_lower_corner      = value                          -> Domain limits for regional partposit files.
lat_left_lower_corner      = value
lon_right_upper_corner     = value
lat_right_upper_corner     = value

==========================================================================================
    - MASK FILE DETAILS -
==========================================================================================
file_mask                  = 'path'                         -> Path to mask file (netcdf format)  
maskname                   = 'mask'                         -> Name of mask variable in the mask file 
maskvar_lat                = 'lat'                          -> Latitude variable name in the mask file 
maskvar_lon                = 'lon'                          -> Longitude variable name in the mask file
mask_value                 = value                          -> Mask value for filterirng parcels in the target region

==========================================================================================
    - OUTPUT DOMAIN RESOLUTION -
==========================================================================================
resolution                 = 0.5                              -> Output resolutiom
numPdX                     = 720                              -> Number of grid points in x-direction
numPdY                     = 360                              -> Number of grid points in y-direction
lon_lower_left             = -180                             -> Longitude in lower left corner
lat_lower_left             = -90                              -> Latitude in lower left corner



```
