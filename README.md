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
||                                   LAgrangian moisTure and heaT trackINg                                ||
===========================================================================================================
```
LATTIN is a Python-based tool for Lagrangian moisture and heat tracking

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

# LATTIN namelist file configuration
```
===========================================================================================================
                                        - OPEN CONFIGURATION FILE -
===========================================================================================================
Parameter                        Value                      Description
verbose                    =  'True'/'False'             -> Print run logs
runID                      =  "Experimental"             -> Run name

==========================================================================================
      - PATHS -
==========================================================================================
raw_partposit_path         = "path"                      -> Directory containing input data
file_gz                    = 'True' / 'False'            -> Checking if model data is compresssed in gz format
output_path                = "path"                      -> Directory to save LATTIN outputs


==========================================================================================
    - MODEL DETAILS -
==========================================================================================
model                      = 'FLEXPART'/'FLEXPART-WRF'  -> Lagrangian model
total_emited_mass          = value                      -> Total emited mass in model simulation
total_release_parcels      = value                      -> Total number of released parcels in model simulation

==========================================================================================
    - LATTIN RUN CONFIGURATION-
==========================================================================================
mode                       = "backward"                  -> Run mode  
year                       = int value  or list          -> Start year. E.g. year=2015 or year=[2015,2016,2017] 
month                      = int value or list           -> Start month. E.g. month=1 or month=[1,7,3]. 
day                        = int value or list           -> Start day. E.g. day=1 or day=[1,2,3]. 
hour                       = int value or list           -> Start hour. E.g. hour=0 or day=[0,6,12]. 
minutes                    = int value or list           -> Start minutes. E.g. minutes=0 or minutes=[0,10,20].
ndays                      = int value                   -> Number of continuos days to start the simulation.
time_step                  = int value                   -> Temporal resolution of input data [minutes]
tracking_time              = int value                   -> Total simulation time for tracking [minutes]

lon_left_lower_corner      = value                       -> Domain limits for regional partposit files.
lat_left_lower_corner      = value
lon_right_upper_corner     = value
lat_right_upper_corner     = value

==========================================================================================
    - MASK FILE DETAILS -
==========================================================================================
file_mask                  = 'path'                      -> Path to mask file (netcdf format)  
maskname                   = 'mask'                      -> Name of mask variable in the mask file 
maskvar_lat                = 'lat'                       -> Latitude variable name in the mask file 
maskvar_lon                = 'lon'                       -> Longitude variable name in the mask file
mask_value                 = value                       -> Mask value for filterirng parcels in the target region

==========================================================================================
    - OUTPUT DOMAIN RESOLUTION -
==========================================================================================
resolution                 = 0.5                         -> Output resolutiom
numPdX                     = 720                         -> Number of grid points in x-direction
numPdY                     = 360                         -> Number of grid points in y-direction
lon_lower_left             = -180                        -> Longitude in lower left corner
lat_lower_left             = -90                         -> Latitude in lower left corner

==========================================================================================
    - SPECIFIC FOR HEAT TRACKING -
==========================================================================================
tracking_heat              = 'True'/'False'              -> Activate heat tracking

heat_tracking_method       = 'SCH19'                     -> Heat tracking method [SCH19, SCH20, JK22, CUSTOM].
                                                            If you select one of this method [SCH19, SCH20, JK22],
                                                            you do not to specify the next parameters.
                                                            WARNING: The default values inly work for
                                                            time_step=360 minutes

var_heat_track             = 'potTemp'/'dse'             -> Variable for heat tracking
                                                            dse: Dry static energy
                                                            potTemp: Potential Temperature  

dvarheatthreshold          = value                       -> Minimun change in tracking var to be considered an uptake
                                                            If tracking var is potential temperature,
                                                            dvarheatthreshold is in Kelvin
                                                            If tracking var is dry static energy,
                                                            dvarheatthreshold is in kJ

filter_pbl_parcels         = 'True'/'False'              -> Filter parcels within the target region within the PBL
heat_custom_limits_highs   = [lower_limit, upper_limit]  -> Custom limits for filtering parcel within the target region [m]
                                                            Set heat_custom_limits_highs = [0,0] to use PBL highs for filtering.
                                                            Only it works if filter_pbl_parcels=True

pblcheck                   =  int value                   -> checking PBL condition along the parcels trajectories
                                                             0: no PBL check, use everything
                                                             1: at least one location within the PBL
                                                             2: both locations within the PBL

pbl_method                 = "maxval"                     -> PBL method for PBL check. [maxval, meanval, actualval] 
trk_rh_check               = 'True'/'False'               -> Check relative humidity
rh_threshold               = value                        -> Allowed relative humidity changes.
                                                             Only needed if trk_rh_check=True

dqcheck                    = 'True'/'False'               -> Checking changes in specific humidty along the parcels trajectory.
dqthreshold                = value                        -> Allowed changes in specific humidity.
                                                             Only needed if dqcheck=True

heat_linear_adjustment     = 'True'/'False'               -> Apply linear adjusment to detected uptakes


==========================================================================================
    - SPECIFIC FOR MOISTURE TRACKING -
==========================================================================================
tracking_moisture           = 'True' /'False'             -> Activate moisture tracking
moisture_tracking_method    =  "SOD08"                    -> Misture tracking method [SOD08, SJ05, FAS19, JK22, APA22, CUSTOM]
                                                             If you select one of this method[SOD08, SJ05, FAS19, JK22],
                                                             you do not to specify the next parameters.
                                                             WARNING: The default values only work for
                                                             time_step=360 minutes


filter_dqdt_parcels         = 'True' /'True'              -> Only track precipitating parcesl within the target region
filter_pbl_dq_parcels       = 'True'/'False'              -> Filter parcels within the target region within the PBL
moist_custom_limits_highs   = [lower_limit, upper_limit]  -> Custom limits for filtering parcel within the target region [m]
                                                             Set moist_custom_limits_highs = [0,0] to use PBL highs for filtering.
                                                             Only it works if filter_pbl_dq_parcels=True

dqdt_threshold              = value                        -> Change in specific humidity for considering that a precipitation
                                                              event occurred within the target region.
                                                              Only needed if filter_dqdt_parcels=True

precip_minrh                = 80                           -> Minumim relative humidity to account for precipitation  [%]             
                                                              Set precip_minrh=0 to do not apply

dqpblcheck                  = value                        -> checking PBL condition along the parcels trajectories
                                                             0: no PBL check, use everything
                                                             1: at least one location within the PBL
                                                             2: both locations within the PBL

dqpbl_method                = 'maxval'                     -> PBL method for PBL check [maxval, meanval, actualval] 
trkdq_rh_check              = 'True'/'False'               -> Check relative humidity 
dqrh_threshold              = value                        -> Allowed relative humidity changes
                                                              Only needed if trkdq_rh_check=True
mindq_gain                  = value                        -> Minimun change in specific humidity to be considered an uptake 
moisture_linear_adjustment  = 'True'/'False'               -> Apply linear adjusment to detected uptakes

===========================================================================================================
                                        - CLOSE CONFIGURATION FILE -
===========================================================================================================
```
Please note that we provide an example of LATTIN namelist input file. You can modifiy this or create a new one based on the description above

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


# Input data

LATTIN can read files from  FLEXPARTv9+ (Piso et al., 2019) and  FLEXPART-WRFv3.3.2 (Brioude et al., 2013) outputs in binary file format.

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

1 - Create a bash script (run_lattin.sh). This example is valid for FINESTARRAE III cluster on Galician Supercomputing Center.
  
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


# Contact and Support

* This code is not bug-free. Please report any bugs through 'Issues': https://github.com/apalarcon/LATTIN/issues

or

* Contact to Albenis Pérez Alarcón:
  
  apalarcon1991[a]gmail.com
  
  albenis.perez.alarcon[a]uvigo.es


# LICENSE
Copyright 2023 Albenis Pérez-Alarcón, Patricia Coll-Hidalgo, José C. Fernández-Alvarez, Raquel Nieto and Luis Gimeno

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
* Brioude, J., Arnold, D., Stohl, A., Cassiani, M., Morton, D., Seibert, P., et al. 2013. The Lagrangian particle dispersion model FLEXPART-WRF version 3.1. Geosci. Model Dev., 6(6), 1889-1904. https://doi.org/10.5194/gmd-6-1889-2013.
* Keune, J., Schumacher, D.L., Miralles, D.G. 2022. A unified framework to estimate the origins of atmospheric moisture and heat using Lagrangian models. Geoscientific Model Development, 15(5), 1875-1898. Geosci. Model Dev., 15, 1875–1898. https://doi.org/10.5194/gmd-15-1875-2022
* Pérez-Alarcón, A., Sorí, R., Fernández-Alvarez, J.C., Nieto, R., Gimeno, L. 2022. Where does the moisture for North Atlantic tropical cyclones come from? J. Hydrometeorol., 23(3), 457–472. https://doi.org/10.1175/JHM-D-21-0117.1 
* Pisso, I., Sollum, E., Grythe, H., Kristiansen, N. I., Cassiani, M., Eckhardt, S., et al. 2019. The Lagrangian particle dispersion model FLEXPART version 10.4. Geosci. Model Dev., 12(12), 4955-4997. https://doi.org/10.5194/gmd-12-4955-2019
* Schumacher, D.L., Keune, J., Van Heerwaarden, C.C., Vilà-Guerau de Arellano, J., Teuling, A.J.,  Miralles, D.G. 2019. Amplification of mega-heatwaves through heat torrents fuelled by upwind drought. Nat. Geosci., 12, 712–717. https://doi.org/10.1038/s41561-019-0431-6.
* Schumacher, D. L., Keune, J.,  Miralles, D. G. 2020. Atmospheric heat and moisture transport to energy‐and water‐limited ecosystems. Ann. NY Acad. Sci., 1472, 123–138. https://doi.org/10.1111/nyas.14357
* Sodemann H, Schwierz C, Wernli H. Interannual variability of Greenland winter precipitation sources: Lagrangian moisture diagnostic and North Atlantic Oscillation influence. J. Geophys. Res.-Atmos. 2008; 113:D03107. https://doi.org/10.1029/2007JD008503. 
* Stohl, A., James, P. 2005. A Lagrangian analysis of the atmospheric branch of the global water cycle. Part II: Moisture transports between Earth’s ocean basins and river catchments. J. Hydrometeorol. 6(6), 961-984. https://doi.org/10.1175/JHM470.1

