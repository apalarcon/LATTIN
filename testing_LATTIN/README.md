# Testing LATTIN
## Steps

1 - You can check the reproducible LATTIN capsule in Code Ocean

[![Code Ocean: 10.24433/CO.4396135.v1](https://img.shields.io/badge/CodeOcean-10.24433/CO.4396135.v1-blue)](https://doi.org/10.24433/CO.4396135.v1)

OR

1 - After install LATTIN, create a directory on your local PC: Example: ```~/testing_LATTIN```
```
mkdir ~/testing_LATTIN
```

2 - Download raw partposit data for atmospheric moisture and heat tracking from the Zenodo repository

[![Zenodo: 10.5281/zenodo.10114851](https://img.shields.io/badge/Zenodo-10.5281/zenodo.10114851-blue)](https://doi.org/10.5281/zenodo.10114851)

3 - Unzip and copy data to  ```~/testing_LATTIN/DATA/```
```
cp -r DATA ~/testing_LATTIN/DATA/
```
5 - Download the IP_mask.nc file and copy it into ```~/testing_LATTIN/mask/```

6 - Download ``` input_heat.cfg ``` and ``` input_moisture.cfg ``` files and copy them into ``` ~/testing_LATTIN ```

<b>NOTE:</b> You can create your directory tree for testing LATTIN but you will need to modify ```input_heat.cfg``` and ```input_moisture.cfg``` files.


## EXAMPLE: Atmospheric moisture tracking
### Case: Extreme precipitation event in the Iberian Peninsula from 1 to 3 March 2001
Here we provided an example of LATTIN runs for 2 March 2001 using the SOD08 approach.
* Using the run_testing_latting.py file
  
```
mpirun -n 10 python run_testing_latting.py input_moisture.cfg
  ```
* Results
* 
  You should obtain four netCDF files in ```~/testing_LATTIN/LATTIN_outs/TEST/moisture_SOD08/```
  ```
  lattin_backward_200103020000.nc
  lattin_backward_200103020600.nc
  lattin_backward_200103021200.nc
  lattin_backward_200103021800.nc
  ```


## EXAMPLE: Atmospheric heat tracking
### Case: Extreme heat wave in the Iberian Peninsula from 27 June to 22 July 2015
Here we provided an example of LATTIN runs for 5 July 2015 using the SCH20 approach.
* Using the run_testing_latting.py file
  
```
mpirun -n 10 python run_testing_latting.py input_heat.cfg
  ```
* Results
  
  You should obtain  four netCDF files in ```~/testing_LATTIN/LATTIN_outs/TEST/heat_SCH20/```
  ```
  lattin_backward_201507050000.nc
  lattin_backward_201507050600.nc
  lattin_backward_201507051200.nc
  lattin_backward_201507051800.nc
  ```

### NOTE
The full analysis for the moisture and heat sources for these extreme events can be found in 

Pérez-Alarcón, A.; Coll-Hidalgo, P.; Fernández-Alvarez, J.C.; Nieto, R.; Gimeno, L. (2023). LATTIN: A Python-based tool for Lagrangian atmospheric moisture and heat
tracking. Environmental Modelling and Software (Under Review).


## EXAMPLE: Reading ```gz``` compressed  parposit files

1 - Download  raw input data from Zenodo repository

[![Zenodo: 10.5281/zenodo.6490365](https://img.shields.io/badge/Zenodo-10.5281/zenodo.6490365-blue)](https://doi.org/10.5281/zenodo.6490365)

2 - Copy *.gz files to ```~/testing_LATTIN/DATA/gzfiles/```
```
cp *.gz ~/testing_LATTIN/DATA/gzfiles/
```
3 - Using the run_testing_latting.py file
```
mpirun -n 10 python run_testing_latting.py input_heat_moisture_gz.cfg
  ```

* Result
  
  This example was performed for 17 October 2014 at 1800 UTC.  You should obtain a single netCDF in ```~/testing_LATTIN/LATTIN_outs/TEST/moisture_heat_gz/```
  ```
  lattin_backward_201410171800.nc
  ```

*<b>NOTE</b>: This test also shows an example of how to perform an atmospheric moisture and heat tracking in the same simulation.
