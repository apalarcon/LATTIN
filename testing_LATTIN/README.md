# Testing LATTIN
## Steps

1 - After install LATTIN, create a directory on your local PC: Example: ```~/testing_LATTIN```
```
mkdir ~/testing_LATTIN
```
2 - Create directory ```~/testing_LATTIN/DATA/```
```
mkdir -p ~/testing_LATTIN/DATA/
```
3 - Download raw partposit data for atmospheric moisture and heat tracking from the Zenodo repository

[![Zenodo: 10.5281/zenodo.10114851](https://img.shields.io/badge/Zenodo-10.5281/zenodo.10114851-blue)]([https://git-scm.com/](https://doi.org/10.5281/zenodo.10114851))

4 - Unzip and copy data to  ```~/testing_LATTIN/DATA/```
```
cp -r heat ~/testing_LATTIN/DATA/
cp -r moisture ~/testing_LATTIN/DATA/
```
5 - Download the IP_mask.nc file and copy it to ```~/testing_LATTIN/mask/```

6 - Download ``` input_heat.cfg ``` and ``` input_moisture.cfg ``` files and copy them to ``` ~/testing_LATTIN ```

<b>NOTE:</b> You can create your directory tree for testing LATTIN but you will need to modify ```input_heat.cfg``` and ```input_moisture.cfg``` files.


## EXAMPLE: Atmospheric moisture tracking
### Case: Extreme precipitation event in the Iberian Peninsula from 1 to 3 March 2001
Here we provided an example of LATTIN runs for 2 March 2001 using the SOD08 approach.
* Using the run_testing_latting.py file
  
```
mpirun -n 10 python run_testing_latting.py input_moisture.cfg
  ```
* Results
  You should obtain in ```~/testing_LATTIN/LATTIN_outs/TEST/moisture_SOD08/``` four netCDF files
  ```
  lattin_backward_200103020000.nc
  lattin_backward_200103020600.nc
  lattin_backward_200103021200.nc
  lattin_backward_200103021800.nc
  ```
* The accumulated moisture source pattern for this example should look like the following Figure.


## EXAMPLE: Atmospheric heat tracking
### Case: Extreme heat wave in the Iberian Peninsula from 27 June to 22 July 2015
Here we provided an example of LATTIN runs for 5 July 2015 using the SCH20 approach.
* Using the run_testing_latting.py file
  
```
mpirun -n 10 python run_testing_latting.py input_heat.cfg
  ```
* Results
  You should obtain in ```~/testing_LATTIN/LATTIN_outs/TEST/heat_SCH20/``` four netCDF files
  ```
  lattin_backward_201507050000.nc
  lattin_backward_201507050600.nc
  lattin_backward_201507051200.nc
  lattin_backward_201507051800.nc
  ```
* The accumulated heat source pattern for this example should look like the following Figure.
