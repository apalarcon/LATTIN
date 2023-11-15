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
5 - Download the IP_mask.nc file and copy it to ```~/testing_LATTIN```

-------------------------------------------------------------------------------------------------------------
<b> You can also create your directory tree for testing LATTIN but you will need to modify input_heat.cfg and input_moisture.cfg files </b>
------------------------------------------------------------------------------------------------------------

## EXAMPLE: Atmospheric moisture tracking
### Case: Extreme precipitation event in the Iberian Peninsula from 1 to 3 March 2001
Here we provided an example of LATTIN runs for 2 March 2001
