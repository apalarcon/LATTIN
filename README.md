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

## LATTIN requirements

### To run LATTIN, you need

   * python 3
   * git

and

  *  anaconda (or similar to manage python packages)

or

  *  python 3 and the required modules on a cluster

## The main packages required to run LATTIN are:
  
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

* First Method
  
1 - Clone LATTIN repository

  ```
https://github.com/apalarcon/LATTIN.git
  ```
2 - Verify you have installed all packages requiered for LATTIN (see LATTIN requirements section). If you use an Anaconda environment, please be sure you have activate the environment

3 - Enter the lattin/ directory and run  buid_lattin_so.sh,  which will compile LATTIN subroutines in FORTRAN 90.

4 - Copy the lattin directory to your Anaconda instalation

* Second Method

1 - Clone LATTIN repository.

2 - Verify you have installed all packages requiered for LATTIN (see LATTIN requirements section). If you use an Anaconda environment, please be sure you have activate the environment.

3 - run install_lattin.sh.
