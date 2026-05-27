
Prerequisites and Installation
=================================

Prerequisites
----------------

To run LATTIN, you need

- `Python 3.8+ <https://www.python.org/downloads/release/python-380/>`__ 
- `Git <https://git-scm.com/>`__ 
- `Anaconda 3 <https://www.anaconda.com/>`__ 
- `Linux <https://www.linux.org/>`__ 
- `Fortran <https://fortran-lang.org/>`__ 


The main Python packages that must be installed are the following:

- `time`
- `struct`
- `datetime`
- `functools`
- `pathlib`
- `gzip`
- `shutil`
- `math`
- `fnmatch`
- `sys`
- `os`
- `imp`
- `aiohttp==3.13.5`
- `dask==2023.9.2`
- `mpi4py==3.1.4`
- `netcdf4==1.7.2`
- `numpy==1.26.4`
- `pandas==2.2.3`
- `psutil==7.0.0`
- `requests==2.34.2`
- `scipy==1.12.0`
- `setuptools==82.0.1`
- `tenacity==9.1.2`
- `tqdm==4.67.1`
- `xarray==2023.6.0` 
- `zarr==3.0.5`
- `matplotlib==3.10.3`
- `FORTRAN 90 Compiler`


Installation
------------------

**1 - First Method**

.. code-block:: bash

    conda install -c tramo-ephyslab lattin

If you install LATTIN using conda, the first time you import lattin, it automatically will compile the FORTRAN subroutines. Therefore, we strongly recommend typing the following sentences in terminal after the LATTIN installation.

.. code-block:: bash
    
    python
    import lattin

If you have problems compiling FORTRAN subroutines when importing `lattin`` the first time, you can manually compile these subroutines. Go to your Anaconda installation and run the build_lattin_so.sh script.

.. code-block:: bash
  
    cd  patho_to_anaconda_installation/.../site-packages/lattin
    sh build_lattin_so.sh

To knnow the exactly path of your Anaconda installation (`patho_to_anaconda_installation/.../site-packages/`), you can follow these instructions:

.. code-block:: bash

   python
   from distutils.sysconfig import get_python_lib
   print(get_python_lib())

**2- Second Method**
  
You must check that all the packages are installed and that there is no error message when they are imported.

- Clone the repository:

.. code-block:: bash

    git git clone https://github.com/apalarcon/LATTIN.git


- Verify you have installed all packages requiered for LATTIN (see LATTIN requirements section). If you use an Anaconda environment, please be sure you have activated the environment.

.. code-block:: bash
    cd src
    sh run install_lattin.sh

.. note::
    From now on it should be installed in the python environment and can be used like any other Python package.

Possible problems with python packages:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have a problem with the `mpi4py` library, try these steps:

- Remove the `mpi4py` library conda remove `mpi4py`
- Install the `openmpi` library `conda install conda-forge::openmpi`
- Install again the `mpi4py` library conda install `mpi4py`
- If the problem continue (the problem is frequently related with the `libmpi.so.12`  or similar), you can also try

Search the mising library on your system and link it to your Anaconda lib path.

.. code-block:: bash

    ln -s path_to_missing_library/libmpi.so.12 patho_to_anaconda_installation/lib/

or

Contact your system administrator


If you have a problem compiling the FORTRAN package in Python 3.12+, the error is caused by the shift in how Python and NumPy handle compilation in newer versions. 

The Problem:

- Python 3.12+ removed the distutils package.

- f2py (part of NumPy) now defaults to using Meson as the build backend instead of distutils.

Try these steps to solve it:

* Install `meson` and  `ninja`   

.. code-block:: bash

    conda install meson ninja

* Install `gcc_linux-64` and `gfortran_linux-64`

.. code-block:: bash
    
    conda install -c conda-forge gcc_linux-64 gfortran_linux-64