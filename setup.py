from __future__ import print_function

from distutils.dep_util import newer
import os, os.path
import setuptools
import subprocess
import sysconfig
import sys

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("lattin/VERSION", "r") as fh:
    version_ = fh.read()


# Build the functions extension
ext_suffix = sysconfig.get_config_var('EXT_SUFFIX') or '.so'
lattin_so = os.path.join('functions' + ext_suffix)


try:
	print(subprocess.check_output(
		"cd lattin/; "
		"sh buid_lattin_so.sh", shell=True))
except subprocess.CalledProcessError as e:
	print(e.output)
	print("Problems compiling the function module.  "
			"Will continue using a slower fallback...")
	sys.exit()
else:

	print()

	print("-----------------------------------------------------------------------")
	print("FORTRAN functions extension has been created as {}".format(lattin_so))
	print("-----------------------------------------------------------------------")


setuptools.setup(
    name="lattin",
    version=version_,
    developer="Albenis Pérez-Alarcón",
    CoDevelopers ="Patricia Coll-Hidalgo, José C. Fernández-Alvarez, Raquel Nieto, and Luis Gimeno",
    author_email="albenis.pérez.alarcon@uvigo.es",
    description="LATTIN is a Python-based tool for Lagrangian atmospheric moisture and heat tracking ",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/apalarcon/LATTIN/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Development Status :: Alpha",
        "Programming Language :: Python :: 3+",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    #install_requires=["numpy","mpi4py","time","struct","datetime","netCDF4","scipy","functools","pathlib","gzip","shutil","imp", "matplotlib", "sys", "os", "fnmatch","math"],
    install_requires=["numpy", "mpi4py", "netCDF4", "scipy", "matplotlib"],
    include_package_data=True,
    package_data={"":['*.so','VERSION', "constants.py",'*.f90',"_version.py","lattin_functions.py", "LAST_UPDATE"]},

)
