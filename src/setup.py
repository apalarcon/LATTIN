from __future__ import print_function

from distutils.dep_util import newer
import os, os.path
import setuptools
import subprocess
import sysconfig
import sys



def parse_requirements(filename):
    """Load requirements from a pip requirements file."""
    with open(filename, 'r') as f:
        lines = f.readlines()
    # Remove whitespace, comments, and empty lines
    requirements = [
        line.strip() for line in lines
        if line.strip() and not line.startswith('#')
    ]
    return requirements


with open("../README.md", "r") as fh:
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




# Read the requirements from the file
install_requires_ = parse_requirements('requirements.txt')



setuptools.setup(
    name="lattin",
    version=version_,
    developer="Albenis Pérez-Alarcón",
    CoDevelopers ="Raquel Nieto, and Luis Gimeno",
    author_email="albenis.pérez.alarcon@uvigo.es",
    description="LATTIN is a Python-based tool for Lagrangian atmospheric moisture and heat tracking ",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/apalarcon/LATTIN/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Development Status :: Alpha",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
    install_requires=install_requires_,
    include_package_data=True,
    package_data={"":['*.so','VERSION', "constants.py",'*.f90',"_version.py","lattin_functions.py", "LAST_UPDATE", "dry_intrusion_functions.py","lattin_rp_functions.py","lara_reanalysis.py","functions_forward.py","functions_flex11"]},

)
