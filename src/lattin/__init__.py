import logging
import warnings
import os
import subprocess
import sysconfig

try:
	from .lattin import *
except:
	
	ext_suffix = sysconfig.get_config_var('EXT_SUFFIX') or '.so'
	lattin_so = os.path.join('functions' + ext_suffix)

	pathpkg = os.path.dirname(__file__)

	try:
		subprocess.check_output(
			"cd  " + pathpkg +";"
			"sh buid_lattin_so.sh", shell=True)
		from .lattin import *
	except subprocess.CalledProcessError as e:
		print(e.output)
		raise SystemExit("Problems compiling the function module.  "
				"Will continue using a slower fallback...")

from .fmodules import *

from .constants import *
from ._version import get_versions 
from .lattin_functions import *

	

__version__ = get_versions()['version']
__author__ = get_versions()['author']
__contact__ = get_versions()['contact']
__last_update__ = get_versions()['last_update']
del get_versions

try:
    	logging.lastResort
except AttributeError:
        logging.getLogger(__name__).addHandler(logging.StreamHandler())
