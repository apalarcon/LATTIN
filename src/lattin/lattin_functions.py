import numpy as np
import fnmatch
from netCDF4 import Dataset, date2num
from numpy import dtype
import matplotlib
import matplotlib.pylab as plt
import time
import struct
from pathlib import Path as check_PATH
import scipy.interpolate as interp
from scipy.interpolate import griddata
import time
import math
import gzip
import shutil
from mpi4py import MPI
import functools
import warnings
import sys
import os
import gc
from datetime import datetime, timedelta
from lattin.fmodules import read_binary_file as RBF
from lattin.fmodules import readbinid as RBFid
from lattin.fmodules import determined_id as D_id


from lattin.fmodules import len_file as lf
from lattin.fmodules import compute_grid_integrated_heat as compute_grid_integrated_heat
from lattin.fmodules import (
    compute_grid_integrated_moist as compute_grid_integrated_moist,
)
import lattin.constants as lc

warnings.filterwarnings("ignore")
print = functools.partial(print, flush=True)


##### HEADER FUNCTIONS #####################################
def program_name():
    """
    Returns the short name of the program.

    This function provides the abbreviated name of the Lagrangian Atmospheric moisture
    and heat tracking software, commonly referred to as LATTIN.

    Returns:
    str: The short name of the program, "LATTIN".
    """

    return "LATTIN"


def program_fullname():
    """
    Returns the full name of the program.

    This function provides the complete name of the Lagrangian Atmospheric moisture
    and heat tracking software, commonly referred to as LATTIN.

    Returns:
            str: The full name of the program, "Lagrangian Atmospheric moisTure and heaT trackINg".
    """

    return "Lagrangian Atmospheric moisTure and heaT trackINg"


def get_currentversion():
    """
    Reads the current version number from the VERSION file.

    The VERSION file is a simple text file located in the same directory
    as this Python module. It contains a single line with the version number.
    This function reads the file and returns the version number as a string.

    Returns:
            str: The current version number as a string.
    """
    pathpkg = os.path.dirname(__file__)
    version_file = pathpkg + "/VERSION"
    with open(version_file) as vfile:
        version = vfile.readlines()[0].strip()
    return version


def get_lsatupdate():
	"""
	Reads the last update date from the LAST_UPDATE file.

	The LAST_UPDATE file is a simple text file located in the same directory
	as this Python module. It contains a single line with the date of the last
	update of the LATTIN code.

	Returns:
		str: The last update date as a string.
	"""
	lupathpkg = os.path.dirname(__file__)
	version_upd = lupathpkg + "/LAST_UPDATE"
	with open(version_upd) as ufile:
		uversion = ufile.readlines()[0].strip()
	return uversion


def disclaimer():
	print(
		"\n============================================================================================================"
	)
	print("||                    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  ~ ~ ~ _                                      ||")
	print("||                    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~                                     ||")
	print("||                    ++           ++++    +++++++++++ ++++++++++  [~~]  ++++++   ++++                    ||")
	print("||                    ++          ++  ++   +   ++    + +   ++   +  [~~]   ++ ++    ++                     ||")
	print("||                    ++         ++    ++      ++          ++      [~~]   ++  ++   ++                     ||")
	print("||                    ++        ++++++++++     ++          ++      [~~]   ++   ++  ++                     ||")
	print("||                    ++       ++        ++    ++          ++      [~~]   ++    ++ ++                     ||")
	print("||                    +++++++ ++          ++   ++          ++      [~~]  ++++   ++++++                    ||")
	print("||               <<======================================================================>>               ||")
	print("||                         " + program_fullname() +" (v" +str(get_currentversion())+")                    ||")
	print("||                                        Last Update:  " +  get_lsatupdate() +"                                        ||" )
	print("||                                             Copyright 2025                                             ||")
	print("||                                                                                                        ||")
	print("||                                                                                                        ||")
	print("||                " +program_name() + " Version " +str(get_currentversion())+ " is free under the terms of GNU General Public license V3          ||")
	print("||                                  EphysLab (Environmental Physics Laboratory)                           ||")
	print("||                                           Universidade de Vigo                                         ||")
	print("||                                 contact: albenis.perez.alarcon@uvigo.es                                ||")
	print("||                            " +program_name() +" is distributed WITHOUT ANY WARRANTY! (see LICENSE)                   ||")
	print(
		"============================================================================================================\n")


def ending_credits(runtime):
	"""
	Prints ending credits and runtime information.

	Parameters
	----------
	runtime : float
		Time of execution of the program in seconds.

	"""
	print(
		"\n\n============================================================================================================"
	)
	print(
		program_name()
		+ " Version "
		+ str(get_currentversion())
		+ " has successfully finished"
	)
	print("Runtime: %.2f seconds." % np.round(runtime, 2))
	print("Bye :)")
	print(
		"============================================================================================================"
	)


def print_error_message(message):
    """
    Prints an error message and stops program execution.

    Parameters
    ----------
    message : str
            The error message to be printed.

    """
    print(
        "\n============================================================================================================="
    )
    print("ERROR: " + message)
    print(
        "============================================================================================================="
    )
    raise SystemExit("Bye :)")


##### GENERAL FUNCTIONS #############################


def check_paths(pfile, path):
    """
    Returns the value of the attribute given by `path` in the `pfile` object.
    If the attribute does not exist, returns an empty string.
    """
    try:
        fpath = getattr(pfile, path)
    except:
        fpath = ""
    return fpath


def str2boolean(arg):
    """
    Convert a string to a boolean value.

    Parameters
    ----------
    arg : str
            The string to be converted to a boolean value.

    Returns
    -------
    bool
            The boolean value corresponding to the string.

    Raises
    ------
    ArgumentTypeError
            If the string is not recognized as a boolean value.

    """
    if isinstance(arg, bool):
        return arg
    if arg.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif arg.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise arg.ArgumentTypeError("Boolean value expected.")


def get_date_format(format_string, variable):
	"""
	Convert a date format string to a Python date format string.

	Parameters
	----------
	format_string : str
		A string representing the date format.
	variable : str
		The variable name for which the date format is being requested.

	Returns
	-------
	tuple
		A tuple containing the Python date format string, a boolean indicating
		whether the format string was recognized, and a message string.

	Raises
	------
	ValueError
		If the format string is not recognized.
	"""
	formats = {
		"yyyy-mm-dd H:M": "%Y-%m-%d %H:%M",
		"yyyymmddHM": "%Y%m%d%H%M",
		"yyyy-mm-dd_H": "%Y-%m-%d_%H",
		"yyyymmdd_H": "%Y%m%d_%H",
		"mmdd_H": "%m%d_%H",
		"mmddH": "%m%d%H",
		"mmdd_HM": "%m%d_%H%M",
		"mmddHM": "%m%d%H%M",
		"yyyymmdd_HM": "%Y%m%d_%H%M",
	}
	if format_string not in formats:
		return (
			None,
			True,
			f"Unsupported format: '{format_string}' for variable '{variable}'. Available formats are 'yyyy-mm-dd H:M', 'yyyymmddHM', 'yyyy-mm-dd_H', 'yyyymmdd_H', 'yyyymmdd_HM', 'mmdd_H', 'mmddH','mmdd_HM','mmddHM'",
		)
	return formats[format_string], False, "No Message"


def check_init_parms(
    verbose,
    file_gz,
    save_full_parts_position,
    moist_bias_correction,
    tracking_Tanom,
    DI_analysis,
    tracking_moisture,
    tracking_heat,
    check_lara_files,
    lara_from_https,
    only_non_precip
):
    """
    Checks and sets default values for the initialization parameters.

    Parameters
    ----------
    verbose : bool
            Verbose mode.
    file_gz : bool
            Compress output files.
    save_full_parts_position : bool
            Save parcels' positions in full resolution.
    moist_bias_correction : bool
            Apply bias correction to moisture tracking.

    Returns
    -------
    verbose : bool
            Verbose mode.
    file_gz : bool
            Compress output files.
    save_full_parts_position : bool
            Save parcels' positions in full resolution.
    moist_bias_correction : bool
            Apply bias correction to moisture tracking.
    """
    if save_full_parts_position == "" or save_full_parts_position == None:
        save_full_parts_position = False

    if verbose == "" or verbose == None:
        verbose = True

    if file_gz == "" or file_gz == None:
        file_gz = False

    if moist_bias_correction == "" or moist_bias_correction == None:
        moist_bias_correction = False

    if tracking_Tanom == "" or tracking_Tanom == None:
        tracking_Tanom = False

    if DI_analysis == "" or DI_analysis == None:
        DI_analysis = False

    if not DI_analysis:
        DI_moisture_tracking = False

    if tracking_moisture == "" or tracking_moisture == None:
        tracking_moisture = False

    if tracking_heat == "" or tracking_heat == None:
        tracking_heat = False

    if check_lara_files == "" or check_lara_files == None:
        check_lara_files = True

    if lara_from_https == "" or lara_from_https == None:
        lara_from_https = False

    if lara_from_https:
        check_lara_files = False

    if only_non_precip=="" or only_non_precip == None:
        only_non_precip=False

    return (
        verbose,
        file_gz,
        save_full_parts_position,
        moist_bias_correction,
        tracking_Tanom,
        DI_analysis,
        tracking_moisture,
        tracking_heat,
        check_lara_files,
        lara_from_https,
        DI_moisture_tracking,
        only_non_precip
    )


def heat_tracking_parms(
    heat_tracking_method,
    trk_rh_check,
    rh_threshold,
    pblcheck,
    dqcheck,
    pbl_method,
    var_heat_track,
    dqthreshold,
    filter_pbl_parcels,
    dvarheatthreshold,
    heat_linear_adjustment,
    dtime,
    heat_custom_limits_highs,
    save_heat_parts_position,
):
    errors_found = False
    errors = ""
    if heat_tracking_method.upper() == "SCH19" and dtime == 360:
        trk_rh_check = False
        filter_pbl_parcels = True
        rh_threshold = 0
        pblcheck = 2
        dqcheck = True
        pbl_method = "maxval"
        var_heat_track = "dse"
        dqthreshold = 0.1
        filter_pbl_parcels = True
        dvarheatthreshold = 1
        heat_linear_adjustment = True
        heat_tracking_method = heat_tracking_method
        heat_custom_limits_highs = [0, 0]
    elif heat_tracking_method.upper() == "SCH20" and dtime == 360:
        trk_rh_check = False
        filter_pbl_parcels = True
        rh_threshold = 0
        pblcheck = 1
        dqcheck = False
        pbl_method = "maxval"
        var_heat_track = "dse"
        dqthreshold = 1
        filter_pbl_parcels = True
        dvarheatthreshold = 1
        heat_linear_adjustment = True
        heat_tracking_method = heat_tracking_method
        heat_custom_limits_highs = [0, 0]
    elif heat_tracking_method.upper() == "JK22" and dtime == 360:
        trk_rh_check = True
        filter_pbl_parcels = True
        rh_threshold = 10
        pblcheck = 1
        dqcheck = False
        pbl_method = "maxval"
        var_heat_track = "potTemp"
        dqthreshold = 1
        filter_pbl_parcels = True
        dvarheatthreshold = 0
        heat_linear_adjustment = True
        heat_tracking_method = heat_tracking_method
        heat_custom_limits_highs = [0, 0]

    else:
        trk_rh_check = trk_rh_check
        filter_pbl_parcels = filter_pbl_parcels
        rh_threshold = rh_threshold
        pblcheck = pblcheck
        dqcheck = dqcheck
        pbl_method = pbl_method
        var_heat_track = var_heat_track
        dqthreshold = dqthreshold
        filter_pbl_parcels = filter_pbl_parcels
        dvarheatthreshold = dvarheatthreshold
        heat_linear_adjustment = heat_linear_adjustment
        heat_tracking_method = heat_tracking_method
        heat_custom_limits_highs = heat_custom_limits_highs
        if heat_tracking_method in ("SCH19", "SCH20", "JK22"):
            print("")
            print(
                "RUN WARNING FOR HEAT TRACKING!!!!: Default values for "
                + heat_tracking_method
                + " only work for time_step = 360 minutes. Using CUSTOM tracking instead"
            )
            print("")
            time.sleep(2)
            heat_tracking_method = "CUSTOM"

        heat_message = (
            " parameter is missing in the HEAT TRACKING module of the input file \n"
        )
        heat_parms = [
            "trk_rh_check",
            "rh_threshold",
            "filter_pbl_parcels",
            "pblcheck",
            "pbl_method",
            "var_heat_track",
            "dvarheatthreshold",
            "dqcheck",
            "dqthreshold",
            "heat_linear_adjustment",
        ]

        if not heat_tracking_method in ("SCH19", "SCH20", "JK22"):
            for param in heat_parms:
                if vars()[param] == "" or vars()[param] == None:
                    errors = errors + "(**) ERROR: " + param + heat_message
                    errors_found = True

        if (
            heat_custom_limits_highs == None
            or heat_custom_limits_highs == ""
            or heat_custom_limits_highs[0] > heat_custom_limits_highs[1]
        ):
            errors = (
                errors
                + "(**) ERROR: heat_custom_limits_highs is missing or lower_high in heat_custom_limits_highs parameter in input file is higher than upper_high"
            )
            errors_found = True

    if save_heat_parts_position == "" or save_heat_parts_position == None:
        save_heat_parts_position = True

    return (
        trk_rh_check,
        rh_threshold,
        pblcheck,
        dqcheck,
        pbl_method,
        var_heat_track,
        dqthreshold,
        filter_pbl_parcels,
        dvarheatthreshold,
        heat_linear_adjustment,
        heat_tracking_method,
        heat_custom_limits_highs,
        errors,
        errors_found,
        save_heat_parts_position,
    )


def moisture_tracking_parms(
    moisture_tracking_method,
    filter_dqdt_parcels,
    dqdt_threshold,
    mindq_gain,
    mindq_loss,
    dqpbl_method,
    moisture_linear_adjustment,
    dqpblcheck,
    trkdq_rh_check,
    dqrh_threshold,
    dtime,
    filter_pbl_dq_parcels,
    moist_custom_limits_highs,
    precip_minrh,
    save_moist_parts_position,
    check_RH_route_precip,
    precip_minrh_en_route,
    precip_date_format,
    moist_bias_correction,
    mode
):

    errors_found = False
    errors = ""
    if moisture_tracking_method == "SOD08" and dtime == 360:
        filter_dqdt_parcels = True
        dqdt_threshold = -0.0002
        mindq_gain = 0.0002
        mindq_loss = -0.0002
        dqpbl_method = "meanval"
        moisture_linear_adjustment = True
        dqpblcheck = 1
        trkdq_rh_check = False
        dqrh_threshold = 100
        moisture_tracking_method = moisture_tracking_method
        filter_pbl_dq_parcels = False
        moist_custom_limits_highs = [0, 0]
        precip_minrh = 80
        check_RH_route_precip = False
        precip_minrh_en_route = 80
    elif moisture_tracking_method == "FAS19" and dtime == 360:
        filter_dqdt_parcels = True
        dqdt_threshold = -0.0001
        mindq_gain = 0.0001
        mindq_loss = -0.0001
        dqpbl_method = "meanval"
        moisture_linear_adjustment = True
        dqpblcheck = 0
        trkdq_rh_check = False
        dqrh_threshold = 100
        moisture_tracking_method = moisture_tracking_method
        filter_pbl_dq_parcels = False
        moist_custom_limits_highs = [0, 0]
        precip_minrh = 80
        check_RH_route_precip = False
        precip_minrh_en_route = 80
    elif moisture_tracking_method == "SJ05" and dtime == 360:
        filter_dqdt_parcels = False
        dqdt_threshold = 0
        if mode=="backward":
            mindq_gain = -10000000
            mindq_loss = -10000000
        elif mode=="forward":
            mindq_gain = 10000000
            mindq_loss = 10000000
        dqpbl_method = "maxval"
        moisture_linear_adjustment = False
        dqpblcheck = 0
        trkdq_rh_check = False
        dqrh_threshold = 100
        moisture_tracking_method = moisture_tracking_method
        filter_pbl_dq_parcels = False
        moist_custom_limits_highs = [0, 0]
        precip_minrh = 0
        check_RH_route_precip = False
        precip_minrh_en_route = 0
    elif moisture_tracking_method == "JK22" and dtime == 360:
        filter_dqdt_parcels = True
        dqdt_threshold = 0
        mindq_gain = 0
        mindq_loss = 0
        dqpbl_method = "maxval"
        moisture_linear_adjustment = True
        dqpblcheck = 1
        trkdq_rh_check = True
        dqrh_threshold = 20
        moisture_tracking_method = moisture_tracking_method
        filter_pbl_dq_parcels = False
        moist_custom_limits_highs = [0, 0]
        precip_minrh = 80
        check_RH_route_precip = False
        precip_minrh_en_route = 80
    elif moisture_tracking_method == "APA22" and dtime == 360:
        filter_dqdt_parcels = True
        dqdt_threshold = -0.0001
        mindq_gain = 0
        mindq_loss = 0
        dqpbl_method = "maxval"
        moisture_linear_adjustment = True
        dqpblcheck = 0
        trkdq_rh_check = False
        dqrh_threshold = 20
        moisture_tracking_method = moisture_tracking_method
        filter_pbl_dq_parcels = False
        moist_custom_limits_highs = [0, 0]
        precip_minrh = 0
        check_RH_route_precip = False
        precip_minrh_en_route = 0
    elif moisture_tracking_method == "APA25" and dtime == 360:
        filter_dqdt_parcels = True
        dqdt_threshold = 0
        mindq_gain = 0
        mindq_loss = 0
        dqpbl_method = "maxval"
        moisture_linear_adjustment = True
        dqpblcheck = 0
        trkdq_rh_check = True
        dqrh_threshold = 20
        moisture_tracking_method = moisture_tracking_method
        filter_pbl_dq_parcels = False
        moist_custom_limits_highs = [0, 0]
        precip_minrh = 65
        check_RH_route_precip = True
        precip_minrh_en_route = 65
    else:
        filter_dqdt_parcels = filter_dqdt_parcels
        dqdt_threshold = dqdt_threshold
        mindq_gain = mindq_gain
        mindq_loss = mindq_loss
        dqpbl_method = dqpbl_method
        moisture_linear_adjustment = moisture_linear_adjustment
        dqpblcheck = dqpblcheck
        trkdq_rh_check = trkdq_rh_check
        dqrh_threshold = dqrh_threshold
        moisture_tracking_method = moisture_tracking_method
        filter_pbl_dq_parcels = filter_pbl_dq_parcels
        moist_custom_limits_highs = moist_custom_limits_highs

        check_RH_route_precip = check_RH_route_precip
        precip_minrh_en_route = precip_minrh_en_route

        if moisture_tracking_method in ("SOD08", "FAS19", "SJ05", "JK22", "APA22"):
            print("")
            print(
                "RUN WARNING FOR MOISTURE TRACKING!!!!: Default values for "
                + moisture_tracking_method
                + " only work for time_step = 360 minutes. Using CUSTOM tracking instead"
            )
            print("")
            time.sleep(2)
            moisture_tracking_method = "CUSTOM"

        moisture_message = (
            " parameter is missing in the MOISTURE TRACKING module of the input file\n"
        )
        moist_parms = [
            "filter_dqdt_parcels",
            "dqdt_threshold",
            "dqpblcheck",
            "dqpbl_method",
            "trkdq_rh_check",
            "dqrh_threshold",
            "moisture_linear_adjustment",
            "filter_pbl_dq_parcels",
            "precip_minrh",
            "check_RH_route_precip",
            "precip_minrh_en_route",
            "mindq_loss",
            "mindq_gain",
        ]
        if not moisture_tracking_method in ("SOD08", "FAS19", "SJ05", "JK22", "APA22"):
            for param in moist_parms:
                if vars()[param] == "" or vars()[param] == None:
                    errors = errors + "(**) ERROR: " + param + moisture_message
                    errors_found = True
        if (
            moist_custom_limits_highs == None
            or moist_custom_limits_highs == ""
            or moist_custom_limits_highs[0] > moist_custom_limits_highs[1]
        ):
            errors = (
                errors
                + "(**) ERROR: moist_custom_limits_highs is missing or lower_high in moist_custom_limits_highs parameter in input file is higher than upper_high"
            )
            errors_found = True

    if save_moist_parts_position == "" or save_moist_parts_position == None:
        save_moist_parts_position = True

    if moist_bias_correction:
        if not precip_date_format == "" or not precip_date_format == None:
            precip_date_format, etemp, emessage = get_date_format(
                precip_date_format, "precip_date_format"
            )
            if etemp:
                errors = errors + "(**) ERROR: " + emessage
                errors_found = True

    return (
        filter_dqdt_parcels,
        dqdt_threshold,
        mindq_gain,
        mindq_loss,
        dqpbl_method,
        moisture_linear_adjustment,
        dqpblcheck,
        trkdq_rh_check,
        dqrh_threshold,
        moisture_tracking_method,
        filter_pbl_dq_parcels,
        moist_custom_limits_highs,
        precip_minrh,
        errors,
        errors_found,
        save_moist_parts_position,
        check_RH_route_precip,
        precip_minrh_en_route,
        precip_date_format,
    )


def checking_input_parameters(
    size,
    totalfiles,
    totaltime,
    dtime,
    file_mask,
    model,
    mode,
    numPdx,
    numPdy,
    resolution,
    latgrid,
    longrid,
    ndays,
    cenlon,
    DI_analysis,
    tracking_moisture,
    tracking_Tanom,
    tracking_heat,
    moisture_tracking_method,
):
    errors = ""
    errors_found = False
    if (totalfiles) % size != 0:
        errors = (
            "(**) ERROR: The number of processors must exactly divide the number of partposit files to process (tracking_time/time_step).\n     From your input file, the recommended number of processors is "
            + str(int(totaltime / dtime))
            + "\n"
        )
        errors_found = True

    try:
        nc_ = Dataset(file_mask)
    except:
        errors = errors + "(**) ERROR: No such file or directory: " + file_mask + "\n"
        errors_found = True

    if not model in ["FLEXPART", "FLEXPART-WRF", "LARA", "FLEXPART11", "HPT",'HYSPLIT']:
        errors = errors + "(**) ERROR: Model must be FLEXPART, FLEXPART11, FLEXPART-WRF, LARA, HYSPLIT or HPT\n"
        errors_found = True

    if not mode in ("backward","forward"):
        errors = errors + "(**) ERROR: mode must be 'backward' or 'forward' in time\n"
        errors_found = True
    #elif mode in ("forward") and tracking_moisture==True and not moisture_tracking_method in ("SJ05","APA22","FAS19","SOD08","JK22",):
    #    errors = errors + "(**) ERROR: mode must be 'forward' in time only for moisture_tracking=True and moisture_tracking_method='SJ05' \n"
    #    errors_found = True
    elif mode == "forward" and (DI_analysis or tracking_Tanom or tracking_heat):
        errors = errors + "(**) ERROR: mode must be 'forward' in time only for moisture_tracking=True and moisture_tracking_method='SJ05' \n"
        errors_found = True




    if numPdx <= 0 or numPdy <= 0:
        errors = errors + "(**) ERROR: numPdX  and  numPdY  must be greater than zero\n"
        errors_found = True
    if resolution <= 0:
        errors = errors + "(**) ERROR:  resolution must be greater than zero \n"
        errors_found = True

    if latgrid.min() < -90 or latgrid.max() > 90:
        errors = errors + "(**) ERROR:  Latitudes for output grid are incorrect \n"
        errors_found = True

    if model == "LARA" and not dtime % 60 == 0:
        errors = (
            errors
            + "(**) ERROR:  LARA reanalysis is available for a time step multiple of 60  \n"
        )
        errors_found = True

    if cenlon == "0":
        if longrid.min() < -180 or longrid.max() > 180:
            errors = errors + "(**) ERROR:  Longitudes for output grid are incorrect \n"
            errors_found = True
    elif cenlon == "180":
        if longrid.min() < -0 or longrid.max() > 360:
            errors = errors + "(**) ERROR:  Longitudes for output grid are incorrect \n"
            errors_found = True

    if ndays < 1:
        errors = (
            errors
            + "(**) ERROR:  ndays = "
            + str(int(ndays))
            + " is incorrect. Parameter ndays must be higher than 0  \n"
        )
        errors_found = True

    return errors, errors_found


def checking_raw_partposti_files(partpositfiles):
    errors = ""
    find_error = False
    for i in partpositfiles:
        my_file = check_PATH(i)
        if not my_file.is_file():
            errors = errors + "     (**) ERROR: No such file or directory: " + i + "\n"
            find_error = True

    return errors, find_error


def checking_precip_files(
    precip_path,
    precip_fname_prefix,
    list_year,
    list_month,
    list_day,
    list_hour,
    list_min,
    precip_date_format,
):
    errors = ""
    find_error = False

    for year, month, day, hour, minn in zip(
        list_year, list_month, list_day, list_hour, list_min
    ):
        itime = datetime.strptime(
            year + "-" + month + "-" + day + " " + hour + ":" + minn + ":00",
            "%Y-%m-%d %H:%M:%S",
        )

        precipfile = f"{precip_path}/{precip_fname_prefix}_{itime.strftime(precip_date_format)}.nc"

        try:
            nc = Dataset(precipfile)
        except:
            errors = (
                errors
                + "     (**) ERROR: No such file or directory: "
                + precipfile
                + "\n"
            )
            find_error = True

    return errors, find_error


def get_lara_config():

    base_url = "https://data.eodc.eu/collections/LARA"
    periods = """
    1940-1947
    1947-1954
    1954-1961
    1961-1968
    1968-1975
    1975-1982
    1982-1989
    1989-1996
    1996-2003
    2003-2010
    2010-2017
    2017-2024
    """

    nparcels = 5914326
    atmos_mass = 5.09256513e18

    variables = ["lat", "lon", "T", "z", "sh", "rho", "prs", "hmix"]
    gridded_variables = ["hmix"]

    return base_url, periods, nparcels, atmos_mass, variables, gridded_variables


def printting_run_information(
    verbose,
    raw_partposit_path,
    output_path,
    list_year,
    list_month,
    list_day,
    list_hour,
    list_min,
    ndays,
    model,
    mode,
    total_emited_mass,
    total_release_parcels,
    totaltime,
    dtime,
    var_heat_track,
    file_mask,
    heat_tracking_method,
    moisture_tracking_method,
    tracking_heat,
    tracking_moisture,
    density,
    dvarheatthreshold,
    pbl_method,
    pblcheck,
    dqcheck,
    dqthreshold,
    trk_rh_check,
    rh_threshold,
    filter_pbl_parcels,
    filter_dqdt_parcels,
    dqdt_threshold,
    dqpblcheck,
    mindq_gain,
    dqpbl_method,
    dqtrk_rh_check,
    dqrh_threshold,
    heat_linear_adjustment,
    moisture_linear_adjustment,
    heat_custom_limits_highs,
    filter_pbl_dq_parcels,
    moist_custom_limits_highs,
    precip_minrh,
    calendar,
    runID,
    save_heat_parts_position,
    save_moist_parts_position,
    save_full_parts_position,
    cenlon,
    moist_bias_correction,
    check_RH_route_precip,
    precip_minrh_en_route,
    mindq_loss,
    tracking_Tanom,
    Tanom_tracking_method,
    save_Tanom_parts_position,
    analysis_levels,
    dp_sfc,
    dp_upper,
    Tanom_linear_adjustment,
    Tanom_threshold,
    DI_dt_thresh,
    DI_step,
    DI_dp_change,
    DI_start_pressure_level,
    DI_end_pressure_level,
    DI_time_limits_checking,
    DI_moisture_tracking,
    DI_analysis,
    DI_moist_bias_correction,
    DI_save_raw_parcels,
    area_conserving_grid,
    basedir_lara,
    lara_from_https,
    check_lara_files,
    interpolate_parcel_temperature,
    path_clim_temperature,
    only_non_precip
):

    #print("-------------------------------------------------------------------------------------------------------------\n")
	if verbose:
		print("RUNNING PARAMETERS")
		print("-------------------------------------------------------------------------------------------------------------")
		print("->  | RunID                                      :", runID)
		print("-------------------------------------------------------------------------------------------------------------")

		if model=="LARA":
			print("   -> Model                                      :", f"FLEXPART ({model} Reanalysis)")
		else:
			print("   -> Model                                      :", model)
		print("   -> Mask file                                  :", file_mask)

		if model=="LARA":
			print("   -> Read LARA from zarr https                  :", lara_from_https)
			if not lara_from_https:
				print("   -> LARA basedir                               :", basedir_lara)
				print("   -> Check LARA local files before starting     :", check_lara_files)
			else:
				url_zarr, _, _, _,_, _ = get_lara_config()
				print("   -> LARA zarr                                  :", url_zarr)
		else:
			print("   -> Raw partposit data                         :", raw_partposit_path)
		print("   -> Output directory                           :", output_path)
		print("   -> Tracking mode                              :", mode, "in time")
		print("   -> Time step                                  :", dtime, "minutes")
		print("   -> Tracking time                              :", totaltime, "minutes", "(" +str(totaltime/1440) + " days)")
		print("   -> Total emitted mass                         :", total_emited_mass, "kg")
		print("   -> Total released parcels                     :", total_release_parcels,"parcels")
		print("   -> Heat tracking                              :", tracking_heat)
		print("   -> Moisture tracking                          :", tracking_moisture)
		print("   -> Temperature anomaly tracking               :", tracking_Tanom)
		print("   -> Dry intrusion analysis                     :", DI_analysis)
		print("   -> Simulation starts at                       :", list_year[0]+list_month[0]+list_day[0]+" "+list_hour[0]+":"+list_min[0]+":00")
		print("   -> Simulation ends at                         :", list_year[-1]+list_month[-1]+list_day[-1]+" "+list_hour[-1]+":"+list_min[-1]+":00")
		print("   -> Calendar                                   :", calendar)
		print("   -> Saving full parcels data                   :", save_full_parts_position)
		print("   -> Central longitude of the out grid          :", cenlon)
		print("   -> Area conserving grid                       :",area_conserving_grid)

		if tracking_heat:
			
			if var_heat_track=="dse":
				varunits="kJ/Kg"
				varname="Dry static energy"
			elif var_heat_track=="potTemp":
				varunits="K"
				varname="Potential Temperature"
				
			if pblcheck==0:
				pblcheck_="no PBL check, use everything"
			elif pblcheck==1:
				pblcheck_="at least one location within the PBL"
			elif pblcheck==2:
				pblcheck_="both locations within the PBL"
				
			if heat_custom_limits_highs[0]==0 and heat_custom_limits_highs[1]==0:
				filterhigh="(Not Applied)"
			else:
				filterhigh=""

			print("\n   -------------    Heat Tracking Information    -------------")
			print("    + Tracking method                            :", heat_tracking_method)
			print("    + Heat tracking using                        :", varname)
			print("    + Uptake threshold                           :", ">", dvarheatthreshold, varunits)
			print("    + Filter parcels within PBL                  :", filter_pbl_parcels)
			print("    + PBL check                                  :", pblcheck_)
			if pblcheck!=0:
				print("    + PBL method                                 :", pbl_method)
			print("    + Check specific humidity                    :", dqcheck)
			if dqcheck:
				print("    + Specific humidity chage                    :", "<", dqthreshold)
			print("    + Check relative humidity                    :", trk_rh_check)
			if trk_rh_check:
				print("    + Relative humidity change                  :", "<=", rh_threshold, "%")
			print("    + Linear adjustment                          :", heat_linear_adjustment)
			print("    + Lower and upper limits for filter parcels  :", heat_custom_limits_highs, 'meters', filterhigh)
			print("    + Saving heat parcels data                   :", save_heat_parts_position)
		
		if tracking_Tanom:

	
			print("\n   -------------    Temperature anomaly Tracking Information    -------------")
			print("    + Tracking method                            :", Tanom_tracking_method)
			print("    + Path to temperature data                   :",path_clim_temperature)
			print("    + Analysis levels                            :", analysis_levels)
			print("    + Recalculating air parcel temperature       :", interpolate_parcel_temperature)
			if "sfc" in analysis_levels:
				print("    + dP to filter parcels at surface            :", dp_sfc, "hPa") 
			print("    + dP to filter parcels at upper levels       :", dp_upper, "hPa")
			print("    + Temperature anomaly threshold              :", Tanom_threshold, "K")
			if Tanom_tracking_method=="RP23":
				print("    + Linear adjustment                          :", "not required")
			else:
				print("    + Linear adjustment                          :", Tanom_linear_adjustment)
			print("    + Saving processed parcels data              :", save_Tanom_parts_position)


		if tracking_moisture or DI_moisture_tracking:
			print("\n   -------------  Moisture Tracking Information  -------------")
			print("    + Tracking method                            :", moisture_tracking_method)
			print("    + Moisture tracking using                    :", "Specific Humidity")

			if mode=="backward":

				print("    + Filter precipitating parcels               :", filter_dqdt_parcels)
				if filter_dqdt_parcels:
					print("    + Tracking only non-precip parcels           :", only_non_precip)
					print("    + dq/dt threshold                            :", "<", dqdt_threshold, "kg/kg")

			print("    + Filter parcels within PBL                  :", filter_pbl_dq_parcels)
			
			if moist_custom_limits_highs[0]==0 and moist_custom_limits_highs[1]==0:
				mfilterhigh="(Not Applied)"
			else:
				mfilterhigh=""
			if filter_pbl_dq_parcels:
				print("    + Lower and upper limits for filter parcels  :", moist_custom_limits_highs, 'meters', mfilterhigh )	
				
			if dqpblcheck==0:
				dqpblcheck_="no PBL check, use everything"
			elif dqpblcheck==1:
				dqpblcheck_="at least one location within the PBL"
			elif dqpblcheck==2:
				dqpblcheck_="both locations within the PBL"
			print("    + PBL check                                  :", dqpblcheck_)
			if dqpblcheck!=0:
				print("    + PBL method                                 :", dqpbl_method)
			
			if moisture_tracking_method!="SJ05":
				print("    + Uptake threshold                           :", ">", mindq_gain, "kg/kg")
				print("    + Losses threshold                           :","<", mindq_loss, "kg/kg" )
			print("    + Check RH for moisture uptake               :", dqtrk_rh_check)
			if dqtrk_rh_check:
				print("    + Relative humidity change                   :", "<=", dqrh_threshold, "%")
			print("    + Linear adjustment                          :", moisture_linear_adjustment)
			print("    + Minimum RH to account for precipitation    :", ">",precip_minrh, "%")
			
			print("    + Check RH for precipitation en route        :", check_RH_route_precip)
			if check_RH_route_precip:
				print("    + Minimum RH for precipitation en route      :", ">",precip_minrh_en_route, "%")
				
			print("    + Lower and upper limits for filter parcels  :", moist_custom_limits_highs, 'meters', mfilterhigh)
			print("    + Bias-correct moisture sources              :", moist_bias_correction)


			print("    + Saving moisture parcels data               :", save_moist_parts_position)

		if DI_analysis:
			print("\n   -------------  Dry Intrusion analysis  --------------------")
			print("    + DI analysis method                         :", 'Following Raveh-Rubin, S (2017) -> DOI:10.1175/JCLI-D-16-0782.1')
			print("    + Time step for checking pressure change     :", DI_dt_thresh, 'minutes')
			print("    + Time step for evaluating DI                :", DI_step,'minutes')
			print("    + Pressure change                            :", DI_dp_change,'Pa')
			print("    + Start pressure lower than                  :", DI_start_pressure_level,'Pa')
			print("    + End pressure higher than                   :", DI_end_pressure_level,'Pa')
			print("    + Range of the trajectory to evaluate DI     :", DI_time_limits_checking,'minutes')

			if DI_moisture_tracking:
				print("    + Compute moisture changes                   :", DI_moisture_tracking,'(See Moisture Tracking Information)')
				print("    + Bias-correct moisture sources              :", DI_moist_bias_correction,'(No moisture sources bias-correction is applied in DI analysis)')

			else:
				print("    + Compute moisture changes                   :", DI_moisture_tracking)
			print("    + Saving raw parcel trajectories with DI     :",DI_save_raw_parcels)

		print("-------------------------------------------------------------------------------------------------------------")


def writing_netcdf(
    latitude,
    longitude,
    tracking_heat,
    var_heat_track,
    array_heat_day,
    tracking_moisture,
    array_moist_day,
    tensor_org,
    tensor_heat,
    tensor_moist,
    CR,
    area,
    track_days,
    mode,
    run_date,
    listdates,
    meantimes,
    dtime,
    moisture_tracking_method,
    heat_tracking_method,
    moisture_linear_adjustment,
    filename,
    save_heat_parts_position,
    save_moist_parts_position,
    save_full_parts_position,
    precip_matrix,
    filter_dqdt_parcels,
    bias_cor_precip_matrix,
    array_moistday_corected,
    moist_bias_correction,
    save_Tanom_parts_position,
    Tanom_tracking_method,
    tracking_Tanom,
    analysis_levels,
    dict_output,
    DI_analysis,
    DI_dict_output,
    DI_moisture_tracking,
    DI_save_raw_parcels,
):

    try:
        os.remove(filename + ".nc")
    except OSError:
        pass

    ncout = Dataset(filename + ".nc", "w", format="NETCDF4")
    # define axis size
    ncout.createDimension("lat", len(latitude))
    ncout.createDimension("lon", len(longitude))
    ncout.createDimension("tracking_days", len(track_days))
    ncout.createDimension("time", 1)

    ncout.createDimension("ntime", tensor_org.shape[0])
    ncout.createDimension("nparcels", tensor_org.shape[1])
    ncout.createDimension("nvars", tensor_org.shape[2])

    time = ncout.createVariable("time", "f8", "time")
    time.units = "hours since 1900-01-01 00:00:00"
    time.calendar = "Standard"
    time.standard_name = "Run started at " + str(run_date)

    tracking_days = ncout.createVariable("tracking_days", "f4", "tracking_days")
    tracking_days.units = "days"

    ntime = ncout.createVariable("ntime", "f8", "ntime")
    ntime.units = "hours since 1900-01-01 00:00:00"
    ntime.calendar = "Standard"
    ntime.standard_name = "Run started at " + str(run_date)

    # create latitude axis
    lat = ncout.createVariable("lat", "f8", ("lat"), zlib=True)
    lat.standard_name = "latitude"
    lat.long_name = "latitude"
    lat.units = "degrees"
    lat.axis = "Y"

    # create longitude axis
    lon = ncout.createVariable("lon", "f8", ("lon"), zlib=True)
    lon.standard_name = "longitude"
    lon.long_name = "longitude"
    lon.units = "degrees"
    lon.axis = "X"

    gridarea = ncout.createVariable("grid_area", "f8", ("lat", "lon"), zlib=True)
    gridarea.standard_name = "Grid cell area"
    gridarea.units = "m2"

    if save_full_parts_position:
        parcels_position = ncout.createVariable(
            "raw_parcels_position",
            "f8",
            ("time", "ntime", "nparcels", "nvars"),
            zlib=True,
        )
        parcels_position.standard_name = "Raw parcels position at each time step. Only parcels within the target region at time t0"

    if DI_analysis:
        ncout.createDimension("back_times", len(DI_dict_output["DI_back_times"]))

        back_times = ncout.createVariable("back_times", "f8", "back_times")
        back_times.units = f"Minutes from  {str(run_date)}"

        di_outgrid = ncout.createVariable(
            "gridded_DI", "f8", ("time", "back_times", "lat", "lon"), zlib=True
        )
        di_outgrid.standard_name = "DI Intrusion per square meter at the analysis time"
        di_outgrid.units = "trajectories/m2"

        di_footprint = ncout.createVariable(
            "DI_footprint", "f8", ("time", "back_times", "lat", "lon"), zlib=True
        )
        di_footprint.standard_name = "Trajectories with DI Intrusion per square meter during the time step for checking DI"
        di_footprint.units = "trajectories/m2"

        if DI_save_raw_parcels:
            ncout.createDimension("DI_parcels", int(DI_dict_output[f"parcels with DI"]))

            DI_parcels_position = ncout.createVariable(
                "DI_raw_parcels",
                "f8",
                ("time", "ntime", "DI_parcels", "nvars"),
                zlib=True,
            )
            DI_parcels_position.standard_name = "Raw parcels with dry intrusion at any time of their trajectory. Only parcels within the target region at time t0"

        if DI_moisture_tracking:
            DImoistd = ncout.createVariable(
                "DI_moisture_changes",
                "f8",
                ("time", "back_times", "tracking_days", "lat", "lon"),
                zlib=True,
            )
            DImoistd.standard_name = (
                "Evaporation minus precipitation by days for each DI time"
            )
            DImoistd.units = "mm/day"

            DImoist = ncout.createVariable(
                "DI_moisture_integrated",
                "f8",
                ("time", "back_times", "lat", "lon"),
                zlib=True,
            )
            DImoist.standard_name = "Integrated moisture uptake for each DI time"
            DImoist.units = "mm/day"

            fullDImoistd = ncout.createVariable(
                "allDI_moisture_changes",
                "f8",
                ("time", "tracking_days", "lat", "lon"),
                zlib=True,
            )
            fullDImoistd.standard_name = (
                "Evaporation minus precipitation by days for all DI parcels"
            )
            fullDImoistd.units = "mm/day"

            fullnonDImoistd = ncout.createVariable(
                "allNonDI_moisture_changes",
                "f8",
                ("time", "tracking_days", "lat", "lon"),
                zlib=True,
            )
            fullnonDImoistd.standard_name = (
                "Evaporation minus precipitation by days for all non DI parcels"
            )
            fullnonDImoistd.units = "mm/day"

            fullDImoist = ncout.createVariable(
                "allDI_moisture_integrated", "f8", ("time", "lat", "lon"), zlib=True
            )
            fullDImoist.standard_name = "Integrated moisture uptake for all DI parcels"
            fullDImoist.units = "mm/day"

            fullnonDImoist = ncout.createVariable(
                "allNoNDI_moisture_integrated", "f8", ("time", "lat", "lon"), zlib=True
            )
            fullnonDImoist.standard_name = (
                "Integrated moisture uptake for all non DI parcels"
            )
            fullnonDImoist.units = "mm/day"

    if tracking_heat:

        ncout.createDimension("heat_time", tensor_heat.shape[0])
        ncout.createDimension("heat_parcels", tensor_heat.shape[1])
        ncout.createDimension("heat_vars", 3)

        heat_time = ncout.createVariable("heat_time", "f8", "heat_time")
        heat_time.units = "hours since 1900-01-01 00:00:00"
        heat_time.calendar = "Standard"

        heatd = ncout.createVariable(
            "Heat_days", "f8", ("time", "tracking_days", "lat", "lon"), zlib=True
        )
        heatd.standard_name = "Surface Sensible Heat Flux by days"
        heatd.units = "W/m2"

        heat = ncout.createVariable(
            "Heat_integrated", "f8", ("time", "lat", "lon"), zlib=True
        )
        heat.standard_name = "Integrated Surface Sensible Heat Flux"
        heat.units = "W/m2"

        if save_heat_parts_position:
            heat_parcels = ncout.createVariable(
                "heat_parcels_position",
                "f8",
                ("time", "heat_time", "heat_parcels", "heat_vars"),
                zlib=True,
            )
            heat_parcels.standard_name = "Parcels position for heat tracking"
            if var_heat_track == "sde":
                heat_parcels.units = "J/kg per " + str(int(dtime)) + " minutes"
            elif var_heat_track == "potTemp":
                heat_parcels.units = "K per " + str(int(dtime)) + " minutes"

    if tracking_Tanom and Tanom_tracking_method.upper() == "RP23":
        ncout.createDimension("trjtime", len(listdates) - 1)
        mtime = ncout.createVariable("trjtime", "f8", "trjtime")
        mtime.units = "hours since 1900-01-01 00:00:00"
        mtime.calendar = "Standard"
        mtime.standard_name = "Backward trajectories time. Run started at " + str(
            run_date
        )
        for level in analysis_levels:

            # creating variables for genesis Temperature anomalies
            vars()[f"vgenlat_{level}"] = ncout.createVariable(
                f"gen_lat_{level}",
                "f8",
                ("time", "lat", "lon"),
                zlib=True,
                fill_value=np.nan,
            )
            vars()[
                f"vgenlat_{level}"
            ].standard_name = f"latitude of genesis location (level:{level})"
            vars()[f"vgenlat_{level}"].units = "degrees north"

            vars()[f"vgenlon_{level}"] = ncout.createVariable(
                f"gen_lon_{level}",
                "f8",
                ("time", "lat", "lon"),
                zlib=True,
                fill_value=np.nan,
            )
            vars()[
                f"vgenlon_{level}"
            ].standard_name = f"longitude of genesis location (level:{level})"
            vars()[f"vgenlon_{level}"].units = "degrees north"

            vars()[f"vage_{level}"] = ncout.createVariable(
                f"age_{level}",
                "f8",
                ("time", "lat", "lon"),
                zlib=True,
                fill_value=np.nan,
            )
            vars()[
                f"vage_{level}"
            ].standard_name = f"Lagrangian age of temperature anomaly (level:{level})"
            vars()[f"vage_{level}"].units = "hours"

            vars()[f"vdist_{level}"] = ncout.createVariable(
                f"lag_dist_{level}",
                "f8",
                ("time", "lat", "lon"),
                zlib=True,
                fill_value=np.nan,
            )
            vars()[
                f"vdist_{level}"
            ].standard_name = (
                f"Lagrangian formation distance of temperature anomaly (level:{level})"
            )
            vars()[f"vdist_{level}"].units = "km"

            vars()[f"vres1_{level}"] = ncout.createVariable(
                f"res1_{level}",
                "f8",
                ("time", "lat", "lon"),
                zlib=True,
                fill_value=np.nan,
            )
            vars()[
                f"vres1_{level}"
            ].standard_name = f"temperature anomaly at genesis time (level:{level})"
            vars()[f"vres1_{level}"].units = "K"

            vars()[f"vres2_{level}"] = ncout.createVariable(
                f"res2_{level}",
                "f8",
                ("time", "lat", "lon"),
                zlib=True,
                fill_value=np.nan,
            )
            vars()[
                f"vres2_{level}"
            ].standard_name = f"numerical residual of temperature anomaly decomposition (level:{level})"
            vars()[f"vres2_{level}"].units = "K"

            vars()[f"vgenp_{level}"] = ncout.createVariable(
                f"gen_p_{level}",
                "f8",
                ("time", "lat", "lon"),
                zlib=True,
                fill_value=np.nan,
            )
            vars()[
                f"vgenp_{level}"
            ].standard_name = f"pressure of genesis location (level:{level})"
            vars()[f"vgenp_{level}"].units = "hPa"

            vars()[f"vdeltap_{level}"] = ncout.createVariable(
                f"delta_p_{level}",
                "f8",
                ("time", "lat", "lon"),
                zlib=True,
                fill_value=np.nan,
            )
            vars()[
                f"vdeltap_{level}"
            ].standard_name = (
                f"vertical air parcel displacement since genesis time (level:{level})"
            )
            vars()[f"vdeltap_{level}"].units = "hPa"

            # creating variables for Temperature anomalies properties
            # proc_properties=["lon_p","lat_p","pressure","dist_traj","T_anom","seas_i","adv_i","adiab1_i","adiab2_i","adiab3_i","diab_i"]

            vars()[f"vTanom_{level}"] = ncout.createVariable(
                f"T_anom_{level}",
                "f8",
                ("time", "trjtime", "lat", "lon"),
                zlib=True,
                fill_value=np.nan,
            )
            vars()[
                f"vTanom_{level}"
            ].standard_name = f"temperature anomaly (level:{level})"
            vars()[f"vTanom_{level}"].units = "K"

            vars()[f"vseas_{level}"] = ncout.createVariable(
                f"seas_{level}",
                "f8",
                ("time", "trjtime", "lat", "lon"),
                zlib=True,
                fill_value=np.nan,
            )
            vars()[
                f"vseas_{level}"
            ].standard_name = f"seasonality temperature anomaly (level:{level})"
            vars()[f"vseas_{level}"].units = "K"

            vars()[f"vadv_{level}"] = ncout.createVariable(
                f"adv_{level}",
                "f8",
                ("time", "trjtime", "lat", "lon"),
                zlib=True,
                fill_value=np.nan,
            )
            vars()[
                f"vadv_{level}"
            ].standard_name = f"advective temperature anomaly (level:{level})"
            vars()[f"vadv_{level}"].units = "K"

            vars()[f"vadb1_{level}"] = ncout.createVariable(
                f"adiab1_{level}",
                "f8",
                ("time", "trjtime", "lat", "lon"),
                zlib=True,
                fill_value=np.nan,
            )
            vars()[
                f"vadb1_{level}"
            ].standard_name = f"adiabatic temperature anomaly 1 (level:{level})"
            vars()[f"vadb1_{level}"].units = "K"

            vars()[f"vadb2_{level}"] = ncout.createVariable(
                f"adiab2_{level}",
                "f8",
                ("time", "trjtime", "lat", "lon"),
                zlib=True,
                fill_value=np.nan,
            )
            vars()[
                f"vadb2_{level}"
            ].standard_name = f"adiabatic temperature anomaly 2 (level:{level})"
            vars()[f"vadb2_{level}"].units = "K"

            vars()[f"vadb3_{level}"] = ncout.createVariable(
                f"adiab3_{level}",
                "f8",
                ("time", "trjtime", "lat", "lon"),
                zlib=True,
                fill_value=np.nan,
            )
            vars()[
                f"vadb3_{level}"
            ].standard_name = f"adiabatic temperature anomaly 3 (level:{level})"
            vars()[f"vadb3_{level}"].units = "K"

            vars()[f"vdiab_{level}"] = ncout.createVariable(
                f"diab_{level}",
                "f8",
                ("time", "trjtime", "lat", "lon"),
                zlib=True,
                fill_value=np.nan,
            )
            vars()[
                f"vdiab_{level}"
            ].standard_name = f"diabatic temperature anomaly _{level}"
            vars()[f"vdiab_{level}"].units = "K"

            vars()[f"vdistrj_{level}"] = ncout.createVariable(
                f"dist_traj_{level}",
                "f8",
                ("time", "trjtime", "lat", "lon"),
                zlib=True,
                fill_value=np.nan,
            )
            vars()[
                f"vdistrj_{level}"
            ].standard_name = f"great circle distance between trajectory position and event (level:{level})"
            vars()[f"vdistrj_{level}"].units = "km"

            vars()[f"avdistrj_{level}"] = ncout.createVariable(
                f"mean_dist_traj_{level}",
                "f8",
                ("time", "lat", "lon"),
                zlib=True,
                fill_value=np.nan,
            )
            vars()[
                f"avdistrj_{level}"
            ].standard_name = f"average great circle distance between trajectory position and TX1day event over the trajectory length (for all time > genesis time) (level:{level})"
            vars()[f"avdistrj_{level}"].units = "km"

            # Tanom_matrix_parts
            if save_Tanom_parts_position:
                ncout.createDimension(
                    f"{level}_parcels",
                    dict_output[f"{level}_Tanom_matrix_parts"].shape[1],
                )
                ncout.createDimension(
                    f"{level}_parcels_properties",
                    len(dict_output[f"{level}_parcels_properties"]),
                )

                vars()[f"trjparcels_{level}"] = ncout.createVariable(
                    f"trjparcels_{level}",
                    "f8",
                    (
                        "time",
                        "trjtime",
                        f"{level}_parcels",
                        f"{level}_parcels_properties",
                    ),
                    zlib=True,
                    fill_value=-999,
                )
                vars()[
                    f"trjparcels_{level}"
                ].standard_name = f"Parcels position for temperature anomaly tracking tracking (level:{level})"

    if tracking_Tanom and Tanom_tracking_method.upper() == "PR23":

        tanom_terms = ["tanom", "seas", "adv", "adiab", "diab"]
        ncout.createDimension("extend_tanom_terms", len(tanom_terms) + 3)
        ncout.createDimension("tanom_terms", len(tanom_terms))
        ncout.createDimension("trjtime", len(listdates) - 1)
        ncout.createDimension("contribs", 3)
        mtime = ncout.createVariable("trjtime", "f8", "trjtime")
        mtime.units = "hours since 1900-01-01 00:00:00"
        mtime.calendar = "Standard"
        mtime.standard_name = "Backward trajectories time. Run started at " + str(
            run_date
        )
        for level in analysis_levels:
            ncout.createDimension(
                f"{level}_parcels", dict_output[f"{level}_matrix_gen"].shape[1]
            )
            ncout.createDimension(
                f"{level}_parcels_properties",
                dict_output[f"{level}_matrix_gen"].shape[2],
            )

            vars()[f"contrib_{level}"] = ncout.createVariable(
                f"contributions_{level}",
                "f8",
                (
                    "time",
                    "contribs",
                    "tanom_terms",
                ),
                zlib=True,
            )
            vars()[
                f"contrib_{level}"
            ].standard_name = f"Contributions from each process to T' (level: {level}). Rows represent overall, positive and negative contributions. Columns represent contributions from tanom, seasonality, advective, adiabatic and diabatic processes"
            vars()[f"contrib_{level}"].units = "%"

            for term in tanom_terms:
                vars()[f"{term}dayspos{level}"] = ncout.createVariable(
                    f"{term}_days_pos_{level}",
                    "f8",
                    ("time", "tracking_days", "lat", "lon"),
                    zlib=True,
                    fill_value=np.nan,
                )
                vars()[
                    f"{term}dayspos{level}"
                ].standard_name = (
                    f"Source positive contribution for {term} (level: {level})"
                )
                vars()[f"{term}dayspos{level}"].units = "%/km2"

                vars()[f"{term}daysneg{level}"] = ncout.createVariable(
                    f"{term}_days_neg_{level}",
                    "f8",
                    ("time", "tracking_days", "lat", "lon"),
                    zlib=True,
                    fill_value=np.nan,
                )
                vars()[
                    f"{term}daysneg{level}"
                ].standard_name = (
                    f"Source positive contribution for {term} (level: {level})"
                )
                vars()[f"{term}daysneg{level}"].units = "%/km2"

                vars()[f"{term}integrated_pos{level}"] = ncout.createVariable(
                    f"{term}_integrated_pos_{level}",
                    "f8",
                    ("time", "lat", "lon"),
                    zlib=True,
                )
                vars()[
                    f"{term}integrated_pos{level}"
                ].standard_name = f"Integrated {term} for level {level}"
                vars()[f"{term}integrated_pos{level}"].units = "%/km2"

                vars()[f"{term}integrated_neg{level}"] = ncout.createVariable(
                    f"{term}_integrated_neg_{level}",
                    "f8",
                    ("time", "lat", "lon"),
                    zlib=True,
                )
                vars()[
                    f"{term}integrated_neg{level}"
                ].standard_name = f"Integrated {term} for level {level}"
                vars()[f"{term}integrated_neg{level}"].units = "%/km2"

                vars()[f"{term}gen{level}"] = ncout.createVariable(
                    f"{term}_genesis_{level}",
                    "f8",
                    ("time", f"{level}_parcels", f"{level}_parcels_properties"),
                    zlib=True,
                )
                vars()[
                    f"{term}gen{level}"
                ].standard_name = f"{term} genesis properties for level {level}, {term}_genesis_{level}[0] - gen_lon (degrees),  {term}_genesis_{level}[1] - gen_lat (degrees), {term}_genesis_{level}[2] - gen_p (Pa), {term}_genesis_{level}[3] - gen_age (hours), {term}_genesis_{level}[4] - gen_dist (km), {term}_genesis_{level}[5] - res (K). if term is Tanom, res is the overall res, if term is each of the physical processes (adiabatic, advective, diabatic or seasonality), the res correspond the pre-existent anomaly before the process genesis."

                vars()[f"{term}anom_target{level}"] = ncout.createVariable(
                    f"{term}_target_region_{level}",
                    "f8",
                    ("time", "lat", "lon"),
                    zlib=True,
                    fill_value=-999,
                )
                vars()[
                    f"{term}anom_target{level}"
                ].standard_name = (
                    f"Anomalies due to {term} at the target region for level {level} "
                )
                vars()[f"{term}anom_target{level}"].units = "K"

            # Tanom_matrix_parts

            if save_Tanom_parts_position:
                ncout.createDimension(
                    f"{level}_nparcels", dict_output[f"{level}_matrix_trajs"].shape[1]
                )
                vars()[f"process_{level}"] = ncout.createVariable(
                    f"trjs_process_{level}",
                    "f8",
                    ("time", "trjtime", f"{level}_nparcels", "extend_tanom_terms"),
                    zlib=True,
                    fill_value=-999,
                )
                vars()[
                    f"process_{level}"
                ].standard_name = f"Physical process leading to the temperature anomaly event (level:{level})"

    if tracking_moisture:

        ncout.createDimension("moisture_time", tensor_moist.shape[0])
        ncout.createDimension("moisture_parcels", tensor_moist.shape[1])
        ncout.createDimension("moisture_vars", 3)

        moisture_time = ncout.createVariable("moisture_time", "f8", "moisture_time")
        moisture_time.units = "hours since 1900-01-01 00:00:00"
        moisture_time.calendar = "Standard"

        moistd = ncout.createVariable(
            "moisture_days", "f8", ("time", "tracking_days", "lat", "lon"), zlib=True
        )
        moistd.standard_name = "Evaporation minus precipitation by days"
        moistd.units = "mm/day"

        moist = ncout.createVariable(
            "moisture_integrated", "f8", ("time", "lat", "lon"), zlib=True
        )
        moist.standard_name = "Integrated moisture uptake"
        moist.units = "mm/day"

        if moist_bias_correction:
            moistdbias = ncout.createVariable(
                "moisture_days_corrected",
                "f8",
                ("time", "tracking_days", "lat", "lon"),
                zlib=True,
            )
            moistdbias.standard_name = (
                "Bias corrected evaporation minus precipitation by days"
            )
            moistd.units = f"mm/{int(dtime/60)}h"

            moistbias = ncout.createVariable(
                "moisture_integrated_corrected", "f8", ("time", "lat", "lon"), zlib=True
            )
            moistbias.standard_name = "Bias corrected integrated moisture uptake"
            moistbias.units = f"mm/{int(dtime/60)}h"

            pmatrixbias = ncout.createVariable(
                "Lag_precip_corrected", "f8", ("time", "lat", "lon"), zlib=True
            )
            pmatrixbias.standard_name = (
                "Bias corrected estimated Lagrangian precipitation in the target region"
            )
            pmatrixbias.units = f"mm/{int(dtime/60)}h"

        if save_moist_parts_position:
            moist_parcels = ncout.createVariable(
                "moisture_parcels_position",
                "f8",
                ("time", "moisture_time", "moisture_parcels", "moisture_vars"),
                zlib=True,
            )
            moist_parcels.long_name = "Parcels position for moisture tracking"
            moist_parcels.units = "kg/kg per " + str(int(dtime)) + " minutes"

       
        if filter_dqdt_parcels:
            if mode=="backward":
                pmatrix = ncout.createVariable(
                    "Lag_precip", "f8", ("time", "lat", "lon"), zlib=True
                )
                pmatrix.standard_name = (
                    "Estimated Lagrangian precipitation in the target region"
                )
                pmatrix.units = f"mm/{int(dtime/60)}h"

    lon[:] = longitude[:]
    lat[:] = latitude[:]
    gridarea[:] = area[:]
    tensor_description = ""
    if save_full_parts_position:

        try:

            if not tensor_org.flags["C_CONTIGUOUS"]:
                tensor_org = np.ascontiguousarray(tensor_org)

            chunk_rows = 5
            num_total_rows = tensor_org.shape[0]

            for i in range(0, num_total_rows, chunk_rows):
                start_row = i
                end_row = min(i + chunk_rows, num_total_rows)

                # Define the slice for the source data
                data_slice = tensor_org[start_row:end_row, :, :]
                parcels_position[0, start_row:end_row, :, :] = data_slice

            # parcels_position[0,:]=tensor_org[:]
            tensor_description = "nvars[0]-particle id, nvars[1]-longitude, nvars[2]-latitude, nvars[3]-specific humidity(kg/kg), nvars[4]-parcel high (m), nvars[5]-topography high [m], nvar[6]-parcel density (kg/m3), nvar[7]-PBL high [m], nvar[8]-tropopause high [m], nvar[9]-parcel Temperature [K], nvar[10]-parcel mass [kg], nvar[11]-heat tracking var, nvar[12]-reltive humidity[%], nvar[13] - Pressure (Pa)"
        except:
            print(
                "\n     => WARNING: Raw partposit data is too large to be saved in the netCDF: Skipping"
            )

    time[:] = date2num(
        datetime(
            int(run_date.split(" ")[0].split("-")[0]),
            int(run_date.split(" ")[0].split("-")[1]),
            int(run_date.split(" ")[0].split("-")[2]),
            int(run_date.split(" ")[1].split(":")[0]),
            int(run_date.split(" ")[1].split(":")[1]),
            0,
        ),
        time.units,
        time.calendar,
    )

    tracking_days[:] = track_days[:]

    ntime[:] = date2num(listdates, ntime.units, ntime.calendar)

    if tracking_heat:
        heatd[0, :] = array_heat_day[:]
        heat[0, :] = np.sum(array_heat_day, axis=0)

        htensor_des = ""
        if save_heat_parts_position:
            heat_parcels[0, :] = tensor_heat[:, :, :3]
            htensor_des = (
                " method. heat_vars[0] = longitude, heat_vars[1] = latitude, heat_vars[2] = d"
                + var_heat_track
                + "/dt"
            )

        heat_time[:] = date2num(meantimes[1:], heat_time.units, heat_time.calendar)

        heat_description = "Heat Tracking using " + heat_tracking_method + htensor_des
    else:
        heat_description = "Heat Tracking no applied"

    if tracking_Tanom and Tanom_tracking_method == "RP23":
        # Filling the variables for genesis
        for level in analysis_levels:

            matrix_gen = dict_output[f"{level}_gen_matrix"]

            matrix_gen[matrix_gen == -999.0] = np.nan

            vars()[f"vgenlat_{level}"][0, :, :] = matrix_gen[1, :, :]
            vars()[f"vgenlon_{level}"][0, :, :] = matrix_gen[0, :, :]

            vars()[f"vage_{level}"][0, :, :] = matrix_gen[5, :, :]
            vars()[f"vdist_{level}"][0, :, :] = matrix_gen[6, :, :]
            vars()[f"vres1_{level}"][0, :, :] = matrix_gen[3, :, :]
            vars()[f"vres2_{level}"][0, :, :] = matrix_gen[4, :, :]
            vars()[f"vgenp_{level}"][0, :, :] = matrix_gen[2, :, :] / 100
            vars()[f"vdeltap_{level}"][0, :, :] = matrix_gen[7, :, :] / 100

            # Filling the variables for temperature anomalies properties
            Tanom_matrix_gridded = dict_output[f"{level}_Tanom_matrix_gridded"]
            Tanom_matrix_gridded[Tanom_matrix_gridded == -999] = np.nan
            # proc_properties=["lon_p","lat_p","pressure","dist_traj","T_anom","seas_i","adv_i","adiab1_i","adiab2_i","adiab3_i","diab_i"]
            mtime[:] = date2num(listdates[1:], mtime.units, mtime.calendar)

            vars()[f"vTanom_{level}"][0, :, :, :] = Tanom_matrix_gridded[4, :, :, :]
            vars()[f"vseas_{level}"][0, :, :, :] = Tanom_matrix_gridded[5, :, :, :]
            vars()[f"vadv_{level}"][0, :, :, :] = Tanom_matrix_gridded[6, :, :, :]
            vars()[f"vadb1_{level}"][0, :, :, :] = Tanom_matrix_gridded[7, :, :, :]
            vars()[f"vadb2_{level}"][0, :, :, :] = Tanom_matrix_gridded[8, :, :, :]
            vars()[f"vadb3_{level}"][0, :, :, :] = Tanom_matrix_gridded[9, :, :, :]
            vars()[f"vdiab_{level}"][0, :, :, :] = Tanom_matrix_gridded[10, :, :, :]
            vars()[f"vdistrj_{level}"][0, :, :, :] = Tanom_matrix_gridded[3, :, :, :]


            txtensor_des = ""
            if save_Tanom_parts_position:
                # mtime[:]=mtimes[:]

                try:

                    if not dict_output[f"{level}_Tanom_matrix_parts"].flags[
                        "C_CONTIGUOUS"
                    ]:
                        tensor_org = np.ascontiguousarray(
                            dict_output[f"{level}_Tanom_matrix_parts"]
                        )

                    chunk_rows = 2
                    num_total_rows = dict_output[f"{level}_Tanom_matrix_parts"].shape[0]

                    for i in range(0, num_total_rows, chunk_rows):
                        start_row = i
                        end_row = min(i + chunk_rows, num_total_rows)

                        # Define the slice for the source data
                        data_slice = dict_output[f"{level}_Tanom_matrix_parts"][
                            start_row:end_row, :
                        ]
                        vars()[f"trjparcels_{level}"][
                            0, start_row:end_row, :
                        ] = data_slice

                    # 	vars()[f'trjparcels_{level}'][0,:]=dict_output[f"{level}_Tanom_matrix_parts"]
                    txtensor_des = " method. trjparcels[0] = longitude (degrees), trjparcels[1] = latitude (degrees), trjparcels[2] = pressure (hPa), trjparcels[3] = great circle distance between trajecory position and temperature anomaly event (km), trjparcels[4] = temperature anomaly (K), trjparcels[5] = seasonality temperature anomaly (K), trjparcels[6] = advective temperature anomaly (K), trjparcels[7] = adiabatic temperature anomaly 1 (K), trjparcels[8] = adiabatic temperature anomaly 12 (K), trjparcels[9] = adiabatic temperature anomaly 3 (K), trjparcels[10] = diabatic temperature anomaly (K)"
                except:
                    print(
                        "\n     => WARNING: Parcel trajectories for Tanom analysis is too large to be saved in the netCDF: Skipping"
                    )

        tx_description = (
            "Temperature anomaly tracking using " + Tanom_tracking_method + txtensor_des
        )

    elif tracking_Tanom and Tanom_tracking_method == "PR23":
        mtime[:] = date2num(listdates[1:], mtime.units, mtime.calendar)
        tanom_terms = ["tanom", "seas", "adv", "adiab", "diab"]
        for level in analysis_levels:
            day_matrix_level = dict_output[f"{level}_matrix_day_integrated"]
            matrix_gen = dict_output[f"{level}_matrix_gen"]
            anoms_target_gridded = dict_output[f"{level}_gridded_process_anoms_mask"]
            iterm = 0
            itemgen = 0

            vars()[f"contrib_{level}"][0, :] = dict_output[
                f"{level}_contribution_matrix"
            ][:]

            for term in tanom_terms:

                vars()[f"{term}dayspos{level}"][0, :] = day_matrix_level[iterm, :]
                vars()[f"{term}daysneg{level}"][0, :] = day_matrix_level[iterm + 1, :]
                vars()[f"{term}integrated_pos{level}"][0, :] = np.sum(
                    day_matrix_level[iterm, :], axis=0
                )
                vars()[f"{term}integrated_neg{level}"][0, :] = np.sum(
                    day_matrix_level[iterm + 1, :], axis=0
                )
                iterm = iterm + 2

                vars()[f"{term}gen{level}"][0, :] = matrix_gen[itemgen, :]

                vars()[f"{term}anom_target{level}"][0, :] = anoms_target_gridded[
                    itemgen, :
                ]

                itemgen = itemgen + 1

            txtensor_des = ""
            if save_Tanom_parts_position:
                vars()[f"process_{level}"][0, :] = dict_output[f"{level}_matrix_trajs"][
                    :
                ]
                txtensor_des = " method. process_level[0] = longitude (degrees), process_level[1] = latitude (degrees), ,process_level[2] = Pressure (Pa),  process_level[3] = Tanom (K),  process_level[4] = seas (K), process_level[5] = adv (K),process_level[6] = adiab (K),process_level[7] = diab (K),"

        tx_description = (
            "Temperature anomaly tracking using " + Tanom_tracking_method + txtensor_des
        )

    else:
        tx_description = "Temperature anomaly tracking no applied"

    if tracking_moisture:
        moistd[0, :] = array_moist_day[:]
        moist[0, :] = np.sum(array_moist_day, axis=0)

        if moist_bias_correction:

            moistdbias[0, :] = array_moistday_corected[:]
            pmatrixbias[0, :] = bias_cor_precip_matrix[:]
            moistbias[0, :] = np.sum(array_moistday_corected, axis=0)

        mtensor_des = ""
        if save_moist_parts_position:
            moist_parcels[0, :] = tensor_moist[:, :, :3]
            mtensor_des = " method. moisture_vars[0] = longitude, moisture_vars[1] = latitude, moisture_vars[2] = dq/dt"

        

        if filter_dqdt_parcels:
            if mode=="backward":
                pmatrix[0, :] = precip_matrix[:]

        moisture_time[:] = date2num(
            meantimes[:-1], moisture_time.units, moisture_time.calendar
        )

        moist_description = (
            "Moisture Tracking using " + moisture_tracking_method + mtensor_des
        )
    else:
        moist_description = "Moisture Tracking no applied"

    if DI_analysis:
        back_times[:] = DI_dict_output["DI_back_times"][:]
        di_outgrid[0, :] = DI_dict_output["matrix_DIoutgrid"][:]
        di_footprint[0, :] = DI_dict_output["matrix_DIfootprint"][:]

        if DI_save_raw_parcels:
            DI_parcels_position[0, :] = DI_dict_output["raw parcels with DI"][:]

        if DI_moisture_tracking:
            DImoistd[0, :] = DI_dict_output["matrix_MUgrid"][:]
            fullDImoistd[0, :] = DI_dict_output["matrix_MUgrid with DI"][:]
            fullnonDImoistd[0, :] = DI_dict_output["matrix_MUgrid without DI"][:]

            fullDImoist[0, :] = np.sum(
                DI_dict_output["matrix_MUgrid with DI"][:], axis=0
            )
            fullnonDImoist[0, :] = np.sum(
                DI_dict_output["matrix_MUgrid without DI"][:], axis=0
            )

            DI_integrated = np.zeros((len(back_times), len(latitude), len(longitude)))
            for i, bt in enumerate(back_times):
                DI_integrated[i, :] = np.sum(
                    DI_dict_output["matrix_MUgrid"][i, :], axis=0
                )

            DImoist[0, :] = DI_integrated

    ncout.title = program_fullname()
    ncout.description = (
        tensor_description
        + ". <====> "
        + heat_description
        + ". <====> "
        + moist_description
        + ". <====> "
        + tx_description
    )
    today = datetime.now()
    ncout.history = (
        "Created " + today.strftime("%d/%m/%Y %H:%M:%S") + " using " + program_name()
    )
    ncout.institution = (
        "Environmental Physics Laboratory (EPhysLab), University of Vigo, Spain"
    )
    ncout.source = program_name() + " " + str(get_currentversion())

    ncout.close()


def saving_data(
    count_parcels,
    heat_parcels,
    no_heat_uptake_parcels,
    precipitating_parcels,
    no_evap_uptakes,
    rundates,
    fname,
    save_moist_stats,
    save_heat_stats,
    ctime,
    mode,
    heat_tracking_method,
    moist_tracking_method,
    lwvrt,
    moisture_linear_adjustment,
    moist_bias_correction,
    no_evap_uptakes_bias,
    lwvrt_bias,
    case_heat,
    case_Tanom,
    Tanom_tracking_method,
    save_Tanom_stats,
    dict_output,
    analysis_levels,
    DI_analysis,
    DI_dict_output,
    DI_moisture_tracking,
    save_DI_stats,
):
    f = open(fname, "w")
    f.write(f"+ Tool:  {program_name()} (version: {get_currentversion()})\n")
    f.write("+ Run Date: " + rundates + "\n")
    f.write("+ Run Mode: " + mode + " in time \n")
    f.write(
        "+ Number of parcels within the target region: "
        + str(int(count_parcels))
        + "\n"
    )
    if save_heat_stats:
        f.write("+ Heat tracking method: " + heat_tracking_method + "\n")
        f.write(
            "Number of filtered parcels at time t0 for heat tracking: "
            + str(int(heat_parcels))
            + " ({:.2f}".format(100 * (heat_parcels) / (count_parcels))
            + "%)\n"
        )
        f.write(
            f"Number of parcels without {case_heat} in the trajectory: "
            + str(int(no_heat_uptake_parcels))
            + " ({:.2f}".format(100 * (no_heat_uptake_parcels) / (heat_parcels))
            + "%)\n"
        )
        f.write(
            "-----------------------------------------------------------------------------------------------------------------------------\n"
        )
    if save_Tanom_stats == True and Tanom_tracking_method == "RP23":
        # print("--------")
        f.write(
            "+ Temperature anomaly tracking method: " + Tanom_tracking_method + "\n"
        )
        for level in analysis_levels:
            f.write(f"-->ANALYSIS LEVEL: {level}\n")
            f.write(
                "Number of filtered parcels at time t0 for Tanom analysis: "
                + str(int(dict_output[f"{level}_counter_part_Tanom"]))
                + " ({:.2f}".format(
                    100 * (dict_output[f"{level}_counter_part_Tanom"]) / (count_parcels)
                )
                + "%)\n"
            )
            f.write(
                f"Number of parcels without {case_Tanom}: "
                + str(int(dict_output[f"{level}_no_Tanom_parts"]))
                + " ({:.2f}".format(
                    100
                    * (dict_output[f"{level}_no_Tanom_parts"])
                    / (dict_output[f"{level}_counter_part_Tanom"])
                )
                + "%)\n"
            )
        f.write(
            "-----------------------------------------------------------------------------------------------------------------------------\n"
        )

    if save_Tanom_stats == True and Tanom_tracking_method == "PR23":
        f.write(
            "+ Temperature anomaly tracking method: " + Tanom_tracking_method + "\n"
        )
        for level in analysis_levels:
            f.write(f"-->ANALYSIS LEVEL: {level}\n")
            f.write(
                "Number of filter parcels at time t0 within the level: "
                + str(int(dict_output[f"{level}_filtered_parcels_within_level"]))
                + " ({:.2f}".format(
                    100
                    * (dict_output[f"{level}_filtered_parcels_within_level"])
                    / (count_parcels)
                )
                + "%)\n"
            )
            f.write(
                "Number of filter parcels at time t0 after apply temperature anomaly threshold: "
                + str(int(dict_output[f"{level}_filtered_parcels"]))
                + " ({:.2f}".format(
                    100
                    * (dict_output[f"{level}_filtered_parcels"])
                    / (dict_output[f"{level}_filtered_parcels_within_level"])
                )
                + "%)\n"
            )
            f.write(
                f"Number of parcels without Tanom contribution to temperature anonaly: "
                + str(int(dict_output[f"{level}_no_part_tanom"]))
                + " ({:.2f}".format(
                    100
                    * (dict_output[f"{level}_no_part_tanom"])
                    / ((dict_output[f"{level}_filtered_parcels"]))
                )
                + "%)\n"
            )
            f.write(
                f"Number of parcels without seas contribution to temperature anonaly: "
                + str(int(dict_output[f"{level}_no_part_seas"]))
                + " ({:.2f}".format(
                    100
                    * (dict_output[f"{level}_no_part_seas"])
                    / ((dict_output[f"{level}_filtered_parcels"]))
                )
                + "%)\n"
            )
            f.write(
                f"Number of parcels without adiab contribution to temperature anonaly: "
                + str(int(dict_output[f"{level}_no_part_adiab"]))
                + " ({:.2f}".format(
                    100
                    * (dict_output[f"{level}_no_part_adiab"])
                    / ((dict_output[f"{level}_filtered_parcels"]))
                )
                + "%)\n"
            )
            f.write(
                f"Number of parcels without adv contribution to temperature anonaly: "
                + str(int(dict_output[f"{level}_no_part_adv"]))
                + " ({:.2f}".format(
                    100
                    * (dict_output[f"{level}_no_part_adv"])
                    / ((dict_output[f"{level}_filtered_parcels"]))
                )
                + "%)\n"
            )
            f.write(
                f"Number of parcels without diab contribution to temperature anonaly: "
                + str(int(dict_output[f"{level}_no_part_diab"]))
                + "({:.2f}".format(
                    100
                    * (dict_output[f"{level}_no_part_diab"])
                    / ((dict_output[f"{level}_filtered_parcels"]))
                )
                + "%)\n"
            )
            f.write(f"!Contributions (% T')\n")
            f.write(f"!Term      Overall    Positive     Negative\n")
            f.write(f"!------------------------------------------\n")
            for iterm, term in enumerate(dict_output[f"{level}_tanom_terms"]):
                f.write(
                    f'!{term.ljust(5)}   {dict_output[f"{level}_overall_contributions"][iterm]:.2f}  {dict_output[f"{level}_positive_contributions"][iterm]:.2f}   {dict_output[f"{level}_negative_contributions"][iterm]:.2f}\n'
                )
        f.write(
            "-----------------------------------------------------------------------------------------------------------------------------\n"
        )

    if save_moist_stats:
        f.write("+ Moisture tracking method: " + moist_tracking_method + "\n")
        if mode=="backward":
            f.write(
                "Number of filtered parcels within the target region at time t0: "
                + str(int(precipitating_parcels))
                + " ({:.2f}".format(100 * (precipitating_parcels) / (count_parcels))
                + "%)\n"
            )
            f.write(
                "Number of parcels without moisture uptake in the trajectory: "
                + str(int(no_evap_uptakes))
                + " ({:.2f}".format(100 * (no_evap_uptakes) / (precipitating_parcels))
                + "%)\n"
            )
        if mode=="forward":
            f.write(
                "Number of parcels tracked forward after filtering in the source region: "
                + str(int(no_evap_uptakes))
                + " ({:.2f}".format(100 * (no_evap_uptakes) / (count_parcels))
                + "%)\n"
            )

        if moist_bias_correction:
            f.write(
                "!Number of parcels without moisture contribution to precipitation in the target region after the bias correction: "
                + str(int(no_evap_uptakes_bias))
                + " ({:.2f}".format(
                    100 * (no_evap_uptakes_bias) / (precipitating_parcels)
                )
                + "%)\n"
            )

        if moisture_linear_adjustment:
            f.write(
                "Lagrangian mean water vapour residence time: "
                + str(lwvrt)[0:5]
                + " days\n"
            )

        if moist_bias_correction:
            f.write(
                "Lagrangian mean water vapour residence time after the bias correction: "
                + str(lwvrt_bias)[0:5]
                + " days\n"
            )

        f.write(
            "-----------------------------------------------------------------------------------------------------------------------------\n"
        )

    if save_DI_stats:
        f.write("+ Dry Intrusion analysis\n")
        f.write("» Moisture tracking method: " + moist_tracking_method + "\n")
        back_times = DI_dict_output["DI_back_times"]
        DI_parcels_by_step = DI_dict_output["DI_parcels_by_step"]
        precip_parcels_by_step = DI_dict_output["precip_parcels_by_step"]
        no_uptake_parcels_by_step = DI_dict_output["no_uptake_parcels_by_step"]
        lmwvrt_by_step = DI_dict_output["lmwvrt_by_step"]
        for i, dt in enumerate(back_times):
            f.write(f"» Backward time: {np.abs(dt)} minutes ({np.abs(dt)/1440} days)\n")
            f.write(
                f"  Number of trajectories with DI: {DI_parcels_by_step[i]}"
                + " ({:.2f}".format(100 * (DI_parcels_by_step[i]) / count_parcels)
                + "%)\n"
            )
            if DI_moisture_tracking:
                f.write(
                    "  Number of precipitating parcels within the target region at time t0: "
                    + str(int(precip_parcels_by_step[i]))
                    + " ({:.2f}".format(
                        100 * (int(precip_parcels_by_step[i])) / (DI_parcels_by_step[i])
                    )
                    + "%)\n"
                )
                f.write(
                    "  Number of parcels without moisture uptake in the trajectory: "
                    + str(int(no_uptake_parcels_by_step[i]))
                    + " ({:.2f}".format(
                        100
                        * (no_uptake_parcels_by_step[i])
                        / (precip_parcels_by_step[i])
                    )
                    + "%)\n"
                )
                if moisture_linear_adjustment:
                    f.write(
                        "  Lagrangian mean water vapour residence time: "
                        + str(lmwvrt_by_step[i])[0:5]
                        + " days\n"
                    )

        des = ["with DI", "without DI"]

        for ds in des:

            f.write(
                f"» Computing moisture changes for all air parcel trajectories {ds}\n"
            )
            f.write(
                f"  Number of trajectories with DI: {DI_dict_output[f'parcels {ds}']}"
                + " ({:.2f}".format(
                    100 * (DI_dict_output[f"parcels {ds}"]) / count_parcels
                )
                + "%)\n"
            )

            if DI_moisture_tracking:

                f.write(
                    "  Number of precipitating parcels within the target region at time t0: "
                    + str(int(DI_dict_output[f"precip_parcels {ds}"]))
                    + " ({:.2f}".format(
                        100
                        * (DI_dict_output[f"precip_parcels {ds}"])
                        / DI_dict_output[f"parcels {ds}"]
                    )
                    + "%)\n"
                )
                f.write(
                    "  Number of parcels without moisture uptake in the trajectory: "
                    + str(int(DI_dict_output[f"no_uptake_parcels {ds}"]))
                    + " ({:.2f}".format(
                        100
                        * (DI_dict_output[f"no_uptake_parcels {ds}"])
                        / (DI_dict_output[f"precip_parcels {ds}"])
                    )
                    + "%)\n"
                )
                if moisture_linear_adjustment:
                    f.write(
                        "  Lagrangian mean water vapour residence time: "
                        + str(DI_dict_output[f"lmwvrt {ds}"])[0:5]
                        + " days\n"
                    )

        f.write(
            "-----------------------------------------------------------------------------------------------------------------------------\n"
        )

    f.write("+ Computing time: " + str(ctime) + " seconds\n")


def calc_A(resolution, lat, lon):
    """
    Calculate the area of each grid cell in m^2.

    Parameters
    ----------
    resolution : float
            The resolution of the grid in degrees.
    lat : 2D array
            The latitudes of the grid points.
    lon : 2D array
            The longitudes of the grid points.

    Returns
    -------
    area : 2D array
            The area of each grid cell in m^2.
    """
    rt = lc.earth_radius
    gr = np.pi / 180.0
    a, b = lat.shape
    area = np.empty((a - 1, b - 1))
    area[:, :] = 0
    for j in range(len(lat[0, :]) - 1):
        for i in range(len(lat[:, 0]) - 1):
            area[i, j] = np.abs(
                (gr * rt**2) * (np.sin(gr * lat[i, j]) - np.sin(gr * lat[i + 1, j]))
            ) * np.abs(resolution)
    return area


def get_area_conserving_grid(lon_start, lon_end, lat_start, lat_end, resolution):
    """
    Create an area-conserving grid with specified resolution and start/end points.

    Parameters:
    ----------
    lon_start : float
            The starting longitude.
    lon_end : float
            The ending longitude.
    lat_start : float
            The starting latitude.
    lat_end : float
            The ending latitude.
    resolution : float
            The resolution of the grid in degrees (for both latitude and longitude).

    Returns:
    -------
    lon_grid : 1D array
            The longitudes of the grid points.
    lat_grid : 1D array
            The latitudes of the grid points.
    """
    # Create longitude grid from lon_start to lon_end
    lon_grid = np.arange(lon_start, lon_end + resolution, resolution)

    # Generate latitudes using sine spacing (area-conserving)
    num_lat_points = int((lat_end - lat_start) / resolution) + 1
    lat_edges = np.linspace(
        np.sin(np.radians(lat_start)), np.sin(np.radians(lat_end)), num_lat_points
    )
    lat_grid = np.degrees(np.arcsin(lat_edges))  # Convert back to degrees

    lonn, latt = np.meshgrid(lon_grid, lat_grid)

    if lon_grid.max() > 180 + resolution:
        cenlon = "180"

    else:
        cenlon = "0"

    return latt, lonn, cenlon


def grid_point(resolution, numPdX, numPdY, lon_lower_left, lat_lower_left):
    """
    Create a grid of points given a resolution, number of points in the x and y direction, and the lower left coordinates.

    Parameters
    ----------
    resolution : float
            The resolution of the grid in degrees.
    numPdX : int
            Number of points in the x direction.
    numPdY : int
            Number of points in the y direction.
    lon_lower_left : float
            The longitude of the lower left corner of the grid.
    lat_lower_left : float
            The latitude of the lower left corner of the grid.

    Returns
    -------
    lat : 2D array
            The latitudes of the grid points.
    lon : 2D array
            The longitudes of the grid points.
    cenlon : str
            The longitude of the center of the grid, either "0" or "180".
    """
    lat_new = []
    lon_new = []
    lat_min = lat_lower_left
    lon_min = lon_lower_left
    lat_new = np.append(lat_new, lat_min)
    lon_new = np.append(lon_new, lon_min)

    for i in range(numPdY):
        lat_min = lat_min + resolution
        lat_new = np.append(lat_new, lat_min)
    for j in range(numPdX):
        lon_min = lon_min + resolution
        lon_new = np.append(lon_new, lon_min)
    lon, lat = np.meshgrid(lon_new, lat_new)

    if lon.max() > 180 + resolution:
        cenlon = "180"

    else:
        cenlon = "0"

    return lat, lon, cenlon


def grid_plot_final(lat, lon):
    """
    Generate a grid of midpoint latitudes and longitudes for plotting.

    Parameters
    ----------
    lat : 2D array
            The latitudes of the grid points.
    lon : 2D array
            The longitudes of the grid points.

    Returns
    -------
    lat_plot : 2D array
            The midpoint latitudes for the grid.
    lon_plot : 2D array
            The midpoint longitudes for the grid.
    """

    lat_new = []
    lon_new = []
    for i in range(len(lat[:, 0]) - 1):
        lat_new = np.append(lat_new, (lat[i + 1, 0] + lat[i, 0]) / 2.0)
    for j in range(len(lon[0, :]) - 1):
        lon_new = np.append(lon_new, (lon[0, j + 1] + lon[0, j]) / 2.0)
    lon_plot, lat_plot = np.meshgrid(lon_new, lat_new)
    return lat_plot, lon_plot


def generate_simulation_dates(
    ndays, cyear, cmonth, cday, chours, cminutes, calendar="366d"
):

    if not isinstance(chours, list):
        chours = [chours]
    if not isinstance(cminutes, list):
        cminutes = [cminutes]
    if not isinstance(cyear, list):
        cyear = [cyear]
    if not isinstance(cmonth, list):
        cmonth = [cmonth]
    if not isinstance(cday, list):
        cday = [cday]

    if calendar == "365d":

        date_start = datetime(int(cyear[0]), int(cmonth[0]), int(cday[0]))

        end_date = str(
            time_calc(
                str(int(cyear[0]))
                + "-"
                + str(int(cmonth[0])).zfill(2)
                + "-"
                + str(int(cday[0])).zfill(2)
                + " "
                + str(int(chours[0])).zfill(2)
                + ":"
                + str(int(cminutes[0])).zfill(2)
                + ":00",
                int(ndays * 24),
            )
        ).split(" ")[0]

        date_end = datetime(
            int(end_date.split("-")[0]),
            int(end_date.split("-")[1]),
            int(end_date.split("-")[2]),
        )

        date_list = [date_start.date() + timedelta(days=x) for x in range(ndays + 1)]
        leapyears = 0
        for tdate in date_list:

            if tdate.year % 4 == 0 and tdate.month == 2 and tdate.day == 29:
                leapyears += 1

        nhour = int((ndays + leapyears) * 24)
    else:
        nhour = int(ndays * 24)
    year = []
    mes = []
    dia = []
    hora = []
    mins = []
    array = np.arange(0, nhour, 24)
    for yy in cyear:
        yy = str(int(yy)).zfill(4)
        for i in array:
            for mm in cmonth:
                mm = str(int(mm)).zfill(2)
                for dd in cday:
                    dd = str(int(dd)).zfill(2)
                    for hh in chours:
                        for mmin in cminutes:
                            fecha = (
                                yy
                                + "-"
                                + mm
                                + "-"
                                + dd
                                + " "
                                + str(int(hh)).zfill(2)
                                + ":"
                                + str(int(mmin)).zfill(2)
                                + ":00"
                            )
                            a = str(time_calc(fecha, float(i)))
                            var1 = a.split(" ")
                            var11 = var1[0].split("-")
                            var12 = var1[1].split(":")
                            year_ = str(var11[0])
                            mes_ = str(var11[1])
                            dia_ = str(var11[2])
                            hora_ = str(var12[0])
                            minn_ = str(var12[1])

                            if (
                                calendar == "365d"
                                and int(year_) % 4 == 0
                                and int(mes_) == 2
                                and int(dia_) == 29
                            ):
                                msg = "nothing to do"
                            else:
                                year.append(year_)
                                mes.append(mes_)
                                dia.append(dia_)
                                hora.append(hora_)
                                mins.append(minn_)

    return year, mes, dia, hora, mins


def read_binaryFile_fortran(
    filename,
    type_file,
    lon_left_lower_corner,
    lat_left_lower_corner,
    lon_right_upper_corner,
    lat_right_upper_corner,
):

    if type_file == 1:
        with open(filename, "rb") as inputfile:
            a = b"".join([line for line in inputfile])
        npart = struct.unpack("iiii", a[0:16])
        npart = npart[2]
        data = RBF(
            filename,
            npart,
            lon_left_lower_corner,
            lat_left_lower_corner,
            lon_right_upper_corner,
            lat_right_upper_corner,
        )

    if type_file == 2:
        len_a = lf(filename)
        npart = ((len_a - 12) / 60) - 1
        data = RBF(
            filename,
            npart,
            lon_left_lower_corner,
            lat_left_lower_corner,
            lon_right_upper_corner,
            lat_right_upper_corner,
        )

    mask_to_keep = data[:, 0] != -999
    # ind=np.where(data[:, 0]==-999)
    # data=data[:int(ind[0][0]), :]
    # print(len(ind[0]))

    data = data[mask_to_keep]

    return data


def search_row_fortran(lista, matrix):

    matrix_ = np.array(sRow(matrix, lista, len(lista), len(matrix[:, 0])), np.float64)

    return matrix_


def read_binaryFile_fortranID(
    filename,
    type_file,
    lon_left_lower_corner,
    lat_left_lower_corner,
    lon_right_upper_corner,
    lat_right_upper_corner,
    idparts,
    fparts,
):

    if type_file == 1:
        with open(filename, "rb") as inputfile:
            a = b"".join([line for line in inputfile])
        npart = struct.unpack("iiii", a[0:16])
        npart = npart[2]
        # data= RBFid(filename,npart,lon_left_lower_corner,lat_left_lower_corner, lon_right_upper_corner,lat_right_upper_corner,  idparts, fparts)
        data = RBF(
            filename,
            npart,
            lon_left_lower_corner,
            lat_left_lower_corner,
            lon_right_upper_corner,
            lat_right_upper_corner,
        )

    if type_file == 2:
        len_a = lf(filename)
        npart = ((len_a - 12) / 60) - 1

        # data= RBFid(filename,npart, lon_left_lower_corner,lat_left_lower_corner, lon_right_upper_corner,lat_right_upper_corner,  idparts, fparts)
        data = RBF(
            filename,
            npart,
            lon_left_lower_corner,
            lat_left_lower_corner,
            lon_right_upper_corner,
            lat_right_upper_corner,
        )

    # Sort data by first column
    data = data[np.argsort(data[:, 0])]

    # Remove duplicates - keep first occurrence only
    _, unique_indices = np.unique(data[:, 0], return_index=True)

    data = data[unique_indices]

    # Filter existing data
    mask_ids = np.isin(data[:, 0], idparts)
    filtered_data = data[mask_ids]

    # Find missing IDs
    existing_ids = filtered_data[:, 0]
    missing_ids = np.setdiff1d(idparts, existing_ids)

    if len(missing_ids) > 0:
        # Create rows for missing IDs filled with -999
        missing_rows = np.full((len(missing_ids), data.shape[1]), -999.9)
        missing_rows[:, 0] = missing_ids

        # Combine existing and missing data
        data = np.vstack([filtered_data, missing_rows])

        # Sort by ID to maintain order matching idparts
        data = data[np.argsort(data[:, 0])]
    else:
        data = filtered_data

    return data


def load_mask_grid_NR(filename, maskname, maskvar_lon, maskvar_lat):

    wrfile = Dataset(filename)
    lat = wrfile.variables[maskvar_lat][:]
    lon = wrfile.variables[maskvar_lon][:]
    mask = wrfile.variables[maskname][:]

    if len(lon.shape) < 2:
        lon, lat = np.meshgrid(lon, lat)

    for i in range(0, lon.shape[0]):
        for j in range(0, lon.shape[1]):
            if lon[i, j] > 180:
                lon[i, j] = lon[i, j] - 360

    return lat, lon, mask


def funtion_interpol_mascara_old(maskvar_lat, maskvar_lon, mascara, data):

    lat_lon = np.empty((len(maskvar_lat), 2))
    lat_lon[:, 0] = maskvar_lon
    lat_lon[:, 1] = maskvar_lat

    prsInterpu = interp.NearestNDInterpolator(lat_lon, mascara)
    si = np.empty((data[:, 1].size, 2))
    si[:, 0] = data[:, 1]
    si[:, 1] = data[:, 2]
    result = prsInterpu(si)

    return result


def funtion_interpol_mascara(maskvar_lat, maskvar_lon, mascara, data):
    # Direct array stacking - eliminates intermediate array creation
    lat_lon = np.column_stack((maskvar_lon, maskvar_lat))

    # Create interpolator
    prsInterpu = interp.NearestNDInterpolator(lat_lon, mascara)

    # Direct slicing - eliminates intermediate array creation
    query_points = data[:, 1:3]

    return prsInterpu(query_points)


def desc_gz(name_file):

    with gzip.open(name_file, "rb") as f_in:
        with open(name_file[:-3], "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


def determine_id_binary_grid_NR_fortran(
    data, maskvar_lat, maskvar_lon, value_mascara, value_mask
):

    part_inmask = funtion_interpol_mascara(
        maskvar_lat, maskvar_lon, value_mascara, data
    )

    id_vector = np.array(D_id(part_inmask, value_mask, len(part_inmask)), dtype=int)

    submatrix = []
    ind = []
    for ii in id_vector:
        if ii != -999:
            submatrix = np.append(submatrix, data[ii, :])
            ind.append(ii)
    submatrix = np.reshape(submatrix, (len(ind), 11))
    return submatrix


def time_calc(init_time, h_diff):

    formatted_time = datetime.strptime(init_time, "%Y-%m-%d %H:%M:%S")
    calculated_time = formatted_time + timedelta(hours=h_diff)
    return calculated_time


def time_calcminutes(init_time, h_diff):

    formatted_time = datetime.strptime(init_time, "%Y-%m-%d %H:%M:%S")
    calculated_time = formatted_time + timedelta(minutes=h_diff)
    return calculated_time


def generate_file(mode, dtime, totaltime, fecha, path, key_gz, calendar):

    list_fecha = []
    listdates = []
    if mode == "backward":

        leapyearc = 0
        if calendar == "365d":
            array = np.arange(int(totaltime) + dtime, 0, -dtime)

            for i in array:
                a = str(time_calcminutes(fecha, float(i) * (-1))).split(" ")[0]
                if int(a.split("-")[1]) == 2 and int(a.split("-")[2]) == 29:
                    leapyearc += 1

        nhour = int(totaltime) + dtime + dtime * leapyearc

        array = np.arange(nhour, 0, -dtime)
        for i in array:
            a = str(time_calcminutes(fecha, float(i) * (-1)))
            var1 = a.split(" ")
            var11 = var1[0].split("-")
            var12 = var1[1].split(":")
            fecha_dia = str(
                var11[0] + var11[1] + var11[2] + var12[0] + var12[1] + var12[2]
            )
            name = path + "partposit_" + fecha_dia
            if key_gz:
                name = path + "partposit_" + fecha_dia + ".gz"
            else:
                name = path + "partposit_" + fecha_dia

            if (
                calendar == "365d"
                and int(var11[0]) % 4 == 0
                and int(var11[1]) == 2
                and int(var11[2]) == 29
            ):
                msg = "nothong to do"
            else:

                list_fecha = np.append(list_fecha, name)

                listdates = np.append(listdates, int(fecha_dia))
        fecha_ = fecha.split(" ")
        var11 = fecha_[0].split("-")
        var12 = fecha_[1].split(":")
        fecha_dia = str(var11[0] + var11[1] + var11[2] + var12[0] + var12[1] + var12[2])
        if key_gz:
            name = path + "partposit_" + fecha_dia + ".gz"
        else:
            name = path + "partposit_" + fecha_dia

        list_fecha = np.append(list_fecha, name)
        listdates = np.append(listdates, int(fecha_dia))

    if mode == "forward":

        leapyearc = 0
        if calendar == "365d":
            array = np.arange(0, int(totaltime) + dtime, dtime)

            for i in array:
                a = str(time_calcminutes(fecha, float(i) * (1))).split(" ")[0]
                if int(a.split("-")[1]) == 2 and int(a.split("-")[2]) == 29:
                    leapyearc += 1

        nhour = int(totaltime) + dtime + dtime * leapyearc
        array = np.arange(0, nhour + dtime, dtime)

        for i in array:
            a = str(time_calcminutes(fecha, float(i)))


            var1 = a.split(" ")
            var11 = var1[0].split("-")
            var12 = var1[1].split(":")
            fecha_dia = str(
                var11[0] + var11[1] + var11[2] + var12[0] + var12[1] + var12[2]
            )
            name = path + "partposit_" + fecha_dia
            if key_gz:
                name = path + "partposit_" + fecha_dia + ".gz"
            else:
                name = path + "partposit_" + fecha_dia


            if (
                calendar == "365d"
                and int(var11[0]) % 4 == 0
                and int(var11[1]) == 2
                and int(var11[2]) == 29
            ):
                msg = "nothong to do"
            else:

                list_fecha = np.append(list_fecha, name)
                listdates = np.append(listdates, int(fecha_dia))

    return list_fecha, listdates


def read_proccesor(
    verbose,
    partpositfiles,
    submatrix,
    rank,
    lon_left_lower_corner,
    lat_left_lower_corner,
    lon_right_upper_corner,
    lat_right_upper_corner,
    model,
    key_gz,
    type_file,
):

    a1 = np.arange(len(partpositfiles))
    dx, dy = submatrix.shape
    tensor_local = np.ones((len(partpositfiles), dx, dy)) * (-999.9)
    for i in a1:
        if verbose:
            print("Reading | " + model + " -> ", partpositfiles[i])
        if key_gz:
            desc_gz(partpositfiles[i])
            matrix_i = read_binaryFile_fortranID(
                partpositfiles[i][:-3],
                type_file,
                lon_left_lower_corner,
                lat_left_lower_corner,
                lon_right_upper_corner,
                lat_right_upper_corner,
                submatrix[:, 0],
                len(submatrix[:, 0]),
            )
            cmd_rm = "rm -rf " + partpositfiles[i][:-3]
            os.system(cmd_rm)
        else:
            matrix_i = read_binaryFile_fortranID(
                partpositfiles[i],
                type_file,
                lon_left_lower_corner,
                lat_left_lower_corner,
                lon_right_upper_corner,
                lat_right_upper_corner,
                submatrix[:, 0],
                len(submatrix[:, 0]),
            )

            # part_post_i=read_binaryFile_fortran(partpositfiles[-2], type_file, lon_left_lower_corner,lat_left_lower_corner,
            # 	lon_right_upper_corner,lat_right_upper_corner)

        # print((submatrix[:,0] - matrix_i[:,0]).min(), (submatrix[:,0] - matrix_i[:,0]).max())
        # matrix_i = matrix_i[np.argsort(matrix_i[:, 0])]
        # print(submatrix[:,0] - matrix_i[:,0])

        tensor_local[i, :, :] = matrix_i
    return tensor_local


def determine_id_target_region(
    data, maskvar_lat, maskvar_lon, value_mascara, value_mask
):
    """
    Final optimized version that completely replaces the Fortran D_id call
    This should be significantly faster than the original implementation
    """
    # Get interpolated mask values
    part_inmask = funtion_interpol_mascara(
        maskvar_lat, maskvar_lon, value_mascara, data
    )

    # Direct boolean indexing - most efficient approach
    # This replaces both D_id and the subsequent filtering in one step
    valid_mask = part_inmask == value_mask

    # Extract valid rows directly
    submatrix = data[valid_mask, :]

    return submatrix


def get_vars_from_partposit(
    verbose,
    partpositfiles,
    file_mask,
    maskname,
    maskvar_lon,
    maskvar_lat,
    lat_f,
    lon_f,
    rank,
    size,
    comm,
    type_file,
    lon_left_lower_corner,
    lat_left_lower_corner,
    lon_right_upper_corner,
    lat_right_upper_corner,
    model,
    value_mask,
    key_gz,
    var_heat_track,
    mode,
):
    if mode=="backward":
        name_file = partpositfiles[-1]
    elif mode=="forward":
        name_file = partpositfiles[0]
    if rank == 0:
        if verbose:
            print("Reading | " + model + " -> ", name_file)
        name_txt_part = name_file.split("/")
    if key_gz:
        desc_gz(name_file)
        part_post = read_binaryFile_fortran(
            name_file[:-3],
            type_file,
            lon_left_lower_corner,
            lat_left_lower_corner,
            lon_right_upper_corner,
            lat_right_upper_corner,
        )
        cmd_rm = "rm -rf " + name_file[:-3]
        os.system(cmd_rm)
    else:
        part_post = read_binaryFile_fortran(
            name_file,
            type_file,
            lon_left_lower_corner,
            lat_left_lower_corner,
            lon_right_upper_corner,
            lat_right_upper_corner,
        )

    lat_masked, lon_masked, mascara = load_mask_grid_NR(
        file_mask, maskname, maskvar_lon, maskvar_lat
    )

    # submatrix=determine_id_binary_grid_NR_fortran(part_post, lat_masked.flatten(), lon_masked.flatten(), mascara.flatten(), value_mask)
    submatrix = determine_id_target_region(
        part_post,
        lat_masked.flatten(),
        lon_masked.flatten(),
        mascara.flatten(),
        value_mask,
    )

    submatrix = submatrix[np.argsort(submatrix[:, 0])]

    _, unique_indices = np.unique(submatrix[:, 0], return_index=True)
    submatrix = submatrix[unique_indices]

    # print(submatrix.shape)
    # quit()

    # print(submatrix.shape)

    if mode=="backward":
        name_file2 = partpositfiles[-2]
    elif mode=="forward":
        name_file2 = partpositfiles[1]


    if rank == 0:
        if verbose:
            print("Reading | " + model + " -> ", name_file2)

    if key_gz:
        desc_gz(name_file2)
        matrix_i = read_binaryFile_fortranID(
            name_file2[:-3],
            type_file,
            lon_left_lower_corner,
            lat_left_lower_corner,
            lon_right_upper_corner,
            lat_right_upper_corner,
            submatrix[:, 0],
            len(submatrix[:, 0]),
        )
        cmd_rm = "rm -rf " + name_file2[:-3]
        os.system(cmd_rm)
    else:
        matrix_i = read_binaryFile_fortranID(
            name_file2,
            type_file,
            lon_left_lower_corner,
            lat_left_lower_corner,
            lon_right_upper_corner,
            lat_right_upper_corner,
            submatrix[:, 0],
            len(submatrix[:, 0]),
        )

        # part_post_i=read_binaryFile_fortran(partpositfiles[-2], type_file, lon_left_lower_corner,lat_left_lower_corner,
        # 		lon_right_upper_corner,lat_right_upper_corner)

    # matrix_i = matrix_i[np.argsort(matrix_i[:, 0])]

    # matrix_i=search_row_fortran(submatrix[:,0],part_post_i)

    # print("khhgfhgdgfdfgdf")
    # print(submatrix[:,0] - matrix_ii[:,0])
    # print((submatrix[:,0] - matrix_i[:,0]).min(), (submatrix[:,0] - matrix_i[:,0]).max())
    # tt2=matrix_i[:,0]
    # print("==================", len(tt2[tt2>0]))
    # print()

    # matrix_i[matrix_i==-999]=-999.9

    # print()
    # print(matrix_i)
    # print()
    # print("....................................................")
    # print()
    # print(matrix_ii)

    # print()

    # tt=matrix_ii[:,0]
    # print("-------------------",len(tt[tt>0]))
    # tt2=matrix_i[:,0]
    # print("==================", len(tt2[tt2>0]))

    # a = matrix_i - matrix_ii

    # print(a)

    # print(a[:,0].min(), a[:,0].max())

    # print(matrix_i.shape)
    # print(matrix_ii.shape)
    # quit()
    dimX, dimY = submatrix.shape

    tensor_org = np.ones((len(partpositfiles), dimX, dimY)) * (-999.9)
    tensor_org[-1, :, :] = submatrix
    tensor_org[-2, :, :] = matrix_i

    n = len(partpositfiles) - 2
    count = n // size
    remainder = n % size

    if rank < remainder:
        start = rank * (count + 1)
        stop = start + count + 1
    else:
        start = rank * count + remainder
        stop = start + count

    if mode=="forward":
       partpositfiles=partpositfiles[2:]
       partpositfiles=partpositfiles[::-1]
    local_list = partpositfiles[start:stop]
    # local_list=partpositfiles[:-2][rank::size]
    local_results = np.empty((len(local_list), dimX, dimY))
    local_results = read_proccesor(
        verbose,
        local_list,
        submatrix,
        rank,
        lon_left_lower_corner,
        lat_left_lower_corner,
        lon_right_upper_corner,
        lat_right_upper_corner,
        model,
        key_gz,
        type_file,
    )

    if rank > 0:
        comm.Send(local_results, dest=0, tag=14)
    else:
        i_start = []
        i_stop = []
        for i in range(size):
            if i < remainder:
                i_start = np.append(i_start, i * (count + 1))
                i_stop = np.append(i_stop, i_start + count + 1)
            else:
                ii_start = i * count + remainder
                ii_stop = ii_start + count
                i_start = np.append(i_start, ii_start)
                i_stop = np.append(i_stop, ii_stop)
        final_results = np.copy(local_results)

        tensor_org[int(i_start[0]) : int(i_stop[0]), :, :] = final_results
        for i in range(1, size):
            if i < remainder:
                rank_size = count + 1
            else:
                rank_size = count
            tmp = np.empty(
                (rank_size, final_results.shape[1], final_results.shape[2]),
                dtype=np.float64,
            )

            comm.Recv(tmp, source=i, tag=14)
            tensor_org[int(i_start[i]) : int(i_stop[i]), :, :] = tmp
    comm.Bcast(tensor_org, root=0)

    tensor_final = np.empty(
        (tensor_org.shape[0], tensor_org.shape[1], tensor_org.shape[2] + 3)
    )
    tensor_final[:] = -999.9
    tensor_final[:, :, :-3] = tensor_org

    if var_heat_track == "dse":

        tensor_final[:, :, -3] = compute_dry_static_energy(
            tensor_org[:, :, 2], tensor_org[:, :, 9], tensor_org[:, :, 4]
        )
        tensor_final[:, :, -3] = tensor_final[:, :, -3] / 1000

    elif var_heat_track == "potTemp":

        tensor_final[:, :, -3] = compute_theta(
            tensor_org[:, :, 6], tensor_org[:, :, 3], tensor_org[:, :, 9]
        )
    else:
        tensor_final[:, :, -3] = compute_theta(
            tensor_org[:, :, 6], tensor_org[:, :, 3], tensor_org[:, :, 9]
        )

    tensor_final[:, :, -2] = compute_rh(
        tensor_org[:, :, 6], tensor_org[:, :, 3], tensor_org[:, :, 9]
    )
    tensor_final[:, :, -1] = calc_pres(
        tensor_org[:, :, 6], tensor_org[:, :, 3], tensor_org[:, :, 9]
    )

    tensor_final[:, :, -3][tensor_final[:, :, -3] < 0] = -999.9
    tensor_final[:, :, -2][tensor_final[:, :, -2] < 0] = -999.9
    tensor_final[:, :, -1][tensor_final[:, :, -1] < 0] = -999.9

    comm.Barrier()

    if rank == 0:
        print(f"\n     => Data reading: Finished")
    # print(tensor_final[-1,0,:])
    # del tensor_org
    gc.collect()
    return tensor_final


def compute_rh(rho_kgm3, q_kgkg, T_K, press_data=False, p_Pa=np.array([None])):

    if not press_data:
        p_Pa = calc_pres(rho_kgm3, q_kgkg, T_K)

    e = q_kgkg * p_Pa / (0.622 + 0.378 * q_kgkg)

    es = 611.2 * np.exp(17.67 * (T_K - lc.TREF) / (T_K - lc.TREF + 243.5))

    return 1e2 * e / es


def calc_pres(rho_kgm3, q_kgkg, T_K):
    r_kgkg = -q_kgkg / (q_kgkg - 1)
    Tv_K = T_K * (1 + r_kgkg / lc.EPSILON) / (1 + r_kgkg)
    return rho_kgm3 * lc.RSPECIFIC * Tv_K


def calc_pottemp(p_Pa, q_kgkg, T_K):
    basic_theta = False
    if np.all(q_kgkg < 0) or np.any(np.isnan(q_kgkg)):
        basic_theta = True

    if basic_theta:
        return T_K * ((lc.p0 / p_Pa) ** lc.kappa)
    else:

        r_kgkg = -q_kgkg / (q_kgkg - 1)
        r_gkg = r_kgkg * 1e3

        p_hPa = p_Pa / 1e2
        return T_K * (1000 / p_hPa) ** (0.2854 * (1 - 0.00028 * r_gkg))


def ajust_units(array, area, density, dtime, varid, var_heat_track):
    if varid == 0:
        if var_heat_track == "dse":
            array = 1000 * array * density / (area * dtime * 60)
        elif var_heat_track == "potTemp":
            array = lc.cp * array * density / (area * dtime * 60)
    elif varid == 1:
        array = array * density / (area)
    return array


def cal_track_diff(var_vals, case):

    var_vals[var_vals < -500] = np.nan
    # difference
    if case in ["diff"]:
        dvalue = var_vals[1:] - var_vals[:-1]
    # mean
    if case in ["mean"]:
        dvalue = (var_vals[1:] + var_vals[:-1]) / 2
    # max
    if case in ["max"]:
        dvalue = np.amax(var_vals[1:], var_vals[:-1])
    return dvalue


############ SPECIFIC FUNCTIONS FOR HEAT TRACKING  ################################


def compute_dry_static_energy(lat, T_K, heigh):
    dse = (
        lc.cp * T_K
        + lc.g0
        * (
            1
            + 0.0053024 * (np.sin(lat * np.pi / 180)) ** 2
            - 0.0000058 * (np.sin(2 * lat * np.pi / 180)) ** 2
        )
        * heigh
    )

    return dse


def compute_theta(rho_kgm3, q_kgkg, T_K, press_data=False, parpres=np.array([None])):

    if not press_data:
        parpres = calc_pres(rho_kgm3, q_kgkg, T_K)

    theta = calc_pottemp(parpres, q_kgkg, T_K)

    return theta


def compute_var_integarated_day_heat(
    array,
    t,
    area,
    density,
    dtime,
    var_heat_track,
    lon,
    lat,
    numPdY,
    numPdX,
    cenlon,
    varid,
):
    dimX, dimY = lat.shape
    array_day = np.empty((len(t) - 1, dimX - 1, dimY - 1))
    ndb = np.arange(len(t) - 1, 0, -1)

    for ii in range(len(t) - 1):
        heatd = np.array(
            compute_grid_integrated_heat(
                array[t[ii] : t[ii + 1], :, :],
                lon,
                lat,
                numPdY,
                numPdX,
                len(array[t[ii] : t[ii + 1], 0, 0]),
                len(array[0, :, 0]),
            ),
            dtype=np.float64,
        )
        heatd[np.isnan(heatd)] = 0
        array_day[int(ndb[ii] - 1), :, :] = ajust_units(
            heatd, area, density, dtime, varid, var_heat_track
        )

    return array_day


def is_withinpbl(parts_high, pblhigh, lower_limit, upper_limit, method):
    check_pblupper = np.empty(len(parts_high))
    check_pbllower = np.empty(len(parts_high))

    pblhighs = []

    for i in range(0, pblhigh.shape[1]):
        pbl_vars = pblhigh[:, i]

        if method == "meanval":
            pblval = np.mean(pbl_vars[-4:])
            pblhighs = np.append(pblhighs, pblval)
        elif method == "maxval":
            pblval = np.max(pbl_vars[-4:])

            pblhighs = np.append(pblhighs, pblval)
        elif method == "actualval":
            pblval = pbl_vars[-1]
            pblhighs = np.append(pblhighs, pblval)
        else:
            print_error_message(
                " pbl_method or dqpbl_method is not valid. This parameter must be equal to  ['maxval'/'meanval'/'actualval']"
            )

    if upper_limit != 0:
        pblhighs[pblhighs <= upper_limit] = upper_limit
        pblhighs[pblhighs >= upper_limit] = upper_limit

    check_pblupper[parts_high <= pblhighs] = True
    check_pblupper[parts_high > pblhighs] = False

    check_pbllower[parts_high >= lower_limit] = True
    check_pbllower[parts_high < lower_limit] = False

    check_pbl = np.logical_and(check_pbllower, check_pblupper)

    return check_pbl


def maxval(x, n=2):

    if len(x) == 2:
        return max(x[0], x[1])
    else:
        return np.array([np.max(x[i : i + n]) for i in range(len(x) - (n - 1))])


def meanval(x, n=2):
    if len(x) == 2:
        return np.mean([x[0], x[1]])
    else:
        return np.array([np.mean(x[i : i + n]) for i in range(len(x) - (n - 1))])


def is_in_pbl(pblcheck, par_vals, pbl_vars, method, lendvar):

    if not pblcheck in (0, 1, 2):
        print_error_message(
            "pblcheck is not valid. This parameter must be equal to  0: no PBL check is applied, use everything, 1:  at least one location is within the PBL\n       2: both locations are within the max PBL "
        )

    if pblcheck == 0:
        check_pbl = np.ones((lendvar), dtype=bool)
        check_pbl[:] = True
        return check_pbl

    else:
        if method == "meanval":
            before = par_vals[:-1] <= meanval(pbl_vars, n=2)
            after = par_vals[1:] <= meanval(pbl_vars, n=2)
        elif method == "maxval":
            before = par_vals[:-1] <= maxval(pbl_vars, n=2)
            after = par_vals[1:] <= maxval(pbl_vars, n=2)

        elif method == "actualval":
            before = par_vals[:-1] <= pbl_vars[:-1]
            after = par_vals[1:] <= pbl_vars[1:]
        else:
            print_error_message(
                " pbl_method or dqpbl_method is not valid. This parameter must be equal to  ['maxval'/'meanval'/'actualval']"
            )
        if pblcheck == 2:
            return np.logical_and(before, after)
        elif pblcheck == 1:
            return np.logical_or(before, after)


def processing_heat_track_backward(
    tensor_org,
    pblcheck,
    filter_pbl_parcels,
    pbl_method,
    heat_custom_limits,
    trk_rh_check,
    rh_threshold,
    dqcheck,
    dqthreshold,
    dvar_threshold,
    var_heat_track,
    lag_times,
    area,
    density,
    dtime,
    heat_linear_adjustment,
    lon,
    lat,
    numPdY,
    numPdX,
    cenlon,
    varid,
    rank,
    size,
    comm,
):

    if filter_pbl_parcels:

        check_pbl = is_withinpbl(
            tensor_org[-1, :, 4],
            tensor_org[1:, :, 7],
            heat_custom_limits[0],
            heat_custom_limits[1],
            "maxval",
        )

        if len(check_pbl[check_pbl == True]) == 0:
            tensor_heat = np.empty((tensor_org.shape[0] - 1, 1, tensor_org.shape[2]))
            tensor_heat[:, :, :] = 0

        else:

            tensor_heat = np.empty(
                (
                    tensor_org.shape[0] - 1,
                    len(check_pbl[check_pbl == True]),
                    tensor_org.shape[2],
                )
            )

            for i in range(0, tensor_heat.shape[0]):
                tensor_heat[i, :] = tensor_org[i + 1, :, :][check_pbl == True]
    else:
        tensor_heat = np.copy(tensor_org[1:, :, :])

    matrix_heat_files = ["lon", "lat", "dq", "ds", "rh", "FH"]
    dmatrix = np.empty(
        (tensor_heat.shape[0] - 1, tensor_heat.shape[1], len(matrix_heat_files) + 1)
    )
    dmatrix[:] = 0

    # paralell_heat_process(dmatrix, tensor_heat, pblcheck, pbl_method, trk_rh_check, rh_threshold, dqcheck, dqthreshold, dvar_threshold, var_heat_track, heat_linear_adjustment)

    n = tensor_heat.shape[1]
    count = n // size
    remainder_ = n % size
    remainder = 0

    if rank < remainder:
        start = rank * (count + 1)
        stop = start + count + 1
    else:
        start = rank * count + remainder
        stop = start + count

    local_parts = tensor_heat[0, start:stop, 0]
    local_results = np.empty((dmatrix.shape[0], len(local_parts), dmatrix.shape[2]))
    local_results[:, :, :] = 0

    local_results = parallel_heat_process_backward(
        local_results,
        tensor_heat[:, start:stop, :],
        pblcheck,
        pbl_method,
        trk_rh_check,
        rh_threshold,
        dqcheck,
        dqthreshold,
        dvar_threshold,
        var_heat_track,
        heat_linear_adjustment,
        cenlon,
    )

    if rank > 0:
        comm.Send(local_results, dest=0, tag=14)
    else:
        i_start = []
        i_stop = []
        for i in range(size):
            if i < remainder:
                i_start = np.append(i_start, i * (count + 1))
                i_stop = np.append(i_stop, i_start + count + 1)
            else:
                ii_start = i * count + remainder
                ii_stop = ii_start + count
                i_start = np.append(i_start, ii_start)
                i_stop = np.append(i_stop, ii_stop)

        final_results = np.copy(local_results)

        dmatrix[:, int(i_start[0]) : int(i_stop[0]), :] = final_results

        for i in range(1, size):

            if i < remainder:
                rank_size = count + 1
            else:
                rank_size = count

            tmp = np.empty(
                (final_results.shape[0], rank_size, final_results.shape[2]),
                dtype=np.float64,
            )

            comm.Recv(tmp, source=i, tag=14)

            dmatrix[:, int(i_start[i]) : int(i_stop[i]), :] = tmp

    comm.Bcast(dmatrix, root=0)

    if remainder_ > 0:
        dmatrix[:, count * size :] = parallel_heat_process_backward(
            dmatrix[:, count * size :],
            tensor_heat[:, count * size :, :],
            pblcheck,
            pbl_method,
            trk_rh_check,
            rh_threshold,
            dqcheck,
            dqthreshold,
            dvar_threshold,
            var_heat_track,
            heat_linear_adjustment,
            cenlon,
        )

    matrix_heat = np.copy(dmatrix[:, :, :-1])
    matrix_heat[:, :, 3][matrix_heat[:, :, 3] == -999.9] = 0

    parts_uptakes = dmatrix[0, :, 6]

    array_day = compute_var_integarated_day_heat(
        matrix_heat,
        lag_times,
        area,
        density,
        dtime,
        var_heat_track,
        lon,
        lat,
        numPdY,
        numPdX,
        cenlon,
        varid,
    )

    counter_part = matrix_heat.shape[1]
    no_uptakes_parts = counter_part - int(np.sum(parts_uptakes))

    # return array_day/(lag_times[1]-lag_times[0]), counter_part
    return array_day, matrix_heat, counter_part, int(no_uptakes_parts)


def parallel_heat_process_backward(
    dmatrix,
    tensor_heat,
    pblcheck,
    pbl_method,
    trk_rh_check,
    rh_threshold,
    dqcheck,
    dqthreshold,
    dvar_threshold,
    var_heat_track,
    heat_linear_adjustment,
    cenlon,
):

    for i in range(0, tensor_heat.shape[1]):

        lons = tensor_heat[:, i, 1]

        dlon = compute_mean_lon(lons, cenlon)
        if cenlon == "180":
            dlon = (dlon + 360) % 360

            dlon[dlon < 0] = -999.9

        lats = tensor_heat[:, i, 2]
        dlat = cal_track_diff(lats, "mean")

        # dlon, dlat = compute_mean_lon_v2(lons, lats, cenlon)

        qs = tensor_heat[:, i, 3]
        dq = cal_track_diff(qs, "diff")

        var = tensor_heat[:, i, 11]
        dvar = cal_track_diff(var, "diff")

        rh = tensor_heat[:, i, 12]
        drh = cal_track_diff(rh, "diff")

        dmatrix[:, i, 0] = dlon
        dmatrix[:, i, 1] = dlat
        dmatrix[:, i, 2] = dq
        dmatrix[:, i, 3] = dvar
        dmatrix[:, i, 4] = drh
        dmatrix[:, i, 5] = 1

        trk_drh = np.ones((len(dlon)), dtype=bool)

        trk_check_pbl = is_in_pbl(
            pblcheck, tensor_heat[:, i, 4], tensor_heat[:, i, 7], pbl_method, len(dlon)
        )

        if trk_rh_check:

            trk_drh[np.abs(drh) <= rh_threshold] = True
            trk_drh[np.abs(drh) > rh_threshold] = False
        else:
            trk_drh[:] = True

        change_dq = np.ones((len(dlon)), dtype=bool)
        if dqcheck:

            changedq = np.abs(dq) / qs[:-1]
            change_dq[changedq <= dqthreshold] = True
            change_dq[changedq > dqthreshold] = False

        else:
            change_dq[:] = True

        check_dvar = np.ones((len(dlon)), dtype=bool)

        check_dvar[dvar > dvar_threshold] = True
        check_dvar[dvar <= dvar_threshold] = False

        valid = trk_check_pbl & trk_drh & change_dq & check_dvar

        dmatrix[:, i, 5][valid == True] = 1
        dmatrix[:, i, 5][valid == False] = 0

        # print(dmatrix.max())

        dmatrix[np.isnan(dmatrix)] = -999.9

        if heat_linear_adjustment:

            dmatrix[:, i, 3] = compute_linear_discounted(dmatrix[:, i, 3], valid)

            if np.sum((dmatrix[:, i, 3])) > 0:
                dmatrix[:, i, 6] = 1
        else:
            if len(valid[valid == True]) >= 1:
                # uptakes_parts = uptakes_parts + 1
                dmatrix[:, i, 6] = 1

    return dmatrix


############ SPECIFIC FUNCTIONS FOR MOISTURE TRACKING  ################################


def is_precipitating_parcel(parts_dq, dqdt_threshold, parts_rh, minrh,is_non_precip):
    check_precipdq = np.empty(len(parts_dq))
    check_preciprh = np.empty(len(parts_dq))

    if is_non_precip:
        check_precipdq[parts_dq < dqdt_threshold] = False
        check_precipdq[parts_dq >= dqdt_threshold] = True

        check_preciprh[parts_rh < minrh] = True
        check_preciprh[parts_rh >= minrh] = False
    else:
        check_precipdq[parts_dq < dqdt_threshold] = True
        check_precipdq[parts_dq >= dqdt_threshold] = False

        check_preciprh[parts_rh < minrh] = False
        check_preciprh[parts_rh >= minrh] = True

    check_precip = np.logical_and(check_precipdq, check_preciprh)

    return check_precip


def compute_var_integarated_day_moist(
    array, t, area, density, dtime, lon, lat, numPdY, numPdX, cenlon, varid
):
    dimX, dimY = lat.shape

    array_day = np.empty((len(t) - 1, dimX - 1, dimY - 1))

    array_day_cr = np.empty((len(t) - 1, dimX - 1, dimY - 1))

    ndb = np.arange(len(t) - 1, 0, -1)
    for ii in range(len(t) - 1):

        moistd = np.array(
            compute_grid_integrated_moist(
                array[t[ii] : t[ii + 1], :, :],
                lon,
                lat,
                numPdY,
                numPdX,
                len(array[t[ii] : t[ii + 1], 0, 0]),
                len(array[0, :, 0]),
            ),
            dtype=np.float64,
        )
        moistd[np.isnan(moistd)] = 0
        array_day[int(ndb[ii] - 1), :, :] = ajust_units(
            moistd, area, density, dtime, varid, "None"
        )
        array_day_cr[int(ndb[ii] - 1), :, :] = moistd
    return array_day, array_day_cr


def compute_precipitation_target_region(
    lat, lon, tensor, dqdt_threshold, cenlon, area, density
):
    # Following Kenue et al. (2022)   - implemented
    # DOI: https://doi.org/10.5194/gmd-15-1875-2022
    dq = tensor[-1, :, 3] - tensor[-2, :, 3]
    pres = tensor[-1, :, 13]
    plons = tensor[-1, :, 1]
    plats = tensor[-1, :, 2]

    if cenlon == "180":
        plons = (plons + 360) % 360

    fmatrix = np.zeros((lat.shape[0] - 1, lon.shape[1] - 1))

    # Creating boolean masks for the conditions
    dq_mask = (-999 < dq) & (dq <= dqdt_threshold)

    for k in range(lat.shape[0] - 1):
        for l in range(lon.shape[1] - 1):

            lat_mask = (lat[k, l] <= plats) & (plats < lat[k + 1, l])
            lon_mask = (lon[k, l] <= plons) & (plons < lon[k, l + 1])

            mask = dq_mask & lat_mask & lon_mask

            if np.any(mask):
                # fmatrix[k, l] = (1 / lc.g0)* np.sum((-1) * dq[mask] * pres[mask])
                fmatrix[k, l] = np.sum((-1) * dq[mask])

    fmatrix = fmatrix * density / area

    return fmatrix, dq


def correct_sink_precip(
    LagP,
    tensor,
    density,
    area,
    lat,
    lon,
    cenlon,
    precipfile,
    precip_var,
    precip_lat,
    precip_lon,
    maskfile,
    maskvar,
    maskval,
    masklat,
    masklon,
):

    dq = tensor[-1, :, 3] - tensor[-2, :, 3]
    plons = tensor[-1, :, 1]
    plats = tensor[-1, :, 2]

    lat_plot, lon_plot = grid_plot_final(lat, lon)
    points_output = np.array([lon_plot.flatten(), lat_plot.flatten()])

    # Reading mask file
    mfile = Dataset(maskfile)
    latm = mfile.variables[masklat][:]
    lonm = mfile.variables[masklon][:]
    mask = mfile.variables[maskvar][:]

    if len(mask.shape) > 2:
        mask = mas[0, :]

    mask[mask != maskval] = 0

    if len(lonm.shape) < 2:
        lonm, latm = np.meshgrid(lonm, latm)

    if cenlon == "0":
        lonm = np.where(lonm >= 180, lonm - 360, lonm)
    elif cenlon == "180":
        lonm = np.where(lonm < 0, lonm + 360, lonm)

    #####Regridding the mask file to the output grid

    regmask = griddata(
        (lonm.flatten(), latm.flatten()),
        mask.flatten(),
        (lon_plot.flatten(), lat_plot.flatten()),
        method="nearest",
    )
    regmask = regmask.reshape(lon_plot.shape[0], lon_plot.shape[1])
    regmask[regmask != maskval] = 0

    ###reading precipitation data
    nc = Dataset(precipfile)
    latp = nc.variables[precip_lat][:]
    lonp = nc.variables[precip_lon][:]

    Oprecip = nc.variables[precip_var][:]

    if len(lonp.shape) < 2:
        lonp, latp = np.meshgrid(lonp, latp)

    if len(Oprecip.shape) > 2:
        Oprecip = Oprecip[0, :]

    if cenlon == "0":
        lonp = np.where(lonp >= 180, lonp - 360, lonp)
    elif cenlon == "180":
        lonp = np.where(lonp < 0, lonp + 360, lonp)

    regprecip = griddata(
        (lonp.flatten(), latp.flatten()),
        Oprecip.flatten(),
        (lon_plot.flatten(), lat_plot.flatten()),
        method="nearest",
    )
    regprecip = regprecip.reshape(lon_plot.shape[0], lon_plot.shape[1])

    maskprecip = regprecip * regmask

    bias = LagP - maskprecip

    biasdq = bias * (area / density)

    dq_mask = -999 < dq

    cdq = np.empty_like(dq)
    cdq[:] = 0

    fmatrix = np.zeros((lat.shape[0] - 1, lon.shape[1] - 1))

    for k in range(lat.shape[0] - 1):
        for l in range(lon.shape[1] - 1):

            lat_mask = (lat[k, l] <= plats) & (plats < lat[k + 1, l])
            lon_mask = (lon[k, l] <= plons) & (plons < lon[k, l + 1])

            mask = dq_mask & lat_mask & lon_mask

            if np.any(mask):
                # fmatrix[k, l] = (1 / lc.g0)* np.sum((-1) * dq[mask] * pres[mask])

                dqweight = np.abs(dq[mask]) / (np.sum(np.abs(dq[mask])))

                cdq[mask] = dq[mask] * (-1) - dqweight * biasdq[k, l]
                fmatrix[k, l] = np.sum(cdq[mask])
    fmatrix[fmatrix < 0] = 0
    fmatrix = fmatrix * density / area

    cdq[cdq < 0] = 0

    return cdq, fmatrix


def processing_moisture_track_backward(
    tensor_org,
    filter_dqdt_parcels,
    dqdt_threshold,
    filter_pbl_dq_parcels,
    moist_custom_limits_highs,
    dqpblcheck,
    dqpbl_method,
    trkdq_rh_check,
    dqrh_threshold,
    mindq_gain,
    mindq_loss,
    check_RH_route_precip,
    precip_minrh_en_route,
    lag_times,
    area,
    density,
    dtime,
    moisture_linear_adjustment,
    precip_minrh,
    lon,
    lat,
    numPdY,
    numPdX,
    cenlon,
    varid,
    moisture_tracking_method,
    trackingtime_steps,
    moist_bias_correction,
    precipfile,
    precip_var,
    precip_lat,
    precip_lon,
    maskfile,
    maskvar,
    maskval,
    masklat,
    masklon,
    only_non_precip,
    rank,
    size,
    comm,
):
    # filtering parcels using t=0  and t-6#
    # tensor_moist_=np.copy(tensor_org)



    if filter_dqdt_parcels:

        tmp_matrix = tensor_org[-1, :, 3] - tensor_org[-2, :, 3]

        check_precip = is_precipitating_parcel(
            tmp_matrix,
            dqdt_threshold,
            (tensor_org[-1, :, 12] + tensor_org[-2, :, 12]) / 2,
            precip_minrh,
            only_non_precip
        )
        # check_precip=is_precipitating_parcel(tmp_matrix, dqdt_threshold,  tensor_org[-1,:,12], precip_minrh)

        if len(check_precip[check_precip == True]) == 0:
            tensor_moist_ = np.empty((tensor_org.shape[0], 1, tensor_org.shape[2]))
            tensor_moist_[:, :, :] = 0
            counter_precip_partdq = 0

        else:
            tensor_moist_ = np.empty(
                (
                    tensor_org.shape[0],
                    len(check_precip[check_precip == True]),
                    tensor_org.shape[2],
                )
            )

            for i in range(0, tensor_moist_.shape[0]):
                tensor_moist_[i, :] = tensor_org[i, :, :][check_precip == True]

            counter_precip_partdq = tensor_moist_.shape[1]
    else:
        tensor_moist_ = np.copy(tensor_org[:, :, :])
        counter_precip_partdq = None

    if filter_pbl_dq_parcels:
        check_pbl = is_withinpbl(
            tensor_moist_[-1, :, 4],
            tensor_moist_[1:, :, 7],
            moist_custom_limits_highs[0],
            moist_custom_limits_highs[1],
            "maxval",
        )

        if len(check_pbl[check_pbl == True]) == 0:
            tensor_moist = np.empty((tensor_moist_.shape[0], 1, tensor_moist_.shape[2]))
            tensor_moist[:, :, :] = 0
            counter_precip_part_pbl = 0
        else:
            tensor_moist = np.empty(
                (
                    tensor_moist_.shape[0],
                    len(check_pbl[check_pbl == True]),
                    tensor_moist_.shape[2],
                )
            )

            for i in range(0, tensor_moist.shape[0]):
                tensor_moist[i, :] = tensor_moist_[i, :, :][check_pbl == True]

            counter_precip_part_pbl = tensor_moist.shape[1]
    else:
        tensor_moist = np.copy(tensor_moist_[:, :, :])
        counter_precip_part_pbl = None

    if filter_dqdt_parcels:
        precip_matrix, dqfrack = compute_precipitation_target_region(
            lat, lon, tensor_moist, dqdt_threshold, cenlon, area, density
        )

        if moist_bias_correction:
            corrected_dq0, bias_cor_precip_matrix = correct_sink_precip(
                precip_matrix,
                tensor_moist,
                density,
                area,
                lat,
                lon,
                cenlon,
                precipfile,
                precip_var,
                precip_lat,
                precip_lon,
                maskfile,
                maskvar,
                maskval,
                masklat,
                masklon,
            )
        else:
            bias_cor_precip_matrix = np.empty((lat.shape[0] - 1, lon.shape[1] - 1))
            bias_cor_precip_matrix[:] = 0
            corrected_dq0 = np.empty_like(tensor_moist[-1, :, 3])
            corrected_dq0[:] = 0
    else:
        precip_matrix = np.empty((lat.shape[0] - 1, lon.shape[1] - 1))
        precip_matrix[:] = 0
        bias_cor_precip_matrix = np.empty((lat.shape[0] - 1, lon.shape[1] - 1))
        bias_cor_precip_matrix[:] = 0
        corrected_dq0 = np.empty_like(tensor_moist[-1, :, 3])
        corrected_dq0[:] = 0

    precipvals = tensor_moist[-1, :, 3] - tensor_moist[-2, :, 3]
    qt0 = tensor_moist[-1, :, 3]
    qt0_ = tensor_moist[-2, :, 3]
    Qt0_ = np.sum(tensor_moist[-2, :, 3])

    tensor_moist = tensor_moist[:-1, :, :]
    Qt0_ = np.sum(tensor_moist[-1, :, 3])

    matrix_moist_files = ["lon", "lat", "dq", "ds", "rh", "FH"]
    dmatrix = np.empty(
        (
            len(tensor_moist[:, 0, 0]) - 1,
            len(tensor_moist[0, :, 0]),
            len(matrix_moist_files) + 8,
        )
    )
    dmatrix[:, :, :] = 0

    if moisture_tracking_method == "testing":

        dmatrix = process_SJ05_backward(dmatrix, tensor_moist, qt0, precipvals, cenlon)

    else:

        n = tensor_moist.shape[1]
        count = n // size
        remainder_ = n % size
        remainder = 0

        if rank < remainder:
            start = rank * (count + 1)
            stop = start + count + 1
        else:
            start = rank * count + remainder
            stop = start + count

        local_parts = tensor_moist[0, start:stop, 0]
        local_results = np.empty((dmatrix.shape[0], len(local_parts), dmatrix.shape[2]))
        local_results[:, :, :] = 0

        local_results = parallel_moisture_process_backward(
            local_results,
            tensor_moist[:, start:stop, :],
            dqpblcheck,
            dqpbl_method,
            trkdq_rh_check,
            dqrh_threshold,
            mindq_gain,
            mindq_loss,
            moisture_linear_adjustment,
            precip_minrh,
            qt0[start:stop],
            qt0[start:stop],
            precipvals[start:stop],
            trackingtime_steps,
            cenlon,
            corrected_dq0[start:stop],
            moist_bias_correction,
            check_RH_route_precip,
            precip_minrh_en_route,
        )

        if rank > 0:
            comm.Send(local_results, dest=0, tag=14)

        else:
            i_start = []
            i_stop = []
            for i in range(size):
                if i < remainder:
                    i_start = np.append(i_start, i * (count + 1))
                    i_stop = np.append(i_stop, i_start + count + 1)
                else:
                    ii_start = i * count + remainder
                    ii_stop = ii_start + count
                    i_start = np.append(i_start, ii_start)
                    i_stop = np.append(i_stop, ii_stop)

            final_results = np.copy(local_results)

            dmatrix[:, int(i_start[0]) : int(i_stop[0]), :] = final_results

            for i in range(1, size):

                if i < remainder:
                    rank_size = count + 1
                else:
                    rank_size = count

                tmp = np.empty(
                    (final_results.shape[0], rank_size, final_results.shape[2]),
                    dtype=np.float64,
                )

                comm.Recv(tmp, source=i, tag=14)

                dmatrix[:, int(i_start[i]) : int(i_stop[i]), :] = tmp

        comm.Bcast(dmatrix, root=0)

        if remainder_ > 0:
            dmatrix[:, count * size :] = parallel_moisture_process_backward(
                dmatrix[:, count * size :],
                tensor_moist[:, count * size :, :],
                dqpblcheck,
                dqpbl_method,
                trkdq_rh_check,
                dqrh_threshold,
                mindq_gain,
                mindq_loss,
                moisture_linear_adjustment,
                precip_minrh,
                qt0[count * size :],
                qt0[count * size :],
                precipvals[count * size :],
                trackingtime_steps,
                cenlon,
                corrected_dq0[count * size :],
                moist_bias_correction,
                check_RH_route_precip,
                precip_minrh_en_route,
            )

    matrix_moist = np.copy(dmatrix[:, :, :6])

    matrix_moist[:, :, 2][matrix_moist[:, :, 2] == -999.9] = 0

    matrix_moist_corrected = np.copy(dmatrix[:, :, :6])
    matrix_moist_corrected[:, :, 2] = dmatrix[:, :, 11]
    matrix_moist_corrected[:, :, 2][matrix_moist_corrected[:, :, 2] == -999.9] = 0
    matrix_moist_corrected[:, :, 2][matrix_moist_corrected[:, :, 2] < 0] = 0

    sum_prec = dmatrix[0, :, 6]
    part_uptakes = dmatrix[0, :, 7]
    partt = dmatrix[0, :, 8]
    noprecip = dmatrix[0, :, 9]

    lwvrt_ = dmatrix[0, :, 10]

    part_uptakes_bias = dmatrix[0, :, 12]

    if moisture_linear_adjustment:

        aux_lwvrt = lwvrt_[lwvrt_ > 0]

        lmrt = aux_lwvrt / (1440)

        lwvrt = np.mean(lmrt[np.isfinite(lmrt)])
    else:
        lwvrt = "not computed"

    array_day, moistd = compute_var_integarated_day_moist(
        matrix_moist,
        lag_times,
        area,
        density,
        dtime,
        lon,
        lat,
        numPdY,
        numPdX,
        cenlon,
        varid,
    )

    array_day_corrected, moistd_corrected = compute_var_integarated_day_moist(
        matrix_moist_corrected,
        lag_times,
        area,
        density,
        dtime,
        lon,
        lat,
        numPdY,
        numPdX,
        cenlon,
        varid,
    )

    array_day_corrected[array_day_corrected < 0] = 0

    if moist_bias_correction:

        totalmu = np.sum(array_day_corrected, axis=0)

        bias = np.sum(totalmu) - np.sum(bias_cor_precip_matrix)

        weights = totalmu / np.sum(totalmu)

        bias_weights = bias * weights

        corrected_mu = totalmu - bias_weights

        daysum = []
        newdaysum = []

        adjusted_moist_day = np.empty_like(array_day_corrected)
        adjusted_moist_day[:, :, :] = 0
        for iday in range(0, array_day_corrected.shape[0]):
            sumday = np.sum(array_day_corrected[iday, :])

            adjust_sumday = (sumday / np.sum(totalmu)) * np.sum(corrected_mu)

            daysum = np.append(daysum, sumday)
            newdaysum = np.append(newdaysum, adjust_sumday)

            dayweights = array_day_corrected[iday, :] / sumday

            adjusted_moist_day[iday, :] = (
                array_day_corrected[iday, :] - (sumday - adjust_sumday) * dayweights
            )

        lwvrt_bias_ = dmatrix[0, :, 13]

        aux_lwvrt_bias = lwvrt_bias_[lwvrt_ > 0]

        lmrt_bias = aux_lwvrt_bias / (1440)

        lwvrt_bias = np.mean(lmrt_bias[np.isfinite(lmrt_bias)])

    else:
        adjusted_moist_day = None
        lwvrt_bias = None

    counter_precip_part = matrix_moist.shape[1]

    no_evap_uptakes = matrix_moist.shape[1] - int(np.sum(part_uptakes))
    no_evap_uptakes_bias = matrix_moist.shape[1] - int(np.sum(part_uptakes_bias))
    if moisture_linear_adjustment:
        CR = (moistd / Qt0_) * 100

    else:
        CR = None

    sum_prec_ = np.sum(sum_prec)
    if sum_prec_ == 0:
        sum_prec_ = 1

    return (
        array_day,
        matrix_moist,
        counter_precip_part,
        no_evap_uptakes,
        np.sum(partt) / sum_prec_,
        CR,
        lwvrt,
        precip_matrix,
        bias_cor_precip_matrix,
        adjusted_moist_day,
        moistd_corrected,
        no_evap_uptakes_bias,
        lwvrt_bias,
    )


def compute_linear_discounted(var, valid):

    var[var == -999.9] = 0

    result_var = np.empty_like(var)
    result_var[:] = 0

    for i in range(0, len(var)):
        if var[i] > 0:
            result_var[i] = var[i]

        else:
            suma = np.sum(result_var[:i])

            for j in range(0, i):
                if suma > 0:
                    aux = result_var[j] - ((result_var[j] / suma) * abs(var[i]))
                else:
                    aux = 0

                if aux <= 0:
                    result_var[j] = 0
                else:
                    result_var[j] = aux

    return result_var


def compute_mean_lon(lons, cenlon):
    lons1 = np.array(lons[1:])
    lons2 = np.array(lons[:-1])
    mean = np.empty(len(lons1))

    # Normalize longitudes to the range [-180, 180]
    lons1 = np.where(lons1 >= 180, lons1 - 360, lons1)
    lons2 = np.where(lons2 >= 180, lons2 - 360, lons2)

    for i in range(len(lons1)):
        if lons1[i] >= -180 and lons2[i] >= -180:
            check = (np.abs(lons1[i]) + np.abs(lons2[i])) / 2

            if lons1[i] > 0 and lons2[i] < 0 and np.abs(check - 180) < 90:
                lons2[i] += 360
            elif lons1[i] < 0 and lons2[i] > 0 and np.abs(check - 180) < 90:
                lons1[i] += 360

            val = (lons1[i] + lons2[i]) / 2
            val = val - 360 if val >= 180 else val

            mean[i] = val
        elif lons1[i] < -180 <= lons2[i] <= 180:
            mean[i] = lons2[i]
        elif lons2[i] < -180 <= lons1[i] <= 180:
            mean[i] = lons1[i]
        else:
            mean[i] = -999.9

    return mean


def midpoint(lat1, lon1, lat2, lon2):
    # from following http://www.movable-type.co.uk/scripts/latlong.html

    # Convert latitude and longitude from degrees to radians
    lat1 = math.radians(lat1)
    lon1 = math.radians(lon1)
    lat2 = math.radians(lat2)
    lon2 = math.radians(lon2)

    # Calculate the midpoint
    dlon = lon2 - lon1

    Bx = math.cos(lat2) * math.cos(dlon)
    By = math.cos(lat2) * math.sin(dlon)

    lat3 = math.atan2(
        math.sin(lat1) + math.sin(lat2), math.sqrt((math.cos(lat1) + Bx) ** 2 + By**2)
    )
    lon3 = lon1 + math.atan2(By, math.cos(lat1) + Bx)

    # Convert the midpoint back to degrees
    lat3 = math.degrees(lat3)
    lon3 = math.degrees(lon3)

    return lat3, lon3


def compute_mean_lon_v2(lons, lats, cenlon):
    # Input values as degrees
    lon1 = np.array(lons[1:])
    lon2 = np.array(lons[:-1])

    lat1 = np.array(lats[1:])
    lat2 = np.array(lats[:-1])

    mlats = []
    mlons = []
    for la1, lo1, la2, lo2 in zip(lat1, lon1, lat2, lon2):
        mid_lat, mid_lon = midpoint(la1, lo1, la2, lo2)

        mlats = np.append(mlats, mid_lat)
        mlons = np.append(mlons, mid_lon)

    if cenlon == "180":
        mlons = (mlons + 360) % 360
    else:
        mlons = (mlons + 180) % 360 - 180

    return mlons, mlats


def process_SJ05_backward(dmatrix, tensor_moist, qt0, precipvals, cenlon):

    tensor_moist[tensor_moist < -500] = np.nan

    for i in range(0, tensor_moist.shape[1], cenlon):
        dlon = compute_mean_lon(tensor_moist[:, i, 1])
        dmatrix[:, i, 0] = dlon

    dmatrix[:, :, 1] = (tensor_moist[1:, :, 2] + tensor_moist[:-1, :, 2]) / 2
    dmatrix[:, :, 2] = tensor_moist[1:, :, 3] - tensor_moist[:-1, :, 3]
    dmatrix[:, :, 3] = tensor_moist[1:, :, 11] - tensor_moist[:-1, :, 11]
    dmatrix[:, :, 4] = tensor_moist[1:, :, 12] - tensor_moist[:-1, :, 12]

    dmatrix[:, :, 5] = 1
    dmatrix[:, :, 6] = 1
    dmatrix[:, :, 7] = 1
    dmatrix[:, :, 8] = 1
    dmatrix[:, :, 9] = 1
    dmatrix[:, :, 10] = 0

    return dmatrix


def compute_lag_mwvrt(mtime, fc, qt0):

    fc[fc == -999.9] = 0

    if np.sum(fc) > 0:
        wvrt = np.sum(mtime * fc / np.sum(fc))
    else:
        wvrt = -999.9
    return np.abs(wvrt)


def parallel_moisture_process_backward(
    dmatrix,
    tensor_moist,
    dqpblcheck,
    dqpbl_method,
    trkdq_rh_check,
    dqrh_threshold,
    mindq_gain,
    mindq_loss,
    moisture_linear_adjustment,
    precip_minrh,
    qt0,
    Qt0_,
    precipvals,
    trackingtime_steps,
    cenlon,
    corrected_dq0,
    moist_bias_correction,
    check_RH_route_precip,
    precip_minrh_en_route,
):
    # uptakes_parts_precip=0
    # sum_prec=0
    # partt=0
    # noprecip=0

    trackingtime_steps = trackingtime_steps * (-1)
    trackingtime_steps = trackingtime_steps[::-1]
    mtime = (trackingtime_steps[1:] + trackingtime_steps[:-1]) / 2

    for i in range(0, tensor_moist.shape[1]):

        lons = tensor_moist[:, i, 1]

        dlon = compute_mean_lon(lons, cenlon)

        if cenlon == "180":
            dlon = (dlon + 360) % 360

            dlon[dlon < 0] = -999.9

        lats = tensor_moist[:, i, 2]

        # dlon, dlat = compute_mean_lon_v2(lons, lats, cenlon)

        dlat = cal_track_diff(lats, "mean")

        qs = tensor_moist[:, i, 3]
        dq = cal_track_diff(qs, "diff")

        # sum_prec= sum_prec + np.abs(precipvals[i])
        dmatrix[:, i, 6] = np.abs(precipvals[i])

        var = tensor_moist[:, i, 11]
        dvar = cal_track_diff(var, "diff")

        rh = tensor_moist[:, i, 12]
        prh = tensor_moist[1:, i, 12]
        drh = cal_track_diff(rh, "diff")

        l1 = (dq > mindq_loss) & (dq < mindq_gain)

        dq[l1] = 0

        if check_RH_route_precip:

            dq_mask = (dq < 0) & (prh < precip_minrh_en_route)
            dq[dq_mask] = 0

        dmatrix[:, i, 0] = dlon
        dmatrix[:, i, 1] = dlat
        dmatrix[:, i, 2] = dq
        dmatrix[:, i, 3] = dvar
        dmatrix[:, i, 4] = drh
        dmatrix[:, i, 5] = 0

        # trk_check_pbl=is_pbl_check(dqpblcheck, tensor_moist[:,i,4], tensor_moist[:,i, 7], len(dlon))

        trk_check_pbl = is_in_pbl(
            dqpblcheck,
            tensor_moist[:, i, 4],
            tensor_moist[:, i, 7],
            dqpbl_method,
            len(dlon),
        )

        check_dvar = np.ones((len(dlon)), dtype=bool)
        check_dvar[dq >= mindq_gain] = True
        check_dvar[dq < mindq_gain] = False

        trk_drh = np.ones((len(dlon)), dtype=bool)
        if trkdq_rh_check:

            trk_drh[np.abs(drh) <= dqrh_threshold] = True
            trk_drh[np.abs(drh) > dqrh_threshold] = False

        valid = trk_check_pbl & trk_drh & check_dvar

        dmatrix[:, i, 5][valid == True] = 1
        dmatrix[:, i, 5][valid == False] = 0

        dmatrix[np.isnan(dmatrix)] = -999.9

        if moisture_linear_adjustment:

            dmatrix[:, i, 2] = compute_linear_discounted(dmatrix[:, i, 2], valid)

            if corrected_dq0[i] <= 0:
                dmatrix[:, i, 11] = 0

            else:

                dqi = dmatrix[:, i, 2]

                new_dqi = np.empty_like(dqi)
                new_dqi[:] = 0

                dq_mask = dqi > 0

                sumdq = np.sum(dqi[dq_mask])

                dqfrac = dqi / sumdq

                biasdq = sumdq - (corrected_dq0[i])

                new_dqi[dq_mask] = dqi[dq_mask] - dqfrac[dq_mask] * biasdq

                dmatrix[:, i, 11] = new_dqi

                dmatrix[:, i, 13] = compute_lag_mwvrt(
                    mtime, dmatrix[:, i, 11] * dmatrix[:, i, 5], qt0[i]
                )

            # print(np.sum(dmatrix[:,i,2]), np.sum(dmatrix[:,i,11]), corrected_dq0[i])

            # adj_fac = dmatrix[:,i,2] / qs[-1]

            if np.sum((dmatrix[:, i, 2])) > 0:
                dmatrix[:, i, 7] = 1

            if np.sum((dmatrix[:, i, 11])) > 0:
                dmatrix[:, i, 12] = 1

            dmatrix[:, i, 2] = dmatrix[:, i, 2]
            dmatrix[:, i, 8] = np.sum(dmatrix[:, i, 2])

            dmatrix[:, i, 10] = compute_lag_mwvrt(
                mtime, dmatrix[:, i, 2] * dmatrix[:, i, 5], qt0[i]
            )
        else:
            if len(valid[valid == True]) >= 1:
                # uptakes_parts_precip = uptakes_parts_precip + 1
                dmatrix[:, i, 7] = 1

        # noprecip+=qt0[i]
        dmatrix[:, i, 9] = qt0[i]
    return dmatrix
