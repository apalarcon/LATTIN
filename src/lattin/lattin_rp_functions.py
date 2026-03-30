import numpy as np
from netCDF4 import Dataset, num2date, date2num
import os
import datetime
from datetime import datetime, timedelta
import scipy.stats as scps
from scipy import ndimage
import pandas as pd
from mpi4py import MPI
import time
import gc
from threading import Thread, Event
from lattin.constants import *
from lattin.lattin_functions import (
    compute_mean_lon,
    is_withinpbl,
    compute_theta,
    str2boolean,
    calc_pottemp,
    get_date_format,
)
from scipy.interpolate import RegularGridInterpolator
from lattin.fmodules import (
    compute_grid_integrated_temperature_anom as compute_grid_integrated_temperature_anom,
)
from lattin.fmodules import (
    compute_temperature_anom_genesis_properties as compute_temperature_anom_genesis_properties,
)
from lattin.fmodules import (
    compute_temperature_anomalies_sources as compute_temperature_anomalies_sources,
)


def deg2rad(inval):
    """
    helper function to convert lat/lon in degrees to radians
    """
    outval = inval / 360 * 2 * np.pi
    return outval


def gc_dist(lat1, lon1, lat2, lon2, r):
    """
    lat1, lon1 locate the grid points from which distances should be computed (lat/lon both in degrees)
    lat2, lon2 define arrays of grid points for which to compute the distances
    r denotes the radius of the earth in meters (6371000 is usually appropriate)
    """
    lat1_rad = deg2rad(lat1)
    lon1_rad = deg2rad(lon1)
    lat2_rad = deg2rad(lat2)
    lon2_rad = deg2rad(lon2)

    dist = r * np.arccos(
        np.sin(lat1_rad) * np.sin(lat2_rad)
        + np.cos(lat1_rad) * np.cos(lat2_rad) * np.cos(lon2_rad - lon1_rad)
    )
    dist[np.isnan(dist)] = 0
    return dist


def Tanom_tracking_parms(
    Tanom_tracking_method,
    save_Tanom_parts_position,
    path_clim_temperature,
    climT_fname_prefix,
    climT_date_format,
    Tlat_var_name,
    Tlon_var_name,
    climTvar_name,
    dTdp_var_name,
    dTdt_var_name,
    Tplves_var_name,
    analysis_levels,
    dp_sfc,
    dp_upper,
    Tanom_linear_adjustment,
    psfc_var_name,
    Tanom_threshold,
    interpolate_parcel_temperature,
    interpolate_sfc,
    Tvar_name,
    path_to_meteodata,
    meteodata_fname_prefix,
    meteodata_date_format,
    meteolat_var_name,
    meteolon_var_name,
    meteo_plves_var_name
):
    """
    This function checks the parameters of the TEMPERATURE ANOMALY TRACKING module

    Parameters
    ----------
    Tanom_tracking_method : str
        Temperature anomaly tracking method: RP23 or PR23
    save_Tanom_parts_position : bool
        To save processed parcels trajectories (temperature anomaly tracking)
    path_clim_temperature : str
        Path to netCDF files containing climatological temperature and other variables if needed
    climT_fname_prefix : str
        Prefix of the name of the climatological netCDf files
    climT_date_format : str
        Format of the date in the name of the climatological netCDf files
    Tlat_var_name : str
        Latitude variable name in the climatological temperature files
    Tlon_var_name : str
        Longitude variable name in the climatological temperature files
    climTvar_name : str
        Climatological air temperature variable name
    dTdp_var_name : str
        dTdp variable name
    dTdt_var_name : str
        dTdt variable name
    Tplves_var_name : str
        Pressure levels variable name
    analysis_levels : list
        Atmospheric levels to perform the analysis
    dp_sfc : float
        To filter air parcels close to surface [hPa]
    dp_upper : float
        To filter air parcels at specific upper level [hPa]
    Tanom_linear_adjustment : bool
        Apply linear adjusment to temperature changes along parcel trajectories
    psfc_var_name : str
        Surface pressure variable name
    Tanom_threshold : list
        Temperature anomaly threshold for filtering air parcels trajectories [K]
    interpolate_parcel_temperature : bool
        To interpotate air temperature to parcel trajectories
    interpolate_sfc : bool
        To interpotate surface pressure to parcel trajectories
    Tvar_name : str
        Air temperature variable name

    Returns
    -------
    save_Tanom_parts_position : bool
        To save processed parcels trajectories (temperature anomaly tracking)
    analysis_levels : list
        Atmospheric levels to perform the analysis
    dp_sfc : float
        To filter air parcels close to surface [hPa]
    dp_upper : float
        To filter air parcels at specific upper level [hPa]
    Tanom_linear_adjustment : bool
        Apply linear adjusment to temperature changes along parcel trajectories
    Tanom_threshold : list
        Temperature anomaly threshold for filtering air parcels trajectories [K]
    climT_date_format : str
        Format of the date in the name of the climatological netCDf files
    interpolate_parcel_temperature : bool
        To interpotate air temperature to parcel trajectories
    interpolate_sfc : bool
        To interpotate surface pressure to parcel trajectories
    errors : str
        Errors found in the parameters
    errors_found : bool
        If errors were found
    """
    errors_found = False
    errors = ""

    Tanom_tracking_method = Tanom_tracking_method
    save_Tanom_parts_position = save_Tanom_parts_position
    path_clim_temperature = path_clim_temperature
    climT_fname_prefix = climT_fname_prefix
    climT_date_format = climT_date_format
    Tlat_var_name = Tlat_var_name
    Tlon_var_name = Tlon_var_name
    climTvar_name = climTvar_name
    dTdp_var_name = dTdp_var_name
    dTdt_var_name = dTdt_var_name
    Tplves_var_name = Tplves_var_name
    psfc_var_name = psfc_var_name
    analysis_levels = analysis_levels
    dp_sfc = dp_sfc
    dp_upper = dp_upper
    Tanom_linear_adjustment = Tanom_linear_adjustment
    Tanom_threshold = Tanom_threshold
    path_to_meteodata=path_to_meteodata
    meteodata_fname_prefix=meteodata_fname_prefix
    meteodata_date_format=meteodata_date_format
    meteolat_var_name=meteolat_var_name
    meteolon_var_name=meteolon_var_name
    meteo_plves_var_name=meteo_plves_var_name

    if not analysis_levels:
        analysis_levels = ["pbl", 500]

    if dp_sfc == "" or dp_sfc == None:
        dp_sfc = 50

    if dp_upper == "" or dp_upper == None:
        dp_upper = 50

    if not Tanom_threshold:
        Tanom_threshold = [0, 0]

    if save_Tanom_parts_position == "" or save_Tanom_parts_position == None:
        save_Tanom_parts_position = True

    if Tanom_tracking_method == "RP23":
        Tanom_linear_adjustment = False
    if Tanom_tracking_method == "PR23":
        Tanom_linear_adjustment = True

    if not Tanom_tracking_method in ["RP23", "PR23"]:
        errors = errors + "(**) ERROR: Tanom_tracking_method must be RP23 or PR23\n"
        errors_found = True

    Tanom_message = " parameter is missing in the TEMPERATURE ANOMALY TRACKING module of the input file\n"
    if interpolate_parcel_temperature == None or interpolate_parcel_temperature == "":
        interpolate_parcel_temperature = False
    else:
        interpolate_parcel_temperature = str2boolean(interpolate_parcel_temperature)



    if interpolate_sfc == None or interpolate_sfc == "":
        interpolate_sfc = False
    else:
        interpolate_sfc = str2boolean(interpolate_sfc)



    Tanom_parms = [
        "path_clim_temperature",
        "climT_fname_prefix",
        "climT_date_format",
        "Tlat_var_name",
        "Tlon_var_name",
        "climTvar_name",
        "dTdp_var_name",
        "dTdt_var_name",
        "Tplves_var_name",
        "Tanom_linear_adjustment",
    ]

    for param in Tanom_parms:
        if vars()[param] == "" or vars()[param] == None:
            errors = errors + "(**) ERROR: " + param + Tanom_message
            errors_found = True

    if not climT_date_format == "" or not climT_date_format == None:
        climT_date_format, etemp, emessage = get_date_format(
            climT_date_format, "climT_date_format"
        )
        if etemp:
            errors = errors + "(**) ERROR: " + emessage
            errors_found = True

    if not analysis_levels:
        errors = (
            errors
            + "(**) ERROR: analysis_levels parameter is missing or is empty in the TEMPERATURE ANOMALY TRACKING module of the input file\n"
        )
        errors_found = True
    Tlparams = ["dp_sfc", "dp_upper"]
    for param in Tlparams:
        if vars()[param] == "" or vars()[param] == None:
            if any(item != "pbl" for item in analysis_levels):
                errors = errors + "(**) ERROR: " + param + Tanom_message
                errors_found = True

    if not Tanom_threshold:
        errors = (
            errors
            + "(**) ERROR: Tanom_threshold parameter is missing or is empty in the TEMPERATURE ANOMALY TRACKING module of the input file\n"
        )
        errors_found = True

    if len(Tanom_threshold) != len(analysis_levels):
        errors = (
            errors
            + "(**) ERROR: Tanom_threshold and analysis_levels parameters must have the same number of elements\n"
        )
        errors_found = True

    if save_Tanom_parts_position == "" or save_Tanom_parts_position == None:
        save_Tanom_parts_position = True

    if any(item == "sfc" for item in analysis_levels):
        if not interpolate_sfc:
            errors = (
                errors
                + "(**) ERROR: 'sfc' is in 'analysis_levels' but 'interpolate_sfc=False'. Please remove 'sfc' from 'analysis_levels' or set 'interpolate_sfc=True'\n"
            )
            errors_found = True

    if interpolate_sfc or interpolate_parcel_temperature:

        for param in ['path_to_meteodata','meteodata_fname_prefix','meteolat_var_name','meteolon_var_name']:
            if vars()[param] == "" or vars()[param] == None:
                errors = errors + "(**) ERROR: " + param + Tanom_message
                errors_found = True

        if interpolate_sfc:
            if psfc_var_name == "" or psfc_var_name == None:
                errors = errors + "(**) ERROR: psfc_var_name " + Tanom_message
                errors_found = True

        if  interpolate_parcel_temperature:
            if Tvar_name == "" or Tvar_name == None:
                errors = errors + "(**) ERROR: Tvar_name" + Tanom_message
                errors_found = True

            if meteo_plves_var_name == "" or meteo_plves_var_name == None:
                errors = errors + "(**) ERROR: meteo_plves_var_name" + Tanom_message
                errors_found = True



        if not meteodata_date_format == "" or not meteodata_date_format == None:
            meteodata_date_format, epsfc, emessage = get_date_format(
                meteodata_date_format, "meteodata_date_format"
            )
            if epsfc:
                errors = errors + "(**) ERROR: " + emessage
                errors_found = True


    return (
        save_Tanom_parts_position,
        analysis_levels,
        dp_sfc,
        dp_upper,
        Tanom_linear_adjustment,
        Tanom_threshold,
        climT_date_format,
        interpolate_parcel_temperature,
        interpolate_sfc,
        meteodata_date_format,
        errors,
        errors_found,
    )


def cal_track_diff_forward(var_vals, case):
    """
    Compute differences, means, or maxima of adjacent elements in an array,
    depending on the value of the string argument "case". This function is
    used to compute differences, means, or maxima of parcel properties along
    trajectories. The input array is first masked to exclude values below -500
    (presumably missing data). The function then takes the difference, mean,
    or maximum of adjacent elements in the array, depending on the value of
    "case". The output is an array of the same length as the input, but with
    the first element removed (since the first element will not have a
    difference, mean, or maximum with the "next" element, which does not exist).

    Parameters
    ----------
    var_vals : ndarray
        Input array of parcel properties
    case : str
        String indicating whether to compute differences, means, or maxima.
        Options are "diff", "mean", and "max".

    Returns
    -------
    dvalue : ndarray
        Array of differences, means, or maxima of adjacent elements in the input
        array.
    """
    var_vals[var_vals < -500] = np.nan
    # difference
    if case in ["diff"]:
        dvalue = var_vals[:-1] - var_vals[1:]
    # mean
    if case in ["mean"]:
        dvalue = (var_vals[:-1] + var_vals[1:]) / 2
    # max
    if case in ["max"]:
        dvalue = np.amax(var_vals[:-1], var_vals[1:])
    return dvalue


def interpolate_temperature_field_regulargrid(
    parcel_positions,
    Tlatitude,
    Tlongitude,
    Tpressure_levels,
    climTemperature_field,
    dT_dp,
    dT_dt,
    surface_pressure,
    meteo_latitude,
    meteo_longitude,
    meteo_plves,
    interpolate_parcel_temperature,
    interpolate_sfc,
    temperature_field,
    interpo_method="linear",
    chunk_size=250000,
):
    """
    Interpolates the temperature field to the positions of air parcels.

    Parameters:
        parcel_positions (ndarray): Array of shape (N, 3) with [lat, lon, pressure] of air parcels.
        Tlatitude (ndarray): 2D Latitude array matching the spatial grid of Temperature_field.
        Tlongitude (ndarray): 2D Longitude array matching the spatial grid of Temperature_field.
        Tpressure_levels (ndarray): 1D array of pressure levels.
        Temperature_field (ndarray): 3D temperature field of shape (pressure, lat, lon).
        interpo_method (str): Interpolation method ('linear', 'nearest', 'cubic').

    Returns:
        tuple: Interpolated temperatures, dT/dp, and dT/dt at the positions of air parcels.
    """
    # Ensure latitude, longitude, and pressure levels are sorted in ascending order
    lat_sort_idx = np.argsort(Tlatitude)
    Tlatitude = Tlatitude[lat_sort_idx]

    lon_sort_idx = np.argsort(Tlongitude)
    Tlongitude = Tlongitude[lon_sort_idx]

    meteo_lat_sort_idx = np.argsort(meteo_latitude)
    meteo_latitude = meteo_latitude[meteo_lat_sort_idx]

    meteo_lon_sort_idx = np.argsort(meteo_longitude)
    meteo_longitude = meteo_longitude[meteo_lon_sort_idx]


    pres_sort_idx = np.argsort(Tpressure_levels)
    Tpressure_levels = Tpressure_levels[pres_sort_idx]

    if interpolate_parcel_temperature:
        meteo_plves_sort_idx = np.argsort(meteo_plves)
        meteo_plves = meteo_plves[meteo_plves_sort_idx]


    # --- 2. Sort the large 3D data fields axis by axis to avoid copies ---
    # This is the most critical change to prevent OOM errors.
    # print("      -> Sorting large data fields in a memory-safe way...")

    # List of all large 3D fields to be sorted
    fields_to_sort = [climTemperature_field, dT_dp, dT_dt]
    #if interpolate_parcel_temperature:
    #    fields_to_sort.append(temperature_field)

    for field in fields_to_sort:
        # Sort along the last axis (longitude)
        field[:] = field[:, :, lon_sort_idx]
        # Sort along the middle axis (latitude)
        field[:] = field[:, lat_sort_idx, :]
        # Sort along the first axis (pressure)
        field[:] = field[pres_sort_idx, :, :]
        # The [:] is crucial - it forces an in-place modification.

    # Sort the smaller 2D surface pressure field
    if interpolate_sfc:
        surface_pressure[:] = surface_pressure[:, meteo_lon_sort_idx]
        surface_pressure[:] = surface_pressure[meteo_lat_sort_idx, :]

    # print("      -> Sorting complete.")
    gc.collect()

    # --- 3. Create Interpolators (low memory usage) ---
    # print("      -> Creating interpolators...")
    climtemperature_interpolator = RegularGridInterpolator(
        (Tpressure_levels, Tlatitude, Tlongitude),
        climTemperature_field,
        method=interpo_method,
        bounds_error=False,
        fill_value=np.nan,
    )
    dTdp_interpolator = RegularGridInterpolator(
        (Tpressure_levels, Tlatitude, Tlongitude),
        dT_dp,
        method=interpo_method,
        bounds_error=False,
        fill_value=np.nan,
    )
    dTdt_interpolator = RegularGridInterpolator(
        (Tpressure_levels, Tlatitude, Tlongitude),
        dT_dt,
        method=interpo_method,
        bounds_error=False,
        fill_value=np.nan,
    )

    # We can now delete the original fields as they are copied inside the interpolator
    del climTemperature_field, dT_dp, dT_dt
    gc.collect()

    if interpolate_sfc:
        psfc_interpolator = RegularGridInterpolator(
            (meteo_latitude, meteo_longitude),
            surface_pressure,
            method=interpo_method,
            bounds_error=False,
            fill_value=np.nan,
        )
        del surface_pressure
        gc.collect()

    if interpolate_parcel_temperature:
        # Sort along the last axis (longitude)
        temperature_field[:] = temperature_field[:, :, meteo_lon_sort_idx]
        # Sort along the middle axis (latitude)
        temperature_field[:] = temperature_field[:, meteo_lat_sort_idx, :]
        # Sort along the first axis (pressure)
        temperature_field[:] = temperature_field[meteo_plves_sort_idx, :, :]


        temperature_interpolator = RegularGridInterpolator(
            (Tpressure_levels, meteo_latitude, meteo_longitude),
            temperature_field,
            method=interpo_method,
            bounds_error=False,
            fill_value=np.nan,
        )
        del temperature_field
        gc.collect()

    # --- 4. Chunked Interpolation (Already memory-efficient) ---
    num_points = parcel_positions.shape[0]
    final_climtemperatures = np.empty(num_points, dtype=np.float64)
    final_dTdp = np.empty(num_points, dtype=np.float64)
    final_dTdt = np.empty(num_points, dtype=np.float64)
    final_psfc = np.full(num_points, -999.9, dtype=np.float64)
    final_temperatures = np.full(num_points, -999.9, dtype=np.float64)

    # print(f"      -> Processing {num_points} points in chunks of {chunk_size}...")

    for start_idx in range(0, num_points, chunk_size):
        end_idx = min(start_idx + chunk_size, num_points)
        positions_chunk = parcel_positions[start_idx:end_idx]
        coords_chunk_3d = positions_chunk[:, [2, 0, 1]]
        coords_chunk_2d = positions_chunk[:, [0, 1]]

        final_climtemperatures[start_idx:end_idx] = climtemperature_interpolator(
            coords_chunk_3d
        )
        final_dTdp[start_idx:end_idx] = dTdp_interpolator(coords_chunk_3d)
        final_dTdt[start_idx:end_idx] = dTdt_interpolator(coords_chunk_3d)

        if interpolate_sfc:
            final_psfc[start_idx:end_idx] = psfc_interpolator(coords_chunk_2d)
        if interpolate_parcel_temperature:
            final_temperatures[start_idx:end_idx] = temperature_interpolator(
                coords_chunk_3d
            )

    # print("      -> Chunked interpolation finished.")
    # quit()
    return (
        final_climtemperatures,
        final_dTdp,
        final_dTdt,
        final_psfc,
        final_temperatures,
    )


def generate_climT_filenames(climT_fname_prefix, climT_date_format, track_times):
    """
    Generates a list of filenames for the climate temperature files.

    Parameters:
        climT_fname_prefix (str): Prefix of the filename.
        climT_date_format (str): Date format string for strftime.
        track_times (list): List of datetime objects.

    Returns:
        list: List of filenames for the climate temperature files.
    """
    fnames = []
    for itime in track_times:
        fnames.append(
            f"{climT_fname_prefix}" + itime.strftime(climT_date_format) + ".nc"
        )
    return fnames



def  generate_meteo_filenames(meteodata_fname_prefix, meteodata_date_format, track_times):
    """
    Generates a list of filenames for the climate temperature files.

    Parameters:
        meteodata_fname_prefix (str): Prefix of the filename.
        psfc_date_format (str): Date format string for strftime.
        track_times (list): List of datetime objects.

    Returns:
        list: List of filenames for the climate temperature files.
    """
    fnames = []
    for itime in track_times:
        fnames.append(
            f"{meteodata_fname_prefix}" + itime.strftime(meteodata_date_format) + ".nc"
        )
    return fnames


def checking_Tclimatological_files(
    path_clim_temperature,
    climT_filenames,
    Tlat_var_name,
    Tlon_var_name,
    climTvar_name,
    dTdp_var_name,
    dTdt_var_name,
    psfc_var_name,
    Tvar_name,
    interpolate_parcel_temperature,
    interpolate_sfc,
):

    """
    Checks the presence of specified variables in climatological NetCDF files.

    Parameters:
        path_clim_temperature (str): Directory path containing the climate temperature files.
        climT_filenames (list): List of climate temperature NetCDF filenames.
        Tlat_var_name (str): Latitude variable name in the files.
        Tlon_var_name (str): Longitude variable name in the files.
        climTvar_name (str): Climatological temperature variable name in the files.
        dTdp_var_name (str): Temperature gradient with respect to pressure variable name in the files.
        dTdt_var_name (str): Temperature gradient with respect to time variable name in the files.
        psfc_var_name (str): Surface pressure variable name, included if interpolate_sfc is True.
        Tvar_name (str): Air temperature variable name, included if interpolate_parcel_temperature is True.
        interpolate_parcel_temperature (bool): Indicates if air temperature should be interpolated.
        interpolate_sfc (bool): Indicates if surface pressure should be interpolated.

    Returns:
        tuple: A string containing error messages (if any), and a boolean indicating if any errors were found.
    """

    variable_names = [
        Tlat_var_name,
        Tlon_var_name,
        climTvar_name,
        dTdp_var_name,
        dTdt_var_name,
    ]

    if interpolate_sfc:
        variable_names = np.append(variable_names, psfc_var_name)

    if interpolate_parcel_temperature:
        variable_names = np.append(variable_names, Tvar_name)

    errors = ""
    find_error = False
    for Tfilename in climT_filenames:
        try:
            nc = Dataset(f"{path_clim_temperature}/{Tfilename}")
            for varname in variable_names:
                if not varname in nc.variables.keys():
                    errors = (
                        errors
                        + f"     (**) ERROR: variable {varname} is not in file: "
                        + Tfilename
                        + "\n"
                    )
                    find_error = True

            nc.close()
        except:
            errors = (
                errors
                + "     (**) ERROR: No such file or directory: "
                + Tfilename
                + "\n"
            )
            find_error = True

    return errors, find_error




def checking_meteo_files(
    path_to_meteodata,
    meteo_filenames,
    meteolat_var_name,
    meteolon_var_name,
    psfc_var_name,
    Tvar_name,
    interpolate_sfc,
    interpolate_parcel_temperature,
    meteo_plves_var_name
):

    """
    Checks the presence of specified variables in climatological NetCDF files.

    Parameters:
        path_psfc (str): Directory path containing the files.
        climT_filenames (list): List of climate temperature NetCDF filenames.
        plat_var_name (str): Latitude variable name in the files.
        plon_var_name (str): Longitude variable name in the files.
        psfc_var_name (str): Surface pressure variable name, included if interpolate_sfc is True.
        interpolate_sfc (bool): Indicates if surface pressure should be interpolated.

    Returns:
        tuple: A string containing error messages (if any), and a boolean indicating if any errors were found.
    """

    variable_names = [
        meteolat_var_name,
        meteolon_var_name,
    ]

    if interpolate_sfc:
        variable_names = np.append(variable_names, psfc_var_name)
    if interpolate_parcel_temperature:
        variable_names = np.append(variable_names, Tvar_name)
        variable_names = np.append(variable_names,meteo_plves_var_name)

    errors = ""
    find_error = False
    for pfilename in meteo_filenames:
        try:
            nc = Dataset(f"{path_to_meteodata}/{pfilename}")
            for varname in variable_names:
                if not varname in nc.variables.keys():
                    errors = (
                        errors
                        + f"     (**) ERROR: variable {varname} is not in file: "
                        + pfilename
                        + "\n"
                    )
                    find_error = True

            nc.close()
        except:
            errors = (
                errors
                + "     (**) ERROR: No such file or directory: "
                + pfilename
                + "\n"
            )
            find_error = True

    return errors, find_error




def read_climT_file(
    path_clim_temperature,
    climT_filename,
    Tlat_var_name,
    Tlon_var_name,
    climTvar_name,
    Tvar_name,
    dTdp_var_name,
    dTdt_var_name,
    Tplves_var_name,
    start_date,
):
    """
    Reads climate temperature data from a NetCDF file and returns temperature-related variables.

    Parameters:
        path_clim_temperature (str): Path to the directory containing the climate temperature files.
        climT_filename (str): Name of the climate temperature NetCDF file.
        Tlat_var_name (str): Variable name for latitude in the NetCDF file.
        Tlon_var_name (str): Variable name for longitude in the NetCDF file.
        climTvar_name (str): Variable name for temperature in the NetCDF file.
        dTdp_var_name (str): Variable name for temperature gradient with respect to pressure in the NetCDF file.
        dTdt_var_name (str): Variable name for temperature gradient with respect to time in the NetCDF file.
        Tplves_var_name (str): Variable name for pressure levels in the NetCDF file.
        start_date (datetime): Start date for the data extraction.

    Returns:
        tuple: Contains temperature (T), temperature gradient dT/dp (dT_dp), temperature gradient dT/dt (dT_dt),
               latitude (lat), longitude (lon), and pressure levels (plv) arrays.
    """

    nc = Dataset(f"{path_clim_temperature}/{climT_filename}", "r")
    climT = nc.variables[climTvar_name][:]
    dT_dp = nc.variables[dTdp_var_name][:]
    dT_dt = nc.variables[dTdt_var_name][:]
    lat = nc.variables[Tlat_var_name][:]
    lon = nc.variables[Tlon_var_name][:]
    plv = nc.variables[Tplves_var_name][:]

    #if interpolate_sfc:
    #    psfc = nc.variables[psfc_var_name][:]
    #else:
    #    psfc = 0

    #if interpolate_parcel_temperature:
    #    T = nc.variables[Tvar_name][:]
    #else:
    #    T = 0

    nc.close()
    lon = np.where(lon >= 180, lon - 360, lon)

    if len(str(int(plv.max()))) < 5:
        plv = plv * 100

    #if len(str(int(psfc.max()))) < 5:
    #    psfc = psfc * 100

    return climT, dT_dp, dT_dt, lat, lon, plv



def read_meteo_files(psfc_var_name, Tvar_name,
                        path_to_meteodata,
                        meteolat_var_name,
                        meteolon_var_name,
                        meteo_plves_var_name,
                        meteo_filename,
                        interpolate_parcel_temperature,
                        interpolate_sfc):

    nc = Dataset(f"{path_to_meteodata}/{meteo_filename}", "r")
    lat = nc.variables[meteolat_var_name][:]
    lon = nc.variables[meteolon_var_name][:]



    if interpolate_sfc:
        psfc = nc.variables[psfc_var_name][:]
    else:
        psfc = 0

    if interpolate_parcel_temperature:
        T = nc.variables[Tvar_name][:]
        plevs = nc.variables[meteo_plves_var_name][:]
    else:
        T = 0
        plevs=[0]
    nc.close()

    lon = np.where(lon >= 180, lon - 360, lon)

    if len(str(int(psfc.max()))) < 5:
        psfc = psfc * 100

    if len(psfc.shape)>2:
        psfc=psfc[0,:]

    return psfc, T, lat, lon, plevs


def parallel_interpolation(
    local_airtrajs,
    path_clim_temperature,
    localfnames,
    Tlat_var_name,
    Tlon_var_name,
    climTvar_name,
    Tvar_name,
    dTdp_var_name,
    dTdt_var_name,
    Tplves_var_name,
    psfc_var_name,
    path_to_meteodata,
    meteo_filenames,
    meteolat_var_name,
    meteolon_var_name,
    meteo_plves_var_name,
    start_date,
    track_times,
    interpolate_parcel_temperature,
    interpolate_sfc,
    comm,
    rank,
    size,
):
    """
    Interpolates temperature fields from climate temperature files in parallel.

    Parameters:
        local_airtrajs (ndarray): Array with shape (n_times, n_points, n_variables) containing the air parcel trajectories.
        path_clim_temperature (str): Path to the directory containing the climate temperature files.
        localfnames (list): List of strings with the filenames of the climate temperature files.
        Tlat_var_name (str): Variable name for latitude in the NetCDF files.
        Tlon_var_name (str): Variable name for longitude in the NetCDF files.
        climTvar_name (str): Variable name for temperature in the NetCDF files.
        dTdp_var_name (str): Variable name for temperature gradient with respect to pressure in the NetCDF files.
        dTdt_var_name (str): Variable name for temperature gradient with respect to time in the NetCDF files.
        Tplves_var_name (str): Variable name for pressure levels in the NetCDF files.
        start_date (datetime): Start date of the backward analysis.
        track_times (list): List of datetime objects with the dates of the postion of the trajectories.
        comm (MPI communicator): MPI communicator object.
        rank (int): Rank of the current process.
        size (int): Total number of processes.

    Returns:
        ndarray: Array with shape (n_times, n_points, n_variables) containing the interpolated temperature fields for each trajectory.
    """
    n_times = local_airtrajs.shape[0]

    # Divide the work among ranks
    local_indices = np.array_split(np.arange(n_times), size)[rank]

    local_results = np.zeros_like(local_airtrajs[:, :, :])
    # local_results = local_airtrajs.copy()

    chunk_size = local_airtrajs.shape[1]
    num_chunks = 1

    for i in local_indices:

        # Read climate data for the current trajectory
        (
            climTemperature_field,
            dT_dp,
            dT_dt,
            Tlatitude,
            Tlongitude,
            Tpressure_levels,
        ) = read_climT_file(
            path_clim_temperature,
            localfnames[i],
            Tlat_var_name,
            Tlon_var_name,
            climTvar_name,
            Tvar_name,
            dTdp_var_name,
            dTdt_var_name,
            Tplves_var_name,
            start_date,

        )

        if interpolate_sfc or interpolate_parcel_temperature:
            (surface_pressure,
            temperature_field,
            meteo_latitude,
            meteo_longitude,
            meteo_plves) = read_meteo_files(psfc_var_name, Tvar_name,
                        path_to_meteodata,
                        meteolat_var_name,
                        meteolon_var_name,
                        meteo_plves_var_name,
                        meteo_filenames[i],
                        interpolate_parcel_temperature,
                        interpolate_sfc)

        else:
            surface_pressure=0
            temperature_field=0
            meteo_latitude=Tlatitude
            meteo_longitude=Tlongitude
            meteo_plves=Tpressure_levels
        # Perform interpolation
        # print(f"             .... Rank {rank}: Interpolating track {localfnames[i]} -> ({track_times[i]})")
        for chunk_idx in range(num_chunks):
            start_idx = chunk_idx * chunk_size
            end_idx = min(start_idx + chunk_size, local_airtrajs.shape[1])

            # Prepare parcel positions
            # parcel_positions = np.zeros((local_airtrajs.shape[1], 3))
            parcel_positions = np.zeros(
                (len(local_airtrajs[i, start_idx:end_idx, 0]), 3)
            )
            parcel_positions[:, 0] = local_airtrajs[i, start_idx:end_idx, 2]  # Latitude
            parcel_positions[:, 1] = local_airtrajs[
                i, start_idx:end_idx, 1
            ]  # Longitude
            parcel_positions[:, 2] = local_airtrajs[
                i, start_idx:end_idx, 13
            ]  # Pressure

            interpTclimVals, interpdTdp, interpdTdt, intpsfc, intTemperature = (
                interpolate_temperature_field_regulargrid(
                    parcel_positions,
                    Tlatitude,
                    Tlongitude,
                    Tpressure_levels,
                    climTemperature_field,
                    dT_dp,
                    dT_dt,
                    surface_pressure,
                    meteo_latitude,
                    meteo_longitude,
                    meteo_plves,
                    interpolate_parcel_temperature,
                    interpolate_sfc,
                    temperature_field,
                    interpo_method="linear",
                    chunk_size=500000,
                )
            )

            local_results[i, start_idx:end_idx, -4] = interpTclimVals
            local_results[i, start_idx:end_idx, -3] = interpdTdp
            local_results[i, start_idx:end_idx, -2] = interpdTdt

            if interpolate_parcel_temperature:
                local_results[i, start_idx:end_idx, 9] = intTemperature
            else:
                local_results[i, start_idx:end_idx, 9] = local_airtrajs[
                    i, start_idx:end_idx, 9
                ]

            local_results[i, start_idx:end_idx, 11] = compute_theta(
                local_airtrajs[i, start_idx:end_idx, 6],
                local_airtrajs[i, start_idx:end_idx, 3],
                local_results[i, start_idx:end_idx, 9],
                press_data=True,
                parpres=local_airtrajs[i, start_idx:end_idx, 13],
            )

            if interpolate_sfc:
                local_results[i, start_idx:end_idx, -1] = intpsfc

            del parcel_positions
            if "interpTclimVals" in locals():
                del interpTclimVals, interpdTdp, interpdTdt, intpsfc, intTemperature
            gc.collect()

        # Clean up climate data after processing trajectory
        del climTemperature_field, temperature_field, dT_dp, dT_dt, surface_pressure
        del Tlatitude, Tlongitude, Tpressure_levels
        gc.collect()

    # Gather results from all ranks at rank 0
    global_results = None
    if rank == 0:
        global_results = np.empty_like(local_results)

    # Perform reduce operation: sum the local results across all processors into global_results
    comm.Reduce(local_results, global_results, op=MPI.SUM, root=0)

    # Ensure all ranks have the same global_results shape for Bcast
    if rank == 0:
        # Flatten the results to a 1D array (buffer)
        global_results_flat = global_results.flatten()
    else:
        global_results_flat = np.empty_like(
            local_results.flatten()
        )  # Create an empty array with the correct size

    # Broadcast the flattened array from rank 0 to all other ranks
    comm.Bcast(global_results_flat, root=0)

    # Reshape the results back to the original shape
    global_results = global_results_flat.reshape(local_results.shape)
    comm.Barrier()

    return global_results


def compute_var_integarated_day_heatRP(array, t, lon, lat, numPdY, numPdX):
    """
    Compute integrated heat for daily averages in RP simulation.

    Parameters
    ----------
    array : 3D numpy array
        Array with shape (t, x, y) containing parcel positions and properties.
    t : int
        Number of time steps.
    lon : 2D numpy array
        Array with shape (numPdY+1, numPdX+1) containing longitudes.
    lat : 2D numpy array
        Array with shape (numPdY+1, numPdX+1) containing latitudes.
    numPdY : int
        Number of grid points in y-direction.
    numPdX : int
        Number of grid points in x-direction.

    Returns
    -------
    array_day : 3D numpy array
        Array with shape (t, numPdY, numPdX) containing integrated heat for daily averages.
    """
    dimX, dimY = lat.shape
    array_day = np.zeros((t, dimX - 1, dimY - 1))
    array[np.isnan(array)] = -999.0
    for ii in range(t):
        tmp_array = np.array(
            compute_grid_integrated_temperature_anom(
                array[ii, :, :], lon, lat, numPdY, numPdX, len(array[0, :, 0])
            ),
            dtype=np.float64,
        )
        array_day[ii, :, :] = tmp_array
    return array_day


def process_by_batches(tensor_org, batch_size=10):
    """Process in small batches to avoid memory explosion"""
    # Replace sentinel values in-place first
    tensor_org[tensor_org == -999.9] = -999.0

    # Initialize output list to collect results
    result_batches = []

    for i in range(0, tensor_org.shape[0], batch_size):
        end_idx = min(i + batch_size, tensor_org.shape[0])

        # Process small batch
        batch = tensor_org[i:end_idx]
        additional_cols = np.full(
            (batch.shape[0], batch.shape[1], 5), -999.0, dtype=batch.dtype
        )
        batch_result = np.concatenate([batch, additional_cols], axis=2)
        result_batches.append(batch_result)

        # Force garbage collection
        del batch, additional_cols

    # Combine all batches
    airtrajs_org = np.concatenate(result_batches, axis=0)
    del result_batches

    return airtrajs_org


def interpolate_to_parcel_trajectories(
    path_clim_temperature,
    climT_filenames,
    Tlat_var_name,
    Tlon_var_name,
    climTvar_name,
    Tvar_name,
    dTdp_var_name,
    dTdt_var_name,
    Tplves_var_name,
    psfc_var_name,
    path_to_meteodata,
    meteo_filenames,
    meteolat_var_name,
    meteolon_var_name,
    meteo_plves_var_name,
    comm,
    rank,
    size,
    start_date,
    track_times,
    tensor_org,
    interpolate_parcel_temperature,
    interpolate_sfc,
):







    # tensor_org[tensor_org==-999.9]=-999.
    # airtrajs_org = np.full((tensor_org.shape[0],tensor_org.shape[1],tensor_org.shape[2]+5), -999.)
    # airtrajs_org[:,:,:-5] = tensor_org[:,:,:]

    """
    Interpolates climatological temperature and related fields to air parcel trajectories.

    This function processes air parcel trajectory data by interpolating climatological 
    temperature, temperature gradients, and surface pressure from climate temperature files.
    It operates in parallel using MPI for distributed processing, and includes a real-time 
    timer to display elapsed time during interpolation.

    Parameters:
        path_clim_temperature (str): Path to the directory containing climate temperature files.
        climT_filenames (list): List of filenames for the climate temperature files.
        Tlat_var_name (str): Variable name for latitude in the NetCDF files.
        Tlon_var_name (str): Variable name for longitude in the NetCDF files.
        climTvar_name (str): Variable name for climatological temperature in the NetCDF files.
        Tvar_name (str): Variable name for temperature in the NetCDF files.
        dTdp_var_name (str): Variable name for temperature gradient with respect to pressure.
        dTdt_var_name (str): Variable name for temperature gradient with respect to time.
        Tplves_var_name (str): Variable name for pressure levels in the NetCDF files.
        psfc_var_name (str): Variable name for surface pressure in the NetCDF files.
        comm (MPI communicator): MPI communicator for parallel processing.
        rank (int): Rank of the current MPI process.
        size (int): Total number of MPI processes.
        start_date (datetime): Start date of the backward analysis.
        track_times (list): List of datetime objects for trajectory positions.
        tensor_org (ndarray): Original tensor containing air parcel trajectory data.
        interpolate_parcel_temperature (bool): Whether to interpolate parcel temperature.
        interpolate_sfc (bool): Whether to interpolate surface pressure.

    Returns:
        ndarray: Updated tensor with interpolated temperature and related fields.
    """

    airtrajs_org = process_by_batches(tensor_org, batch_size=2)  # Very small batches

    if rank == 0:
        print(
            "      --- Interpolating climatological Temperature, dTclim_dp, dTclim_dt and surface pressure to air parcels trajectories"
        )
    start_time_int = time.time()
    stop_event = Event()

    def print_timer():
        while not stop_event.is_set():
            elapsed = int(time.time() - start_time_int)
            if rank == 0:
                secunits = "seconds"
                print(f"\r          ... Elapsed Time: {elapsed} {secunits}", end="")
            time.sleep(1)

        # Print completion message after the timer ends
        if rank == 0:
            print(
                f"\r          ... Interpolation Done: (Elapsed Time: {elapsed} {secunits}).                 "
            )

    # Run the interpolation function and print a timer simultaneously

    timer_thread = Thread(target=print_timer, daemon=True)
    timer_thread.start()

    tmp_airtrajs = parallel_interpolation(
        airtrajs_org,
        path_clim_temperature,
        climT_filenames,
        Tlat_var_name,
        Tlon_var_name,
        climTvar_name,
        Tvar_name,
        dTdp_var_name,
        dTdt_var_name,
        Tplves_var_name,
        psfc_var_name,
        path_to_meteodata,
        meteo_filenames,
        meteolat_var_name,
        meteolon_var_name,
        meteo_plves_var_name,
        start_date,
        track_times,
        interpolate_parcel_temperature,
        interpolate_sfc,
        comm,
        rank,
        size,
    )
    airtrajs_org[:, :, -4:] = tmp_airtrajs[:, :, -4:]

    if interpolate_parcel_temperature:
        airtrajs_org[:, :, 9] = tmp_airtrajs[:, :, 9]
    airtrajs_org[:, :, 11] = tmp_airtrajs[:, :, 11]

    stop_event.set()
    timer_thread.join()

    return airtrajs_org


def gettting_temperature_anom_budget(
    airtrajs_org,
    cenlon,
    dtime,
    lat_vals,
    lon_vals,
    analysis_levels,
    dp_sfc,
    dp_upper,
    Tanom_threshold,
    ntimes,
    start_date,
    numPdX,
    numPdY,
    rank,
):
   
    """
    Computes temperature anomaly budget for air trajectories over specified analysis levels.

    Args:
        airtrajs_org (np.ndarray): Original air trajectory data with dimensions [time, parcels, variables].
        cenlon (float or str): Central longitude for mean longitude computation or '180'.
        dtime (int): Time step in minutes for trajectory data.
        lat_vals (np.ndarray): Latitude grid values for spatial analysis.
        lon_vals (np.ndarray): Longitude grid values for spatial analysis.
        analysis_levels (list): List of analysis levels ('pbl', 'sfc', or pressure levels as strings).
        dp_sfc (float): Pressure difference threshold for surface analysis in hPa.
        dp_upper (float): Pressure difference threshold for upper levels in hPa.
        Tanom_threshold (list): Temperature anomaly thresholds for each analysis level.
        ntimes (int): Number of time steps for trajectory analysis.
        start_date (str): Start date of the trajectory analysis period.
        numPdX (int): Number of grid points in the x-direction for trajectory integration.
        numPdY (int): Number of grid points in the y-direction for trajectory integration.
        rank (int): Rank of the process in a parallel computing environment.

        
    This function is adapted from Röthlisberger & Papritz (2023) and Papritz and Röthlisberger (2023)

    Röthlisberger, M., & Papritz, L. (2023). Quantifying the physical processes leading to atmospheric
    hot extremes at a global scale. Nature Geoscience, 16(3), 210-216. DOI: 10.1038/s41561-023-01126-1

    Returns:
        dict: A dictionary containing gridded temperature anomaly matrix, parts matrix,
              genesis matrix, parcel properties, and counters for each analysis level.
    """

    dt = dtime * 60
    nts = airtrajs_org.shape[0] - 2

    dict_output = {}
    for ilevel, level in enumerate(analysis_levels):
        if rank == 0:
            print(f"      --- Processing level: {level}              ")
        tanom0 = Tanom_threshold[ilevel]
        if level == "pbl":

            check_condition = is_withinpbl(
                airtrajs_org[-1, :, 4], airtrajs_org[1:, :, 7], 0, 0, "maxval"
            )

        elif level == "sfc":

            diff = np.abs(airtrajs_org[-1, :, 18] - airtrajs_org[-1, :, 13])
            check_condition = diff <= dp_sfc * 100
            check_condition[:]=True #delete

        else:
            check_condition = (
                float(level) * 100 - dp_upper * 100 <= airtrajs_org[-1, :, 13]
            ) & (airtrajs_org[-1, :, 13] <= float(level) * 100 + dp_upper * 100)

        airtrajs = check_filtering_condition(airtrajs_org, check_condition)

        tanom_vector = np.expand_dims(airtrajs[:, :, 9] - airtrajs[:, :, 15], axis=-1)

        # Stack the new column along the last axis
        airtrajs_with_Tanom_tmp = np.concatenate((airtrajs, tanom_vector), axis=-1)


        tanom_condition = (airtrajs_with_Tanom_tmp[-1, :, -1] > tanom0) & (
            airtrajs_with_Tanom_tmp[-1, :, 9] != -999.0
        )
        tanom_condition[:]=True #delete

        airtrajs_with_Tanom = check_filtering_condition(
            airtrajs_with_Tanom_tmp, tanom_condition
        )

        aux_airtrajsT_ = airtrajs_with_Tanom[::-1, :, :]

        aux_airtrajsT = aux_airtrajsT_[:-1, :, :].copy()

        aux_airtrajsT[:, :, 14] = cal_track_diff_forward(
            aux_airtrajsT_[:, :, 13], "diff"
        ) / (dt)

        proc_properties = [
            "lon",
            "lat",
            "lon_p",
            "lat_p",
            "pressure",
            "dist_traj",
            "T_anom",
            "seas_i",
            "adv_i",
            "adiab1_i",
            "adiab2_i",
            "adiab3_i",
            "diab_i",
        ]

        # matrix_heat = np.zeros((aux_airtrajsT.shape[0],aux_airtrajsT.shape[1], len(proc_properties)))-999.
        matrix_heat = (
            np.zeros(
                (aux_airtrajsT.shape[0], aux_airtrajsT.shape[1], len(proc_properties))
            )
            - 999.0
        )

        gen_properties = [
            "lon",
            "lat",
            "gen_lon",
            "gen_lat",
            "gen_p",
            "res1",
            "res2",
            "age",
            "dist",
            "delta_p_tmp",
            "cont",
        ]
        matrix_gen = np.zeros((aux_airtrajsT.shape[1], len(gen_properties))) - 999.0

        back_hours = get_hours_from_start_date(start_date, ntimes)[::-1]

        
        for i in range(0, aux_airtrajsT.shape[1]):
            # temperature anomaly index is 19
            T_anom = aux_airtrajsT[:, i, 19]

            # Same sign of T' for the entire trajectory...
            if len(np.unique(np.sign(T_anom))) == 1:
                genesis_i = nts - 1
            # T' changes sign at least oonce during the trajectory time.
            else:
                genesis_i = np.where(np.sign(T_anom[0]) != np.sign(T_anom))[0][0] - 1
                if genesis_i < 0:
                    genesis_i = 0

            # shorten trajectory to time t > t_g
            short_traj = aux_airtrajsT[: genesis_i + 1, i, :]
            T_anom = T_anom[: genesis_i + 1]

            # compute "genesis characteristics"
            gen_lat = short_traj[genesis_i, 2]
            gen_lon = short_traj[genesis_i, 1]
            gen_p = short_traj[genesis_i, 13]
            
            delta_p_tmp = short_traj[0, 13] - short_traj[genesis_i, 13]
            age = np.abs(back_hours[genesis_i]) 

            full_lons = short_traj[:, 1]
            full_lats = short_traj[:, 2]

            plon = compute_mean_lon(short_traj[:, 1], cenlon)
            if cenlon == "180":
                plon = (plon + 360) % 360

                plon[plon < 0] = -999.0

                full_lons = (full_lons + 360) % 360

                full_lons[full_lons < 0] = -999.0

            # quit()
            p_longitude = plon
            p_latitude = cal_track_diff_forward(short_traj[:, 2], "mean")
            if genesis_i == 0:
                p_latitude = np.array([short_traj[0, 2]])
                p_longitude = np.array([short_traj[0, 1]])

            P_presure_m = cal_track_diff_forward(short_traj[:, 13], "mean")  # (Pa)
            P_presure = short_traj[0 : (genesis_i) + 1, 13]
            p_Temp = cal_track_diff_forward(short_traj[:, 9], "mean")  # (K)

            p_Theta = cal_track_diff_forward(short_traj[:, 11], "mean")  # (K)
            p_climT = cal_track_diff_forward(short_traj[:, 15], "mean")  # (K)
            p_ddp_T_clim = short_traj[0:(genesis_i), 16]
            p_ddt_T_clim = short_traj[0:(genesis_i), 17]

            p_climT_m = (
                short_traj[1 : (genesis_i + 1), 15] + short_traj[0:(genesis_i), 15]
            ) / 2
            p_omega_m = short_traj[0:(genesis_i), 14]
            p_ddt_T_clim_m = (
                short_traj[1 : (genesis_i + 1), 17] + short_traj[0:(genesis_i), 17]
            ) / 2
            p_ddp_T_clim_m = (
                short_traj[1 : (genesis_i + 1), 16] + short_traj[0:(genesis_i), 16]
            ) / 2

            DT_anom_Dt_m = (T_anom[0:(genesis_i)] - T_anom[1 : (genesis_i + 1)]) / dt
            DT_clim_Dt_m = (
                short_traj[0:(genesis_i), 15] - short_traj[1 : (genesis_i + 1), 15]
            ) / dt
            DTH_Dt_m = (
                short_traj[0:(genesis_i), 11] - short_traj[1 : (genesis_i + 1), 11]
            ) / dt  # (K s-1)
            T_anom_m = cal_track_diff_forward(T_anom[: genesis_i + 1], "mean")

            # compute terms in Eq. (2) of Röthlisberger & Papritz (2023)

            #   res1, T' at the (discrete) genesis time.
            res1 = T_anom[-1]

            #   seasonality T'
            seas_m = -p_ddt_T_clim_m

            #   adiab 1
            adiab1_m = np.array(-p_ddp_T_clim_m * p_omega_m)

            #   adv, computed from DT_clim_DT and adiab1
            adv_m = -(DT_clim_Dt_m) + p_ddt_T_clim_m - adiab1_m

            #   adiab 2
            adiab2_m = (kappa / P_presure_m) * T_anom_m * p_omega_m

            #   adiab 3
            adiab3_m = (kappa / P_presure_m) * p_climT_m * p_omega_m

            #   diab
            diab_m = ((P_presure_m / p0) ** kappa) * DTH_Dt_m

            # integrating these quantities along the trajectory in the forward
            #   direction. Note that here the time dimension flips to forward directed.
            if genesis_i > 0:
                seas_i = np.array(
                    [seas_m[-1]]
                    + [np.sum(seas_m[-k:] * dt) for k in range(1, genesis_i + 1)]
                )
                adv_i = np.array(
                    [adv_m[-1]]
                    + [np.sum(adv_m[-k:] * dt) for k in range(1, genesis_i + 1)]
                )
                adiab1_i = np.array(
                    [adiab1_m[-1]]
                    + [np.sum(adiab1_m[-k:] * dt) for k in range(1, genesis_i + 1)]
                )
                adiab2_i = np.array(
                    [adiab2_m[-1]]
                    + [np.sum(adiab2_m[-k:] * dt) for k in range(1, genesis_i + 1)]
                )
                adiab3_i = np.array(
                    [adiab3_m[-1]]
                    + [np.sum(adiab3_m[-k:] * dt) for k in range(1, genesis_i + 1)]
                )
                diab_i = np.array(
                    [diab_m[-1]]
                    + [np.sum(diab_m[-k:] * dt) for k in range(1, genesis_i + 1)]
                )

                seas_i = seas_i[::-1]
                adv_i = adv_i[::-1]
                adiab1_i = adiab1_i[::-1]
                adiab2_i = adiab2_i[::-1]
                adiab3_i = adiab3_i[::-1]
                diab_i = diab_i[::-1]

                res2 = T_anom[0] - (
                    res1
                    + seas_i[0]
                    + adv_i[0]
                    + adiab1_i[0]
                    + adiab2_i[0]
                    + adiab3_i[0]
                    + diab_i[0]
                )

                cont = 1
            #   special case: the anomaly is generated at the last trajectory time step,
            #   i.e., has age 0 hours  ...

            else:
                seas_i = [0]
                adv_i = [0]
                adiab1_i = [0]
                adiab2_i = [0]
                adiab3_i = [0]
                diab_i = [0]
                res2 = 0
                cont = 0

            # compute (gc)dist for all trajectoory time steps
            dist_traj = gc_dist(
                full_lats, full_lons, full_lats[0], full_lons[0], 6371000
            )
            # the formation distance is the distance between genesis and event.
            dist = dist_traj[-1]

            matrix_heat[0 : len(full_lons), i, 0] = full_lons[0]
            matrix_heat[0 : len(full_lats), i, 1] = full_lats[0]
            matrix_heat[0 : len(full_lons), i, 2] = full_lons
            matrix_heat[0 : len(full_lats), i, 3] = full_lats
            matrix_heat[0 : len(P_presure), i, 4] = P_presure
            matrix_heat[0 : len(dist_traj), i, 5] = dist_traj / 1000
            matrix_heat[0 : len(T_anom), i, 6] = T_anom
            save_prop = [
                "seas_i",
                "adv_i",
                "adiab1_i",
                "adiab2_i",
                "adiab3_i",
                "diab_i",
            ]

            xval = 0
            for prop in save_prop:
                index = proc_properties.index(prop)

                if len(np.unique(np.sign(vars()[prop]))) == 1:
                    genesis_pi = len(vars()[prop]) - 1
                # Tp' changes sign at least oonce during the trajectory time.
                else:
                    genesis_pi = (
                        np.where(np.sign(vars()[prop][0]) != np.sign(vars()[prop]))[0][
                            0
                        ]
                        - 1
                    )
                    if genesis_pi < 0:
                        genesis_pi = 0

                if genesis_pi > 0:
                    vars()[prop][genesis_pi + 1 :] = -999.0

                matrix_heat[0 : len(vars()[prop]), i, index] = vars()[prop]

            matrix_gen[i, 0] = full_lons[0]
            matrix_gen[i, 1] = full_lats[0]
            # gen_properties=["lon","lat","gen_lon","gen_lat","gen_p","res1","res2","age","dist","delta_p_tmp","cont"]
            for prop in gen_properties[2:]:
                index = gen_properties.index(prop)
                matrix_gen[i, index] = vars()[prop]

        ####Processing 3D properties (proc_properties)
        array_proc_3d = proc_properties[2:]
        dimX, dimY = lat_vals.shape
        final_matrix_Tanom = np.zeros((len(array_proc_3d), nts + 1, dimX - 1, dimY - 1))

        temp_matrix = np.copy(matrix_heat[:, :, :3])

        for i, prop in enumerate(array_proc_3d):

            temp_matrix[:, :, 2] = matrix_heat[:, :, i + 2]

            aux_data = compute_var_integarated_day_heatRP(
                temp_matrix, nts + 1, lon_vals, lat_vals, numPdY, numPdX
            )
            final_matrix_Tanom[i, :, :, :] = aux_data

        ####Processing 2D properties (gen_properties)
        array_proc_2d = gen_properties[2:]

        final_matrix_Tanom_2d = np.zeros((len(array_proc_2d), dimX - 1, dimY - 1))
        temp_matrix_gen = np.copy(matrix_gen[:, :3])

        # quit()
        for i, prop in enumerate(array_proc_2d):
            temp_matrix_gen[:, 2] = matrix_gen[:, i + 2]
            final_matrix_Tanom_2d[i, :, :] = np.array(
                compute_temperature_anom_genesis_properties(
                    temp_matrix_gen,
                    lon_vals,
                    lat_vals,
                    numPdY,
                    numPdX,
                    temp_matrix_gen.shape[0],
                ),
                dtype=np.float64,
            )

        contributing_trjs = len(matrix_gen[:, -1][matrix_gen[:, -1] == 1])

        save_matrix_heat = matrix_heat[:, :, 2:]

        dict_output[f"{level}_Tanom_matrix_gridded"] = final_matrix_Tanom[:, ::-1, :, :]
        dict_output[f"{level}_Tanom_matrix_parts"] = save_matrix_heat[::-1, :, :]
        dict_output[f"{level}_gen_matrix"] = final_matrix_Tanom_2d
        dict_output[f"{level}_parcels_properties"] = array_proc_3d
        dict_output[f"{level}_counter_part_Tanom"] = temp_matrix_gen.shape[0]
        dict_output[f"{level}_no_Tanom_parts"] = (
            temp_matrix_gen.shape[0] - contributing_trjs
        )

    return dict_output


def check_filtering_condition(matrix_org, check_condition):
    """
    Filters the input matrix based on a boolean condition array.

    Parameters
    ----------
    matrix_org : np.ndarray
        Original 3D matrix of shape (X, Y, Z) to be filtered.
    check_condition : np.ndarray
        Boolean array of shape (Y,) representing the condition for filtering.
        Only elements with True value in this array will be retained.

    Returns
    -------
    matrix : np.ndarray
        Filtered 3D matrix of shape (X, Y_filtered, Z) where Y_filtered
        is the number of True values in `check_condition`. If no condition
        is True, returns a matrix of shape (X, 1, Z) filled with zeros.

    Notes
    -----
    - This function assumes that the second dimension of `matrix_org` corresponds
      to the dimension being filtered by `check_condition`.
    - The function uses numpy's advanced indexing to efficiently filter
      the matrix.
    """

    if len(check_condition[check_condition == True]) == 0:
        matrix = np.empty((matrix_org.shape[0], 1, matrix_org.shape[2]))
        matrix[:, :, :] = 0

    else:
        matrix = np.empty(
            (
                matrix_org.shape[0],
                len(check_condition[check_condition == True]),
                matrix_org.shape[2],
            )
        )

        for i in range(0, matrix.shape[0]):
            matrix[i, :] = matrix_org[i, :, :][check_condition == True]
    return matrix


def get_hours_from_start_date(start_date, track_times):
    """
    Compute the time difference between the start_date and a list of track_times
    in hours.

    Parameters
    ----------
    start_date : datetime
        The start date of the analysis
    track_times : list of datetime
        The dates of the forecast times

    Returns
    -------
    dt : list of float
        The time differences in hours
    """
    dt = [(itime - start_date) / timedelta(hours=1) for itime in track_times]

    return dt


def gettting_temperature_anom_budget_sources(
    airtrajs_org,
    cenlon,
    dtime,
    lat_vals,
    lon_vals,
    analysis_levels,
    dp_sfc,
    dp_upper,
    Tanom_threshold,
    ntimes,
    start_date,
    Tanom_linear_adjustment,
    area,
    lag_times,
    numPdX,
    numPdY,
    rank,
):
    

    
    """
    Computes temperature anomaly budget sources from air trajectory data.

    This function processes air trajectory data to compute various temperature anomaly
    terms and their contributions, specifically targeting different atmospheric levels
    (e.g., pbl, sfc, and others). It performs filtering based on specified conditions
    and integrates various terms over time to evaluate their contributions to temperature
    anomalies.

    Parameters
    ----------
    airtrajs_org : np.ndarray
        Original air trajectory data of shape (time, parcels, variables).
    cenlon : str
        Central longitude reference point ("0" or "180").
    dtime : int
        Time step in minutes.
    lat_vals : np.ndarray
        Latitude values for gridding.
    lon_vals : np.ndarray
        Longitude values for gridding.
    analysis_levels : list of str
        List of atmospheric levels to analyze (e.g., ["pbl", "sfc"]).
    dp_sfc : float
        Pressure difference threshold for surface level analysis.
    dp_upper : float
        Pressure difference threshold for upper levels.
    Tanom_threshold : list of float
        Threshold values for temperature anomaly at each level.
    ntimes : list of datetime
        Time points corresponding to each trajectory step.
    start_date : datetime
        Starting date for the analysis.
    Tanom_linear_adjustment : bool
        Flag indicating whether to apply linear adjustment to temperature anomaly.
    area : float
        Area of the grid cell for normalization.
    lag_times : list of int
        Time lags for integration.
    numPdX : int
        Number of grid points in X direction.
    numPdY : int
        Number of grid points in Y direction.
    rank : int
        Rank of the current process (0 for master process).

        


    Returns
    -------
    dict_output : dict
        Dictionary containing processed matrices and contribution metrics for each 
        specified atmospheric level. Includes integrated day matrices, trajectory 
        matrices, gridded process anomaly masks, and contribution statistics.

    Notes
    -----
    - This function is designed for parallel processing, with the `rank` parameter
      indicating whether the current process is master (rank 0) or a worker.
    - This function is adapted from Röthlisberger & Papritz (2023) and Papritz and Röthlisberger (2023)

    Papritz, L., & Röthlisberger, M. (2023). A novel temperature anomaly source diagnostic: Method and application to the 2021 heatwave in the Pacific Northwest. Geophysical Research Letters, 50, e2023GL105641. https://doi.org/10.1029/2023GL105641
    """

    dt = dtime * 60

    nts = airtrajs_org.shape[0] - 2

    dict_output = {}
    for ilevel, level in enumerate(analysis_levels):
        if rank == 0:
            print(f"      --- Processing level: {level}              ")
        tanom0 = Tanom_threshold[ilevel]
        if level == "pbl":

            check_condition = is_withinpbl(
                airtrajs_org[-1, :, 4], airtrajs_org[1:, :, 7], 0, 0, "maxval"
            )

        elif level == "sfc":

            diff = np.abs(airtrajs_org[-1, :, 18] - airtrajs_org[-1, :, 13])
            check_condition = diff <= dp_sfc * 100
           

        else:
            check_condition = (
                float(level) * 100 - dp_upper * 100 <= airtrajs_org[-1, :, 13]
            ) & (airtrajs_org[-1, :, 13] <= float(level) * 100 + dp_upper * 100)

        airtrajs = check_filtering_condition(airtrajs_org, check_condition)

        tanom_vector = np.expand_dims(airtrajs[:, :, 9] - airtrajs[:, :, 15], axis=-1)

        # Stack the new column along the last axis
        airtrajs_with_Tanom_tmp = np.concatenate((airtrajs, tanom_vector), axis=-1)
        filtered_parcels_level = airtrajs_with_Tanom_tmp.shape[1]

        tanom_condition = (airtrajs_with_Tanom_tmp[-1, :, -1] > tanom0) & (
            airtrajs_with_Tanom_tmp[-1, :, 9] != -999.0
        )

        airtrajs_with_Tanom = check_filtering_condition(
            airtrajs_with_Tanom_tmp, tanom_condition
        )

        aux_airtrajsT_ = airtrajs_with_Tanom[::-1, :, :]

        aux_airtrajsT = aux_airtrajsT_[:-1, :, :].copy()

        aux_airtrajsT[:, :, 14] = cal_track_diff_forward(
            aux_airtrajsT_[:, :, 13], "diff"
        ) / (
            dt
        )  # (computing omega as dp/dt, P/s )

        tanom_terms = ["tanom", "seas", "adv", "adiab", "diab"]

        matrix_trajs = (
            np.zeros(
                (aux_airtrajsT.shape[0], aux_airtrajsT.shape[1], len(tanom_terms) + 3)
            )
            - 999.0
        )
        matrix_trajs_ = (
            np.zeros((aux_airtrajsT.shape[0], aux_airtrajsT.shape[1], 2)) - 999.0
        )

        proc_properties = ["lon_p", "lat_p", "variable_pos"]
        matrix_heat = np.zeros(
            (
                len(tanom_terms),
                aux_airtrajsT.shape[0] - 1,
                aux_airtrajsT.shape[1],
                len(proc_properties) + 1,
            )
        )
        matrix_heat_neg = np.zeros(
            (
                len(tanom_terms),
                aux_airtrajsT.shape[0] - 1,
                aux_airtrajsT.shape[1],
                len(proc_properties) + 1,
            )
        )

        gen_properties = ["gen_lon", "gen_lat", "gen_p", "age", "dist", "res", "cont"]
        matrix_gen = np.zeros(
            (len(tanom_terms), aux_airtrajsT.shape[1], len(gen_properties))
        )

        back_hours = get_hours_from_start_date(start_date, ntimes)[::-1]
        Tanomomaly0 = []
        nconts = 0
        contrib_tanom_pos = []
        contrib_tanom_neg = []

        for pterm in tanom_terms[1:]:
            vars()[f"contrib_{pterm}_pos"] = []
            vars()[f"contrib_{pterm}_neg"] = []
            vars()[f"contrib_{pterm}"] = []
        vars()[f"contrib_tanom"] = []
        contrib_res = []

        
        for i in range(0, aux_airtrajsT.shape[1]):

            # temperature anomaly index is 19
            # temperature anomaly index is 19
            T_anom = aux_airtrajsT[:, i, 19]

            
            # Same sign of T' for the entire trajectory...
            if len(np.unique(np.sign(T_anom))) == 1:
                genesis_i = nts - 1
            # T' changes sign at least oonce during the trajectory time.
            else:
                genesis_i = np.where(np.sign(T_anom[0]) != np.sign(T_anom))[0][0] - 1
                if genesis_i < 0:
                    genesis_i = 0

            # shorten trajectory to time t > t_g
            short_traj = aux_airtrajsT[: genesis_i + 1, i, :]
            T_anom = T_anom[: genesis_i + 1]

            # print(T_anom)

            # compute "genesis characteristics"
            gen_lat = short_traj[genesis_i, 2]
            gen_lon = short_traj[genesis_i, 1]
            gen_p = short_traj[genesis_i, 13]
           
            delta_p_tmp = short_traj[0, 13] - short_traj[genesis_i, 13]
            age = np.abs(back_hours[genesis_i])  # ((genesis_i) * dt)

            full_lons = short_traj[:, 1]
            full_lats = short_traj[:, 2]

            plon = compute_mean_lon(short_traj[:, 1], cenlon)
            if cenlon == "180":
                plon = (plon + 360) % 360

                plon[plon < 0] = -999.0
                full_lons = (full_lons + 360) % 360

                full_lons[full_lons < 0] = -999.0

            plat = cal_track_diff_forward(short_traj[:, 2], "mean")
            if genesis_i == 0:
                plat = np.array([short_traj[0, 2]])
                plon = np.array([short_traj[0, 1]])

            P_presure_m = cal_track_diff_forward(short_traj[:, 13], "mean")  # (Pa)
            P_presure = short_traj[0 : (genesis_i) + 1, 13]
            p_Temp = cal_track_diff_forward(short_traj[:, 9], "mean")  # (K)

            p_Theta = cal_track_diff_forward(short_traj[:, 11], "mean")  # (K)
            p_climT = cal_track_diff_forward(short_traj[:, 15], "mean")  # (K)
            p_ddp_T_clim = short_traj[0:(genesis_i), 16]
            p_ddt_T_clim = short_traj[0:(genesis_i), 17]

            p_climT_m = (
                short_traj[1 : (genesis_i + 1), 15] + short_traj[0:(genesis_i), 15]
            ) / 2
            p_omega_m = short_traj[0:(genesis_i), 14]
            p_ddt_T_clim_m = (
                short_traj[1 : (genesis_i + 1), 17] + short_traj[0:(genesis_i), 17]
            ) / 2
            p_ddp_T_clim_m = (
                short_traj[1 : (genesis_i + 1), 16] + short_traj[0:(genesis_i), 16]
            ) / 2

            DT_anom_Dt_m = (T_anom[0:(genesis_i)] - T_anom[1 : (genesis_i + 1)]) / dt
            DT_clim_Dt_m = (
                short_traj[0:(genesis_i), 15] - short_traj[1 : (genesis_i + 1), 15]
            ) / dt
            DTH_Dt_m = (
                short_traj[0:(genesis_i), 11] - short_traj[1 : (genesis_i + 1), 11]
            ) / dt  # (K s-1)
            T_anom_m = cal_track_diff_forward(T_anom[: genesis_i + 1], "mean")

            # compute terms in Eq. (2) of Röthlisberger & Papritz (2023)

            #   res1, T' at the (discrete) genesis time.
            res1 = T_anom[-1]

            #   seasonality T'
            seas_m = -p_ddt_T_clim_m

            #   adiab 1
            adiab1_m = np.array(-p_ddp_T_clim_m * p_omega_m)

            #   adv, computed from DT_clim_DT and adiab1
            adv_m = -(DT_clim_Dt_m) + p_ddt_T_clim_m - adiab1_m

            #   adiab 2
            adiab2_m = (kappa / P_presure_m) * T_anom_m * p_omega_m

            #   adiab 3
            adiab3_m = (kappa / P_presure_m) * p_climT_m * p_omega_m

            adiab_m = adiab1_m + adiab2_m + adiab3_m
            #   diab
            diab_m = ((P_presure_m / p0) ** kappa) * DTH_Dt_m

            # integrating these quantities along the trajectory in the forward
            #   direction. Note that here the time dimension flips to forward directed.
            if genesis_i > 0:

                Tanomomaly0 = np.append(Tanomomaly0, T_anom[0])

                seas_i = np.array(
                    [seas_m[-1]]
                    + [np.sum(seas_m[-k:] * dt) for k in range(1, genesis_i + 1)]
                )
                adv_i = np.array(
                    [adv_m[-1]]
                    + [np.sum(adv_m[-k:] * dt) for k in range(1, genesis_i + 1)]
                )
                adiab_i = np.array(
                    [adiab_m[-1]]
                    + [np.sum(adiab_m[-k:] * dt) for k in range(1, genesis_i + 1)]
                )
                diab_i = np.array(
                    [diab_m[-1]]
                    + [np.sum(diab_m[-k:] * dt) for k in range(1, genesis_i + 1)]
                )

                seas_i = seas_i[::-1]
                adv_i = adv_i[::-1]
                adiab_i = adiab_i[::-1]
                diab_i = diab_i[::-1]

                res2 = T_anom[0] - (
                    res1 + seas_i[0] + adv_i[0] + adiab_i[0] + diab_i[0]
                )
                res = res1 + res2 + seas_i[0]
                cont = 1
                nconts = nconts + 1
                #   special case: the anomaly is generated at the last trajectory time step,
                #   i.e., has age 0 hours  ...

                contrib_res = np.append(contrib_res, ((res1 + res2) / T_anom[0]) * 100)

                vars()["contrib_tanom"] = np.append(
                    vars()["contrib_tanom"], (T_anom[0] / T_anom[0]) * 100
                )
                vars()["contrib_seas"] = np.append(
                    vars()["contrib_seas"], (seas_i[0] / T_anom[0]) * 100
                )
                vars()["contrib_adv"] = np.append(
                    vars()["contrib_adv"], (adv_i[0] / T_anom[0]) * 100
                )
                vars()["contrib_adiab"] = np.append(
                    vars()["contrib_adiab"], (adiab_i[0] / T_anom[0]) * 100
                )
                vars()["contrib_diab"] = np.append(
                    vars()["contrib_diab"], (diab_i[0] / T_anom[0]) * 100
                )

                matrix_trajs[: len(full_lons), i, 0] = full_lons
                matrix_trajs[: len(full_lats), i, 1] = full_lats
                matrix_trajs[: len(P_presure), i, 2] = P_presure
                matrix_trajs[: len(T_anom), i, 3] = T_anom
                matrix_trajs[: len(seas_i), i, 4] = seas_i
                matrix_trajs[: len(adv_i), i, 5] = adv_i
                matrix_trajs[: len(adiab_i), i, 6] = adiab_i
                matrix_trajs[: len(diab_i), i, 7] = diab_i

                matrix_trajs_[: len(full_lons), i, 0] = full_lons[0]
                matrix_trajs_[: len(full_lons), i, 1] = full_lats[0]
                # compute (gc)dist for all trajectoory time steps
                # the formation distance is the distance between genesis and event.
                dist_traj = gc_dist(
                    full_lats, full_lons, full_lats[0], full_lons[0], 6371000
                )
                dist = dist_traj[-1] / 1000

                # gen_properties=["gen_lon","gen_lat","gen_p","age","dist","cont"]
                # tanom_terms=["tanom","seas","adv","adiab","diab"]

                for prop in gen_properties:
                    index = gen_properties.index(prop)
                    matrix_gen[0, i, index] = vars()[prop]

                save_prop = ["seas_i", "adv_i", "adiab_i", "diab_i"]
                matrix_heat[0, 0 : len(plon), i, 0] = plon
                matrix_heat[0, 0 : len(plat), i, 1] = plat
                matrix_heat_neg[0, 0 : len(plon), i, 0] = plon
                matrix_heat_neg[0, 0 : len(plat), i, 1] = plat

                delta_Tanom = cal_track_diff_forward(T_anom, "diff")
                if np.all(T_anom >= 0):

                    if Tanom_linear_adjustment:
                        adjusted_Tanom = compute_linear_discounted_Tanom(
                            delta_Tanom[::-1]
                        )
                    else:
                        adjusted_Tanom = delta_Tanom.copy()

                    matrix_heat[0, 0 : len(adjusted_Tanom), i, 2] = (
                        adjusted_Tanom / T_anom[0]
                    )
                    matrix_heat[
                        0, 0 : len(adjusted_Tanom[adjusted_Tanom > 0]), i, 3
                    ] = 1

                    contrib_tanom_pos = np.append(
                        contrib_tanom_pos, (np.sum(adjusted_Tanom) / T_anom[0]) * 100
                    )

                elif np.all(T_anom <= 0):

                    if Tanom_linear_adjustment:
                        adjusted_Tanom_negative = (
                            compute_linear_discounted_Tanom_negative(delta_Tanom[::-1])
                        )
                    else:
                        adjusted_Tanom_negative = delta_Tanom.copy()
                    matrix_heat_neg[0, 0 : len(adjusted_Tanom_negative), i, 2] = (
                        adjusted_Tanom_negative / T_anom[0]
                    )
                    matrix_heat_neg[
                        0,
                        0 : len(adjusted_Tanom_negative[adjusted_Tanom_negative < 0]),
                        i,
                        3,
                    ] = 1

                    contrib_tanom_neg = np.append(
                        contrib_tanom_neg,
                        (np.sum(adjusted_Tanom_negative) / T_anom[0]) * 100,
                    )

                for iterm, pterm in enumerate(tanom_terms[1:]):

                    if len(np.unique(np.sign(vars()[f"{pterm}_i"]))) == 1:
                        genesis_pi = len(vars()[f"{pterm}_i"]) - 1
                    # Tp' changes sign at least oonce during the trajectory time.
                    else:
                        genesis_pi = (
                            np.where(
                                np.sign(vars()[f"{pterm}_i"][0])
                                != np.sign(vars()[f"{pterm}_i"])
                            )[0][0]
                            - 1
                        )
                        if genesis_pi < 0:
                            genesis_pi = 0

                    vars()[f"{pterm}_i"][genesis_pi + 1 :] = 0
                    pfull_lats = full_lats[0 : genesis_pi + 1]
                    pfull_lons = full_lons[0 : genesis_pi + 1]
                    pP_presure = P_presure[0 : genesis_pi + 1]
                    pgen_lat = full_lats[genesis_pi]
                    pgen_lon = full_lons[genesis_pi]
                    pgen_p = P_presure[genesis_pi]
                    page = np.abs(back_hours[genesis_pi])  #
                    pdist_traj = gc_dist(
                        pfull_lats, pfull_lons, pfull_lats[0], pfull_lons[0], 6371000
                    )
                    pdist = pdist_traj[-1] / 1000
                    pres = vars()[f"{pterm}_i"][-1]
                    if len(pfull_lons) > 0:
                        pcont = 1

                    for prop in gen_properties:
                        index = gen_properties.index(prop)
                        matrix_gen[iterm + 1, i, index] = vars()[f"p{prop}"]

                    plon = compute_mean_lon(pfull_lons, cenlon)
                    if cenlon == "180":
                        plon = (plon + 360) % 360

                        plon[plon < 0] = -999.0

                    plat = cal_track_diff_forward(pfull_lats, "mean")

                    matrix_heat[iterm + 1, 0 : len(plon), i, 0] = plon
                    matrix_heat[iterm + 1, 0 : len(plat), i, 1] = plat

                    matrix_heat_neg[iterm + 1, 0 : len(plon), i, 0] = plon
                    matrix_heat_neg[iterm + 1, 0 : len(plat), i, 1] = plat

                    delta_iterm = cal_track_diff_forward(vars()[f"{pterm}_i"], "diff")

                    if np.all(vars()[f"{pterm}_i"] >= 0):

                        if Tanom_linear_adjustment:
                            adjusted_iterm = compute_linear_discounted_Tanom(
                                delta_iterm[::-1]
                            )  ###the adjusted function receives the data from the genesis to the end
                        else:
                            adjusted_iterm = delta_iterm.copy()

                        matrix_heat[iterm + 1, 0 : len(adjusted_iterm), i, 2] = (
                            adjusted_iterm / T_anom[0]
                        )
                        matrix_heat[
                            iterm + 1, 0 : len(adjusted_iterm[adjusted_iterm > 0]), i, 3
                        ] = 1

                        vars()[f"contrib_{pterm.split('_')[0]}_pos"] = np.append(
                            vars()[f"contrib_{pterm.split('_')[0]}_pos"],
                            (np.sum(adjusted_iterm) / T_anom[0]) * 100,
                        )

                    elif np.all(vars()[f"{pterm}_i"] <= 0):

                        if Tanom_linear_adjustment:
                            adjusted_iterm_negative = (
                                compute_linear_discounted_Tanom_negative(
                                    delta_iterm[::-1]
                                )
                            )
                        else:
                            adjusted_iterm_negative = delta_iterm.copy()
                        matrix_heat_neg[
                            iterm + 1, 0 : len(adjusted_iterm_negative), i, 2
                        ] = (adjusted_iterm_negative / T_anom[0])
                        matrix_heat_neg[
                            iterm + 1,
                            0 : len(
                                adjusted_iterm_negative[adjusted_iterm_negative < 0]
                            ),
                            i,
                            3,
                        ] = 1

                        vars()[f"contrib_{pterm.split('_')[0]}_neg"] = np.append(
                            vars()[f"contrib_{pterm.split('_')[0]}_neg"],
                            (np.sum(adjusted_iterm_negative) / T_anom[0]) * 100,
                        )

            else:
                matrix_trajs[: len(full_lons), i, 0] = full_lons
                matrix_trajs[: len(full_lats), i, 1] = full_lats
                matrix_trajs[: len(P_presure), i, 2] = P_presure
                matrix_trajs[: len(T_anom), i, 3] = T_anom
                matrix_trajs[:, i, 4] = 0
                matrix_trajs[:, i, 5] = 0
                matrix_trajs[:, i, 6] = 0
                matrix_trajs[:, i, 7] = 0

                matrix_trajs_[: len(full_lons), i, 0] = full_lons[0]
                matrix_trajs_[: len(full_lons), i, 1] = full_lats[0]

        # tanom_terms=["tanom","seas","adv","adiab","diab"]

        matrix_day_integrated = np.zeros(
            (len(tanom_terms) * 2, len(lag_times) - 1, numPdY, numPdX)
        )

        mean_Tanom0 = np.mean(Tanomomaly0)

        iloc = 0
        positive_contributions = []
        negative_contributions = []
        overall_contributions = []
        for iterm, pterm in enumerate(tanom_terms):

            array_days = compute_var_integarated_day_Tanom(
                matrix_heat[iterm, ::-1, :, :],
                lag_times,
                area,
                lon_vals,
                lat_vals,
                numPdY,
                numPdX,
            )
            matrix_day_integrated[iloc, :, :] = (
                (array_days) * 10e6 * 100
            )  # mulpitply * 10e6 to convert from fractional K/m2 to fractional K/km2 and multiply *100 to convert from  fractional  K/km2 to %/km2

            array_days_neg = compute_var_integarated_day_Tanom(
                matrix_heat_neg[iterm, ::-1, :, :],
                lag_times,
                area,
                lon_vals,
                lat_vals,
                numPdY,
                numPdX,
            )
            matrix_day_integrated[iloc + 1, :, :] = (array_days_neg) * 10e6 * 100

            iloc = iloc + 2

            if len(vars()[f"contrib_{pterm}_pos"]) >= 1:
                positive_contributions = np.append(
                    positive_contributions, np.nanmedian(vars()[f"contrib_{pterm}_pos"])
                )
            else:
                positive_contributions = np.append(positive_contributions, np.nan)

            if len(vars()[f"contrib_{pterm}_neg"]) >= 1:
                negative_contributions = np.append(
                    negative_contributions, np.nanmedian(vars()[f"contrib_{pterm}_neg"])
                )
            else:
                negative_contributions = np.append(negative_contributions, np.nan)

            overall_contributions = np.append(
                overall_contributions, np.nanmedian(vars()[f"contrib_{pterm}"])
            )

        parcels_properties = ["T_anom", "seas", "adv", "adiab", "diab"]
        full_matrix_Tanom = np.zeros(
            (len(parcels_properties), len(lag_times), numPdY, numPdX)
        )

        matrix_trajs[: len(full_lons), i, 0] = full_lons

        temp_matrix = np.copy(matrix_trajs[:, :, :3])
        temp_matrix[:, :, 0] = matrix_trajs_[:, :, 0]
        temp_matrix[:, :, 1] = matrix_trajs_[:, :, 1]
        for i, prop in enumerate(parcels_properties):

            temp_matrix[:, :, 2] = matrix_trajs[:, :, i + 3]

            aux_data = compute_var_integarated_day_heatRP(
                temp_matrix, len(lag_times), lon_vals, lat_vals, numPdY, numPdX
            )
            full_matrix_Tanom[i, :, :, :] = aux_data

        gridded_process_anoms_mask = full_matrix_Tanom[:, 0, :, :]

        contrib_matrix = np.zeros((3, len(tanom_terms)))
        contrib_matrix[0, :] = overall_contributions
        contrib_matrix[1, :] = positive_contributions
        contrib_matrix[2, :] = negative_contributions

        dict_output[f"{level}_matrix_day_integrated"] = matrix_day_integrated
        dict_output[f"{level}_matrix_gen"] = matrix_gen[:, :, :-1]
        dict_output[f"{level}_matrix_trajs"] = matrix_trajs[::-1, :, :]
        dict_output[f"{level}_gridded_process_anoms_mask"] = gridded_process_anoms_mask

        dict_output[f"{level}_positive_contributions"] = positive_contributions
        dict_output[f"{level}_negative_contributions"] = negative_contributions
        dict_output[f"{level}_overall_contributions"] = overall_contributions
        dict_output[f"{level}_tanom_terms"] = tanom_terms
        dict_output[f"{level}_contribution_matrix"] = contrib_matrix

        filtered_parcels = aux_airtrajsT.shape[1]
        no_part_tanom = aux_airtrajsT.shape[1] - len(
            matrix_gen[0, :, -1][matrix_gen[0, :, -1] == 1]
        )
        no_part_seas = aux_airtrajsT.shape[1] - len(
            matrix_gen[1, :, -1][matrix_gen[1, :, -1] == 1]
        )
        no_part_adiab = aux_airtrajsT.shape[1] - len(
            matrix_gen[2, :, -1][matrix_gen[2, :, -1] == 1]
        )
        no_part_adv = aux_airtrajsT.shape[1] - len(
            matrix_gen[3, :, -1][matrix_gen[3, :, -1] == 1]
        )
        no_part_diab = aux_airtrajsT.shape[1] - len(
            matrix_gen[4, :, -1][matrix_gen[4, :, -1] == 1]
        )
        dict_output[f"{level}_filtered_parcels_within_level"] = filtered_parcels_level
        dict_output[f"{level}_filtered_parcels"] = filtered_parcels
        dict_output[f"{level}_no_part_tanom"] = no_part_tanom
        dict_output[f"{level}_no_part_seas"] = no_part_seas
        dict_output[f"{level}_no_part_adiab"] = no_part_adiab
        dict_output[f"{level}_no_part_adv"] = no_part_adv
        dict_output[f"{level}_no_part_diab"] = no_part_diab

    return dict_output


def compute_var_integarated_day_Tanom(array, t, area, lon, lat, numPdY, numPdX):
    """
    Compute integrated temperature anomalies for daily averages.

    Parameters
    ----------
    array : 3D numpy array
        Array with shape (t, x, y) containing parcel positions and properties.
    t : array-like
        Array of time indices.
    area : float
        Area over which the integration is performed.
    lon : 2D numpy array
        Array with shape (numPdY+1, numPdX+1) containing longitudes.
    lat : 2D numpy array
        Array with shape (numPdY+1, numPdX+1) containing latitudes.
    numPdY : int
        Number of grid points in y-direction.
    numPdX : int
        Number of grid points in x-direction.

    Returns
    -------
    array_day : 3D numpy array
        Array with shape (len(t)-1, numPdY, numPdX) containing integrated temperature anomalies for daily averages.
    """

    dimX, dimY = lat.shape
    array_day = np.empty((len(t) - 1, dimX - 1, dimY - 1))

    ndb = np.arange(len(t) - 1, 0, -1)
    for ii in range(len(t) - 1):

        moistd = np.array(
            compute_temperature_anomalies_sources(
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
        # array_day[int (ndb[ii]-1), :,:] = ajust_units(moistd, area, density, dtime, varid, "None")
        array_day[int(ndb[ii] - 1), :, :] = moistd / area
    return array_day


def compute_linear_discounted_Tanom(var):

    """
    Applies a linear discount to a time series of temperature anomalies.

    Parameters
    ----------
    var : 1D numpy array
        Array containing the time series of temperature anomalies.

    Returns
    -------
    result_var : 1D numpy array
        Array containing the time series with the linear discount applied.
    """
    var[var == -999.0] = 0
    var[np.isnan(var)] = 0

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

    return result_var[::-1]


def compute_linear_discounted_Tanom_negative(var):

    """
    Applies a linear discount to the temperature anomalies (Tanom) time series, considering only negative values.
    The discount is done by subtracting the absolute value of the current negative Tanom from the previous positive Tanom values.
    If the positive Tanom is smaller than the absolute value of the current negative Tanom, the positive Tanom is set to 0.
    The result is returned in reverse order.

    Parameters
    ----------
    var : 1D numpy array
        Array containing the time series of temperature anomalies.

    Returns
    -------
    result_var : 1D numpy array
        Array containing the time series with the linear discount applied.
    """
    var[var == -999.0] = 0
    var[np.isnan(var)] = 0
    result_var = np.empty_like(var)
    result_var[:] = 0

    for i in range(0, len(var)):
        if var[i] < 0:
            result_var[i] = var[i]

        else:
            suma = np.sum(result_var[:i])

            for j in range(0, i):
                if suma < 0:
                    aux = result_var[j] + ((result_var[j] / suma) * abs(var[i]))
                else:
                    aux = 0

                if aux >= 0:
                    result_var[j] = 0
                else:
                    result_var[j] = aux

    return result_var[::-1]
