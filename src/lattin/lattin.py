import fnmatch
from netCDF4 import Dataset
from numpy import dtype
import matplotlib
import matplotlib.pyplot as plt
import time
from pathlib import Path as check_PATH
import math
from mpi4py import MPI
import functools
import warnings
import sys
import numpy as np
import importlib as imp
import lattin.constants as lc
from lattin.lattin_functions import *
from lattin.lattin_rp_functions import *
from lattin.lara_reanalysis import (
    get_lara_variables_from_repo,
    get_date_list_lara,
    get_lara_config,
    get_vars_from_lara_zarr,
)
from lattin.dry_intrusion_functions import (
    getting_dry_intrusion_main,
    get_minutes_from_start_date,
    DI_analysis_parms,
)

from lattin.functions_forward import processing_moisture_track_forward
from lattin.functions_flex11 import generate_flexpart11_filename, checking_raw_flex11_partposti_files, get_vars_from_partoutput
import tempfile
import shutil
import gc


def lattin_main(pathfile):
    """
    Main function for executing the LATTIN program.

    This function initializes the MPI environment, loads configuration parameters from a file, and
    performs various atmospheric tracking analyses based on the configuration. It handles different
    types of tracking such as heat, moisture, temperature anomaly, and dry intrusion analysis, and
    generates output files and statistics for each tracking mode.

    Parameters:
    pathfile (str): Path to the configuration file containing parameters for the LATTIN run.

    Raises:
    ImportError: If the configuration module cannot be imported.
    SystemExit: If there are errors in configuration parameters or required files are missing.

    Note:
    The function makes extensive use of MPI for parallel computation and assumes a suitable MPI
    environment is available.
    """

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    if rank == 0:
        start_time = time.time()
        disclaimer()
        print("Starting " + program_name() + " --->\n")
        print("Using parameters from: " + pathfile)
        print(
            "------------------------------------------------------------------------------------------------------------\n"
        )
    # content = imp.load_source("", pathfile)
    temp_file = tempfile.NamedTemporaryFile(suffix=".py", delete=False)
    shutil.copyfile(pathfile, temp_file.name)
    temp_file.close()

    spec = imp.util.spec_from_file_location("config_module", temp_file.name)
    if spec is None:
        raise ImportError(f"Cannot create a spec for the file: {pathfile}")

    content = imp.util.module_from_spec(spec)
    spec.loader.exec_module(content)

    verbose = check_paths(content, "verbose")
    runID = check_paths(content, "runID")
    # Reading paths
    raw_partposit_path = check_paths(content, "raw_partposit_path")
    basedir_lara = check_paths(content, "basedir_lara")
    outputpath = check_paths(content, "output_path")
    file_gz = check_paths(content, "file_gz")
    check_lara_files = check_paths(content, "check_lara_files")
    lara_from_https = check_paths(content, "lara_from_https")

    # Reading model details
    model = check_paths(content, "model")
    total_emited_mass = check_paths(content, "total_emited_mass")
    total_release_parcels = check_paths(content, "total_release_parcels")

    # Reading Lattin run configuration
    mode = check_paths(content, "mode")
    year = check_paths(content, "year")
    month = check_paths(content, "month")
    day = check_paths(content, "day")
    hour = check_paths(content, "hour")
    minutes = check_paths(content, "minutes")
    ndays = check_paths(content, "ndays")
    dtime = check_paths(content, "time_step")
    totaltime = check_paths(content, "tracking_time")
    calendar = check_paths(content, "calendar")
    save_full_parts_position = check_paths(content, "save_full_parts_position")

    lon_left_lower_corner = check_paths(content, "lon_left_lower_corner")
    lat_left_lower_corner = check_paths(content, "lat_left_lower_corner")
    lon_right_upper_corner = check_paths(content, "lon_right_upper_corner")
    lat_right_upper_corner = check_paths(content, "lat_right_upper_corner")

    # Reading mask details
    file_mask = check_paths(content, "file_mask")
    maskname = check_paths(content, "maskname")
    maskvar_lat = check_paths(content, "maskvar_lat")
    maskvar_lon = check_paths(content, "maskvar_lon")
    mask_value = check_paths(content, "mask_value")

    # Reading details for outpur domain
    resolution = check_paths(content, "resolution")
    numPdX = check_paths(content, "numPdX")
    numPdY = check_paths(content, "numPdY")
    lon_lower_left = check_paths(content, "lon_lower_left")
    lat_lower_left = check_paths(content, "lat_lower_left")
    area_conserving_grid = check_paths(content, "area_conserving_grid")

    # Reading details for heat tracking
    tracking_heat = check_paths(content, "tracking_heat")
    heat_tracking_method = check_paths(content, "heat_tracking_method")
    var_heat_track = check_paths(content, "var_heat_track")
    dvarheatthreshold = check_paths(content, "dvarheatthreshold")
    filter_pbl_parcels = check_paths(content, "filter_pbl_parcels")
    pblcheck = check_paths(content, "pblcheck")
    pbl_method = check_paths(content, "pbl_method")
    trk_rh_check = check_paths(content, "trk_rh_check")
    rh_threshold = check_paths(content, "rh_threshold")
    dqcheck = check_paths(content, "dqcheck")
    dqthreshold = check_paths(content, "dqthreshold")
    heat_linear_adjustment = check_paths(content, "heat_linear_adjustment")
    heat_custom_limits_highs = check_paths(content, "heat_custom_limits_highs")
    save_heat_parts_position = check_paths(content, "save_heat_parts_position")

    # Reading details for temperature anomaly decomposition
    tracking_Tanom = check_paths(content, "tracking_Tanom")
    Tanom_tracking_method = check_paths(content, "Tanom_tracking_method")

    save_Tanom_parts_position = check_paths(content, "save_Tanom_parts_position")
    path_clim_temperature = check_paths(content, "path_clim_temperature")
    climT_fname_prefix = check_paths(content, "climT_fname_prefix")
    climT_date_format = check_paths(content, "climT_date_format")
    Tlat_var_name = check_paths(content, "Tlat_var_name")
    Tlon_var_name = check_paths(content, "Tlon_var_name")
    climTvar_name = check_paths(content, "climTvar_name")
    dTdp_var_name = check_paths(content, "dTdp_var_name")
    dTdt_var_name = check_paths(content, "dTdt_var_name")
    Tplves_var_name = check_paths(content, "Tplves_var_name")
    analysis_levels = check_paths(content, "analysis_levels")
    dp_sfc = check_paths(content, "dp_sfc")
    dp_upper = check_paths(content, "dp_upper")
    Tanom_linear_adjustment = check_paths(content, "Tanom_linear_adjustment")
    psfc_var_name = check_paths(content, "psfc_var_name")
    Tanom_threshold = check_paths(content, "Tanom_threshold")
    interpolate_parcel_temperature = check_paths(
        content, "interpolate_parcel_temperature"
    )
    interpolate_sfc = check_paths(content, "interpolate_sfc")
    Tvar_name = check_paths(content, "Tvar_name")

    # Reading details for moisture tracking
    tracking_moisture = check_paths(content, "tracking_moisture")
    moisture_tracking_method = check_paths(content, "moisture_tracking_method")
    filter_dqdt_parcels = check_paths(content, "filter_dqdt_parcels")
    dqdt_threshold = check_paths(content, "dqdt_threshold")
    dqpblcheck = check_paths(content, "dqpblcheck")
    dqpbl_method = check_paths(content, "dqpbl_method")
    trkdq_rh_check = check_paths(content, "trkdq_rh_check")
    dqrh_threshold = check_paths(content, "dqrh_threshold")
    mindq_gain = check_paths(content, "mindq_gain")
    mindq_loss = check_paths(content, "mindq_loss")
    moisture_linear_adjustment = check_paths(content, "moisture_linear_adjustment")
    filter_pbl_dq_parcels = check_paths(content, "filter_pbl_dq_parcels")
    moist_custom_limits_highs = check_paths(content, "moist_custom_limits_highs")
    precip_minrh = check_paths(content, "precip_minrh")
    check_RH_route_precip = check_paths(content, "check_RH_route_precip")
    precip_minrh_en_route = check_paths(content, "precip_minrh_en_route")
    save_moist_parts_position = check_paths(content, "save_moisture_parts_position")

    ####Reading parameters for bias correcting moisture sources
    moist_bias_correction = check_paths(content, "moist_bias_correction")
    precip_fname_prefix = check_paths(content, "precip_fname_prefix")
    precip_path = check_paths(content, "precip_path")
    precip_lat = check_paths(content, "precip_lat")
    precip_lon = check_paths(content, "precip_lon")
    precip_var = check_paths(content, "precip_var")
    precip_date_format = check_paths(content, "precip_date_format")

    ####Reading parameters for dry intrusion analysi
    DI_analysis = check_paths(content, "DI_analysis")
    DI_moisture_tracking = check_paths(content, "DI_moisture_tracking")
    DI_dt_thresh = check_paths(content, "DI_dt_thresh")
    DI_step = check_paths(content, "DI_step")
    DI_dp_change = check_paths(content, "DI_dp_change")
    DI_start_pressure_level = check_paths(content, "DI_start_pressure_level")
    DI_end_pressure_level = check_paths(content, "DI_end_pressure_level")
    DI_time_limits_checking = check_paths(content, "DI_time_limits_checking")
    DI_save_raw_parcels = check_paths(content, "DI_save_raw_parcels")

    ##########################

    # lara_from_https=False

    (
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
    ) = check_init_parms(
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
    )
    if not isinstance(analysis_levels, list):
        analysis_levels = [analysis_levels]
    if not isinstance(Tanom_threshold, list):
        Tanom_threshold = [Tanom_threshold]

    if lara_from_https and model=="LARA":

        if rank == 0:
            print(
                f"WARNING: Getting LARA from zarr https using {size} CPUS is not available yet\n"
            )
            time.sleep(0.5)
            print(" -> Submit the LATTIN run using 1 CPU")
            print(f" -> If you want to run LATTIN using {size} CPUs, you must:")
            print("    - Set lara_from_https=False")
            print(
                "    - Set check_lara_files=True. LARA will download the necessary files if they are not locally stored."
            )
            print("    - Define your basedir_lara")
            print(
                "============================================================================================================"
            )
            print(program_name() + " fatal error")
            print("Bye:)")
            print(
                "============================================================================================================"
            )
        raise SystemExit()

    verbose = str2boolean(verbose)
    file_gz = str2boolean(file_gz)
    save_full_parts_position = str2boolean(save_full_parts_position)

    if model == "LARA":
        _, lara_periods, total_release_parcels, total_emited_mass, _, _ = (
            get_lara_config()
        )

    filesperday = int(1440 / dtime)
    density = total_emited_mass / total_release_parcels
    lag_times = np.arange(0, int((totaltime / dtime)) + filesperday, filesperday)

    if runID == "" or runID == None:
        runID = program_name() + "_outputs/"

    # if tracking_Tanom:
    # 	area_conserving_grid=True
    if area_conserving_grid == "" or area_conserving_grid == None:
        area_conserving_grid = False
    else:
        area_conserving_grid = area_conserving_grid

    # lon_lower_left=lon_lower_left-resolution/2
    # lat_lower_left=lat_lower_left - resolution/2

    lat, lon, cenlon = grid_point(
        resolution, numPdX, numPdY, lon_lower_left, lat_lower_left
    )

    if area_conserving_grid:
        lat, lon, cenlon = get_area_conserving_grid(
            lon_lower_left, lon.max(), lat_lower_left, lat.max(), resolution
        )
        numPdX = lon.shape[1] - 1
        numPdY = lon.shape[0] - 1

    area = calc_A(resolution, lat, lon)
    lat_plot, lon_plot = grid_plot_final(lat, lon)

    if calendar != "365d":
        calendar = "366d"

    if tracking_heat:
        (
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
        ) = heat_tracking_parms(
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
        )

        if errors_found:
            if rank == 0:
                print(
                    "Checking parameters for HEAT TRACKING: Errors detected. PLEASE TAKE ACTION!!!!\n"
                )
                time.sleep(0.5)
                print(errors)
                print(
                    "============================================================================================================"
                )
                print(program_name() + " fatal error")
                print("Bye:)")
                print(
                    "============================================================================================================"
                )
            raise SystemExit()
        else:
            if rank == 0:
                print("Checking parameters for HEAT TRACKING: PASSED\n")
    else:
        filter_pbl_parcels = "False"
        trk_rh_check = "False"
        heat_linear_adjustment = "False"
        dqcheck = "False"
        save_heat_parts_position = "False"

    if tracking_Tanom:
        (
            save_Tanom_parts_position,
            analysis_levels,
            dp_sfc,
            dp_upper,
            Tanom_linear_adjustment,
            Tanom_threshold,
            climT_date_format,
            interpolate_parcel_temperature,
            interpolate_sfc,
            errors,
            errors_found,
        ) = Tanom_tracking_parms(
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
        )

        if errors_found:
            if rank == 0:
                print(
                    "Checking parameters for TEMPERATURE ANOMALY TRACKING: Errors detected. PLEASE TAKE ACTION!!!!\n"
                )
                time.sleep(0.5)
                print(errors)
                print(
                    "============================================================================================================"
                )
                print(program_name() + " fatal error")
                print("Bye:)")
                print(
                    "============================================================================================================"
                )
            raise SystemExit()
        else:
            if rank == 0:
                print("Checking parameters for TEMPERATURE ANOMALY TRACKING: PASSED\n")

    else:
        save_Tanom_parts_position = False
        Tanom_linear_adjustment = False
        interpolate_parcel_temperature = False
        interpolate_sfc = False

    if tracking_moisture or DI_moisture_tracking:
        (
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
        ) = moisture_tracking_parms(
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
        )

        if errors_found:
            if rank == 0:
                print(
                    "Checking parameters for MOISTURE TRACKING: Errors detected. PLEASE TAKE ACTION!!!!\n"
                )
                time.sleep(0.5)
                print(errors)
                print(
                    "============================================================================================================"
                )
                print(program_name() + " fatal error")
                print("Bye:)")
                print(
                    "============================================================================================================"
                )
            raise SystemExit()
        else:

            if moist_bias_correction:
                if not filter_dqdt_parcels:
                    if rank == 0:
                        print(
                            "(**) WARNING moist_bias_correction="
                            + str(moist_bias_correction)
                            + " but filter_dqdt_parcels = "
                            + str(filter_dqdt_parcels)
                            + ". Setting moist_bias_correction = False\n"
                        )
                        moist_bias_correction = False
                elif not moisture_linear_adjustment:
                    if rank == 0:
                        print(
                            "(**) WARNING moist_bias_correction="
                            + str(moist_bias_correction)
                            + " but moisture_linear_adjustment = "
                            + str(moisture_linear_adjustment)
                            + ". Setting moist_bias_correction = False\n"
                        )
                        moist_bias_correction = False
                else:
                    if rank == 0:
                        print("Checking parameters for MOISTURE TRACKING: PASSED\n")
    else:
        filter_dqdt_parcels = "False"
        trkdq_rh_check = "False"
        moisture_linear_adjustment = "False"
        save_moist_parts_position = "False"
        check_RH_route_precip = "False"

    if DI_analysis:
        (
            errors,
            errors_found,
            DI_dt_thresh,
            DI_step,
            DI_dp_change,
            DI_start_pressure_level,
            DI_end_pressure_level,
            DI_time_limits_checking,
            DI_moisture_tracking,
            DI_moist_bias_correction,
            DI_save_raw_parcels,
        ) = DI_analysis_parms(
            DI_dt_thresh,
            DI_step,
            DI_dp_change,
            DI_start_pressure_level,
            DI_end_pressure_level,
            DI_time_limits_checking,
            DI_moisture_tracking,
            totaltime,
            DI_save_raw_parcels,
        )
        if errors_found:
            if rank == 0:
                print(
                    "Checking parameters for DRY INTRUSION ANALYSIS: Errors detected. PLEASE TAKE ACTION!!!!\n"
                )
                time.sleep(0.5)
                print(errors)
                print(
                    "============================================================================================================"
                )
                print(program_name() + " fatal error")
                print("Bye:)")
                print(
                    "============================================================================================================"
                )
            raise SystemExit()
        else:
            if rank == 0:
                print("Checking parameters for DRY INTRUSION ANALYSIS: PASSED\n")

    else:
        DI_moisture_tracking = False
        DI_save_raw_parcels = False
        DI_moist_bias_correction = False

    if (
        not tracking_heat
        and not tracking_moisture
        and not tracking_Tanom
        and not DI_analysis
    ):
        if rank == 0:
            print(
                "(**) ERROR: No tracking method is defined in the input file. Set at least heat_tracking=True or  moisture_tracking=True or tracking_Tanom = True or DI_analysis=True \n"
            )
        raise SystemExit()

    filter_pbl_parcels = str2boolean(filter_pbl_parcels)
    trk_rh_check = str2boolean(trk_rh_check)
    filter_dqdt_parcels = str2boolean(filter_dqdt_parcels)
    trkdq_rh_check = str2boolean(trkdq_rh_check)
    dqcheck = str2boolean(dqcheck)
    heat_linear_adjustment = str2boolean(heat_linear_adjustment)
    moisture_linear_adjustment = str2boolean(moisture_linear_adjustment)
    moist_bias_correction = str2boolean(moist_bias_correction)
    save_heat_parts_position = str2boolean(save_heat_parts_position)
    save_moist_parts_position = str2boolean(save_moist_parts_position)
    check_RH_route_precip = str2boolean(check_RH_route_precip)
    DI_moisture_tracking = str2boolean(DI_moisture_tracking)
    DI_save_raw_parcels = str2boolean(DI_save_raw_parcels)

    errors, errors_found = checking_input_parameters(
        size,
        lag_times[-1],
        totaltime,
        dtime,
        file_mask,
        model,
        mode,
        numPdX,
        numPdY,
        resolution,
        lat_plot[:, 0],
        lon_plot[0, :],
        ndays,
        cenlon,
        DI_analysis,
        tracking_moisture,
        tracking_Tanom,
        tracking_heat,
        moisture_tracking_method,


    )
    if errors_found:
        if rank == 0:
            print(
                "Checking config parameters: Errors detected. PLEASE TAKE ACTION!!!!\n"
            )
            time.sleep(0.5)
            print(errors)
            print(
                "============================================================================================================"
            )
            print(program_name() + " fatal error")
            print("Bye:)")
            print(
                "============================================================================================================"
            )
        raise SystemExit()
    else:
        if rank == 0:
            print("Checking config parameters: PASSED\n")

    if model in ["FLEXPART"]:
        type_file = 2
    elif model in ("FLEXPART-WRF"):
        type_file = 1

    if mode == "backward":
        track_days = np.arange(-1, -1 * len(lag_times), -1)
    elif mode == "forward":
        track_days = np.arange(1, len(lag_times), 1)

    list_year, list_month, list_day, list_hour, list_min = generate_simulation_dates(
        ndays, year, month, day, hour, minutes, calendar
    )

    output_path = outputpath + "/" + runID + "/"
    if rank == 0:
        # if not os.path.exists(output_path): os.makedirs(output_path, exist_ok=True)
        os.makedirs(output_path, exist_ok=True)  # Solo el proceso 0 crea el directorio

    if rank == 0:
        printting_run_information(
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
            trkdq_rh_check,
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
        )

    if moist_bias_correction:
        pferrors_p, find_error_p = checking_precip_files(
            precip_path,
            precip_fname_prefix,
            list_year,
            list_month,
            list_day,
            list_hour,
            list_min,
            precip_date_format,
        )

        if find_error_p:
            if rank == 0:
                print(
                    "\n     Checking precipitation files for bias-correct moisture sources: Files not found. PLEASE TAKE ACTION!!!!\n"
                )
                time.sleep(0.5)
                print(pferrors_p)
                print(
                    "======================================================================================================="
                )
                print(program_name() + " fatal error")
                print("Bye:)")
                print(
                    "======================================================================================================="
                )
            raise SystemExit()
        else:
            if rank == 0:
                print(
                    "\n     Checking precipitation files for bias-correct moisture sources: PASSED\n"
                )

    if rank == 0:
        print("\n     -> " + program_name() + " is running with", size, "CPUs\n")

    lara_tmp = list([])
    for year, month, day, hour, minn in zip(
        list_year, list_month, list_day, list_hour, list_min
    ):

        date = year + "-" + month + "-" + day + " " + hour + ":" + minn + ":00"
        filename_out = "lattin_" + mode + "_" + year + month + day + hour + minn
        date_format = "%Y-%m-%d %H:%M:%S"
        start_date = datetime.strptime(date, date_format)

        stats_fame = (
            "lattin_" + mode + "_stats_" + year + month + day + hour + minn + ".dat"
        )
        if rank == 0:
            partial_start_time = time.time()
            print("")
            print("===> STARTING TRACKING FOR ", date)
            print(
                "     -------------------------------------------------------------------------------------------------------"
            )

        partpositfiles, listdates = generate_file(
            mode, dtime, totaltime, date, raw_partposit_path, file_gz, calendar
        )


        if model in ["FLEXPART", "FLEXPART-WRF"]:
            pferrors, find_error = checking_raw_partposti_files(partpositfiles)

            if find_error:
                if rank == 0:
                    print(
                        "     Checking raw partposit files: Files not found. PLEASE TAKE ACTION!!!!\n"
                    )
                    time.sleep(0.5)
                    print(pferrors)
                    print(
                        "======================================================================================================="
                    )
                    print(program_name() + " fatal error")
                    print("Bye:)")
                    print(
                        "======================================================================================================="
                    )
                raise SystemExit()
            else:
                if rank == 0:
                    print("     Checking raw partposit files: PASSED\n")

        elif model in ['FLEXPART11']:
            partpositfiles, listdates = generate_flexpart11_filename(
            mode, dtime, totaltime, date, raw_partposit_path, file_gz, calendar
            )


            pferrors, find_error = checking_raw_flex11_partposti_files(partpositfiles, listdates)

            if find_error:
                if rank == 0:
                    print(
                        "     Checking raw partoutput files: Files not found. PLEASE TAKE ACTION!!!!\n"
                    )
                    time.sleep(0.5)
                    print(pferrors)
                    print(
                        "======================================================================================================="
                    )
                    print(program_name() + " fatal error")
                    print("Bye:)")
                    print(
                        "======================================================================================================="
                    )
                raise SystemExit()
            else:
                if rank == 0:
                    print("     Checking raw partoutput files: PASSED\n")


        elif model == "LARA":
            lara_list_of_dates_case = get_date_list_lara(
                year,
                month,
                day,
                hour,
                minn,
                year,
                month,
                day,
                hour,
                minn,
                totaltime,
                dtime,
                mode=mode,
                calendar=calendar,
            )
            lara_list_of_dates_case_ = [date_ for date_ in lara_list_of_dates_case]
            #
            if len(lara_tmp) >= 1:
                #
                lara_dates_to_check = [
                    file_date
                    for file_date in lara_list_of_dates_case_
                    if file_date not in lara_tmp
                ]
            else:
                lara_dates_to_check = lara_list_of_dates_case_

                lara_tmp = np.append(lara_tmp, lara_dates_to_check)

            if len(lara_dates_to_check) >= 1:
                if rank == 0:
                    if check_lara_files:
                        print("\n     -> CHECKING LARA FILES\n")
                        get_lara_variables_from_repo(
                            basedir_lara, pd.DatetimeIndex(lara_dates_to_check)
                        )
                    else:
                        if not lara_from_https:
                            print(
                                "\n     -> WARNING: LARA local files will not be checked\n        Set check_lara_files=True to check files"
                            )

                comm.Barrier()
            else:
                if rank == 0:
                    if check_lara_files:
                        print("\n     -> CHECKING LARA FILES: PASSED\n")

        ntimes = []
        for i in range(0, len(listdates)):
            auxt = datetime(
                int(str(int(listdates[i]))[0:4]),
                int(str(int(listdates[i]))[4:6]),
                int(str(int(listdates[i]))[6:8]),
                int(str(int(listdates[i]))[8:10]),
                int(str(int(listdates[i]))[10:12]),
                int(str(int(listdates[i]))[12:14]),
            )
            ntimes = np.append(ntimes, auxt)

        meantimes = []
        for i in range(0, len(ntimes) - 1):
            meantimes = np.append(
                meantimes, ntimes[i] + (ntimes[i + 1] - ntimes[i]) / 2
            )

        trackingtime_steps = np.arange(dtime, totaltime + 2 * dtime, dtime)

        climT_filenames = generate_climT_filenames(
            climT_fname_prefix, climT_date_format, ntimes
        )
        if tracking_Tanom:
            pferrors_T, find_error_T = checking_Tclimatological_files(
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
            )

            if find_error_T:
                if rank == 0:
                    print(
                        "     Checking climatological 3D temperature files temperature anomaly analysis: Files not found. PLEASE TAKE ACTION!!!!\n"
                    )
                    time.sleep(0.5)
                    print(pferrors_T)
                    print(
                        "======================================================================================================="
                    )
                    print(program_name() + " fatal error")
                    print("Bye:)")
                    print(
                        "======================================================================================================="
                    )
                raise SystemExit()
            else:
                if rank == 0:
                    print(
                        "\n     Checking climatological 3D temperature files for temperature anomaly analysis: PASSED\n"
                    )

        if model in ["FLEXPART", "FLEXPART-WRF"]:
            if rank == 0:
                print("\n     !Getting data from raw partposit files")
            tensor_org = get_vars_from_partposit(
                verbose,
                partpositfiles,
                file_mask,
                maskname,
                maskvar_lon,
                maskvar_lat,
                lat,
                lon,
                rank,
                size,
                comm,
                type_file,
                lon_left_lower_corner,
                lat_left_lower_corner,
                lon_right_upper_corner,
                lat_right_upper_corner,
                model,
                mask_value,
                file_gz,
                var_heat_track,
                mode,

            )
        elif model in ["FLEXPART11"]:
            if rank == 0:
                print("\n     !Getting data from raw partoutput files")
            tensor_org = get_vars_from_partoutput(
                verbose,
                partpositfiles,
                file_mask,
                maskname,
                maskvar_lon,
                maskvar_lat,
                lat,
                lon,
                rank,
                size,
                comm,
                lon_left_lower_corner,
                lat_left_lower_corner,
                lon_right_upper_corner,
                lat_right_upper_corner,
                model,
                mask_value,
                var_heat_track,
                mode,
                listdates,
            )
            comm.Barrier()
            #quit()
        elif model == "LARA":
            if rank == 0:
                if lara_from_https:
                    print("\n     !Getting data from LARA zarr https")
                else:
                    print("\n     !Getting data from LARA local files")

            tensor_org = get_vars_from_lara_zarr(
                 verbose,
                 lara_list_of_dates_case,
                 basedir_lara,
                 file_mask,
                 maskname,
                 maskvar_lon,
                 maskvar_lat,
                 lat,
                 lon,
                 rank,
                 size,
                 comm,
                 lon_left_lower_corner,
                 lat_left_lower_corner,
                 lon_right_upper_corner,
                 lat_right_upper_corner,
                 mask_value,
                 var_heat_track,
                 start_date,
                 from_https=lara_from_https,
                 mode=mode,
            )

        # if rank==0:
        #
        # 	np.save(f"tensor_flex_test.npy", tensor_org) #file_mask, maskvar_lon, maskvar_lat, mask_value, maskname
        # quit()
        # print(tensor_tx.shape)
        # quit()
        # tensor_org=np.load("test_new_tensor.npy")
        # testfilename =  f"/mnt/lustre/hsm/nlsas/notape/home/uvi/fi/exmerisk/albenis/Cuba_HW_May24/TEST_NWP_HW21/test_orgiginal_code/lagranto_arrays_org_mod/10/traj_{start_date}.npy"
        # tensor_org=np.load(testfilename)
        # print(testfilename)

        parcels_count = len(tensor_org[0, :, 0])

        precipfile = f"{precip_path}/{precip_fname_prefix}_{start_date.strftime(precip_date_format)}.nc"

        if rank == 0:
            print(
                "\n     => Number of parcels within the target region:", parcels_count
            )

        # Processing heat tracking #######################################
        if tracking_heat:
            if verbose and rank == 0:
                print("\n     + PROCESSING HEAT")
            if mode == "backward":

                array_day_heat, tensor_heat, counter_part_heat, no_heatuptakes_parts = (
                    processing_heat_track_backward(
                        tensor_org,
                        pblcheck,
                        filter_pbl_parcels,
                        pbl_method,
                        heat_custom_limits_highs,
                        trk_rh_check,
                        rh_threshold,
                        dqcheck,
                        dqthreshold,
                        dvarheatthreshold,
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
                        0,
                        rank,
                        size,
                        comm,
                    )
                )
                case_heat = "heat uptake"

                if filter_pbl_parcels:
                    if verbose and rank == 0:
                        print(
                            "      !Number of filtered parcels within the PBL or custom highs at time t0:",
                            counter_part_heat,
                            "({:.2f}".format(
                                100 * (counter_part_heat) / (parcels_count)
                            )
                            + "%)",
                        )
                        print(
                            f"      !Number of parcels without {case_heat} in the trajectory:",
                            no_heatuptakes_parts,
                            "({:.2f}".format(
                                100 * (no_heatuptakes_parts) / (counter_part_heat)
                            )
                            + "%)",
                        )

                    save_heat_stats = True
                else:
                    save_heat_stats = True
                    if rank == 0:
                        print(
                            f"      !Number of parcels without {case_heat} in the trajectory:",
                            no_heatuptakes_parts,
                            "({:.2f}".format(
                                100 * (no_heatuptakes_parts) / (counter_part_heat)
                            )
                            + "%)",
                        )

        else:
            case_heat = "None"
            save_heat_stats = False
            counter_part_heat = None
            no_heatuptakes_parts = None
            array_day_heat = np.empty((len(track_days), lat.shape[0], lat.shape[1]))
            tensor_heat = np.empty_like(tensor_org[1:, :, :])

        if tracking_Tanom:
            if verbose and rank == 0:
                print("\n     + PROCESSING TEMPERATURE ANOMALY")
            if mode == "backward":

                # interpolating climatological variables to air parcel trajectories

                gc.collect()
                airtrajs_org = interpolate_to_parcel_trajectories(
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
                    comm,
                    rank,
                    size,
                    start_date,
                    ntimes,
                    tensor_org,
                    interpolate_parcel_temperature,
                    interpolate_sfc,
                )
                gc.collect()

                # airtrajs_org=np.load(f"/mnt/lustre/hsm/nlsas/notape/home/uvi/fi/exmerisk/albenis/Cuba_HW_May24/TEST_NWP_HW21/test_orgiginal_code/lagranto_arrays_org/10/traj_{start_date}.npy")

                # np.save("interpolated_airtrajs.npy", airtrajs_org)
                # quit()
                # airtrajs_org=np.load("interpolated_airtrajs.npy")

                # quit()
                case_Tanom = "Temperature anomaly contribution"
                if Tanom_tracking_method == "RP23":
                    dict_output = gettting_temperature_anom_budget(
                        airtrajs_org,
                        cenlon,
                        dtime,
                        lat,
                        lon,
                        analysis_levels,
                        dp_sfc,
                        dp_upper,
                        Tanom_threshold,
                        ntimes,
                        start_date,
                        numPdX,
                        numPdY,
                        rank,
                    )

                    if verbose and rank == 0:
                        for level in analysis_levels:
                            print("      - ANALYSIS LEVEL:", level)
                            print(
                                "      !Number of filter parcels at time t0 :",
                                dict_output[f"{level}_counter_part_Tanom"],
                                "({:.2f}".format(
                                    100
                                    * (dict_output[f"{level}_counter_part_Tanom"])
                                    / (parcels_count)
                                )
                                + "%)",
                            )
                            print(
                                f"      !Number of parcels without {case_Tanom}:",
                                dict_output[f"{level}_no_Tanom_parts"],
                                "({:.2f}".format(
                                    100
                                    * (dict_output[f"{level}_no_Tanom_parts"])
                                    / ((dict_output[f"{level}_counter_part_Tanom"]))
                                )
                                + "%)\n",
                            )

                    save_Tanom_stats = True
                elif Tanom_tracking_method == "PR23":
                    dict_output = gettting_temperature_anom_budget_sources(
                        airtrajs_org,
                        cenlon,
                        dtime,
                        lat,
                        lon,
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
                    )
                    if verbose and rank == 0:
                        for level in analysis_levels:
                            print("\n      - ANALYSIS LEVEL:", level)
                            print(
                                "      !Number of filter parcels at time t0 within the level:",
                                dict_output[f"{level}_filtered_parcels_within_level"],
                                "({:.2f}".format(
                                    100
                                    * (
                                        dict_output[
                                            f"{level}_filtered_parcels_within_level"
                                        ]
                                    )
                                    / (parcels_count)
                                )
                                + "%)",
                            )
                            print(
                                "      !Number of filter parcels at time t0 after apply temperature anomaly threshold:",
                                dict_output[f"{level}_filtered_parcels"],
                                "({:.2f}".format(
                                    100
                                    * (dict_output[f"{level}_filtered_parcels"])
                                    / (
                                        dict_output[
                                            f"{level}_filtered_parcels_within_level"
                                        ]
                                    )
                                )
                                + "%)",
                            )
                            print(
                                f"      !Number of parcels without Tanom contribution to temperature anonaly:",
                                dict_output[f"{level}_no_part_tanom"],
                                "({:.2f}".format(
                                    100
                                    * (dict_output[f"{level}_no_part_tanom"])
                                    / ((dict_output[f"{level}_filtered_parcels"]))
                                )
                                + "%)",
                            )
                            print(
                                f"      !Number of parcels without seas contribution to temperature anonaly:",
                                dict_output[f"{level}_no_part_seas"],
                                "({:.2f}".format(
                                    100
                                    * (dict_output[f"{level}_no_part_seas"])
                                    / ((dict_output[f"{level}_filtered_parcels"]))
                                )
                                + "%)",
                            )
                            print(
                                f"      !Number of parcels without adiab contribution to temperature anonaly:",
                                dict_output[f"{level}_no_part_adiab"],
                                "({:.2f}".format(
                                    100
                                    * (dict_output[f"{level}_no_part_adiab"])
                                    / ((dict_output[f"{level}_filtered_parcels"]))
                                )
                                + "%)",
                            )
                            print(
                                f"      !Number of parcels without adv contribution to temperature anonaly:",
                                dict_output[f"{level}_no_part_adv"],
                                "({:.2f}".format(
                                    100
                                    * (dict_output[f"{level}_no_part_adv"])
                                    / ((dict_output[f"{level}_filtered_parcels"]))
                                )
                                + "%)",
                            )
                            print(
                                f"      !Number of parcels without diab contribution to temperature anonaly:",
                                dict_output[f"{level}_no_part_diab"],
                                "({:.2f}".format(
                                    100
                                    * (dict_output[f"{level}_no_part_diab"])
                                    / ((dict_output[f"{level}_filtered_parcels"]))
                                )
                                + "%)\n",
                            )
                            print(f"      !Contributions (% T')")
                            print(f"      Term      Overall    Positive     Negative")
                            print(f"      ------------------------------------------")
                            for iterm, term in enumerate(
                                dict_output[f"{level}_tanom_terms"]
                            ):
                                print(
                                    f'      {term.ljust(5)}     {dict_output[f"{level}_overall_contributions"][iterm]:.2f}    {dict_output[f"{level}_positive_contributions"][iterm]:.2f}     {dict_output[f"{level}_negative_contributions"][iterm]:.2f}'
                                )

                    save_Tanom_stats = True
                else:
                    print("No tracking method available")

        else:
            dict_output = {}
            case_Tanom = "None"
            save_Tanom_stats = False

        # Processing moisture tracking #######################################
        if tracking_moisture:

            if verbose and rank == 0:
                print("\n     + PROCESSING MOISTURE")
            if mode == "backward":
                (
                    array_day_moist,
                    tensor_moist,
                    counter_precip_part,
                    no_evap_uptakes,
                    attributed_precip,
                    CR,
                    lwvrt,
                    precip_matrix,
                    bias_cor_precip_matrix,
                    array_moistday_corected,
                    moistd_corrected,
                    no_evap_uptakes_bias,
                    lwvrt_bias,
                ) = processing_moisture_track_backward(
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
                    1,
                    moisture_tracking_method,
                    trackingtime_steps,
                    moist_bias_correction,
                    precipfile,
                    precip_var,
                    precip_lat,
                    precip_lon,
                    file_mask,
                    maskname,
                    mask_value,
                    maskvar_lat,
                    maskvar_lon,
                    rank,
                    size,
                    comm,
                )

                if filter_dqdt_parcels == True or filter_pbl_dq_parcels == True:
                    if verbose and rank == 0:
                        print(
                            "      !Number of precipitating parcels within the target region at time t0:",
                            counter_precip_part,
                            "({:.2f}".format(
                                100 * (counter_precip_part) / (parcels_count)
                            )
                            + "%)",
                        )

                        print(
                            "      !Number of parcels without moisture uptake in the trajectory:",
                            no_evap_uptakes,
                            "({:.2f}".format(
                                100 * (no_evap_uptakes) / (counter_precip_part)
                            )
                            + "%)",
                        )

                        if moist_bias_correction:
                            print(
                                "      !Number of parcels without moisture contribution to precipitation in the target region after the bias correction:",
                                no_evap_uptakes_bias,
                                "({:.2f}".format(
                                    100 * (no_evap_uptakes_bias) / (counter_precip_part)
                                )
                                + "%)",
                            )

                        if moisture_linear_adjustment:
                            print(
                                "      !Lagrangian mean water vapour residence time: ",
                                str(lwvrt)[0:5],
                                "days",
                            )

                        if moist_bias_correction:
                            print(
                                "      !Lagrangian mean water vapour residence time after the bias correction: ",
                                str(lwvrt_bias)[0:5],
                                "days",
                            )

                    save_moist_stats = True
                else:
                    save_moist_stats = True
                    if verbose and rank == 0:

                        print(
                            "      !Number of precipitating parcels within the target region at time t0:",
                            counter_precip_part,
                            "({:.2f}".format(
                                100 * (counter_precip_part) / (parcels_count)
                            )
                            + "%)",
                        )


            elif mode=="forward":
                tensor_org=tensor_org[::-1,:]
                (
                    array_day_moist,
                    tensor_moist,
                    counter_precip_part,
                    no_evap_uptakes,
                    attributed_precip,
                    CR,
                    lwvrt,
                    precip_matrix,
                    bias_cor_precip_matrix,
                    array_moistday_corected,
                    moistd_corrected,
                    no_evap_uptakes_bias,
                    lwvrt_bias,
                )= processing_moisture_track_forward(
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
                    1,
                    moisture_tracking_method,
                    trackingtime_steps,
                    moist_bias_correction,
                    precipfile,
                    precip_var,
                    precip_lat,
                    precip_lon,
                    file_mask,
                    maskname,
                    mask_value,
                    maskvar_lat,
                    maskvar_lon,
                    rank,
                    size,
                    comm,
                )
                save_moist_stats = True
                if verbose and rank==0:
                    print(
                            "      !Number of parcels tracked forward after filtering in the source region:",
                            no_evap_uptakes,
                            "({:.2f}".format(
                                100 * (no_evap_uptakes) / (parcels_count)
                            )
                            + "%)",
                        )


                if moisture_linear_adjustment:
                    if verbose and rank == 0:
                        print(
                            "      !Lagrangian mean water vapour residence time: ",
                            str(lwvrt)[0:5],
                            "days",
                        )


        else:
            array_day_moist = np.empty((len(track_days), lat.shape[0], lat.shape[1]))
            tensor_moist = np.empty_like(tensor_org[:-1, :, :])
            CR = np.empty((len(track_days), lat.shape[0], lat.shape[1]))
            save_moist_stats = False
            save_attrib = False
            counter_precip_part = None
            no_evap_uptakes = None
            lwvrt = None
            precip_matrix = None
            bias_cor_precip_matrix = None
            array_moistday_corected = None
            moistd_corrected = None
            no_evap_uptakes_bias = None
            lwvrt_bias = None

        if DI_analysis:
            if verbose and rank == 0:
                print("\n     + PROCESSING DRY INTRUSION")
            if mode == "backward":
                DI_dict_output = getting_dry_intrusion_main(
                    verbose,
                    tensor_org,
                    start_date,
                    ntimes,
                    dtime,
                    DI_dt_thresh,
                    DI_step,
                    DI_start_pressure_level,
                    DI_end_pressure_level,
                    DI_dp_change,
                    DI_time_limits_checking,
                    lat_plot,
                    lon_plot,
                    DI_moisture_tracking,
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
                    1,
                    moisture_tracking_method,
                    trackingtime_steps,
                    DI_moist_bias_correction,
                    precipfile,
                    precip_var,
                    precip_lat,
                    precip_lon,
                    file_mask,
                    maskname,
                    mask_value,
                    maskvar_lat,
                    maskvar_lon,
                    rank,
                    size,
                    comm,
                )

                save_DI_stats = True
        else:
            save_DI_stats = False
            DI_dict_output = {}

        if rank == 0:
            if verbose:
                print("\n     + SAVING TO")
                print("       !Output File:", output_path + "/" + filename_out + ".nc")
            writing_netcdf(
                latitude=lat_plot[:, 0],
                longitude=lon_plot[0, :],
                tracking_heat=tracking_heat,
                var_heat_track=var_heat_track,
                array_heat_day=array_day_heat,
                tracking_moisture=tracking_moisture,
                array_moist_day=array_day_moist,
                tensor_org=tensor_org,
                tensor_heat=tensor_heat,
                tensor_moist=tensor_moist,
                CR=CR,
                area=area,
                track_days=track_days,
                mode=mode,
                run_date=date,
                listdates=ntimes,
                meantimes=meantimes,
                dtime=dtime,
                moisture_tracking_method=moisture_tracking_method,
                heat_tracking_method=heat_tracking_method,
                moisture_linear_adjustment=moisture_linear_adjustment,
                filename=output_path + "/" + filename_out,
                save_heat_parts_position=save_heat_parts_position,
                save_moist_parts_position=save_moist_parts_position,
                save_full_parts_position=save_full_parts_position,
                precip_matrix=precip_matrix,
                filter_dqdt_parcels=filter_dqdt_parcels,
                bias_cor_precip_matrix=bias_cor_precip_matrix,
                array_moistday_corected=array_moistday_corected,
                moist_bias_correction=moist_bias_correction,
                save_Tanom_parts_position=save_Tanom_parts_position,
                Tanom_tracking_method=Tanom_tracking_method,
                tracking_Tanom=tracking_Tanom,
                analysis_levels=analysis_levels,
                dict_output=dict_output,
                DI_analysis=DI_analysis,
                DI_dict_output=DI_dict_output,
                DI_moisture_tracking=DI_moisture_tracking,
                DI_save_raw_parcels=DI_save_raw_parcels,
            )

            if verbose:
                print("\n     + SAVING STATS TO")
                print("       !Stats File:", output_path + "/" + stats_fame)
            partial_runtime = time.time() - partial_start_time
            saving_data(
                parcels_count,
                counter_part_heat,
                no_heatuptakes_parts,
                counter_precip_part,
                no_evap_uptakes,
                date,
                output_path + "/" + stats_fame,
                save_moist_stats,
                save_heat_stats,
                np.round(partial_runtime, 2),
                mode,
                heat_tracking_method,
                moisture_tracking_method,
                lwvrt,
                moisture_linear_adjustment,
                moist_bias_correction,
                no_evap_uptakes_bias,
                lwvrt_bias,
                case_heat,
                case_Tanom,
                Tanom_tracking_method,
                save_Tanom_stats,
                dict_output=dict_output,
                analysis_levels=analysis_levels,
                DI_analysis=DI_analysis,
                DI_dict_output=DI_dict_output,
                DI_moisture_tracking=DI_moisture_tracking,
                save_DI_stats=save_DI_stats,
            )

            if verbose:
                print(
                    "\n     --------------------------------------------------------- "
                )
                print(
                    "     "
                    + program_name()
                    + " Version "
                    + str(get_currentversion())
                    + " has successfully finished this date"
                )
                print(
                    "     Partial Runtime: %.2f seconds." % np.round(partial_runtime, 2)
                )
                print("     --------------------------------------------------------- ")
        gc.collect()
        comm.Barrier()
        time.sleep(0.05)
    if rank == 0:

        runtime = time.time() - start_time
        ending_credits(runtime)
    MPI.Finalize()
