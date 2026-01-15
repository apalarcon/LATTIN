import numpy as np
import pandas as pd
import os
import time
from datetime import datetime, timedelta
import operator
import functools
from lattin.lattin_functions import processing_moisture_track_backward

def DI_analysis_parms(DI_dt_thresh,
                    DI_step,
                    DI_dp_change,
                    DI_start_pressure_level,
                    DI_end_pressure_level,
                    DI_time_limits_checking,
                    DI_moisture_tracking,
                    totaltime,
                    DI_save_raw_parcels
                    ):


    errors_found=False
    errors=""

    DI_dt_thresh=DI_dt_thresh
    DI_step=DI_step
    DI_dp_change=DI_dp_change
    DI_start_pressure_level=DI_start_pressure_level
    DI_end_pressure_level=DI_end_pressure_level
    DI_time_limits_checking=DI_time_limits_checking
    DI_moisture_tracking=DI_moisture_tracking
    DI_moist_bias_correction=False
    DI_save_raw_parcels=DI_save_raw_parcels

    if DI_dt_thresh=="" or DI_dt_thresh==None:
        DI_dt_thresh=2880
    if DI_step=="" or DI_step==None:
        DI_step=1440

    if DI_dp_change=="" or DI_dp_change==None:
        DI_dp_change=40000
    if DI_start_pressure_level =="" or DI_start_pressure_level==None:
        DI_start_pressure_level=60000

    if DI_end_pressure_level =="" or DI_end_pressure_level==None:
        DI_end_pressure_level=70000

    if DI_moisture_tracking=="" or DI_moisture_tracking==None:
        DI_moisture_tracking=False

    if not isinstance(DI_time_limits_checking, list):
        DI_time_limits_checking=[DI_time_limits_checking]

    if DI_save_raw_parcels==""  or DI_save_raw_parcels==None:
        DI_save_raw_parcels=True

    if len(DI_time_limits_checking)<2:
        errors_found=True
        errors = errors + "(**) ERROR: DI_time_limits_checking should have two elements, the lower and upper limits of the backward time\n"
    else:
        if DI_time_limits_checking[0]<-1*totaltime:
            errors_found=True
            errors = errors + f"(**) ERROR: The lower limit should be higher than the backward trajectory length of {totaltime} minutes\n"
        elif DI_time_limits_checking[1]<-1*totaltime:
            errors_found=True
            errors = errors + f"(**) ERROR: The upper limit should be higher than the backward trajectory length {totaltime} minutes\n"

    return errors, errors_found, DI_dt_thresh, DI_step,DI_dp_change,DI_start_pressure_level,DI_end_pressure_level,DI_time_limits_checking,DI_moisture_tracking,DI_moist_bias_correction, DI_save_raw_parcels


def compute_grid_to_grid_DIfreq(xlons, xlats, lonmap, latmap, outgrid):
    """
    Compute a 2D grid of Dry Intrusion frequencies from a list of
    longitude and latitude points.

    Parameters
    ----------
    xlons : array_like
        Longitudes of the points
    xlats : array_like
        Latitudes of the points
    lonmap : array_like
        Longitudes of the grid points
    latmap : array_like
        Latitudes of the grid points
    outgrid : array_like
        Initial grid of zeros to be modified

    Returns
    -------
    outgrid_mod : array_like
        A 2D grid with the same shape as outgrid, where each cell has
        been incremented by one for each point that falls within that
        cell.
    """
    outgrid_mod = np.copy(outgrid)

    # Find indices for xlons and xlats that correspond to the grid cells
    lon_indices = np.searchsorted(lonmap[0, :], xlons, side='right') - 1
    lat_indices = np.searchsorted(latmap[:, 0], xlats, side='right') - 1

    # Clip indices to ensure they are within valid ranges
    lon_indices = np.clip(lon_indices, 0, lonmap.shape[1] - 2)
    lat_indices = np.clip(lat_indices, 0, latmap.shape[0] - 2)

    if not isinstance(lon_indices, list):
        lon_indices=[lon_indices]

    if not isinstance(lat_indices, list):
        lat_indices=[lat_indices]

    # Increment the corresponding grid cells in outgrid_mod
    for lat_idx, lon_idx in zip(lat_indices, lon_indices):
        outgrid_mod[lat_idx, lon_idx] += 1

    return outgrid_mod


def getting_all_DI_trajectories(DI_parcels_input):

    """
    Parameters
    ----------
    DI_parcels_input : array_like
        Array with air parcel IDs that are associated with dry intrusions

    Returns
    -------
    DIparcelsIDs : array_like
        Unique air parcel IDs associated with dry intrusions
    """
    DIparcelsIDs  = np.unique(DI_parcels_input)

    return DIparcelsIDs


def getting_DI_back_trajectories_IDs(Ptracks, dP, DI_dt_thresh, index_start, index_end,
                                     start_pressure_level, end_pressure_level,
                                     outgridDI, DIfootprint, lat, lon):
    """
    Identifies and processes DI back trajectories based on pressure levels and thresholds.

    Parameters:
    ----------
    Ptracks : np.ndarray
        3D array containing trajectory data with shape (time, parcels, variables).
        Variable at index 13 represents pressure, and index 0 represents parcel ID.
    dP : float
        Maximum allowable pressure change for DI parcels.
    DI_dt_thresh : int or float
        Time difference threshold used for identifying DI parcels (not directly used in this function).
    index_start : int
        Index representing the start of the period in Ptracks.
    index_end : int
        Index representing the end of the period in Ptracks.
    start_pressure_level : float
        Minimum pressure level for a parcel to be considered at the start of the period.
    end_pressure_level : float
        Maximum pressure level for a parcel to be considered at the end of the period.
    outgridDI : np.ndarray
        2D array representing the DI frequency at analysis time to be updated.
    DIfootprint : np.ndarray
        2D array representing the DI frequency grid to be updated.
    lat : np.ndarray
        Array of latitude coordinates for the grid.
    lon : np.ndarray
        Array of longitude coordinates for the grid.

    Returns:
    -------
    DIparcelsIDs : np.ndarray
        Array containing the IDs of DI parcels satisfying the criteria.
    outgridDI : np.ndarray
        Updated DI frequency grid.
    DIfootprint : np.ndarray
        Updated DI frequency grid.

    Notes:
    -----
    - Parcels are considered DI if they move from a pressure lower than start_pressure_level
      to a pressure higher than end_pressure_level, with a pressure change smaller than dP.
    - The function assumes Ptracks contains sufficient dimensions and data.

    Example:
    -------
    >>> Ptracks = np.random.rand(10, 100, 14)  # Example shape
    >>> outgridDI = np.zeros((180, 360))
    >>> DIfootprint = np.zeros((180, 360))
    >>> getting_DI_back_trajectories_IDs(Ptracks, dP=5, DI_dt_thresh=60,
                                         index_start=0, index_end=9,
                                         start_pressure_level=900, end_pressure_level=1000,
                                         outgridDI=outgridDI, DIfootprint=DIfootprint, lat=np.linspace(-90, 90, 180),
                                         lon=np.linspace(-180, 180, 360))
    (array([1, 5, 9]), outgridDI)
    """
    DIparcelsIDs = []

    for i in range(Ptracks.shape[1]):
        start_pressure = Ptracks[index_start, i, 13]
        end_pressure = Ptracks[index_end, i, 13]

        if (start_pressure < start_pressure_level and
            end_pressure > end_pressure_level and
            np.abs(start_pressure - end_pressure) > dP):

            DIparcelsIDs.append(Ptracks[index_start, i, 0])

            DIfootprint = compute_grid_to_grid_DIfreq(
                Ptracks[index_start:index_end + 1, i, 1],
                Ptracks[index_start:index_end + 1, i, 2],
                lon, lat, DIfootprint
            )

            outgridDI = compute_grid_to_grid_DIfreq(
                Ptracks[index_end, i, 1],
                Ptracks[index_end, i, 2],
                lon, lat, outgridDI
            )

    return np.array(DIparcelsIDs), outgridDI, DIfootprint



def filtering_DI_back_trajectories(Ptracks, parcelsDIs):
    """
    Parameters
    ----------
    Ptracks : array_like
        Array with the air parcels trajectories from the LATTIN ouputs.
    parcelsDIs : array_like
        Array with air parcel IDs associated with dry intrusions.

    Returns
    -------
    nPtracks, nDIPtracks : array_like
        Tensors with the Dry Intrusion and non-Dry Intrusion air parcels trajectories, respectively.
    """

    # # Remove duplicates and filter only valid IDs
    # valid_parcelsDIs = set(Ptracks[0, :, 0])  # Get all valid parcel IDs
    # parcelsDIs = set(parcelsDIs).intersection(valid_parcelsDIs)  # Keep only valid IDs
    #
    # # Compute correct sizes
    # n_DI_count = sum(Ptracks[0, i, 0] in parcelsDIs for i in range(Ptracks.shape[1]))
    # n_non_DI_count = Ptracks.shape[1] - n_DI_count
    #
    # # Allocate arrays with correct dimensions
    # nPtracks = np.empty((Ptracks.shape[0], n_DI_count, Ptracks.shape[2]))
    # nDIPtracks = np.empty((Ptracks.shape[0], n_non_DI_count, Ptracks.shape[2]))


    # Convert parcelsDIs to a NumPy array for efficiency
    parcelsDIs = np.array(parcelsDIs)

    # Use NumPy's `isin()` function for fast element checking
    mask_DI = np.isin(Ptracks[0, :, 0], parcelsDIs)  # Boolean mask for DI parcels

    # Compute counts efficiently
    n_DI_count = np.count_nonzero(mask_DI)  # Count of DI parcels
    n_non_DI_count = Ptracks.shape[1] - n_DI_count  # Count of non-DI parcels

    # Allocate output arrays
    nPtracks = np.empty((Ptracks.shape[0], n_DI_count, Ptracks.shape[2]))
    nDIPtracks = np.empty((Ptracks.shape[0], n_non_DI_count, Ptracks.shape[2]))



    # Fill arrays
    j = 0  # Dry Intrusion counter
    k = 0  # Non-Dry Intrusion counter
    for i in range(Ptracks.shape[1]):
        if Ptracks[0, i, 0] in parcelsDIs:
            if j < n_DI_count:
                nPtracks[:, j, :] = Ptracks[:, i, :]
                j += 1
        else:
            if k < n_non_DI_count:
                nDIPtracks[:, k, :] = Ptracks[:, i, :]
                k += 1

    return nPtracks, nDIPtracks



def get_minutes_from_start_date(start_date, track_times):

    """
    Compute the time difference in hours from a given start date
    for a list of dates.

    Parameters
    ----------
    start_date : datetime
        Starting date
    track_times : list of datetime
        List of dates

    Returns
    -------
    dt : list of float
        Time difference in minutes for each date in track_times
    """
    dt =[(itime-start_date)/timedelta(minutes=1) for itime in track_times]

    return  dt


def get_days_index(dt_mins, DI_dt_thresh, DI_step, DI_time_limits_checking):
    """
    Computes the indices for the start and end of periods based on thresholds and steps.

    Parameters:
    ----------
    dt_mins : list or np.ndarray
        List or array of time values in minutes, ordered from most recent to oldest.
    DI_dt_thresh : int or float
        Time difference threshold in minutes used to determine the start of a period.
    DI_step : int or float
        Step value in minutes for which periods are evaluated.

    DI_time_limits_checking: list:
        Interval of the backward trajectory to check DI

    Returns:
    -------
    indices_DI_end : list of int
        Indices marking the end of each period in dt_mins.
    indices_DI_start : list of int
        Indices marking the start of each period in dt_mins.

    Notes:
    -----
    - The function assumes dt_mins is sorted from most recent to oldest.
    - Periods are identified where the current value is a multiple of DI_step
      and the difference with DI_dt_thresh also exists in the list.

    Example:
    -------
    >>> dt_mins = [0, 30, 60, 90, 120, 150]
    >>> get_days_index(dt_mins, DI_dt_thresh=60, DI_step=30)
    ([5, 3], [3, 1])
    """
    dt_mins = np.array(dt_mins)

    indices_DI_end = []
    indices_DI_start = []
    DI_back_times=[]


    min_time=DI_time_limits_checking[0]
    max_time=DI_time_limits_checking[1]


    for i, dt in enumerate(dt_mins[::-1]):
        if dt % DI_step == 0 and (dt - DI_dt_thresh) in dt_mins and min_time<=dt<=max_time:
            id_end = len(dt_mins) - 1 - i
            id_start = np.where(dt_mins == dt - DI_dt_thresh)[0][0]
            indices_DI_end.append(id_end)
            indices_DI_start.append(id_start)
            DI_back_times.append(int(dt))


    return indices_DI_end[::-1], indices_DI_start[::-1],DI_back_times[::-1]



def getting_dry_intrusion_main(verbose, tensor_org, start_date, tracking_times, time_step, DI_dt_thresh, DI_step, start_pressure_level,end_pressure_level,dp_change, DI_time_limits_checking, lat_plot, lon_plot, moisture_tracking, filter_dqdt_parcels, dqdt_threshold, filter_pbl_dq_parcels,moist_custom_limits_highs, dqpblcheck, dqpbl_method, trkdq_rh_check, dqrh_threshold, mindq_gain, mindq_loss, check_RH_route_precip, precip_minrh_en_route, lag_times, area, density, dtime,moisture_linear_adjustment,  precip_minrh, lon,lat, numPdY,numPdX, cenlon, varid, moisture_tracking_method, trackingtime_steps, moist_bias_correction, precipfile, precip_var, precip_lat, precip_lon, file_mask, maskname, mask_value, maskvar_lat, maskvar_lon, rank,size, comm):

    dt_mins=get_minutes_from_start_date(start_date,  tracking_times)

    indices_DI_end, indices_DI_start,DI_back_times=get_days_index(dt_mins, DI_dt_thresh, DI_step, DI_time_limits_checking)


    matrix_DIoutgrid=np.zeros((len(DI_back_times), lat_plot.shape[0], lat_plot.shape[1]))
    matrix_DIfootprint=np.zeros((len(DI_back_times), lat_plot.shape[0], lat_plot.shape[1]))
    matrix_MUgrid=np.zeros((len(DI_back_times), len(lag_times)-1, lat_plot.shape[0], lat_plot.shape[1]))

    DI_parcels_by_step=np.zeros((len(DI_back_times)))
    precip_parcels_by_step=np.zeros((len(DI_back_times)))
    no_uptake_parcels_by_step=np.zeros((len(DI_back_times)))
    lmwvrt_by_step=np.zeros((len(DI_back_times)))

    dict_output={}
    dict_output["DI_back_times"]=DI_back_times
    mergingDI_parcels=[]
    for i, (i_start, i_end) in enumerate(zip(indices_DI_start, indices_DI_end)):

        outgridDI=np.zeros((lat_plot.shape[0], lat_plot.shape[1]))
        DIfootprint = np.zeros((lat_plot.shape[0], lat_plot.shape[1]))
        if verbose and rank==0:
            print(f"\n       » Backward time: {np.abs(dt_mins[i_end])} minutes ({np.abs(dt_mins[i_end])/1440} days)")

        parcelsDIs, outgridDI, DIfootprint =getting_DI_back_trajectories_IDs(Ptracks=tensor_org,
                                                dP=dp_change,
                                                DI_dt_thresh=DI_dt_thresh,
                                                index_start=i_start,
                                                index_end=i_end,
                                                start_pressure_level=start_pressure_level,
                                                end_pressure_level=end_pressure_level,
                                                outgridDI=outgridDI,
                                                DIfootprint=DIfootprint,
                                                lat=lat,
                                                lon=lon,
                                                )

        outgridDI=outgridDI/area
        DIfootprint=DIfootprint/area

        matrix_DIoutgrid[i,:]=outgridDI
        matrix_DIfootprint[i,:]=DIfootprint
        parcelsDIs=np.setdiff1d(parcelsDIs, mergingDI_parcels)

        mergingDI_parcels=np.append(mergingDI_parcels, parcelsDIs)



        if len(parcelsDIs)>0:

            DI_Ptracks, nDI_Ptracks = filtering_DI_back_trajectories(Ptracks=tensor_org,
                                                parcelsDIs=parcelsDIs)

            DI_parcels_by_step[i] = len(parcelsDIs)
            if verbose and rank==0:
                print(f"         - Number of trajectories with DI: {len(parcelsDIs)}", "({:.2f}".format(100 * (len(parcelsDIs)) / tensor_org.shape[1])+"%)")


            if moisture_tracking:
                array_day_moist_DI, tensor_moist, counter_precip_part, no_evap_uptakes, attributed_precip, CR, lwvrt, precip_matrix_DI, bias_cor_precip_matrix, array_moistday_corected,moistd_corrected, no_evap_uptakes_bias, lwvrt_bias = processing_moisture_track_backward(DI_Ptracks, filter_dqdt_parcels, dqdt_threshold, filter_pbl_dq_parcels,moist_custom_limits_highs, dqpblcheck, 	dqpbl_method, trkdq_rh_check, dqrh_threshold, mindq_gain, mindq_loss, check_RH_route_precip, precip_minrh_en_route, lag_times, area, density, dtime,moisture_linear_adjustment,  precip_minrh, lon,lat, numPdY,numPdX, cenlon, 1, moisture_tracking_method, trackingtime_steps, moist_bias_correction, precipfile, precip_var, precip_lat, precip_lon, file_mask, maskname, mask_value, maskvar_lat, maskvar_lon, rank,size, comm)

                matrix_MUgrid[i,:]=array_day_moist_DI
                precip_parcels_by_step[i]=counter_precip_part
                no_uptake_parcels_by_step[i]=no_evap_uptakes
                lmwvrt_by_step[i]=lwvrt
                if verbose and rank==0:
                    print("         - Computing moisture changes for air parcel trajectories with DI")
                    print("            .Number of precipitating parcels within the target region at time t0:", counter_precip_part, "({:.2f}".format(100 * (counter_precip_part) / len(parcelsDIs))+"%)")
                    print("            .Number of parcels without moisture uptake in the trajectory:", no_evap_uptakes, "({:.2f}".format(100 * (no_evap_uptakes) / (counter_precip_part))+"%)")
                    if moisture_linear_adjustment:
                        print("            .Lagrangian mean water vapour residence time: ", str(lwvrt)[0:5], "days")
        else:
            if verbose and rank==0:
                print(f"         - Number of trajectories with DI: {len(parcelsDIs)}")

    dict_output['matrix_DIoutgrid']=matrix_DIoutgrid
    dict_output['matrix_MUgrid']=matrix_MUgrid
    dict_output["DI_parcels_by_step"]=DI_parcels_by_step
    dict_output["precip_parcels_by_step"]=precip_parcels_by_step
    dict_output["no_uptake_parcels_by_step"]=no_uptake_parcels_by_step
    dict_output["lmwvrt_by_step"]=lmwvrt_by_step
    dict_output["matrix_DIfootprint"]=matrix_DIfootprint


    parcelsDIs=getting_all_DI_trajectories(mergingDI_parcels)

    if len(parcelsDIs)>0:
        DI_Ptracks, nDI_Ptracks = filtering_DI_back_trajectories(Ptracks=tensor_org,
                                        parcelsDIs=parcelsDIs)

        dict_output[f'parcels with DI']=DI_Ptracks.shape[1]
        dict_output[f'parcels without DI']=nDI_Ptracks.shape[1]
        dict_output[f"raw parcels with DI"]=DI_Ptracks
    else:
        dict_output[f'parcels with DI']=0
        dict_output[f'parcels without DI']=0

        no_DI_tmp=tensor_org.copy()
        no_DI_tmp[:,:,:]=np.nan
        dict_output[f"raw parcels with DI"]=no_DI_tmp


    if moisture_tracking and len(parcelsDIs)>0 :


        variables=["DI_Ptracks","nDI_Ptracks"]
        des=["with DI","without DI"]

        for ivar, var in enumerate(variables):

            if verbose and rank==0:
                print(f"\n       » Computing moisture changes for all air parcel trajectories {des[ivar]}")
                print(f"         - Number of trajectories {des[ivar]}: {vars()[var].shape[1]}", "({:.2f}".format(100 * (vars()[var].shape[1]) / tensor_org.shape[1])+"%)")
            array_day_moist_temp, tensor_moist, counter_precip_part, no_evap_uptakes, attributed_precip, CR, lwvrt, precip_matrix_DI, bias_cor_precip_matrix, array_moistday_corected,moistd_corrected, no_evap_uptakes_bias, lwvrt_bias = processing_moisture_track_backward(vars()[var]                                                                                                                                                                                                            , filter_dqdt_parcels, dqdt_threshold, filter_pbl_dq_parcels,moist_custom_limits_highs, dqpblcheck, dqpbl_method, trkdq_rh_check, dqrh_threshold, mindq_gain, mindq_loss, check_RH_route_precip, precip_minrh_en_route, lag_times, area, density, dtime,moisture_linear_adjustment,  precip_minrh, lon,lat, numPdY,numPdX, cenlon, 1, moisture_tracking_method, trackingtime_steps, moist_bias_correction, precipfile, precip_var, precip_lat, precip_lon, file_mask, maskname, mask_value, maskvar_lat, maskvar_lon, rank,size, comm)


            dict_output[f'matrix_MUgrid {des[ivar]}']=array_day_moist_temp
            dict_output[f"precip_parcels {des[ivar]}"]=counter_precip_part
            dict_output[f"no_uptake_parcels {des[ivar]}"]=no_evap_uptakes
            dict_output[f"lmwvrt {des[ivar]}"]=lwvrt


            if verbose and rank==0:
                print("         - Number of precipitating parcels within the target region at time t0:", counter_precip_part, "({:.2f}".format(100 * (counter_precip_part) / vars()[var].shape[1])+"%)")
                print("         - Number of parcels without moisture uptake in the trajectory:", no_evap_uptakes, "({:.2f}".format(100 * (no_evap_uptakes) / (counter_precip_part))+"%)")
                if moisture_linear_adjustment:
                    print("         - Lagrangian mean water vapour residence time: ", str(lwvrt)[0:5], "days")



    return dict_output

