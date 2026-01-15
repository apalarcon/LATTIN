import numpy as np
from lattin.lattin_functions import compute_mean_lon, is_in_pbl,compute_var_integarated_day_moist,cal_track_diff, compute_lag_mwvrt





def processing_moisture_track_forward(
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
                    file_mask,
                    maskname,
                    mask_value,
                    maskvar_lat,
                    maskvar_lon,
                    rank,
                    size,
                    comm,
                ):


    tensor_moist_ = np.copy(tensor_org[:, :, :])
    counter_precip_partdq = None


    tensor_moist = np.copy(tensor_moist_[:, :, :])
    counter_precip_part_pbl = None



    precip_matrix = np.empty((lat.shape[0] - 1, lon.shape[1] - 1))
    precip_matrix[:] = 0
    bias_cor_precip_matrix = np.empty((lat.shape[0] - 1, lon.shape[1] - 1))
    bias_cor_precip_matrix[:] = 0

    qt0 = tensor_moist[0, :, 3]

    #tensor_moist = tensor_moist[:-1, :, :]


    matrix_moist_files = ["lon", "lat", "dq", "ds", "rh", "FH"]
    dmatrix = np.empty(
        (
            len(tensor_moist[:, 0, 0]) - 2,
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

        local_results = parallel_moisture_process_forward(
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
            trackingtime_steps,
            cenlon,
            moist_bias_correction,
            check_RH_route_precip,
            precip_minrh_en_route,
            moisture_tracking_method,
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
            dmatrix[:, count * size :] = parallel_moisture_process_forward(
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
                trackingtime_steps,
                cenlon,
                moist_bias_correction,
                check_RH_route_precip,
                precip_minrh_en_route,
                moisture_tracking_method,
            )

    matrix_moist = np.copy(dmatrix[:, :, :6])

    matrix_moist[:, :, 2][matrix_moist[:, :, 2] == -999.9] = 0


    sum_prec = dmatrix[0, :, 6]
    part_uptakes = dmatrix[0, :, 7]
    evap_uptakes = dmatrix[0, :, 8]




    no_attributed_precip = np.median(dmatrix[0, :, 9][dmatrix[0, :, 9]!=-999.9])*100

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

    array_day_corrected = np.empty_like(array_day)
    moistd_corrected = np.empty_like(moistd)


    adjusted_moist_day = None
    lwvrt_bias = None

    counter_precip_part = matrix_moist.shape[1]

    no_evap_uptakes_bias = None
    CR = None


    return (
        array_day[::-1,:],
        matrix_moist,
        counter_precip_part,
        np.sum(evap_uptakes),
        None,
        CR,
        lwvrt,
        precip_matrix,
        bias_cor_precip_matrix,
        adjusted_moist_day,
        moistd_corrected,
        no_evap_uptakes_bias,
        lwvrt_bias,
    )



def parallel_moisture_process_forward(
    dmatrix,
    tensor_moist_,
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
    trackingtime_steps,
    cenlon,
    moist_bias_correction,
    check_RH_route_precip,
    precip_minrh_en_route,
    moisture_tracking_method,
):
    # uptakes_parts_precip=0
    # sum_prec=0
    # partt=0
    # noprecip=0

    trackingtime_steps = trackingtime_steps * (-1)
    #trackingtime_steps = trackingtime_steps[::-1]
    mtime = (trackingtime_steps[1:] + trackingtime_steps[:-1]) / 2

    matrix_q0 = tensor_moist_[:2,:,:]
    tensor_moist=tensor_moist_[1:,:,:].copy()
    for i in range(0, tensor_moist_.shape[1]):

        dq_source = matrix_q0[1,i,3] - matrix_q0[0,i,3]
        dq_source_cont=1

        if not moisture_tracking_method=="SJ05":
            if dq_source < mindq_gain:
               dq_source = 0
               dq_source_cont=0

        trk_check_pblq0 = is_in_pbl(
                dqpblcheck,
                matrix_q0[:, i, 4],
                matrix_q0[:, i, 7],
                dqpbl_method,
                1,
            )

        drhq0 = matrix_q0[1,i,12] - matrix_q0[0,i,12]
        trk_drhq0 = np.ones((1), dtype=bool)
        if trkdq_rh_check:

            trk_drhq0[np.abs(drhq0) <= dqrh_threshold] = True
            trk_drhq0[np.abs(drhq0) > dqrh_threshold] = False

        validq0 = trk_check_pblq0 & trk_drhq0

        if not validq0[0]:
            dq_source=0
            dq_source_cont=0


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
        dmatrix[:, i, 6] = 0

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
        check_dvar[dq <= mindq_loss] = True
        check_dvar[dq > mindq_loss] = False

        trk_drh = np.ones((len(dlon)), dtype=bool)
        if trkdq_rh_check:

            trk_drh[np.abs(drh) <= dqrh_threshold] = True
            trk_drh[np.abs(drh) > dqrh_threshold] = False

        valid = trk_check_pbl & trk_drh & check_dvar

        dmatrix[:, i, 5][valid == True] = 1
        dmatrix[:, i, 5][valid == False] = 0

        dmatrix[np.isnan(dmatrix)] = -999.9

        if moisture_linear_adjustment:

            precip_from_source,precip_from_backgroud = compute_linear_discounted_forward(dmatrix[:, i, 2], dq_source, qt0[i])

            dmatrix[:, i, 2] =  precip_from_source + precip_from_backgroud
            #print(dmatrix[:, i, 2])

            # if corrected_dq0[i] <= 0:
            #     dmatrix[:, i, 11] = 0
            #
            # else:
            #
            #     dqi = dmatrix[:, i, 2]
            #
            #     new_dqi = np.empty_like(dqi)
            #     new_dqi[:] = 0
            #
            #     dq_mask = dqi > 0
            #
            #     sumdq = np.sum(dqi[dq_mask])
            #
            #     dqfrac = dqi / sumdq
            #
            #     biasdq = sumdq - (corrected_dq0[i])
            #
            #     new_dqi[dq_mask] = dqi[dq_mask] - dqfrac[dq_mask] * biasdq
            #
            #     dmatrix[:, i, 11] = new_dqi
            #
            #     dmatrix[:, i, 13] = compute_lag_mwvrt(
            #         mtime, dmatrix[:, i, 11] * dmatrix[:, i, 5], qt0[i]
            #     )
            #
            # # print(np.sum(dmatrix[:,i,2]), np.sum(dmatrix[:,i,11]), corrected_dq0[i])
            #
            # # adj_fac = dmatrix[:,i,2] / qs[-1]
            #
            if np.sum((dmatrix[:, i, 2])) > 0:
                 dmatrix[:, i, 7] = 1
            #
            if np.sum((dmatrix[:, i, 11])) > 0:
                 dmatrix[:, i, 12] = 1
            #
            # dmatrix[:, i, 2] = dmatrix[:, i, 2]

            #
            dmatrix[:, i, 10] = compute_lag_mwvrt(
                 mtime, precip_from_source, dq_source
            )
            if dq_source>0:
                dmatrix[:,i,9]=np.sum(precip_from_backgroud)/(np.sum(dmatrix[:, i, 2]))
            else:
                dmatrix[:,i,9]=-999.9
        else:
            if len(valid[valid == True]) >= 1:
                # uptakes_parts_precip = uptakes_parts_precip + 1
                dmatrix[:, i, 7] = 1

            # noprecip+=qt0[i]
            dmatrix[:, i, 9] = 0
        dmatrix[:, i, 8] = dq_source_cont
    return dmatrix


def compute_linear_discounted_forward(dq, dq_source, q_initial):
    """
    Tracks the fate of 'dq_source' AND 'q_initial' forward in time.

    Args:
        dq: Array of moisture changes.
        dq_source: The specific evaporation amount we want to track.
        q_initial: The background moisture.

    Returns:
        precip_from_source (np.array): Negative values showing rain derived from dq_source.
        precip_from_background (np.array): Negative values showing rain derived from q_initial.
    """

    # Initialize output arrays
    precip_from_source = np.zeros_like(dq)
    precip_from_background = np.zeros_like(dq)

    # 1. Validation
    if dq_source <= 0:
        return precip_from_source, precip_from_background

    # 2. Initialize Trackers for specific components
    remaining_source_mass = dq_source
    remaining_background_mass = q_initial

    # Initialize Total Mass (Source + Background)
    # Note: New evaporation from dq[0] will be added inside the loop if dq[0]>0
    total_parcel_mass = q_initial + dq_source

    # 3. Iterate Forward (starting from T=0)
    for i in range(0, len(dq)):
        val = dq[i]

        if val == -999.9: continue

        if val > 0:
            # --- EVAPORATION (Dilution) ---
            # Add new moisture to total.
            # This dilutes both Source and Background fractions.
            total_parcel_mass += val

        elif val < 0:
            # --- PRECIPITATION (Loss) ---
            rain_amount = abs(val)

            if total_parcel_mass > 0:
                # 1. Cap rain at physical limit
                if rain_amount > total_parcel_mass:
                    rain_amount = total_parcel_mass

                # 2. Calculate Fractions for each component
                # Note: These fractions might not sum to 1.0 if there was 'New Evaporation'
                # added during the trajectory (which is neither source nor background).
                frac_source = remaining_source_mass / total_parcel_mass
                frac_background = remaining_background_mass / total_parcel_mass

                # 3. Calculate specific rain amounts
                rain_from_source = rain_amount * frac_source
                rain_from_background = rain_amount * frac_background

                # 4. Update Mass Trackers
                remaining_source_mass -= rain_from_source
                remaining_background_mass -= rain_from_background
                total_parcel_mass -= rain_amount

                # 5. Record Results
                precip_from_source[i] = -rain_from_source
                precip_from_background[i] = -rain_from_background

                # Optimization: If both tracked components are gone, we can stop
                if remaining_source_mass <= 1e-9 and remaining_background_mass <= 1e-9:
                    break
            else:
                remaining_source_mass = 0.0
                remaining_background_mass = 0.0

    return np.abs(precip_from_source), np.abs(precip_from_background)

