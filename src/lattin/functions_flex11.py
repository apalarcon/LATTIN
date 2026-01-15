import numpy as np
from lattin.lattin_functions import (
    time_calcminutes,
    compute_dry_static_energy,
    compute_theta,
    compute_rh,
    calc_pres,
    load_mask_grid_NR,
    generate_file,
    funtion_interpol_mascara,
    determine_id_target_region)
import xarray as xr
import pandas as pd
import scipy.interpolate as interp
import gc
def generate_flexpart11_filename(mode, dtime, totaltime, fecha, path, key_gz, calendar):
    #partoutput_20210608180000.nc
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
            name = path + "partoutput_" + fecha_dia+".nc"
            if key_gz:
                name = path + "partoutput_" + fecha_dia + ".gz"
            else:
                name = path + "partoutput_" + fecha_dia+".nc"

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
            name = path + "partoutput_" + fecha_dia + ".gz"
        else:
            name = path + "partoutput_" + fecha_dia+".nc"

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
            name = path + "partoutput_" + fecha_dia+".nc"
            if key_gz:
                name = path + "partoutput_" + fecha_dia + ".gz"
            else:
                name = path + "partoutput_" + fecha_dia+".nc"


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



def checking_raw_flex11_partposti_files(partpositfiles, listdates):
    errors = ""
    find_error = False
    for i, fdate in zip(partpositfiles, listdates):

        ts_fdate = ts = pd.to_datetime(str(int(fdate)), format='%Y%m%d%H%M%S')

        #try:
        ds = xr.open_dataset(i)

        timeds = ds['time'].values
        ds.close()

        if ts_fdate in timeds:
            pass
        else:
             errors = errors + "     (**) ERROR: No such file or directory: " + i + "\n"
             find_error = True

    return errors, find_error


def interpol_gridded_to_parcels(grid_lat, grid_lon, grid_variable, parcels):

    """
    Interpolates a gridded field to air parcel positions

    Parameters
    ----------
    grid_lat : 1D array
        Latitudes of the grid points
    grid_lon : 1D array
        Longitudes of the grid points
    grid_variable : 2D array
        Gridded field to be interpolated
    parcels : 2D array
        Array of shape (N, 3) with columns [lat, lon, pressure]

    Returns
    -------
    1D array
        Interpolated values at the parcel positions
    """

    grid_lon = (grid_lon + 180) % 360 - 180

    lon_, lat_ = np.meshgrid(grid_lon, grid_lat)

    # Direct array stacking - eliminates intermediate array creation
    flat_lon = lon_.flatten()
    flat_lat = lat_.flatten()
    flat_var = grid_variable.flatten()

    # 2. Create a mask: Keep only points where data AND coordinates are finite (not NaN)
    # This effectively ignores land masks or missing data points
    mask = np.isfinite(flat_var) & np.isfinite(flat_lon) & np.isfinite(flat_lat)

    # 3. Apply the mask to create clean inputs
    # We only stack the coordinates that have valid data
    lat_lon_clean = np.column_stack((flat_lon[mask], flat_lat[mask]))
    var_clean = flat_var[mask]

    # 4. Check if we have data left (prevent crash if everything was NaN)
    if len(var_clean) == 0:
        raise ValueError("All data in grid_variable or coordinates are NaN/Inf.")

    # 5. Create interpolator using only valid points
    prsInterpu = interp.NearestNDInterpolator(lat_lon_clean, var_clean)

    return prsInterpu(parcels)


def get_array_indices(variable):

    """
    Return the index of the given variable in the array

    Parameters
    ----------
    variable: str
        Name of the variable, e.g. "parcel", "lon", "sh", etc.

    Returns
    -------
    int
        Index of the given variable in the array

    Raises
    ------
    KeyError
        If the given variable is not in the dictionary
    """
    array_indices = {
        "parcel": 0,
        "lon": 1,
        "lat": 2,
        "sh": 3,
        "z": 4,
        "to": 5,
        "rho": 6,
        "hmix": 7,
        "tro": 8,
        "T": 9,
        "m": 10,
        "pottemp": 11,
        "prs": 13,
    }

    return array_indices[variable]




def get_parcels_over_target_region2(
    variables, name_file, file_mask, maskname, maskvar_lon, maskvar_lat, mask_value,
):

    """
    Get parcels over target region.

    Parameters
    ----------
    variables : list of str
        list of variables to read from LARA
    name_file : str
        path to the filename containing start date
    file_mask : str
        path to mask file
    maskname : str
        name of mask variable in mask file
    maskvar_lat : str
        latitude variable name in mask file
    maskvar_lon : str
        longitude variable name in mask file
    mask_value : int
        mask value for filtering parcels in the target region

    Returns
    -------
    filter_matrix : 2D array
        2D array with the parcel IDs in the target region
    """

    ds = xr.open_dataset(name_file)
    parcels = ds["particle"].values
    for var_name in variables:



        data_values = ds[var_name].values

        if var_name == "lon":
            data_values = (data_values + 180) % 360 - 180

        vars()[var_name] = data_values




    ds.close()

    gc.collect()

    data_temp = np.column_stack(
        (parcels, vars()["lon"], vars()["lat"])
    )

    mask = ~np.isnan(data_temp[:, 1]) & ~np.isnan(data_temp[:, 2])

    # Apply the mask to keep only the valid rows
    data_temp = data_temp[mask]

    lat_masked, lon_masked, mascara = load_mask_grid_NR(
        file_mask, maskname, maskvar_lon, maskvar_lat
    )

    filter_matrix = determine_id_target_region(
        data_temp,
        lat_masked.flatten(),
        lon_masked.flatten(),
        mascara.flatten(),
        mask_value,
    )

    del lat_masked, lon_masked, mascara
    gc.collect()

    return filter_matrix[:, 0]

def get_parcels_over_target_region(variables,
    name_file, file_mask, maskname, maskvar_lon, maskvar_lat, mask_value,fdate,
):
    """
    Optimized function to get parcels over target region.
    Removed 'variables' list because only lon/lat/id are needed for filtering.
    """

    # --- 1. Load Particle Data (Only what is needed) ---
    # Using 'with' ensures the file is closed automatically
    ts_fdate = ts = pd.to_datetime(str(int(fdate)), format='%Y%m%d%H%M%S')
    with xr.open_dataset(name_file, chunks={}) as ds:
        # Load directly into numpy arrays.

        ds_t = ds.sel(time=ts_fdate)


        part_id = ds["particle"].values.flatten()
        part_lon = ds["lon"].values.flatten()
        part_lat = ds["lat"].values.flatten()

    # --- 2. Normalize Longitude ---
    # In-place operation is slightly faster and saves memory
    part_lon = (part_lon + 180) % 360 - 180

    # --- 3. Filter NaNs (Vectorized) ---
    # Create a boolean mask for valid coordinates
    valid_mask = np.isfinite(part_lon) & np.isfinite(part_lat)

    # Apply filter immediately to reduce array size
    part_id = part_id[valid_mask]
    part_lon = part_lon[valid_mask]
    part_lat = part_lat[valid_mask]

    data_temp_ = np.column_stack(
        (part_id, part_lon, part_lat)
    )

    mask = ~np.isnan(data_temp_[:, 1]) & ~np.isnan(data_temp_[:, 2])

    # Apply the mask to keep only the valid rows
    data_temp = data_temp_[mask]

    lat_masked, lon_masked, mascara = load_mask_grid_NR(
        file_mask, maskname, maskvar_lon, maskvar_lat
    )

    filter_matrix = determine_id_target_region(
        data_temp,
        lat_masked.flatten(),
        lon_masked.flatten(),
        mascara.flatten(),
        mask_value,
    )

    del lat_masked, lon_masked, mascara
    gc.collect()

    return filter_matrix[:, 0]

def read_by_proccesor(
    verbose,
    partpositfiles,
    local_dates,
    parcel_ids,
    submatrix,
    rank,
    lon_left_lower_corner,
    lat_left_lower_corner,
    lon_right_upper_corner,
    lat_right_upper_corner,
    model,
    ):

    a1 = np.arange(len(partpositfiles))
    dx, dy = submatrix.shape
    tensor_local = np.ones((len(partpositfiles), dx, dy)) * (-999.9)
    for i in a1:
        if verbose:
            print("Reading | " + model + " -> ", partpositfiles[i])
        matrix_i = read_flex11_file(partpositfiles[i],
                        parcel_ids,
                        local_dates[i],
                        lon_left_lower_corner,
                        lat_left_lower_corner,
                        lon_right_upper_corner,
                        lat_right_upper_corner,)

        tensor_local[i, :, :] = matrix_i
    return tensor_local



def read_flex11_file(name_file,
            parcel_ids,
            fdate,
            lon_left_lower_corner,
            lat_left_lower_corner,
            lon_right_upper_corner,
            lat_right_upper_corner):

    ds = xr.open_dataset(name_file, chunks={})
    variables = ["lat", "lon", "T", "z", "sh", "rho"]
    ds_subset_ = ds[variables]

    ts_fdate = pd.to_datetime(str(int(fdate)), format='%Y%m%d%H%M%S')
    ds_subset=ds_subset_.sel(time=ts_fdate)

    sort_indices = np.argsort(parcel_ids)

    parcel_ids = parcel_ids.astype(int)
    target_ds = ds_subset.sel(particle=parcel_ids).load()


    gridded_variables = ["hmix"]
    tmp_tensor = np.full((len(parcel_ids), 11), -999.0)
    tmp_tensor[:, 0]=parcel_ids[sort_indices]


    for var_name in variables:

        indice = get_array_indices(var_name)


        # Load raw data

        raw_data = target_ds[var_name].values

        # Apply longitude correction if needed
        if var_name == "lon":
            raw_data = (raw_data + 180) % 360 - 180


        # Apply the pre-calculated sort order (Fast indexing)
        sorted_data = raw_data[sort_indices]

        # Assign to tensor
        tmp_tensor[:, indice] = sorted_data

    target_ds.close()

    target_ds_ = ds[gridded_variables]
    target_ds = target_ds_.sel(time=ts_fdate)
    for var_name in gridded_variables:
        indice = get_array_indices(var_name)
        aux_data = target_ds[var_name].values[:,:]



        values_data = interpol_gridded_to_parcels(
                            target_ds["latitude"].values,
                            target_ds["longitude"].values,
                            aux_data,
                            tmp_tensor[:, 1:3],
                        )
        tmp_tensor[:, indice] = values_data


    ds.close()
    if (
        lon_left_lower_corner > -180
        or lon_right_upper_corner < 180
        or lat_left_lower_corner > -90
        or lat_right_upper_corner < 90
    ):

        # Check for out of bounds coordinates
        lon_out_of_bounds = (tmp_tensor[:, 1] <= lon_left_lower_corner) | (
            tmp_tensor[:, 1] >= lon_right_upper_corner
        )
        lat_out_of_bounds = (tmp_tensor[:, 2] <= lat_left_lower_corner) | (
            tmp_tensor[:, 2] >= lat_right_upper_corner
        )

        # Combine masks for points that are out of bounds
        bad_points_mask = lon_out_of_bounds | lat_out_of_bounds

        # Set all elements from index 1 onwards to -999 for out-of-bounds points
        tmp_tensor[bad_points_mask, 1:] = -999.0

    rows_to_keep = np.any(tmp_tensor[:, 0] != -999.0, axis=0)
    tmp_tensor = tmp_tensor[rows_to_keep]


    return tmp_tensor[0,:,:]

def get_vars_from_partoutput(
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
            ):

    if mode=="backward":
        name_file = partpositfiles[-1]
        fdate=listdates[-1]
    elif mode=="forward":
        name_file = partpositfiles[0]
        fdate=listdates[0]

    if rank == 0:
        if verbose:
            print("\n      --> Getting IDs of parcels over the target region")

    parcel_ids=get_parcels_over_target_region( variables=['lon','lat'], name_file=name_file, file_mask=file_mask, maskname=maskname, maskvar_lon=maskvar_lon, maskvar_lat=maskvar_lat, mask_value=mask_value,fdate=fdate)

    if rank == 0:
        if verbose:
            print("Reading | " + model + " -> ", name_file)

    submatrix = read_flex11_file(name_file,
            parcel_ids,
            fdate,
            lon_left_lower_corner,
            lat_left_lower_corner,
            lon_right_upper_corner,
            lat_right_upper_corner,)


    if mode=="backward":
        name_file2 = partpositfiles[-2]
        fdate=listdates[-2]
    elif mode=="forward":
        name_file2 = partpositfiles[1]
        fdate=listdates[1]


    if rank == 0:
        if verbose:
            print("Reading | " + model + " -> ", name_file2)

    matrix_i = read_flex11_file(name_file2,
            parcel_ids,
            fdate,
            lon_left_lower_corner,
            lat_left_lower_corner,
            lon_right_upper_corner,
            lat_right_upper_corner,)


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
       listdates = listdates[2:]
       listdates=listdates[::-1]
    local_list = partpositfiles[start:stop]
    local_dates=listdates[start:stop]
    # local_list=partpositfiles[:-2][rank::size]
    local_results = np.empty((len(local_list), dimX, dimY))
    local_results = read_by_proccesor(
                    verbose,
                    local_list,
                    local_dates,
                    parcel_ids,
                    submatrix,
                    rank,
                    lon_left_lower_corner,
                    lat_left_lower_corner,
                    lon_right_upper_corner,
                    lat_right_upper_corner,
                    model,
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
