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
import datetime
import os
import functools

print = functools.partial(print, flush=True)

def add_specific_humidity(df):
    """
    Calculates specific humidity (g/kg) from HYSPLIT data columns.
    Assumes AIR_TEMP is in Kelvin, PRESSURE is in hPa, and RELHUMID is a %.
    """
    # 1. Convert Temp to Celsius
    temp_c = df['AIR_TEMP'] - 273.15

    # 2. Saturation Vapor Pressure (Tetens formula)
    e_s = 6.112 * np.exp((17.67 * temp_c) / (temp_c + 243.5))

    # 3. Actual Vapor Pressure
    e = e_s * (df['RELHUMID'] / 100.0)

    # 4. Specific Humidity in kg/kg, then multiplied by 1000 for g/kg
    q_kg_kg = (0.622 * e) / (df['PRESSURE'] - (0.378 * e))
    df['SPCHUMID'] = q_kg_kg

    return df

def is_number(s):
    """Helper function to properly identify negative numbers and decimals."""
    try:
        float(s)
        return True
    except ValueError:
        return False

def is_new_row(parts):
    """
    HYSPLIT trajectory rows always start with exactly these integers:
    TRAJ_ID, MET_GRID, YEAR, MONTH, DAY.
    Wrapped variable text (like wind floats '1.9 -0.9') will fail this test.
    """
    if len(parts) < 5:
        return False
    # isdigit() rejects decimals, meaning it correctly identifies pure integers
    return all(p.lstrip('-').isdigit() for p in parts[:5])

def hysplit_to_dataframe(input_file):
    base_columns = [
        'TRAJ_ID', 'MET_GRID', 'YEAR', 'MONTH', 'DAY',
        'HOUR', 'MINUTE', 'FORECAST_HR', 'AGE',
        'LAT', 'LON', 'HEIGHT'
    ]

    if not os.path.isfile(input_file):
        return pd.DataFrame(columns=base_columns)

    extra_columns = []
    data_start_line = 0

    # 1. Parse the text file into memory
    with open(input_file, 'r') as file:
        lines = file.readlines()

    # 2. Safely find the diagnostic variables header
    for i, line in enumerate(lines):
        parts = line.split()
        if not parts:
            continue

        if len(parts) > 1 and is_number(parts[0]):
            num_vars = int(float(parts[0]))

            if num_vars < 50 and not is_number(parts[1]):
                extra_columns.extend(parts[1:])

                j = i + 1
                while len(extra_columns) < num_vars and j < len(lines):
                    next_parts = lines[j].split()

                    # Safety Break: If we hit trajectory integers, the header lied. Stop looking.
                    if next_parts and is_new_row(next_parts):
                        break

                    extra_columns.extend(next_parts)
                    j += 1

                data_start_line = j
                break

    if not extra_columns:
        print("Error: Could not find the diagnostic variables header.")
        return pd.DataFrame(columns=base_columns)

    # 3. Read Data with Self-Correcting Buffer
    expected_row_len = len(base_columns) + len(extra_columns)
    data_rows = []
    current_row = []

    for line in lines[data_start_line:]:
        parts = line.split()
        if not parts:
            continue

        # If we detect the start of a brand new trajectory row...
        if is_new_row(parts):
            # ...we discard any leftover mangled data to prevent cascading errors
            current_row = parts.copy()
        else:
            # ...otherwise, it is wrapped data (like wind floats) belonging to the current row
            current_row.extend(parts)

        # Check if our row buffer is now perfectly full
        if len(current_row) == expected_row_len:
            data_rows.append(current_row)
            current_row = [] # Reset the buffer for the next row

    if not data_rows:
        print("Error: Header found, but no valid formatted trajectory data was extracted.")
        return pd.DataFrame(columns=base_columns)

    # 4. Build and return the clean DataFrame
    final_columns = base_columns + extra_columns



    df = pd.DataFrame(data_rows, columns=final_columns)

    # 4. Convert the string data into numeric values (floats/integers)
    df = df.apply(pd.to_numeric, errors='ignore')


    df['YEAR'] = df['YEAR'].apply(lambda y: y + 2000 if y < 50 else (y + 1900 if y < 100 else y))

    # 2. Add the temporary SECOND column filled with zeros
    df['SECOND'] = 0

    # 3. Combine the columns into a true Pandas datetime object
    datetime_col = pd.to_datetime(df[['YEAR', 'MONTH', 'DAY', 'HOUR', 'MINUTE', 'SECOND']])

    # 4. Format the datetime objects into the specific string layout
    df['DATES'] = datetime_col.dt.strftime('%Y%m%d %H%M%S')

    # Optional: Drop the temporary SECOND column
    df = df.drop(columns=['SECOND'])
    if not "SPCHUMID" in df.columns:
        if "RELHUMID" in df.columns:
            df = add_specific_humidity(df)

    return df


def generate_hysplit_filename(mode, dtime, totaltime, fecha, path, calendar, file_extension):

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
        dt_obj = datetime.datetime.strptime(fecha, "%Y-%m-%d %H:%M:%S")
        formatted_string = dt_obj.strftime("%Y%m%d_%H%M%S")


        name = path + "hysplit_particles_" + formatted_string.split('_')[0]+ "_"+ formatted_string.split('_')[1]
        array = np.arange(nhour, 0, -dtime)
        for i in array:
            a = str(time_calcminutes(fecha, float(i) * (-1)))
            var1 = a.split(" ")
            var11 = var1[0].split("-")
            var12 = var1[1].split(":")
            fecha_dia = str(
                var11[0] + var11[1] + var11[2] )

            fecha_hour = str(var12[0] + var12[1] + var12[2])




            if (
                calendar == "365d"
                and int(var11[0]) % 4 == 0
                and int(var11[1]) == 2
                and int(var11[2]) == 29
            ):
                msg = "nothing to do"
            else:

                #list_fecha = np.append(list_fecha, name)

                listdates = np.append(listdates, int(fecha_dia+fecha_hour))
        fecha_ = fecha.split(" ")
        var11 = fecha_[0].split("-")
        var12 = fecha_[1].split(":")
        fecha_dia = str(
                var11[0] + var11[1] + var11[2]   )

        fecha_hour = str(var12[0] + var12[1] + var12[2])


        #list_fecha = np.append(list_fecha, name)
        listdates = np.append(listdates, int(fecha_dia+fecha_hour))

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
                var11[0] + var11[1] + var11[2])

            fecha_hour = str(var12[0] + var12[1] + var12[2])




            if (
                calendar == "365d"
                and int(var11[0]) % 4 == 0
                and int(var11[1]) == 2
                and int(var11[2]) == 29
            ):
                msg = "nothong to do"
            else:


                listdates = np.append(listdates, int(fecha_dia+fecha_hour))

    return [name], listdates



def checking_raw_hysplit_partposti_files(partpositfiles, listdates, file_extension):
    errors = ""
    find_error = False


    df = hysplit_to_dataframe(partpositfiles[0])

    if df.empty:
        errors = errors + "     (**) ERROR: No such file or directory or it is corrupt: " + partpositfiles[0] + "\n"
        find_error = True
        return errors, find_error,None


    df['DATES']=pd.to_datetime(df['DATES'], format='%Y%m%d %H%M%S')
    string_array = listdates.astype(np.int64).astype(str)
    # 2. Parse the strings as datetimes using the 14-digit format
    dates_dt = pd.to_datetime(string_array, format='%Y%m%d%H%M%S')



    subset_df = df[df['DATES'].isin(dates_dt)]

    # If you want to reset the index of the new DataFrame so it starts from 0:
    subset_df = subset_df.reset_index(drop=True)


    missing_dates = np.setdiff1d(dates_dt, subset_df['DATES'].unique())

    for i, date_check in enumerate(missing_dates):
        errors = errors + f"     (**) ERROR: This date {date_check} is missing in the HYSPLIT : " + partpositfiles[0] + " file \n"
        find_error = True



    return errors, find_error, subset_df


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
        "LON": 1,
        "LAT": 2,
        "SPCHUMID": 3,
        "HEIGHT": 4,
        "to": 5,
        "rho": 6,
        "hmix": 7,
        "tro": 8,
        "AIR_TEMP": 9,
        "m": 10,
        "THETA": 11,
        "PRESSURE": 13,
    }

    return array_indices[variable]




def get_parcels_over_target_region(variables,
    df, file_mask, maskname, maskvar_lon, maskvar_lat, mask_value,fdate
):
    """
    Optimized function to get parcels over target region.
    Removed 'variables' list because only lon/lat/id are needed for filtering.
    """

    # --- 1. Load Particle Data (Only what is needed) ---
    # Using 'with' ensures the file is closed automatically

    part_id = df["TRAJ_ID"].values.flatten()
    part_lon = df["LON"].values.flatten()
    part_lat = df["LAT"].values.flatten()
    part_pres = df["PRESSURE"].values.flatten()

    # --- 2. Normalize Longitude ---
    # In-place operation is slightly faster and saves memory
    part_lon = (part_lon + 180) % 360 - 180

    # --- 3. Filter NaNs (Vectorized) ---
    # Create a boolean mask for valid coordinates
    valid_mask = np.isfinite(part_lon) & np.isfinite(part_lat)& np.isfinite(part_pres)

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
    df_subset,
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

    a1 = np.arange(len(local_dates))
    dx, dy = submatrix.shape
    tensor_local = np.ones((len(local_dates), dx, dy)) * (-999.9)
    for i in a1:
        df = df_subset[df_subset['DATES']==local_dates[i]]
        if verbose:
            print("Reading | " + model + " -> ", local_dates[i])
        matrix_i = read_file(df,
                        parcel_ids,
                        local_dates[i],
                        lon_left_lower_corner,
                        lat_left_lower_corner,
                        lon_right_upper_corner,
                        lat_right_upper_corner,
       )

        tensor_local[i, :, :] = matrix_i
    return tensor_local



def read_file(target_ds,
            parcel_ids,
            fdate,
            lon_left_lower_corner,
            lat_left_lower_corner,
            lon_right_upper_corner,
            lat_right_upper_corner,
       ):


    variables = ['LAT', 'LON', 'HEIGHT', 'PRESSURE', 'THETA', 'AIR_TEMP']
    sort_indices = np.argsort(parcel_ids)







    tmp_tensor = np.full((len(parcel_ids), 14), -999.0)
    tmp_tensor[:, 0]=parcel_ids[sort_indices]


    for var_name in variables:

        indice = get_array_indices(var_name)


        # Load raw data

        if var_name in target_ds.columns:
            raw_data = target_ds[var_name].values

            # Apply longitude correction if needed
            if var_name == "LON":
                raw_data = (raw_data + 180) % 360 - 180


            # Apply the pre-calculated sort order (Fast indexing)
            sorted_data = raw_data[sort_indices]

            # Assign to tensor
            tmp_tensor[:, indice] = sorted_data




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

def get_vars_from_hysplit(
                verbose,
                df_subset,
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


    string_array = listdates.astype(np.int64).astype(str)
    # 2. Parse the strings as datetimes using the 14-digit format
    dates_dt = pd.to_datetime(string_array, format='%Y%m%d%H%M%S')
    if mode=="backward":
        df = df_subset[df_subset['DATES']==dates_dt[-1]]
        fdate=dates_dt[-1]
    elif mode=="forward":
        df = df_subset[df_subset['DATES']==dates_dt[0]]
        fdate=dates_dt[0]

    if rank == 0:
        if verbose:
            print("\n      --> Getting IDs of parcels over the target region")




    parcel_ids=get_parcels_over_target_region( variables=['LON','LAT'], df=df, file_mask=file_mask, maskname=maskname, maskvar_lon=maskvar_lon, maskvar_lat=maskvar_lat, mask_value=mask_value,fdate=fdate)


    if rank == 0:
        if verbose:
            print("Reading | " + model + " -> ", fdate)

    submatrix = read_file(df,
            parcel_ids,
            fdate,
            lon_left_lower_corner,
            lat_left_lower_corner,
            lon_right_upper_corner,
            lat_right_upper_corner,
         )


    if mode=="backward":
        fdate=dates_dt[-2]
    elif mode=="forward":
        fdate=dates_dt[1]


    if rank == 0:
        if verbose:
            print("Reading | " + model + " -> ", fdate)

    matrix_i = read_file(df_subset,
            parcel_ids,
            fdate,
            lon_left_lower_corner,
            lat_left_lower_corner,
            lon_right_upper_corner,
            lat_right_upper_corner,
           )


    dimX, dimY = submatrix.shape

    tensor_org = np.ones((len(dates_dt), dimX, dimY)) * (-999.9)
    tensor_org[-1, :, :] = submatrix
    tensor_org[-2, :, :] = matrix_i

    n = len(dates_dt) - 2
    count = n // size
    remainder = n % size

    if rank < remainder:
        start = rank * (count + 1)
        stop = start + count + 1
    else:
        start = rank * count + remainder
        stop = start + count

    if mode=="forward":
       dates_dt = dates_dt[2:]
       dates_dt=dates_dt[::-1]

    local_dates=dates_dt[start:stop]
    # local_list=partpositfiles[:-2][rank::size]
    local_results = np.empty((len(local_dates), dimX, dimY))
    local_results = read_by_proccesor(
                    verbose,
                    df_subset,
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

    #tensor_final = np.empty(
    #    (tensor_org.shape[0], tensor_org.shape[1], tensor_org.shape[2])
    #)
    #tensor_final[:] = -999.9
    tensor_final = tensor_org.copy()
    tensor_final[tensor_final==-999]=-999.9
    tensor_final[:, :, -1]=tensor_final[:, :, -1]*100

    if not "PRESSURE" in df_subset.columns:
        tensor_final[:, :, -1] = calc_pres(
            tensor_org[:, :, 6], tensor_org[:, :, 3], tensor_org[:, :, 9])



    if var_heat_track == "dse":

        tensor_final[:, :, -3] = compute_dry_static_energy(
            tensor_org[:, :, 2], tensor_org[:, :, 9], tensor_org[:, :, 4]
        )
        tensor_final[:, :, -3] = tensor_final[:, :, -3] / 1000

    elif var_heat_track == "potTemp":
        if not "THETA" in df_subset.columns:
            tensor_final[:, :, -3] = compute_theta(
            tensor_org[:, :, 6], tensor_org[:, :, 3], tensor_org[:, :, 9], press_data=True, parpres=tensor_final[:, :, -1])

    else:
        tensor_final[:, :, -3] = compute_theta(
            tensor_org[:, :, 6], tensor_org[:, :, 3], tensor_org[:, :, 9], press_data=True, parpres=tensor_final[:, :, -1])


    tensor_final[:, :, -2] = compute_rh(
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

