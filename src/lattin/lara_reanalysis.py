import xarray as xr
import numpy as np
import pandas as pd
import dask
import time
import os
from pathlib import Path
import logging
from datetime import datetime, timedelta
from tenacity import (
    retry,
    stop_after_attempt,
    wait_exponential,
    retry_if_exception_type,
)
import aiohttp
import scipy.interpolate as interp

from lattin.lattin_functions import (
    compute_dry_static_energy,
    compute_theta,
    compute_rh,
    calc_pres,
    load_mask_grid_NR,
    generate_file,
    funtion_interpol_mascara,
    determine_id_target_region,
    get_lara_config,
)
import gc
import functools
import warnings

warnings.filterwarnings("ignore")
print = functools.partial(print, flush=True)

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def configure_dask_for_stability():
    """Configure Dask for more stable network operations."""
    dask.config.set(
        {
            "distributed.comm.timeouts.connect": "120s",
            "distributed.comm.timeouts.tcp": "120s",
            "distributed.comm.retry.count": 5,
            "array.chunk-size": "128MB",
            "array.slicing.split_large_chunks": True,
        }
    )


@retry(
    stop=stop_after_attempt(5),
    wait=wait_exponential(multiplier=2, min=4, max=120),
    retry=retry_if_exception_type(
        (
            aiohttp.client_exceptions.ClientPayloadError,
            aiohttp.client_exceptions.ClientConnectorError,
            ConnectionError,
            TimeoutError,
            OSError,
        )
    ),
)
def download_single_day(
    zarr_url,
    target_date,
    output_dir,
    var_name="T",
    gridded_variables=["hmix"],
    particle_chunk_size=500000,
):
    """
    Download data for a single day with retry logic.

    Parameters:
    -----------
    zarr_url : str
        URL to the Zarr dataset
    target_date : datetime or str
        Date to download (YYYY-MM-DD)
    output_dir : Path
        Directory to save the daily file
    var_name : str
        Variable name to download
    particle_chunk_size : int
        Number of particles per chunk
    """

    if isinstance(target_date, str):
        target_date = pd.to_datetime(target_date)

    date_str = target_date.strftime("%Y%m%d")
    output_file = output_dir / f"LARA_{date_str}_{var_name}.nc"

    # Skip if file already exists and is complete

    expected_timestamps_pd = pd.date_range(start=date_str, periods=24, freq="H")

    expected_timestamps_np = expected_timestamps_pd.to_numpy()

    if output_file.exists():
        try:
            # Quick check if file is complete (expecting 24 hours for hourly data)
            with xr.open_dataset(output_file) as test_ds:
                expected_timesteps = 24  # For hourly data
                if "time" in test_ds.dims and len(test_ds.time) >= expected_timesteps:
                    actual_times_np = test_ds.time.values

                    present_mask = np.isin(expected_timestamps_np, actual_times_np)
                    if present_mask.all():
                        # print(f"     File {output_file.name} already exists and appears complete ({len(test_ds.time)} timesteps), skipping")
                        return output_file

                    else:
                        print(
                            f"         File {output_file.name} exists but appears incomplete ({len(test_ds.time)} timesteps, expected {expected_timesteps}), re-downloading"
                        )
                    output_file.unlink()

                else:
                    print(
                        f"         File {output_file.name} exists but appears incomplete ({len(test_ds.time)} timesteps, expected {expected_timesteps}), re-downloading"
                    )
                    output_file.unlink()
        except Exception as e:
           
            output_file.unlink()

    print(f"         -> Downloading {var_name} data for {date_str}")

    # Configure storage options for stability - SIMPLIFIED VERSION
    storage_options = {
        "timeout": 600,  # 10 minute timeout in seconds
    }

    try:
        # Open dataset
        ds = xr.open_zarr(zarr_url, storage_options=storage_options)

        # Get time range for the target date (precise to avoid next day)
        start_time = target_date.strftime("%Y-%m-%d 00:00:00")
        end_time = target_date.strftime("%Y-%m-%d 23:59:59")

        # Select time range for this day only
        day_data = ds.sel(time=slice(start_time, end_time))

        # Debug: Check what we actually got
        # logger.info(f"Selected time range: {start_time} to {end_time}")
        if len(day_data.time) > 0:

            # Check if it's really hourly data
            if len(day_data.time) > 1:
                time_diff = pd.to_datetime(day_data.time.values[1]) - pd.to_datetime(
                    day_data.time.values[0]
                )
                

            if len(day_data.time) > 24:
                print(
                    f"           Got {len(day_data.time)} timesteps - expected 24 for hourly data!"
                )
                
                if len(day_data.time) > 5:
                    # logger.info("  ...")
                    for i in range(max(5, len(day_data.time) - 3), len(day_data.time)):
                        logger.info(f"  {i}: {day_data.time.values[i]}")
            elif len(day_data.time) < 24:
                print(
                    f"           Got only {len(day_data.time)} timesteps - expected 24 for hourly data!"
                )

        # Check if we got any data
        if len(day_data.time) == 0:
            print(f"           No data found for date {date_str}")

        
        # Rechunk for efficient processing

        if not var_name in gridded_variables:

            day_data_chunked = day_data.chunk(
                {
                    "particle": particle_chunk_size,
                    "time": min(12, len(day_data.time)),  # 12-hour chunks max
                }
            )
        else:
            day_data_chunked = day_data.chunk(
                {"time": min(12, len(day_data.time))}  # 12-hour chunks max
            )

        # Load data progressively
        
        day_data_loaded = day_data_chunked.load()

        # Save with compression

        encoding = {}
        for var in day_data_loaded.data_vars:
            encoding[var] = {
                "zlib": True,
                "complevel": 4,
                "shuffle": True,
                "fletcher32": True,
            }

        # Add attributes for tracking
        day_data_loaded.attrs["download_date"] = datetime.now().isoformat()
        day_data_loaded.attrs["source_url"] = zarr_url
        day_data_loaded.attrs["date_downloaded"] = date_str

        day_data_loaded.to_netcdf(output_file, engine="netcdf4", encoding=encoding)

        return output_file

    except Exception as e:
        print(f"           Failed to download {date_str}: {e}")
        # Clean up partial file if it exists
        if output_file.exists():
            output_file.unlink()
        raise


def download_monthly_data_by_days_simple(
    zarr_url,
    year,
    month,
    output_dir,
    var_name="T",
    target_days=[1],
    gridded_variables=["hmix"],
    max_workers=2,
):
    """
    Simplified version using the simple download function.
    Use this if the main function has issues with storage options.
    """
    configure_dask_for_stability()

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate list of dates for the month
   

    # Find last day of month
    if month == 12:
        end_date = datetime(year + 1, 1, 1) - timedelta(days=1)
    else:
        end_date = datetime(year, month + 1, 1) - timedelta(days=1)

    dates = []
    for day in target_days:

        current_date = datetime(year, month, day)
        if current_date <= end_date:

            dates = np.append(dates, current_date)

    
    successful_downloads = []
    failed_downloads = []

    # Download each day sequentially using simple function
    for i, date in enumerate(dates):
        try:
            # logger.info(f"Processing day {i+1}/{len(dates)}: {date.strftime('%Y-%m-%d')}")

            output_file = download_single_day_simple(
                zarr_url,
                date,
                output_dir,
                var_name,
                gridded_variables,
            )

            successful_downloads.append(output_file)
            

            # Small delay between downloads to be nice to the server
            if i < len(dates) - 1:  # Don't sleep after last download
                time.sleep(2)

        except Exception as e:
            print(f"           Failed to download {date.strftime('%Y-%m-%d')}: {e}")
            failed_downloads.append((date, str(e)))

            # Continue with next day rather than stopping
            continue

    return successful_downloads, failed_downloads


def download_monthly_data_by_days(
    zarr_url,
    year,
    month,
    output_dir,
    var_name="T",
    target_days=[1],
    gridded_variables=["hmix"],
    max_workers=2,
):
    """
    Download all days in a month as separate files.

    Parameters:
    -----------
    zarr_url : str
        URL to the Zarr dataset
    year : int
        Year to download
    month : int
        Month to download (1-12)
    output_dir : str or Path
        Directory to save daily files
    var_name : str
        Variable name to download
    max_workers : int
        Maximum number of concurrent downloads (keep low to avoid overwhelming server)
    """

    configure_dask_for_stability()

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    
    # Find last day of month
    if month == 12:
        end_date = datetime(year + 1, 1, 1) - timedelta(days=1)
    else:
        end_date = datetime(year, month + 1, 1) - timedelta(days=1)
    dates = []
    for day in target_days:

        current_date = datetime(year, month, day)
        if current_date <= end_date:

            dates = np.append(dates, current_date)

    
    successful_downloads = []
    failed_downloads = []

    # Download each day sequentially (safer than parallel for large files)
    for i, date in enumerate(dates):
        try:
            # logger.info(f"Processing day {i+1}/{len(dates)}: {date.strftime('%Y-%m-%d')}")

            output_file = download_single_day(
                zarr_url,
                date,
                output_dir,
                var_name,
                gridded_variables,
            )

            successful_downloads.append(output_file)
            
            # Small delay between downloads to be nice to the server
            if i < len(dates) - 1:  # Don't sleep after last download
                time.sleep(2)

        except Exception as e:
            
            failed_downloads.append((date, str(e)))

            # Continue with next day rather than stopping
            continue

    return successful_downloads, failed_downloads


def create_monthly_index(daily_files, output_dir):
    """
    Create an index file listing all daily files for easy reference.
    """

    index_file = Path(output_dir) / "daily_files_index.txt"

    with open(index_file, "w") as f:
        f.write("Daily LARA files index\n")
        f.write("=" * 50 + "\n")
        f.write(f"Generated: {datetime.now().isoformat()}\n")
        f.write(f"Total files: {len(daily_files)}\n\n")

        for daily_file in sorted(daily_files):
            size_mb = daily_file.stat().st_size / (1024 * 1024)
            f.write(f"{daily_file.name} ({size_mb:.1f} MB)\n")

    print(f"Created index file: {index_file}")
    return index_file


def retry_failed_downloads(
    zarr_url, failed_downloads, output_dir, var_name="T", gridded_variables=["hmix"]
):
    """
    Retry downloading failed days.
    """

   

    print(f"         -> Retrying {len(failed_downloads)} failed downloads...")

    # Ensure output_dir is a Path object
    output_dir = Path(output_dir)

    successful_retries = []
    still_failed = []

    for date, original_error in failed_downloads:
        try:
            print(f"           Retrying {date.strftime('%Y-%m-%d')}")
            output_file = download_single_day(
                zarr_url, date, output_dir, var_name, gridded_variables
            )
            successful_retries.append(output_file)
            print(f"           Retry successful for {date.strftime('%Y-%m-%d')}")

        except Exception as e:
            print(f"           Retry failed for {date.strftime('%Y-%m-%d')}: {e}")
            still_failed.append((date, str(e)))

    print(
        f"           Retry results: {len(successful_retries)} successful, {len(still_failed)} still failed"
    )

    return successful_retries, still_failed


# Alternative function with simplified storage options if the above still doesn't work
def download_single_day_simple(
    zarr_url,
    target_date,
    output_dir,
    var_name="T",
    gridded_variables=["hmix"],
    particle_chunk_size=500000,
):
    """
    Alternative download function with simplified storage options.
    Use this if the main function still has issues.
    """

    if isinstance(target_date, str):
        target_date = pd.to_datetime(target_date)

    date_str = target_date.strftime("%Y%m%d")
    output_file = output_dir / f"LARA_{date_str}_{var_name}.nc"

    expected_timestamps_pd = pd.date_range(start=date_str, periods=24, freq="H")

    
    expected_timestamps_np = expected_timestamps_pd.to_numpy()

    if output_file.exists():
        try:
            # Quick check if file is complete (expecting 24 hours for hourly data)
            with xr.open_dataset(output_file) as test_ds:
                expected_timesteps = 24  # For hourly data
                if "time" in test_ds.dims and len(test_ds.time) >= expected_timesteps:
                    actual_times_np = test_ds.time.values

                    present_mask = np.isin(expected_timestamps_np, actual_times_np)
                    if present_mask.all():
                        
                        return output_file

                    else:
                        print(
                            f"         File {output_file.name} exists but appears incomplete ({len(test_ds.time)} timesteps, expected {expected_timesteps}), re-downloading"
                        )
                    output_file.unlink()

                else:
                    print(
                        f"         File {output_file.name} exists but appears incomplete ({len(test_ds.time)} timesteps, expected {expected_timesteps}), re-downloading"
                    )
                    output_file.unlink()
        except Exception as e:
            
            output_file.unlink()

    print(f"         -> Downloading {var_name} data for {date_str}")

    try:
        # Open dataset with minimal storage options
        ds = xr.open_zarr(zarr_url)

        # Get time range for the target date (more precise)
        start_time = target_date.strftime("%Y-%m-%d 00:00:00")
        end_time = target_date.strftime("%Y-%m-%d 23:59:59")

        # Select time range for this day
        day_data = ds.sel(time=slice(start_time, end_time))

       
        if len(day_data.time) > 0:
            
            # Check if it's really hourly data
            if len(day_data.time) > 1:
                time_diff = pd.to_datetime(day_data.time.values[1]) - pd.to_datetime(
                    day_data.time.values[0]
                )
                

            if len(day_data.time) > 24:
                print(
                    f"       Got {len(day_data.time)} timesteps - expected 24 for hourly data!"
                )
                
                for i in range(min(5, len(day_data.time))):
                    logger.info(f"  {i}: {day_data.time.values[i]}")
                if len(day_data.time) > 5:
                    logger.info("  ...")
                    for i in range(max(5, len(day_data.time) - 3), len(day_data.time)):
                        logger.info(f"  {i}: {day_data.time.values[i]}")
            elif len(day_data.time) < 24:
                print(
                    f"         Got only {len(day_data.time)} timesteps - expected 24 for hourly data!"
                )

        # Check if we got any data
        if len(day_data.time) == 0:
            print(f"           No data found for date {date_str}")

        

        # Rechunk for efficient processing
        if not var_name in gridded_variables:
            day_data_chunked = day_data.chunk(
                {
                    "particle": particle_chunk_size,
                    "time": min(12, len(day_data.time)),  # 12-hour chunks max
                }
            )
        else:
            day_data_chunked = day_data.chunk(
                {"time": min(12, len(day_data.time))}  # 12-hour chunks max
            )

        # Load data progressively
        
        day_data_loaded = day_data_chunked.load()

        # Save with compression
        
        encoding = {}
        for var in day_data_loaded.data_vars:
            encoding[var] = {
                "zlib": True,
                "complevel": 4,
                "shuffle": True,
                "fletcher32": True,
            }

        # Add attributes for tracking
        day_data_loaded.attrs["download_date"] = datetime.now().isoformat()
        day_data_loaded.attrs["source_url"] = zarr_url
        day_data_loaded.attrs["date_downloaded"] = date_str

        day_data_loaded.to_netcdf(output_file, engine="netcdf4", encoding=encoding)

        

        return output_file

    except Exception as e:
        print(f"           Failed to download {date_str}: {e}")
        # Clean up partial file if it exists
        if output_file.exists():
            output_file.unlink()
        raise


def get_lara_url(year_to_check, month):
    """
    LARA VARIABLEs:
    - lon: Longitude of Particle
    - lat: Latitude of Particle
    - z: Height of Particle
    - pv: Potential Vorticity at Particle location
    - sh: Specific Humidity at Particle location
    - rho: Density at Particle location
    - T: Temperature a Particle location
    - prs: Pressure at Particle location

    Averaged VARIABLEs:
    - lon_av: Longitude
    - lat_av: Latitude
    - z_av: Height
    - pv_av: Potential Vorticity at Particle location
    - sh_av: Specific Humidity at Particle location
    - rho_av: Density at Particle location
    - T_av: Temperature a Particle location
    - prs_av: Pressure at Particle location

    gridded VARIABLEs:
    - hmix: Atmospheric boundary layer height
    - tro: Tropopause height
    - to: Topography (fixed in time)

    """

    # 1. Parse the periods string into a list of tuples (start_year, end_year, original_string)
    # We assume a period S-E includes years from S up to E, including, E.
    # So, 1940-1947 covers 1940, 1941, ..., 1946, 1947. The year  1947 falls also into 1947-1954.

    lara_base_url, periods_str, _, _, _, _ = get_lara_config()

    parsed_periods = []
    for line in periods_str.strip().splitlines():
        if not line.strip():  # Skip any empty lines that might result from strip()
            continue
        parts = line.strip().split("-")
        start_year = int(parts[0])
        end_year = int(parts[1])
        parsed_periods.append((start_year, end_year, line.strip()))

    # Sort periods by start year, just in case they weren't already
    parsed_periods.sort(key=lambda x: x[0])

    find_period = None
    for start_year, end_year, period_str in parsed_periods:

        if start_year == 1940:
            start_year = 1939
        # A year belongs to the period [start_year, end_year)
        # This means start_year <= year_to_check <=end_year
        if start_year <= int(year_to_check) <= end_year:
            find_period = period_str

    if not find_period == None:
        return f"{lara_base_url}/{find_period}/{int(year_to_check)}/{str(int(month)).zfill(2)}"

    else:
        return None


def get_date_list_lara(
    start_year,
    start_month,
    start_day,
    start_hour,
    start_min,
    end_year,
    end_month,
    end_day,
    end_hour,
    end_min,
    total_tracking_time,
    time_step,
    mode,
    calendar,
):

    """
    Generate a list of unique dates from start and end date.

    Parameters
    ----------
    start_year : int
        start year
    start_month : int
        start month
    start_day : int
        start day
    start_hour : int
        start hour
    start_min : int
        start minute
    end_year : int
        end year
    end_month : int
        end month
    end_day : int
        end day
    end_hour : int
        end hour
    end_min : int
        end minute
    total_tracking_time : int
        total tracking time [minutes]
    time_step : int
        time step [minutes]
    mode : str
        "backward" or "forward"
    calendar : str
        calendar type "365d" or "366d"

    Returns
    -------
    datetime_series : pandas.Series
        A pandas series of unique datetime objects
    """
    start_date = datetime(
        int(start_year),
        int(start_month),
        int(start_day),
        int(start_hour),
        int(start_min),
        0,
    )
    end_date = datetime(
        int(end_year), int(end_month), int(end_day), int(end_hour), int(end_min), 0
    )

    _, list_dates_start = generate_file(
        mode,
        time_step,
        total_tracking_time,
        start_date.strftime("%Y-%m-%d %H:%M:%S"),
        "./",
        False,
        calendar,
    )
    _, list_dates_end = generate_file(
        mode,
        time_step,
        total_tracking_time,
        end_date.strftime("%Y-%m-%d %H:%M:%S"),
        "./",
        False,
        calendar,
    )

    combined_list_dates = np.concatenate((list_dates_start, list_dates_end))

    unique_dates = np.unique(combined_list_dates).astype(np.int64)

    datetime_series = pd.to_datetime(unique_dates, format="%Y%m%d%H%M%S")

    return datetime_series


# Main execution function
def get_lara_variables_from_repo(output_basedir, lara_list_of_dates):
    """
    Main function to download monthly data as daily files.
    Modify the parameters below for your specific needs.
    """

    years = np.unique(lara_list_of_dates.year)
    _, _, _, _, variables, gridded_variables = get_lara_config()

    
    for year in years:
        year_dates = lara_list_of_dates[lara_list_of_dates.year == year]

        months = np.unique(year_dates.month)

        for month in months:
            month_dates = year_dates[year_dates.month == month]

            target_days = np.unique(month_dates.day)

            temp_url = get_lara_url(year, month)

            for var_name in variables:
                #
                zarr_url = f"{temp_url}/{var_name}"

                output_dir = f"{output_basedir}/{var_name}/{int(year)}/{str(int(month)).zfill(2)}"
                try:
                    
                    print(f"     --- Checking files for {var_name}")
                    try:

                        successful_files, failed_downloads = (
                            download_monthly_data_by_days_simple(
                                zarr_url=zarr_url,
                                year=year,
                                month=month,
                                output_dir=output_dir,
                                var_name=var_name,
                                target_days=target_days,
                                gridded_variables=gridded_variables,
                            )
                        )

                    except Exception as e:
                        
                        print(f"     --- Checking files for {var_name}")
                        successful_files, failed_downloads = (
                            download_monthly_data_by_days(
                                zarr_url=zarr_url,
                                year=year,
                                month=month,
                                output_dir=output_dir,
                                var_name=var_name,
                                target_days=target_days,
                                gridded_variables=gridded_variables,
                            )
                        )

                    # Retry failed downloads
                    if failed_downloads:
                        retry_successful, still_failed = retry_failed_downloads(
                            zarr_url, failed_downloads, Path(output_dir), var_name
                        )
                        successful_files.extend(retry_successful)

                    if len(successful_files) == len(target_days):
                        print(f"         PASSED")

                except Exception as e:
                    print(f"     -- Download failed: {e} for {year} {month} {var_name}")
                    raise


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
        "xmass": 10,
        "pottemp": 11,
        "prs": 13,
    }

    return array_indices[variable]


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
    lat_lon = np.column_stack((lon_.flatten(), lat_.flatten()))

    # Create interpolator
    prsInterpu = interp.NearestNDInterpolator(lat_lon, grid_variable.flatten())

    return prsInterpu(parcels)


def D_id_python_exact(value_mascara, value_mask, len_value_mascara):
    
    """
    Parameters
    ----------
    value_mascara : 1D array
        1D array of mask values
    value_mask : int
        Value to be searched in value_mascara
    len_value_mascara : int
        Length of value_mascara

    Returns
    -------
    1D array
        Array with indices of value_mascara where value_mask is found, -999.9 elsewhere
    """
    return np.where(value_mascara == value_mask, np.arange(len_value_mascara), -999.9)


def extract_valid_ids(data, idparts):
    # Sort data by first column

    # Remove duplicates - keep first occurrence only
    """
    Filters and processes data based on given IDs, ensuring all specified IDs are present in the output.

    Parameters
    ----------
    data : np.ndarray
        2D array with data, where the first column contains IDs.
    idparts : array_like
        Array of IDs that need to be present in the returned data.

    Returns
    -------
    np.ndarray
        A 2D array with rows corresponding to the IDs in `idparts`. If an ID in `idparts` is not found
        in `data`, a row with that ID and other values set to -999 is included. The result is sorted
        by ID to match the order of `idparts`.

    Notes
    -----
    - The function removes duplicate IDs, keeping the first occurrence.
    - Rows with IDs not present in `idparts` are filtered out.
    - Missing IDs in the result are filled with rows of -999 values, except for the ID column.
    """

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
        missing_rows = np.full((len(missing_ids), data.shape[1]), -999)
        missing_rows[:, 0] = missing_ids

        # Combine existing and missing data
        data = np.vstack([filtered_data, missing_rows])

        # Sort by ID to maintain order matching idparts
        data = data[np.argsort(data[:, 0])]
    else:
        data = filtered_data

    return data


#####################################################################################################################################################
def extract_valid_date_indices(
    zarr_url, month_dates, parcel_ids, var_name, date_str, from_https=False
):

    """
    Extract valid date indices from a given Zarr URL.

    Parameters
    ----------
    zarr_url : str
        URL of the Zarr file.
    month_dates : array_like
        Array of datetime64 dates.
    parcel_ids : array_like
        Array of parcel IDs.
    var_name : str
        Name of the variable to access.
    date_str : str
        String representation of the date.
    from_https : bool, optional
        Whether to use HTTPS to access the Zarr file. Default is False.

    Returns
    -------
    target_ds : xr.Dataset
        Extracted dataset with the specified parcel IDs and time range.

    Raises
    ------
    ValueError
        If there is an error accessing the Zarr file or if no parcel IDs are found.
    """
    if from_https:
        ds = xr.open_zarr(zarr_url)
        target_dates = month_dates
    else:

        ds = xr.open_dataset(zarr_url)
        target_dates = date_str
    try:
        target_ds_ = ds.sel(time=target_dates.strftime("%Y-%m-%d %H:%M:%S"))
    except Exception as e:
        raise ValueError(f"    Error in file for {var_name} on date")

    ds.close()

    try:
        target_ds = target_ds_.sel(particle=parcel_ids)
    except Exception as e:
        raise ValueError(
            f"    No parcels IDs found. Something is wrong with your maskfile or the accessed input data"
        )

    del target_ds_
    gc.collect()

    return target_ds


def reading_lara_variables_from_zarr(
    verbose,
    listdates,
    lara_basedir,
    variables,
    gridded_variables,
    lon_left_lower_corner,
    lat_left_lower_corner,
    lon_right_upper_corner,
    lat_right_upper_corner,
    start_date,
    parcel_ids,
    allproc=True,
    rank=0,
    from_https=False,
):

    """
    Reads variables from LARA reanalysis data for a given set of dates and parcel IDs.

    Parameters
    ----------
    verbose : bool
        Whether to print progress information.
    listdates : array_like
        Array of dates to read.
    lara_basedir : str
        Base directory for LARA reanalysis data.
    variables : list
        List of variables to read.
    gridded_variables : list
        List of gridded variables to read.
    lon_left_lower_corner : float
        Lower left longitude of the target region.
    lat_left_lower_corner : float
        Lower left latitude of the target region.
    lon_right_upper_corner : float
        Upper right longitude of the target region.
    lat_right_upper_corner : float
        Upper right latitude of the target region.
    start_date : datetime
        Starting date of the data.
    parcel_ids : array_like
        Array of parcel IDs to read.
    allproc : bool
        Whether all processes should print progress information.
    rank : int
        Rank of the current process.
    from_https : bool, optional
        Whether to use HTTPS to access the Zarr file. Default is False.

    Returns
    -------
    tmp_tensor : array_like
        Array of shape (len(listdates), len(parcel_ids), 14) containing the data for the specified variables.

    Raises
    ------
    ValueError
        If there is an error accessing the Zarr file or if no parcel IDs are found.
    """
    
    tmp_gridded_vars = set(gridded_variables)
    tmp_variables = [item for item in variables if item not in tmp_gridded_vars]

    tmp_tensor = np.full((len(listdates), len(parcel_ids), 14), -999.0)

    sp = "        "

    

    # index_target_date = listdates.tolist().index(start_date)
    listdates_list = listdates.tolist()
    years = np.unique(listdates.year)
    for year in years:
        year_dates = listdates[listdates.year == year]

        months = np.unique(year_dates.month)

        for month in months:
            month_dates = listdates[
                (listdates.year == year) & (listdates.month == month)
            ]

            if from_https:
                date_ranges = [month_dates[0]]
            else:
                date_ranges = month_dates

            for date_str_ in date_ranges:

                date_str = date_str_.strftime("%Y%m%d")

                temp_url = get_lara_url(year, month)
                for var_name in tmp_variables:

                    if from_https:
                        zarr_url = f"{temp_url}/{var_name}"
                        fact = 1
                    else:
                        zarr_url = f"{lara_basedir}/{var_name}/{int(year)}/{str(int(month)).zfill(2)}/LARA_{date_str}_{var_name}.nc"
                        fact = 2

                    target_ds = extract_valid_date_indices(
                        zarr_url,
                        month_dates,
                        parcel_ids,
                        var_name,
                        date_str_,
                        from_https,
                    )

                    indice = get_array_indices(var_name)

                    ds_times = target_ds["time"]

                    if ds_times.size == 1:
                        ds_times = [ds_times]

                    for trk_date in ds_times:

                        if verbose:
                            if allproc:
                                if rank == 0:
                                    print(
                                        f'Reading | LARA reanalysis data -> {zarr_url}{sp[:len(sp)-fact*len(var_name)]} -> {pd.Timestamp(trk_date.values.item()).strftime("%Y-%m-%d %H:%M:%S")}'
                                    )
                            else:
                                print(
                                    f'Reading | LARA reanalysis data -> {zarr_url}{sp[:len(sp)-fact*len(var_name)]} -> {pd.Timestamp(trk_date.values.item()).strftime("%Y-%m-%d %H:%M:%S")}'
                                )
                        index_target_date = listdates_list.index(trk_date)

                        if from_https:
                            aux_data_ds = target_ds.sel(time=trk_date)
                        else:
                            aux_data_ds = target_ds

                        aux_data = aux_data_ds[var_name].values

                        if var_name == "lon":
                            aux_data = (aux_data + 180) % 360 - 180

                        tmp_data = np.column_stack(
                            (aux_data_ds["particle"].values, aux_data)
                        )
                        tmp_data = tmp_data[np.argsort(tmp_data[:, 0])]

                        tmp_data = extract_valid_ids(tmp_data, parcel_ids)

                        tmp_tensor[index_target_date, :, indice] = tmp_data[:, 1]
                        tmp_tensor[index_target_date, :, 0] = parcel_ids

                        del aux_data

                    del target_ds
                    gc.collect()

                for var_name in gridded_variables:

                    if from_https:
                        zarr_url = f"{temp_url}/{var_name}"
                        target_dates = month_dates
                        ds = xr.open_zarr(zarr_url)
                        fact = 1
                    else:
                        zarr_url = f"{lara_basedir}/{var_name}/{int(year)}/{str(int(month)).zfill(2)}/LARA_{date_str}_{var_name}.nc"
                        target_dates = date_str_
                        ds = xr.open_dataset(zarr_url)
                        fact = 2

                    try:
                        target_ds = ds.sel(
                            time=target_dates.strftime("%Y-%m-%d %H:%M:%S")
                        )
                    except Exception as e:
                        raise ValueError(f"    Error in file for {var_name} on date")

                    indice = get_array_indices(var_name)

                    ds_times = target_ds["time"]

                    if ds_times.size == 1:
                        ds_times = [ds_times]

                    for trk_date in ds_times:
                        if verbose:
                            if allproc:
                                if rank == 0:
                                    print(
                                        f'Reading | LARA reanalysis data -> {zarr_url}{sp[:len(sp)-fact*len(var_name)]} -> {pd.Timestamp(trk_date.values.item()).strftime("%Y-%m-%d %H:%M:%S")}'
                                    )
                            else:
                                print(
                                    f'Reading | LARA reanalysis data -> {zarr_url}{sp[:len(sp)-fact*len(var_name)]} -> {pd.Timestamp(trk_date.values.item()).strftime("%Y-%m-%d %H:%M:%S")}'
                                )

                        index_target_date = listdates_list.index(trk_date)

                        if from_https:
                            aux_data = target_ds.sel(time=trk_date)
                        else:
                            aux_data = target_ds

                        values_data = interpol_gridded_to_parcels(
                            aux_data["latitude"].values,
                            aux_data["longitude"].values,
                            aux_data[var_name].values,
                            tmp_tensor[index_target_date, :, 1:3],
                        )

                        tmp_tensor[index_target_date, :, indice] = values_data

                    del target_ds
                    gc.collect()

    if (
        lon_left_lower_corner > -180
        or lon_right_upper_corner < 180
        or lat_left_lower_corner > -90
        or lat_right_upper_corner < 90
    ):

        # Check for out of bounds coordinates
        lon_out_of_bounds = (tmp_tensor[:, :, 1] <= lon_left_lower_corner) | (
            tmp_tensor[:, :, 1] >= lon_right_upper_corner
        )
        lat_out_of_bounds = (tmp_tensor[:, :, 2] <= lat_left_lower_corner) | (
            tmp_tensor[:, :, 2] >= lat_right_upper_corner
        )

        # Combine masks for points that are out of bounds
        bad_points_mask = lon_out_of_bounds | lat_out_of_bounds

        # Set all elements from index 1 onwards to -999 for out-of-bounds points
        tmp_tensor[bad_points_mask, 1:] = -999.0

    rows_to_keep = np.any(tmp_tensor[:, :, 0] != -999.0, axis=1)
    tmp_tensor = tmp_tensor[rows_to_keep]

    return tmp_tensor


def get_lara_variables_from_zarr_backward(
    verbose,
    listdates,
    lara_basedir,
    variables,
    gridded_variables,
    lon_left_lower_corner,
    lat_left_lower_corner,
    lon_right_upper_corner,
    lat_right_upper_corner,
    start_date,
    parcel_ids,
    var_heat_track,
    from_https,
    rank,
    size,
    comm,
):

    """
    Reads LARA variables from Zarr files for a list of dates (from the most recent to the oldest) 
    and for a given set of parcel IDs. The results are stored in a 3D array and returned to the root process.
    
    Parameters
    ----------
    verbose : boolean
        print information about the process
    listdates : list
        list of dates for which the data is read
    lara_basedir : str
        path to the LARA data repository
    variables : list
        list of variables to read
    gridded_variables : list
        list of variables to read from the gridded data
    lon_left_lower_corner : float
        longitude of the lower left corner of the target region
    lat_left_lower_corner : float
        latitude of the lower left corner of the target region
    lon_right_upper_corner : float
        longitude of the upper right corner of the target region
    lat_right_upper_corner : float
        latitude of the upper right corner of the target region
    start_date : datetime
        start date of the simulation
    parcel_ids : array
        array of parcel IDs
    var_heat_track : str
        variable to track heat (e.g. 'theta', 'T')
    from_https : boolean
        read data from https or local directory
    rank : int
        rank of the MPI process
    size : int
        number of MPI processes
    comm : MPI communicator
        MPI communicator
    
    Returns
    -------
    tensor_final : 3D array
        3D array with the data read from the Zarr files
    """
    tmp_tensor = reading_lara_variables_from_zarr(
        verbose,
        listdates[-2:],
        lara_basedir,
        variables,
        gridded_variables,
        lon_left_lower_corner,
        lat_left_lower_corner,
        lon_right_upper_corner,
        lat_right_upper_corner,
        start_date,
        parcel_ids,
        allproc=True,
        rank=rank,
        from_https=from_https,
    )

    _, dimX, dimY = tmp_tensor.shape
    tensor_final = np.ones((len(listdates), dimX, dimY)) * (-999.9)
    tensor_final[-1, :, :] = tmp_tensor[1, :, :]
    tensor_final[-2, :, :] = tmp_tensor[0, :, :]

    tmp_listdates = listdates[:-2]
    n = len(tmp_listdates)
    count = n // size
    remainder = n % size

    if rank < remainder:
        start = rank * (count + 1)
        stop = start + count + 1
    else:
        start = rank * count + remainder
        stop = start + count

    local_list = tmp_listdates[start:stop]

    local_results = np.empty((len(local_list), dimX, dimY))
    local_results = reading_lara_variables_from_zarr(
        verbose,
        local_list,
        lara_basedir,
        variables,
        gridded_variables,
        lon_left_lower_corner,
        lat_left_lower_corner,
        lon_right_upper_corner,
        lat_right_upper_corner,
        start_date,
        parcel_ids,
        allproc=False,
        rank=rank,
        from_https=from_https,
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

        tensor_final[int(i_start[0]) : int(i_stop[0]), :, :] = final_results
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
            tensor_final[int(i_start[i]) : int(i_stop[i]), :, :] = tmp
    comm.Bcast(tensor_final, root=0)

    tensor_final = compute_other_quantities(tensor_final, var_heat_track)

    if rank == 0:
        print(f"\n     => Data reading: Finished")

    gc.collect()
    return tensor_final




def get_lara_variables_from_zarr_forward(
    verbose,
    listdates,
    lara_basedir,
    variables,
    gridded_variables,
    lon_left_lower_corner,
    lat_left_lower_corner,
    lon_right_upper_corner,
    lat_right_upper_corner,
    start_date,
    parcel_ids,
    var_heat_track,
    from_https,
    rank,
    size,
    comm,
):

    """
    Reads LARA variables from Zarr files for a list of dates (from the most recent to the oldest)
    and for a given set of parcel IDs. The results are stored in a 3D array and returned to the root process.

    Parameters
    ----------
    verbose : boolean
        print information about the process
    listdates : list
        list of dates for which the data is read
    lara_basedir : str
        path to the LARA data repository
    variables : list
        list of variables to read
    gridded_variables : list
        list of variables to read from the gridded data
    lon_left_lower_corner : float
        longitude of the lower left corner of the target region
    lat_left_lower_corner : float
        latitude of the lower left corner of the target region
    lon_right_upper_corner : float
        longitude of the upper right corner of the target region
    lat_right_upper_corner : float
        latitude of the upper right corner of the target region
    start_date : datetime
        start date of the simulation
    parcel_ids : array
        array of parcel IDs
    var_heat_track : str
        variable to track heat (e.g. 'theta', 'T')
    from_https : boolean
        read data from https or local directory
    rank : int
        rank of the MPI process
    size : int
        number of MPI processes
    comm : MPI communicator
        MPI communicator

    Returns
    -------
    tensor_final : 3D array
        3D array with the data read from the Zarr files
    """


    tmp_tensor = reading_lara_variables_from_zarr(
        verbose,
        listdates[:2],
        lara_basedir,
        variables,
        gridded_variables,
        lon_left_lower_corner,
        lat_left_lower_corner,
        lon_right_upper_corner,
        lat_right_upper_corner,
        start_date,
        parcel_ids,
        allproc=True,
        rank=rank,
        from_https=from_https,
    )

    _, dimX, dimY = tmp_tensor.shape
    tensor_final = np.ones((len(listdates), dimX, dimY)) * (-999.9)
    tensor_final[-1, :, :] = tmp_tensor[0, :, :]
    tensor_final[-2, :, :] = tmp_tensor[1, :, :]

    listdates=listdates[2:]

    tmp_listdates = listdates[::-1]

    n = len(tmp_listdates)
    count = n // size
    remainder = n % size

    if rank < remainder:
        start = rank * (count + 1)
        stop = start + count + 1
    else:
        start = rank * count + remainder
        stop = start + count

    local_list = tmp_listdates[start:stop]

    local_results = np.empty((len(local_list), dimX, dimY))
    local_results = reading_lara_variables_from_zarr(
        verbose,
        local_list,
        lara_basedir,
        variables,
        gridded_variables,
        lon_left_lower_corner,
        lat_left_lower_corner,
        lon_right_upper_corner,
        lat_right_upper_corner,
        start_date,
        parcel_ids,
        allproc=False,
        rank=rank,
        from_https=from_https,
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

        tensor_final[int(i_start[0]) : int(i_stop[0]), :, :] = final_results
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
            tensor_final[int(i_start[i]) : int(i_stop[i]), :, :] = tmp
    comm.Bcast(tensor_final, root=0)

    tensor_final = compute_other_quantities(tensor_final, var_heat_track)

    if rank == 0:
        print(f"\n     => Data reading: Finished")

    gc.collect()
    return tensor_final





def get_parcels_over_target_region(
    variables, start_date, file_mask, maskname, maskvar_lon, maskvar_lat, mask_value, from_https, lara_basedir,
):

    """
    Get parcels over target region.

    Parameters
    ----------
    variables : list of str
        list of variables to read from LARA
    start_date : datetime
        start date of the simulation
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
    date_str = start_date.strftime("%Y-%m-%d %H:%M:%S")
    date_str_ = start_date.strftime("%Y%m%d")
    year = start_date.year
    month = start_date.month

    temp_url = get_lara_url(year, month)

    for var_name in variables:

        if from_https:
            zarr_url = f"{temp_url}/{var_name}"
            ds = xr.open_zarr(zarr_url)
        else:
            zarr_url = f"{lara_basedir}/{var_name}/{int(year)}/{str(int(month)).zfill(2)}/LARA_{date_str_}_{var_name}.nc"
            ds = xr.open_dataset(zarr_url)



        try:
            target_ds = ds.sel(time=date_str)
        except Exception as e:
            raise ValueError(f"    Error in file for {var_name} on date {date_str}")

        ds.close()

        parcels = target_ds["particle"].values
        data_values = target_ds[var_name].values

        if var_name == "lon":
            data_values = (data_values + 180) % 360 - 180

        vars()[var_name] = np.column_stack((parcels, data_values))
        vars()[var_name] = vars()[var_name][np.argsort(vars()[var_name][:, 0])]

    del target_ds
    del data_values
    gc.collect()

    data_temp = np.column_stack(
        (vars()["lon"][:, 0], vars()["lon"][:, 1], vars()["lat"][:, 1])
    )

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


def compute_other_quantities(tensor_final, var_heat_track):
    if var_heat_track == "dse":

        tensor_final[:, :, 11] = compute_dry_static_energy(
            tensor_final[:, :, 2], tensor_final[:, :, 9], tensor_final[:, :, 4]
        )
        tensor_final[:, :, 11] = tensor_final[:, :, 11] / 1000.0

    elif var_heat_track == "potTemp":

        tensor_final[:, :, 11] = compute_theta(
            tensor_final[:, :, 6],
            tensor_final[:, :, 3],
            tensor_final[:, :, 9],
            press_data=True,
            parpres=tensor_final[:, :, 13],
        )
    else:
        tensor_final[:, :, 11] = compute_theta(
            tensor_final[:, :, 6],
            tensor_final[:, :, 3],
            tensor_final[:, :, 9],
            press_data=True,
            parpres=tensor_final[:, :, 13],
        )

    tensor_final[:, :, 12] = compute_rh(
        tensor_final[:, :, 6],
        tensor_final[:, :, 3],
        tensor_final[:, :, 9],
        press_data=True,
        p_Pa=tensor_final[:, :, 13],
    )
    # tensor_final[:,:,13] = calc_pres(tensor_final[:,:,6], tensor_final[:,:,3], tensor_final[:,:,9])

    tensor_final[:, :, 11][tensor_final[:, :, 11] < 0] = -999.9
    tensor_final[:, :, 12][tensor_final[:, :, 12] < 0] = -999.9
    tensor_final[:, :, 13][tensor_final[:, :, 13] < 0] = -999.9

    tensor_final[tensor_final == -999.0] = -999.9

    return tensor_final


def get_vars_from_lara_zarr(
    verbose,
    listdates,
    lara_basedir,
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
    from_https=False,
    mode="backward",
):

    """
    Reads LARA data from Zarr files for a given list of dates and extracts variables 
    for parcels over a target region.

    Parameters
    ----------
    verbose : bool
        Print information about the process.
    listdates : list of datetime objects
        List of dates to read data from.
    lara_basedir : str
        Path to the LARA data directory.
    file_mask : str
        Path to the mask file.
    maskname : str
        Name of the mask variable in the mask file.
    maskvar_lon : str
        Name of the longitude variable in the mask file.
    maskvar_lat : str
        Name of the latitude variable in the mask file.
    lat : array like
        Array of latitudes of the target region.
    lon : array like
        Array of longitudes of the target region.
    rank : int
        Rank of the process.
    size : int
        Number of processes.
    comm : MPI communicator
        MPI communicator.
    lon_left_lower_corner : float
        Longitude of the lower left corner of the target region.
    lat_left_lower_corner : float
        Latitude of the lower left corner of the target region.
    lon_right_upper_corner : float
        Longitude of the upper right corner of the target region.
    lat_right_upper_corner : float
        Latitude of the upper right corner of the target region.
    mask_value : int
        Mask value for filtering parcels in the target region.
    var_heat_track : str
        Variable to track heat from.
    start_date : datetime object
        Start date of the simulation.
    from_https : bool
        If True, read data from the HTTPS server.
    mode : str
        Mode of the simulation. Can be "backward" or "forward".

    Returns
    -------
    tensor_final : array like
        Array of shape (num_parcels, num_steps, num_variables) with the variables
        for each parcel.

    """
    if rank == 0:
        if verbose:
            print("\n      --> Getting IDs of parcels over the target region")

    parcel_ids = get_parcels_over_target_region(
        variables=["lon", "lat"],
        start_date=start_date,
        file_mask=file_mask,
        maskname=maskname,
        maskvar_lon=maskvar_lon,
        maskvar_lat=maskvar_lat,
        mask_value=mask_value,
        from_https=from_https,
        lara_basedir=lara_basedir,
    )

    
    _, _, _, _, variables, gridded_variables = get_lara_config()

    comm.Barrier()

    if mode == "backward":
        tensor_final = get_lara_variables_from_zarr_backward(
            verbose,
            listdates,
            lara_basedir,
            variables,
            gridded_variables,
            lon_left_lower_corner,
            lat_left_lower_corner,
            lon_right_upper_corner,
            lat_right_upper_corner,
            start_date,
            parcel_ids,
            var_heat_track,
            from_https,
            rank,
            size,
            comm,
        )

    elif mode == "forward":
        tensor_final = get_lara_variables_from_zarr_forward(
            verbose,
            listdates,
            lara_basedir,
            variables,
            gridded_variables,
            lon_left_lower_corner,
            lat_left_lower_corner,
            lon_right_upper_corner,
            lat_right_upper_corner,
            start_date,
            parcel_ids,
            var_heat_track,
            from_https,
            rank,
            size,
            comm,
        )

    gc.collect()
    return tensor_final
