Input configuration file 
========================

.. code-block:: bash

    ===========================================================================================================
                                        - OPEN CONFIGURATION FILE -
    ===========================================================================================================
    Parameter                        Value                      Description
    verbose                       =  'True'/'False'             -> Print run logs
    runID                         =  "Experimental"             -> Run name

    ==========================================================================================
        - PATHS -
    ==========================================================================================
    raw_partposit_path            = "path"                      -> Directory containing input data
    basedir_lara                  = 'path'                      -> Directory to read or save LARA reanalysis files
    file_gz                       = 'True' / 'False'            -> Checking if model data is compresssed in gz format
    output_path                   = "path"                      -> Directory to save LATTIN outputs


    ==========================================================================================
        - MODEL DETAILS -
    ==========================================================================================
    model                           = 'FLEXPART'/'FLEXPART-WRF'  -> Lagrangian model
                                      'FLEXPART11'/'LARA' 
    check_lara_files                = 'True'/'False'             -> To check if LARA files exists
    lara_from_https                 = 'True'/'False'             -> To download LARA files if they are missing.
    total_emited_mass               = value                      -> Total emited mass in model simulation
    total_release_parcels           = value                      -> Total number of released parcels in model simulation

    ==========================================================================================
        - LATTIN RUN CONFIGURATION-
    ==========================================================================================
    mode                            = "backward"                  -> Run mode  
    year                            = int value  or list          -> Start year. E.g. year=2015 or year=[2015,2016,2017] 
    month                           = int value or list           -> Start month. E.g. month=1 or month=[1,7,3]. 
    day                             = int value or list           -> Start day. E.g. day=1 or day=[1,2,3]. 
    hour                            = int value or list           -> Start hour. E.g. hour=0 or day=[0,6,12]. 
    minutes                         = int value or list           -> Start minutes. E.g. minutes=0 or minutes=[0,10,20].
    ndays                           = int value                   -> Number of continuos days to start the simulation.
    time_step                       = int value                   -> Temporal resolution of input data [minutes]
    tracking_time                   = int value                   -> Total simulation time for tracking [minutes]
    calendar                        = '365d'/'366d'               -> Calendar type.
                                                                    Use calendar="365d" to discard Feb 29 in leap years

    lon_left_lower_corner           = value                       -> Domain limits for regional partposit files.
    lat_left_lower_corner           = value
    lon_right_upper_corner          = value
    lat_right_upper_corner          = value
    save_full_parts_position        = True / False                -> To save raw parcel trajectories

    ==========================================================================================
        - MASK FILE DETAILS -
    ==========================================================================================
    file_mask                       = 'path'                      -> Path to mask file (netcdf format)  
    maskname                        = 'mask'                      -> Name of mask variable in the mask file 
    maskvar_lat                     = 'lat'                       -> Latitude variable name in the mask file 
    maskvar_lon                     = 'lon'                       -> Longitude variable name in the mask file
    mask_value                      = value                       -> Mask value for filterirng parcels in the target region

    ==========================================================================================
        - OUTPUT DOMAIN RESOLUTION -
    ==========================================================================================
    resolution                      = 0.5                         -> Output resolutiom
    numPdX                          = 720                         -> Number of grid points in x-direction
    numPdY                          = 360                         -> Number of grid points in y-direction
    lon_lower_left                  = -180                        -> Longitude in lower left corner
    lat_lower_left                  = -90                         -> Latitude in lower left corner
    area_conserving_grid            = 'True'/'False'              -> To generate an area-conserving grid output

    ==========================================================================================
        - SPECIFIC FOR HEAT TRACKING -
    ==========================================================================================
    tracking_heat                   = 'True'/'False'              -> Activate heat tracking

    heat_tracking_method            = 'SCH19'                     -> Heat tracking method [SCH19, SCH20, JK22, CUSTOM].
                                                                    If you select one of this method [SCH19, SCH20, JK22],
                                                                    you do not to specify the next parameters.
                                                                    WARNING: The default values inly work for
                                                                    time_step=360 minutes

    var_heat_track                  = 'potTemp'/'dse'             -> Variable for heat tracking
                                                                    dse: Dry static energy
                                                                    potTemp: Potential Temperature  

    dvarheatthreshold               = value                       -> Minimun change in tracking var to be considered an uptake
                                                                    If tracking var is potential temperature,
                                                                    dvarheatthreshold is in Kelvin
                                                                    If tracking var is dry static energy,
                                                                    dvarheatthreshold is in kJ

    filter_pbl_parcels              = 'True'/'False'              -> Filter parcels within the target region within the PBL
    heat_custom_limits_highs        = [lower_limit, upper_limit]  -> Custom limits for filtering parcel within the target region [m]
                                                                    Set heat_custom_limits_highs = [0,0] to use PBL highs for filtering.
                                                                    Only it works if filter_pbl_parcels=True

    pblcheck                        = int value                    -> checking PBL condition along the parcels trajectories
                                                                    0: no PBL check, use everything
                                                                    1: at least one location within the PBL
                                                                    2: both locations within the PBL

    pbl_method                      = "maxval"                     -> PBL method for PBL check. [maxval, meanval, actualval] 
    trk_rh_check                    = 'True'/'False'               -> Check relative humidity
    rh_threshold                    = value                        -> Allowed relative humidity changes.
                                                                    Only needed if trk_rh_check=True

    dqcheck                         = 'True'/'False'               -> Checking changes in specific humidty along the parcels trajectory.
    dqthreshold                     = value                        -> Allowed changes in specific humidity.
                                                                    Only needed if dqcheck=True

    heat_linear_adjustment          = 'True'/'False'               -> Apply linear adjusment to detected uptakes

    save_heat_parts_position        = True / False                 -> To save processed parcels trajectories (heat)

    ==========================================================================================
        - SPECIFIC FOR MOISTURE TRACKING -
    ==========================================================================================
    tracking_moisture               = 'True' /'False'             -> Activate moisture tracking
    moisture_tracking_method        =  "SOD08"                    -> Misture tracking method [SOD08, SJ05, FAS19, JK22, APA22, APA25,CUSTOM]
                                                                    If you select one of this method[SOD08, SJ05, FAS19, JK22, APA22, APA25],
                                                                    you do not to specify the next parameters.
                                                                    WARNING: The default values only work for
                                                                    time_step=360 minutes


    filter_dqdt_parcels             = 'True' /'False'              -> Only track precipitating parcesl within the target region
    filter_pbl_dq_parcels           = 'True'/'False'               -> Filter parcels within the target region within the PBL
    moist_custom_limits_highs       = [lower_limit, upper_limit]   -> Custom limits for filtering parcel within the target region [m]
                                                                    Set moist_custom_limits_highs = [0,0] to use PBL highs for filtering.
                                                                    Only it works if filter_pbl_dq_parcels=True

    dqdt_threshold                  = value                        -> Change in specific humidity for considering that a precipitation
                                                                    event occurred within the target region.
                                                                    Only needed if filter_dqdt_parcels=True

    precip_minrh                    = value [float]                -> Minumim relative humidity to account for precipitation  [%]             
                                                                    Set precip_minrh=0 to do not apply

    dqpblcheck                      = value                        -> checking PBL condition along the parcels trajectories
                                                                    0: no PBL check, use everything
                                                                    1: at least one location within the PBL
                                                                    2: both locations within the PBL

    dqpbl_method                    = 'maxval'                     -> PBL method for PBL check [maxval, meanval, actualval] 
    trkdq_rh_check                  = 'True'/'False'               -> Check relative humidity 
    dqrh_threshold                  = value                        -> Allowed relative humidity changes
                                                                    Only needed if trkdq_rh_check=True
    mindq_gain                      = value [float]                -> Minimun change in specific humidity to be considered an uptake 
    check_RH_route_precip           = 'True' /'False'              -> To evaluate RH during moisture diagnostic                                                           
    precip_minrh_en_route           = value [float]                -> Minumim relative humidity to account for precipitation in route [%]  
    mindq_loss                      = value [float]                -> Minimum change in specific humidity to account for moisture loss [kg/kg]


    moisture_linear_adjustment      = 'True'/'False'               -> Apply linear adjusment to detected uptakes
    save_moisture_parts_position    = True / False                 -> To save processed parcels trajectories (moisture)

    ------------------------ Bias correction approach------------------------------------
    moist_bias_correction           = True / False                 -> To apply bias correction to Lagrangian precipitation estimate and moiture sources
    precip_fname_prefix             = value [str]                  -> Prefix of the precipitation filename
    precip_date_format              = value [str]                  -> Format of the date in the name of the precipitation file
                                                                    * "yyyymmddHM"  
                                                                    * 'yyyy-mm-dd H:M'
                                                                    * 'yyyymmddHM'
                                                                    * 'yyyy-mm-dd_H'
                                                                    * 'yyyymmdd_H'
                                                                    * "mmdd_H"
                                                                    * "mmddH"
                                                                    * "yyyymmdd_HM"

    precip_path                     = value                         -> Path to the precipitation file (netcdf format).
                                                                    It should contains the accumulated preciptation over a time period equal to `time_step`
    precip_lat                      = value                         -> Name of the latitude variable in the precipitation file
    precip_lon                      = value                         -> Name of the longitude variable in the precipitation file
    precip_var                      = value                         -> Name of the precipitation variable in the precipitation file

    ==========================================================================================
        - SPECIFIC FOR DRY INTRUSION ANALYSIS -
    ==========================================================================================
    DI_analysis                     = 'True'/'False'                -> Activate dry intrusion (DI) analysis
    DI_dt_thresh                    = value [float]                 -> Time period to check for DI [minutes]
    DI_step                         = value [float]                 -> Time step to define start backward time to checking DI. [minutes]
    DI_dp_change                    = value [float]                 -> Pressure change to account for DI  [Pa] 
    DI_start_pressure_level         = value [float]                 -> Parcel should be above this level at the start point to be considered for a Dry Intrusion [Pa] 
    DI_end_pressure_level           = value [float]                 -> Parcel should do downn this level to be considered for a Dry Intrusion (e.g., to PBL)  [Pa] 
    DI_time_limits_checking         = list [backward minutes]       -> Minimum and maximum value of track length to check for dry intrusion. 
                                                                    Negative value for backward tracking
    DI_moisture_tracking            = 'True' / 'False'              -> To account for moisture uptake of parcels with DI
    DI_save_raw_parcels             = 'True' / 'False'              -> To save processed parcels trajectories (DI) 

    ==========================================================================================
        - SPECIFIC FOR TEMPERATURE ANOMALY TRACKING -
    ==========================================================================================
    tracking_Tanom                    = 'True' / 'False'              -> Activating temperature anomaly tracking
    Tanom_tracking_method             = value [str]                   -> Temperature anomaly tracking method: ['RP23','PR23']
    path_clim_temperature             = Path [str]                    -> Path to netCDF files containing climatological temperature and other variables if needed
    climT_fname_prefix                = value [str]                   -> Prefix of the name of the climatological netCDf files
    climT_date_format                 = value [str]                   -> Format of the date in the name of the precipitation file
                                                                        * "yyyymmddHM"  
                                                                        * 'yyyy-mm-dd H:M'
                                                                        * 'yyyymmddHM'
                                                                        * 'yyyy-mm-dd_H'
                                                                        * 'yyyymmdd_H'
                                                                        * "mmdd_H"
                                                                        * "mmddH"
                                                                        * "yyyymmdd_HM"

    Tlat_var_name                     = value [str]                   -> Latitude variable name in the climatological temperature files
    Tlon_var_name                     = value [str]                   -> Longitude variable name in the climatological temperature files
    climTvar_name                     = value [str]                   -> Climatological air temperature variable name
    dTdp_var_name                     = value [str]                   -> dTdp variable name
    dTdt_var_name                     = value [str]                   -> dTdt variable name
    psfc_var_name                     = value [str]                   -> Surface pressure variable name
    Tvar_name                         = value [str]                   -> Air temperature variable name 
    Tplves_var_name                   = value [str]                   -> Pressure levels variable name 
    interpolate_parcel_temperature    = 'True'/ 'False'               -> To interpotate air temperature to parcel trajectories
    interpolate_sfc                   = 'True'/ 'False'               -> To interpotate surface pressure to parcel trajectories
    Tanom_threshold                   = value list[float]             -> Temperature anomaly threshold for filtering air parcels trajectories [K]
    analysis_levels                   = value list                    -> Atmopsheric levels to perform the analysis
                                                                        sfc: for surface
                                                                        pbl: for PBL
                                                                        Except for 'sfc' and 'pbl', levels should be in [hPa]
    dp_sfc                            = value [float]                 -> To filter air parcels close to surface [hPa]. 
                                                                        It will retain all parcels betweenn psfc and psfc - dp_sfc

    dp_upper                          = value [float]                 -> To filter air parcels at specific upper level. 
                                                                        It will retain all parcels between upper_level-dp_upper and upper_level+dp_upper
    Tanom_linear_adjustment           = 'True' / 'False'              -> Apply linear adjusment to temperature changes along parcel trajectories
                                                                        It is only required for PR23 method.
    save_Tanom_parts_position         = 'True' / 'False'              -> To save processed parcels trajectories (temperature anomaly tracking) 

    ===========================================================================================================
                                            - CLOSE CONFIGURATION FILE -
    ===========================================================================================================



.. note::
    
    Precipitation and climatological temperature files should be one for each time step (`time_step`) of the tracking period