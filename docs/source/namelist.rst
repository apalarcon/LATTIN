Input configuration file 
========================

.. code-block:: bash

    ===========================================================================================================
                                            - OPEN CONFIGURATION FILE -
    ===========================================================================================================
    Parameter                        Value                      Description
    verbose                    =  'True'/'False'             -> Print run logs
    runID                      =  "Experimental"             -> Run name

    ==========================================================================================
        - PATHS -
    ==========================================================================================
    raw_partposit_path         = "path"                      -> Directory containing input data
    file_gz                    = 'True' / 'False'            -> Checking if model data is compresssed in gz format
    output_path                = "path"                      -> Directory to save LATTIN outputs


    ==========================================================================================
        - MODEL DETAILS -
    ==========================================================================================
    model                      = 'FLEXPART'/'FLEXPART-WRF'  -> Lagrangian model
    total_emited_mass          = value                      -> Total emited mass in model simulation
    total_release_parcels      = value                      -> Total number of released parcels in model simulation

    ==========================================================================================
        - LATTIN RUN CONFIGURATION-
    ==========================================================================================
    mode                       = "backward"                  -> Run mode  
    year                       = int value  or list          -> Start year. E.g. year=2015 or year=[2015,2016,2017] 
    month                      = int value or list           -> Start month. E.g. month=1 or month=[1,7,3]. 
    day                        = int value or list           -> Start day. E.g. day=1 or day=[1,2,3]. 
    hour                       = int value or list           -> Start hour. E.g. hour=0 or day=[0,6,12]. 
    minutes                    = int value or list           -> Start minutes. E.g. minutes=0 or minutes=[0,10,20].
    ndays                      = int value                   -> Number of continuos days to start the simulation.
    time_step                  = int value                   -> Temporal resolution of input data [minutes]
    tracking_time              = int value                   -> Total simulation time for tracking [minutes]
    calendar                   = '365d'/'366d'               -> Calendar type.
                                                                Use calendar="365d" to discard Feb 29 in leap years

    lon_left_lower_corner      = value                       -> Domain limits for regional partposit files.
    lat_left_lower_corner      = value
    lon_right_upper_corner     = value
    lat_right_upper_corner     = value
    save_full_parts_position   = True / False                -> To save raw parcel trajectories
    ==========================================================================================
        - MASK FILE DETAILS -
    ==========================================================================================
    file_mask                  = 'path'                      -> Path to mask file (netcdf format)  
    maskname                   = 'mask'                      -> Name of mask variable in the mask file 
    maskvar_lat                = 'lat'                       -> Latitude variable name in the mask file 
    maskvar_lon                = 'lon'                       -> Longitude variable name in the mask file
    mask_value                 = value                       -> Mask value for filterirng parcels in the target region

    ==========================================================================================
        - OUTPUT DOMAIN RESOLUTION -
    ==========================================================================================
    resolution                 = 0.5                         -> Output resolutiom
    numPdX                     = 720                         -> Number of grid points in x-direction
    numPdY                     = 360                         -> Number of grid points in y-direction
    lon_lower_left             = -180                        -> Longitude in lower left corner
    lat_lower_left             = -90                         -> Latitude in lower left corner

    ==========================================================================================
        - SPECIFIC FOR HEAT TRACKING -
    ==========================================================================================
    tracking_heat              = 'True'/'False'              -> Activate heat tracking

    heat_tracking_method       = 'SCH19'                     -> Heat tracking method [SCH19, SCH20, JK22, CUSTOM].
                                                                If you select one of this method [SCH19, SCH20, JK22],
                                                                you do not to specify the next parameters.
                                                                WARNING: The default values inly work for
                                                                time_step=360 minutes

    var_heat_track             = 'potTemp'/'dse'             -> Variable for heat tracking
                                                                dse: Dry static energy
                                                                potTemp: Potential Temperature  

    dvarheatthreshold          = value                       -> Minimun change in tracking var to be considered an uptake
                                                                If tracking var is potential temperature,
                                                                dvarheatthreshold is in Kelvin
                                                                If tracking var is dry static energy,
                                                                dvarheatthreshold is in kJ

    filter_pbl_parcels         = 'True'/'False'              -> Filter parcels within the target region within the PBL
    heat_custom_limits_highs   = [lower_limit, upper_limit]  -> Custom limits for filtering parcel within the target region [m]
                                                                Set heat_custom_limits_highs = [0,0] to use PBL highs for filtering.
                                                                Only it works if filter_pbl_parcels=True

    pblcheck                   =  int value                   -> checking PBL condition along the parcels trajectories
                                                                0: no PBL check, use everything
                                                                1: at least one location within the PBL
                                                                2: both locations within the PBL

    pbl_method                 = "maxval"                     -> PBL method for PBL check. [maxval, meanval, actualval] 
    trk_rh_check               = 'True'/'False'               -> Check relative humidity
    rh_threshold               = value                        -> Allowed relative humidity changes.
                                                                Only needed if trk_rh_check=True

    dqcheck                    = 'True'/'False'               -> Checking changes in specific humidty along the parcels trajectory.
    dqthreshold                = value                        -> Allowed changes in specific humidity.
                                                                Only needed if dqcheck=True

    heat_linear_adjustment     = 'True'/'False'               -> Apply linear adjusment to detected uptakes

    save_heat_parts_position   = True / False                 -> To save processed parcels trajectories (heat)
    ==========================================================================================
        - SPECIFIC FOR MOISTURE TRACKING -
    ==========================================================================================
    tracking_moisture           = 'True' /'False'             -> Activate moisture tracking
    moisture_tracking_method    =  "SOD08"                    -> Misture tracking method [SOD08, SJ05, FAS19, JK22, APA22, CUSTOM]
                                                                If you select one of this method[SOD08, SJ05, FAS19, JK22, APA22],
                                                                you do not to specify the next parameters.
                                                                WARNING: The default values only work for
                                                                time_step=360 minutes


    filter_dqdt_parcels         = 'True' /'True'              -> Only track precipitating parcesl within the target region
    filter_pbl_dq_parcels       = 'True'/'False'              -> Filter parcels within the target region within the PBL
    moist_custom_limits_highs   = [lower_limit, upper_limit]  -> Custom limits for filtering parcel within the target region [m]
                                                                Set moist_custom_limits_highs = [0,0] to use PBL highs for filtering.
                                                                Only it works if filter_pbl_dq_parcels=True

    dqdt_threshold              = value                        -> Change in specific humidity for considering that a precipitation
                                                                event occurred within the target region.
                                                                Only needed if filter_dqdt_parcels=True

    precip_minrh                = 80                           -> Minumim relative humidity to account for precipitation  [%]             
                                                                Set precip_minrh=0 to do not apply

    dqpblcheck                  = value                        -> checking PBL condition along the parcels trajectories
                                                                0: no PBL check, use everything
                                                                1: at least one location within the PBL
                                                                2: both locations within the PBL

    dqpbl_method                = 'maxval'                     -> PBL method for PBL check [maxval, meanval, actualval] 
    trkdq_rh_check              = 'True'/'False'               -> Check relative humidity 
    dqrh_threshold              = value                        -> Allowed relative humidity changes
                                                                Only needed if trkdq_rh_check=True
    mindq_gain                  = value                        -> Minimun change in specific humidity to be considered an uptake 
    moisture_linear_adjustment  = 'True'/'False'               -> Apply linear adjusment to detected uptakes
    save_moisture_parts_position  = True / False               -> To save processed parcels trajectories (moisture)
    ===========================================================================================================
                                            - CLOSE CONFIGURATION FILE -
    ===========================================================================================================