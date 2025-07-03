import fnmatch
from netCDF4 import Dataset
from numpy import dtype
import matplotlib
import matplotlib.pylab as plt
import time
from pathlib import Path as check_PATH
import math
from mpi4py import MPI
import functools
import warnings
import sys
import numpy as np
import imp
import lattin.constants as lc
from lattin.lattin_functions import *

def lattin_main(pathfile):
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	size = comm.Get_size()
	
	
	
	if rank==0:
		start_time = time.time()
		disclaimer()
		print("Starting " + program_name() +" --->\n")
		print("Using parameters from: " + pathfile)
		print("------------------------------------------------------------------------------------------------------------\n")
	content = imp.load_source("", pathfile)
	verbose = check_paths(content, "verbose")
	runID = check_paths(content, "runID")
	#Reading paths 
	raw_partposit_path = check_paths(content, "raw_partposit_path")
	outputpath = check_paths(content, "output_path")
	file_gz = check_paths(content, "file_gz")
	
	#Reading model details
	model = check_paths(content, "model")
	total_emited_mass = check_paths(content, "total_emited_mass")
	total_release_parcels = check_paths(content, "total_release_parcels")
	
	#Reading Lattin run configuration
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
	
	#Reading mask details
	file_mask = check_paths(content, "file_mask")
	maskname = check_paths(content, "maskname")
	maskvar_lat = check_paths(content, "maskvar_lat")
	maskvar_lon = check_paths(content, "maskvar_lon")
	mask_value = check_paths(content, "mask_value")
	
	#Reading details for outpur domain		
	resolution = check_paths(content, "resolution")
	numPdX = check_paths(content, "numPdX")
	numPdY = check_paths(content, "numPdY")
	lon_lower_left = check_paths(content, "lon_lower_left")
	lat_lower_left = check_paths(content, "lat_lower_left")
	
	#Reading details for heat tracking
	tracking_heat= check_paths(content, "tracking_heat")
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
	heat_linear_adjustment=check_paths(content, "heat_linear_adjustment")
	heat_custom_limits_highs = check_paths(content, "heat_custom_limits_highs")
	save_heat_parts_position = check_paths(content, "save_heat_parts_position")
	
	#Reading details for moisture tracking
	tracking_moisture= check_paths(content, "tracking_moisture")
	moisture_tracking_method = check_paths(content, "moisture_tracking_method") 
	filter_dqdt_parcels = check_paths(content, "filter_dqdt_parcels")
	dqdt_threshold = check_paths(content, "dqdt_threshold")
	dqpblcheck = check_paths(content, "dqpblcheck")
	dqpbl_method = check_paths(content, "dqpbl_method")
	trkdq_rh_check = check_paths(content, "trkdq_rh_check")
	dqrh_threshold = check_paths(content, "dqrh_threshold")
	mindq_gain = check_paths(content, "mindq_gain")
	moisture_linear_adjustment=check_paths(content, "moisture_linear_adjustment")
	filter_pbl_dq_parcels=check_paths(content, "filter_pbl_dq_parcels")
	moist_custom_limits_highs=check_paths(content, "moist_custom_limits_highs")
	precip_minrh = check_paths(content, "precip_minrh")
	save_moist_parts_position = check_paths(content, "save_moisture_parts_position")


	verbose, file_gz, save_full_parts_position = check_init_parms(verbose, file_gz, save_full_parts_position,
																)

	verbose=str2boolean(verbose)
	file_gz=str2boolean(file_gz)
	save_full_parts_position = str2boolean(save_full_parts_position)
	
	
	filesperday=int(1440/dtime)
	density=total_emited_mass/total_release_parcels
	lag_times=np.arange(0,int((totaltime/dtime))+filesperday, filesperday)
	
	if runID=="" or runID==None:
		runID=program_name()+"_outputs/"
	
	
		
	#lon_lower_left=lon_lower_left-resolution/2
	#lat_lower_left=lat_lower_left - resolution/2
	
	
	lat,lon, cenlon= grid_point (resolution, numPdX, numPdY,lon_lower_left,lat_lower_left)
	area=calc_A(resolution, lat, lon)
	lat_plot,lon_plot= grid_plot_final(lat, lon)

	
	if calendar!="365d":
		calendar="366d"


	
	if tracking_heat:
		trk_rh_check,rh_threshold, pblcheck, dqcheck, pbl_method,var_heat_track,dqthreshold,filter_pbl_parcels,dvarheatthreshold,heat_linear_adjustment,heat_tracking_method, heat_custom_limits_highs, errors, errors_found, save_heat_parts_position = heat_tracking_parms(heat_tracking_method,
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
			if rank==0:
				print("Checking parameters for HEAT TRACKING: Errors detected. PLEASE TAKE ACTION!!!!\n" )
				time.sleep(0.5)
				print(errors)
				print("============================================================================================================")
				print(program_name() + " fatal error")
				print("Bye:)")
				print("============================================================================================================")
			raise SystemExit()
		else:
			if rank==0:
				print("Checking parameters for HEAT TRACKING: PASSED\n" )
	else:
		filter_pbl_parcels="False"
		trk_rh_check="False"
		heat_linear_adjustment="False"
		dqcheck="False"
		save_heat_parts_position = "False"
			
	if tracking_moisture:
		filter_dqdt_parcels, dqdt_threshold, mindq_gain, dqpbl_method, moisture_linear_adjustment, dqpblcheck, trkdq_rh_check, dqrh_threshold, moisture_tracking_method, filter_pbl_dq_parcels, 	moist_custom_limits_highs, precip_minrh,  errors, errors_found, save_moist_parts_position = moisture_tracking_parms(moisture_tracking_method, filter_dqdt_parcels, dqdt_threshold, mindq_gain, dqpbl_method, moisture_linear_adjustment,dqpblcheck, trkdq_rh_check, dqrh_threshold, dtime, filter_pbl_dq_parcels, moist_custom_limits_highs, precip_minrh, save_moist_parts_position)
	
	
		if errors_found:
			if rank==0:
				print("Checking parameters for MOISTURE TRACKING: Errors detected. PLEASE TAKE ACTION!!!!\n" )
				time.sleep(0.5)
				print(errors)
				print("============================================================================================================")
				print(program_name() + " fatal error")
				print("Bye:)")
				print("============================================================================================================")
			raise SystemExit()
		else:
			if rank==0:
				print("Checking parameters for MOISTURE TRACKING: PASSED\n" )
	else:
		filter_dqdt_parcels="False"
		trkdq_rh_check="False"
		moisture_linear_adjustment="False"
		save_moist_parts_position = "False"
	
	if not tracking_heat and not tracking_moisture:
		if rank==0:
			print("(**) ERROR: No tracking is active in the input file. Set heat_tracking=True or moisture_tracking=True or both = True\n")
		raise SystemExit()	
	
	
	filter_pbl_parcels=str2boolean(filter_pbl_parcels)
	trk_rh_check=str2boolean(trk_rh_check)
	filter_dqdt_parcels=str2boolean(filter_dqdt_parcels)
	trkdq_rh_check=str2boolean(trkdq_rh_check)
	dqcheck=str2boolean(dqcheck)
	heat_linear_adjustment=str2boolean(heat_linear_adjustment)
	moisture_linear_adjustment=str2boolean(moisture_linear_adjustment)
	
	save_heat_parts_position = str2boolean(save_heat_parts_position)
	save_moist_parts_position = str2boolean(save_moist_parts_position)
	
	errors, errors_found = checking_input_parameters(size,
												  lag_times[-1],
												  totaltime, 
												  dtime, 
												  file_mask,
												  model, 
												  mode, 
												  numPdX, 
												  numPdY, 
												  resolution, 
												  lat_plot[:,0],
												  lon_plot[0,:],
												  ndays,
												  cenlon,
												  )
	if errors_found:
		if rank==0:
			print("Checking config parameters: Errors detected. PLEASE TAKE ACTION!!!!\n" )
			time.sleep(0.5)
			print(errors)
			print("============================================================================================================")
			print(program_name() + " fatal error")
			print("Bye:)")
			print("============================================================================================================")
		raise SystemExit()
	else:
		if rank==0:
			print("Checking config parameters: PASSED\n" )
	
	
	if model in ["FLEXPART"]:
		type_file=2
	elif model in ("FLEXPART-WRF"):
		type_file=1


	if mode=="backward":
		track_days=np.arange(-1,-1*len(lag_times),-1)
	elif mode=="forward":
		track_days=np.arange(1,len(lag_times),1)



	list_year, list_month, list_day, list_hour, list_min=generate_simulation_dates(ndays, year, month, day, hour, minutes, calendar)
	
	
	output_path=outputpath+"/"+runID+"/"
	if rank==0:
		if not os.path.exists(output_path): os.makedirs(output_path)

	
	if rank==0:
		printting_run_information(verbose, raw_partposit_path, output_path, list_year, list_month, list_day, list_hour, list_min, ndays, model, mode, totaltime, dtime, var_heat_track,file_mask, heat_tracking_method, moisture_tracking_method, tracking_heat, tracking_moisture, density, dvarheatthreshold, pbl_method,pblcheck, dqcheck,dqthreshold, trk_rh_check, rh_threshold, filter_pbl_parcels, filter_dqdt_parcels, dqdt_threshold, dqpblcheck,mindq_gain, dqpbl_method, trkdq_rh_check, dqrh_threshold, heat_linear_adjustment, moisture_linear_adjustment,heat_custom_limits_highs, filter_pbl_dq_parcels, moist_custom_limits_highs,precip_minrh, calendar, runID, save_heat_parts_position, save_moist_parts_position, save_full_parts_position, cenlon)
	
	if rank==0:
		print("\n   -> "+program_name() + " is running with", size, "CPUs\n")
		
	
	

	
	
	
	for year, month, day, hour, minn in zip(list_year, list_month, list_day, list_hour, list_min):
		
		date=year+"-"+month+"-"+day+" "+hour+":"+minn+":00"
		filename_out="lattin_"+mode+"_"+year+month+day+hour+minn
		
		stats_fame="lattin_"+mode+"_stats_"+year+month+day+hour+minn+".dat"
		if rank==0:
			partial_start_time = time.time()
			print("")
			print("===> STARTING TRACKING FOR ", date)
			print("     -------------------------------------------------------------------------------------------------------")

		
		partpositfiles, listdates=generate_file(mode, dtime, totaltime, date, raw_partposit_path, file_gz, calendar)
		
		
		pferrors, find_error=checking_raw_partposti_files( partpositfiles)
		if find_error:
			if rank==0:
				print("     Checking raw partposit files: Files not found. PLEASE TAKE ACTION!!!!\n" )
				time.sleep(0.5)
				print(pferrors)
				print("=======================================================================================================")
				print(program_name() + " fatal error")
				print("Bye:)")
				print("=======================================================================================================")
			raise SystemExit()
		else:
			if rank==0:
				print("    Checking raw partposit files: PASSED\n" )
		
		
		
		ntimes=[]
		for i in range(0, len(listdates)):
			auxt=datetime(int(str(int(listdates[i]))[0:4]), int(str(int(listdates[i]))[4:6]), int(str(int(listdates[i]))[6:8]), int(str(int(listdates[i]))[8:10]), int(str(int(listdates[i]))[10:12]), int(str(int(listdates[i]))[12:14]) )
			ntimes=np.append(ntimes, auxt)
		
		meantimes=[]
		for i in range(0, len(ntimes)-1):
			meantimes=np.append(meantimes, ntimes[i] + (ntimes[i+1]-ntimes[i])/2)
		
		
		trackingtime_steps=np.arange(dtime, totaltime + 2*dtime, dtime)

		if rank==0:
			print("     !Getting data from raw partposit files")
		tensor_org = get_vars_from_partposit(verbose,partpositfiles ,file_mask, maskname,maskvar_lon, maskvar_lat,lat, lon, rank,size, comm, type_file,
               lon_left_lower_corner,lat_left_lower_corner, lon_right_upper_corner,lat_right_upper_corner, model, mask_value, file_gz, var_heat_track)
		
		
		#np.save("tensor_test.npy", tensor_org)
		#tensor_org=np.load("tensor_test.npy")
				
		parcels_count=len(tensor_org[0,:,0])
	
		
		if rank==0:
			print("\n     => Number of parcels within the target region:", parcels_count)
		
		
		#if rank==0:
		#Processing heat tracking #######################################
		if tracking_heat:
			if verbose and rank==0:
				print("\n     + PROCCESSING HEAT") 
			if mode=="backward":
				array_day_heat, tensor_heat, counter_part_heat, no_heatuptakes_parts=processing_heat_track_backward(tensor_org, pblcheck, filter_pbl_parcels, pbl_method, heat_custom_limits_highs, trk_rh_check, rh_threshold, dqcheck, dqthreshold, dvarheatthreshold, var_heat_track, lag_times, area, density, dtime, heat_linear_adjustment, lon,lat, numPdY,numPdX, cenlon, 0, rank,size, comm)
			
				if filter_pbl_parcels:
					if verbose and rank==0:
						print("      !Number of filter parcels within the PBL or custom highs at time t0:", counter_part_heat, "({:.2f}".format(100 * (counter_part_heat) / (parcels_count))+"%)")
						print("      !Number of parcels without heat uptake in the trajectory:", no_heatuptakes_parts, "({:.2f}".format(100 * (no_heatuptakes_parts) / (counter_part_heat)) + "%)" )
				
					save_heat_stats=True
				else:
					save_heat_stats=True
					if rank==0:
						print("      !Number of parcels without heat uptake in the trajectory:", no_heatuptakes_parts, "({:.2f}".format(100 * (no_heatuptakes_parts) / (counter_part_heat)) + "%)" )
		else:
			save_heat_stats=False
			counter_part_heat=None
			no_heatuptakes_parts=None
			array_day_heat=np.empty((len(track_days), lat.shape[0], lat.shape[1]))
			tensor_heat=np.empty_like(tensor_org[1:,:,:])
		#Processing moisture tracking #######################################
		if tracking_moisture:
			if verbose and rank==0:
				print("\n     + PROCCESSING MOISTURE")
			if mode=="backward":
				array_day_moist, tensor_moist, counter_precip_part, no_evap_uptakes, attributed_precip, CR, lwvrt = processing_moisture_track_backward(tensor_org, filter_dqdt_parcels, dqdt_threshold, 	filter_pbl_dq_parcels,moist_custom_limits_highs, dqpblcheck, 	dqpbl_method, trkdq_rh_check, dqrh_threshold, mindq_gain, lag_times, area, density, dtime,moisture_linear_adjustment,  precip_minrh, lon,lat, numPdY,numPdX, cenlon, 1, moisture_tracking_method, trackingtime_steps, rank,size, comm)
				if filter_dqdt_parcels==True or filter_pbl_dq_parcels==True:
					if verbose and rank==0:
						print("      !Number of precipitating parcels within the target region at time t0:", counter_precip_part, "({:.2f}".format(100 * (counter_precip_part) / (parcels_count))+"%)")
						
						print("      !Number of parcels without moisture uptake in the trajectory:", no_evap_uptakes, "({:.2f}".format(100 * (no_evap_uptakes) / (counter_precip_part))+"%)")
						if moisture_linear_adjustment:
							print("      !Lagrangian mean water vapour residence time: ", str(lwvrt)[0:5], "days")
						
					save_moist_stats=True
				else:
					save_moist_stats=True
					if rank==0:
						print("      !Number of parcels without moisture uptake in the trajectory:", no_evap_uptakes, "({:.2f}".format(100 * (no_evap_uptakes) / (counter_precip_part))+"%)")
			
		else:
			array_day_moist=np.empty((len(track_days), lat.shape[0], lat.shape[1]))
			tensor_moist=np.empty_like(tensor_org[:-1,:,:])
			CR = np.empty((len(track_days), lat.shape[0], lat.shape[1]))
			save_moist_stats=False
			save_attrib=False
			counter_precip_part=None
			no_evap_uptakes=None
			lwvrt=None
				
		if rank==0:		
			if verbose:
				print("\n     + SAVING TO")
				print("       !Output File:", output_path+"/"+filename_out+".nc")
			writing_netcdf(latitude=lat_plot[:,0],
							longitude=lon_plot[0,:],  
							tracking_heat=tracking_heat,  
							var_heat_track=var_heat_track, 
							array_heat_day=array_day_heat,
							tracking_moisture=tracking_moisture,
							array_moist_day=array_day_moist, 
							tensor_org=tensor_org, 
							tensor_heat=tensor_heat, 
							tensor_moist=tensor_moist,
							CR=CR,
							track_days=track_days,
							mode=mode,
							run_date=date,
							listdates=ntimes,
							meantimes=meantimes,
							dtime=dtime,
							moisture_tracking_method=moisture_tracking_method,
							heat_tracking_method=heat_tracking_method,
							moisture_linear_adjustment=moisture_linear_adjustment,
							filename=output_path+"/"+filename_out,
							save_heat_parts_position = save_heat_parts_position,
							save_moist_parts_position = save_moist_parts_position,
							save_full_parts_position = save_full_parts_position,
							)
			
			if verbose:
				print("\n     + SAVING STATS TO")
				print("       !Stats File:", output_path+"/"+stats_fame)
			partial_runtime = time.time() - partial_start_time
			saving_data(parcels_count,
				counter_part_heat,
				no_heatuptakes_parts,
				counter_precip_part,
				no_evap_uptakes,
				date,
				output_path+"/"+stats_fame,
				save_moist_stats,
				save_heat_stats,
				np.round(partial_runtime, 2),
				mode,
				heat_tracking_method,
				moisture_tracking_method,
				lwvrt,
				moisture_linear_adjustment
				)
			
			if verbose:
				print("\n     --------------------------------------------------------- ")
				print("     "+program_name() +" Version " +str(get_currentversion()) + " has successfully finished this date")
				print("     Partial Runtime: %.2f seconds." % np.round(partial_runtime, 2))
				print("     --------------------------------------------------------- ")
			
		time.sleep(0.05)
	if rank==0:
			
		runtime = time.time() - start_time	
		ending_credits(runtime)
