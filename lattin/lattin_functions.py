import numpy as np
import fnmatch
from netCDF4 import Dataset, date2num
from numpy import dtype
import matplotlib
import matplotlib.pylab as plt
import time
import struct
from pathlib import Path as check_PATH
import scipy.interpolate as interp
import time
import math
import gzip
import shutil
from mpi4py import MPI
import functools
import warnings
import sys
import os
from datetime import datetime, timedelta
from lattin.fmodules import read_binary_file as RBF
from lattin.fmodules import determined_id as D_id
from lattin.fmodules import search_row as sRow
from lattin.fmodules import len_file as lf
from lattin.fmodules import compute_grid_integrated_heat as compute_grid_integrated_heat
from lattin.fmodules import compute_grid_integrated_moist as compute_grid_integrated_moist
import lattin.constants as lc

warnings.filterwarnings("ignore")
print = functools.partial(print, flush=True)

##### HEADER FUNCTIONS #####################################
def program_name():
	return "LATTIN"

def program_fullname():
	return "Lagrangian Atmospheric moisTure and heaT trackINg"

def get_currentversion():
	pathpkg = os.path.dirname(__file__)
	version_file = pathpkg+"/VERSION"
	with open(version_file) as vfile:
		version = vfile.readlines()[0].strip()
	return(version)


def get_lsatupdate():
        lupathpkg = os.path.dirname(__file__)
        version_upd = lupathpkg+"/LAST_UPDATE"
        with open(version_upd) as ufile:
                uversion = ufile.readlines()[0].strip()
        return(uversion)


def disclaimer():
	print(
		"\n============================================================================================================"
	)
	print("||                    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  ~ ~ ~ _                                      ||")
	print("||                    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~                                     ||")
	print("||                    ++           ++++    +++++++++++ ++++++++++  [~~]  ++++++   ++++                    ||")
	print("||                    ++          ++  ++   +   ++    + +   ++   +  [~~]   ++ ++    ++                     ||")
	print("||                    ++         ++    ++      ++          ++      [~~]   ++  ++   ++                     ||")
	print("||                    ++        ++++++++++     ++          ++      [~~]   ++   ++  ++                     ||")
	print("||                    ++       ++        ++    ++          ++      [~~]   ++    ++ ++                     ||")
	print("||                    +++++++ ++          ++   ++          ++      [~~]  ++++   ++++++                    ||")
	print("||               <<======================================================================>>               ||")
	print("||                         " + program_fullname() +" (v" +str(get_currentversion())+")                     ||")
	print("||                                        Last Update :  " +  get_lsatupdate() +"                                       ||" )
	print("||                                             Copyright 2023                                             ||")
	print("||                                                                                                        ||")
	print("||                                                                                                        ||")
	print("||                " +program_name() + " Version " +str(get_currentversion())+ " is free under the terms of GNU General Public license V3           ||")
	print("||                                  EphysLab (Environmental Physics Laboratory)                           ||")
	print("||                                           Universidade de Vigo                                         ||")
	print("||                                 contact: albenis.perez.alarcon@uvigo.es                                ||")
	print("||                            " +program_name() +" is distributed WITHOUT ANY WARRANTY! (see LICENSE)                   ||")
	print(
		"============================================================================================================\n")



def ending_credits(runtime):
	print("\n\n============================================================================================================")
	print(program_name() +" Version " +str(get_currentversion()) + " has successfully finished")
	print("Runtime: %.2f seconds." % np.round(runtime, 2))
	print("Bye :)")
	print("============================================================================================================")


def print_error_message(message):
	print("\n=============================================================================================================")
	print ("ERROR: "+message)
	print("=============================================================================================================")
	raise SystemExit("Bye :)")
##### GENERAL FUNCTIONS #############################

def check_paths(pfile, path):
	try: 
		fpath = getattr(pfile, path)
	except:
		fpath = ""
	return fpath

def str2boolean(arg):
	if isinstance(arg, bool):
		return arg
	if arg.lower() in ("yes", "true", "t", "y", "1"):
		return True
	elif arg.lower() in ("no", "false", "f", "n", "0"):
		return False
	else:
		raise argparse.ArgumentTypeError("Boolean value expected.")


def heat_tracking_parms(heat_tracking_method,
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
						heat_custom_limits_highs
						):
	errors_found=False
	errors=""
	if heat_tracking_method.upper()=="SCH19" and dtime==360:
		trk_rh_check=False
		filter_pbl_parcels=True
		rh_threshold=0
		pblcheck=2
		dqcheck=True
		pbl_method="maxval"
		var_heat_track="dse"
		dqthreshold=0.1
		filter_pbl_parcels=True
		dvarheatthreshold=1
		heat_linear_adjustment=True
		heat_tracking_method=heat_tracking_method
		heat_custom_limits_highs=[0,0]
	elif heat_tracking_method.upper()=="SCH20"  and dtime==360:
		trk_rh_check=False
		filter_pbl_parcels=True
		rh_threshold=0
		pblcheck=1
		dqcheck=False
		pbl_method="maxval"
		var_heat_track="dse"	
		dqthreshold=1
		filter_pbl_parcels=True
		dvarheatthreshold=1
		heat_linear_adjustment=True
		heat_tracking_method=heat_tracking_method
		heat_custom_limits_highs=[0,0]
	elif heat_tracking_method.upper()=="JK22"  and dtime==360:
		trk_rh_check=True
		filter_pbl_parcels=True
		rh_threshold=10
		pblcheck=1
		dqcheck=False
		pbl_method="maxval"
		var_heat_track="potTemp"	
		dqthreshold=1
		filter_pbl_parcels=True
		dvarheatthreshold=0
		heat_linear_adjustment=True
		heat_tracking_method=heat_tracking_method
		heat_custom_limits_highs=[0,0]
	else:
		trk_rh_check=trk_rh_check
		filter_pbl_parcels=filter_pbl_parcels
		rh_threshold=rh_threshold
		pblcheck=pblcheck
		dqcheck=dqcheck
		pbl_method=pbl_method
		var_heat_track=var_heat_track	
		dqthreshold=dqthreshold
		filter_pbl_parcels=filter_pbl_parcels
		dvarheatthreshold=dvarheatthreshold
		heat_linear_adjustment=heat_linear_adjustment
		heat_tracking_method=heat_tracking_method
		heat_custom_limits_highs=heat_custom_limits_highs
		if heat_tracking_method in ("SCH19", "SCH20", "JK22"):
			print("")
			print ("RUN WARNING FOR HEAT TRACKING!!!!: Default values for "+ heat_tracking_method  + " only work for time_step = 360 minutes. Using CUSTOM tracking instead")
			print("")
			time.sleep(2)
			heat_tracking_method="CUSTOM"
		
		heat_message=" parameter is missing in the HEAT TRACKING module of the input file\n"
		heat_parms=["trk_rh_check", "rh_threshold", "filter_pbl_parcels", "pblcheck", "pbl_method", "var_heat_track", "dvarheatthreshold", "dqcheck", "dqthreshold", "heat_linear_adjustment"]
		
		if not heat_tracking_method in ("SCH19", "SCH20", "JK22"):
			for param in heat_parms:
				if vars()[param]=="" or vars()[param]==None:
					errors = errors + "(**) ERROR: " + param + heat_message
					errors_found=True

		if heat_custom_limits_highs==None or heat_custom_limits_highs=="" or heat_custom_limits_highs[0]>heat_custom_limits_highs[1]:
			errors = errors + "(**) ERROR: heat_custom_limits_highs is missing or lower_high in heat_custom_limits_highs parameter in input file is higher than upper_high"
			errors_found=True	
			
		
	return trk_rh_check,rh_threshold,pblcheck, dqcheck, pbl_method,var_heat_track,dqthreshold, filter_pbl_parcels,dvarheatthreshold,heat_linear_adjustment,heat_tracking_method,heat_custom_limits_highs, errors, errors_found


def moisture_tracking_parms(moisture_tracking_method, filter_dqdt_parcels, dqdt_threshold, mindq_gain, dqpbl_method, moisture_linear_adjustment,dqpblcheck, trkdq_rh_check, dqrh_threshold, dtime, filter_pbl_dq_parcels, 	moist_custom_limits_highs, precip_minrh):
	
	errors_found=False
	errors=""
	if moisture_tracking_method=="SOD08" and dtime==360:
		filter_dqdt_parcels=True
		dqdt_threshold=-0.0002
		mindq_gain=0.0002
		dqpbl_method="maxval"
		moisture_linear_adjustment=True
		dqpblcheck=1
		trkdq_rh_check=False
		dqrh_threshold=100
		moisture_tracking_method=moisture_tracking_method
		filter_pbl_dq_parcels=False
		moist_custom_limits_highs=[0,0]
		precip_minrh= 80
	elif moisture_tracking_method=="FAS19" and dtime==360:
		filter_dqdt_parcels=True
		dqdt_threshold=-0.0001
		mindq_gain=0.0001
		dqpbl_method="meanval"
		moisture_linear_adjustment=True
		dqpblcheck=0
		trkdq_rh_check=False
		dqrh_threshold=100
		moisture_tracking_method=moisture_tracking_method
		filter_pbl_dq_parcels=False
		moist_custom_limits_highs=[0,0]
		precip_minrh=80
	elif moisture_tracking_method=="SJ05" and dtime==360 :
		filter_dqdt_parcels=False
		dqdt_threshold=0
		mindq_gain=-10000000
		dqpbl_method="maxval"
		moisture_linear_adjustment=False
		dqpblcheck=0
		trkdq_rh_check=False
		dqrh_threshold=100
		moisture_tracking_method=moisture_tracking_method
		filter_pbl_dq_parcels=False
		moist_custom_limits_highs=[0,0]
		precip_minrh=0
	elif moisture_tracking_method=="JK22" and dtime==360:
		filter_dqdt_parcels=True
		dqdt_threshold=0
		mindq_gain=0
		dqpbl_method="maxval"
		moisture_linear_adjustment=True
		dqpblcheck=1
		trkdq_rh_check=True
		dqrh_threshold=20
		moisture_tracking_method=moisture_tracking_method
		filter_pbl_dq_parcels=False
		moist_custom_limits_highs=[0,0]
		precip_minrh=80
	elif moisture_tracking_method=="APA22" and dtime==360:
		filter_dqdt_parcels=True
		dqdt_threshold=-0.0001
		mindq_gain=0
		dqpbl_method="maxval"
		moisture_linear_adjustment=True
		dqpblcheck=0
		trkdq_rh_check=False
		dqrh_threshold=20
		moisture_tracking_method=moisture_tracking_method
		filter_pbl_dq_parcels=False
		moist_custom_limits_highs=[0,0]
		precip_minrh=0
	else:
		filter_dqdt_parcels=filter_dqdt_parcels
		dqdt_threshold=dqdt_threshold
		mindq_gain=mindq_gain
		dqpbl_method=dqpbl_method
		moisture_linear_adjustment=moisture_linear_adjustment
		dqpblcheck=dqpblcheck
		trkdq_rh_check=trkdq_rh_check
		dqrh_threshold=dqrh_threshold
		moisture_tracking_method=moisture_tracking_method
		filter_pbl_dq_parcels=filter_pbl_dq_parcels
		moist_custom_limits_highs=moist_custom_limits_highs
		
		if moisture_tracking_method in ("SOD08", "FAS19", "SJ05", "JK22", "APA22"):
			print("")
			print ("RUN WARNING FOR MOISTURE TRACKING!!!!: Default values for "+ moisture_tracking_method  + " only work for time_step = 360 minutes. Using CUSTOM tracking instead")
			print("")
			time.sleep(2)
			moisture_tracking_method="CUSTOM"

		moisture_message=" parameter is missing in the MOISTURE TRACKING module of the input file\n"
		moist_parms=["filter_dqdt_parcels", "dqdt_threshold", "dqpblcheck", "dqpbl_method", "trkdq_rh_check", "dqrh_threshold", "moisture_linear_adjustment","filter_pbl_dq_parcels","precip_minrh"]
		if not moisture_tracking_method in ("SOD08", "FAS19", "SJ05", "JK22", "APA22"):
			for param in moist_parms:
				if vars()[param]=="" or vars()[param]==None:
					errors = errors + "(**) ERROR: " + param + moisture_message
					errors_found=True
		if moist_custom_limits_highs==None or moist_custom_limits_highs=="" or moist_custom_limits_highs[0]>moist_custom_limits_highs[1]:
			errors = errors + "(**) ERROR: moist_custom_limits_highs is missing or lower_high in moist_custom_limits_highs parameter in input file is higher than upper_high"
			errors_found=True

	
	return filter_dqdt_parcels, dqdt_threshold, mindq_gain, dqpbl_method, moisture_linear_adjustment,dqpblcheck, trkdq_rh_check, dqrh_threshold,moisture_tracking_method, filter_pbl_dq_parcels, 	moist_custom_limits_highs, precip_minrh, errors, errors_found



def checking_input_parameters(size,
							totalfiles,
							totaltime,
							dtime, 
							file_mask,
							model, 
							mode, 
							numPdx,
							numPdy, 
							resolution,
							latgrid, 
							longrid,
							ndays,
							):
	errors=""
	errors_found=False
	if (totalfiles)%size!=0:
		errors="(**) ERROR: The number of processors must exactly divide the number of partposit files to process (trackin_time/time_step).\n     From your input file, the recommended number of processors is " + str(int(totaltime/dtime))+"\n"
		errors_found=True
	
		
	try:
		nc_=Dataset(file_mask)
	except:
		errors = errors + "(**) ERROR: No such file or directory: "+ file_mask +"\n"
		errors_found=True
	
	if not model in ("FLEXPART", "FLEXPART-WRF"):
		errors = errors + "(**) ERROR: Model must be FLEXPART or FLEXPART-WRF\n"
		errors_found=True
	if not mode in ("backward"):
		errors = errors + "(**) ERROR: mode must be 'backward' in time\n"
		errors_found=True
	
	if numPdx <=0  or  numPdy<=0:
		errors = errors + "(**) ERROR: numPdX  and  numPdY  must be greater than zero\n"
		errors_found=True
	if resolution<=0:
		errors = errors + "(**) ERROR:  resolution must be greater than zero \n"
		errors_found=True
	
	if latgrid.min()<-90 or latgrid.max()>90:
		errors = errors + "(**) ERROR:  Latitudes for output grid are incorrect \n"
		errors_found=True
	if longrid.min()<-180 or longrid.max()>180:
		errors = errors + "(**) ERROR:  Longitudes for output grid are incorrect \n"
		errors_found=True
	if ndays<1:
		errors = errors + "(**) ERROR:  ndays = " + str(int(ndays)) + " is incorrect. Parameter ndays must be higher than 0  \n"
		errors_found=True
			
	return errors, errors_found


def checking_raw_partposti_files( partpositfiles):
	errors=""
	find_error=False
	for i in partpositfiles:
		my_file=check_PATH(i)
		if not my_file.is_file():
			errors = errors + "     (**) ERROR: No such file or directory: " +i +"\n"
			find_error=True

	return errors, find_error



def printting_run_information(verbose,
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
							  dqtrk_rh_check, 
							  dqrh_threshold,
							  heat_linear_adjustment, 
							  moisture_linear_adjustment,
							  heat_custom_limits_highs,
							  filter_pbl_dq_parcels, 
							  moist_custom_limits_highs,
							  precip_minrh,
							  calendar,
							  ):
	
	#print("-------------------------------------------------------------------------------------------------------------\n")
	if verbose:
		print("RUNNING PARAMETERS")
		print("-------------------------------------------------------------------------------------------------------------")
		print("   -> Model                                      :", model)
		print("   -> Mask file                                  :", file_mask)
		print("   -> Raw partposit data                         :", raw_partposit_path)
		print("   -> Output directory                           :", output_path)
		print("   -> Tracking mode                              :", mode, "in time")
		print("   -> Time step                                  :", dtime, "minutes")
		print("   -> Tracking time                              :", totaltime, "minutes", "(" +str(totaltime/1440) + " days)")
		print("   -> Heat tracking                              :", tracking_heat)
		print("   -> Moisture tracking                          :", tracking_moisture)
		print("   -> Simulation starts at                       :", list_year[0]+list_month[0]+list_day[0]+" "+list_hour[0]+":"+list_min[0]+":00")
		print("   -> Simulation ends at                         :", list_year[-1]+list_month[-1]+list_day[-1]+" "+list_hour[-1]+":"+list_min[-1]+":00")
		print("   -> Calendar                                   :", calendar)

		if tracking_heat:
			
			if var_heat_track=="dse":
				varunits="kJ/Kg"
				varname="Dry static energy"
			elif var_heat_track=="potTemp":
				varunits="K"
				varname="Potential Temperature"
				
			if pblcheck==0:
				pblcheck_="no PBL check, use everything"
			elif pblcheck==1:
				pblcheck_="at least one location within the PBL"
			elif pblcheck==2:
				pblcheck_="both locations within the PBL"
				
			if heat_custom_limits_highs[0]==0 and heat_custom_limits_highs[1]==0:
				filterhigh="(Not Apply)"
			else:
				filterhigh=""

			print("\n   -------------    Heat Tracking Information    -------------")
			print("    + Tracking method                            :", heat_tracking_method)
			print("    + Heat tracking using                        :", varname)
			print("    + Uptake threshold                           :", ">", dvarheatthreshold, varunits)
			print("    + Filter parcels within PBL                  :", filter_pbl_parcels)
			print("    + PBL check                                  :", pblcheck_)
			if pblcheck!=0:
				print("    + PBL method                                 :", pbl_method)
			print("    + Check specific humidity                    :", dqcheck)
			if dqcheck:
				print("    + Specific humidity chage                    :", "<", dqthreshold)
			print("    + Check relative humidity                    :", trk_rh_check)
			if trk_rh_check:
				print("    + Relative humidity change                   :", "<=", rh_threshold, "%")
			print("    + Linear adjustment                          :", heat_linear_adjustment)
			print("    + Lower and upper limits for filter parcels  :", heat_custom_limits_highs, 'meters', filterhigh)
				
		if tracking_moisture:		
			print("\n   -------------  Moisture Tracking Information  -------------")
			print("    + Tracking method                            :", moisture_tracking_method)
			print("    + Moisture tracking using                    :", "Specific Humidty")
			print("    + Filter precipitating parcels               :", filter_dqdt_parcels)
			if filter_dqdt_parcels:
				print("    + dq/dt threshold                            :", dqdt_threshold, "kg/kg")
			print("    + Filter parcels within PBL                  :", filter_pbl_dq_parcels) 
			
			if moist_custom_limits_highs[0]==0 and moist_custom_limits_highs[1]==0:
				mfilterhigh="(Not Apply)"
			else:
				mfilterhigh=""
			if filter_pbl_dq_parcels:
				print("    + Lower and upper limits for filter parcels  :", moist_custom_limits_highs, 'meters', mfilterhigh )	
				
			if dqpblcheck==0:
				dqpblcheck_="no PBL check, use everything"
			elif dqpblcheck==1:
				dqpblcheck_="at least one location within the PBL"
			elif dqpblcheck==2:
				dqpblcheck_="both locations within the PBL"
			print("    + PBL check                                  :", dqpblcheck_)
			if dqpblcheck!=0:
				print("    + PBL method                                 :", dqpbl_method)
			
			print("    + Uptake threshold                           :", ">", mindq_gain)
			print("    + Check relative humidity                    :", dqtrk_rh_check)
			if dqtrk_rh_check:
				print("    + Relative humidity change                    :", "<=", dqrh_threshold, "%")
			print("    + Linear adjustment                          :", moisture_linear_adjustment)
			print("    + Minimum RH to account for precipitation    :", ">",precip_minrh, "%")
		print("-------------------------------------------------------------------------------------------------------------")
		
		
def writing_netcdf(latitude,
				   longitude, 
				   tracking_heat, 
				   var_heat_track, 
				   array_heat_day,
				   tracking_moisture, 
				   array_moist_day, 
				   tensor_org, 
				   tensor_heat, 
				   tensor_moist,
				   CR,
				   track_days,
				   mode, 
				   run_date,
				   listdates,
				   meantimes,
				   dtime,
				   moisture_tracking_method,
					heat_tracking_method,
					moisture_linear_adjustment,
				   filename
				   ):
	
	try:
		os.remove(filename+".nc")
	except OSError:
		pass
	
	ncout = Dataset(filename+".nc", 'w', format='NETCDF4')
	# define axis size
	ncout.createDimension('lat', len(latitude))
	ncout.createDimension('lon', len(longitude))
	ncout.createDimension('tracking_days', len(track_days))
	ncout.createDimension("time", 1)
	
	ncout.createDimension('ntime', tensor_org.shape[0])
	ncout.createDimension('nparcels', tensor_org.shape[1])
	ncout.createDimension('nvars', tensor_org.shape[2])


	time = ncout.createVariable("time", "f8", "time")
	time.units = "hours since 1900-01-01 00:00:00"
	time.calendar = "Standard"
	time.standard_name = "Run started at " + str(run_date)
	
	
	tracking_days = ncout.createVariable("tracking_days", "f4", "tracking_days")
	tracking_days.units = "days"
	

	ntime = ncout.createVariable("ntime", "f8", "ntime")
	ntime.units = "hours since 1900-01-01 00:00:00"
	ntime.calendar = "Standard"
	ntime.standard_name = "Run started at " + str(run_date)


	# create latitude axis
	lat = ncout.createVariable('lat', 'f8', ('lat'),zlib=True)
	lat.standard_name = 'latitude'
	lat.long_name = 'latitude'
	lat.units = 'degrees'
	lat.axis = 'Y'

	# create longitude axis
	lon = ncout.createVariable('lon', 'f8', ('lon'),zlib=True)
	lon.standard_name = 'longitude'
	lon.long_name = 'longitude'
	lon.units = 'degrees'
	lon.axis = 'X'


	parcels_position = ncout.createVariable('raw_parcels_postion', 'f8', ("time","ntime",'nparcels', "nvars"),zlib=True)
	parcels_position.standard_name = 'Raw parcels position at each time step. Only parcels within the target region at time t0'




	if tracking_heat:
		
		ncout.createDimension('heat_time', tensor_heat.shape[0])
		ncout.createDimension('heat_parcels', tensor_heat.shape[1])
		ncout.createDimension('heat_vars', 3)
		
		
		heat_time = ncout.createVariable("heat_time", "f8", "heat_time")
		heat_time.units = "hours since 1900-01-01 00:00:00"
		heat_time.calendar = "Standard"
	
		
		
		heatd = ncout.createVariable('Heat_days', 'f8', ("time",'tracking_days','lat', "lon"),zlib=True)
		heatd.standard_name = 'Surface Sensible Heat Flux by days'
		heatd.units ="W/m2"
		


		heat = ncout.createVariable('Heat_integrated', 'f8', ("time",'lat', "lon"),zlib=True)
		heat.standard_name = 'Integrated Surface Sensible Heat Flux'
		heat.units ="W/m2"
		
		heat_parcels = ncout.createVariable('heat_parcels_position', 'f8', ("time","heat_time",'heat_parcels', "heat_vars"),zlib=True)
		heat_parcels.standard_name = 'Parcels position for heat tracking'
		if var_heat_track=="sde":
			heat_parcels.units ="J/kg per " + str(int(dtime)) + " minutes"
		elif var_heat_track=="potTemp":
			heat_parcels.units ="K per " + str(int(dtime)) + " minutes" 
	
	if tracking_moisture:
		
		ncout.createDimension('moisture_time', tensor_moist.shape[0])
		ncout.createDimension('moisture_parcels', tensor_moist.shape[1])
		ncout.createDimension('moisture_vars', 3)
		
		moisture_time = ncout.createVariable("moisture_time", "f8", "moisture_time")
		moisture_time.units = "hours since 1900-01-01 00:00:00"
		moisture_time.calendar = "Standard"
		
		
		moistd = ncout.createVariable('moisture_days', 'f8', ("time",'tracking_days','lat', "lon"),zlib=True)
		moistd.standard_name = 'Evaporation minus precipitation by days'
		moistd.units ="mm/day"

		moist = ncout.createVariable('moisture_integrated', 'f8', ("time",'lat', "lon"),zlib=True)
		moist.standard_name = 'Integrated Surface Sensible Heat Flux'
		moist.units ="mm/day"
		
		moist_parcels = ncout.createVariable('moisture_parcels_position', 'f8', ("time","moisture_time",'moisture_parcels', "moisture_vars"),zlib=True)
		moist_parcels.long_name = 'Parcels position for moisture tracking'
		moist_parcels.units ="kg/kg per " + str(int(dtime)) + " minutes"
	
		#if moisture_linear_adjustment:
		#	moistCRd = ncout.createVariable('moisture_contribution_days', 'f8', ("time",'tracking_days','lat', "lon"),zlib=True)
		#	moistCRd.standard_name = 'Moisture contribution in for each tracking day'
		#	moistCRd.units ="%"

		#	moistCR = ncout.createVariable('integrated_moisture_contribution', 'f8', ("time",'lat', "lon"),zlib=True)
		#	moistCR.standard_name = 'Integrated moisture contribution in for each tracking day'
		#	moistCR.units ="%"


	
	lon[:] = longitude[:]
	lat[:]= latitude[:]
	parcels_position[0,:]=tensor_org[:]
	time[:]=date2num(datetime(int(run_date.split(" ")[0].split("-")[0]),int(run_date.split(" ")[0].split("-")[1]),int(run_date.split(" ")[0].split("-")[2]),int(run_date.split(" ")[1].split(":")[0]),int(run_date.split(" ")[1].split(":")[1]),0   ), time.units, time.calendar)
	
	tracking_days[:]=track_days[:]
	
	
	ntime[:]=date2num(listdates, ntime.units, ntime.calendar)
	
	tensor_description="nvars[0]-particle id, nvars[1]-longitude, nvars[2], latitude, nvars[3]-specific humidity(kg/kg), nvars[4]-parcel high (m), nvars[5]-topography high [m], nvar[6]-parcel density (kg/m3), nvar[7]-PBL high [m], nvar[8]-tropopause high [m], nvar[9]-parcel Temperature [K], nvar[10]-parcel mass [kg], nvar[11]-heat tracking var, nvar[12]-reltive humidity[%]"
	
	if tracking_heat:
		heatd[0,:]=array_heat_day[:]
		heat[0,:]=np.sum(array_heat_day, axis=0)
		heat_parcels[0,:]=tensor_heat[:,:,:3]
		heat_time[:]=date2num(meantimes[1:], heat_time.units, heat_time.calendar)
		
		heat_description="Heat Tracking using " + heat_tracking_method + " method. heat_vars[0] = longitude, heat_vars[1] = latitude, heat_vars[2] = d"+var_heat_track+"/dt" 
	else:
		heat_description="Heat Tracking no applied"
		
		
	if tracking_moisture:
		moistd[0,:]=array_moist_day[:]
		moist[0,:]=np.sum(array_moist_day, axis=0)
		moist_parcels[0,:]=tensor_moist[:,:,:3]

		#if moisture_linear_adjustment:
		#	moistCRd[0,:] = CR[:]
		#	moistCR[0,:] = np.sum(CR, axis=0)



		moisture_time=date2num(meantimes[:-1], moisture_time.units, moisture_time.calendar)

		moist_description="Moisture Tracking using " + heat_tracking_method + " method. moisture_vars[0] = longitude, moisture_vars[1] = latitude, moisture_vars[2] = dq/dt"
	else:
		moist_description="Moisture Tracking no applied"
	
	
	
	ncout.title = program_fullname()
	ncout.description =tensor_description + ". <====> " +  heat_description + ". <====> " + moist_description
	today = datetime.now()
	ncout.history = "Created " + today.strftime("%d/%m/%Y %H:%M:%S") + " using " + program_name()
	ncout.institution = (
		"Environmental Physics Laboratory (EPhysLab), University of Vigo, Spain"
	)
	ncout.source = (
		program_name()+" "
		+ str(get_currentversion())
	)
		
	ncout.close()


def saving_data(count_parcels,
		heat_parcels,
		no_heat_uptake_parcels,
		precipitating_parcels,
		no_evap_uptakes,
		rundates,
		fname,
		save_moist_stats,
		save_heat_stats,
		):
	f=open(fname, "w")
	f.write("Run Date: " + rundates+"\n")
	f.write("Number of parcels within the target region: " + str(int(count_parcels))+"\n")
	if save_heat_stats:
		f.write("Number of filter parcels within the PBL at time t0 for heat tracking: " + str(int(heat_parcels))+ " ({:.2f}".format(100 * (heat_parcels) / (count_parcels))+"%)\n")	
		f.write("Number of parcels without heat uptake in the trajectory: " + str(int(no_heat_uptake_parcels))+ " ({:.2f}".format(100 * (no_heat_uptake_parcels) / (heat_parcels)) + "%)\n")	
	if save_moist_stats:
		f.write("Number of precipitating parcels within the target region at time t0: " + str(int(precipitating_parcels))+ " ({:.2f}".format(100 * (precipitating_parcels) / (count_parcels)) + "%)\n")
		f.write("!Number of parcels without moisture uptake in the trajectory: " + str(int(no_evap_uptakes)) + " ({:.2f}".format(100 * (no_evap_uptakes) / (precipitating_parcels))+"%)\n")
	
		
def calc_A(resolution, lat, lon):

	rt = lc.earth_radius
	gr = np.pi/180.
	a,b=lat.shape
	area=np.empty((a-1,b-1))
	area[:,:]=0
	for j in range(len(lat[0,:])-1):
		for i in range(len(lat[:,0])-1):
			area[i,j]=np.abs((gr*rt**2)*( np.sin(gr*lat[i,j]) - np.sin(gr*lat[i+1,j])))*np.abs(resolution)
	return area

def grid_point (resolution, numPdX, numPdY, lon_lower_left,lat_lower_left):

	lat_new=[]
	lon_new=[]
	lat_min=lat_lower_left
	lon_min=lon_lower_left
	lat_new=np.append(lat_new, lat_min)
	lon_new=np.append(lon_new, lon_min)

	for i in range(numPdY):
		lat_min=lat_min+resolution
		lat_new=np.append(lat_new, lat_min)
	for j in range(numPdX):
		lon_min=lon_min+resolution
		lon_new= np.append(lon_new, lon_min)
	lon, lat=np.meshgrid(lon_new, lat_new)
	return lat, lon

def grid_plot_final(lat, lon):

	lat_new=[]
	lon_new=[]
	for i in range(len(lat[:,0])-1):
		lat_new= np.append(lat_new, (lat[i+1,0]+ lat[i,0])/2.)
	for j in range(len(lon[0,:])-1):
		lon_new= np.append(lon_new, (lon[0,j+1]+ lon[0,j])/2.)
	lon_plot, lat_plot=np.meshgrid(lon_new, lat_new)
	return lat_plot,lon_plot


def generate_simulation_dates(ndays, cyear, cmonth, cday, chours, cminutes, calendar="366d" ):

	if not isinstance(chours, list):
		chours=[chours]
	if not isinstance(cminutes, list):
		cminutes=[cminutes]
	if not isinstance(cyear, list):
		cyear=[cyear]
	if not isinstance(cmonth, list):
		cmonth=[cmonth]       
	if not isinstance(cday, list):
		cday=[cday]  

	
	if calendar=="365d":
	
		date_start=datetime(int(cyear[0]), int(cmonth[0]), int(cday[0]))
		
		end_date = str(time_calc(str(int(cyear[0]))+"-"+str(int(cmonth[0])).zfill(2)+"-"+str(int(cday[0])).zfill(2) + " " + str(int(chours[0])).zfill(2)+":"+str(int(cminutes[0])).zfill(2)+":00",int(ndays*24))).split(" ")[0]


		date_end = datetime(int(end_date.split("-")[0]), int(end_date.split("-")[1]), int(end_date.split("-")[2]))

		date_list = [date_start.date() + timedelta(days=x) for x in range(ndays+1)]
		leapyears=0
		for tdate in date_list:
					
			if tdate.year%4==0 and tdate.month==2 and tdate.day==29:
				leapyears+=1
		
		
		nhour=int((ndays+leapyears)*24)
	else:
		nhour=int(ndays*24)
	year=[]
	mes=[]
	dia=[]
	hora=[]
	mins=[]
	array =np.arange(0,nhour,24)

	  


	
		
	for yy in cyear:
		yy=str(int(yy)).zfill(4)
		for i in array:
			for mm in cmonth:
				mm=str(int(mm)).zfill(2)
				for dd in cday:
					dd=str(int(dd)).zfill(2)
					for hh in chours:
						for mmin in cminutes:
							fecha=yy+"-"+mm+"-"+dd+" "+str(int(hh)).zfill(2)+":"+str(int(mmin)).zfill(2)+":00"
							a=str(time_calc(fecha,float(i)))
							var1=a.split(" ")
							var11=var1[0].split("-")
							var12=var1[1].split(":")
							year_=str(var11[0])
							mes_=str(var11[1])
							dia_=str(var11[2])
							hora_=str(var12[0])
							minn_=str(var12[1])
							
							if calendar=="365d" and int(year_)%4==0 and int(mes_)==2 and int(dia_)==29:
								msg="nothing to do"
							else:
								year.append(year_)
								mes.append(mes_)
								dia.append(dia_)
								hora.append(hora_)
								mins.append(minn_)

	return year, mes, dia, hora,mins

def read_binaryFile_fortran(filename, type_file,lon_left_lower_corner,lat_left_lower_corner,lon_right_upper_corner,lat_right_upper_corner):

    if type_file==1:
        with open(filename,'rb') as inputfile:
            a=b''.join([line for line in inputfile])
        npart=struct.unpack('iiii', a[0:16])
        npart=npart[2]
        data= RBF(filename,npart,lon_left_lower_corner,lat_left_lower_corner, lon_right_upper_corner,lat_right_upper_corner)

    if type_file==2:
        len_a=lf(filename)
        npart=((len_a-12)/60)-1
        data= RBF(filename,npart, lon_left_lower_corner,lat_left_lower_corner, lon_right_upper_corner,lat_right_upper_corner)
    ind=np.where(data[:, 0]==-999)
    data=data[:int(ind[0][0]), :]

    return data

def load_mask_grid_NR(filename, maskname,maskvar_lon, maskvar_lat):

    wrfile = Dataset(filename)
    lat  = wrfile.variables[maskvar_lat][:]
    lon  = wrfile.variables[maskvar_lon][:]
    mask  = wrfile.variables[maskname][:]

    if len(lon.shape)<2:
        lon, lat=np.meshgrid(lon,lat)
    
    
    
    for i in range(0,lon.shape[0]):
       for j in range(0,lon.shape[1]):
          if lon[i,j]>180:
                lon[i,j]=lon[i,j]-360

    return lat, lon,mask

def funtion_interpol_mascara (maskvar_lat, maskvar_lon, mascara, data):

    lat_lon=np.empty((len(maskvar_lat), 2))
    lat_lon[:,0]=maskvar_lon
    lat_lon[:,1]=maskvar_lat

    prsInterpu = interp.NearestNDInterpolator(lat_lon,mascara)
    si = np.empty((data[:,1].size, 2))
    si[:,0] = data[:,1]
    si[:,1] = data[:,2]
    result=prsInterpu(si)
    
    return result

def desc_gz(name_file):

    with gzip.open(name_file, 'rb') as f_in:
        with open(name_file[:-3], 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
            
            
def determine_id_binary_grid_NR_fortran(data, maskvar_lat, maskvar_lon, value_mascara, value_mask):

    part_inmask=funtion_interpol_mascara (maskvar_lat, maskvar_lon, value_mascara, data)
    
       
    id_vector=np.array(D_id( part_inmask, value_mask, len( part_inmask)),dtype=int)
    
    submatrix=[]
    ind=[]
    for ii in id_vector:
        if ii !=-999:
            submatrix=np.append(submatrix, data[ii,:])
            ind.append(ii)
    submatrix=np.reshape(submatrix,(len(ind), 11))
    return submatrix

def search_row_fortran(lista, matrix):

    matrix_=np.array(sRow(matrix, lista, len(lista), len(matrix[:,0])), np.float64)

    return matrix_


def time_calc(init_time,h_diff):

	formatted_time = datetime.strptime(init_time, "%Y-%m-%d %H:%M:%S")
	calculated_time=formatted_time+timedelta(hours=h_diff)
	return calculated_time

def time_calcminutes(init_time,h_diff):

	formatted_time = datetime.strptime(init_time, "%Y-%m-%d %H:%M:%S")
	calculated_time=formatted_time+timedelta(minutes=h_diff)
	return calculated_time

def generate_file(mode, dtime, totaltime, fecha, path, key_gz, calendar):
  
	list_fecha=[]
	listdates=[]
	if mode == "backward":
		
		leapyearc=0
		if calendar=="365d":
			array =np.arange(int(totaltime)+dtime,0, -dtime)
			
			for i in array:
				a=str(time_calcminutes(fecha,float(i)*(-1))).split(" ")[0]
				if int(a.split("-")[1])==2 and int(a.split("-")[2])==29: 
					leapyearc+=1
			
		nhour=int(totaltime)+dtime +dtime*leapyearc
		
		
		array =np.arange(nhour,0, -dtime)
		for i in array:
			a=str(time_calcminutes(fecha,float(i)*(-1)))
			var1=a.split(" ")
			var11=var1[0].split("-")
			var12=var1[1].split(":")
			fecha_dia=str(var11[0]+var11[1]+var11[2]+var12[0]+var12[1]+var12[2])
			name=path+"partposit_"+fecha_dia
			if key_gz:
				name=path+"partposit_"+fecha_dia+".gz"
			else:
				name=path+"partposit_"+fecha_dia
			
			
			if calendar=="365d" and int(var11[0])%4==0 and int(var11[1])==2 and  int(var11[2])==29:
				msg='nothong to do'
			else:
			
				list_fecha=np.append(list_fecha, name)
			
				listdates=np.append(listdates, int(fecha_dia))
		fecha_=fecha.split(" ")
		var11=fecha_[0].split("-")
		var12=fecha_[1].split(":")
		fecha_dia=str(var11[0]+var11[1]+var11[2]+var12[0]+var12[1]+var12[2])
		if key_gz:
			name=path+"partposit_"+fecha_dia+".gz"
		else:
			name=path+"partposit_"+fecha_dia
			
					
		list_fecha=np.append(list_fecha, name)
		listdates=np.append(listdates, int(fecha_dia))
		

		
	if mode == "forward":
		
		leapyearc=0
		if calendar=="365d":
			array =np.arange(0,int(totaltime)+dtime, dtime)
			
			for i in array:
				a=str(time_calcminutes(fecha,float(i)*(1))).split(" ")[0]
				if int(a.split("-")[1])==2 and int(a.split("-")[2])==29: 
					leapyearc+=1
		
		
	
		nhour=int(totaltime)+dtime +dtime*leapyearc
		array =np.arange(0,nhour+dtime, dtime)
	
		for i in array:
			a=str(time_calcminutes(fecha,float(i)))
			var1=a.split(" ")
			var11=var1[0].split("-")
			var12=var1[1].split(":")
			fecha_dia=str(var11[0]+var11[1]+var11[2]+var12[0]+var12[1]+var12[2])
			name=path+"partposit_"+fecha_dia
			if key_gz:
				name=path+"partposit_"+fecha_dia+".gz"
			else:
				name=path+"partposit_"+fecha_dia
				
				
			if calendar=="365d" and int(var11[0])%4==0 and int(var11[1])==2 and  int(var11[2])==29:
				msg='nothong to do'
			else:	
				
				
				list_fecha=np.append(list_fecha, name)
				listdates=np.append(listdates, int(fecha_dia))
    
	return list_fecha, listdates

def read_proccesor(verbose, partpositfiles,submatrix, rank, lon_left_lower_corner,lat_left_lower_corner,
               lon_right_upper_corner,lat_right_upper_corner, model, key_gz, type_file):

	a1=np.arange(len(partpositfiles))
	dx,dy =submatrix.shape
	tensor_local=np.ones((len(partpositfiles),dx,dy))*(-999.9)
	for i in a1:
		if verbose:
			print ("Reading | " + model+" -> ",  partpositfiles[i])
		if key_gz:
			desc_gz(partpositfiles[i])
			part_post_i=read_binaryFile_fortran(partpositfiles[i][:-3], type_file,lon_left_lower_corner,lat_left_lower_corner,
				lon_right_upper_corner,lat_right_upper_corner)
			cmd_rm= "rm -rf "+partpositfiles[i][:-3]
			os.system(cmd_rm)
		else:
			part_post_i=read_binaryFile_fortran(partpositfiles[i], type_file,lon_left_lower_corner,lat_left_lower_corner,
				lon_right_upper_corner,lat_right_upper_corner)
		matrix_i=search_row_fortran(submatrix[:,0],part_post_i)
			
		tensor_local[i,:,:]=matrix_i
	return tensor_local


def get_vars_from_partposit(verbose, partpositfiles ,file_mask, maskname,maskvar_lon, maskvar_lat,lat_f, lon_f,rank,size, comm, type_file,
                 lon_left_lower_corner,lat_left_lower_corner, lon_right_upper_corner,lat_right_upper_corner, model, value_mask, key_gz, var_heat_track):
	
	name_file=partpositfiles[-1]
	if rank==0:
		if verbose:
			print ("Reading | " + model+" -> ",  name_file)
		name_txt_part=name_file.split("/")
	if key_gz:
		desc_gz(name_file)
		part_post=read_binaryFile_fortran(name_file[:-3], type_file, lon_left_lower_corner,lat_left_lower_corner,
				lon_right_upper_corner,lat_right_upper_corner)
		cmd_rm= "rm -rf "+name_file[:-3]
		os.system(cmd_rm)
	else:
		part_post=read_binaryFile_fortran(name_file, type_file, lon_left_lower_corner,lat_left_lower_corner,
				lon_right_upper_corner,lat_right_upper_corner)

	


	lat_masked, lon_masked, mascara=load_mask_grid_NR(file_mask, maskname,maskvar_lon, maskvar_lat)

	submatrix=determine_id_binary_grid_NR_fortran(part_post, lat_masked.flatten(), lon_masked.flatten(), mascara.flatten(), value_mask)

	submatrix=submatrix[np.argsort(submatrix[:, 0])]


	if rank==0:
		if verbose:
			print ("Reading | " + model+" -> ",  partpositfiles[-2])
		

	if key_gz:
		desc_gz(partpositfiles[-2])
		part_post_i=read_binaryFile_fortran(partpositfiles[-2][:-3], type_file, lon_left_lower_corner,lat_left_lower_corner,lon_right_upper_corner,lat_right_upper_corner)
		cmd_rm= "rm -rf "+partpositfiles[-2][:-3]
		os.system(cmd_rm)
	else:
		part_post_i=read_binaryFile_fortran(partpositfiles[-2], type_file, lon_left_lower_corner,lat_left_lower_corner,lon_right_upper_corner,lat_right_upper_corner)
	matrix_i=search_row_fortran(submatrix[:,0],part_post_i)

	


	dimX, dimY=matrix_i.shape

	tensor_org=np.ones((len(partpositfiles),dimX ,dimY ))*(-999.9)
	tensor_org[-1,:,:]=submatrix
	tensor_org[-2,:,:]=matrix_i


	n = len(partpositfiles)-2
	count = n // size
	remainder = n % size


	if rank < remainder:
		start = rank * (count + 1)
		stop = start + count + 1
	else:
		start = rank * count + remainder
		stop = start + count

	local_list=partpositfiles[start:stop]
	local_results = np.empty((len(local_list), dimX, dimY))
	local_results= read_proccesor(verbose,local_list, submatrix, rank,lon_left_lower_corner,lat_left_lower_corner,
				lon_right_upper_corner,lat_right_upper_corner, model, key_gz, type_file)
	
	if rank > 0:
		comm.Send(local_results, dest=0, tag=14)
	else:
		i_start=[]
		i_stop=[]
		for i in range(size):
			if i < remainder:
				i_start = np.append(i_start,i * (count + 1))
				i_stop = np.append(i_stop,i_start + count + 1)
			else:
				ii_start=i * count + remainder
				ii_stop=ii_start + count
				i_start = np.append(i_start,ii_start)
				i_stop = np.append(i_stop, ii_stop)
		final_results = np.copy(local_results)
		tensor_org[int(i_start[0]):int(i_stop[0]),:,:]= final_results
		for i in range(1, size):
			if i < remainder:
				rank_size = count + 1
			else:
				rank_size = count
			tmp = np.empty((rank_size, final_results.shape[1],final_results.shape[2]), dtype=np.float64)
			comm.Recv(tmp, source=i, tag=14)
			tensor_org[int(i_start[i]):int(i_stop[i]),:,:]=tmp
	comm.Bcast(tensor_org, root=0)

	
	tensor_final=np.empty((tensor_org.shape[0],tensor_org.shape[1],tensor_org.shape[2]+2))
	tensor_final[:]=-999.9
	tensor_final[:,:,:-2]=tensor_org
	
	if var_heat_track=="dse":

		tensor_final[:,:,-2] = compute_dry_static_energy(tensor_org[:,:,2], tensor_org[:,:,9], tensor_org[:,:,4])
		tensor_final[:,:,-2]=tensor_final[:,:,-2]/1000
		
	elif var_heat_track=="potTemp":
		
		tensor_final[:,:,-2] = compute_theta(tensor_org[:,:,6], tensor_org[:,:,3], tensor_org[:,:,9])



	tensor_final[:,:,-1] = compute_rh(tensor_org[:,:,6], tensor_org[:,:,3], tensor_org[:,:,9])

	tensor_final[:,:,-2][tensor_final[:,:,-2]<0]=-999.9
	tensor_final[:,:,-1][tensor_final[:,:,-1]<0]=-999.9
	
	#print(tensor_final[-1,0,:])
	
	return tensor_final
	

def compute_rh(rho_kgm3, q_kgkg, T_K):
    
    p_Pa=calc_pres(rho_kgm3, q_kgkg, T_K)
    
    
    e = q_kgkg * p_Pa / (0.622 + 0.378 * q_kgkg)

   
    es = 611.2 * np.exp(17.67 * (T_K - lc.TREF) / (T_K - lc.TREF + 243.5))

    return 1e2 * e / es


def calc_pres(rho_kgm3, q_kgkg, T_K):
	r_kgkg = -q_kgkg / (q_kgkg - 1)
	Tv_K = T_K * (1 + r_kgkg / lc.EPSILON) / (1 + r_kgkg)
	return rho_kgm3 * lc.RSPECIFIC * Tv_K



def calc_pottemp(p_Pa, q_kgkg, T_K):
    r_kgkg = -q_kgkg / (q_kgkg - 1) 
    r_gkg = r_kgkg * 1e3

    p_hPa = p_Pa / 1e2
    return T_K * (1000 / p_hPa) ** (0.2854 * (1 - 0.00028 * r_gkg))


def ajust_units(array, area, density, dtime, varid, var_heat_track):
	if varid==0:
		if var_heat_track=="dse":
			array = 1000*array*density/(area*dtime*60)
		elif var_heat_track=="potTemp":
			array = lc.cp*array*density/(area*dtime*60)
	elif varid==1:
		array = array*density/(area)
	return array


def cal_track_diff(var_vals, case):
	
	
	var_vals[var_vals<-500]=np.nan
	# difference
	if case in ["diff"]:
		dvalue = var_vals[1:] - var_vals[:-1]
	# mean
	if case in ["mean"]:
		dvalue = (var_vals[1:] + var_vals[:-1]) / 2
	# max
	if case in ["max"]:
		dvalue = np.amax(var_vals[1:], var_vals[:-1])
	return dvalue




############ SPECIFIC FUNCTIONS FOR HEAT TRACKING  ################################ 

def compute_dry_static_energy(lat, T_K, heigh):
	dse=lc.cp*T_K + lc.g0*(1 + 0.0053024 * (np.sin(lat*np.pi/180 ))**2 - 0.0000058*(np.sin( 2*lat*np.pi/180))**2)*heigh

	return dse



def compute_theta(rho_kgm3, q_kgkg, T_K):
	
	parpres =  calc_pres(rho_kgm3, q_kgkg, T_K)
	
	theta = calc_pottemp(parpres,  q_kgkg, T_K)
	
	return theta

def compute_var_integarated_day_heat(array, t, area, density, dtime, var_heat_track, lon,lat,numPdY,numPdX, varid):
	dimX, dimY=lat.shape
	array_day=np.empty((len(t)-1,dimX-1, dimY-1))
	ndb=np.arange(len(t)-1,0,-1)
	for ii in range(len(t)-1):
		heatd=np.array(compute_grid_integrated_heat(array[t[ii]:t[ii+1],:,:],lon,lat,numPdY,numPdX,len(array[t[ii]:t[ii+1],0,0]),len(array[0,:,0])),dtype=np.float64)

		array_day[int (ndb[ii]-1), :,:] = ajust_units(heatd, area, density, dtime, varid, var_heat_track)
		
	return array_day



def is_withinpbl(parts_high, pblhigh, lower_limit, upper_limit, method):
	check_pblupper=np.empty(len(parts_high))
	check_pbllower=np.empty(len(parts_high))

	pblhighs=[]

	for i in range(0, pblhigh.shape[1]):
		pbl_vars=pblhigh[:,i]
	
		if method == "meanval":
			pblval = np.mean(pbl_vars[-4:])	
			pblhighs=np.append(pblhighs, pblval)
		elif method == "maxval":
			pblval = np.max(pbl_vars[-4:])
			
			pblhighs=np.append(pblhighs, pblval)
		elif method == "actualval":
			pblval = pbl_vars[-1]
			pblhighs=np.append(pblhighs, pblval)
		else:
			print_error_message(" pbl_method or dqpbl_method is not valid. This parameter must be equal to  ['maxval'/'meanval'/'actualval']")
	
	
	pblhighs[pblhighs<=upper_limit]=upper_limit
	
	check_pblupper[parts_high<=pblhighs]=True
	check_pblupper[parts_high>pblhighs]=False

	check_pbllower[parts_high>=lower_limit]=True
	check_pbllower[parts_high<lower_limit]=False
	
	
	check_pbl=np.logical_and(check_pbllower,check_pblupper)

	return check_pbl


def maxval(x, n=2):
	
	if len(x) == 2:
		return max(x[0], x[1])
	else:
		return np.array([np.max(x[i : i + n]) for i in range(len(x) - (n - 1))])


def meanval(x, n=2):
	if len(x) == 2:
		return np.mean([x[0], x[1]])
	else:
		return np.array([np.mean(x[i : i + n]) for i in range(len(x) - (n - 1))])


def is_in_pbl(pblcheck, par_vals, pbl_vars, method, lendvar):
	
	if not pblcheck in (0,1,2):
		print_error_message("pblcheck is not valid. This parameter must be equal to  0: no PBL check is applied, use everything, 1:  at least one location is within the PBL\n       2: both locations are within the max PBL ")
	
	if pblcheck == 0:
		check_pbl=np.ones((lendvar), dtype=bool)
		check_pbl[:]=True
		return check_pbl
		
	else:
		if method == "meanval":
			before = par_vals[:-1] <= meanval(pbl_vars, n=2)
			after = par_vals[1:] <= meanval(pbl_vars, n=2)
		elif method == "maxval":
			before = par_vals[:-1] <= maxval(pbl_vars, n=2)
			after = par_vals[1:] <= maxval(pbl_vars, n=2)
		
		elif method == "actualval":
			before = par_vals[:-1] <= pbl_vars[:-1]
			after = par_vals[1:] <= pbl_vars[1:]
		else:
			print_error_message(" pbl_method or dqpbl_method is not valid. This parameter must be equal to  ['maxval'/'meanval'/'actualval']")
		if pblcheck == 2:
			return np.logical_and(before, after)
		elif pblcheck == 1:
			return np.logical_or(before, after)




	
def processing_heat_track(tensor_org, pblcheck, filter_pbl_parcels, pbl_method, heat_custom_limits, trk_rh_check, rh_threshold, dqcheck, dqthreshold, dvar_threshold, var_heat_track, lag_times, area, density, dtime,heat_linear_adjustment, lon,lat, numPdY,numPdX, varid):
	
	
	if filter_pbl_parcels:
	
		check_pbl=is_withinpbl(tensor_org[-1,:,4], tensor_org[1:,:,7], heat_custom_limits[0], heat_custom_limits[1], "maxval")
		
		if len(check_pbl[check_pbl==True])==0:
			tensor_heat=np.empty((tensor_org.shape[0]-1, 1, tensor_org.shape[2]))
			tensor_heat[:,:,:]=0
		
		else:
		
			tensor_heat=np.empty((tensor_org.shape[0]-1, len(check_pbl[check_pbl==True]), tensor_org.shape[2]))

			for i in range(0, tensor_heat.shape[0]):
				tensor_heat[i,:]=tensor_org[i+1,:,:][check_pbl==True]
	else:
		tensor_heat=np.copy(tensor_org[1:,:,:])
	

	matrix_heat_files=[ "lon", "lat", "dq", "ds", "rh", "FH"]
	dmatrix=np.empty((tensor_heat.shape[0]-1, tensor_heat.shape[1], len(matrix_heat_files)))
	uptakes_parts=0
	for i in range(0, tensor_heat.shape[1]):
		
		lons=tensor_heat[:,i, 1]
		dlon=cal_track_diff(lons, "mean")
			
		lats=tensor_heat[:,i, 2]
		dlat=cal_track_diff(lats, "mean")
		
		qs=tensor_heat[:,i, 3]
		dq=cal_track_diff(qs, "diff")
		
		var=tensor_heat[:,i, 11]
		dvar=cal_track_diff(var, "diff")

		rh=tensor_heat[:,i, 12]
		drh=cal_track_diff(rh, "diff")

		dmatrix[:,i,0]=dlon
		dmatrix[:,i,1]=dlat
		dmatrix[:,i,2]=dq
		dmatrix[:,i,3]=dvar
		dmatrix[:,i,4]=drh
		dmatrix[:,i,5]=1
		

		trk_drh=np.ones((len(dlon)), dtype=bool)
		
		#trk_check_pbl=is_pbl_check(pblcheck, tensor_heat[:,i,4], tensor_heat[:,i,4], len(dlon))
		trk_check_pbl=is_in_pbl(pblcheck, tensor_heat[:,i,4], tensor_heat[:,i,7], pbl_method, len(dlon))
		
			
		if trk_rh_check:
			
			trk_drh[np.abs(drh)<=rh_threshold]=True
			trk_drh[np.abs(drh)>rh_threshold]=False
		else:
			trk_drh[:]=True
		
		change_dq=np.ones((len(dlon)), dtype=bool)
		if dqcheck:
			
			changedq=np.abs(dq)/qs[:-1]		
			change_dq[changedq<=dqthreshold]=True
			change_dq[changedq>dqthreshold]=False
					
		else:
			change_dq[:]=True
		
		
		
		
		check_dvar=np.ones((len(dlon)), dtype=bool)
		
		
		check_dvar[dvar>dvar_threshold]=True
		check_dvar[dvar<=dvar_threshold]=False
	
		valid=trk_check_pbl&trk_drh&change_dq&check_dvar

		dmatrix[:,i,5][valid==True]=1
		dmatrix[:,i,5][valid==False]=0
	
	
		#print(dmatrix.max())
	
		dmatrix[np.isnan(dmatrix)]=-999
	
		if len(valid[valid==True])>=1:
			uptakes_parts = uptakes_parts + 1
	
	
		if heat_linear_adjustment:
		
			dmatrix[:,i,3] = compute_linear_discounted(dmatrix[:,i,3], valid)
	
	matrix_heat=np.copy(dmatrix)
	matrix_heat[:,:,3][matrix_heat[:,:,3]==-999]=0
	

	array_day=compute_var_integarated_day_heat(matrix_heat, lag_times, area, density, dtime, var_heat_track, lon,lat,numPdY,numPdX, varid)
		
	counter_part=matrix_heat.shape[1]
	no_uptakes_parts=counter_part-uptakes_parts
	#return array_day/(lag_times[1]-lag_times[0]), counter_part
	return array_day, matrix_heat, counter_part, no_uptakes_parts
		

############ SPECIFIC FUNCTIONS FOR MOISTURE TRACKING  ################################

def is_precipitating_parcel(parts_dq, dqdt_threshold, parts_rh, minrh):
	check_precipdq=np.empty(len(parts_dq))
	check_preciprh=np.empty(len(parts_dq))
	
	
	check_precipdq[parts_dq<dqdt_threshold]=True
	check_precipdq[parts_dq>=dqdt_threshold]=False

	check_preciprh[parts_rh<minrh]=False
	check_preciprh[parts_rh>=minrh]=True

	check_precip=np.logical_and(check_precipdq,check_preciprh)


	return check_precip


def	calc_dvar_moisture(matrixresult, tensorvar, slices, npp,dqdt_threshold):
	
	for i in slices[::-1]:
					
		matrix=get_dqdt(tensorvar[i,:,:], tensorvar[i+1,:,:], npp, len(tensorvar[0,:,0]), 13)
		
		matrixresult[i,:,2]=matrix[:,2]
		matrixresult[i,:,1]=matrix[:,1]
		matrixresult[i,:,0]=matrix[:,0]
		matrixresult[i,:,3]=matrix[:,3]
		matrixresult[i,:,4]=matrix[:,4]
		matrixresult[i,:,5]=matrix[:,5]
	return matrixresult

def compute_var_integarated_day_moist(array, t, area, density, dtime, lon,lat,numPdY,numPdX, varid):
	dimX, dimY=lat.shape
	array_day=np.empty((len(t)-1,dimX-1, dimY-1))

	array_day_cr=np.empty((len(t)-1,dimX-1, dimY-1))

	
	ndb=np.arange(len(t)-1,0,-1)
	for ii in range(len(t)-1):
		moistd=np.array(compute_grid_integrated_moist(array[t[ii]:t[ii+1],:,:],lon,lat,numPdY,numPdX,len(array[t[ii]:t[ii+1],0,0]),len(array[0,:,0])),dtype=np.float64)

		array_day[int (ndb[ii]-1), :,:] = ajust_units(moistd, area, density, dtime, varid, "None")
		array_day_cr[int (ndb[ii]-1), :,:] = moistd
	return array_day, array_day_cr


def processing_moisture_track(tensor_org, filter_dqdt_parcels, dqdt_threshold, 	filter_pbl_dq_parcels, moist_custom_limits_highs, dqpblcheck, dqpbl_method, trkdq_rh_check, dqrh_threshold, mindq_gain, lag_times, area, density, dtime, moisture_linear_adjustment, precip_minrh, lon,lat, numPdY,numPdX, varid):
	#filtering parcels using t=0  and t-6#
	#tensor_moist_=np.copy(tensor_org)
	if filter_dqdt_parcels:

		tmp_matrix=tensor_org[-1,:,3] - tensor_org[-2,:,3]

		check_precip=is_precipitating_parcel(tmp_matrix, dqdt_threshold,  (tensor_org[-1,:,12] + tensor_org[-2,:,12])/2, precip_minrh)
				
		if len(check_precip[check_precip==True])==0:
			tensor_moist_=np.empty((tensor_org.shape[0], 1, tensor_org.shape[2]))
			tensor_moist_[:,:,:]=0
			counter_precip_partdq=0
			
		else:
			tensor_moist_=np.empty((tensor_org.shape[0], len(check_precip[check_precip==True]), tensor_org.shape[2]))

			for i in range(0, tensor_moist_.shape[0]):
				tensor_moist_[i,:]=tensor_org[i,:,:][check_precip==True]


			counter_precip_partdq=tensor_moist_.shape[1]
	else:
		tensor_moist_=np.copy(tensor_org[:,:,:])
		counter_precip_partdq=None
	
	
	if filter_pbl_dq_parcels:
		check_pbl=is_withinpbl(tensor_moist_[-1,:,4], tensor_moist_[1:,:,7], moist_custom_limits_highs[0], moist_custom_limits_highs[1], "maxval")
		
		if len(check_pbl[check_pbl==True])==0:
			tensor_moist=np.empty((tensor_moist_.shape[0], 1, tensor_moist_.shape[2]))
			tensor_moist[:,:,:]=0
			counter_precip_part_pbl=0
		else:
			tensor_moist=np.empty((tensor_moist_.shape[0], len(check_pbl[check_pbl==True]), tensor_moist_.shape[2]))

			for i in range(0, tensor_moist.shape[0]):
				tensor_moist[i,:]=tensor_moist_[i,:,:][check_pbl==True]
		
			counter_precip_part_pbl=tensor_moist.shape[1]
	else:
		tensor_moist=np.copy(tensor_moist_[:,:,:])
		counter_precip_part_pbl=None
	



	precipvals=tensor_moist[-1,:,3]-tensor_moist[-2,:,3]
	qt0=tensor_moist[-1,:,3]
	Qt0_=np.sum(tensor_moist[-2,:,3])
	
	tensor_moist=tensor_moist[:-1,:,:]
	Qt0_=np.sum(tensor_moist[-1,:,3])
	
	matrix_moist_files=[ "lon", "lat", "dq", "ds",'rh', "FH"]
	dmatrix=np.ones((len(tensor_moist[:,0,0])-1, len(tensor_moist[0,:,0]),len(matrix_moist_files)))*(-999.9)
	
	uptakes_parts_precip=0
	sum_prec=0
	partt=0
	noprecip=0
	for i in range(0, tensor_moist.shape[1]):
		
		
		
			
		lons=tensor_moist[:,i, 1]
		dlon=cal_track_diff(lons, "mean")
		
		lats=tensor_moist[:,i, 2]
		dlat=cal_track_diff(lats, "mean")
		
		qs=tensor_moist[:,i, 3]
		dq=cal_track_diff(qs, "diff")
		
		sum_prec= sum_prec + np.abs(precipvals[i])
		
		
		var=tensor_moist[:,i, 11]
		dvar=cal_track_diff(var, "diff")

		rh=tensor_moist[:,i, 12]
		drh=cal_track_diff(rh, "diff")

		dmatrix[:,i,0]=dlon
		dmatrix[:,i,1]=dlat
		dmatrix[:,i,2]=dq
		dmatrix[:,i,3]=dvar
		dmatrix[:,i,4]=drh
		dmatrix[:,i,5]=0
		
		#trk_check_pbl=is_pbl_check(dqpblcheck, tensor_moist[:,i,4], tensor_moist[:,i, 7], len(dlon))
		
				
		trk_check_pbl=is_in_pbl(dqpblcheck, tensor_moist[:,i,4], tensor_moist[:,i,7], dqpbl_method, len(dlon))
		
		check_dvar=np.ones((len(dlon)), dtype=bool)
		check_dvar[dq>=mindq_gain]=True
		check_dvar[dq<mindq_gain]=False
		
		trk_drh=np.ones((len(dlon)), dtype=bool)
		if trkdq_rh_check:
			
			trk_drh[np.abs(drh)<=dqrh_threshold]=True
			trk_drh[np.abs(drh)>dqrh_threshold]=False
		
		valid=trk_check_pbl&trk_drh&check_dvar
		
		dmatrix[:,i,5][valid==True]=1
		dmatrix[:,i,5][valid==False]=0
		
		if len(valid[valid==True])>=1:
			uptakes_parts_precip = uptakes_parts_precip + 1
	
		
		dmatrix[np.isnan(dmatrix)]=-999
		
		if moisture_linear_adjustment:
		
			dmatrix[:,i,2] = compute_linear_discounted(dmatrix[:,i,2], valid)
				
			#adj_fac = dmatrix[:,i,2] / qs[-1]
							
			dmatrix[:,i,2]=dmatrix[:,i,2]
			
			partt+=np.sum(dmatrix[:,i,2])
	
		
		noprecip+=qt0[i]
	
	
	matrix_moist=np.copy(dmatrix)
	
	matrix_moist[:,:,2][matrix_moist[:,:,2]==-999]=0
	

	array_day,moistd=compute_var_integarated_day_moist(matrix_moist, lag_times, area, density, dtime, lon,lat,numPdY,numPdX, varid)
	

	counter_precip_part=matrix_moist.shape[1]
		
	no_evap_uptakes=matrix_moist.shape[1]-uptakes_parts_precip
	
	if moisture_linear_adjustment:
		CR = (moistd/Qt0_)*100
				
	else:
		CR=None

	return array_day, matrix_moist, counter_precip_part, no_evap_uptakes, partt/sum_prec, CR



def compute_linear_discounted(var, valid):
	
	var[var==-999]=0
		
	result_var=np.empty_like(var)
	result_var[:]=0
	
	for i in range(0, len(var)):
		if var[i]>0:
			result_var[i]=var[i]
			
		else:
			suma=np.sum(result_var[:i])
			
			for j in range(0,i):
				if suma>0:
					aux=result_var[j]-((result_var[j]/suma)*abs(var[i]))
				else:
					aux=0
				
				if aux<=0:
					result_var[j]=0
				else:
					result_var[j]=aux
		

	return result_var
	
