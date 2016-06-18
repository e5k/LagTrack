#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer

# Time of dataset
date_start  = 'var_date_start'
date_end 	= 'var_date_end'
# Area
north       = var_north
south       = var_south
west        = var_west
east        = var_east
# Output folder, i.e. replace by your project name
out_path    = 'var_out'


##################################################
#os.mkdir(out_path)
#os.mkdir(out_path+"nc_output_files")
#os.mkdir(out_path+"txt_output_files")
server = ECMWFDataServer()

server.retrieve({
	'dataset'   : "interim",
	'date'      : "%s/to/%s"%(date_start,date_end),
	'time'      : "00/06/12/18",
	'step'      : "0",
	'stream'    : "oper",
	'levtype'   : "pl",
	'levelist'  : "all",
	'type'      : "an",
	'class'     : "ei",
	'grid'      : "0.25/0.25",
	'param'     : "129/130/131/132/135/156/157",
	'area'      : "%d/%d/%d/%d"%(north, west, south, east),
	'format'	: 'netcdf',
	'target'    : "%s"%(out_path)
}) 

