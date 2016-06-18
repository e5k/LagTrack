#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer

# Time of dataset
date_start  = '20060201'
date_end 	= '20060228'
# Area
north       = 47
south       = 44
west        = 4
east        = 6
# Output folder, i.e. replace by your project name
out_path    = 'input/wind/gva/gva.nc'


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

