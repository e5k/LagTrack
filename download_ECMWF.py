#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer

# Time of dataset
date_start  = '19960601'
date_end 	= '19960630'
# Area
north       = -37.5
south       = -41
west        = 173.8
east        = 177.5
# Output folder, i.e. replace by your project name
out_path    = 'input/wind/ruapehu_R2/ruapehu_R2.nc'


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

