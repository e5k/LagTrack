#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
import calendar
import os
## Enter the required data below

# Time of dataset
year_start  = 1980
year_end    = 1980
month_start = 5
month_end   = 5
# Area
north       = 48
south       = 41
west        = 237
east        = 260
# Output folder, i.e. replace by your project name
out_path    = 'input/wind/MSH_ECMWF/MSH_ECMWF.nc'


## Time of dataset
#year_start  = 2010
#year_end    = 2013
#month_start = 1
#month_end   = 12
## Area
#north       = -37
#south       = -40
#west        = -40
#east        = -38
## Output folder, i.e. replace by your project name
#out_path    = 'OUT/'


##################################################
#os.mkdir(out_path)
#os.mkdir(out_path+"nc_output_files")
#os.mkdir(out_path+"txt_output_files")
server = ECMWFDataServer()
count  = 1
for year in range(year_start, year_end+1):
    print 'YEAR ',year
    for month in range(month_start, month_end+1):
        lastday1=calendar.monthrange(year,month)
        lastday=lastday1[1]
        bdate="%s%02d01"%(year,month)
        edate="%s%02d%s"%(year,month,lastday)

        print "######### ERA-interim  #########"
        print 'Accessing wind data from ', bdate,' to ',edate,' (YYYYMMDD)'
        print "################################"
        
        server.retrieve({
            'dataset'   : "interim",
            'date'      : "%s/to/%s"%(bdate,edate),
            'time'      : "00/06/12/18",
            'step'      : "0",
            'stream'    : "oper",
            'levtype'   : "pl",
            'levelist'  : "all",
            'type'      : "an",
            'class'     : "ei",
            'grid'      : "0.25/0.25",
            'param'     : "129/131/132/156",
            'area'      : "%d/%d/%d/%d"%(north, west, south, east),
            'format'	: 'netcdf',
            'target'    : "%s%05d_%s_%04d.nc"%(out_path+"nc_output_files/", count, calendar.month_abbr[month],year)
        }) 
        
        count = count + 1