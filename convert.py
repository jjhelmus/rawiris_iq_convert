#!/usr/bin/env python
# Convert a RAW IRIS IQ file to a netCDF file
# usage: convert.py iris_file netcdf_file

import sys

import iris_ts_file

if __name__ == '__main__':

    # check command line arguments
    if len(sys.argv) != 3:
        print "convert.py: Convert a RAW IRIS IQ file to netCDF file"
        print "usage: convert.py iris_file netcdf_file"
        sys.exit()

    # unpack command line arguments
    iris_file = sys.argv[1]
    netcdf_file = sys.argv[2]

    print "Reading in file:", iris_file
    iristsfile = iris_ts_file.IrisTSFile(iris_file)

    print "Writing NetCDF file:", netcdf_file
    iris_ts_file.write_iristsfile_to_netcdf(netcdf_file, iristsfile)
