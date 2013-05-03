rawiris_iq_convert
==================

Convert a RAW IRIS IQ file to netcdf

Requirements
============
* Python 2.7
* Numpy 1.6.1
* Cython ?

Building
========
Run python setup.py build_ext -i to compile the decode_highsnr.pyx module.

Running
=======
From within this directory run
convert.py iris_file netcdf_file
