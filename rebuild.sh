# remove and rebuild _rsl_interface.so and _fourdd_interface.so modules
rm  decode_highsnr.so
cython decode_highsnr.pyx
python setup.py build_ext -i
