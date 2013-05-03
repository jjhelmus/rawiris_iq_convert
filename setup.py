from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from numpy import get_include

ext_modules = [Extension(
    "decode_highsnr",
    ["decode_highsnr.pyx"],
    include_dirs=[get_include()])]

setup(
    name='foobar',
    cmdclass={'build_ext': build_ext},
    ext_modules=ext_modules
)
