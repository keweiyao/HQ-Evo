from setuptools import setup, Extension
from Cython.Build import cythonize
from glob import glob

# Remove the "-Wstrict-prototypes" compiler option, which isn't valid for C++.
import distutils.sysconfig
cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")

#-------------HqEvo Module------------------
filelist1 = list(set(['cython/HqEvo.pyx'] + glob('src/*.cpp'))-set(["src/main.cpp"]))
filelist3 = list(set(['cython/HqLGV.pyx'] + glob('src/*.cpp'))-set(["src/main.cpp"]))

extensions = [
        Extension('HqEvo', filelist1, language="c++", extra_compile_args=["-std=c++11", '-march=native'],libraries=["m", "gsl", "gslcblas", "boost_filesystem", "hdf5", "hdf5_cpp"]),
        Extension('HqLGV', filelist3, language="c++", extra_compile_args=["-std=c++11", '-march=native'],libraries=["m", "gsl", "gslcblas", "boost_filesystem", "hdf5", "hdf5_cpp"])
]


setup(
        ext_modules=cythonize(extensions)
)
