from setuptools import setup, Extension
from Cython.Build import cythonize
from glob import glob

extensions = [
    Extension(
        'update',
	    ['cython/update.pyx']+
		glob('src/*.cpp'),
		 language="c++",
        extra_compile_args=["-std=c++11"],
		libraries=["m", "gsl", "gslcblas"])
]

setup(
    ext_modules=cythonize(extensions)
	)
