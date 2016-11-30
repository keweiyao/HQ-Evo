from setuptools import setup, Extension
from Cython.Build import cythonize
from glob import glob

#-------------HqEvo Module------------------
filelist1 = list(set(['cython/HqEvo.pyx'] + glob('src/*.cpp'))-set(["src/main.cpp", "src/utility.cpp"]))
extensions = [
    Extension(
        'HqEvo',
		filelist1,
		 language="c++",
        extra_compile_args=["-std=c++11"],
		libraries=["m", "gsl", "gslcblas"])
]

setup(
    ext_modules=cythonize(extensions)
	)

#-------------Transform Module------------------
filelist2 = ["cython/Transform.pyx", "src/utility.cpp"]
extensions = [
    Extension(
        'Transform',
		filelist2,
		language="c++",
        extra_compile_args=["-std=c++11"],
		libraries=["m"])
]

setup(
    ext_modules=cythonize(extensions)
	)
