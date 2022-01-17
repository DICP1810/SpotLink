import numpy as np
from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules=cythonize("MSOperatorMass.pyx"),
    include_dirs=[np.get_include()]
)