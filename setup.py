from distutils.core import setup
from Cython.Build import cythonize
import numpy as np

setup(
    name = "AmpliDiff",
    version = "0.0.1",
    ext_modules  = cythonize('amplicon_generation.pyx'), 
    include_dirs=np.get_include()
    )