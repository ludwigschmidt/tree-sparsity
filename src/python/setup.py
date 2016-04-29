from setuptools import Extension, find_packages, setup
import numpy as np

module = Extension(
    'treesparsity',
    sources=['python_wrapper.cc'],
    extra_compile_args=['-march=native', '-O3', '-std=c++11', '-Wall', '-Wextra', '-pedantic'],
    include_dirs=[np.get_include()])

setup(name = 'Tree sparsity',
      version = 0.2,
      packages = find_packages(),
      ext_modules = [module])
