#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy
import glob

ext_modules=[ Extension("flipflop_model",
              ["flipflop/flipflop_model.pyx"],
              extra_compile_args = ["-ffast-math"])]

scripts = glob.glob("scripts/*.py")

setup(
    name="flipflop",
    url="https://github.com/CalumGabbutt/flipflop",
    version=1.0,
    author="Calum Gabbutt",
    author_email="c.gabbutt@qmul.ac.uk",
    packages=["flipflop"],
    license="MIT",
    scripts=scripts,
    ext_modules = cythonize(ext_modules, annotate=False),
    include_dirs=[numpy.get_include()],
    description=("A Bayesian pipeline to infer stem cell"
                 "dynamics from methylation array data."),
    install_requires=["numpy", "scipy", "matplotlib", "pandas", 
                    "dynesty", "joblib", "seaborn", "arviz", "cython"],
    package_data={
        "flipflop": ["files/*"],
    }
)