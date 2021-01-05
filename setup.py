#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import re
import glob
from setuptools import setup


    
setup(
    name="ticktockclock",
    url="https://github.com/CalumGabbutt/ticktockclock",
    version=1.0,
    author="Calum Gabbutt",
    author_email="c.gabbutt@qmul.ac.uk",
    packages=["ticktockclock"],
    license="MIT",
    description=("A Bayesian pipeline to infer stem cell"
                 "dynamics from methylation array data."),
    install_requires=["numpy", "scipy", "matplotlib", "pandas", 
                    "dynesty", "joblib", "seaborn", "arviz", "cython"]
)