#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import re
import glob
try:
    from setuptools import setup
    setup
except ImportError:
    from distutils.core import setup
    setup

# dir_path = os.path.dirname(os.path.realpath(__file__))

# init_string = open(os.path.join(dir_path, 'dynesty', '__init__.py')).read()
# VERS = r"^__version__ = ['\"]([^'\"]*)['\"]"
# mo = re.search(VERS, init_string, re.M)
# __version__ = mo.group(1)

# try:
#     import pypandoc
#     with open('README.md', 'r') as f:
#         txt = f.read()
#     txt = re.sub('<[^<]+>', '', txt)
#     long_description = pypandoc.convert(txt, 'rst', 'md')
# except ImportError:
#     long_description = open('README.md').read()

    
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
    # package_data={"": ["README.md", "LICENSE", "AUTHORS.md"]},
    # include_package_data=True,
    install_requires=["numpy", "scipy", "matplotlib", "pandas", 
                    "dynesty", "joblib", "seaborn", "arviz"]
)