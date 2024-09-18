#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
import fnmatch
from setuptools import setup, find_packages


def find_files(directory):
    matches = []

    for root, dirnames, filenames in os.walk(directory):
        for filename in fnmatch.filter(filenames, '*'):
            matches.append(os.path.join(root, filename))

    return matches


setup(name="bottchscore",
      version='0.2.0',
      description="Calculate the information complexity of a molecular structure",
      author="Forli Lab",
      author_email="forli@scripps.edu",
      url="https://github.com/forlilab/bottchscore",
      packages=find_packages(exclude=['docs']),
      include_package_data=True,
      zip_safe=False,
      license="LGPL-2.1",
      keywords=["cheminformatics", "drug design"],
      classifiers=["Programming Language :: Python",
                   "Operating System :: Unix",
                   "Operating System :: MacOS",
                   "Topic :: Scientific/Engineering"],
      entry_points={
          'console_scripts': [
              'bottchscore_mols=bottchscore.cli.bottchscore_mols:main',
          ]
      }
      
)