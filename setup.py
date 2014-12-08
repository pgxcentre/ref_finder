#!/usr/bin/env python

# How to build source distribution
# python setup.py sdist --format bztar
# python setup.py sdist --format gztar
# python setup.py sdist --format zip

import os
from setuptools import setup

from ref_finder import __version__ as ref_finder_version

setup(name="ref-allele-finder",
      version=ref_finder_version,
      description="Finds the reference allele in the reference genome.",
      author="Louis-Philippe Lemieux Perreault",
      author_email="louis-philippe.lemieux.perreault@statgen.org",
      url="http://www.statgen.org",
      license="GPL",
      py_modules=["ref_finder"],
      entry_points={"console_scripts": ["ref-finder=ref_finder:main"]},
      install_requires=["pyfaidx >= 0.3.0"],
      classifiers=['Operating System :: Linux',
                   'Programming Language :: Python',
                   'Programming Language :: Python :: 3.4'])
