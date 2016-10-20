# coding: utf-8
import os
import sys
from setuptools import setup
from setuptools.command.test import test as TestCommand

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='georasters',
      version='0.5',
      description='Tools for working with Geographical Information System Rasters',
      url='http://github.com/ozak/georasters',
      author='Ömer Özak',
      author_email='omer@omerozak.com',
      license='GPLv3',
      #package_dir={'': 'src'},
      packages=['georasters'],
      install_requires=[
          'numpy',
          'GDAL',
          'docopt',
          'pandas', 
          'pyproj',
          'scikit-image',
          'matplotlib',
      ],
      classifiers=[
          "Development Status :: 1 - Planning",
          "Topic :: Utilities",
          "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'Operating System :: OS Independent',
          'Programming Language :: Python',
          'Topic :: Scientific/Engineering :: GIS',
      ],
      zip_safe=False,
      long_description=readme(),
      )
