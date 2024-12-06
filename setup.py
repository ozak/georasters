# coding: utf-8
import os
import sys
from setuptools import setup
from setuptools.command.test import test as TestCommand

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def readme():
    with open('README.md') as f:
        return f.read()

class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest
        errno = pytest.main(self.test_args)
        sys.exit(errno)

setup(name='georasters',
      version='0.5.29',
      description='Tools for working with Geographical Information System Rasters',
      url='http://github.com/ozak/georasters',
      keywords="gis geospatial geographic raster vector zonal statistics spatial analysis",
      author='Ömer Özak',
      author_email='omer@omerozak.com',
      license='GPLv3',
      #package_dir={'': 'src'},
      packages=['georasters'],
      long_description=read('README.md'),
      install_requires=[
                'numpy',
                'pandas',
                'docopt',
                'GDAL',
                'pyproj',
                'scikit-image',
                'matplotlib',
                'coverage',
                'fiona',
                'geopandas',
                'pysal',
                'affine',
                'rasterstats'
            ],
      #install_requires=read('requirements.txt').splitlines(),
      tests_require=['pytest', 'pytest-cov>=2.2.0', 'pyshp>=1.1.4',
                     'coverage', 'simplejson'],
      cmdclass={'test': PyTest},
      classifiers=[
          "Development Status :: 3 - Alpha",
          "Topic :: Utilities",
          "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'Operating System :: OS Independent',
          'Programming Language :: Python',
          'Topic :: Scientific/Engineering :: GIS',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Programming Language :: Python :: 3.9',
          'Programming Language :: Python :: 3.10',
          'Programming Language :: Python :: 3.11',
          'Programming Language :: Python :: 3.12',
      ],
      zip_safe=False,
      )
