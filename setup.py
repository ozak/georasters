# coding: utf-8
import os
import sys
from setuptools import setup
from setuptools.command.test import test as TestCommand

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def readme():
    with open('README.rst') as f:
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
      version='0.5.4',
      description='Tools for working with Geographical Information System Rasters',
      url='http://github.com/ozak/georasters',
      keywords="gis geospatial geographic raster vector zonal statistics spatial analysis",
      author='Ömer Özak',
      author_email='omer@omerozak.com',
      license='GPLv3',
      #package_dir={'': 'src'},
      packages=['georasters'],
      long_description=read('README.rst'),
      install_requires=[
                'numpy',
                'GDAL',
                'docopt',
                'pandas', 
                'pyproj',
                'scikit-image',
                'matplotlib',
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
          'Programming Language :: Python :: 2.6',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
      ],
      zip_safe=False,
      )
