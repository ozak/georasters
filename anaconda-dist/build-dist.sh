#!/usr/bin/bash
source activate GeoPython2env

VERSION="0.5.5"

# Build georasters
conda skeleton pypi georasters
conda build -c conda-forge -c anaconda -c ozak --python 2.7 --skip-existing georasters
conda build -c conda-forge -c anaconda -c ozak --python 3.5 --skip-existing georasters
conda build purge

## if it fails, you need to run these
#conda skeleton pypi geopandas
#conda build -c conda-forge --python 2.7 geopandas
#conda build -c conda-forge --python 3.5 geopandas
#conda build purge
#
#conda skeleton pypi dask
#conda build -c conda-forge --python 2.7 dask
#conda build -c conda-forge --python 3.5 dask
#conda build purge

# Convert to other OS
conda convert --platform all ~/Anaconda3/conda-bld/osx-64/georasters-$VERSION-py27_0.tar.bz2 -o ~/GitHub/georasters/anaconda-dist/
conda convert --platform all ~/Anaconda3/conda-bld/osx-64/georasters-$VERSION-py35_0.tar.bz2 -o ~/GitHub/georasters/anaconda-dist/

# Upload 
#2.7
anaconda upload ~/GitHub/georasters/anaconda-dist/linux-32/georasters-$VERSION-py27_0.tar.bz2
anaconda upload ~/GitHub/georasters/anaconda-dist/linux-64/georasters-$VERSION-py27_0.tar.bz2
anaconda upload ~/GitHub/georasters/anaconda-dist/win-32/georasters-$VERSION-py27_0.tar.bz2
anaconda upload ~/GitHub/georasters/anaconda-dist/win-64/georasters-$VERSION-py27_0.tar.bz2
# 3.5
anaconda upload ~/GitHub/georasters/anaconda-dist/linux-32/georasters-$VERSION-py35_0.tar.bz2
anaconda upload ~/GitHub/georasters/anaconda-dist/linux-64/georasters-$VERSION-py35_0.tar.bz2
anaconda upload ~/GitHub/georasters/anaconda-dist/win-32/georasters-$VERSION-py35_0.tar.bz2
anaconda upload ~/GitHub/georasters/anaconda-dist/win-64/georasters-$VERSION-py35_0.tar.bz2
