#!/bin/bash
$PYTHON setup.py install 

# Add more build steps here, if they are necessary.

# See
# http://docs.continuum.io/conda/build.html
# for a list of environment variables that are set during the build process.
# what real Python executable to use
#PYTHON=$(which python)
#
## find the root of the virtualenv, it should be the parent of the dir this script is in
#ENV=`$PYTHON -c "import os; print(os.path.abspath(os.path.join(os.path.dirname(\"$0\"), '..')))"`
#
## now run Python with the virtualenv set as Python's HOME
#export PYTHONHOME=$ENV
#exec $PYTHON "$@"
