# Test
python setup.py test

# To upload to PyPy
#python setup.py register
#python setup.py sdist bdist_wheel bdist_wininst upload
python setup.py sdist
python setup.py sdist bdist_wheel --universal
twine upload dist/*

# To get sha256
wget https://pypi.python.org/packages/28/fc/561ea81499587e592fc7311d457d7228eb8ed84af372d040bb68cd948a73/georasters-0.5.9b.tar.gz#md5=3d4ba6cb4b3487f9a0074ac176d7df7c
openssl dgst -sha256 ./georasters-0.5.9b.tar.gz