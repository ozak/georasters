# Test
python setup.py test

# To upload to PyPy
#python setup.py register
#python setup.py sdist bdist_wheel bdist_wininst upload
python setup.py sdist
python setup.py sdist bdist_wheel --universal
twine upload dist/*

# To get sha256
wget https://files.pythonhosted.org/packages/c5/53/898fc08433a3fc5c5cac490edfeb61e621f9970c2d8647711037bae1190e/georasters-0.5.13.tar.gz
openssl dgst -sha256 ./georasters-0.5.13.tar.gz