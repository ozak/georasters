# Test
python setup.py test

# To upload to PyPy
#python setup.py register
#python setup.py sdist bdist_wheel bdist_wininst upload
python setup.py sdist
python setup.py sdist bdist_wheel --universal
twine upload dist/*

# To get sha256
wget https://files.pythonhosted.org/packages/43/9d/bb33917c788f9283e0adf5158ab6ff070f67c971943cca1d4424645a287a/georasters-0.5.16.tar.gz
openssl dgst -sha256 ./georasters-0.5.16.tar.gz
