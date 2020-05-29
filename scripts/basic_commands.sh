# Test
python setup.py test

# To upload to PyPy
#python setup.py register
#python setup.py sdist bdist_wheel bdist_wininst upload
python setup.py sdist
python setup.py sdist bdist_wheel --universal
twine upload dist/*

# To get sha256
wget https://files.pythonhosted.org/packages/88/f7/a0b045408b6447c2e64b6ccc3cd75ebdf4e689fc02fbb522966c86f34bfc/georasters-0.5.17.tar.gz
openssl dgst -sha256 ./georasters-0.5.17.tar.gz
