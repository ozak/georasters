# Test
python setup.py test

# To upload to PyPy
#python setup.py register
#python setup.py sdist bdist_wheel bdist_wininst upload
python setup.py sdist
python setup.py sdist bdist_wheel --universal
twine upload dist/*

# To get sha256
wget https://files.pythonhosted.org/packages/a0/e9/63d9e6a03d4214a15c3b87f49e80a9d49613015396702677ec1aa013e443/georasters-0.5.12.tar.gz
openssl dgst -sha256 ./georasters-0.5.12.tar.gz