# Test
python setup.py test

# To upload to PyPy
#python setup.py register
#python setup.py sdist bdist_wheel bdist_wininst upload
python setup.py sdist
python setup.py sdist bdist_wheel --universal
twine upload dist/*

# To get sha256
wget https://files.pythonhosted.org/packages/bb/26/ba5d2443facf71d720d78c418348d05e3c6fe1c0510c556f0ff0d169f81d/georasters-0.5.16-py2.py3-none-any.whl
openssl dgst -sha256 ./georasters-0.5.16-py2.py3-none-any.whl
