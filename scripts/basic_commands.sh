# Test
python setup.py test

# To upload to PyPy
#python setup.py register
#python setup.py sdist bdist_wheel bdist_wininst upload
python setup.py sdist
python setup.py sdist bdist_wheel --universal
twine upload dist/*

# To get sha256
wget https://pypi.python.org/packages/24/b2/5aa7e9f9e389017396d2e1a45032decd775376b74c4f0391d461389f55a6/georasters-0.5.11.tar.gz#md5=cd6ebfc5d143df753ab42be6d127c4fd
openssl dgst -sha256 ./georasters-0.5.11.tar.gz