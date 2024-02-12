# Test
python setup.py test

# To upload to PyPy
#python setup.py register
#python setup.py sdist bdist_wheel bdist_wininst upload
python setup.py sdist
python setup.py sdist bdist_wheel --universal
twine upload --repository georasters dist/*

# To get sha256
wget https://files.pythonhosted.org/packages/62/c3/d00378a2a1f294a97c38024041a45489a8d22fbf79ed25421fdacccf784a/georasters-0.5.18.tar.gz
openssl dgst -sha256 ./georasters-0.5.18.tar.gz
