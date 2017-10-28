# Test
python setup.py test

# To upload to PyPy
#python setup.py register
#python setup.py sdist bdist_wheel bdist_wininst upload
python setup.py sdist bdist_wheel --universal
twine upload dist/*

