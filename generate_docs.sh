cd docs/
sphinx-apidoc -force -o ./api/. ../src
# Recall to change the maxdepth in the pyforce.rst file
rm api/modules.rst
make clean
make html
