cd docs/
sphinx-apidoc -force -o ./api/. ../src
rm api/modules.rst
make clean
make html
