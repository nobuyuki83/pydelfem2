#################################
# download & build submodules

git submodule update --init --recursive
git submodule foreach git pull origin master

PATH_PYTHON=$(which python3)
echo ${PATH_PYTHON}

cd src_pybind/core
mkdir buildMake
cd buildMake
cmake -DPYTHON_EXECUTABLE:PATH=${PATH_PYTHON}  ..
cd ../../../
cd src_pybind/gl

mkdir buildMake
cd buildMake
cmake -DPYTHON_EXECUTABLE:PATH=${PATH_PYTHON}  ..
cd ../../../

pip3 uninstall PyDelFEM2 -y
pip3 uninstall PyDelFEM2 -y
pip3 install -e .
python3 setup.py test
