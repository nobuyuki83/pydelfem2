#################################
# download & build submodules

git submodule update --init 
git submodule foreach git checkout master
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

echo "###################"
echo "# setup environment"
python3 -m venv myenv
source myenv/bin/activate
python3 -m pip install --upgrade pip
# 3D utility
python3 -m pip install -U PyOpenGL 
python3 -m pip install -U glfw

pip3 uninstall PyDelFEM2 -y
pip3 uninstall PyDelFEM2 -y
pip3 install -e .
python3 setup.py test
