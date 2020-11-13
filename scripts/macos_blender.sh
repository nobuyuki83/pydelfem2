echo "###################"
echo "# blender"

version_blender=$(/Applications/blender.app/Contents/MacOS/blender --version)
version_blender=$(cut -d' ' -f 2 <<< ${version_blender})
version_blender=${version_blender%.*}
echo "blender version: ${version_blender}"

path_blender=/Applications/blender.app/Contents/MacOS/Blender
path_python=$(dirname ${path_blender})/../Resources/${version_blender}/python/bin
path_python=$(find ${path_python} -maxdepth 1 -type f -name 'python*')
echo "path python: ${path_python}"

version_python=$(${path_python} --version)
version_python=$(cut -d' ' -f 2 <<< ${version_python})
version_python=${version_python%.*}
echo "version python: ${version_python}"

brew install "python@${version_python}"
path_include_from=/usr/local/opt/python@${version_python}/Frameworks/Python.framework/Headers
path_include_to=/Applications/Blender.app/Contents/Resources/${version_blender}/python/include/python${version_python}m

echo ${path_include_from}
echo ${path_include_to}
cp -r ${path_include_from}/*.h ${path_include_to}
cp -r ${path_include_from}/internal ${path_include_to}

${path_python} -m ensurepip     # make sure pip is installed
${path_python} -m pip install --upgrade pip # upgrade pip

${path_python} -m pip install -e . # install delfem2


${path_blender} -b -P examples_blender/01_simple.py
open examples_blender/out/01_out.png
${path_blender} -b -P examples_blender/02_envtex.py
open examples_blender/out/02_out.png
${path_blender} -b -P examples_blender/10_mesh.py
open examples_blender/out/10_out.png