# PyDelFEM2
A python binding of DelFEM2

<a href="http://doge.mit-license.org"><img src="http://img.shields.io/:license-mit-blue.svg"></a> 

| Ubuntu | 
| --- |
| ![CI_ubuntu](https://github.com/nobuyuki83/pydelfem2/workflows/CI_ubuntu/badge.svg) |


## About

PyDelFEM2 is a python binding of DelFEM2. PyDelFEM2 can run various types of FEM simulation just a 10-20 lines of codes. Here is the example of solving the Poisson's equation in a square domain.

```
import PyDelFEM2 as dfm2
import PyDelFEM2.gl.glfw

cad = dfm2.Cad2D()
cad.add_polygon(list_xy=[-1,-1, +1,-1, +1,0, +0,+0, 0,+1, -1,+1.0])
mesh,map_cad2mesh = cad.mesh(0.05)
fem = dfm2.FEM_Poisson(source=1.0)
fem.updated_topology(mesh)
npIdP = cad.points_edge([0,1,2,3], mesh.np_pos)
fem.ls.bc[npIdP] = 1
fem.solve()
field = dfm2.VisFEM_ColorContour(fem,"value")
dfm2.gl._glfw.winDraw3d([field])
```
The result of this code woud be the following window

![Poisson](docs/imgs/poisson.png)


## Files

+ Python examples
  + [delfem2/examples_py](examples_py) : examples using python
  + [delfem2/examples_pyqt](examples_pyqt) examples using PyQt gui library
  + [delfem2/examples_jupyter](examples_jupyter) : examples using Jupyter
  + [delfem2/examples_blender](examples_blender) : examples with Blender python scripting
+  Python tests
  + [delfem2/test_py](test_py) : test using python



## Installation 

The most recommended way to install PyDelFEM2, which is the python binding of DelFEM2, is to build it from the source code. The following command down load the source code and its C++ dependencies and build python modules and download its python dependencies.

```
git clone https://github.com/nobuyuki83/delfem2.git
git submodle update --init --recursive
pip3 install -e .
```

Here are some trouble shooting tips: 
- For Ubuntu if you don't have git install it with ```sudo apt-get install git```

- For Ubuntu, if you don't have pip installed, get it with ```sudo apt-get install python3-pip```

- The installation fails if OpenGL packages are missing. For Ubuntu, install them ```sudo apt-get install freeglut3-dev libglfw3-dev libglew-dev```



Alternatively, PyDelFEM2 can be also installed from the GitHub repository can be done with the command:
```
pip3 install git+https://github.com/nobuyuki83/delfem2
```

PyDelFEM2 can be installed from PyPL simply with the following command, however this is not recommended as the version in PyPL is not updated frequently.

```
pip3 install PyDelFEM2
```


PyDelFEM runs on Python3. Python2 is not supported. PyDelFEM2 depends on following awesome python packages:
- numpy
- glfw
- PyOpenGL
- PySide2

These dependencies are written in ```REQUIRED_PACKAGES``` in the setup.py, so they are automatically installed when installing the PyDelFEM2 pakage using the ```setup.py``` or ```pip3```.  