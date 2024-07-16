==================================
hermite_prolate_spheroid_mechanics
==================================

Hermite prolate spheroid 

Building the example
====================

Instructions on how to configure and build with CMake::

  git clone https://github.com/OpenCMISS-Examples/hermite_prolate_spheroid_mechanics.git
  cd hermite_prolate_spheroid_mechanics
  mkdir build
  cd build
  cmake -DOpenCMISS_INSTALL_ROOT=/path/to/opencmiss/install ../.
  make  # cmake --build . will also work here and is much more platform agnostic.

Running the example
===================

Explain how the example is run::

  ./src/fortran/hermite_prolate_spheroid_mechanics.F90

or maybe it is a Python only example::

  source /path/to/opencmisslibs/install/virtaul_environments/oclibs_venv_pyXY_release/bin/activate
  python src/python/hermite_prolate_spheroid_mechanics.py

where the XY in the path are the Python major and minor versions respectively.

Prerequisites
=============

None

License
=======

Apache 2.0 License
