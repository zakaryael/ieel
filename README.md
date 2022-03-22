#  i-eel

This C++ code integrates the dynamics of a thin flexible fiber and implement various machine-learning algorithms in order to optimise locomotion.

## Requirements

To compile :
* CMake (version 3.1 or higher)
* C++11 supportive compiler

Some libraries are also needed:
* HDF5(cxx) and NetCDF (version 4 or higher) for I/O
* SuperLU (linear algebra library for sparse systems) together with the wrapper Armadillo

## Installation

A Makefile to compile and install the code can be generated using the CMake build scripts.
CMake allows two build types: "release" and "debug".
The build type decides on compiler options and whether assertions are compiled.
Use "release" for optimal performance (e.g. for production runs) and
"debug" for optimal diagnostics (e.g. for code development).
Out of source builds are recommended.

```
mkdir build
cd build
cmake PATH_TO_SOURCE -DCMAKE_BUILD_TYPE=release (/debug) -DCMAKE_INSTALL_PREFIX=PATH_TO_INSTALL_DIR
make -j
make install
```
## Usage

After installation, several example programs are available in the directory bin/ of the folder 'PATH_TO_INSTALL_DIR'. The corresponding source codes are in the program/ directory.
