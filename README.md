[![Build Status](https://travis-ci.org/abgandar/dace.svg?branch=master)](https://travis-ci.org/abgandar/dace)

# DACE
The Differential Algebra Computational Toolbox.

## Requirements
To use the DACE library, you need a C++11 compatible C++ compiler. Almost all modern C++ compilers (GNU, Clang, MS Visual Studio) meet that requirement nowadays. We also highly recommend using CMake as a build system, and obviously Git as your source code management system.

On Windows, you need to download Microsoft Visual Studio (the free Community Edition is fine). On MacOS X, you need to install the Xcode package from Apple from the App Store. On Linux, you can install the compiler toolchain through your package system if they don't already come pre-installed.

To build Windows DACE library installer packages, you additionally need the NSIS installer compiler.

## Getting started
To get started, either download one of the pre-built release packages for your platform (see the Releases tab above), or just clone the entire repository and build the library from scratch (which really is very easy).

To use the DACE, we recommend using CMake. New DACE releases past version 2.0.1 come with a CMake package system plugin that allows very simple use of the DACE within the CMake framework. To get started, have a look at the Tutorials directories. Both of the Tutorials there are independent projects with their own stand-alone CMake file you can use as a starting point for your own programs.

## Building the DACE
To build your own version of the DACE library, simply clone this repository:
```
git clone "https://github.com/dacelib/dace.git"
```
Then create a build directory and run cmake, and then make to compile everything:
```
cd dace
mkdir build
cd build
cmake ..
make
```
After some compiling, you should have a sparkling new dace library ready for use. To install it directly into your system, you can just use:
```
sudo make install
```
The ```sudo``` is there to give you the required permissions to install into your system directories (usually ```/usr/local```).

## Running the Tutorials
To test if your installation was successful, you can try to build one of the Tutorials:
```
cd ../Tutorials/Tutorial1
mkdir build
cd build
cmake ..
make
```
This should automatically compile all tutorials, ready to be run by you using e.g. the command ```./Example1```.
