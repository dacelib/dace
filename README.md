[![Build Status](https://travis-ci.org/dacelib/dace.svg?branch=master)](https://travis-ci.org/dacelib/dace)

# DACE
The Differential Algebra Computational Toolbox.

## Requirements
To use the DACE library, you need a C++11 compatible C++ compiler. Almost all modern C++ compilers (GNU, Clang, MS Visual Studio) meet that requirement nowadays. We also highly recommend using CMake as a build system, and obviously Git as your source code management system.

On Linux, you can install the compiler toolchain through your package system if they don't already come pre-installed. On Windows it is suggested to use the WSL environment and follow the linux procedure, or download Microsoft Visual Studio (the free Community Edition is fine). On MacOS X, you need to install the Xcode package from Apple from the App Store. 

To build Windows DACE library installer packages, you additionally need the NSIS installer compiler. (Deprecated, will be removed soon)

## Getting started
To get started just clone the entire repository and build the library from scratch (which really is very easy). You can use the current development version if you want the latests improvement and bug fixes or point to the Releases page for official released versions.

## Building the DACE
To build the DACE library, simply clone this repository:
```
git clone "https://github.com/dacelib/dace.git" dace
```
Then create a build directory, run cmake, and then cmake --build to compile everything:
```
mkdir dace-build
cmake -S dace/ -B dace-build/
cmake --build dace-build/
```
After some compiling, you should have a sparkling new dace library ready for use in the dace-build folder. 

Optionally you can install the freshly built DACE library directly into your system:
```
sudo cmake --install dace-build/
```
The ```sudo``` is there to give you the required permissions to install into your system directories (usually ```/usr/local```).

We have moved the tutorials to a different repository https://github.com/dacelib/dace-tutorials to keep the main repository clean. You can clone the tutorials repository and follow the instructions there to get started with DACE.

Also have a look at the [DACE Wiki pages](https://github.com/dacelib/dace/wiki) if you have further questions.

## Embedding the DACE in other projects
To use the DACE library without installing it locally, use the [FetchContent](https://cmake.org/cmake/help/v3.19/module/FetchContent.html) mechanism of CMake.

Include this code at the start of your CMakeList.txt file to have CMake automatically clone and build the DACE within your project and make the resulting libraries available as the `dace::dace_s` (static) and `dace::dace` (dynamic) targets.
```
include(FetchContent)
FetchContent_Declare(
  DACE
  GIT_REPOSITORY https://github.com/dacelib/dace.git
  GIT_TAG v2.1.0
)
FetchContent_MakeAvailable(DACE)
add_library(dace::dace ALIAS dace)
add_library(dace::dace_s ALIAS dace_s)
```

To build against the DACE, simply add your own executable and link it with one of these targets:
```
add_executable(my-code my-code.cpp)
target_link_libraries(my-code PUBLIC dace::dace_s)
```
