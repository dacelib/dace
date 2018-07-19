cmake_minimum_required (VERSION 2.8.4)
include(CMakeDependentOption)

project(DACE C CXX)

set(DACE_VERSION_MAJOR "2")
set(DACE_VERSION_MINOR "0")
set(DACE_VERSION_PATCH "0")


# global OS specific settings
if(WIN32)
  cmake_policy(SET CMP0054 NEW)
elseif(APPLE)
  cmake_policy(SET CMP0042 NEW)
  cmake_policy(SET CMP0068 NEW)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
  set(CMAKE_INSTALL_NAME_DIR "@loader_path")
endif(WIN32)


# Set platform dependent C compiler flags for various compilers and build types
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -std=c99")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O2 -frounding-math -fsignaling-nans")
elseif(CMAKE_C_COMPILER_ID MATCHES ".*Clang")        # Also catch IDs like "AppleClang"
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -std=c99")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O2")
elseif(CMAKE_C_COMPILER_ID STREQUAL "Intel")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O2")
elseif(CMAKE_C_COMPILER_ID STREQUAL "MSVC")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /fp:strict")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /O2")
endif(CMAKE_C_COMPILER_ID STREQUAL "GNU")

# Set platform dependent C++ compiler flags for various compilers and build types
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -frounding-math -fsignaling-nans")
elseif(CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")        # Also catch IDs like "AppleClang"
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /fp:strict")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")
endif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")


# Boolean options
option(WITH_ALGEBRAICMATRIX "Include the AlgebraicMatrix type (experimental)" OFF)
option(WITH_PTHREAD "Enable threading via pthread library" OFF)

# Tri-state option for memory model and the associated choices
set(DACE_MEMORY_MODEL HYBRID CACHE STRING "DACE memory model")
set_property(CACHE DACE_MEMORY_MODEL PROPERTY STRINGS DYNAMIC STATIC HYBRID)
string(TOUPPER "${DACE_MEMORY_MODEL}" DACE_MEMORY_MODEL_U)

# Default build type option if not set already
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
       "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
       FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)  # options shown in CMake GUI
endif(NOT CMAKE_BUILD_TYPE)

# Debug flag
string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_U)
if(CMAKE_BUILD_TYPE_U STREQUAL "DEBUG")
  set(WITH_DEBUG 1)
endif(CMAKE_BUILD_TYPE_U STREQUAL "DEBUG")

# Pass options to source by configuring config.h file in core
configure_file(core/include/dace/config.h.in core/include/dace/config.h)


# globally set the include directory for the C core headers
include_directories("core/include")
include_directories("${CMAKE_BINARY_DIR}/core/include")
# globally set the include directory for the C++ interface headers
include_directories("interfaces/cxx/include")

# now process the subdirectories
add_subdirectory(core)
add_subdirectory(interfaces/cxx)


# Package generation setting
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "DACE Library")
set(CPACK_PACKAGE_NAME "DACE")
set(CPACK_PACKAGE_VENDOR "DACE Development Group")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION "${DACE_VERSION_MAJOR}.${DACE_VERSION_MINOR}.${DACE_VERSION_PATCH}")
set(CPACK_PACKAGE_VERSION_MAJOR "${DACE_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${DACE_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${DACE_VERSION_PATCH}")
set(CPACK_PACKAGE_CONTACT "dace@dinamicatech.com")
set(CPACK_PACKAGING_INSTALL_PREFIX /usr/local)
set(CPACK_PACKAGE_INSTALL_DIRECTORY "DACE ${DACE_VERSION_MAJOR}.${DACE_VERSION_MINOR}")
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}")
set(CPACK_COMPONENTS_ALL libraries headers)
set(CPACK_COMPONENT_LIBRARIES_DISPLAY_NAME "DACE C++ Library")
set(CPACK_COMPONENT_HEADERS_DISPLAY_NAME "DACE C++ Headers")
set(CPACK_COMPONENT_LIBRARIES_DESCRIPTION "DACE C++ shared library")
set(CPACK_COMPONENT_HEADERS_DESCRIPTION "DACE C++ headers files")
set(CPACK_COMPONENT_HEADERS_DEPENDS libraries)
#set(CPACK_MONOLITHIC_INSTALL 1)

# OS specific packager settings
if(WIN32)
  set(CPACK_GENERATOR "NSIS")
  set(CPACK_NSIS_EXTRA_PREINSTALL_COMMANDS "!include \\\"x64.nsh\\\"")
  set (CPACK_NSIS_EXTRA_INSTALL_COMMANDS "

    \\\${DisableX64FSRedirection}

    CopyFiles \\\"$INSTDIR\\\\lib\\\\dace.dll\\\" \\\"$SYSDIR\\\"

     ")
    set(CPACK_NSIS_EXTRA_UNINSTALL_COMMANDS "

    \\\${DisableX64FSRedirection}

    Delete \\\"$SYSDIR\\\\dace.dll\\\"

  ")
  unset(CPACK_PACKAGING_INSTALL_PREFIX)	# replace by whatever it needs to be
elseif(APPLE)
  set(CPACK_GENERATOR "TGZ;productbuild")
  set(CPACK_PACKAGE_VENDOR "DACEDevelopmentGroup")    # no spaces allowed for packagebuild
  set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README-MAC.rtf")
  file(COPY "${CMAKE_SOURCE_DIR}/LICENSE" DESTINATION "${CMAKE_BINARY_DIR}")
  file(RENAME "${CMAKE_BINARY_DIR}/LICENSE" "${CMAKE_BINARY_DIR}/LICENSE.txt")
  set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_BINARY_DIR}/LICENSE.txt")
elseif(UNIX)
  set(CPACK_GENERATOR "TGZ;DEB")
  set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
  # dependencies are set to default libraries on trusty (Ubuntu 14.04 LTS). Should be set to actual versions on build system.
  set(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6 (>= 2.19), libstdc++ (>= 4.8.4)")
  set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README-LINUX")
  install(FILES README-LINUX DESTINATION share/libdace/ )
endif(WIN32)

include(CPack)
