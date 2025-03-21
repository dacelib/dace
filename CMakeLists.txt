cmake_minimum_required (VERSION 3.6.0)
include(CMakeDependentOption)

# since MacOS 10.14 (XCode 10.0), default includes are no longer installed to /usr/include or /usr/local/include
# the SDK to use must be specified explicitly in CMAKE_OSX_SYSROOT as it seems at least cmake 3.12.4 does not find it automatically if left empty
# this must happen before project()
if(APPLE AND NOT DEFINED ENV{SDKROOT})
  execute_process(COMMAND xcodebuild -sdk macosx -version Path OUTPUT_VARIABLE RES OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(CMAKE_OSX_SYSROOT "${RES}" CACHE PATH "The product will be built against the headers and libraries located inside the indicated SDK.")
endif(APPLE AND NOT DEFINED ENV{SDKROOT})

project(DACE C CXX)

set(DACE_MAJOR_VERSION 2)
set(DACE_MINOR_VERSION 0)
set(DACE_PATCH_VERSION 1)


# add our own CMake modules to path
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")


# global OS specific settings
if(WIN32)
  cmake_policy(SET CMP0054 NEW)
elseif(APPLE)
  cmake_policy(SET CMP0042 NEW)
  cmake_policy(SET CMP0068 NEW)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
  set(CMAKE_INSTALL_NAME_DIR "@loader_path")
endif(WIN32)


# Set C/C++ standard (adapted https://github.com/ceres-solver/ceres-solver/blob/master/CMakeLists.txt)
# Respect user-specified CMAKE_C/C++_STANDARD or default to C99/C++11.
# Standard is a hard requirement and that option is hidden from CMake GUI
if(NOT CMAKE_C_STANDARD)
  set(CMAKE_C_STANDARD 99 CACHE STRING "C standard (minimum 99)" FORCE)
  set_property(CACHE CMAKE_C_STANDARD PROPERTY STRINGS 99 11 17)  # options shown in CMake GUI
endif(NOT CMAKE_C_STANDARD)
set(CMAKE_C_STANDARD_REQUIRED ON CACHE BOOL "")
mark_as_advanced(CMAKE_C_STANDARD_REQUIRED)

if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 14 CACHE STRING "C++ standard (minimum 14)" FORCE)
  set_property(CACHE CMAKE_CXX_STANDARD PROPERTY STRINGS 14 17 20)  # options shown in CMake GUI
endif(NOT CMAKE_CXX_STANDARD)
set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "")
mark_as_advanced(CMAKE_CXX_STANDARD_REQUIRED)

#message("-- Building with C${CMAKE_C_STANDARD} and C++${CMAKE_CXX_STANDARD}")


# Set generally understood compiler flags (enforcing -O2)
add_compile_options("$<$<CONFIG:Release>:-O2>" "$<$<CONFIG:RelWithDebInfo>:-O2>")

# Set compiler specific flags
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
  add_compile_options("-Wall" "$<$<CONFIG:Release>:-frounding-math>" "$<$<CONFIG:Release>:-fsignaling-nans>" "$<$<CONFIG:RelWithDebInfo>:-frounding-math>" "$<$<CONFIG:RelWithDebInfo>:-fsignaling-nans>")
elseif(CMAKE_C_COMPILER_ID MATCHES ".*Clang")        # Also catch IDs like "AppleClang"
  add_compile_options("-Wall")
elseif(CMAKE_C_COMPILER_ID STREQUAL "MSVC")
  add_compile_options("/W4" "$<$<CONFIG:Release>:/fp:strict>" "$<$<CONFIG:RelWithDebInfo>:/fp:strict>")
endif(CMAKE_C_COMPILER_ID STREQUAL "GNU")


# Boolean options
option(WITH_ALGEBRAICMATRIX "Include the AlgebraicMatrix type (experimental)" OFF)
option(WITH_PTHREAD "Enable threading via pthread library" OFF)

# Tri-state option for memory model and the associated choices
set(DACE_MEMORY_MODEL DYNAMIC CACHE STRING "DACE memory model")
set_property(CACHE DACE_MEMORY_MODEL PROPERTY STRINGS DYNAMIC STATIC HYBRID)
string(TOUPPER "${DACE_MEMORY_MODEL}" DACE_MEMORY_MODEL_U)

# Default build type option if not set already
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
       "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
       FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)  # options shown in CMake GUI
endif(NOT CMAKE_BUILD_TYPE)

# Debug flag
string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_U)
if(CMAKE_BUILD_TYPE_U STREQUAL "DEBUG")
  set(WITH_DEBUG 1)
endif(CMAKE_BUILD_TYPE_U STREQUAL "DEBUG")

# Static memory management default sizes
# If you do not know what these values are, you probably do not want to touch them!
if(DACE_MEMORY_MODEL_U STREQUAL "STATIC")
  # Pick any maximum NO, NV combination to support, then set these values like this:
  # NMMAX = (NO+NV)!/NO!/NV!         LIAMAX = (NO+1)^((NV+1)/2)
  # For example, for NO=10 and NV=10 set DACE_STATIC_NMMAX=184756 and DACE_STATIC_LIAMAX=161051
  # The following default values are sized for NO=10 and NV=7
  set(DACE_STATIC_NOMAX 10 CACHE STRING "(EXPERT) Maximum order." FORCE)
  set(DACE_STATIC_NVMAX 10 CACHE STRING "(EXPERT) Maximum number of variables." FORCE)
  set(DACE_STATIC_NMMAX 19448 CACHE STRING "(EXPERT) Maximum number of monomials [(NO+NV)!/NO!/NV!]." FORCE)
  set(DACE_STATIC_LIAMAX 161051 CACHE STRING "(EXPERT) Maximum addressing array size [(NO+1)^((NV+1)/2)]." FORCE)
  set(DACE_STATIC_VAR_SIZE 100 CACHE STRING "(EXPERT) Maximum number of DA variables." FORCE)
  set(DACE_STATIC_MEM_SIZE 1944800 CACHE STRING "(EXPERT) Maximum stack size [NMMAX*VAR_SIZE]." FORCE)
endif(DACE_MEMORY_MODEL_U STREQUAL "STATIC")


# Run checks for platform-specific functions
include(CheckSafeStrings)
check_safe_strings(HAVE_SAFE_STRINGS)


# Pass options to source by configuring config.h file in core
configure_file(core/include/dace/config.h.in core/include/dace/config.h)


# globally set the include directory for the C core headers
include_directories("core/include")
include_directories("${PROJECT_BINARY_DIR}/core/include")
# globally set the include directory for the C++ interface headers
include_directories("interfaces/cxx/include")

# now process the subdirectories
add_subdirectory(core/contrib)
add_subdirectory(core)
add_subdirectory(interfaces/cxx)


# Package generation setting
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "DACE Library")
set(CPACK_PACKAGE_NAME "DACE")
set(CPACK_PACKAGE_VENDOR "DACE Development Group")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION "${DACE_MAJOR_VERSION}.${DACE_MINOR_VERSION}.${DACE_PATCH_VERSION}")
set(CPACK_PACKAGE_VERSION_MAJOR "${DACE_MAJOR_VERSION}")
set(CPACK_PACKAGE_VERSION_MINOR "${DACE_MINOR_VERSION}")
set(CPACK_PACKAGE_VERSION_PATCH "${DACE_PATCH_VERSION}")
set(CPACK_PACKAGE_CONTACT "dace@dinamicatech.com")
set(CPACK_PACKAGING_INSTALL_PREFIX /usr/local)
set(CPACK_PACKAGE_INSTALL_DIRECTORY "DACE ${DACE_MAJOR_VERSION}.${DACE_MINOR_VERSION}")
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(CPACK_COMPONENTS_ALL libraries headers cmake)
set(CPACK_COMPONENT_LIBRARIES_DISPLAY_NAME "DACE C++ Library")
set(CPACK_COMPONENT_HEADERS_DISPLAY_NAME "DACE C++ Headers")
set(CPACK_COMPONENT_CMAKE_DISPLAY_NAME "DACE CMake Package")
set(CPACK_COMPONENT_LIBRARIES_DESCRIPTION "DACE C++ shared library")
set(CPACK_COMPONENT_HEADERS_DESCRIPTION "DACE C++ headers files")
set(CPACK_COMPONENT_CMAKE_DESCRIPTION "DACE CMake package for easy integration into CMake infrastructure")
set(CPACK_COMPONENT_HEADERS_DEPENDS libraries)
set(CPACK_COMPONENT_CMAKE_DEPENDS headers)

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
  unset(CPACK_PACKAGING_INSTALL_PREFIX)    # replace by whatever it needs to be
elseif(APPLE)
  set(CPACK_GENERATOR "TGZ;productbuild")
  set(CPACK_PACKAGE_VENDOR "DACEDevelopmentGroup")    # no spaces allowed for packagebuild
  set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README-MAC.rtf")
  file(COPY "${PROJECT_SOURCE_DIR}/LICENSE" DESTINATION "${PROJECT_BINARY_DIR}")
  file(RENAME "${PROJECT_BINARY_DIR}/LICENSE" "${PROJECT_BINARY_DIR}/LICENSE.txt")
  set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_BINARY_DIR}/LICENSE.txt")
elseif(UNIX)
  set(CPACK_GENERATOR "TGZ;DEB")
  set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
  # dependencies are set to default libraries on trusty (Ubuntu 14.04 LTS). Should be set to actual versions on build system.
  set(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6 (>= 2.19), libstdc++6 (>= 4.8.4)")
  set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README-LINUX")
  install(FILES README-LINUX DESTINATION share/libdace/ )
endif(WIN32)

include(CPack)
