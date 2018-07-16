cmake_minimum_required (VERSION 2.8.4)
include(CMakeDependentOption)

if(WIN32)
  cmake_policy(SET CMP0054 NEW)
endif(WIN32)

project(DACE C CXX)

set(DACE_VERSION_MAJOR "2")
set(DACE_VERSION_MINOR "0")
set(DACE_VERSION_PATCH "0")

if(APPLE)
  cmake_policy(SET CMP0042 NEW)
  cmake_policy(SET CMP0068 NEW)
  SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
  set(CMAKE_INSTALL_NAME_DIR "@loader_path")
endif(APPLE)

option(ALGEBRAICMATRIX "Include the AlgebraicMatrix type" OFF)
option(BUILDDEBUG "Build with debugging enabled" OFF)
option(DACEPTHREAD "Enable Pthread library" OFF)
option(DACESTATIC "Static memory management version" OFF)
option(DACEDYNAMIC "Dynamic memory management version" OFF)

# Handle Options
if(DACESTATIC AND DACEDYNAMIC)
    message(FATAL_ERROR "It is not possible to select both static and dynamic memory management at the same time")
endif(DACESTATIC AND DACEDYNAMIC)

# Set a default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
       "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
       FORCE)
endif(NOT CMAKE_BUILD_TYPE)

set(DACE_OPT_FLAGS "")

if(DACEPTHREAD)
  # must distinguish Windows style /options and Unix style -options
  if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(DACE_OPT_FLAGS "${DACE_OPT_FLAGS} /D DACE_PTHREAD")
  else(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(DACE_OPT_FLAGS "${DACE_OPT_FLAGS} -D DACE_PTHREAD")
  endif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
endif(DACEPTHREAD)

if(DACESTATIC)
  # must distinguish Windows style /options and Unix style -options
  if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(DACE_OPT_FLAGS "${DACE_OPT_FLAGS} /D DACE_STATIC_MEMORY")
  else(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(DACE_OPT_FLAGS "${DACE_OPT_FLAGS} -D DACE_STATIC_MEMORY")
  endif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
endif(DACESTATIC)

if(DACEDYNAMIC)
  # must distinguish Windows style /options and Unix style -options
  if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(DACE_OPT_FLAGS "${DACE_OPT_FLAGS} /D DACE_DYNAMIC_MEMORY")
  else(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(DACE_OPT_FLAGS "${DACE_OPT_FLAGS} -D DACE_DYNAMIC_MEMORY")
  endif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
endif(DACEDYNAMIC)

if(ALGEBRAICMATRIX)
  # must distinguish Windows style /options and Unix style -options
  if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(DACE_OPT_FLAGS "${DACE_OPT_FLAGS} /D WITH_ALGEBRAICMATRIX")
  else(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(DACE_OPT_FLAGS "${DACE_OPT_FLAGS} -D WITH_ALGEBRAICMATRIX")
  endif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
endif(ALGEBRAICMATRIX)

# Set C compiler flags for various compilers and build types
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${DACE_OPT_FLAGS} -Wall -std=c99 -D WITH_DEBUG")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${DACE_OPT_FLAGS} -std=c99 -O2 -frounding-math -fsignaling-nans")
elseif(CMAKE_C_COMPILER_ID MATCHES ".*Clang")		# Also catch IDs like "AppleClang"
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${DACE_OPT_FLAGS} -std=c99 -Wall -D WITH_DEBUG")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${DACE_OPT_FLAGS} -std=c99 -O2")
elseif(CMAKE_C_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${DACE_OPT_FLAGS} -D WITH_DEBUG")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${DACE_OPT_FLAGS} -O2")
elseif(CMAKE_C_COMPILER_ID STREQUAL "MSVC")
# XXX: set the debug flags correctly
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${DACE_OPT_FLAGS} /DNDEBUG /fp:strict /D MSDOS /D WITH_DEBUG")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${DACE_OPT_FLAGS} /O2 /DNDEBUG /fp:strict /D MSDOS")
endif(CMAKE_C_COMPILER_ID STREQUAL "GNU")

# Set C++ compiler flags for various compilers and build types
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DACE_OPT_FLAGS} -Wall -std=c++11 -D WITH_DEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${DACE_OPT_FLAGS} -std=c++11 -O2 -frounding-math -fsignaling-nans")
elseif(CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")		# Also catch IDs like "AppleClang"
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DACE_OPT_FLAGS} -std=c++11 -Wall -D WITH_DEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${DACE_OPT_FLAGS} -std=c++11 -O2")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DACE_OPT_FLAGS} -D WITH_DEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${DACE_OPT_FLAGS} -O2")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
# XXX: set the debug flags correctly
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DACE_OPT_FLAGS} /DNDEBUG /fp:strict /D MSDOS /D WITH_DEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${DACE_OPT_FLAGS} /O2 /DNDEBUG /fp:strict /D MSDOS")
endif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

# globally set the include directory for the C core headers
include_directories("core/include")
# globally set the include directory for the C++ interface headers
include_directories("interfaces/cxx/include")

# now process the subdirectories
add_subdirectory(core)
add_subdirectory(interfaces/cxx)

# probably can be removed now:
# include directories for the DACE (used in SF CMake file)
set(${PROJECT_NAME}_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/interfaces/cxx/include CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)


#install(DIRECTORY ${PROJECT_SOURCE_DIR}/interfaces/cxx/include/DA DESTINATION include COMPONENT headers)

# do we need these?
#install(FILES LICENCE.txt DESTINATION share/dace/)
#install(FILES install.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE DESTINATION ./)

# Defaults for package generation
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "DACE Library")
SET(CPACK_PACKAGE_NAME "DACE")
SET(CPACK_PACKAGE_VENDOR "DACE Development Group")
#SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
SET(CPACK_PACKAGE_VERSION "${DACE_VERSION_MAJOR}.${DACE_VERSION_MINOR}.${DACE_VERSION_PATCH}")
SET(CPACK_PACKAGE_VERSION_MAJOR "${DACE_VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${DACE_VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${DACE_VERSION_PATCH}")
SET(CPACK_PACKAGE_CONTACT "dace@dinamicatech.com")
SET(CPACK_PACKAGING_INSTALL_PREFIX /usr/local)
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "DACE ${DACE_VERSION_MAJOR}.${DACE_VERSION_MINOR}")
SET(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}")
SET(CPACK_COMPONENTS_ALL libraries headers)
SET(CPACK_COMPONENT_LIBRARIES_DISPLAY_NAME "DACE C++ Library")
SET(CPACK_COMPONENT_HEADERS_DISPLAY_NAME "DACE C++ Headers")
SET(CPACK_COMPONENT_LIBRARIES_DESCRIPTION "DACE C++ shared library")
SET(CPACK_COMPONENT_HEADERS_DESCRIPTION "DACE C++ headers files") 
SET(CPACK_COMPONENT_HEADERS_DEPENDS libraries)
#SET(CPACK_MONOLITHIC_INSTALL 1)

# OS specific settings
IF(WIN32)
    SET(CPACK_GENERATOR "NSIS")
ELSEIF(APPLE)
    SET(CPACK_GENERATOR "productbuild")
ELSEIF(UNIX)
    SET(CPACK_GENERATOR "TGZ;DEB")
    install(FILES README-LINUX DESTINATION ./ )
ENDIF(WIN32)

# Package generation time overrides for some CPack settings
configure_file("${CMAKE_SOURCE_DIR}/CPackConf.cmake.in" "${CMAKE_BINARY_DIR}/CPackConf.cmake" @ONLY)
SET(CPACK_PROJECT_CONFIG_FILE "${CMAKE_BINARY_DIR}/CPackConf.cmake")

include(CPack)
