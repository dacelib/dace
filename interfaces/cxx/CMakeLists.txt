# C++ interface source files
set(SRCCXXIF AlgebraicVector.cpp DA.cpp DACEException.cpp DAFormatter.cpp MathExtension.cpp Monomial.cpp compiledDA.cpp)

if(WITH_ALGEBRAICMATRIX)
  list(APPEND SRCCXXIF "AlgebraicMatrix.cpp")
endif(WITH_ALGEBRAICMATRIX)

# Build object files for static and dynamic libraries of C++ interface
add_library(dacecxx_s OBJECT ${SRCCXXIF})
add_library(dacecxx OBJECT ${SRCCXXIF})
set_target_properties(dacecxx PROPERTIES POSITION_INDEPENDENT_CODE True)
if(WIN32)
  target_compile_definitions(dacecxx PRIVATE "DACE_API=__declspec(dllexport)")
  target_compile_definitions(dacecxx_s PRIVATE "DACE_API=")
endif(WIN32)

# statically and dynamically linked C++ interface library with the corresponding c core object files compiled in
add_library(dace SHARED $<TARGET_OBJECTS:c_core> $<TARGET_OBJECTS:c_core_contrib> $<TARGET_OBJECTS:dacecxx>)
add_library(dace_s STATIC $<TARGET_OBJECTS:c_core_s> $<TARGET_OBJECTS:c_core_contrib_s> $<TARGET_OBJECTS:dacecxx_s>)

target_include_directories(dace INTERFACE $<INSTALL_INTERFACE:include>)
target_include_directories(dace_s INTERFACE $<INSTALL_INTERFACE:include>)

if(WITH_PTHREAD)
  target_link_libraries(dace PUBLIC pthread)
  target_link_libraries(dace_s PUBLIC pthread)
endif(WITH_PTHREAD)

# library versioning
set_property(TARGET dace PROPERTY VERSION "${DACE_MAJOR_VERSION}.${DACE_MINOR_VERSION}.${DACE_PATCH_VERSION}")
set_property(TARGET dace PROPERTY SOVERSION ${DACE_MAJOR_VERSION})
set_property(TARGET dace PROPERTY INTERFACE_dace_MAJOR_VERSION ${DACE_MAJOR_VERSION})
set_property(TARGET dace APPEND PROPERTY COMPATIBLE_INTERFACE_STRING dace_MAJOR_VERSION)
set_property(TARGET dace_s PROPERTY VERSION "${DACE_MAJOR_VERSION}.${DACE_MINOR_VERSION}.${DACE_PATCH_VERSION}")
set_property(TARGET dace_s PROPERTY INTERFACE_dace_s_MAJOR_VERSION ${DACE_MAJOR_VERSION})
set_property(TARGET dace_s APPEND PROPERTY COMPATIBLE_INTERFACE_STRING dace_s_MAJOR_VERSION)

# explicitly set required minimum C++ standard so it carries through to exported targets and consumers
set_property(TARGET dace PROPERTY INTERFACE_COMPILE_FEATURES cxx_std_11)
set_property(TARGET dace_s PROPERTY INTERFACE_COMPILE_FEATURES cxx_std_11)

install(TARGETS dace EXPORT dace
        LIBRARY DESTINATION lib COMPONENT libraries
        RUNTIME DESTINATION lib COMPONENT libraries
        ARCHIVE DESTINATION lib COMPONENT libraries
        PUBLIC_HEADER DESTINATION include COMPONENT headers)
install(TARGETS dace_s EXPORT dace
        LIBRARY DESTINATION lib COMPONENT libraries
        RUNTIME DESTINATION lib COMPONENT libraries
        ARCHIVE DESTINATION lib COMPONENT libraries
        PUBLIC_HEADER DESTINATION include COMPONENT headers)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/interfaces/cxx/include/dace DESTINATION include
        COMPONENT headers FILES_MATCHING PATTERN "*.h")

# Automatically create cmake package files to find the DACE
include(CMakePackageConfigHelpers)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/dace/daceConfigVersion.cmake"
  VERSION "${DACE_MAJOR_VERSION}.${DACE_MINOR_VERSION}.${DACE_PATCH_VERSION}"
  COMPATIBILITY AnyNewerVersion
)

export(EXPORT dace
  FILE "${CMAKE_CURRENT_BINARY_DIR}/dace/dace.cmake"
  NAMESPACE dace::
)

configure_file("../../cmake/daceConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/dace/daceConfig.cmake"
  COPYONLY
)

# install cmake package files
install(EXPORT dace NAMESPACE dace:: DESTINATION lib/cmake/dace COMPONENT cmake)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/dace/daceConfig.cmake" "${CMAKE_CURRENT_BINARY_DIR}/dace/daceConfigVersion.cmake" DESTINATION lib/cmake/dace COMPONENT cmake)
