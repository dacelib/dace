# C++ interface source files
file( GLOB SRCCXXIF *.cpp )
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Build object files for dynamic and static libraries of c core
add_library(dacecxx_s OBJECT ${SRCCXXIF})
add_library(dacecxx OBJECT ${SRCCXXIF})
set_target_properties(dacecxx PROPERTIES POSITION_INDEPENDENT_CODE True)

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  set_target_properties(dacecxx PROPERTIES COMPILE_FLAGS "/D DACE_DLL /D DACE_DLL_EXPORTS")
else(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  set_target_properties(dacecxx PROPERTIES COMPILE_FLAGS "-D DACE_DLL -D DACE_DLL_EXPORTS")
endif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")

# statically and dynamically linked C++ interface library with the corresponding c core object files compiled in
add_library(dace SHARED $<TARGET_OBJECTS:c_core> $<TARGET_OBJECTS:dacecxx>)
add_library(dace_s STATIC $<TARGET_OBJECTS:c_core_s> $<TARGET_OBJECTS:dacecxx_s>)
install(TARGETS dace LIBRARY DESTINATION lib COMPONENT libraries
            RUNTIME DESTINATION lib COMPONENT libraries 
            ARCHIVE DESTINATION lib COMPONENT libraries
            PUBLIC_HEADER DESTINATION include COMPONENT headers)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/interfaces/cxx/include/DA DESTINATION include COMPONENT headers)