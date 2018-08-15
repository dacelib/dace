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
  target_compile_definitions(c_core PRIVATE "DACE_API=__declspec(dllexport)")
endif(WIN32)

# statically and dynamically linked C++ interface library with the corresponding c core object files compiled in
add_library(dace SHARED $<TARGET_OBJECTS:c_core> $<TARGET_OBJECTS:dacecxx>)
add_library(dace_s STATIC $<TARGET_OBJECTS:c_core_s> $<TARGET_OBJECTS:dacecxx_s>)

if(WITH_PTHREAD)
target_link_libraries(dace PUBLIC pthread)
target_link_libraries(dace_s PUBLIC pthread)
endif(WITH_PTHREAD)

install(TARGETS dace LIBRARY DESTINATION lib COMPONENT libraries
        RUNTIME DESTINATION lib COMPONENT libraries
        ARCHIVE DESTINATION lib COMPONENT libraries
        PUBLIC_HEADER DESTINATION include COMPONENT headers)
install(TARGETS dace_s LIBRARY DESTINATION lib COMPONENT libraries
        RUNTIME DESTINATION lib COMPONENT libraries
        ARCHIVE DESTINATION lib COMPONENT libraries
        PUBLIC_HEADER DESTINATION include COMPONENT headers)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/interfaces/cxx/include/dace DESTINATION include COMPONENT headers FILES_MATCHING PATTERN "*.h")
