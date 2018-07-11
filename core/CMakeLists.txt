if(APPLE)
  cmake_policy(SET CMP0042 NEW)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
  set(CMAKE_INSTALL_NAME_DIR "@loader_path")
endif(APPLE)

file( GLOB SRCCCORE daceaux.c dacebasic.c dacecompat.c daceerror.c daceeval.c daceinit.c daceio.c dacemath.c dacememory.c dacenorm.c)
include_directories(include)
add_library(c_core OBJECT ${SRCCCORE})
set_target_properties(c_core PROPERTIES POSITION_INDEPENDENT_CODE True)
add_library(c_core_s OBJECT ${SRCCCORE})

install(DIRECTORY ${PROJECT_SOURCE_DIR}/core/include/ DESTINATION include COMPONENT headers)