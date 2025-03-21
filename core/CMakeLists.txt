set(SRCCCORE daceaux.c dacebasic.c dacecompat.c daceerror.c daceeval.c daceinit.c daceio.c dacemath.c dacememory.c dacenorm.c)
include_directories(include contrib/include)

add_library(c_core_s OBJECT ${SRCCCORE})
add_library(c_core OBJECT ${SRCCCORE})

set_target_properties(c_core PROPERTIES POSITION_INDEPENDENT_CODE True)
if(WIN32)
  target_compile_definitions(c_core PRIVATE "DACE_API=__declspec(dllexport)")
  target_compile_definitions(c_core_s PRIVATE "DACE_API=")
endif(WIN32)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/core/include/ DESTINATION include/ COMPONENT headers FILES_MATCHING PATTERN "*.h")
install(FILES ${PROJECT_BINARY_DIR}/core/include/dace/config.h DESTINATION include/dace COMPONENT headers)
