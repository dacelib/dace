file(GLOB SRCCCORE daceaux.c dacebasic.c dacecompat.c daceerror.c daceeval.c daceinit.c daceio.c dacemath.c dacememory.c dacenorm.c)
include_directories(include)

add_library(c_core OBJECT ${SRCCCORE})
set_target_properties(c_core PROPERTIES POSITION_INDEPENDENT_CODE True)
target_compile_definitions(c_core PRIVATE ${DACE_OPT_FLAGS})

add_library(c_core_s OBJECT ${SRCCCORE})
target_compile_definitions(c_core_s PRIVATE ${DACE_OPT_FLAGS})

install(DIRECTORY ${PROJECT_SOURCE_DIR}/core/include/ DESTINATION include/ COMPONENT headers)
