set(SRCCCORE daceaux.c dacebasic.c dacecompat.c daceerror.c daceeval.c daceinit.c daceio.c dacemath.c dacememory.c dacenorm.c)
include_directories(include)

add_library(c_core OBJECT ${SRCCCORE})
set_target_properties(c_core PROPERTIES POSITION_INDEPENDENT_CODE True)
target_compile_definitions(c_core PRIVATE DACE_DLL=1)
target_compile_definitions(c_core PRIVATE DACE_DLL_EXPORTS=1)

add_library(c_core_s OBJECT ${SRCCCORE})

install(DIRECTORY ${PROJECT_SOURCE_DIR}/core/include/ DESTINATION include/ COMPONENT headers FILES_MATCHING PATTERN "*.h")
