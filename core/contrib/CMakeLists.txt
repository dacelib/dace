# all of these are only used internally by the DACE core and not ever exposed
# hence nothing is installed and no headers are made public
set(SRCCCONTRIB libf2c.c dgamma.c ribesl.c rjbesl.c rkbesl.c rybesl.c psi.c zeta.c)
include_directories(include)

add_library(c_core_contrib_s OBJECT ${SRCCCONTRIB})
add_library(c_core_contrib OBJECT ${SRCCCONTRIB})
set_target_properties(c_core_contrib PROPERTIES POSITION_INDEPENDENT_CODE True)
