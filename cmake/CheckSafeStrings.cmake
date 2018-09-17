# Check if safe C string functions with _s suffix are available
function(check_safe_strings HAVE_SAFE_STRINGS_VAR)
   include(CheckSymbolExists)
   set(CMAKE_REQUIRED_DEFINITIONS "-D__STDC_WANT_LIB_EXT1__=1")
   check_symbol_exists(strncpy_s string.h HAVE_SAFE_STRINGS)
   set(${HAVE_SAFE_STRINGS_VAR} ${HAVE_SAFE_STRINGS} PARENT_SCOPE)
endfunction()
