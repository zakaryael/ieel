# ARMA_INCLUDE_DIR = armadillo.h
# ARMA_LIBRARIES = libarmadillo
# ARMA_FOUND = true if ARMA is found

# To find Armadillo on mac, use the following line
FIND_PATH(ARMA_LIBDIR libarmadillo.dylib HINTS $ENV{LD_LIBRARY_PATH} "/usr/local/lib/" "/opt/local/lib/")
# To find Armadillo on the NEF cluster, comment the above line and uncomment the following:
# FIND_PATH(ARMA_LIBDIR libarmadillo.so HINTS $ENV{LD_LIBRARY_PATH} "/home/jbec/Info/lib64/")

# if(ARMA_LIBDIR)
# else()
#   set(ARMA_LIBDIR "/usr/local/lib/")
# endif()
MESSAGE(STATUS "ARMA_LIB_DIR=${ARMA_LIBDIR}")

FIND_PATH(ARMA_INCLUDE_DIR arma_version.hpp HINTS "${ARMA_LIBDIR}/../include/armadillo_bits/")
# For mac:
set(ARMA_INCLUDE_DIR "/usr/local/include/")
# For NEF:
# set(ARMA_INCLUDE_DIR "/home/jbec/Info/include/")

FIND_LIBRARY(ARMA_LIBRARY armadillo HINTS ${ARMA_LIBDIR})
get_filename_component(ARMA_LIBDIR ${ARMA_LIBRARY} PATH)
include(LibFindMacros)
libfind_process(ARMA)

IF(ARMA_INCLUDE_DIR AND ARMA_LIBRARY)
  MESSAGE(STATUS "ARMA_INCLUDE_DIR=${ARMA_INCLUDE_DIR}")
  MESSAGE(STATUS "ARMA_LIBRARY=${ARMA_LIBRARY}")
  set(ARMA_FOUND TRUE)
ELSE()
  IF(ARMA_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "ARMA required, please specify its location.")
  ELSE()
    MESSAGE(STATUS      "ARMA was not found.")
  ENDIF()
ENDIF()


# MARK_AS_ADVANCED(ARMA_INCLUDE_DIR ARMA_LIBRARY)
