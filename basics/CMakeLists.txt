#
# Define basics as an INTERFACE (i.e. header only) library
#
# This permits to use it to set dependencies on targets, rather than
# setting global include directories
#

set(
    basics_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/Utilities.h
    ${CMAKE_CURRENT_SOURCE_DIR}/RunsIO.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Arglist.h
)

# Define the target with appropriate dependencies
add_library(basics INTERFACE)
target_include_directories(basics SYSTEM INTERFACE ${NETCDF_INCLUDE_DIRS} )
target_link_netcdf(basics INTERFACE)
target_include_directories(basics SYSTEM INTERFACE ${ARMA_INCLUDE_DIR})
target_link_arma(basics INTERFACE)

# Install header files
install(FILES ${basics_HEADERS} DESTINATION include/basics)
