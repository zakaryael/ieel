include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/compiler_helper.cmake)

# Defines a new target library ${name}_${type}, either shared or static.
#
# Installs it in ${CMAKE_INSTALL_PREFIX}/lib
#
# install_ieel_library(name)
#
# Args:
#     name: name of the library to be added as a target. Expects a
#         variable named ${name}_SOURCES to be defined and to contain
#         the list of sources.
#
function(install_ieel_library name)

    # Uncomment to debug calls to this function
    # message(STATUS "Building library: ${name}")

    # Create the target with the required type
    add_library(${name} ${${name}_SOURCES})

    target_include_directories(${name} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
        $<INSTALL_INTERFACE:include>
    )
    target_link_libraries(${name} PUBLIC basics)
    install(TARGETS ${name} DESTINATION lib)
    set_warning_flags(${name})

endfunction()

# Adds NETCDF to the target link libraries
#
function(target_link_netcdf name export_type)
    target_link_libraries(${name} ${export_type} ${NETCDF_LIBRARY})
endfunction()

# Adds ARMADILLO to the target link libraries
#
function(target_link_arma name export_type)
    target_link_libraries(${name} ${export_type} ${ARMA_LIBRARY})
endfunction()

# Adds and install an application that depends on the ieel library.
#
function(install_ieel_application name destination)
    set(app_name "${name}_app")
    add_executable(${app_name} ${name}.cpp)
    target_link_libraries(${app_name} PUBLIC ieel)
    set_warning_flags(${app_name})
    set_target_properties(${app_name} PROPERTIES OUTPUT_NAME ${name})
    if (destination)
        install(TARGETS ${app_name} DESTINATION ${destination})
    endif ()
endfunction()

