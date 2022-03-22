#
# Compile i-eel source files
#
# This file is a part of the i-eel code
#
set(
    ieel_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/Fiber.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Flow.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fiber2D.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Flow2D.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/QLearning.cpp
)

set(
    ieel_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/Fiber.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Flow.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Fiber2D.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Flow2D.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/QLearning.hpp
)

# Define the target with appropriate dependencies
install_ieel_library(ieel)
target_include_directories(ieel SYSTEM INTERFACE ${ARMA_INCLUDE_DIR})
target_link_arma(ieel PUBLIC)

if (NETCDF_FOUND)
    target_include_directories(ieel SYSTEM PRIVATE ${NETCDF_INCLUDE_DIRS})
    target_link_libraries(ieel PRIVATE ${NETCDF_LIBRARIES})
endif ()

if (WITH_HDF5CXX)
    target_include_directories(ieel SYSTEM PRIVATE ${HDF5_INCLUDE_DIRS})
    target_link_libraries(ieel PRIVATE ${HDF5_LIBRARIES})
endif ()

# Install header files
install(FILES ${ieel_HEADERS} DESTINATION include/channelflow)
install(FILES ${PROJECT_BINARY_DIR}/ieel/config.h DESTINATION include/ieel)

