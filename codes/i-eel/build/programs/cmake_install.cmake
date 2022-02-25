# Install script for directory: /Users/zakaryaelkhiyati/swimmers/i-eel/codes/i-eel/programs

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/Users/zakaryaelkhiyati/swimmers/i-eel/codes/i-eel/runs")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/Library/Developer/CommandLineTools/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/zakaryaelkhiyati/swimmers/i-eel/codes/i-eel/build/programs/test_buckling")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_buckling" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_buckling")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/zakaryaelkhiyati/swimmers/i-eel/codes/i-eel/build/src"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_buckling")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/zakaryaelkhiyati/swimmers/i-eel/codes/i-eel/runs/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_buckling")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_buckling")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/zakaryaelkhiyati/swimmers/i-eel/codes/i-eel/build/programs/test_swim")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_swim" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_swim")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/zakaryaelkhiyati/swimmers/i-eel/codes/i-eel/build/src"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_swim")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/zakaryaelkhiyati/swimmers/i-eel/codes/i-eel/runs/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_swim")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_swim")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/zakaryaelkhiyati/swimmers/i-eel/codes/i-eel/build/programs/test_buckling2D")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_buckling2D" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_buckling2D")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/zakaryaelkhiyati/swimmers/i-eel/codes/i-eel/build/src"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_buckling2D")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/zakaryaelkhiyati/swimmers/i-eel/codes/i-eel/runs/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_buckling2D")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_buckling2D")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/zakaryaelkhiyati/swimmers/i-eel/codes/i-eel/build/programs/learn2D_cellflow")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/learn2D_cellflow" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/learn2D_cellflow")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/zakaryaelkhiyati/swimmers/i-eel/codes/i-eel/build/src"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/learn2D_cellflow")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/zakaryaelkhiyati/swimmers/i-eel/codes/i-eel/runs/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/learn2D_cellflow")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/learn2D_cellflow")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/zakaryaelkhiyati/swimmers/i-eel/codes/i-eel/build/programs/learn2D_reduced")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/learn2D_reduced" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/learn2D_reduced")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/zakaryaelkhiyati/swimmers/i-eel/codes/i-eel/build/src"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/learn2D_reduced")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/zakaryaelkhiyati/swimmers/i-eel/codes/i-eel/runs/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/learn2D_reduced")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/learn2D_reduced")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/zakaryaelkhiyati/swimmers/i-eel/codes/i-eel/build/programs/policy_gradient")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/policy_gradient" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/policy_gradient")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/zakaryaelkhiyati/swimmers/i-eel/codes/i-eel/build/src"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/policy_gradient")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/zakaryaelkhiyati/swimmers/i-eel/codes/i-eel/runs/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/policy_gradient")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/policy_gradient")
    endif()
  endif()
endif()

