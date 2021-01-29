# Install script for directory: /home/gautam/sbs_work/vgit4libsbsdg/libsbsdig

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/w/halla-scifs17exp/sbs/thir/vgit4libsbsdg/libsbsdig-install")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
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

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/sbsdig" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/sbsdig")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/sbsdig"
         RPATH "/w/halla-scifs17exp/sbs/thir/vgit4libsbsdg/libsbsdig-install/lib64:/site/12gev_phys/2.3/Linux_CentOS7.7.1908-x86_64-gcc4.8.5/root/6.14.04/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/w/halla-scifs17exp/sbs/thir/vgit4libsbsdg/libsbsdig/sbsdig")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/sbsdig" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/sbsdig")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/sbsdig"
         OLD_RPATH "/site/12gev_phys/2.3/Linux_CentOS7.7.1908-x86_64-gcc4.8.5/root/6.14.04/lib:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
         NEW_RPATH "/w/halla-scifs17exp/sbs/thir/vgit4libsbsdg/libsbsdig-install/lib64:/site/12gev_phys/2.3/Linux_CentOS7.7.1908-x86_64-gcc4.8.5/root/6.14.04/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/sbsdig")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/w/halla-scifs17exp/sbs/thir/vgit4libsbsdg/libsbsdig-install/db")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/w/halla-scifs17exp/sbs/thir/vgit4libsbsdg/libsbsdig-install" TYPE DIRECTORY FILES "/home/gautam/sbs_work/vgit4libsbsdg/libsbsdig/db")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/w/halla-scifs17exp/sbs/thir/vgit4libsbsdg/libsbsdig-install/example")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/w/halla-scifs17exp/sbs/thir/vgit4libsbsdg/libsbsdig-install" TYPE DIRECTORY FILES "/home/gautam/sbs_work/vgit4libsbsdg/libsbsdig/example")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE PROGRAM FILES "/w/halla-scifs17exp/sbs/thir/vgit4libsbsdg/libsbsdig/sbsdigenv.sh")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE PROGRAM FILES "/w/halla-scifs17exp/sbs/thir/vgit4libsbsdg/libsbsdig/sbsdigenv.csh")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/home/gautam/sbs_work/vgit4libsbsdg/libsbsdig/src/SBSDigAuxi.h"
    "/home/gautam/sbs_work/vgit4libsbsdg/libsbsdig/src/SBSDigBkgdGen.h"
    "/home/gautam/sbs_work/vgit4libsbsdg/libsbsdig/src/SBSDigGEMDet.h"
    "/home/gautam/sbs_work/vgit4libsbsdg/libsbsdig/src/SBSDigGEMPlane.h"
    "/home/gautam/sbs_work/vgit4libsbsdg/libsbsdig/src/SBSDigGEMSimDig.h"
    "/home/gautam/sbs_work/vgit4libsbsdg/libsbsdig/src/SBSDigPMTDet.h"
    "/home/gautam/sbs_work/vgit4libsbsdg/libsbsdig/src/SBSDigPMTSignal.h"
    "/home/gautam/sbs_work/vgit4libsbsdg/libsbsdig/src/g4sbs_data.h"
    "/home/gautam/sbs_work/vgit4libsbsdg/libsbsdig/src/g4sbs_tree.h"
    "/home/gautam/sbs_work/vgit4libsbsdg/libsbsdig/src/g4sbs_types.h"
    "/home/gautam/sbs_work/vgit4libsbsdg/libsbsdig/src/gmn_tree.h"
    "/home/gautam/sbs_work/vgit4libsbsdg/libsbsdig/src/sbsdig_LinkDef.h"
    "/home/gautam/sbs_work/vgit4libsbsdg/libsbsdig/src/g4sbs_types.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/w/halla-scifs17exp/sbs/thir/vgit4libsbsdg/libsbsdig/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
