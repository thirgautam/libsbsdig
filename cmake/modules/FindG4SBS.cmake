###############################################################################
## Finds the SBS analysis package (SBS)
##
## Changelog:
##   * Mon Oct 08 2018 Juan Carlos Cornejo <cornejo@jlab.org>
##   - Initial find package module
###############################################################################


## Presently, we only support Podd being compiled in directory, and not
## installed anywhere else.

set(G4SBS_FOUND True)

set(G4SBS_find_library_name "g4sbsroot")
#if(DEFINED SBS_FIND_VERSION)
#  set(SBS_find_library_name ${SBS_find_library_name}.${SBS_FIND_VERSION})
#endif()

#Modifying the CMake build of SBS-offline : It is hereafter assumed that the ${SBS} environment variable points to the
# top-level install directory of Podd, with parallel subdirectories bin, include, and lib

## Find the compiled shared library
find_library(G4SBS_LIBRARY
  NAMES ${G4SBS_find_library_name}
  PATHS $ENV{G4SBS}/lib
  )
if(G4SBS_LIBRARY)
  ## Get the path to the library: under the new cmake build system for Podd,
  ## the include directory is parallel to the "lib" directory, under which the library is found
  get_filename_component(G4SBS_path ${G4SBS_LIBRARY} PATH)
  set(G4SBS_find_filenames cteqpdf.h)
  
  ## Ensure that we have the appropriate include directories
  foreach(G4SBS_filename ${G4SBS_find_filenames})
    find_path(_G4SBS_include_${G4SBS_filename}
      NAMES ${G4SBS_filename}	
      PATHS $ENV{G4SBS}/include
      )	    
    if(_G4SBS_include_${G4SBS_filename})
      list(APPEND G4SBS_INCLUDE_DIR ${_G4SBS_include_${G4SBS_filename}})
    else()
      message(FATAL_ERROR "Missing required header file ${G4SBS_filename} in G4SBS library. Please ensure that found path ${G4SBS_path} points to the correct directory. : ${_G4SBS_include_dir}")
      set(G4SBS_FOUND False)
    endif()
  endforeach()
else()
  set(G4SBS_FOUND False)
  message(FATAL_ERROR "G4SBS library not found. Please set your $G4SBS variable accordingly.")

  ##
endif()


include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(G4SBS DEFAULT_MSG G4SBS_FOUND G4SBS_LIBRARY G4SBS_INCLUDE_DIR)

#if(NOT ANALYZER_FOUND)
#  message(FATAL_ERROR "SBS not found. Set your the environmental variable SBS accordingly.")
#else()
#  set(ANALYZER_INCLUDE_DIR ${SBS_PATH})
#  message(ERROR_FATAL "SBS Analyzer found in ${SBS_PATH}")
#endif()


if(NOT G4SBS_CONFIG_EXEC)
  ## Only execute this if not already done so

  ## Wouldn't it be great if the user actually had an environmental variable
  ## set?
  if(DEFINED ENV{G4SBS})
    set(G4SBS_PATH $ENV{G4SBS})
  endif()
endif(NOT G4SBS_CONFIG_EXEC)

mark_as_advanced(G4SBS_FOUND)
