###############################################################################
## Finds the SBS analysis package (SBS)
##
## Changelog:
##   * Mon Oct 08 2018 Juan Carlos Cornejo <cornejo@jlab.org>
##   - Initial find package module
###############################################################################


## Presently, we only support Podd being compiled in directory, and not
## installed anywhere else.

set(SBS_FOUND True)

set(SBS_find_library_name "sbs")
#if(DEFINED SBS_FIND_VERSION)
#  set(SBS_find_library_name ${SBS_find_library_name}.${SBS_FIND_VERSION})
#endif()

#Modifying the CMake build of SBS-offline : It is hereafter assumed that the ${SBS} environment variable points to the
# top-level install directory of Podd, with parallel subdirectories bin, include, and lib

## Find the compiled shared library
find_library(SBS_LIBRARY
  NAMES ${SBS_find_library_name}
  PATHS $ENV{SBS}/lib
  )
if(SBS_LIBRARY)
  ## Get the path to the library: under the new cmake build system for Podd,
  ## the include directory is parallel to the "lib" directory, under which the library is found
  get_filename_component(SBS_path ${SBS_LIBRARY} PATH)
  set(SBS_find_filenames MPDModule.h)
  
  ## Ensure that we have the appropriate include directories
  foreach(SBS_filename ${SBS_find_filenames})
    find_path(_SBS_include_${SBS_filename}
      NAMES ${SBS_filename}	
      PATHS $ENV{SBS}/include
      )	    
    if(_SBS_include_${SBS_filename})
      list(APPEND SBS_INCLUDE_DIR ${_SBS_include_${SBS_filename}})
    else()
      message(FATAL_ERROR "Missing required header file ${SBS_filename} in SBS library. Please ensure that found path ${SBS_path} points to the correct directory. : ${_SBS_include_dir}")
      set(SBS_FOUND False)
    endif()
  endforeach()
else()
  set(SBS_FOUND False)
  message(FATAL_ERROR "SBS library not found. Please set your $SBS variable accordingly.")

  ##
endif()


include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(SBS DEFAULT_MSG SBS_FOUND SBS_LIBRARY SBS_INCLUDE_DIR)

#if(NOT ANALYZER_FOUND)
#  message(FATAL_ERROR "SBS not found. Set your the environmental variable SBS accordingly.")
#else()
#  set(ANALYZER_INCLUDE_DIR ${SBS_PATH})
#  message(ERROR_FATAL "SBS Analyzer found in ${SBS_PATH}")
#endif()


if(NOT SBS_CONFIG_EXEC)
  ## Only execute this if not already done so

  ## Wouldn't it be great if the user actually had an environmental variable
  ## set?
  if(DEFINED ENV{SBS})
    set(SBS_PATH $ENV{SBS})
  endif()
endif(NOT SBS_CONFIG_EXEC)

mark_as_advanced(SBS_FOUND)
