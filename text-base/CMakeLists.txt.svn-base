cmake_minimum_required(VERSION 1.0)
project(CrossTalkZ CXX)


# <section> some stuff that might be useful.... 
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS_RELEASE "-g -O3 -Wall")
set(CMAKE_CXX_COMPILER "g++")

SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "3")
SET(CPACK_PACKAGE_VERSION_PATCH "3")
set(PACKAGE_VERSION ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH})


set(CPACK_GENERATOR RPM;DEB)
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libxml2 (>= 2.6.16 )")
set(CPACK_RPM_PACKAGE_REQUIRES "libxml2 (>= 2.6.16 )")

set(CPACK_SOURCE_GENERATOR TGZ;ZIP)
set(CPACK_RPM_PACKAGE_ARCHITECTURE "i386")
set(CPACK_DEB_PACKAGE_ARCHITECTURE "i386")
set(CPACK_PACKAGE_CONTACT "Erik.Sonnhammer@sbc.su.se") 
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "CrossTalkZ statistics tool")
set(CPACK_PACKAGE_VENDOR "")
set(CPACK_PACKAGE_NAME "CrossTalkZ")

include(CPack)

#macro(my_find_program progname)
#  find_program(${progname}_EXECUTABLE ${progname})
#  if(NOT ${progname}_EXECUTABLE)
#    message(FATAL_ERROR  "${progname} not found" )
#  endif()
#endmacro(my_find_program)

#my_find_program(gengetopt)
#
#find_program(svnversion_EXECUTABLE svnversion)
#if(svnversion_EXECUTABLE)
#  execute_process( COMMAND ${svnversion_EXECUTABLE} -n ${CMAKE_SOURCE_DIR} OUTPUT_VARIABLE SVN_VERSION_OUTPUT)
#  if(SVN_VERSION_OUTPUT MATCHES "exported")
#    set(SVN_REVISION UNKNOWN)
#  else()
#    set(SVN_REVISION "${SVN_VERSION_OUTPUT}")
#  endif()
#else()
#    set(SVN_REVISION UNKNOWN)
#endif()
#
#if(SVN_REVISION)
#  set(SVN_REVISION_STRING "( subversion revision: ${SVN_REVISION} )")
#endif()
#
#set(VERSION_FOR_HELP_OUTPUT "${PACKAGE_VERSION} ${SVN_REVISION_STRING}")
#
#
## </section>
#
#
## Read more here: http://www.cmake.org/cmake/help/cmake-2-8-docs.html
#
## let cmake do text replacement of cmake variables for all files having a filename ending with ".cmake"
## We use it for putting a version number in main.cpp
#file(GLOB configurefiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cmake )
#foreach(i ${configurefiles})
#  string(REPLACE ".cmake" "" myReplFile "${i}")
#  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${myReplFile}.cmake ${CMAKE_CURRENT_SOURCE_DIR}/${myReplFile} @ONLY)
#endforeach()

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

find_package(LibXml2 REQUIRED)
include_directories( ${LIBXML2_INCLUDE_DIR} )
add_definitions( ${LIBXML2_DEFINITIONS} )
add_executable(CrossTalkZ boostgraphio.cpp main.cpp crosstalkz.cpp )
set(Boost_USE_STATIC_LIBS   ON)
set(Boost_USE_MULTITHREADED ON)
find_package( Boost COMPONENTS program_options graph REQUIRED )
include_directories(${Boost_INCLUDE_DIRS})

#The following is a boost hack if the libraries arent found
#set(Boost_LIBRARIES "optimized;/usr/local/lib/libboost_program_options.a;debug;/usr/local/lib/libboost_program_options.a;optimized;/usr/local/lib/libboost_graph.a;debug;/usr/local/lib/libboost_graph.a")

target_link_libraries(CrossTalkZ ${LIBXML2_LIBRARIES} ${Boost_LIBRARIES})
install(TARGETS CrossTalkZ DESTINATION bin)

#message("${Boost_LIBRARIES}")
#message("${Boost_INCLUDE_DIRS}")
# From this point I include some stuff that might be useful. /Erik Sjolund


