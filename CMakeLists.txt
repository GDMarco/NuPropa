cmake_minimum_required(VERSION 2.8.12)

# Change name of the library to something meaningful!
set(PLUGINNAME "nupropa")
project(${PLUGINNAME} C CXX)

message(STATUS "Building NuPropa...")

find_package(PythonInterp)
find_package(PythonLibs)
find_package(SWIG REQUIRED)

# require C++11
# Note: CRPropa requires C++11
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

# SWIG options
option(ENABLE_SWIG_BUILTIN "Use SWIG builtin option" OFF) 
set(crpropaSwigFile "crpropa.i")
set(swigModeFlag "")
if(ENABLE_SWIG_BUILTIN)
  set(crpropaSwigFile "crpropa-builtin.i")
  set(swigModeFlag "-builtin")
endif(ENABLE_SWIG_BUILTIN)

# get information on crpropa, maybe  CRPROPA_SWIG_PATH and CRPROPA_INSTALL_PREFIX work
execute_process(COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/python/findCRPropa.py" swig_interface OUTPUT_VARIABLE crpropaSwigPath)
execute_process(COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/python/findCRPropa.py" install_prefix OUTPUT_VARIABLE crpropaInstallPrefix)
message(STATUS "CRPropa includes : ${crpropaInstallPrefix}")

find_path(CRPROPA_SWIG_INTERFACE_FILES
  crpropa.i
  HINTS crpropa/ share/crpropa python/ ${crpropaSwigPath} 
  )

find_path(CRPROPA_INCLUDES
  CRPropa.h
  HINTS crpropa/ include/crpropa ${crpropaInstallPrefix}/include/
  )

find_library(CRPROPA_LIB NAMES crpropa libcrpropa HINTS crpropa/ lib/crpropa crpropa/lib ${crpropaInstallPrefix}/lib/)

message(STATUS "CRPropa SWIG interface files: ${CRPROPA_SWIG_INTERFACE_FILES}")
message(STATUS "CRPropa includes : ${CRPROPA_INCLUDES}")
message(STATUS "CRPropa lib: ${CRPROPA_LIB}")

# swig plugin
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}_wrap.cxx PROPERTIES GENERATED true )
add_custom_command( OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}_wrap.cxx
  COMMAND swig -builtin -c++ -python ${PY3} -I${CRPROPA_SWIG_INTERFACE_FILES} -I${CRPROPA_INCLUDES} -I${CMAKE_SOURCE_DIR}/include -I${PYTHON_INCLUDE_DIRS} ${CRPROPA_SWIG_DEFINES} -dirprot -o ${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}_wrap.cxx -outdir ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/python/${PLUGINNAME}.i)

include_directories(include ${CMAKE_CURRENT_BINARY_DIR}/include ${CRPROPA_INCLUDES} ${PYTHON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})

# build interface and source in one library, add more .cc...
add_library(${PLUGINNAME} MODULE 
  ${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}_wrap.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/NeutrinoNeutrinoInteraction.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/NeutrinoAntineutrinoInteraction.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/NeutrinoPhotonInteraction.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Channels.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/ChannelsBundle.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/RelativisticInteraction.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/ParticleData.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/NeutrinoOscillation.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/NeutrinoMixing.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/NeutrinoField.cc
)

set_source_files_properties(
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Channels.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ChannelsBundle.cc
    PROPERTIES COMPILE_FLAGS "-std=c++17"
)

# /Applications/CRPropa/NuNuInteractionv1/CRPropa3/libs/kiss
include_directories(include ${crpropaInstallPrefix}/libs/kiss)

# to change!
include_directories(include /Applications/CRPropa/ALPinist/CRPropa3/libs/eigen3)

target_link_libraries(${PLUGINNAME} ${CRPROPA_LIB} ${PYTHON_LIBRARIES})

set_target_properties(${PLUGINNAME} PROPERTIES PREFIX "")
set_target_properties(${PLUGINNAME} PROPERTIES OUTPUT_NAME "_${PLUGINNAME}")

execute_process(
        COMMAND ${PYTHON_EXECUTABLE} -c "import sys; from distutils import sysconfig; sys.stdout.write(
sysconfig.get_python_lib(1,0,prefix='${CMAKE_INSTALL_PREFIX}'))"
        OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )

# installation 
install(DIRECTORY "${CMAKE_SOURCE_DIR}/python/${PLUGINNAME}" DESTINATION "${PYTHON_SITE_PACKAGES}")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}.py" DESTINATION "${PYTHON_SITE_PACKAGES}/${PLUGINNAME}")
install(TARGETS ${PLUGINNAME} LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/${PLUGINNAME}")


