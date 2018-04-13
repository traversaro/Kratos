set(CMAKE_INCLUDE_CURRENT_DIR ON)

# install the petsc
if(${INSTALL_PETSC} MATCHES ON)

  cmake_minimum_required (VERSION 2.8)

  message("**** compiling PETSc ****")

  set(PETSC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_modules/PETSc/)
  set(PETSC_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/external_libraries/PETSc)
  set(petsc_cmd "./configure")
  set(petsc_arg "--prefix=${PETSC_PREFIX}")
  message(STATUS "petsc cmd: ${petsc_cmd} ${petsc_arg}")
  execute_process(COMMAND ${petsc_cmd} ${petsc_arg}
  WORKING_DIRECTORY ${PETSC_DIR}
  RESULT_VARIABLE petsc_result
  OUTPUT_VARIABLE petsc_output)
  message(STATUS "petsc output[${petsc_result}]:\n${petsc_output}")

  set(petsc_cmd "make")
  set(petsc_arg "-j8")
  message(STATUS "petsc cmd: ${petsc_cmd} ${petsc_arg}")
  execute_process(COMMAND ${petsc_cmd} ${petsc_arg}
  WORKING_DIRECTORY ${PETSC_DIR}
  RESULT_VARIABLE petsc_result
  OUTPUT_VARIABLE petsc_output)
  message(STATUS "petsc output[${petsc_result}]:\n${petsc_output}")

  set(petsc_cmd "make")
  set(petsc_arg "install")
  message(STATUS "petsc cmd: ${petsc_cmd} ${petsc_arg}")
  execute_process(COMMAND ${petsc_cmd} ${petsc_arg}
  WORKING_DIRECTORY ${PETSC_DIR}
  RESULT_VARIABLE petsc_result
  OUTPUT_VARIABLE petsc_output)
  message(STATUS "petsc output[${petsc_result}]: ${petsc_output}")

endif(${INSTALL_PETSC} MATCHES ON)

message("**** configuring KratosPETScApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

# generate variables with the sources
set( KRATOS_PETSC_APPLICATION_CORE
	${CMAKE_CURRENT_SOURCE_DIR}/petsc_application.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/petsc_application_variables.cpp
)

## Generate variables with the sources
set( KRATOS_PETSC_APPLICATION_PYTHON_INTERFACE
        # Custom python
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/petsc_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
)

add_library(KratosPETScCore SHARED ${KRATOS_PETSC_APPLICATION_CORE})
target_link_libraries(KratosPETScCore PUBLIC KratosCore)
set_target_properties(KratosPETScCore PROPERTIES COMPILE_DEFINITIONS "PETSC_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosPETScApplication MODULE ${KRATOS_PETSC_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosPETScApplication PUBLIC KratosPETScCore)
set_target_properties(KratosPETScApplication PROPERTIES PREFIX "")


# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosPETScApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
set_target_properties(KratosPETScApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# install testing files
if(${INSTALL_TESTING_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME} FILES_MATCHING PATTERN "*.py" PATTERN  "*.json" PATTERN "*.mdpa" PATTERN ".svn" EXCLUDE)
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

if(USE_COTIRE MATCHES ON)
    cotire(KratosPETScCore)
    cotire(KratosPETScApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosPETScCore DESTINATION libs )
install(TARGETS KratosPETScApplication DESTINATION libs )

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/PETScApplication.py" DESTINATION KratosMultiphysics )
