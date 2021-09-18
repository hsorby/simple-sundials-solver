

set(SISS_SRCS
  odesolver.cxx
  cvodesolver.cxx
  main.cxx
)

set(SISS_HEADERS
  odesolver.h
  cvodesolver.h
)

set(JSON_CPP_HEADERS
  json.hpp
)

include(${CMAKE_CURRENT_BINARY_DIR}/model_files.cmake)
set(CONFIGURED_MODEL_HEADER model_header.h)
configure_file(model_header.in.h ${CONFIGURED_MODEL_HEADER})
set(APPLICATION_CONFIG_HEADER application_config.h)
configure_file(application_config.in.h ${APPLICATION_CONFIG_HEADER})

add_executable(siss ${SISS_SRCS} ${SISS_HEADERS} ${JSON_CPP_HEADERS} ${MODEL_FILES} ${CONFIGURED_MODEL_HEADER} ${APPLICATION_CONFIG_HEADER})

target_link_libraries(siss PRIVATE SUNDIALS::cvodes)
target_include_directories(siss PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
set_target_properties(siss PROPERTIES
  CXX_STANDARD 17
  CXX_STANDARD_REQUIRED ON
)
