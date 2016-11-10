# compile everything except the main source file into a statis lib to be linked 
# to the main executable
add_library(${LIBRARY_NAME} STATIC
  rates.cpp
  sample_methods.cpp
  Xsection.cpp
  matrix_elements.cpp
)

set_target_properties(${LIBRARY_NAME} PROPERTIES PREFIX "")


# compile the actual executable
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} ${LIBRARY_NAME}  ${GSL_LIBRARIES} ${GSLCALAS_LIBRARIES} -lpthread)

install(TARGETS ${PROJECT_NAME} DESTINATION ${PROJECT_NAME})
