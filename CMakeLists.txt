cmake_minimum_required(VERSION 3.3)
project(ls_reduction)

### c++11 standards
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -std=c++11 -I/usr/local/include")

if(TARGET SCIP::SCIP)
  # find package by SCIP PATH
  find_package(SCIP CONFIG PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()

include_directories(${SCIP_INCLUDE_DIRS})

add_executable(ls_reduction
   src/main.cpp)

target_link_libraries(ls_reduction ${SCIP_LIBRARIES})
