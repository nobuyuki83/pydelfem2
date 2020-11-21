cmake_minimum_required(VERSION 3.0)

option(USE_HEADERONLY "USE_DFM2_HEADERONLY" OFF)

#########################################################

enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
if(MSVC)
  set(CMAKE_CXX_FLAGS "/W4  -O2 \
    /wd4530 /wd4577 /wd4100 /wd4458 /wd4244 /wd4701 \
    /wd4267 /wd4305 /wd4389 /wd4456 /wd4702 /wd4996 \
    /wd4505 /wd4245 /wd4457 /wd4800 /wd4838 /wd4055 \
    /wd4054")
else()
  set(CMAKE_CXX_FLAGS "-O2 -Wall -Wextra \
    -Wno-deprecated-declarations \
    -Wno-unused-parameter \
    -Wno-unused-function \
    -Wno-unused-variable")
endif()

################################

set(PROJ_NAME c_gl)
project(${PROJ_NAME})

set(DELFEM2_INCLUDE_DIR "../../delfem2/include")
add_definitions(-DDFM2_HEADER_ONLY=ON)

set(GLAD_DIR "../../delfem2/3rd_party/glad3")
get_filename_component(GLAD_DIR ${GLAD_DIR} ABSOLUTE)

set(PYBIND11_DIR ../../pybind11)
add_subdirectory(${PYBIND11_DIR} pybind11)

find_package(OpenGL REQUIRED)

################################

set(INPUT_SRC
    ${GLAD_DIR}/include/glad/glad.h     
    ${GLAD_DIR}/src/glad.c
    py_sampler.cpp
    py_texture.cpp
    main.cpp
)

include_directories(
    ${OPENGL_INCLUDE_DIR}
    ${GLAD_DIR}/include
    ${DELFEM2_INCLUDE_DIR}
    ${GLAD_DIR}/include
    ${PYBIND11_DIR}/include)

pybind11_add_module(${PROJ_NAME}
  ${INPUT_SRC} )

target_link_libraries(${PROJ_NAME} PRIVATE
    ${OPENGL_LIBRARY})

file(TO_NATIVE_PATH "${PROJECT_BINARY_DIR}" NATIVE_PROJECT_BINARY_DIR)
file(TO_NATIVE_PATH "${PROJECT_SOURCE_DIR}" NATIVE_PROJECT_SOURCE_DIR)

if(MSVC)
add_custom_command(TARGET ${PROJ_NAME}
POST_BUILD
  COMMAND copy "${NATIVE_PROJECT_BINARY_DIR}\\${CMAKE_CFG_INTDIR}\\c_gl.*.pyd"   ${NATIVE_PROJECT_SOURCE_DIR}\\..\\..\\PyDelFEM2\\gl\\
  COMMENT "Post Event"
)
else()
add_custom_command(TARGET ${PROJ_NAME}
POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${PROJECT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/c_gl.*
  ${PROJECT_SOURCE_DIR}/../../PyDelFEM2/gl
  COMMENT "Post Event"
)
endif()
