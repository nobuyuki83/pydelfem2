cmake_minimum_required(VERSION 3.0)

#########################################################

enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
if(MSVC)
  set(CMAKE_CXX_FLAGS "/W4 -O2 \
    /wd4530 /wd4577 /wd4100 /wd4458 /wd4244 /wd4702 /wd4996 \
    /wd4701 /wd4267 /wd4305 /wd4389 /wd4456 /wd4505 ")
else()
  set(CMAKE_CXX_FLAGS "-O2 -Wall -Wextra \
    -Wno-deprecated-declarations \
    -Wno-unused-parameter \
    -Wno-unused-function")
endif()

##########################################################

set(PROJ_NAME c_core)
project(${PROJ_NAME})

set(DELFEM2_INCLUDE_DIR "../../delfem2/include")
add_definitions(-DDFM2_HEADER_ONLY=ON)

set(PYBIND11_DIR        "../../pybind11")
add_subdirectory(${PYBIND11_DIR} pybind11)


########################

include_directories(
  ${DELFEM2_INCLUDE_DIR}
  ../../
  ../../delfem2/3rd_party/
)

set(DELFEM2_INC         "../../delfem2/include/delfem2")
set(INPUT_SRC
    ${DELFEM2_INC}/evalmathexp.h        
    ${DELFEM2_INC}/evalmathexp.cpp
    ${DELFEM2_INC}/mshtopoio.h          
    ${DELFEM2_INC}/mshtopoio.cpp
    ${DELFEM2_INC}/tinygltf/io_gltf.h
    ${DELFEM2_INC}/tinygltf/io_gltf.cpp
    ../py_funcs.h

    py_mshtopoio.cpp
    py_dynmsh.cpp
    py_primitive.cpp
    py_polyline.cpp
    py_field.cpp
    py_rigging.cpp
    py_fem.cpp
    py_ls.cpp    
    py_sdf.cpp

    main.cpp
)

###########################################

pybind11_add_module(${PROJ_NAME}
  ${INPUT_SRC} )

target_link_libraries(${PROJ_NAME}
  PRIVATE ${INPUT_LIBRARY})

file(TO_NATIVE_PATH "${PROJECT_BINARY_DIR}" NATIVE_PROJECT_BINARY_DIR)
file(TO_NATIVE_PATH "${PROJECT_SOURCE_DIR}" NATIVE_PROJECT_SOURCE_DIR)

if(WIN32)
  add_custom_command(TARGET ${PROJ_NAME}
  POST_BUILD
    COMMAND copy 
      "${NATIVE_PROJECT_BINARY_DIR}\\${CMAKE_CFG_INTDIR}\\c_core.*.pyd"   
      ${NATIVE_PROJECT_SOURCE_DIR}\\..\\..\\PyDelFEM2\\
    COMMENT "Post Event"
  )
else()
  add_custom_command(TARGET ${PROJ_NAME}
  POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy 
      ${PROJECT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/c_core.*
      ${PROJECT_SOURCE_DIR}/../../PyDelFEM2/
    COMMENT "Post Event"
  )
endif()
