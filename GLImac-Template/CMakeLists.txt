cmake_minimum_required(VERSION 2.8)

project(ImacGL)

add_definitions(-std=c++17 -ggdb3)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)

find_package(SDL REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
find_package(Eigen3 REQUIRED)




# Pour gérer un bug a la fac, a supprimer sur machine perso:
# set(OPENGL_LIBRARIES /usr/lib/x86_64-linux-gnu/libGL.so.1)

include_directories(${SDL_INCLUDE_DIR} ${OPENGL_INCLUDE_DIR} ${GLEW_INCLUDE_DIR} glimac/include third-party/include ${EIGEN3_INCLUDE_DIR})

set(ALL_LIBRARIES glimac ${SDL_LIBRARY} ${OPENGL_LIBRARIES} ${GLEW_LIBRARY})

add_subdirectory(glimac)

add_subdirectory(src)
