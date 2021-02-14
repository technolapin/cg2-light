cmake_minimum_required(VERSION 3.5)

project(ga_view)

# set the module path for "findLib.cmake"
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/modules/")

# add a "output" directory in the project binary directory
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output)

# require Eigen
find_package(Eigen3 REQUIRED)
if(${EIGEN3_FOUND})
    message(STATUS "lib EIGEN3 found")
    message(STATUS "  version " ${EIGEN3_VERSION_STRING})
    message(STATUS "  include " ${EIGEN3_INCLUDE_DIR})
endif()

# require lib c3ga
find_package(C3GA)
if(${C3GA_FOUND})
    message(STATUS "lib C3GA found")
    message(STATUS "  include " ${C3GA_INCLUDE_DIRS})
    message(STATUS "  lib     " ${C3GA_LIBRARIES})
else()
    message(STATUS "lib C3GA should be set manually, please uncomment the commented lines in the CMakeLists.txt file")

    # linux users
    #set(C3GA_INCLUDE_DIRS "/home/[full path to]/garamon_c3ga/src/")
    #set(C3GA_LIBRARIES "[full path to]/garamon_c3ga/build/libc3ga.so")

    # mac users
    #set(C3GA_INCLUDE_DIRS "/home/[full path to]/garamon_c3ga/src/")
    #set(C3GA_LIBRARIES "[full path to]/garamon_c3ga/build/libc3ga.dylib")

    # windows users
    #set(C3GA_INCLUDE_DIRS "C:\[full path to]/garamon_c3ga/src/")
    #set(C3GA_LIBRARIES "[full path to]/garamon_c3ga/build/libc3ga.dll")

    # for manual setting, check if valid
    if(EXISTS "${C3GA_INCLUDE_DIRS}")
        message(STATUS "include for lib c3ga found: ${C3GA_INCLUDE_DIRS}" )
    else()
        message(FATAL_ERROR "include for lib c3ga invalid or not found: ${C3GA_INCLUDE_DIRS}" )
    endif()

    # for manual setting, check if valid
    if(EXISTS "${C3GA_LIBRARIES}")
        message(STATUS "library c3ga found: ${C3GA_LIBRARIES}" )
    else()
        message(FATAL_ERROR "library c3ga not found: ${C3GA_LIBRARIES}" )
    endif()
endif()


# compilation flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O2 -std=c++17")
#set(CMAKE_CXX_FLAGS_DEBUG_INIT "-g")

# files to compile
file(GLOB_RECURSE source_files src/*.cpp src/*.hpp)
message(STATUS "  files " ${source_files})









find_package(SDL REQUIRED)
find_package(OpenGL REQUIRED)

if(${SDL_FOUND})
    message(STATUS "lib SDL found")
    message(STATUS "  version " ${SDL_VERSION_STRING})
    message(STATUS "  include " ${SDL_INCLUDE_DIR})
endif()
if(${OpenGL_FOUND})
    message(STATUS "lib OpenGL found")
    message(STATUS "  version " ${OpenGL_VERSION_STRING})
    message(STATUS "  include " ${OpenGL_INCLUDE_DIR})
endif()

find_package(GLEW REQUIRED)
if (GLEW_FOUND)
    message(STATUS "lib GLEW found")
    message(STATUS "  version " ${GLEW_VERSION_STRING})
    message(STATUS "  include " ${GLEW_INCLUDE_DIR})
endif()

# includes
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${C3GA_INCLUDE_DIRS})
include_directories(${SDL_INCLUDE_DIR} ${OPENGL_INCLUDE_DIR} ${GLEW_INCLUDE_DIR} glimac/include third-party/include)



set(ALL_LIBRARIES glimac ${SDL_LIBRARY} ${OPENGL_LIBRARIES} ${GLEW_LIBRARY})


#link_libraries()

add_subdirectory(glimac)















# compile
add_executable(ga_view src/main.cpp)
#add_executable(ga_view src/Directory.cpp src/Directory.hpp src/Entry.hpp src/Entry.cpp src/main.cpp )


# link
target_link_libraries(ga_view ${C3GA_LIBRARIES} ${ALL_LIBRARIES})
