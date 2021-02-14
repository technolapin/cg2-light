file(GLOB HEADER_FILES [_a-zA-Z0-9]*.hpp c2ga/*.hpp)
file(GLOB SRC_FILES main.cpp)
file(GLOB_RECURSE SHADER_FILES [_a-zA-Z0-9]*.glsl)

foreach(SRC_FILE ${SRC_FILES})
    get_filename_component(FILE ${SRC_FILE} NAME_WE)
    get_filename_component(DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    set(OUTPUT ${DIR}_${FILE})
    add_executable(${OUTPUT} ${SRC_FILE} ${HEADER_FILES} ${SHADER_FILES})
    target_link_libraries(${OUTPUT} ${ALL_LIBRARIES})
endforeach()

file(COPY shaders DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
