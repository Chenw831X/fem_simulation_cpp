cmake_minimum_required(VERSION 3.9)
project(fem3d)

# project source files
file(GLOB SRCFILES
    "src/*.cpp"
    "src/LinSysSolver/*.cpp"
)

FOREACH(item ${SRCFILES})
    IF(${item} MATCHES "main.cpp")
        LIST(REMOVE_ITEM SRCFILES ${item})
    ENDIF(${item} MATCHES "main.cpp")
ENDFOREACH(item)

add_library(${PROJECT_NAME}_dev ${SRCFILES})

target_include_directories(${PROJECT_NAME}_dev PUBLIC
    "src"
    "src/LinSysSolver"
)

add_executable(${PROJECT_NAME}_bin "src/main.cpp")
target_link_libraries(${PROJECT_NAME}_bin PUBLIC ${PROJECT_NAME}_dev)
# add path of Eigen header file to include directories
target_include_directories(${PROJECT_NAME}_dev PUBLIC "/usr/local/include/eigen-3.4.0")