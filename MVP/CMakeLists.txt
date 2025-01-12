cmake_minimum_required(VERSION 3.25)
project(MVP VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 14)

include(CTest)
enable_testing()

include_directories(/usr/local/include)
link_directories(/usr/local/lib)

add_executable(MVP
        BGN.cpp
        BGN.h
        main.cpp
        helper.cpp
        helper.h
        SDPC.cpp
        SDPC.h
        SSEDC.cpp
        SSEDC.h)

target_link_libraries(MVP ntl m gmp pthread relic)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)