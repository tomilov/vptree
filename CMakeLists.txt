cmake_minimum_required(VERSION 3.27.1 FATAL_ERROR)

project("vptree" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS YES)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_executable("${PROJECT_NAME}")

target_sources(
    "${PROJECT_NAME}"
    PRIVATE
        "main.cpp")

target_compile_options(
    "${PROJECT_NAME}"
    PRIVATE
        "-march=native"
        #"-fno-strict-aliasing"
        "-W"
        "-Wno-gnu-case-range"
        "-Wall"
        "-Wextra"
        "-Wno-c++98-compat"
        "-Wno-c++98-compat-pedantic"
        "-Wconversion"
        "-Werror")

#find_package(OpenMP REQUIRED)
find_package(fmt REQUIRED)
target_link_libraries(
    "${PROJECT_NAME}"
    PRIVATE
        #OpenMP::OpenMP_CXX
        fmt::fmt)

if(FALSE)
    target_link_options("${PROJECT_NAME}" PRIVATE -fsanitize=undefined)
    target_compile_options("${PROJECT_NAME}" PRIVATE -fsanitize=undefined)
endif()
