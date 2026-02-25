cmake_minimum_required(VERSION 3.20)
project(economic_simulation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

add_executable(economic_simulation ${SOURCES})

target_include_directories(economic_simulation PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  "${CMAKE_CURRENT_SOURCE_DIR}/external"
)

if (MSVC)
  target_compile_options(economic_simulation PRIVATE /W4 /permissive-)
else()
  target_compile_options(economic_simulation PRIVATE -Wall -Wextra -Wpedantic)
endif()