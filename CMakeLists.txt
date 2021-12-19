cmake_minimum_required(VERSION 3.21)
project(fss C)

set(CMAKE_C_STANDARD 11)

include_directories(.)
include_directories(fss32)

add_executable(fss
        fss32/fss32c.c)
