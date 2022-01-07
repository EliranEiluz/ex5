cmake_minimum_required(VERSION 3.21)
project(ex5 C)

set(CMAKE_C_STANDARD 99)

add_executable(ex5 myfunction.c readBMP.c readBMP.h showBMP.c writeBMP.h writeBMP.c)
