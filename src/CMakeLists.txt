cmake_minimum_required(VERSION 2.8)
project(LIGHT)

set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -march=native -Wall")
find_package(Threads)

set(SOURCE_FILES Config.h
        FSEMain.cpp
        Graph.h
        Graph.cpp
        CommandParser.h
        FSECommand.h
        ComputeSetIntersection.h
        ComputeSetIntersection.cpp
        FSE.h
        FSE.cpp
        Utility.h
        Utility.cpp
        TaskSlot.h
        TaskSlot.cpp
        TaskPool.h
        TaskPool.cpp
        )
add_executable(LIGHT ${SOURCE_FILES})
target_link_libraries(LIGHT ${CMAKE_THREAD_LIBS_INIT})

set(CONVERTER_SOURCE_FILES
        Config.h
        ConverterCommand.h
        ConverterMain.cpp
        Graph.h
        Graph.cpp
        CommandParser.h
        Utility.h
        Utility.cpp
        )
add_executable(Converter ${CONVERTER_SOURCE_FILES})