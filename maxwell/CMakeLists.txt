cmake_minimum_required(VERSION 3.0)
project(Maxwell)

set(MAXWELL_INCLUDE_FILES
	include/maxwell.h
       )
set(MAXWELL_SOURCE_FILES
        src/maxwell.cpp
        src/main.cpp
       )

set (SOURCE_FILES ${MAXWELL_INCLUDE_FILES} ${MAXWELL_SOURCE_FILES})
add_executable(Maxwell ${SOURCE_FILES})
target_include_directories(Maxwell PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(Maxwell oops ${EXTRA_LIBS})
