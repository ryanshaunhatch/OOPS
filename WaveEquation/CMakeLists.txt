cmake_minimum_required(VERSION 3.0)
project(WaveEquation)

set(WAVE_INCLUDE_FILES
	include/wave.h
	)
set(WAVE_SOURCE_FILES
	src/wave.cpp
	src/main.cpp
	)

set (SOURCE_FILES ${WAVE_INCLUDE_FILES} ${WAVE_SOURCE_FILES})
add_executable(Wave ${SOURCE_FILES})
target_include_directories(Wave PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(Wave oops ${EXTRA_LIBS})
