cmake_minimum_required (VERSION 3.8)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(Ballistic)

include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/include/general
)

add_library(
	General 

	
	src/Times.cpp
	src/Mathematics.cpp
	  
	src/Quaternion.cpp
	src/Vector.cpp 
	src/Matrix.cpp
	src/Polynomial.cpp)

add_executable(
	TestGeneral
	test/Test.cpp)

target_link_libraries(
	TestGeneral
	General
)