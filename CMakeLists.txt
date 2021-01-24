cmake_minimum_required (VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(Ballistic)

include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/include/general
)

add_library(
	General 

	src/Geometry.cpp
	src/Times.cpp
	src/Mathematics.cpp
	src/Legendre.cpp
	src/Matrix3x3.cpp 
	src/Vector.cpp 
	src/Matrix.cpp 
)

add_executable(
	TestGeneral
	test/Test.cpp)

target_link_libraries(
	TestGeneral
	General
)