cmake_minimum_required (VERSION 3.8)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project("General")

include_directories(
	"${CMAKE_CURRENT_SOURCE_DIR}/include/general"
)

add_library(
	"General" 

	"src/times.cpp"

	"src/mathfuncs.cpp"
	"src/quaternion.cpp"
	"src/vector.cpp" 
	"src/matrix.cpp"
	"src/polynomial.cpp")

add_executable("TestGeneral" "test/test.cpp")

target_link_libraries("TestGeneral" "General")