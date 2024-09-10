cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 20)

# set the project name
project(Vanilla)

find_package(Boost COMPONENTS filesystem graph REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
find_package(CGAL CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
# find_package(Open3D REQUIRED)
# find_package(TBB CONFIG REQUIRED)

# add_subdirectory(src)

# Set the directory containing the .h and .cpp files
set(SOURCE_DIR "src")

# Include all .cpp and .h files from the source directory
file(GLOB SOURCES "${SOURCE_DIR}/*.cpp" "${SOURCE_DIR}/*.h")

# Specify the entry point file
set(ENTRY_POINT "${SOURCE_DIR}/Recon.cpp")

# add the executable
add_executable(Vanilla ${ENTRY_POINT} ${SOURCES})

target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES} Boost::filesystem)
target_link_libraries(${PROJECT_NAME} CGAL::CGAL)
target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)

