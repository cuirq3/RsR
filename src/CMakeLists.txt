cmake_minimum_required(VERSION 3.10)

# set the project name
project(Vanilla)

find_package(Boost COMPONENTS filesystem graph REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
find_package(CGAL CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)
# find_package(Open3D REQUIRED)
# find_package(fmt CONFIG REQUIRED)

# add the executable
add_executable(Vanilla Recon.cpp Recon.h)

target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
target_link_libraries(${PROJECT_NAME} general utils)
target_link_libraries(${PROJECT_NAME} CGAL::CGAL)
target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)
# target_link_libraries(${PROJECT_NAME} fmt::fmt-header-only)
# target_link_libraries(${PROJECT_NAME} Open3D::Open3D)
# target_link_libraries(${PROJECT_NAME} ${Open3D_LIBRARIES})
# target_link_libraries(${PROJECT_NAME} TBB::tbb TBB::tbbmalloc TBB::tbbmalloc_proxy)

