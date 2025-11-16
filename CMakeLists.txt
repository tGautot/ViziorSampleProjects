cmake_minimum_required( VERSION 3.21 )

project( ViziorSampleProjects )

set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )


find_package(Vizior REQUIRED)

add_subdirectory(simple_pendulum)
add_subdirectory(circle_path)
