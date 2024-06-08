cmake_minimum_required(VERSION 3.28.0)
project(bond-pricing VERSION 0.1.0 LANGUAGES C CXX)
set (CMAKE_CXX_STANDARD 20)

include(CTest)
enable_testing()

# INCLUDE JSON LIB
find_package(nlohmann_json REQUIRED)

# DECLARE EXECUTABLE
add_executable(bond-pricer main.cpp)

# LINK JSON LIB
target_link_libraries(bond-pricer PRIVATE nlohmann_json)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)