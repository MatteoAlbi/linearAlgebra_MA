cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

project(linear_algebra_ma VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# check used compiler
set(gcc_like_cxx "$<COMPILE_LANG_AND_ID:CXX,ARMClang,AppleClang,Clang,GNU,LCC>")
set(msvc_cxx "$<COMPILE_LANG_AND_ID:CXX,MSVC>")


#-- matrices
# add library
add_library(matrices 
    "src/linear_algebra_ma/matrices/matrices.cpp"
    "src/linear_algebra_ma/matrices/operators.cpp"
)

# include directories
target_include_directories(matrices PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# turn off warning from #pragma region, #pragma endregion and set default ros2 warnings
target_compile_options(matrices PRIVATE 
    "$<${gcc_like_cxx}:-Wall;-Wextra;-Wshadow;-Wpedantic;-Wunused;-Wno-unknown-pragmas>"
)


#-- geometries
# add library
add_library(geometries
    "src/linear_algebra_ma/geometries/point.cpp"
    "src/linear_algebra_ma/geometries/line.cpp"
)

# include directories
target_include_directories(geometries PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# turn off warning from #pragma region, #pragma endregion and set default ros2 warnings
target_compile_options(geometries PRIVATE 
    "$<${gcc_like_cxx}:-Wall;-Wextra;-Wshadow;-Wpedantic;-Wunused;-Wno-unknown-pragmas>"
)

target_link_libraries(geometries matrices)


################################################################################
# Unit tests
################################################################################

#-- gtest
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
  DOWNLOAD_EXTRACT_TIMESTAMP true
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
if(WIN32)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
endif()
FetchContent_MakeAvailable(googletest)

enable_testing()

# test linear algebra
add_executable(matrices_unit_test
    "tests/test_matrices.cpp"
)
target_link_libraries(matrices_unit_test
    GTest::gtest_main
    matrices
)

#test geometries
add_executable(geometries_unit_test
    "tests/test_geometries.cpp"
)
target_link_libraries(geometries_unit_test
    GTest::gtest_main
    geometries
)

include(GoogleTest)
gtest_discover_tests(
    matrices_unit_test
    geometries_unit_test
)