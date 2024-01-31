# Force CMake to use GCC-13
# set(CMAKE_CXX_COMPILER "/usr/bin/gcc-13")

cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

project(linear_algebra_ma VERSION 0.1 LANGUAGES CXX)


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# check used compiler
set(gcc_like_cxx "$<COMPILE_LANG_AND_ID:CXX,ARMClang,AppleClang,Clang,GNU,LCC>")
set(msvc_cxx "$<COMPILE_LANG_AND_ID:CXX,MSVC>")


#-- complex
add_library(complex 
    "src/complex/complex.cpp"
)
# include directories
target_include_directories(complex PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
# turn off warning from #pragma region, #pragma endregion and set default ros2 warnings
target_compile_options(complex PRIVATE 
    "$<${gcc_like_cxx}:-Wall;-Wextra;-Wshadow;-Wpedantic;-Wunused;-Wno-unknown-pragmas>"
)


#-- matrices
# add library
add_library(matrices INTERFACE)
# set c++20 standard for the library
target_compile_features(matrices INTERFACE cxx_std_20)
# include directories
target_include_directories(matrices INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>

    $<INSTALL_INTERFACE:include>
)
# source files
target_sources(matrices INTERFACE 
    ${CMAKE_CURRENT_SOURCE_DIR}/include/linear_algebra_ma/matrices.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/linear_algebra_ma/matrices/matrices.tpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/linear_algebra_ma/matrices/operators.tpp
)
# turn off warning from #pragma region, #pragma endregion and set default ros2 warnings
target_compile_options(matrices INTERFACE 
    "$<${gcc_like_cxx}:-Wall;-Wextra;-Wshadow;-Wpedantic;-Wunused;-Wno-unknown-pragmas>"
)


# #-- geometries
# # add library
# add_library(geometries
#     "src/linear_algebra_ma/geometries/point.cpp"
#     "src/linear_algebra_ma/geometries/segment.cpp"
# )
# # include directories
# target_include_directories(geometries PUBLIC
#     $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#     $<INSTALL_INTERFACE:include>
# )
# # turn off warning from #pragma region, #pragma endregion and set default ros2 warnings
# target_compile_options(geometries PRIVATE 
#     "$<${gcc_like_cxx}:-Wall;-Wextra;-Wshadow;-Wpedantic;-Wunused;-Wno-unknown-pragmas>"
# )
# target_link_libraries(geometries matrices)


# #-- all
# # add library
# add_library( ${PROJECT_NAME}
#     "src/linear_algebra_ma/matrices/matrices.cpp"
#     "src/linear_algebra_ma/matrices/operators.cpp"
#     "src/linear_algebra_ma/geometries/point.cpp"
#     "src/linear_algebra_ma/geometries/segment.cpp"
# )

# # include directories
# target_include_directories(${PROJECT_NAME} PUBLIC
#     $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#     $<INSTALL_INTERFACE:include>
# )

# # turn off warning from #pragma region, #pragma endregion and set default ros2 warnings
# target_compile_options(${PROJECT_NAME} PRIVATE 
#     "$<${gcc_like_cxx}:-Wall;-Wextra;-Wshadow;-Wpedantic;-Wunused;-Wno-unknown-pragmas>"
# )

# following actions are performed only if this is the root project (not included by other modules)
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)

################################################################################
# Unit tests
################################################################################

add_executable(my_test "tests/test.cpp")
target_link_libraries(my_test PRIVATE
    matrices
)
target_include_directories(my_test PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

#-- gtest
include(FetchContent)
FetchContent_Declare(
  googletest
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
if(WIN32)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
endif()
FetchContent_MakeAvailable(googletest)
add_library(GTest::GTest INTERFACE IMPORTED)
# gtest_main automatically linked together with GTest::GTest
target_link_libraries(GTest::GTest INTERFACE gtest_main)

# enables ctesting
enable_testing()

add_executable(${PROJECT_NAME}_unit_tests
    "tests/test_complex.cpp"
    "tests/test_matrices.cpp"
    # "tests/test_geometries.cpp"
)
target_link_libraries(${PROJECT_NAME}_unit_tests
    PRIVATE
    GTest::GTest
    gtest_main
    complex
    matrices
    # geometries
)

# add the gtest target as a ctest. CTests will be run in block
add_test(${PROJECT_NAME}_gtests ${PROJECT_NAME}_unit_tests)

# make google test register the tests. CTests will run the single defined tests 
include(GoogleTest)
gtest_discover_tests( ${PROJECT_NAME}_unit_tests )


################################################################################
# Create package
################################################################################

# https://www.foonathan.net/2022/06/cmake-fetchcontent/

set(package_files include/ src/ CMakeLists.txt LICENSE)
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-src.zip
    COMMAND ${CMAKE_COMMAND} -E tar c ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-src.zip --format=zip -- ${package_files}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS ${package_files})
add_custom_target(${PROJECT_NAME}_package DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-src.zip)


endif() # CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR