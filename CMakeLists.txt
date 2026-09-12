cmake_minimum_required(VERSION 3.23)

project(haar_thinning
    VERSION 0.1.0
    DESCRIPTION "Online Haar-thinning: (1+eps)n i.i.d. samples -> n low-discrepancy points"
    HOMEPAGE_URL "https://github.com/ekene-e/haar-thinning"
    LANGUAGES CXX)

# ---------------------------------------------------------------- options
if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    set(HAAR_IS_TOP_LEVEL ON)
else()
    set(HAAR_IS_TOP_LEVEL OFF)
endif()

option(HAAR_BUILD_CLI      "Build the haar-thin command-line tool"      ${HAAR_IS_TOP_LEVEL})
option(HAAR_BUILD_TESTS    "Build unit tests"                           ${HAAR_IS_TOP_LEVEL})
option(HAAR_BUILD_BENCH    "Build benchmarks"                           ${HAAR_IS_TOP_LEVEL})
option(HAAR_BUILD_EXAMPLES "Build examples"                             ${HAAR_IS_TOP_LEVEL})
option(HAAR_NATIVE         "Compile executables with -march=native"     OFF)
option(HAAR_INSTALL        "Generate install rules"                     ${HAAR_IS_TOP_LEVEL})

if(HAAR_IS_TOP_LEVEL AND NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ---------------------------------------------------------------- library (header-only)
add_library(haar INTERFACE)
add_library(haar::haar ALIAS haar)
target_compile_features(haar INTERFACE cxx_std_23)
target_include_directories(haar INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

# Strict warnings and reproducible floating point for our own executables.
add_library(haar_dev_flags INTERFACE)
if(MSVC)
    target_compile_options(haar_dev_flags INTERFACE /W4 /permissive- /fp:precise)
else()
    target_compile_options(haar_dev_flags INTERFACE
        -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion
        -ffp-contract=off)
    if(HAAR_NATIVE)
        target_compile_options(haar_dev_flags INTERFACE -march=native)
    endif()
endif()

function(haar_add_executable name)
    add_executable(${name} ${ARGN})
    target_link_libraries(${name} PRIVATE haar::haar haar_dev_flags)
    set_target_properties(${name} PROPERTIES CXX_EXTENSIONS OFF)
endfunction()

# ---------------------------------------------------------------- CLI
if(HAAR_BUILD_CLI)
    haar_add_executable(haar-thin src/haar_thin.cpp)
endif()

# ---------------------------------------------------------------- tests
if(HAAR_BUILD_TESTS)
    enable_testing()
    foreach(t geometry table thinner discrepancy replicability)
        haar_add_executable(test_${t} tests/test_${t}.cpp)
        add_test(NAME ${t} COMMAND test_${t})
    endforeach()
endif()

# ---------------------------------------------------------------- bench / examples
if(HAAR_BUILD_BENCH)
    haar_add_executable(haar-bench bench/bench.cpp)
endif()
if(HAAR_BUILD_EXAMPLES)
    haar_add_executable(example_basic examples/basic.cpp)
    haar_add_executable(example_stream examples/stream.cpp)
endif()

# ---------------------------------------------------------------- install / package
if(HAAR_INSTALL)
    include(GNUInstallDirs)
    include(CMakePackageConfigHelpers)

    install(TARGETS haar EXPORT haarTargets)
    install(DIRECTORY include/haar DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    if(HAAR_BUILD_CLI)
        install(TARGETS haar-thin RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    endif()

    install(EXPORT haarTargets
        NAMESPACE haar::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/haar)

    configure_package_config_file(
        cmake/haarConfig.cmake.in
        ${PROJECT_BINARY_DIR}/haarConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/haar)
    write_basic_package_version_file(
        ${PROJECT_BINARY_DIR}/haarConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
        ARCH_INDEPENDENT)
    install(FILES
        ${PROJECT_BINARY_DIR}/haarConfig.cmake
        ${PROJECT_BINARY_DIR}/haarConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/haar)
endif()
