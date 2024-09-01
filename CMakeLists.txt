cmake_minimum_required(VERSION 3.14)

project(mumosa
    VERSION 0.1.0
    DESCRIPTION "Multi-modal analyses of single-cell data"
    LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Library
add_library(mumosa INTERFACE)
add_library(libscran::mumosa ALIAS mumosa)

target_include_directories(mumosa INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/mumosa>)
target_compile_features(mumosa INTERFACE cxx_std_17)

# Dependencies
option(MUMOSA_FETCH_EXTERN "Automatically fetch mumosa's external dependencies." ON)
if(MUMOSA_FETCH_EXTERN)
    add_subdirectory(extern)
else()
    find_package(tatami_tatami_stats 3.0.0 CONFIG REQUIRED)
    find_package(knncolle_knncolle 2.0.0 CONFIG REQUIRED)
endif()

target_link_libraries(mumosa INTERFACE knncolle::knncolle tatami::tatami_stats)

# Tests
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    option(MUMOSA_TESTS "Build mumosa's test suite." ON)
else()
    option(MUMOSA_TESTS "Build mumosa's test suite." OFF)
endif()

if(MUMOSA_TESTS)
    include(CTest)
    if(BUILD_TESTING)
        add_subdirectory(tests)
    endif() 
endif()

# Install
install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mumosa)

install(TARGETS mumosa
    EXPORT mumosaTargets)

install(EXPORT mumosaTargets
    FILE libscran_mumosaTargets.cmake
    NAMESPACE libscran::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libscran_mumosa)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/libscran_mumosaConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libscran_mumosa)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/libscran_mumosaConfigVersion.cmake"
    COMPATIBILITY SameMajorVersion)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libscran_mumosaConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/libscran_mumosaConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libscran_mumosa)
