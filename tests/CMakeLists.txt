include(FetchContent)
FetchContent_Declare(
  scran_tests
  GIT_REPOSITORY https://github.com/libscran/scran_tests
  GIT_TAG master
)
FetchContent_MakeAvailable(scran_tests)

add_executable(
    libtest 
    src/mumosa.cpp
)

target_link_libraries(
    libtest
    gtest_main
    scran_tests
    mumosa
)

target_compile_options(libtest PRIVATE -Wall -Werror -Wpedantic -Wextra)

option(CODE_COVERAGE "Enable coverage testing" OFF)
if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(libtest PRIVATE -O0 -g --coverage)
    target_link_options(libtest PRIVATE --coverage)
endif()

include(GoogleTest)
gtest_discover_tests(libtest)
