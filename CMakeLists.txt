cmake_minimum_required(VERSION 3.16)
project(FFT_Spectrum_Analyzer VERSION 1.0 LANGUAGES CXX)

# ── C++17 Standard ────────────────────────────────────────────
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)   # Use -std=c++17, not -std=gnu++17

# ── Compiler Warnings ─────────────────────────────────────────
if(MSVC)
    add_compile_options(/W4)
else()
    # GCC / Clang warnings
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# ── Output Directories ─────────────────────────────────────────
# Place executables in build/ (not in nested subdirectories)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# ── Include Directories ────────────────────────────────────────
include_directories(${CMAKE_SOURCE_DIR}/include)

# ── Static Library: FFT + SignalAnalyzer (shared by main + tests) ──
add_library(fft_lib STATIC
    src/fft.cpp
    src/signal_analyzer.cpp
)
target_include_directories(fft_lib PUBLIC ${CMAKE_SOURCE_DIR}/include)

# ── Main Executable: spectrum_analyzer ────────────────────────
add_executable(spectrum_analyzer
    src/main.cpp
)
target_link_libraries(spectrum_analyzer PRIVATE fft_lib)

# ── Test Executable: test_fft ─────────────────────────────────
add_executable(test_fft
    tests/test_fft.cpp
)
target_link_libraries(test_fft PRIVATE fft_lib)

# ── CTest Integration ─────────────────────────────────────────
enable_testing()

add_test(
    NAME FFT_TestSuite
    COMMAND test_fft
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

# ── Installation (optional) ───────────────────────────────────
install(TARGETS spectrum_analyzer DESTINATION bin)

# ── Print build info ─────────────────────────────────────────
message(STATUS "--------------------------------------------")
message(STATUS "Project: ${PROJECT_NAME} v${PROJECT_VERSION}")
message(STATUS "C++ Standard: C++${CMAKE_CXX_STANDARD}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "Build dir: ${CMAKE_BINARY_DIR}")
message(STATUS "--------------------------------------------")
