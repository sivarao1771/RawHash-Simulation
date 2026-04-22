cmake_minimum_required(VERSION 3.10)
project(MemristorRawHash)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(USE_MEMRISTOR "Enable memristor hardware emulation" ON)

find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)
find_package(HDF5 REQUIRED)

# Gather all C/C++ sources from src/ directory
file(GLOB RAWHASH_SOURCES src/*.c src/*.cpp)

# Our new memristor hardware emulation sources
set(MEMRISTOR_SOURCES
    src/memristor_crossbar.cpp
    src/memristor_cam.cpp
)

# Create the executable
add_executable(memristor_rawhash
    ${RAWHASH_SOURCES}
    ${MEMRISTOR_SOURCES}
)

target_compile_definitions(memristor_rawhash PRIVATE
    $<$<BOOL:${USE_MEMRISTOR}>:USE_MEMRISTOR>
    NHDF5RH
    NPOD5RH
    NSLOW5RH
)

# Include directories
target_include_directories(memristor_rawhash PRIVATE
    src
    ${HDF5_INCLUDE_DIRS}
)

# Compile definition to enable memristor code
target_compile_definitions(memristor_rawhash PRIVATE
    $<$<BOOL:${USE_MEMRISTOR}>:USE_MEMRISTOR>
)

# macOS specific: OpenMP via libomp
if(APPLE)
    target_compile_options(memristor_rawhash PRIVATE -Xclang -fopenmp)
    target_include_directories(memristor_rawhash PRIVATE /opt/homebrew/opt/libomp/include /usr/local/opt/libomp/include)
    target_link_directories(memristor_rawhash PRIVATE /opt/homebrew/opt/libomp/lib /usr/local/opt/libomp/lib)
endif()

# Link libraries (all plain, single call)
target_link_libraries(memristor_rawhash
    ${ZLIB_LIBRARIES}
    Threads::Threads
    ${HDF5_LIBRARIES}
    $<$<BOOL:${APPLE}>:omp>
)
