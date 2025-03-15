cmake_minimum_required(VERSION 3.10)
project(Flamme C)

# Set C standard
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Find MPI
find_package(MPI REQUIRED)

# Find HYPRE library and headers
find_path(HYPRE_INCLUDE_DIR NAMES HYPRE.h PATHS /usr/include/hypre /usr/local/include/hypre)
find_library(HYPRE_LIBRARY NAMES HYPRE PATHS /usr/lib /usr/lib64 /usr/local/lib /usr/local/lib64)

# Ensure HYPRE is found
if (NOT HYPRE_INCLUDE_DIR OR NOT HYPRE_LIBRARY)
    message(FATAL_ERROR "HYPRE not found! Please check installation.")
endif()

# Print paths
message(STATUS "HYPRE Include Dir: ${HYPRE_INCLUDE_DIR}")
message(STATUS "HYPRE Library: ${HYPRE_LIBRARY}")

message(STATUS "MPI_INCLUDE_DIRS: ${MPI_C_INCLUDE_DIRS}")
message(STATUS "MPI_LIBRARIES: ${MPI_C_LIBRARIES}")


# Source files
set(SOURCE_FILES
    src/main.c
    src/solver.c
    src/matrix.c
)

# Header files
set(HEADER_FILES
    include/const.h
    include/solver.h
    include/matrix.h
)

# Define executable
add_executable(flamme ${SOURCE_FILES})

# Include directories
target_include_directories(flamme PRIVATE ${PROJECT_SOURCE_DIR}/include ${HYPRE_INCLUDE_DIR} ${MPI_C_INCLUDE_DIRS})

# Link necessary libraries
target_link_libraries(flamme PRIVATE ${HYPRE_LIBRARY} MPI::MPI_C m)
