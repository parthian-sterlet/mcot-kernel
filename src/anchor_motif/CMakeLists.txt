add_executable(anchor_motif mcot_motif.cpp)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
target_include_directories(anchor_pro PRIVATE ../../include)
