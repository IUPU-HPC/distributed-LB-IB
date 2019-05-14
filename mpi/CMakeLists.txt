cmake_minimum_required(VERSION 3.7 FATAL_ERROR)

option(use_itac "Build with itac" OFF)
option(USE_CUDA "Build with CUDA" OFF)

option(USE_STRICT_CHECK "use pedantic flag" OFF)

# include(${CMAKE_SOURCE_DIR}/cmake/common.cmake)

if(CMAKE_BUILD_TYPE MATCHES Debug)
  message("-- ${CMAKE_CURRENT_SOURCE_DIR} > Debug build.")
  set(MPI_CXX_COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
  set(GCC_COVERAGE_COMPILE_FLAGS "-Wall -Wextra -Wconversion -rdynamic")
  if(USE_STRICT_CHECK)
    set(GCC_COVERAGE_COMPILE_FLAGS "-pedantic ${GCC_COVERAGE_COMPILE_FLAGS}")
    set(MPI_CXX_COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
  endif()

  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${GCC_COVERAGE_COMPILE_FLAGS}" )
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${GCC_COVERAGE_COMPILE_FLAGS}" )

elseif(CMAKE_BUILD_TYPE MATCHES Release)
  message("-- ${CMAKE_CURRENT_SOURCE_DIR} > Release build.")
  set(MPI_CXX_COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")

endif()

#output json compile command
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpicxx")

set(CMAKE_CXX_STANDARD 11) #can remove -std=c++11 in CMAKE_CXX_FLAGS

project(distributed-LBIB LANGUAGES CXX)

set(common_library "m" "pthread")

#itac
if(use_itac)
    add_definitions(-DV_T)
    #set(MPI_CXX_COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
    set(VT_ROOT $ENV{VT_ROOT} CACHE PATH "The path vt")
    message("itac enabled! with ${VT_ROOT}")
    include_directories(${VT_ROOT}/include)
    ### see VT_ADD_LIBS
    #set(vt_library ${INTEL_VT}  $ENV{VT_ADD_LIBS})
    set(vt_library "-L$ENV{VT_SLIB_DIR} -lVT $ENV{VT_ADD_LIBS}" CACHE STRING "itac link flags")
    set(common_library ${common_library} ${vt_library})
   #set(common_library ${common_library} "-L/opt/intel/itac/2017.3.030/intel64/lib -lVT -ldwarf -lelf -lvtunwind" )
endif(use_itac)

message("-----current mpicxx flags is ${MPI_CXX_COMPILE_FLAGS}")

# mpi
find_package(MPI REQUIRED)
set(common_library ${common_library} ${MPI_CXX_LIBRARIES})

include(CheckIncludeFiles)

configure_file(
    "${PROJECT_SOURCE_DIR}/config.h.in"
    "${PROJECT_BINARY_DIR}/config.h"
     @ONLY
)

include_directories(${MPI_CXX_INCLUDE_PATH} ${CMAKE_SOURCE_DIR})
include_directories(include ${PROJECT_BINARY_DIR})



add_executable(distributed-lb-ib
  src/main.cpp
  src/cube2thread_and_task.cpp
  src/fiber2thread.cpp
  src/init_gv.cpp
  src/init_df1.cpp
  src/init_eqlbrmdistrfuncDF0.cpp
  src/init_df_inout.cpp
  src/gen_fluid_grid.cpp
  src/gen_fiber_sheet.cpp
  src/do_thread.cpp
  src/compute_bendingforce.cpp
  src/compute_stretchingforce.cpp
  src/compute_elasticforce.cpp
  src/fiber_SpreadForce.cpp
  src/fluid_get_SpreadForce.cpp
  src/compute_eqlbrmdistrfuncDF1.cpp
  src/stream_distrfunc.cpp
  src/bounceback_rigidwalls.cpp
  src/compute_rho_and_u.cpp
  src/fluid_SpreadVelocity.cpp
  src/fiber_get_SpreadVelocity.cpp
  src/replace_old_DF.cpp
  src/copy_inout_to_df2.cpp
  src/print_info.cpp
  src/timer.cpp)

set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
# set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# target_compile_options(distributed-lb-ib ${MPI_CXX_COMPILE_FLAGS} ${CMAKE_CXX_FLAGS})
set_target_properties(distributed-lb-ib PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
# target_link_libraries(distributed-lb-ib m pthread) #remove ${MPI_CXX_LIBRARIES}

target_link_libraries(distributed-lb-ib ${common_library})

#output json compile command
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# if(MPI_COMPILE_FLAGS)
#   set_target_properties(${item} PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
# endif()

# if(MPI_LINK_FLAGS)
#   set_target_properties(${item} PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
# endif()

# add_subdirectory(exe_itac)
