#Minimum requirement of CMake version : 3.0.0
cmake_minimum_required(VERSION 3.0.0)

#Project name and version number
project(${libQc} VERSION 0.1.0)

# set include path
include_directories(${libGnssSrc})
include_directories(${libFusionSrc})

add_definitions()
# for all the header files and source_files
file(GLOB_RECURSE header_files ${libQcSrc}/*.h)
file(GLOB_RECURSE source_files ${libQcSrc}/*.cc)

if (CMAKE_SYSTEM_NAME MATCHES "Windows")
    #    # add files to project
    source_group("CMake Files" FILES CMakeLists.txt)
    source_group("Header Files" FILES header_files)
    source_group("Soruce Files" FILES source_files)
else ()

endif ()

add_library(${PROJECT_NAME} ${header_files} ${source_files})
if (CMAKE_SYSTEM_NAME MATCHES "Windows")
    link_directories(${ROOT}/build/Lib/Debug)
    link_directories(${ROOT}/build/Lib/Release)
    link_directories(${ROOT}/build/Lib/RelWithDebInfo)
    link_directories(${ROOT}/build/Lib/MinSizeRel)
else ()
    link_directories(${ROOT}/build/Lib)
endif ()

target_link_libraries(${PROJECT_NAME} ${libGnss} ${libFusion} ${libUnit})
IF(MATLAB_FOUND)
    target_link_libraries(${PROJECT_NAME} ${MATLAB_LIBRARIES})
ENDIF(MATLAB_FOUND)
