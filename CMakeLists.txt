#Minimum requirement of CMake version : 3.0.0
cmake_minimum_required(VERSION 3.0.0)
cmake_policy(SET CMP0074 NEW)

#Project name and version number
project(PPP_AR VERSION 0.1.0)

if(CMAKE_SYSTEM_NAME MATCHES "Linux")
    add_compile_options(-O0)
    add_compile_options(-pthread)
    add_compile_options(-w)
    add_compile_options(-m64)
    add_compile_options(-lz)
    add_compile_options(-lstdc++)
elseif(CMAKE_SYSTEM_NAME MATCHES "Darwin")
    add_compile_options(-g)
    add_compile_options(-O0)
    add_compile_options(-Wall)
    add_compile_options(-pthread)
endif()

if(WIN32)
    MATH(EXPR stack_size "600*1024*1024")
    set(CMAKE_EXE_LINKER_FLAGS "-Wl,--stack,${stack_size}")
endif()

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -pthread")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

#Set the ROOT and subdirectory, you should put the CMakeList.txt in these file directories
set(ROOT ${PROJECT_SOURCE_DIR})
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${ROOT}/build/Bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${ROOT}/build/Lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${ROOT}/build/Lib)

set(CMAKE_DEBUG_POSTFIX "d")
set(CMAKE_RELEASE_POSTFIX "")
set(CMAKE_RELWITHDEBINFO_POSTFIX "rd")
set(CMAKE_MINSIZEREL_POSTFIX "s")

#Choose different compilation configurations according to VS compilation
if(CMAKE_BUILD_TYPE MATCHES "Release")
    set(CMAKE_BUILD_POSTFIX "${CMAKE_RELEASE_POSTFIX}")
elseif(CMAKE_BUILD_TYPE MATCHES "Debug")
    set(CMAKE_BUILD_POSTFIX "${CMAKE_DEBUG_POSTFIX}")
elseif(CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
    set(CMAKE_BUILD_POSTFIX "${CMAKE_RELWITHDEBINFO_POSTFIX}")
elseif(CMAKE_BUILD_TYPE MATCHES "MinSizeRel")
    set(CMAKE_BUILD_POSTFIX "${CMAKE_MINSIZEREL_POSTFIX}")
else()
    set(CMAKE_BUILD_POSTFIX "")
endif()

# find matlab for GUI APP
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake) # add FindMatlab module
find_package(Matlab REQUIRED)
IF(MATLAB_FOUND)
    add_definitions(-DMATLAB)
    include_directories(${MATLAB_INCLUDE_DIR})
    message(STATUS "MATLAB Found, MATLAB MEX will be compiled.")
ELSE(MATLAB_FOUND)
    MESSAGE("MATLAB not found...nothing will be built.")
ENDIF(MATLAB_FOUND)

include_directories(include)
add_definitions(-DENAGLO -DENAGAL -DENACMP -DENAQZS -DNEXOBS=3 -DNFREQ=5 -DTRACE)

# global value
# set lib name
set(libGnss   libGnss  )
set(libQc     libQc)
set(libKf     libKf)
set(libUnit   libUnit  )
# set lib path
set( libGnssSrc    ${PROJECT_SOURCE_DIR}/src/LibGnss  )
set( libQcSrc      ${PROJECT_SOURCE_DIR}/src/LibQc)
set( libKfSrc      ${PROJECT_SOURCE_DIR}/src/LibKf)
set( libUnitSrc    ${PROJECT_SOURCE_DIR}/src/LibUnit  )

# For lib source files
add_subdirectory(${libGnssSrc}     ${ROOT}/build/${libGnss}  )
add_subdirectory(${libQcSrc}       ${ROOT}/build/${libQc})
add_subdirectory(${libKfSrc}       ${ROOT}/build/${libKf})
add_subdirectory(${libUnitSrc}     ${ROOT}/build/${libUnit}  )

#CUI APP
set(App iRTKLIB)
add_subdirectory(${ROOT}/src/app ${ROOT}/build/${App})

# group
SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)
SET_PROPERTY(TARGET ${libGnss}    PROPERTY FOLDER "LIB")
SET_PROPERTY(TARGET ${libQc}      PROPERTY FOLDER "LIB")
SET_PROPERTY(TARGET ${libKf}      PROPERTY FOLDER "LIB")
SET_PROPERTY(TARGET ${libUnit}    PROPERTY FOLDER "LIB")

#Output Messages for debug the Cmake
message(STATUS "operation system is : ${CMAKE_SYSTEM}")
message(STATUS "current platform is : ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMake version is    : ${CMAKE_SYSTEM_VERSION}")
message(STATUS "C compiler is : ${CMAKE_C_COMPILER}")
message(STATUS "C++ compiler is : ${CMAKE_CXX_COMPILER}")
message(STATUS "The program main directory is : ${ROOT}")