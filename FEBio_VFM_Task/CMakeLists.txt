# CMakeList.txt: FEBio_VFM_Task 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
set(FEBIO_SDK_ROOT "$ENV{FEBIO_SDK_ROOT}" CACHE PATH
    "Path to the FEBio SDK directory containing include/ and lib/")

if(NOT FEBIO_SDK_ROOT)
    if(WIN32)
        set(_febio_sdk_candidates
            "C:/Program Files/FEBioStudio/sdk"
            "C:/Program Files/FEBioStudio2/sdk"
            "C:/Program Files/FEBioStudio3/sdk")
    else()
        set(_febio_sdk_candidates
            "$ENV{HOME}/FEBioStudio/sdk"
            "/opt/FEBioStudio/sdk")
    endif()

    foreach(_candidate IN LISTS _febio_sdk_candidates)
        if(EXISTS "${_candidate}/include/FECore/FECore.h")
            set(FEBIO_SDK_ROOT "${_candidate}" CACHE PATH
                "Path to the FEBio SDK directory containing include/ and lib/" FORCE)
            break()
        endif()
    endforeach()
endif()

if(NOT EXISTS "${FEBIO_SDK_ROOT}/include/FECore/FECore.h")
    message(FATAL_ERROR
        "FEBio SDK not found. Set FEBIO_SDK_ROOT to the SDK directory "
        "that contains include/FECore/FECore.h and lib/.")
endif()

message(STATUS "Using FEBio SDK: ${FEBIO_SDK_ROOT}")

set(FEBIO_LINK_LIBS
    fecore
    febioxml
    febiorve
    febioplot
    febiomix
    febiomech
    febiolib
    febiofluid)

if(WIN32)
    set(_febio_release_lib_dir "${FEBIO_SDK_ROOT}/lib/Release")
    set(_febio_debug_lib_dir "${FEBIO_SDK_ROOT}/lib/Debug")
else()
    set(_febio_release_lib_dir "${FEBIO_SDK_ROOT}/../lib")
    set(_febio_debug_lib_dir "${_febio_release_lib_dir}")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Debug mode")
else()
    message(STATUS "Release mode")
endif()

# 为此项目添加OpenMP支持
find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    message(STATUS "OpenMP found")
endif()

set(FEBIO_VFM_TASK_SOURCES
    src/common_FEBio.cpp
    src/common_pybind.cpp
    src/FEBio_refunction.cpp
    src/FEBio_VFM_Task.cpp
    src/optim.cpp
    src/VFMTask.cpp
    src/VFMTask_callback.cpp
    src/VFMTask_debug_output.cpp
    src/VFMTask_laplace.cpp
    src/VFMTask_pybind.cpp
    src/VFMTask_result_capture.cpp)

add_library(FEBio_VFM_Task SHARED ${FEBIO_VFM_TASK_SOURCES})
target_include_directories(FEBio_VFM_Task
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${FEBIO_SDK_ROOT}/include")
target_link_directories(FEBio_VFM_Task PRIVATE
    "$<$<CONFIG:Debug>:${_febio_debug_lib_dir}>"
    "$<$<NOT:$<CONFIG:Debug>>:${_febio_release_lib_dir}>")
target_link_libraries(FEBio_VFM_Task
    PUBLIC
        ${FEBIO_LINK_LIBS}
        OpenMP::OpenMP_CXX)
target_compile_definitions(FEBio_VFM_Task PRIVATE
    "$<$<CONFIG:Debug>:DEBUG>"
    VFM_ENGINE_VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
    VFM_ENGINE_VERSION_MINOR=${PROJECT_VERSION_MINOR}
    VFM_ENGINE_VERSION_PATCH=${PROJECT_VERSION_PATCH})
target_compile_options(FEBio_VFM_Task PRIVATE
    "$<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:-O3>")

# 为此项目引用vcpkg安装的NLopt
find_package(NLopt CONFIG REQUIRED)
target_link_libraries(FEBio_VFM_Task PRIVATE NLopt::nlopt)

find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(FEBio_VFM_Task PRIVATE Eigen3::Eigen)

# 为此项目添加pybind11支持
message(STATUS "Python_EXECUTABLE: ${Python_EXECUTABLE}")

set(Python_FIND_VIRTUALENV FIRST)
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

if(Python_FOUND)
    message(STATUS "Python found")
    message(STATUS "Python version: ${Python_VERSION}")
    message(STATUS "Python include: ${Python_INCLUDE_DIRS}")
    message(STATUS "Python libraries: ${Python_LIBRARIES}")
    message(STATUS "pybind11 include directories: ${pybind11_INCLUDE_DIRS}")
endif()

target_link_libraries(FEBio_VFM_Task PUBLIC Python::Python)
target_link_libraries(FEBio_VFM_Task PUBLIC pybind11::pybind11 pybind11::embed)

# 为此项目添加mio支持
find_package(mio CONFIG REQUIRED)
target_link_libraries(FEBio_VFM_Task PRIVATE mio::mio)

target_compile_features(FEBio_VFM_Task PRIVATE cxx_std_20)

# TODO: 如有需要，请添加测试并安装目标。
