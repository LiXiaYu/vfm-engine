# CMakeList.txt: 顶层 CMake 项目文件，在此处执行全局配置
# 并包含子项目。
#
cmake_minimum_required (VERSION 3.8)


# if complier is clang
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	set(CMAKE_EXPERIMENTAL_CXX_MODULE_CMAKE_API "2182bf5c-ef0d-489a-91da-49dbc3090d2a")
	set(CMAKE_EXPERIMENTAL_CXX_MODULE_DYNDEP TRUE)
	set(CMAKE_EXPERIMENTAL_CXX_MODULE_MAP_FORMAT "clang")
	set(CMAKE_CXX_EXTENSIONS OFF)
	
	# set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES "/usr/local/include/c++/12.0.0")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem /usr/lib/llvm-17/lib/clang/17/include")
	set(CMAKE_EXE_LINKER_FLAGS "-fuse-ld=lld")
	set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("FEBio_VFM_Task")

# 包含子项目。
add_subdirectory ("FEBio_VFM_Task")
