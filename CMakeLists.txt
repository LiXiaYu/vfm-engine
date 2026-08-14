# CMakeList.txt: 顶层 CMake 项目文件，在此处执行全局配置
# 并包含子项目。
#
cmake_minimum_required (VERSION 3.12)


# if complier is clang
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	set(CMAKE_EXPERIMENTAL_CXX_MODULE_CMAKE_API "2182bf5c-ef0d-489a-91da-49dbc3090d2a")
	set(CMAKE_EXPERIMENTAL_CXX_MODULE_DYNDEP TRUE)
	set(CMAKE_EXPERIMENTAL_CXX_MODULE_MAP_FORMAT "clang")
	set(CMAKE_CXX_EXTENSIONS OFF)
	
	# set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES "/usr/local/include/c++/12.0.0")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem /usr/lib/llvm-20/lib/clang/20/include")
	set(CMAKE_EXE_LINKER_FLAGS "-fuse-ld=lld")
	set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

# 强制 MSVC 生成 PDB 文件（Release/Debug 都生成）
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Zi /DEBUG")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /Zi /DEBUG")
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /DEBUG /OPT:REF /OPT:ICF")
set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} /DEBUG /OPT:REF /OPT:ICF")


project ("FEBio_VFM_Task")

# 包含子项目。
add_subdirectory ("FEBio_VFM_Task")
