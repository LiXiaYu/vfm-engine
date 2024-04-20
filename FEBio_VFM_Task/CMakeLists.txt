# CMakeList.txt: FEBio_VFM_Task 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
if (WIN32)
	set(FEBioSDKPath "C:/Program Files/FEBioStudio2/sdk")
	set(LINK_LIBS
		fecore.lib
		febioxml.lib
		febiorve.lib
		febioplot.lib
		febiomix.lib
		febiomech.lib
		febiolib.lib
		febiofluid.lib
	)
else()
	set(FEBioSDKPath "/home/wangxiaofei/FEBioStudio/sdk")
	set(LINK_LIBS
		libfecore.so
		libfebiorve.so
		libfebiomix.so
		libfebiomech.so
		libfebiolib.so
		libfebiofluid.so
	)
endif()



include_directories ("${FEBioSDKPath}/include")
if (WIN32)
	if(CMAKE_BUILD_TYPE STREQUAL "Release")
		link_directories ("${FEBioSDKPath}/lib/Release")
	else()
		link_directories ("${FEBioSDKPath}/lib/Debug")
	endif()
else()
	# linux上只有release版本
	link_directories ("${FEBioSDKPath}/../lib")
endif()

# 为此项目添加OpenMP支持
find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
	message("OPENMP FOUND")
	
	# if debug not add -O3
	if(CMAKE_BUILD_TYPE STREQUAL "Release")
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -O3")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -O3")
		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} -O3")
	else()
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")		
	endif()

endif()

# 将源代码添加到此项目的可执行文件。
add_library (FEBio_VFM_Task SHARED "FEBio_VFM_Task.cpp" "FEBio_VFM_Task.h" "VFMTask.cpp" "VFMTask.h" "VectorNum.h" "optim.h" "optim.cpp" "common_FEBio.h" "common_FEBio.cpp" "unilang_import.h" "FEBio_refunction.cpp" "FEBio_refunction.h" )
target_link_libraries (FEBio_VFM_Task PUBLIC ${LINK_LIBS})

# 为此项目引用vcpkg安装的NLpot
find_package(NLopt CONFIG REQUIRED)
target_link_libraries (FEBio_VFM_Task PRIVATE NLopt::nlopt)

# 为此项目添加pybind11支持

# set(Python3_FIND_VIRTUALENV FIRST)
find_package(Python COMPONENTS Interpreter Development REQUIRED)

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

find_package(pybind11 CONFIG REQUIRED)

if (Python_FOUND)
	message("Python Found")
	message("Python Version: ${Python_VERSION}")
	message("Python Include: ${Python_INCLUDE_DIRS}")
	message("Python Libraries: ${Python_LIBRARIES}")
	# need libpython and libpython _d

	message("pybind11 INCLUDE_DIRECTORIES: ${pybind11_INCLUDE_DIRS}")
else()
	message("Python Not Found")
endif()

# target_link_libraries (FEBio_VFM_Task PRIVATE ${pybind11_INCLUDE_DIRS})
# target_include_directories (FEBio_VFM_Task PRIVATE "${PYTHON_INCLUDE_DIRS}")
# target_link_directories (FEBio_VFM_Task PRIVATE "${PYTHON_LIBRARY_DIRS}")
# target_link_libraries (FEBio_VFM_Task PRIVATE "${PYTHON_LIBRARIES}")

target_link_libraries (FEBio_VFM_Task PRIVATE Python::Python)
target_link_libraries (FEBio_VFM_Task PRIVATE pybind11::pybind11 pybind11::embed)

# 为此项目引用目录下的dll: unilang.dll
if (WIN32)
	set (UnilangPath "C:/Users/liqi/programming/mona/unilang/unilang/build")

	link_directories("${UnilangPath}")
	target_link_libraries (FEBio_VFM_Task PRIVATE "${UnilangPath}/unilang.lib")
else()
	# set (CMAKE_CXX_FLAGS "-L/lib/x86_64-linux-gnu -lc -lm -ldl -lrt")
	set (UnilangPath "/home/wangxiaofei/LiXiaYu/programming/unilang/unilang_dll/build")

	link_directories("${UnilangPath}")
	target_link_libraries (FEBio_VFM_Task PRIVATE "${UnilangPath}/unilang.dll")

	# set_target_properties(FEBio_VFM_Task PROPERTIES LINK_FLAGS "-Wl,-rpath,/path/to/your/libc.so.6")
endif()

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET FEBio_VFM_Task PROPERTY CXX_STANDARD 20)
endif()

# TODO: 如有需要，请添加测试并安装目标。
