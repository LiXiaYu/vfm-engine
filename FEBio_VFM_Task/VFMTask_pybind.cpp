#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
#include <pybind11/stl.h>

#include "VFMTask.h"

//定义多个embedded模块

PYBIND11_EMBEDDED_MODULE(VFMTask_pybind, m) {
	pybind11::class_<VFMTask_configure>(m, "VFMTask_configure")
		.def(pybind11::init<>())
		.def_readwrite("bpm", &VFMTask_configure::bpm)
		.def_readwrite("isRead_FEMresult_fromsavefile", &VFMTask_configure::isRead_FEMresult_fromsavefile)
		.def_readwrite("isReadfromsaveOptimfunc", &VFMTask_configure::isReadfromsaveOptimfunc)
		.def_readwrite("displacementdataNumber", &VFMTask_configure::displacementdataNumber)
		.def_readwrite("stressdataNumber", &VFMTask_configure::stressdataNumber)
		.def_readwrite("pressurevalueNumber", &VFMTask_configure::pressurevalueNumber)
		.def_readwrite("solution", &VFMTask_configure::solution)
		.def_readwrite("pressure_load", &VFMTask_configure::pressure_load)
		.def_readwrite("fixed", &VFMTask_configure::fixed);
}