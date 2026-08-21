#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

#include "VFMTask.h"
#include "common_FEBio.h"
#include "optim.h"
//定义多个embedded模块

PYBIND11_EMBEDDED_MODULE(VFMTask_pybind, m) {
	pybind11::class_<VFMTask_configure>(m, "VFMTask_configure")
		.def(pybind11::init<>())
		.def_readwrite("bpm", &VFMTask_configure::bpm)
		.def_readwrite("start_index", &VFMTask_configure::start_index)
		.def_readwrite("end_index", &VFMTask_configure::end_index)
		.def_readwrite("run_febio_solve", &VFMTask_configure::run_febio_solve)
		.def_readwrite("reuse_saved_result_buffer", &VFMTask_configure::reuse_saved_result_buffer)
		.def_property("isRead_FEMresult_fromsavefile",
			[](const VFMTask_configure& cfg) {
				return cfg.reuse_saved_result_buffer;
			},
			[](VFMTask_configure& cfg, bool value) {
				if (PyErr_WarnEx(
					PyExc_DeprecationWarning,
					"isRead_FEMresult_fromsavefile is deprecated; use run_febio_solve and reuse_saved_result_buffer",
					1) < 0)
				{
					throw pybind11::error_already_set();
				}
				// Preserve the legacy coupled behavior only for old Python scripts.
				cfg.reuse_saved_result_buffer = value;
				cfg.run_febio_solve = !value;
			})
		.def_readwrite("isReadfromsaveOptimfunc", &VFMTask_configure::isReadfromsaveOptimfunc)
		.def_readwrite("isSetDisplacmentAndPressure", &VFMTask_configure::isSetDisplacmentAndPressure)
		.def_readwrite("allow_failed_solve_postprocessing", &VFMTask_configure::allow_failed_solve_postprocessing)
		.def_readwrite("FEBio_dump_path", &VFMTask_configure::FEBio_dump_path)
		.def_readwrite("Optim_dump_path", &VFMTask_configure::Optim_dump_path)
		.def_readwrite("isReadfromfeblogfile", &VFMTask_configure::isReadfromfeblogfile)
		.def_readwrite("feblog_uxuyuz_path", &VFMTask_configure::feblog_uxuyuz_path)
		.def_readwrite("feblog_sxsyszsxysyzsxz_path", &VFMTask_configure::feblog_sxsyszsxysyzsxz_path)
		.def_readwrite("isLaplaceVFM", &VFMTask_configure::isLaplaceVFM)
		.def_readwrite("InterpolatedDisplacementTimes", &VFMTask_configure::InterpolatedDisplacementTimes)
		.def_readwrite("displacementdataNumber", &VFMTask_configure::displacementdataNumber)
		.def_readwrite("stressdataNumber", &VFMTask_configure::stressdataNumber)
		.def_readwrite("pressurevalueNumber", &VFMTask_configure::pressurevalueNumber)
		.def_readwrite("solution", &VFMTask_configure::solution)
		.def_readwrite("pressure_load", &VFMTask_configure::pressure_load)
		.def_readwrite("constraint_load", &VFMTask_configure::constraint_load)
		.def_readwrite("fixed", &VFMTask_configure::fixed)
		.def_readwrite("fixednode", &VFMTask_configure::fixednode)
		.def_readwrite("initialCoordinate", &VFMTask_configure::initialCoordinate)
		.def_readwrite("vf_u_functions", &VFMTask_configure::vf_u_functions)
		.def_property_readonly("get_timestep", [](VFMTask_configure& cfg) {
			return std::vector<double>(cfg.timestep.begin(), cfg.timestep.end());
				})
		.def("set_timestep", [](VFMTask_configure& cfg, const std::vector<double>& v) {
			cfg.pending_timestep = v; // 仅记录，实际应用由 C++ 端完成
				})
		.def_readwrite("select_solution_element_function", &VFMTask_configure::select_solution_element_function)
		.def_readwrite("setdisplacement_function", &VFMTask_configure::setdisplacement_function)
		.def_readwrite("setstress_function", &VFMTask_configure::setstress_function)
		.def_readwrite("optim_function_iter_max", &VFMTask_configure::optim_function_iter_max)
		.def_readwrite("optim_function", &VFMTask_configure::optim_function)
		.def_readwrite("optim_function_T", &VFMTask_configure::optim_function_T)
		.def_readwrite("optim_function_E__gamma_tau", &VFMTask_configure::optim_function_E__gamma_tau)
		.def_readwrite("optim_function_gamma__E_tau", &VFMTask_configure::optim_function_gamma__E_tau)
		.def_readwrite("optim_function_tau__E_gamma", &VFMTask_configure::optim_function_tau__E_gamma)
		.def_readwrite("optim_function_E_gamma_tau", &VFMTask_configure::optim_function_E_gamma_tau)
		.def_readwrite("optim_function_E_gamma_tau_invF", &VFMTask_configure::optim_function_E_gamma_tau_invF)
		.def_readwrite("optim_function_elastic_E", &VFMTask_configure::optim_function_elastic_E)
		.def_readwrite("optim_function_output_debug_info", &VFMTask_configure::optim_function_output_debug_info)
		.def_readwrite("optim_function_output_iteration_stress", &VFMTask_configure::optim_function_output_iteration_stress)
		.def_readwrite("optim_function_output_iteration_virtual_work", &VFMTask_configure::optim_function_output_iteration_virtual_work)
		.def_readwrite("optim_function_parallel", &VFMTask_configure::optim_function_parallel)
		.def_readwrite("optim_function_num_threads", &VFMTask_configure::optim_function_num_threads)
		.def_readwrite("isSetNewLaplaceVFMs", &VFMTask_configure::isSetNewLaplaceVFMs)
		.def_readwrite("LaplaceVFM_s", &VFMTask_configure::LaplaceVFM_s)
		.def_readwrite("isTalbotLaplaceVFM_s", &VFMTask_configure::isTalbotLaplaceVFM_s);

	// bind function laplace_transform_periodic
	m.def("laplace_transform_periodic", pybind11::overload_cast<const std::vector<double>&,const std::vector<double>&, const std::vector<std::complex<double>>&>(&laplace_transform_periodic), "laplace_transform_periodic");
	m.def("createTalbotPath", &createTalbotPath, "createTalbotPath");
	m.def("talbotInverseLaplaceTransform", &talbotInverseLaplaceTransform, "talbotInverseLaplaceTransform");
	m.def("exit_febio", &exit_febio, "exit_febio");
}
