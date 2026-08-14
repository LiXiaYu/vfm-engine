#include "VFMTask.h"

#include <vector>
#include <set>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <filesystem>
#include <algorithm>
#include <cstring>
#include <stdexcept>

#include <signal.h>

// if in linux
#ifdef __linux__
#include <dlfcn.h>
#endif

#include <omp.h>

#include <nlopt.hpp>

#include "optim.h"
#include "VectorNum.h"

#include "VFMTask_callback.h"
#include "FEBio_refunction.h"

VFMTask::VFMTask(FEModel* pfem) :FECoreTask(pfem)
{
}

VFMTask::~VFMTask()
{

}


int push_back_vector__int(std::uintptr_t p_v, int a1)
{
	((::std::vector<int>*)p_v)->push_back(a1);
	return 0;
}
int push_back_vector__tuple__string_string(std::uintptr_t p_v, char* s1, char* s2)
{
	::std::string ss1(s1);
	::std::string ss2(s2);


	((::std::vector<::std::tuple<::std::string, ::std::string>>*)p_v)->push_back(::std::make_tuple(ss1, ss2));
	return 0;
}

// // openmp 4.5
//#pragma omp declare reduction \
//	(vec_double_plus : ::std::vector<double> : omp_out = ::std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
//	initializer(omp_priv = decltype (omp_orig) (omp_orig.size()))

std::vector<std::vector<double>> createGrid(const std::vector<double>& xl, const std::vector<double>& xu, const std::vector<double>& xd, int now_index, ::std::vector<double>* now_xs) {
	if (now_index <= -1)
	{
		now_index = xl.size() - 1;
	}
	if (now_xs == nullptr)
	{
		now_xs = new ::std::vector<double>(xl.size());
	}

	std::vector<std::vector<double>> grid;
	if (xl.size() == 0 || xu.size() == 0 || xd.size() == 0) {
		return grid; // 如果输入数组为空，则返回空数组
	}

	if (now_index == 0) {
		std::vector<std::vector<double>> grid;

		for (double x = xl[now_index]; x <= xu[now_index]; x += xd[now_index])
		{
			(*now_xs)[now_index] = x;
			grid.push_back(*now_xs);
		}

		return grid;
	}
	else
	{

		for (double x = xl[now_index]; x <= xu[now_index]; x += xd[now_index])
		{
			(*now_xs)[now_index] = x;
			std::vector<std::vector<double>> subGrid = createGrid(xl, xu, xd, now_index - 1, now_xs);
			grid.insert(grid.end(), subGrid.begin(), subGrid.end());
			if (::std::abs(xd[now_index]) <= 1e-7)
			{
				break;
			}
		}

		if (now_index == xl.size() - 1)
		{
			delete now_xs;
		}

		return grid;
	}





	int n = xl.size() - 1; // 当前维度的网格点数
	double dx = xd[n - 1]; // 当前维度的间隔

	for (double x = xl[n - 1]; x <= xu[n - 1]; x += dx) {
		std::vector<double> point;
		for (int i = 0; i < n; i++) {
			point.push_back(x); // 当前维度的值
			x += dx * xu[i] / xl[i]; // 计算下一个维度的值
		}
		grid.push_back(point);
	}

	// 递归处理下一个维度的网格
	std::vector<std::vector<double>> subGrid = createGrid(std::vector<double>(xl.begin(), xl.end() - 2), std::vector<double>(xu.begin(), xu.end() - 2), std::vector<double>(xd.begin(), xd.end() - 2));
	for (const auto& row : subGrid) {
		grid.push_back(row);
	}

	return grid;
}


void sigHandler(int num) {
	std::cout << "sigHandler done\n";
	exit(1);
}


bool VFMTask::Init(const char* szfile)
{
	signal(SIGINT, sigHandler);

	// Read file
	this->szfile = szfile;

	// Parse
		// could import a script language>?
		// normal text file> 

		// python file>
	::std::string sifilename = szfile;
	::std::string suffix = sifilename.substr(sifilename.find_last_of('.') + 1);

	// szfile remove the suffix
	::std::string szfile_no_suffix = sifilename.substr(0, sifilename.find_last_of('.'));

	::std::string python_module_name = szfile_no_suffix;


#ifdef __linux__
	// if debug
//#ifdef __DEBUG__

//#else
//	auto pylibdl = dlopen("/home/wangxiaofei/anaconda3/envs/base-pydebug/lib/libpython3.12.so", RTLD_LAZY | RTLD_GLOBAL);
//#endif
#endif

	// Start python interpreter in Class Constructor (this->guard)
	namespace py = pybind11;
	using namespace py::literals;

	//调用Python函数
	py::print("Hello,World");

	//---------------------------------------------------------------------
	py::print("Python path:");
	//加载Python模块，获取模块变量(py::object)
	auto sys = py::module_::import("sys");
	py::print(sys.attr("path"));

#ifdef __linux__
	auto paths = sys.attr("path").cast<::std::vector<::std::string>>();
	// 根据python路径，找到libpython的路径
	auto libpython2 = paths[2];
	auto pythonname= libpython2.substr(libpython2.find_last_of('/') + 1);

#ifdef DEBUG
	auto tagdebug = "d";
#else
	auto tagdebug = "";
#endif

	auto libpython_path = libpython2.substr(0, libpython2.find_last_of('/')) + "/lib"+ pythonname+ tagdebug +".so";

	std::cout<< "libpython_path: " << libpython_path << std::endl;

	auto pylibdl = dlopen(libpython_path.c_str(), RTLD_LAZY | RTLD_GLOBAL);
#endif
	//---------------------------------------------------------------------

	//加载pybind11嵌入模块，运行命令
	std::cout << "Load python file: " << sifilename << std::endl;
	// 加载外部文件
	try
	{
		// 使用 importlib.util 从文件路径导入 Python 脚本
		py::object importlib = py::module_::import("importlib.util");

		// 创建模块加载 spec
		py::object spec = importlib.attr("spec_from_file_location")("vfm_script", sifilename);

		// 创建模块对象
		py::object module = importlib.attr("module_from_spec")(spec);

		// 执行模块加载
		spec.attr("loader").attr("exec_module")(module);

		// 保存到成员变量
		this->pyfile_module = module;

	}
	catch (pybind11::error_already_set e)
	{
		::std::cout << e.what();
		return false;
	}

	// InitVFMTask
	try
	{
		auto result = pyfile_module.attr("InitVFMTask")(this->configure);
		auto vfm_configure_from_py = result.cast<VFMTask_configure>();
		// 从python文件中获取配置
		this->configure = vfm_configure_from_py;
	}
	catch (pybind11::error_already_set e)
	{
		::std::cout << e.what();
		return false;
	}


	// number of cycles considered at the end
	int numCycle = 2;

	// 获取当前程序运行路径
	::std::string runtimepath = std::filesystem::current_path().string();


	std::filesystem::path filepath(szfile);

	// 获取完整路径
	std::string fullpath = filepath.string();

	// 获取文件名（不含路径）
	std::string filename = filepath.filename().string();

	// 获取不带扩展名的文件名
	std::string stemname = filepath.stem().string();

	// create a save file
	::std::filesystem::path outSavefilePath = std::filesystem::current_path() / (stemname + "_BMP" + ::std::to_string(static_cast<int>(this->configure.bpm)) + "_save.txt");
	this->outSavefile = outSavefilePath.string();
	if (this->configure.isRead_FEMresult_fromsavefile == false)
	{
		// create a save file
		::std::ofstream outFile;
		outFile.open(this->outSavefile, ::std::ios::out);
		if (!outFile.is_open())
		{
			feLogError(("Open file " + this->outSavefile+"\"").c_str());
			feLogError("save.txt");
			feLogError("\" failed.\n");
			return false;
		}

	}

	// create a log file
	::std::filesystem::path outlogfilePath = std::filesystem::current_path() / (stemname + "_BMP" + ::std::to_string(static_cast<int>(this->configure.bpm)) + "_log.txt");
	this->outlogfile = outlogfilePath.string();
	::std::ofstream outFile;
	outFile.open(this->outlogfile, ::std::ios::out);
	if (!outFile.is_open())
	{
		feLogError("Open file \"");
		feLogError("log.txt");
		feLogError("\" failed.\n");
		return false;
	}

	// set dumpfile name
	if (this->configure.FEBio_dump_path.empty() == false)
	{
		this->dumpfile = this->configure.FEBio_dump_path;
	}
	else
	{
		::std::filesystem::path dumpfilePath = std::filesystem::current_path() / (stemname + "_BMP" + ::std::to_string(static_cast<int>(this->configure.bpm)) + "_dump.febdump");
		this->dumpfile = dumpfilePath.string();
	}
	

	// write log
	outFile << endl;
	outFile << "BMP:" << this->configure.bpm << endl;
	outFile << "numCycle:" << numCycle << endl;
	outFile.close();

	FEModel& fem = *GetFEModel();
	// from fem get the model's initial displacement

	fem.AddCallback(read_inited_information, CB_INIT, (void*)this);

	fem.AddCallback(read_stepsolved_information, CB_TIMESTEP_SOLVED, (void*)this);

	//fem.AddCallback(read_solved_information, CB_SOLVED, (void*)this);


	const bool ifsetlc = false;
	if (ifsetlc)
	{
		//find loadcontroller name="BMP120" standard 120
		FELoadController* lc = nullptr;
		for (int j = 0; j < fem.LoadControllers(); j++)
		{
			auto* temp_lc = fem.GetLoadController(j);
			::std::string lcname = temp_lc->GetName();
			if (lcname == "BPM" + ::std::to_string(static_cast<int>(this->configure.bpm)))
			{
				// set the load controller
				lc = temp_lc;
				break;
			}
		}

		// find surface_load name="PressureLoad1" and set new lc by bmp
		for (int j = 0; j < this->configure.pressure_load.size(); j++)
		{
			for (int i = 0; i < fem.ModelLoads(); i++)
			{
				FEModelLoad& load0 = *(fem.ModelLoad(i));
				::std::string loadclassname = load0.GetFactoryClass()->GetClassName();
				::std::string loadname = load0.GetName();

				if (::std::get<0>(this->configure.pressure_load[j]) == loadclassname && ::std::get<1>(this->configure.pressure_load[j]) == loadname)
				{
					if (loadclassname == "FEPressureLoad")
					{
						FEPressureLoad& pressureLoad0 = dynamic_cast<FEPressureLoad&>(load0);

						fem.AttachLoadController(pressureLoad0.GetParameter("pressure"), lc);
					}

					break;
				}
			}
		}
	}
	//// auto adjust FEM total time
	//auto& step = *(fem.GetStep(fem.Steps() - 1));
	//if (static_cast<int>(this->bmp) < 120)
	//{
	//	step.m_ntime *= (120 / static_cast<int>(this->bmp));

	//	step.m_dt0 = step.m_dt0 * (120 / static_cast<int>(this->bmp)); // add time step
	//}

	// continue
	return fem.Init();
}



bool VFMTask::Run()
{
	feLog("start run VFM task.\n");

	// displacement at the initial timestep of the last two cycles
	FEModel& fem = *GetFEModel();
	// from fem get the model's initial displacement

	//fem.AddCallback(read_inited_information, CB_INIT, (void*)this);

	//fem.AddCallback(read_solved_information, CB_SOLVED, (void*)this);


	bool femsolveresult = true;
	if (this->configure.isRead_FEMresult_fromsavefile == false)
	{
		bool femsolveresult = fem.Solve();
	}

	read_solved_information(&fem, CB_SOLVED, (void*)this); // no need callback

	return femsolveresult;
}

//void VFMTask_configure::add_vf_u_function(pybind11::object callback)
//{
//	// callback is a python function
//	std::function<::std::vector<double>(const ::std::vector<double>&, int)> func = [callback](const ::std::vector<double>& x, int index) {
//		//pybind11::gil_scoped_acquire acquire;
//		pybind11::object result = callback(x, index);
//		return result.cast<::std::vector<double>>();
//		};
//	this->vf_u_functions.push_back(func);
//}

void VFMTask::reinit_mmap_and_spans()
{
	if (!configure.pending_timestep.has_value()) return;

	const std::vector<double>& new_ts = *configure.pending_timestep;
	const size_t new_steps = new_ts.size();
	const size_t old_steps = this->total_steps;

	if (new_steps == old_steps) { configure.pending_timestep.reset(); return; }

	// Compute sizes for old/new layout (bytes)
	auto time_bytes = [](size_t steps) { return steps * sizeof(double); };
	auto disp_bytes = [&](size_t steps) { return steps * ndisplacement * 3 * sizeof(double); };
	auto stress_bytes = [&](size_t steps) { return steps * nstress * 6 * sizeof(double); };
	auto nforce_bytes = [&](size_t steps) { return steps * nnodalforce * 3 * sizeof(double); };
	auto cpress_bytes = [&](size_t steps) { return steps * nconstraint * sizeof(double); };
	auto cact_bytes = [&](size_t steps) { return steps * nconstraint * sizeof(uint8_t); };

	// Old layout offsets
	const size_t off_time_old = 0;
	const size_t off_disp_old = off_time_old + time_bytes(old_steps);
	const size_t off_stress_old = off_disp_old + disp_bytes(old_steps);
	const size_t off_nforce_old = off_stress_old + stress_bytes(old_steps);
	const size_t off_cpress_old = off_nforce_old + nforce_bytes(old_steps);
	const size_t off_cact_old = off_cpress_old + cpress_bytes(old_steps);
	const size_t total_old = off_cact_old + cact_bytes(old_steps);

	// New layout offsets
	const size_t off_time_new = 0;
	const size_t off_disp_new = off_time_new + time_bytes(new_steps);
	const size_t off_stress_new = off_disp_new + disp_bytes(new_steps);
	const size_t off_nforce_new = off_stress_new + stress_bytes(new_steps);
	const size_t off_cpress_new = off_nforce_new + nforce_bytes(new_steps);
	const size_t off_cact_new = off_cpress_new + cpress_bytes(new_steps);
	const size_t total_new = off_cact_new + cact_bytes(new_steps);

	// Helper to rebuild spans from mapped base
	auto rebuild_spans = [&](uint8_t* base, size_t steps) {
		uint8_t* p = base;
		configure.timestep = span<double>{ reinterpret_cast<double*>(p), steps };
		p += time_bytes(steps);
		timedisplacement = span3d<double>{ reinterpret_cast<double*>(p), steps, ndisplacement, 3 };
		p += disp_bytes(steps);
		timestress = span3d<double>{ reinterpret_cast<double*>(p), steps, nstress, 6 };
		p += stress_bytes(steps);
		timenodalforce = span3d<double>{ reinterpret_cast<double*>(p), steps, nnodalforce, 3 };
		p += nforce_bytes(steps);
		timeconstraintpressure = span2d<double>{ reinterpret_cast<double*>(p), steps, nconstraint };
		p += cpress_bytes(steps);
		timeconstraintactivate = span2d<uint8_t>{ reinterpret_cast<uint8_t*>(p), steps, nconstraint };
		};

	if (new_steps > old_steps) {
		// GROW: resize up, map to new size, move blocks from tail to head
		this->mmap.unmap();
		std::filesystem::resize_file(this->dumpfile, total_new);
		std::error_code ec;
		this->mmap.map(this->dumpfile, 0, mio::map_entire_file, ec);
		if (ec) {
			throw std::runtime_error("mmap failed: " + ec.message());
		}

		uint8_t* base = (uint8_t*)this->mmap.data();

		// Move blocks forward (use memmove for overlap safety), from last to first
		std::memmove(base + off_cact_new, base + off_cact_old, cact_bytes(old_steps));
		std::memmove(base + off_cpress_new, base + off_cpress_old, cpress_bytes(old_steps));
		std::memmove(base + off_nforce_new, base + off_nforce_old, nforce_bytes(old_steps));
		std::memmove(base + off_stress_new, base + off_stress_old, stress_bytes(old_steps));
		std::memmove(base + off_disp_new, base + off_disp_old, disp_bytes(old_steps));
		// timestep will be overwritten by pending, no need to move

		// Rebuild spans on new layout
		rebuild_spans(base, new_steps);

		// Write pending timestep
		std::copy(new_ts.begin(), new_ts.end(), configure.timestep.begin());

		// Zero-fill newly added time steps for other spans
		const size_t added = new_steps - old_steps;
		if (added > 0) {
			// displacement
			std::memset(reinterpret_cast<uint8_t*>(timedisplacement.base + old_steps * ndisplacement * 3),
				0, disp_bytes(added));
			// stress
			std::memset(reinterpret_cast<uint8_t*>(timestress.base + old_steps * nstress * 6),
				0, stress_bytes(added));
			// nodalforce
			std::memset(reinterpret_cast<uint8_t*>(timenodalforce.base + old_steps * nnodalforce * 3),
				0, nforce_bytes(added));
			// constraint pressure
			std::memset(reinterpret_cast<uint8_t*>(timeconstraintpressure.base + old_steps * nconstraint),
				0, cpress_bytes(added));
			// constraint activate
			std::memset(reinterpret_cast<uint8_t*>(timeconstraintactivate.base + old_steps * nconstraint),
				0, cact_bytes(added));
		}
	}
	else {
		// SHRINK: move blocks forward to new positions while mapped to old size, then truncate
		// Keep current mapping (old size)
		uint8_t* base_old = (uint8_t*)this->mmap.data();
		if (!base_old) throw std::runtime_error("mio mapping missing (shrink)");

		// Move blocks "down" (dest < src). memmove is overlap-safe.
		// Order: from first to last is fine, memmove handles overlap. We skip timestep; we will overwrite.
		std::memmove(base_old + off_disp_new, base_old + off_disp_old, disp_bytes(new_steps));   // only keep first new_steps
		std::memmove(base_old + off_stress_new, base_old + off_stress_old, stress_bytes(new_steps));
		std::memmove(base_old + off_nforce_new, base_old + off_nforce_old, nforce_bytes(new_steps));
		std::memmove(base_old + off_cpress_new, base_old + off_cpress_old, cpress_bytes(new_steps));
		std::memmove(base_old + off_cact_new, base_old + off_cact_old, cact_bytes(new_steps));

		// Unmap old, truncate, remap new size
		this->mmap.unmap();
		std::filesystem::resize_file(this->dumpfile, total_new);
		std::error_code ec;
		this->mmap.map(this->dumpfile, 0, mio::map_entire_file, ec);
		if (ec) {
			throw std::runtime_error("mmap failed: " + ec.message());
		}

		uint8_t* base_new = (uint8_t*)this->mmap.data();

		// Rebuild spans on new layout
		rebuild_spans(base_new, new_steps);

		// Write pending timestep (new length)
		std::copy(new_ts.begin(), new_ts.end(), configure.timestep.begin());
		// No need to clear tail; file is truncated.
	}

	// Commit state
	this->total_steps = new_steps;
	configure.pending_timestep.reset();
}