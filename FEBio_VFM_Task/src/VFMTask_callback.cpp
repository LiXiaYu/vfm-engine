#include "VFMTask_callback.h"
#include "common_pybind.h"

#include <vector>
#include <set>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <filesystem>
#include <tuple>

#include <omp.h>
#include <nlopt.hpp>



bool read_inited_information(FEModel* fem, unsigned int when, void* pd)
{
	VFMTask* task = (VFMTask*)pd;

	::std::ofstream outFile;
	outFile.open(task->outlogfile, ::std::ios::app);
	if (!outFile.is_open())
	{
		write_log(fem, 2, "Open file \"");
		write_log(fem, 2, "log.txt");
		write_log(fem, 2, "\" failed.\n");
		return false;
	}

	int steps_number = fem->Steps();
	task->total_steps = 0;
	for (int i = 0; i < steps_number; i++)
	{
		FEAnalysis* step = fem->GetStep(i);
		FETimeStepController* timecontroller = step->m_timeController;
		
		double tspan = step->m_ntime * step->m_dt;
		double dtmin = timecontroller->m_dtmin;

		// 防御性检查，避免除零
		if (dtmin <= 0) {
			// 如果没设置 dtmin，就退化为用 m_ntime 作为上限
			task->total_steps += step->m_ntime;
		}
		else {
			int maxSteps = static_cast<int>(std::ceil(tspan / dtmin));
			task->total_steps += maxSteps;
		}
	}
	task->total_steps++; // include the initial state
	task->ndisplacement = 0;
	task->nstress = 0;
	task->nnodalforce = 0;
	task->nconstraint = 0; // pressure and activate

	DataStore& datastore = fem->GetDataStore();
	int data_number = datastore.Size();
	for (int j = 0; j < data_number; j++)
	{
		if (j == task->configure.displacementdataNumber)
		{
			DataRecord& datarecord = *(datastore.GetDataRecord(j));
			task->ndisplacement = datarecord.m_item.size();
		}
		else if (j == task->configure.stressdataNumber)
		{
			DataRecord& datarecord = *(datastore.GetDataRecord(j));
			task->nstress = datarecord.m_item.size();
		}
	}

	for (int j = 0; j < task->configure.pressure_load.size(); j++)
	{
		for (int i = 0; i < fem->ModelLoads(); i++)
		{
			FEModelLoad& load = *(fem->ModelLoad(i));
			::std::string loadclassname = load.GetFactoryClass()->GetClassName();
			::std::string loadname = load.GetName();

			if (::std::get<0>(task->configure.pressure_load[j]) == loadclassname && ::std::get<1>(task->configure.pressure_load[j]) == loadname)
			{
				if (loadclassname == "FENodalForce")
				{
					FENodalLoad& fnl = dynamic_cast<FENodalLoad&>(load);

					int dofs = fnl.GetDOFList().Size();
					::std::vector<double> val(dofs, 0.0);

					auto& nset = *(fnl.GetNodeSet());
					task->nnodalforce = nset.Size();
				}
			}
		}
	}

	for (int j = 0; j < task->configure.constraint_load.size(); j++)
	{
		int index_i = 0;
		for (int i = 0; i < fem->NonlinearConstraints(); i++)
		{
			FENLConstraint& fnc = *(fem->NonlinearConstraint(i));
			::std::string fncclassname = fnc.GetFactoryClass()->GetClassName();
			::std::string fncname = fnc.GetName();

			write_to_log_2(fem, ::std::string("fncclassname: ") + fncclassname + " fncname: " + fncname + "\n", outFile);

			if (::std::get<0>(task->configure.constraint_load[j]) == fncclassname && ::std::get<1>(task->configure.constraint_load[j]) == fncname)
			{
				if (fncclassname == "FEVolumeConstraint")
				{
					FEVolumeConstraint& vc = dynamic_cast<FEVolumeConstraint&>(fnc);
					//double pressure = vc.Pressure();


					bool activate = vc.IsActive();

					FESurface& surface = *(vc.GetSurface());

					double pressure = dynamic_cast<FEVolumeSurface&>(surface).m_p;

					task->nconstraint++;
				}
			}
		}
	}

	size_t time_size = task->total_steps * sizeof(double); // task->configure.timestep;
	size_t displacement_size = task->total_steps * task->ndisplacement * 3 * sizeof(double); // task->timedisplacement; //ux, uy, uz
	size_t stress_size = task->total_steps * task->nstress * 6 * sizeof(double); // task->timestress; //sx, sy, sz, sxy, syz, szx
	size_t nodalforce_size = task->total_steps * task->nnodalforce * 3 * sizeof(double); // task->timenodalforce; //fx,fy,fz
	size_t constraintpressure_size = task->total_steps * task->nconstraint * sizeof(double); // task->timeconstraintpressure;
	size_t constraintactivate_size = task->total_steps * sizeof(uint8_t); // task->timeconstraintactivate;

	size_t filesize = time_size + displacement_size + stress_size + nodalforce_size + constraintpressure_size + constraintactivate_size;

	if (task->configure.isRead_FEMresult_fromsavefile == false)
	{
		write_to_log_2(fem, "create new dump file...\n", outFile);
		// 不读取保存结果，创建新文件并覆盖
		if (std::filesystem::exists(task->dumpfile)) {
			write_to_log_2(fem, "dump file exist, remove it: " + task->dumpfile + "\n", outFile);
			std::filesystem::remove(task->dumpfile);
		}
		std::ofstream ofs(task->dumpfile, std::ios::binary | std::ios::trunc);
		ofs.seekp(filesize - 1);
		ofs.write("", 1);
		ofs.close();
	}
	else
	{
		write_to_log_2(fem, "read from savefile dump...\n", outFile);
		// 读取保存结果，检查文件是否存在 //且大小匹配
		if (!std::filesystem::exists(task->dumpfile)) {
			write_to_log_2(fem, "dump file not exist: " + task->dumpfile + "\n", outFile);
			if (task->configure.isSetDisplacmentAndPressure == true)
			{
				// 如果设置了从保存结果读取位移和压力，但 dump 文件不存在，则创建新文件
				write_to_log_2(fem, "since isSetDisplacmentAndPressure is true, create new dump file: " + task->dumpfile + "\n", outFile);
				std::ofstream ofs(task->dumpfile, std::ios::binary | std::ios::trunc);
				ofs.seekp(filesize - 1);
				ofs.write("", 1);
				ofs.close();
			}
			else
			{
				throw std::runtime_error("dump file not exist: " + task->dumpfile);
			}
		}

		//size_t existing_filesize = std::filesystem::file_size(task->dumpfile);
		//if (existing_filesize != filesize) {
		//	std::stringstream ss;
		//	ss << "dump file size mismatch: expected " << filesize << ", got " << existing_filesize;
		//	write_to_log_2(fem, ss.str() + "\n", outFile);
		//	throw std::runtime_error(ss.str());
		//}
	}

	// 映射到mio mmap
	std::error_code ec;
	task->mmap.map(task->dumpfile, 0, mio::map_entire_file, ec);
	if (ec) {
		throw std::runtime_error("mmap failed: " + ec.message());
	}

	char* base = task->mmap.data();

	// 建立视图
	task->configure.timestep = span<double>{ reinterpret_cast<double*>(base), task->total_steps };
	base += time_size;
	task->timedisplacement = span3d<double>{ reinterpret_cast<double*>(base), task->total_steps, task->ndisplacement, 3 };
	base += displacement_size;
	task->timestress = span3d<double>{ reinterpret_cast<double*>(base), task->total_steps, task->nstress, 6 };
	base += stress_size;
	task->timenodalforce = span3d<double>{ reinterpret_cast<double*>(base), task->total_steps, task->nnodalforce, 3 };
	base += nodalforce_size;
	task->timeconstraintpressure = span2d<double>{ reinterpret_cast<double*>(base), task->total_steps, task->nconstraint };
	base += constraintpressure_size;
	task->timeconstraintactivate = span2d<uint8_t>{ reinterpret_cast<uint8_t*>(base), task->total_steps, task->nconstraint };
	base += constraintactivate_size;

	return everytimestep_withinited_savedata(fem, when, pd);
}

bool read_stepsolved_information(FEModel* fem, unsigned int when, void* pd)
{
	return everytimestep_withinited_savedata(fem, when, pd);
}

bool everytimestep_withinited_savedata(FEModel* fem, unsigned int when, void* pd)
{
	VFMTask* task = (VFMTask*)pd;

	::std::ofstream outFile;
	outFile.open(task->outlogfile, ::std::ios::app);
	if (!outFile.is_open())
	{
		write_log(fem, 2, "Open file \"");
		write_log(fem, 2, "log.txt");
		write_log(fem, 2, "\" failed.\n");
		return false;
	}

	if (task->configure.isRead_FEMresult_fromsavefile == false)
	{
		static int INIT_timestep = 0;

		// get timestep
		double currenttime = fem->GetCurrentTime();
		::std::string logs_ct = "currenttime " + to_string(currenttime) + "\n";
		write_log(fem, 0, logs_ct.c_str());
		outFile << logs_ct;
		
		FEAnalysis* currenttimestep = fem->GetCurrentStep();
		int index_timestep = currenttimestep->m_ntimesteps + INIT_timestep;
		task->configure.timestep[index_timestep] = currenttime; 

		if (INIT_timestep == 0)
		{
			INIT_timestep = 1;
		}

		// get displacement
		DataStore& datastore = fem->GetDataStore();
		int data_number = datastore.Size();
		::std::string logs_dn = "data_number " + to_string(data_number) + "\n";
		write_log(fem, 0, logs_dn.c_str());
		outFile << logs_dn;

		//data_number = 1; // only get the first data
		for (int j = 0; j < data_number; j++)
		{
			if (j == task->configure.displacementdataNumber)
			{
				DataRecord& datarecord = *(datastore.GetDataRecord(j));
				::std::string data_name = datarecord.GetName();
				::std::string logs = "data " + to_string(j) + ": " + data_name + "\n";
				write_log(fem, 0, logs.c_str());
				outFile << logs;

				// get log file's data
				for (size_t i = 0; i < datarecord.m_item.size(); ++i)
				{
					int nd = datarecord.Size();
					//::std::vector<double> u(nd);
					//for (int k = 0; k < nd; ++k)
					//{
					//	double val = datarecord.Evaluate(datarecord.m_item[i], k);
					//	
					//	//std::stringstream ss;
					//	//ss << std::setprecision(12) << val;
					//	//::std::string logs = "item " + to_string(i) + " xyz " + to_string(k) + ":" + ss.str() + "\n";
					//	//write_log(fem, 0, logs.c_str());
					//	//outFile << logs;

					//	u[k] = val;
					//}
					//task->timedisplacement[index_timestep][i] = u;
					for (int k = 0; k < nd; ++k)
					{
						double val = datarecord.Evaluate(datarecord.m_item[i], k);
						task->timedisplacement[index_timestep][i][k] = val;
					}
				}
			}
			else if (j == task->configure.stressdataNumber)
			{
				DataRecord& datarecord = *(datastore.GetDataRecord(j));
				::std::string data_name = datarecord.GetName();
				::std::string logs = "data " + to_string(j) + ": " + data_name + "\n";
				write_log(fem, 0, logs.c_str());
				outFile << logs;

				// get log file's data
				for (size_t i = 0; i < datarecord.m_item.size(); ++i)
				{
					int nd = datarecord.Size();
					for (int k = 0; k < nd; ++k)
					{
						double val = datarecord.Evaluate(datarecord.m_item[i], k);
						//std::stringstream ss;
						//ss << std::setprecision(12) << val;
						//::std::string logs = "item " + to_string(i) + " xyz " + to_string(k) + ":" + ss.str() + "\n";
						//write_log(fem, 0, logs.c_str());
						//outFile << logs;

						task->timestress[index_timestep][i][k] = val;
					}
					
				}
			}
		}
	
	
		::std::vector<::std::vector<double>> currentTimeNodalForce;
		for (int j = 0; j < task->configure.pressure_load.size(); j++)
		{
			for (int i = 0; i < fem->ModelLoads(); i++)
			{
				FEModelLoad& load = *(fem->ModelLoad(i));
				::std::string loadclassname = load.GetFactoryClass()->GetClassName();
				::std::string loadname = load.GetName();

				if (::std::get<0>(task->configure.pressure_load[j]) == loadclassname && ::std::get<1>(task->configure.pressure_load[j]) == loadname)
				{
					if (loadclassname == "FENodalForce")
					{
						FENodalLoad& fnl = dynamic_cast<FENodalLoad&>(load);

						int dofs = fnl.GetDOFList().Size();
						::std::vector<double> val(dofs, 0.0);

						auto& nset = *(fnl.GetNodeSet());
						for (int i = 0; i < nset.Size(); ++i)
						{
							auto& node = *nset.Node(i);
							int nodeid = node.GetID();

							auto m_r0 = node.m_r0;
							auto m_rt = node.m_rt;

							
							// get the nodal values
							fnl.GetNodalValues(i, val);

							currentTimeNodalForce.push_back(val);

							// write to log file
							std::stringstream ss;
							ss << std::setprecision(12) << nodeid;
							::std::string logs = "nodeid: " + ss.str() + "\n";
							write_log(fem, 0, logs.c_str());
							outFile << logs;
							for (int k = 0; k < dofs; ++k)
							{
								ss.str("");
								ss << std::setprecision(12) << val[k];
								logs = "val[" + to_string(k) + "]:" + ss.str() + "\n";
								write_log(fem, 0, logs.c_str());
								outFile << logs;
							}

							::std::copy(val.begin(), val.end(), task->timenodalforce[index_timestep][i].begin());
						}
					}
				}
			}
		}


		int nlc = fem->NonlinearConstraints();
		int spc = fem->SurfacePairConstraints();
		write_to_log_2(fem, ::std::string("nlc: ") + ::std::to_string(nlc) + "spc" + ::std::to_string(spc), outFile);
		::std::vector<double> currentPressure;
		::std::vector<bool> currentActivate;
		for (int j = 0; j < task->configure.constraint_load.size(); j++)
		{

			int index_i = 0;
			for (int i = 0; i < fem->NonlinearConstraints(); i++)
			{
				FENLConstraint& fnc = *(fem->NonlinearConstraint(i));
				::std::string fncclassname = fnc.GetFactoryClass()->GetClassName();
				::std::string fncname = fnc.GetName();

				write_to_log_2(fem, ::std::string("fncclassname: ") + fncclassname + " fncname: " + fncname + "\n", outFile);

				if (::std::get<0>(task->configure.constraint_load[j]) == fncclassname && ::std::get<1>(task->configure.constraint_load[j]) == fncname)
				{
					if (fncclassname == "FEVolumeConstraint")
					{
						FEVolumeConstraint& vc = dynamic_cast<FEVolumeConstraint&>(fnc);
						//double pressure = vc.Pressure();
						
						
						bool activate = vc.IsActive();
						
						FESurface& surface = *(vc.GetSurface());

						double pressure = dynamic_cast<FEVolumeSurface&>(surface).m_p;

						task->timeconstraintpressure[index_timestep][index_i] = pressure;
						task->timeconstraintactivate[index_timestep][index_i] = activate ? 1 : 0;

						write_to_log_2(fem, ::std::string("pressure: ") + ::std::to_string(pressure), outFile);
						index_i++;
					}
				}
			}
		}
	}
	else
	{

	}


	outFile.close();

	return true;
}


bool read_solved_information(FEModel* fem, unsigned int when, void* pd)
{
	VFMTask* task = (VFMTask*)pd;

	::std::ofstream outFile;
	outFile.open(task->outlogfile, ::std::ios::app);
	if (!outFile.is_open())
	{
		write_log(fem, 2, "Open file \"");
		write_log(fem, 2, "log.txt");
		write_log(fem, 2, "\" failed.\n");
		return false;
	}

	// get initial initialCoordinate
	write_log(fem, 0, "read_solved_information!!!\n");
	outFile << "read_solved_information!!!\n";
	// not initial initialCoordinate, but mesh

	const bool ONH_sinPulse100s_120_v4 = true;
	const bool BPM_base_v4_BPM = false;

	::std::vector<double> timeArray;
	::std::vector<int> solution_elementsDomainID;
	::std::vector<::std::vector<::std::vector<double>>> trueJArray;
	::std::vector<::std::vector<::std::vector<mat3d>>> truedeformationGradientArray;
	::std::vector<std::vector<std::vector<std::vector<mat3ds>>>> virtualstrainArrayV;
	::std::vector<::std::vector<double>> externalVirtualWork;
	::std::vector<::std::vector<double>> volumeVirtualWork;

	::std::vector<::std::vector<double>> true_internalVirtualWork;

	::std::vector<::std::vector< std::complex<double>>> externalVirtualWork_laplace;
	// [index_vf, index_ElementID, index_gauss, index_laplace_s]
	::std::vector<::std::vector<::std::vector<::std::vector<::std::complex<double>>>>> Sepk2_dotdot_vStrain;

	::std::vector<::std::vector<::std::complex<double>>> fs;
	::std::vector<::std::vector<double>> ft;

	::std::function<double(const ::std::vector<double>&)> function;

	::std::unique_ptr<FunOptimParams> params;

	if (task->configure.isReadfromsaveOptimfunc == false)
	{

#pragma region getInput
		// get mesh's object points
		FEMesh& mesh = fem->GetMesh();

		write_to_log_2(fem, "read initial Coordinate...\n", outFile);

		int nodes_number = mesh.Nodes();

		task->nodes.resize(nodes_number);
		task->configure.initialCoordinate.resize(nodes_number);

		// This parallel for loop is not necessary, because it is not time-consuming
		//#pragma omp parallel for
		for (int j = 0; j < nodes_number; j++)
		{
			vec3d& node = mesh.Node(j).m_r0;
			//node.x;
			//node.y;
			//node.z;
			task->nodes[j] = node;
			task->configure.initialCoordinate[j] = ::std::vector<double>{ node.x, node.y, node.z };
		}

		write_to_log_2(fem, "get solutions element...\n", outFile);
		// get solutions' element
		for (int i = 0; i < task->configure.solution.size(); i++)
		{
			::std::string solution_type = ::std::get<0>(task->configure.solution[i]);
			::std::string solution_name = ::std::get<1>(task->configure.solution[i]);

			if (solution_type == "elements")
			{
				for (int j = 0; j < mesh.ElementSets(); j++)
				{
					FEElementSet& elementset = mesh.ElementSet(j);

					::std::string elementset_name = elementset.GetName();
					if (elementset_name == solution_name)
					{
						::std::vector<int> toadd = elementset.GetElementIDList();
						//// foreach toadd element -=1
						for (int k = 0; k < toadd.size(); k++)
						{
							toadd[k] = mesh.FindElementIndexFromID(toadd[k]);
						}
						// 现在储存真实id，不用减1。甚至不连续，任意编号。
						// 现在改成内部index

						task->solution_elementsID.insert(task->solution_elementsID.end(), toadd.begin(), toadd.end());
					}
				}
			}

		}

		write_to_log_2(fem, "get fixed nodes...\n", outFile);
		// get fixed nodes
		for (int i = 0; i < task->configure.fixed.size(); i++)
		{
			::std::string fixed_type = ::std::get<0>(task->configure.fixed[i]);
			::std::string fixed_name = ::std::get<1>(task->configure.fixed[i]);

			if (fixed_type == "nodeset")
			{
				for (int j = 0; j < mesh.NodeSets(); j++)
				{
					FENodeSet& nodeset = *(mesh.NodeSet(j));
					::std::string nodeset_name = nodeset.GetName();
					if (nodeset_name == fixed_name)
					{
						for (int k = 0; k < nodeset.Size(); k++)
						{
							int node_id = nodeset[k];

							task->configure.fixednode.push_back(node_id);
						}

						break;
					}
				}
			}
			else if (fixed_type == "surface")
			{
				FEFacetSet* facetset = mesh.FindFacetSet(fixed_name);
				if (facetset != nullptr)
				{
					FENodeList&& nodelist = facetset->GetNodeList();
					for (int k = 0; k < nodelist.Size(); k++)
					{
						int node_id = nodelist[k];

						task->configure.fixednode.push_back(node_id);
					}

				}
			}
		}

#pragma endregion

#pragma region getAllInput
		//number of cycles considered at the end
		int numCycle = 1;
		int nNode = nodes_number;


		//::std::vector<int> sidenodes = { 12, 23, 56, 67, 66, 77, 22, 33 };
		//// sidenodes -= 1
		//for (int i = 0; i < sidenodes.size(); i++)
		//{
		//	sidenodes[i] -= 1;
		//}

		//// joint task->configure.fixednode and sidenodes
		//task->configure.fixednode.insert(task->configure.fixednode.end(), sidenodes.begin(), sidenodes.end());



		//%% displacement at the initial timestep of the last two cycles
#pragma region readsavefile_dump

		//write_to_log_2(fem, "read from savefile dump...\n", outFile);
		//if (task->configure.isRead_FEMresult_fromsavefile == true)
		//{
		//	DumpFile dumpfile(*fem);
		//	dumpfile.Open(task->dumpfile.c_str());

		//	dumpfile& task->configure.timestep& task->timedisplacement& task->timestress &task->timenodalforce &task->timeconstraintpressure &task->timeconstraintactivate;

		//	dumpfile.Close();

		//}
		//else
		//{
		//	DumpFile dumpfile(*fem);
		//	dumpfile.Create(task->dumpfile.c_str());

		//	dumpfile& task->configure.timestep& task->timedisplacement& task->timestress &task->timenodalforce& task->timeconstraintpressure& task->timeconstraintactivate;

		//	dumpfile.Close();
		//}

		//// log: read save file end
		//write_to_log_2(fem, "read save file end\n", outFile);

#pragma endregion


		FEAnalysis& laststep = *(fem->GetStep(fem->Steps() - 1));

		int start_index = 0;
		int end_index = 0;

		try
		{
			write_to_log_2(fem, "set timestep from SetPeriodTimeIndex...\n", outFile);
			//task->configure = pyfunctioncall<VFMTask_configure>(fem, "", outFile, task->pyfile_module, "SetPeriodTimeIndex", & task->configure);
			pyfunctioncall(fem, "", outFile, task->pyfile_module, "SetPeriodTimeIndex", &task->configure);
			start_index = task->configure.start_index;
			end_index = task->configure.end_index;
		}
		catch(const std::exception& e)
		{
			// log: calculate timestep for VFM begin...
			write_to_log_2(fem, "calculate timestep for VFM begin...\n", outFile);

			auto [_start_index, _end_index] = find_last_two_cycles(task->configure.timestep, task->configure.bpm);
			start_index = _start_index;
			end_index = _end_index;
		}

		// log: calculate timestep for VFM end...
		write_to_log_2(fem, "calculate timestep for VFM end...\n", outFile);
		// log: start_index, sum_sycle_time
		write_to_log_2(fem, "start_index: " + to_string(start_index) + "\n", outFile);
		// log: end_index, sum_sycle_time
		write_to_log_2(fem, "end_index: " + to_string(end_index) + "\n", outFile);

		//int steps_per_cycle = 10; // 
		//int nStep_total = task->timestep.size();

		//int start_index = nStep_total - steps_per_cycle * numCycle;


		::std::vector<::std::function<::std::vector<double>(const ::std::vector<double>&, int)>> vf_u_functions;

		if constexpr (ONH_sinPulse100s_120_v4)
		{
			double radius = 0.8;
			double height = 0.4;

			double x_r_0 = 0.0;
			double y_r_0 = 0.0;
			double z_r_0 = 0.0;

			double x_r_t = 1.0;
			double y_r_t = 1.0;
			double z_r_t = 1.0;
		}


		// BPM_base_v4_BPM
#pragma region BPM_base_v4_BPM
		if constexpr (BPM_base_v4_BPM)
		{
			double radius = 1.6;
			double height = 1.0;

			double x_r_0 = 142.185;
			double y_r_0 = 141.782;
			double z_r_0 = 0;
			double x_r_t = 1.0;
			double y_r_t = 1.0;
			double z_r_t = 1.0;

			double x_min = 140.185;
			double x_max = 144.185;
			double y_min = 141.289;
			double y_max = 143.289;
			double z_min = -2.0;
			double z_max = 2.0;

			double x_p = 1;
			double y_p = 1;
			double z_p = 1;
		}

		try
		{
			// py call SetSelectSolutionElementFunction
			auto selectsolutionelement_configure_from_py = task->pyfile_module.attr("SetSelectSolutionElementFunction")(task->configure).cast<VFMTask_configure>();
			// 从python文件中获取配置
			task->configure = selectsolutionelement_configure_from_py;

		}
		catch (::std::exception& e)
		{
			::std::string logs = "error: " + ::std::string(e.what()) + "\n";
			write_log(fem, 2, logs.c_str());
			outFile << logs;
		}

		pybind11::gil_scoped_release release;

		int element_number = mesh.Elements();
		// log element nubmer
		write_log(fem, 0, ("element_number: " + to_string(element_number) + "\n").c_str());

		FEElement* pelement0 = mesh.Element(0);
		FEElement* pelement643911 = mesh.Element(643911);
		FEElement* pelement0_ = mesh.FindElementFromID(0);
		FEElement* pelement643911_ = mesh.FindElementFromID(643911);
		mesh.FindElementIndexFromID(643911);
		//log the two element is null or not
		if (pelement0 == nullptr)
		{
			write_log(fem, 2, "error: element 0 is null\n");
		}
		if (pelement643911 == nullptr)
		{
			write_log(fem, 2, "error: element 643911 is null\n");
		}
		if (pelement0_ == nullptr)
		{
			write_log(fem, 2, "error: FindElementFromID 0 is null\n");
		}
		if (pelement643911_ == nullptr)
		{
			write_log(fem, 2, "error: FindElementFromID 643911 is null\n");
		}

		auto it = ::std::remove_if(task->solution_elementsID.begin(), task->solution_elementsID.end(),
			[=, &mesh, &task](int j)->bool {
				FEElement& element = *(mesh.Element(j));
				
				// 检查element 是否为null
				if (mesh.Element(j) == nullptr)
				{
					// log error j, task->solution_elementsID[j]
					::std::stringstream ss;
					ss << "Error: element is null, elementID j: " << j << "\n";
					write_log(fem, 2, ss.str().c_str());
					return true;
				}

				int element_id = element.GetID();
				::std::vector<int> nodeids(element.Nodes());
				::std::vector<::std::vector<double>> node_xyzs(element.Nodes());

				for (int k = 0; k < element.Nodes(); k++)
				{
					int node_id = element.m_node[k];
					nodeids[k] = node_id;

					auto& node = mesh.Node(node_id);

					double x = task->configure.initialCoordinate[node_id][0];
					double y = task->configure.initialCoordinate[node_id][1];
					double z = task->configure.initialCoordinate[node_id][2];

					node_xyzs[k] = { x, y, z };

				}
				pybind11::gil_scoped_acquire acquire;
				::std::vector<int> no_inrange_point_number = task->configure.select_solution_element_function(element_id, nodeids, node_xyzs);


				if (no_inrange_point_number.size() == element.Nodes())
				{
					return true;
				}
				else
				{
					// add node id to fixednode
					for (int k = 0; k < no_inrange_point_number.size(); k++)
					{
						task->configure.fixednode.push_back(element.m_node[no_inrange_point_number[k]]);
					}
					return false;
				}
			});
		task->solution_elementsID.erase(it, task->solution_elementsID.end());
		pybind11::gil_scoped_acquire acquire;

		// 将pressure_load和constraint_load的节点移出fixednode

		// 用集合统计node
		::std::set<int> all_load_node_set;
		for (int j = 0; j < task->configure.pressure_load.size(); j++)
		{
			for (int i = 0; i < fem->ModelLoads(); i++)
			{
				FEModelLoad& load = *(fem->ModelLoad(i));
				::std::string loadclassname = load.GetFactoryClass()->GetClassName();
				::std::string loadname = load.GetName();

				if (::std::get<0>(task->configure.pressure_load[j]) == loadclassname && ::std::get<1>(task->configure.pressure_load[j]) == loadname)
				{
					if (loadclassname == "FEPressureLoad")
					{
						FEPressureLoad& pressureLoad = dynamic_cast<FEPressureLoad&>(load);

						auto& surface = pressureLoad.GetSurface();

						auto nset = surface.GetNodeList();
						for (int i = 0; i < nset.Size(); ++i)
						{
							auto& node = *nset.Node(i);
							int node_index = node.GetID()-1;

							all_load_node_set.insert(node_index);
						}
					}
					else if (loadclassname == "FENodalForce")
					{
						FENodalLoad& fnl = dynamic_cast<FENodalLoad&>(load);

						auto& nset = *(fnl.GetNodeSet());
						for (int i = 0; i < nset.Size(); ++i)
						{
							auto& node = *nset.Node(i);
							int node_index = node.GetID()-1;

							all_load_node_set.insert(node_index);
						}
					}

				}
			}
			
		}
		for (int j = 0; j < task->configure.constraint_load.size(); j++)
		{
			for (int i = 0; i < fem->NonlinearConstraints(); i++)
			{
				FENLConstraint& fnc = *(fem->NonlinearConstraint(i));
				::std::string fncclassname = fnc.GetFactoryClass()->GetClassName();
				::std::string fncname = fnc.GetName();

				if (::std::get<0>(task->configure.constraint_load[j]) == fncclassname && ::std::get<1>(task->configure.constraint_load[j]) == fncname)
				{
					if (fncclassname == "FEVolumeConstraint")
					{
						FEVolumeConstraint& vc = dynamic_cast<FEVolumeConstraint&>(fnc);

						FESurface& surface = *(vc.GetSurface());

						auto nset = surface.GetNodeList();
						for (int i = 0; i < nset.Size(); ++i)
						{
							auto& node = *nset.Node(i);
							int node_index = node.GetID()-1;

							all_load_node_set.insert(node_index);
						}
					}
				}
			}
		}
		// 移除fixednode中所有在all_load_node_set中的节点
		//task->configure.fixednode.erase(std::remove_if(task->configure.fixednode.begin(), task->configure.fixednode.end(), 
		//	[&](int node) { 
		//		return all_load_node_set.find(node) != all_load_node_set.end(); 
		//	}), task->configure.fixednode.end());


		// show all element id in mesh.Domain(0)
		//FEDomain& d0= mesh.Domain(0);
		//int d0en=d0.Elements();
		//for (int i = 0; i < d0en; i++)
		//{
		//	FEElement& element = d0.ElementRef(i);
		//	int element_id = element.GetID();
		//	::std::string logs = "element id: " + to_string(element_id) + "\n";
		//	write_log(fem, 0, logs.c_str());
		//	outFile << logs;
		//}

		solution_elementsDomainID = ::std::vector<int>(task->solution_elementsID.size());
/*    */#pragma omp parallel for
		for (int index_seId = 0; index_seId < task->solution_elementsID.size(); index_seId++)
		{
			int elementid = task->solution_elementsID[index_seId];

			FEElement* element = mesh.Element(elementid);

			int domain_id = -1;
			for (int i = 0; i < mesh.Domains(); i++)
			{
				FEDomain& d = mesh.Domain(i);
				
				for (int j = 0; j < d.Elements(); j++)
				{
					FEElement& el = d.ElementRef(j);
					
					if (element->GetID() == el.GetID())
					{
						domain_id = i;
						break;
					}
				}

				if (domain_id != -1)
				{
					break;
				}

			}
			if (domain_id == -1)
			{
				throw ::std::runtime_error("no domain found for element");
			}

			solution_elementsDomainID[index_seId] = domain_id;
		}

		// string stream
		::std::stringstream elementset_ss;
		int line_number = 0;
		int line_number_max = 7;
		for (int index_seId = 0; index_seId < task->solution_elementsID.size(); index_seId++)
		{
			int elementid = task->solution_elementsID[index_seId];
			elementset_ss << (elementid);
			if (line_number >= line_number_max)
			{
				elementset_ss << ",\n";

				line_number = 0;
			}
			else
			{
				elementset_ss << ", ";

				line_number++;
			}
		}
		::std::string elementset_string = elementset_ss.str();

		::std::stringstream fixednodeset_ss;
		line_number = 0;
		line_number_max = 7;
		for (int index_fnId = 0; index_fnId < task->configure.fixednode.size(); index_fnId++)
		{
			int nodeid = task->configure.fixednode[index_fnId];
			fixednodeset_ss << (nodeid + 1);
			if (line_number >= line_number_max)
			{
				fixednodeset_ss << ",\n";

				line_number = 0;
			}
			else
			{
				fixednodeset_ss << ", ";

				line_number++;
			}
		}
		::std::string fixednodeset_string = fixednodeset_ss.str();


#pragma endregion


		// node_id list in solution_elementsID's elements
		::std::set<int> solution_node_id_set;
		for (int index_seId = 0; index_seId < task->solution_elementsID.size(); index_seId++)
		{
			int j = task->solution_elementsID[index_seId];

			FEElement& element = *(mesh.Element(j));
			for (int k = 0; k < element.Nodes(); k++)
			{
				int node_id = element.m_node[k];
				solution_node_id_set.insert(node_id);
			}
		}
		::std::vector<int> solution_node_id_list(solution_node_id_set.begin(), solution_node_id_set.end()); // should replaced by flat_set?


#ifdef DEBUG
		//// output solution_node_id_list  write to log2
		//::std::stringstream ssxxx;
		//ssxxx << "solution_node_id_list" << "\n";
		//for (int index_node_id = 0; index_node_id < solution_node_id_list.size(); index_node_id++)
		//{
		//	int node_id = solution_node_id_list[index_node_id];
		//	ssxxx << "node_id: " << node_id << "\n";
		//}
		//::std::string ssxxx_string = ssxxx.str();
		//int sss=ssxxx_string.size();
		//outFile << ssxxx_string;
#endif


		// near_node_list in solution_node_id_set;
		::std::vector<::std::set<int>> adjacency_list(solution_node_id_list.size());
		for (int index_seId = 0; index_seId < task->solution_elementsID.size(); index_seId++)
		{
			int j = task->solution_elementsID[index_seId];

			FEElement& element = *(mesh.Element(j));

			::std::vector<::std::vector<int>> element_adjacent;
			switch (element.Type())
			{
			case FE_Element_Type::FE_HEX8G1:
			case FE_Element_Type::FE_HEX8G8:
				element_adjacent = { {1,3,4},{0,2,5},{1,3,6},{0,2,7},
									 {0,5,7},{1,4,6},{2,5,7},{3,4,6} };
				break;
			case FE_Element_Type::FE_PENTA6G6:
				element_adjacent = { {1,2,3},{0,2,4},{0,1,5},
									 {0,4,5},{1,3,5},{2,3,4} };
				break;
			default:
				// all element's node
				// vector from 0 to element.Nodes() - 1
				::std::vector<int> m_lnode(element.Nodes());
				for (int i = 0; i < element.Nodes(); i++)
				{
					m_lnode[i] = i;
				}
				element_adjacent = ::std::vector<::std::vector<int>>(element.Nodes(), m_lnode);
				// remove each index in element_adjacent[i]
				for (int i = 0; i < element.Nodes(); i++)
				{
					element_adjacent[i].erase(
						::std::remove(element_adjacent[i].begin(), element_adjacent[i].end(), i),
						element_adjacent[i].end());
				}
				break;
			}

			auto m_node = element.m_node;
			for (int k = 0; k < element.Nodes(); k++)
			{
				int node_id = element.m_node[k];

				::std::vector<int> adjacent(element_adjacent[k].size());
				::std::transform(element_adjacent[k].begin(), element_adjacent[k].end(), adjacent.begin(), [m_node](int x) {return m_node[x]; });

				auto iter_node_id_solution = ::std::find(solution_node_id_list.begin(), solution_node_id_list.end(), node_id);

				if (iter_node_id_solution != solution_node_id_list.end())
				{
					int index__ = ::std::distance(solution_node_id_list.begin(), iter_node_id_solution);
					adjacency_list[index__].insert(adjacent.begin(), adjacent.end());
				}
				else
				{
					int index__ = index__;
					throw ::std::runtime_error("adjacent: no node_id in solution_node_id_list!");
				}


			}
		}

		if (task->configure.isSetDisplacmentAndPressure)
		{
			try
			{
				// py call SetDisplacmentAndPressureFunction
				auto sdp = task->pyfile_module.attr("SetDisplacmentAndPressureFunction")(task->configure);
				auto setdisplacementandpressure_configure_from_py = sdp.cast<VFMTask_configure>();
				// 从python文件中获取配置
				task->configure = setdisplacementandpressure_configure_from_py;
			}
			catch (::std::exception& e)
			{
				::std::string logs = "error: " + ::std::string(e.what()) + "\n";
				write_log(fem, 2, logs.c_str());
				outFile << logs;
			}

			//task->timedisplacement.resize(task->configure.timestep.size());
			//task->timestress.resize(task->configure.timestep.size());
			task->reinit_mmap_and_spans();

			auto displacements = ::std::vector<::std::vector<::std::vector<double>>>(task->configure.timestep.size());
			auto coordinates = ::std::vector<::std::vector<::std::vector<double>>>(task->configure.timestep.size());

			auto displacments_interpolated = ::std::vector<::std::vector<::std::vector<double>>>(task->configure.timestep.size());
			auto displacments_change_value = ::std::vector<::std::vector<::std::vector<double>>>(task->configure.timestep.size());

			for (int index_timestep = 0; index_timestep < task->configure.timestep.size(); index_timestep++)
			{
				displacements[index_timestep] = ::std::vector<::std::vector<double>>(solution_node_id_list.size());
				coordinates[index_timestep] = ::std::vector<::std::vector<double>>(solution_node_id_list.size());

				displacments_interpolated[index_timestep] = ::std::vector<::std::vector<double>>(solution_node_id_list.size());
				displacments_change_value[index_timestep] = ::std::vector<::std::vector<double>>(solution_node_id_list.size());
				/*            *///#pragma omp parallel for
				for (int index_node_id = 0; index_node_id < solution_node_id_list.size(); index_node_id++)
				{
					int node_id = solution_node_id_list[index_node_id];
					auto& node = mesh.Node(node_id);
					::std::vector<double> xyz_o = task->configure.initialCoordinate[node_id];

					bool debugtag = false;

					if (node_id == 11023)
					{
						debugtag = true;
					}

					::std::vector<double> c_displacement = task->configure.setdisplacement_function(xyz_o, index_timestep, debugtag);

					::std::copy(c_displacement.begin(), c_displacement.end(), task->timedisplacement[index_timestep][node_id].begin());
					//task->timedisplacement[index_timestep][node_id] = c_displacement;

					displacements[index_timestep][index_node_id] = c_displacement;
					coordinates[index_timestep][index_node_id] = { xyz_o[0] + c_displacement[0],xyz_o[1] + c_displacement[1],xyz_o[2] + c_displacement[2] };

				}

				// 平滑次数
				for (int time_InterpolatedDisplacementTime = 0; time_InterpolatedDisplacementTime < task->configure.InterpolatedDisplacementTimes; time_InterpolatedDisplacementTime++)
				{
					// 平滑，保证位移空间连续性
					auto temp_displacements = displacements[index_timestep];
					for (int index_node_id = 0; index_node_id < solution_node_id_list.size(); index_node_id++)
					{
						int node_id = solution_node_id_list[index_node_id];

						::std::vector<double> xyz_o = task->configure.initialCoordinate[node_id];

						::std::vector<int> adjacent_nodeid(adjacency_list[index_node_id].begin(), adjacency_list[index_node_id].end());


						// 反距离加权插值（IDW）
						int power = 2;
						::std::vector<double> interpolated_value = { 0,0,0 };
						::std::vector<double> displacement = span_to_vector(task->timedisplacement[index_timestep][node_id]);
						double weight_sum = 0;
						for (int index_adjacent_nodeid = 0; index_adjacent_nodeid < adjacent_nodeid.size(); index_adjacent_nodeid++)
						{
							int node_id_a = adjacent_nodeid[index_adjacent_nodeid];
							::std::vector<double> xyz_a = task->configure.initialCoordinate[node_id_a];
							::std::vector<double> displacement_a = span_to_vector(task->timedisplacement[index_timestep][node_id_a]);
							double distance = ::std::sqrt(::std::pow(xyz_o[0] - xyz_a[0], 2) + ::std::pow(xyz_o[1] - xyz_a[1], 2) + ::std::pow(xyz_o[2] - xyz_a[2], 2));
							double weight = ::std::pow(1.0 / distance, power);
							weight_sum += weight;
							interpolated_value[0] += weight * displacement_a[0];
							interpolated_value[1] += weight * displacement_a[1];
							interpolated_value[2] += weight * displacement_a[2];
						}
						interpolated_value[0] /= weight_sum;
						interpolated_value[1] /= weight_sum;
						interpolated_value[2] /= weight_sum;

						if (node_id == 10986)
						{
							auto pwokeopkmsoxlckmops = 33;
						}

						temp_displacements[index_node_id] = interpolated_value;
					}
					for (int index_node_id = 0; index_node_id < solution_node_id_list.size(); index_node_id++)
					{
						int node_id = solution_node_id_list[index_node_id];
						::std::copy(temp_displacements[index_node_id].begin(), temp_displacements[index_node_id].end(), task->timedisplacement[index_timestep][node_id].begin());
						//task->timedisplacement[index_timestep][node_id] = temp_displacements[index_node_id];

						displacments_interpolated[index_timestep][index_node_id] = temp_displacements[index_node_id];
						displacments_change_value[index_timestep][index_node_id] = { temp_displacements[index_node_id][0] - displacements[index_timestep][index_node_id][0],temp_displacements[index_node_id][1] - displacements[index_timestep][index_node_id][1] ,temp_displacements[index_node_id][2] - displacements[index_timestep][index_node_id][2] };
					}
				}



				for (int j = 0; j < task->configure.pressure_load.size(); j++)
				{

					for (int i = 0; i < fem->ModelLoads(); i++)
					{
						FEModelLoad& load = *(fem->ModelLoad(i));
						::std::string loadclassname = load.GetFactoryClass()->GetClassName();
						::std::string loadname = load.GetName();

						if (::std::get<0>(task->configure.pressure_load[j]) == loadclassname && ::std::get<1>(task->configure.pressure_load[j]) == loadname)
						{
							if (loadclassname == "FEPressureLoad")
							{
								FEPressureLoad& pressureLoad = static_cast<FEPressureLoad&>(load);

								auto& surface = pressureLoad.GetSurface();

								::std::vector<::std::vector<double>> evw_surface(surface.Elements());
								::std::vector<bool> evw_surface_flag(surface.Elements(), false);

								for (int index_surface_element = 0; index_surface_element < surface.Elements(); index_surface_element++)
								{
									auto& element = surface.Element(index_surface_element);

									auto temp_true_element= element.m_elem[0];

									FEElement* true_element = get_FEElement_p_version(temp_true_element);

									int element_index = mesh.FindElementIndexFromID(true_element->GetID());

									::std::vector<::std::vector<double>> xyz_o; // element center coordinate
									for (int k = 0; k < true_element->Nodes(); k++)
									{
										int node_id = true_element->m_node[k];
										auto& node = mesh.Node(node_id);
										::std::vector<double> xyz = { node.m_r0.x , node.m_r0.y,node.m_r0.z };
										xyz_o.push_back(xyz);

									}

									::std::vector<double> c_stress = task->configure.setstress_function(xyz_o, index_timestep);
									
									::std::copy(c_stress.begin(), c_stress.end(), task->timestress[index_timestep][element_index].begin());
									//task->timestress[index_timestep][elementId] = c_stress;

								}
							}
						}
					}
				}
			}

			//task->pyfile_module.attr("Display_3D_animation")(displacements, "./temp/debug/displacement/displacements.gif");
			//task->pyfile_module.attr("Display_3D_animation")(coordinates, "./temp/debug/displacement/coordinates.gif");

			//task->pyfile_module.attr("Display_3D_animation")(displacments_interpolated, "./temp/debug/displacement/displacements_interpolated.gif");
			//task->pyfile_module.attr("Display_3D_animation")(displacments_change_value, "./temp/debug/displacement/displacments_change_value.gif");

			start_index = 0;
		}

//#ifdef DEBUG
//		// print task->solution_elementsID
//		::std::ostringstream oss;
//		oss << "        <ElementSet name=\"solution_elementsID_add1\">" << ::std::endl << "            ";
//
//		int line_max_e = 8;
//		for (int i = 0; i < task->solution_elementsID.size(); i++)
//		{
//			if (task->solution_elementsID[i] == 57928)
//			{
//				auto pwokeopkmsoxlckmops = i;
//			}
//			oss << (task->solution_elementsID[i] + 1);
//			if (i == task->solution_elementsID.size() - 1)
//			{
//				oss << ::std::endl;
//			}
//			else if (i % line_max_e == line_max_e - 1)
//			{
//				oss << "," << ::std::endl << "            ";
//			}
//			else if (i != task->solution_elementsID.size() - 1)
//			{
//				oss << ", ";
//			}
//		}
//		oss << "        </ElementSet>" << ::std::endl;
//
//		::std::cout << oss.str();
//		//write_to_log_2(fem, oss.str(), outFile);
//
//		// find task->configure.fixednode and solution_node_id_list cross check
//		::std::vector<int> solution_fixednode;
//		for (int i = 0; i < task->configure.fixednode.size(); i++)
//		{
//			if (::std::find(solution_node_id_list.begin(), solution_node_id_list.end(), task->configure.fixednode[i]) != solution_node_id_list.end())
//			{
//				solution_fixednode.push_back(task->configure.fixednode[i]);
//			}
//		}
//
//		// print task->configure.fixednode
//		::std::ostringstream oss1;
//		oss1 << "        <NodeSet name=\"fixednode_add1\">" << ::std::endl << "            ";
//		for (int i = 0; i < solution_fixednode.size(); i++)
//		{
//			oss1 << (solution_fixednode[i] + 1);
//			if (i == solution_fixednode.size() - 1)
//			{
//				oss1 << ::std::endl;
//			}
//			else if (i % line_max_e == line_max_e - 1)
//			{
//				oss1 << "," << ::std::endl << "            ";
//			}
//			else if (i != solution_fixednode.size() - 1)
//			{
//				oss1 << ", ";
//			}
//		}
//		oss1 << "        </NodeSet>" << ::std::endl;
//
//		::std::cout << oss1.str();
//		//write_to_log_2(fem, oss1.str(), outFile);
//
//		//::std::getchar();
//#endif

		// copy ::std::vector task->timedisplacement[nNode*(start_index-1)+1:nNode*start_index] as initialDisp
		::std::vector<::std::vector<double>> initialDisp = span2d_to_vector2d(task->timedisplacement[start_index]);
		::std::vector<::std::vector<::std::vector<double>>> displacementArray = span3d_to_vector3d(task->timedisplacement.begin() + start_index, task->timedisplacement.begin() + end_index + 1);
		timeArray = span_to_vector(task->configure.timestep.begin() + start_index, task->configure.timestep.begin() + end_index+1);


		//convert task->timestress[:][:][0,1,2,3,4,5] -> task->timestress[:][:][0,3,5,1,4,2]
  //      #pragma omp parallel for	
		//for (int i = 0; i < task->timestress.size(); i++)
		//{
		//	for (int j = 0; j < task->timestress[i].size(); j++)
		//	{
		//		::std::vector<double> temp_stress(task->timestress[i][j]);
		//		task->timestress[i][j] = ::std::vector<double>({ temp_stress[0],temp_stress[3],temp_stress[5],temp_stress[1],temp_stress[4],temp_stress[2] }); // [0,3,5,1,2,4]
		//	}
		//}


		::std::vector<::std::vector<double>> initialStress = span2d_to_vector2d(task->timestress[start_index]);
		::std::vector<::std::vector<::std::vector<double>>> stressArray = span3d_to_vector3d(task->timestress.begin() + start_index, task->timestress.begin() + end_index + 1);
		
		::std::vector<double> initialConstraintPressure = span_to_vector(task->timeconstraintpressure[start_index]);
		::std::vector<::std::vector<double>> constraintPressureArray = span2d_to_vector2d(task->timeconstraintpressure.begin() + start_index, task->timeconstraintpressure.begin() + end_index + 1);

		::std::vector<uint8_t> initialConstraintActivate = span_to_vector(task->timeconstraintactivate[start_index]);
		::std::vector<::std::vector<uint8_t>> constraintActivateArray = span2d_to_vector2d(task->timeconstraintactivate.begin() + start_index, task->timeconstraintactivate.begin() + end_index + 1);

		write_to_log_2(fem, ::std::string("task->timeconstraintpressure size: ")+::std::to_string(task->timeconstraintpressure.rows())+"\n", outFile);
		write_to_log_2(fem, ::std::string("task->timeconstraintactivate size: ")+::std::to_string(task->timeconstraintactivate.rows())+"\n", outFile);

		// write_to_log_2 show task->timeconstraintpressure
		for (int i = 0; i < constraintPressureArray.size(); i++)
		{
			write_to_log_2(fem, ::std::string("index i: ") + ::std::to_string(i) + "\n", outFile);
			for (int j = 0; j < constraintPressureArray[i].size(); j++)
			{
				
				::std::ostringstream oss;
				oss <<::std::setprecision(12) << constraintPressureArray[i][j];
				if(j == constraintPressureArray[i].size() - 1 && i != constraintPressureArray.size() - 1)
				{
					oss << ::std::endl;
				}
				else if (j != constraintPressureArray[i].size() - 1)
				{
					oss << ", ";
				}
				write_to_log_2(fem, oss.str(), outFile);
			}
		}

		for(int i = 0; i < constraintActivateArray.size(); i++)
		{
			for (int j = 0; j < constraintActivateArray[i].size(); j++)
			{
				::std::ostringstream oss;
				oss << ::std::setprecision(12) << constraintActivateArray[i][j];
				if (j == constraintActivateArray[i].size() - 1 && i != constraintActivateArray.size() - 1)
				{
					oss << ::std::endl;
				}
				else if (j != constraintActivateArray[i].size() - 1)
				{
					oss << ", ";
				}
				write_to_log_2(fem, oss.str(), outFile);
			}
		}


		if (task->configure.isLaplaceVFM == true)
		{
			if (task->configure.isTalbotLaplaceVFM_s == true)
			{

				int N = 32;
				double shift = CONST_SHIFT;
				auto ss = createTalbotPath(timeArray, N, shift);

				::std::vector<::std::complex<double>> laplace_s;
				for (int i = 0; i < timeArray.size(); i++)
				{
					for (int j = 0; j < N; j++)
					{
						laplace_s.push_back(ss[i][j]);
					}
				}

				task->configure.LaplaceVFM_s = laplace_s;
			}
		}

		::std::vector<::std::vector<double>> initialNodalForce = span2d_to_vector2d(task->timenodalforce[start_index]);
		::std::vector<::std::vector<::std::vector<double>>> nodalforceArray = span3d_to_vector3d(task->timenodalforce.begin() + start_index, task->timenodalforce.begin() + end_index+1);

		::std::vector<::std::vector<::std::vector<double>>> velocityArray(displacementArray);//(task->timedisplacement.begin() + start_index, task->timedisplacement.end());
		::std::vector<::std::vector<::std::vector<double>>> accelerationArray(displacementArray);//(task->timedisplacement.begin() + start_index, task->timedisplacement.end());

		//task->configure.timestep.clear();
		//task->timedisplacement.clear();
		//task->timestress.clear();
		//task->timenodalforce.clear();
		//task->timeconstraintpressure.clear();
		//task->timeconstraintactivate.clear();

		for (int index_velocityArray = 0; index_velocityArray < velocityArray.size(); index_velocityArray++)
		{
#pragma omp parallel for
			for (int index_velocity = 0; index_velocity < velocityArray[index_velocityArray].size(); index_velocity++)
			{
				for (int index_coordinate = 0; index_coordinate < velocityArray[index_velocityArray][index_velocity].size(); index_coordinate++)
				{
					if (index_velocityArray == velocityArray.size() - 1)
					{
						velocityArray[index_velocityArray][index_velocity][index_coordinate] = (displacementArray[0][index_velocity][index_coordinate] - displacementArray[index_velocityArray][index_velocity][index_coordinate]) / (timeArray[index_velocityArray] - timeArray[index_velocityArray - 1]);
					}
					else
					{
						velocityArray[index_velocityArray][index_velocity][index_coordinate] = (displacementArray[index_velocityArray + 1][index_velocity][index_coordinate] - displacementArray[index_velocityArray][index_velocity][index_coordinate]) / (timeArray[index_velocityArray + 1] - timeArray[index_velocityArray]);
					}

				}
			}
		}

		for (int index_accelerationArray = 0; index_accelerationArray < accelerationArray.size(); index_accelerationArray++)
		{
#pragma omp parallel for
			for (int index_acceleration = 0; index_acceleration < accelerationArray[index_accelerationArray].size(); index_acceleration++)
			{
				for (int index_coordinate = 0; index_coordinate < accelerationArray[index_accelerationArray][index_acceleration].size(); index_coordinate++)
				{
					if (index_accelerationArray == accelerationArray.size() - 1)
					{
						accelerationArray[index_accelerationArray][index_acceleration][index_coordinate] = (velocityArray[0][index_acceleration][index_coordinate] - velocityArray[index_accelerationArray][index_acceleration][index_coordinate]) / (timeArray[index_accelerationArray] - timeArray[index_accelerationArray - 1]);
					}
					else
					{
						accelerationArray[index_accelerationArray][index_acceleration][index_coordinate] = (velocityArray[index_accelerationArray + 1][index_acceleration][index_coordinate] - velocityArray[index_accelerationArray][index_acceleration][index_coordinate]) / (timeArray[index_accelerationArray + 1] - timeArray[index_accelerationArray]);
					}
				}
			}
		}

		// task->initialCoordinate = task->initialCoordinate + initialDisp
#pragma omp parallel for
		for (int index_initialCoordinate = 0; index_initialCoordinate < task->configure.initialCoordinate.size(); index_initialCoordinate++)
		{
			for (int index_displacement = 0; index_displacement < task->configure.initialCoordinate[index_initialCoordinate].size(); index_displacement++)
			{
				task->configure.initialCoordinate[index_initialCoordinate][index_displacement] += initialDisp[index_initialCoordinate][index_displacement];
			}
		}

		for (int index_displacementArray = 0; index_displacementArray < displacementArray.size(); index_displacementArray++)
		{
#pragma omp parallel for
			for (int index_displacement = 0; index_displacement < displacementArray[index_displacementArray].size(); index_displacement++)
			{
				for (int index_coordinate = 0; index_coordinate < displacementArray[index_displacementArray][index_displacement].size(); index_coordinate++)
				{
					displacementArray[index_displacementArray][index_displacement][index_coordinate] -= initialDisp[index_displacement][index_coordinate];
				}
			}
		}

		// === Tecplot 导出：节点绝对坐标 + 位移向量（数据源：task->timedisplacement）===
		{
			std::filesystem::create_directories("./temp/debug/tecplot");

			const std::string tecpath = "./temp/debug/tecplot/nodal_pos_disp.dat";

            const bool ok = write_tecplot_nodal_position_and_displacement(
				tecpath,
				task->configure.timestep,
				task->configure.initialCoordinate,
				task->timedisplacement,
				solution_node_id_list,
                mesh,
				task->solution_elementsID,
				start_index,
				end_index);

			if (!ok)
			{
				write_to_log_2(fem, "Failed to write Tecplot file: " + tecpath + "\n", outFile);
			}
			else
			{
				write_to_log_2(fem, "Tecplot written: " + tecpath + "\n", outFile);
			}
		}

#ifdef DEBUG
		{
			int j = 0;
			FEElement& element = *(mesh.Element(j));
			int nodes_number = element.Nodes();


			::std::vector<::std::vector<::std::vector<double>>> node_xyzs_all(timeArray.size());
			for (int index_timestep = 0; index_timestep < timeArray.size(); index_timestep++)
			{
				::std::vector<::std::vector<double>> node_xyzs_r(nodes_number);
				for (int k = 0; k < nodes_number; k++)
				{
					int node_id = element.m_node[k];
					auto& node = mesh.Node(node_id);
					node_xyzs_r[k] = { task->configure.initialCoordinate[node_id][0] + displacementArray[index_timestep][node_id][0], task->configure.initialCoordinate[node_id][1] + displacementArray[index_timestep][node_id][1], task->configure.initialCoordinate[node_id][2] + displacementArray[index_timestep][node_id][2] };
				}
				node_xyzs_all[index_timestep] = node_xyzs_r;
			}
		
			pyfunctioncall(fem, "", outFile, task->pyfile_module, "Display_3D_Element_animation", node_xyzs_all, "./temp/debug/coordinate/Element " + ::std::to_string(j) + ".gif");
		}


#endif // DEBUG


		task->nVirtualFields = 4;

#pragma region getAllVirtualFields

		::std::string logs_vf = "Start get all Virtual Fields.\n";
		write_log(fem, 0, logs_vf.c_str());
		outFile << logs_vf;




#pragma region initialDeformationGradient
		// set initial displacement to nodes
		write_to_log_2(fem, "Set initial displacement to nodes.\n", outFile);

/*    */#pragma omp parallel for
		for (int index_node_id = 0; index_node_id < solution_node_id_list.size(); index_node_id++)
		{
			int node_id = solution_node_id_list[index_node_id];
			auto& node = mesh.Node(node_id);
			node.m_r0 = { task->configure.initialCoordinate[node_id][0], task->configure.initialCoordinate[node_id][1], task->configure.initialCoordinate[node_id][2] };
			node.m_rt = { task->configure.initialCoordinate[node_id][0], task->configure.initialCoordinate[node_id][1], task->configure.initialCoordinate[node_id][2] };
			node.m_rp = { 0,0,0 };
			node.m_d0 = { 0,0,0 };
			node.m_dt = { 0,0,0 };
			node.m_dp = { 0,0,0 };
		}

/*    */#pragma omp parallel for
		for (int index_seId = 0; index_seId < task->solution_elementsID.size(); index_seId++)
		{
			int j = task->solution_elementsID[index_seId];

			FEElement& element = *(mesh.Element(j));

			FESolidDomain& domain = static_cast<FESolidDomain&>(mesh.Domain(solution_elementsDomainID[index_seId]));

			try
			{
				domain_init(mesh, static_cast<FESolidElement&>(element), domain);
			}
			catch (NegativeJacobian& e)
			{
				throw e;
			}

		}

		trueJArray = ::std::vector<::std::vector<::std::vector<double>>>(timeArray.size());
		truedeformationGradientArray = ::std::vector<::std::vector<::std::vector<mat3d>>>(timeArray.size());

		for (int index_timestep = 0; index_timestep < timeArray.size(); index_timestep++)
		{
			::std::vector<::std::vector<double>> currentCoordinate(task->configure.initialCoordinate);
/*        */#pragma omp parallel for
			for (int index_coordinate_i = 0; index_coordinate_i < currentCoordinate.size(); index_coordinate_i++)
			{
				for (int index_coordinate_j = 0; index_coordinate_j < currentCoordinate[index_coordinate_i].size(); index_coordinate_j++)
				{
					currentCoordinate[index_coordinate_i][index_coordinate_j] += displacementArray[index_timestep][index_coordinate_i][index_coordinate_j];
				}
			}

/*        */#pragma omp parallel for
			for (int index_node_id = 0; index_node_id < solution_node_id_list.size(); index_node_id++)
			{
				int node_id = solution_node_id_list[index_node_id];
				auto& node = mesh.Node(node_id);
				node.m_r0 = { task->configure.initialCoordinate[node_id][0], task->configure.initialCoordinate[node_id][1], task->configure.initialCoordinate[node_id][2] }; // true displacement as initial displacement.
				node.m_d0 = { 0,0,0 };
				node.m_rt = { currentCoordinate[node_id][0], currentCoordinate[node_id][1], currentCoordinate[node_id][2] }; // move to virtual displacement
				node.m_dt = { 0,0,0 };
				node.m_rp = { 0,0,0 };
				node.m_dp = { 0,0,0 };
			}

/*        */#pragma omp parallel for
			for (int index_seId = 0; index_seId < task->solution_elementsID.size(); index_seId++)
			{
				int j = task->solution_elementsID[index_seId];

				FEElement& element = *(mesh.Element(j));

				FESolidDomain& domain = static_cast<FESolidDomain&>(mesh.Domain(solution_elementsDomainID[index_seId]));

				domain_init(mesh, static_cast<FESolidElement&>(element), domain);
			}

			trueJArray[index_timestep] = ::std::vector<::std::vector<double>>(task->solution_elementsID.size());
			truedeformationGradientArray[index_timestep] = ::std::vector<::std::vector<mat3d>>(task->solution_elementsID.size());

			/*        *///#pragma omp parallel for
			for (int index_seId = 0; index_seId < task->solution_elementsID.size(); index_seId++)
			{
				int j = task->solution_elementsID[index_seId];

				FEElement& element = *(mesh.Element(j));
				int nodes_number = element.Nodes();

				::std::vector<::std::vector<double>> node_xyzs_r0(nodes_number);
				::std::vector<::std::vector<double>> node_xyzs_rt(nodes_number);
				for (int k = 0; k < nodes_number; k++)
				{
					int node_id = element.m_node[k];
					auto& node = mesh.Node(node_id);
					::std::vector<double> xyz_r0 = { node.m_r0.x , node.m_r0.y,node.m_r0.z };
					node_xyzs_r0[k] = xyz_r0;
					::std::vector<double> xyz_rt = { node.m_rt.x , node.m_rt.y,node.m_rt.z };
					node_xyzs_rt[k] = xyz_rt;
				}
				::std::vector<::std::vector<::std::vector<double>>> node_xyzs_all{ node_xyzs_r0,node_xyzs_rt };

				// convert element to solid element
				FESolidElement& solidElement = static_cast<FESolidElement&>(element);

				// get domain
				FESolidDomain& domain = static_cast<FESolidDomain&>(mesh.Domain(solution_elementsDomainID[index_seId]));

				trueJArray[index_timestep][index_seId] = ::std::vector<double>(solidElement.GaussPoints());
				truedeformationGradientArray[index_timestep][index_seId] = ::std::vector<mat3d>(solidElement.GaussPoints());
				for (int n = 0; n < solidElement.GaussPoints(); n++)
				{
					// get element's stress
					FEMaterialPoint& mp = *(solidElement.GetMaterialPoint(n));
					FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

					double Jc;
					mat3d Ft, Fp;

					// need recalculate
					// calculate inverse jacobian
					try
					{
						Jc = domain_defgrad_GJ(domain, solidElement, Ft, n);
						//domain.defgradp(solidElement, Fp, n);
					}
					catch (NegativeJacobian& e)
					{
						double Jt = Ft.det();
						if (::std::abs(Jt) < 1e-10)
						{

						}
						else
						{
#ifdef DEBUG
							task->pyfile_module.attr("Display_3D_Element_animation")(node_xyzs_all, "./temp/debug/negative/Timestamp " + ::std::to_string(index_timestep) + " NegativeJacobian Element " + ::std::to_string(j) + ".gif");

							// "j" write to log2
							write_to_log_2(fem, "NegativeJacobian Element:" + ::std::to_string(j) + "\n", outFile);
							// write node_xyzs_r0 to log2
							write_to_log_2(fem, "node_xyzs_r0:\n", outFile);
							for (int k = 0; k < node_xyzs_r0.size(); k++)
							{
								write_to_log_2(fem, "Node id:" + ::std::to_string(element.m_node[k]) + " : " + ::std::to_string(node_xyzs_r0[k][0]) + " " + ::std::to_string(node_xyzs_r0[k][1]) + " " + ::std::to_string(node_xyzs_r0[k][2]) + "\n", outFile);
							}
							// write node_xyzs_rt to log2
							write_to_log_2(fem, "node_xyzs_rt:\n", outFile);
							for (int k = 0; k < node_xyzs_rt.size(); k++)
							{
								write_to_log_2(fem, "Node id:" + ::std::to_string(element.m_node[k]) + " : " + ::std::to_string(node_xyzs_rt[k][0]) + " " + ::std::to_string(node_xyzs_rt[k][1]) + " " + ::std::to_string(node_xyzs_rt[k][2]) + "\n", outFile);
							}

							// DEBUG, output and continue...
							// shouldn't break
							// break;

#else


							// Excessive deformation, mesh distortion!!!!
							write_to_log_2(fem, "Shouldn't continue!!! Excessive deformation, mesh distortion at Element:" + ::std::to_string(j) + ",timestep:" + ::std::to_string(index_timestep) + "\n", outFile);
							throw e; // don't continue execution!!!!


							// For debug, output and continue...
							//task->pyfile_module.attr("Display_3D_Element_animation")(node_xyzs_all, "./temp/debug/negative/Timestamp " + ::std::to_string(index_timestep) + " NegativeJacobian Element " + ::std::to_string(j) + ".gif");
#endif

						}

					}

					if (j == 35072)
					{
						int awer23dsfcsafwerf = 0;
					}

					pt.m_F = Ft;
					pt.m_J = Ft.det();

					trueJArray[index_timestep][index_seId][n] = Jc;
					truedeformationGradientArray[index_timestep][index_seId][n] = Ft;
				}


			}
		}


#pragma endregion

		auto& task_fixednode = task->configure.fixednode;

		// py call SetVirtualDisplacementFunction
		auto vfm_configure_from_py = task->pyfile_module.attr("SetVirtualDisplacementFunction")(task->configure).cast<VFMTask_configure>();
		// 从python文件中获取配置
		task->configure = vfm_configure_from_py;
		// 将虚位移计算函数赋值给vf_u_functions
		vf_u_functions = task->configure.vf_u_functions;

#pragma endregion


#pragma region deformationGradient

		// VirtualWork
		externalVirtualWork = ::std::vector<::std::vector<double>>(timeArray.size());

		true_internalVirtualWork = ::std::vector<::std::vector<double>>(timeArray.size());

		volumeVirtualWork = ::std::vector<::std::vector<double>>(timeArray.size());

		write_to_log_2(fem, "Attemp to calculate virtual work.\n", outFile);

		virtualstrainArrayV = ::std::vector<::std::vector<::std::vector<::std::vector<mat3ds>>>>(timeArray.size()); // virtual strain : timestep, vf, element

		// to output node virtural displacement
		::std::vector<::std::vector<::std::vector<::std::vector<double>>>> node_disp_array(timeArray.size(), ::std::vector<::std::vector<::std::vector<double>>>(vf_u_functions.size()));
		// to output surfaceelement node virtual displacement
		::std::vector<::std::vector<::std::vector<::std::vector<double>>>> surfaceelement_node_disp_array(timeArray.size(), ::std::vector<::std::vector<::std::vector<double>>>(vf_u_functions.size()));
		// to output surfaceelement pressure
		::std::vector<::std::vector<double>> surfaceelement_pressure(timeArray.size());
		// to output surfaceelement area
		::std::vector<::std::vector<double>> surfaceelement_area(timeArray.size(), ::std::vector<double>(vf_u_functions.size()));
		// to output srurfaceelement surfae virtual displacement
		::std::vector<::std::vector<::std::vector<double>>> surfaceelement_virtualstrain(timeArray.size(), ::std::vector<::std::vector<double>>(vf_u_functions.size()));
		// set element id to output
		int surfaceelement_output_element_id = 35072;


		// foreach timestep
		for (int index_timestep = 0; index_timestep < timeArray.size(); index_timestep++)
		{
			write_to_log_2(fem, "timestep:" + ::std::to_string(index_timestep) + "\n", outFile);

			::std::vector<::std::vector<double>> currentCoordinate(task->configure.initialCoordinate);
/*        */#pragma omp parallel for
			for (int index_coordinate_i = 0; index_coordinate_i < currentCoordinate.size(); index_coordinate_i++)
			{
				for (int index_coordinate_j = 0; index_coordinate_j < currentCoordinate[index_coordinate_i].size(); index_coordinate_j++)
				{
					currentCoordinate[index_coordinate_i][index_coordinate_j] += displacementArray[index_timestep][index_coordinate_i][index_coordinate_j];
				}
			}

			// true stress for parameters
			write_to_log_2(fem, "get true stress.\n", outFile);

			::std::vector<::std::vector<double>>& currentStress = stressArray[index_timestep];

			::std::vector<double> ivw(vf_u_functions.size()); // internalVirtualWork
			::std::vector<double> evw(vf_u_functions.size()); // externalVirtualWork



			// virtual field
			write_to_log_2(fem, "Cal virtual field.\n", outFile);

			virtualstrainArrayV[index_timestep] = ::std::vector<::std::vector<::std::vector<mat3ds>>>(vf_u_functions.size());
			externalVirtualWork[index_timestep] = ::std::vector<double>(vf_u_functions.size());

			true_internalVirtualWork[index_timestep] = ::std::vector<double>(vf_u_functions.size());

			volumeVirtualWork[index_timestep] = ::std::vector<double>(vf_u_functions.size());
			for (int index_vf = 0; index_vf < vf_u_functions.size(); index_vf++)
			{
				write_to_log_2(fem, "virtual field " + ::std::to_string(index_vf) + "\n", outFile);

				write_to_log_2(fem, "set virtual displacement.\n", outFile);
				pybind11::gil_scoped_release release;

/*            */#pragma omp parallel for
				for (int index_node_id = 0; index_node_id < solution_node_id_list.size(); index_node_id++)
				{
					int node_id = solution_node_id_list[index_node_id];
					auto& node = mesh.Node(node_id);

					pybind11::gil_scoped_acquire acquire;

					::std::vector<double> vf_u;
					if (task->configure.isLaplaceVFM == false)
					{
						vf_u = vf_u_functions[index_vf](currentCoordinate[node_id], node_id);
					}
					else
					{
						vf_u = vf_u_functions[index_vf](task->configure.initialCoordinate[node_id], node_id);
					}

					node.m_r0 = { currentCoordinate[node_id][0], currentCoordinate[node_id][1], currentCoordinate[node_id][2] }; // true displacement as initial displacement.
					node.m_d0 = { 0,0,0 };
					node.m_rt = { currentCoordinate[node_id][0]+vf_u[0],currentCoordinate[node_id][1]+vf_u[1], currentCoordinate[node_id][2]+vf_u[2] }; // move to virtual displacement, only displacement
					node.m_dt = { 0,0,0 };
					node.m_rp = { 0,0,0 };
					node.m_dp = { 0,0,0 };
				}

				pybind11::gil_scoped_acquire acquire;
				/*            *///#pragma omp parallel for
				for (int index_seId = 0; index_seId < task->solution_elementsID.size(); index_seId++)
				{
					int j = task->solution_elementsID[index_seId];

					FEElement& element = *(mesh.Element(j));

					FESolidDomain& domain = static_cast<FESolidDomain&>(mesh.Domain(solution_elementsDomainID[index_seId]));

					try
					{
						domain_init(mesh, static_cast<FESolidElement&>(element), domain);
					}
					catch (NegativeJacobian& e)
					{
						// Warning! virtual displacement is too large, mesh distortion!!!!
#ifdef DEBUG
						int nodes_number = element.Nodes();

						::std::vector<::std::vector<double>> node_xyzs_r0(nodes_number);
						::std::vector<::std::vector<double>> node_xyzs_rt(nodes_number);
						for (int k = 0; k < nodes_number; k++)
						{
							int node_id = element.m_node[k];
							auto& node = mesh.Node(node_id);
							::std::vector<double> xyz_r0 = { node.m_r0.x , node.m_r0.y,node.m_r0.z };
							node_xyzs_r0[k] = xyz_r0;
							::std::vector<double> xyz_rt = { node.m_rt.x , node.m_rt.y,node.m_rt.z };
							node_xyzs_rt[k] = xyz_rt;
						}
						::std::vector<::std::vector<::std::vector<double>>> node_xyzs_all{ node_xyzs_r0,node_xyzs_rt };
						for (int k = 0; k < node_xyzs_r0.size(); k++)
						{
							node_xyzs_all[1][k][0] = node_xyzs_r0[k][0] + node_xyzs_rt[k][0];
							node_xyzs_all[1][k][1] = node_xyzs_r0[k][1] + node_xyzs_rt[k][1];
							node_xyzs_all[1][k][2] = node_xyzs_r0[k][2] + node_xyzs_rt[k][2];
						}

						task->pyfile_module.attr("Display_3D_Element_animation")(node_xyzs_all, "./temp/debug/virtual negative/Timestamp " + ::std::to_string(index_timestep) + " Virtual NegativeJacobian Element " + ::std::to_string(j) + ".gif");

						// "j" write to log2
						write_to_log_2(fem, "NegativeJacobian Element:" + ::std::to_string(j) + "\n", outFile);
						// write node_xyzs_r0 to log2
						write_to_log_2(fem, "node_xyzs_r0:\n", outFile);
						for (int k = 0; k < node_xyzs_r0.size(); k++)
						{
							write_to_log_2(fem, ::std::to_string(node_xyzs_r0[k][0]) + " " + ::std::to_string(node_xyzs_r0[k][1]) + " " + ::std::to_string(node_xyzs_r0[k][2]) + "\n", outFile);
						}
						// write node_xyzs_rt to log2
						write_to_log_2(fem, "node_xyzs_rt:\n", outFile);
						for (int k = 0; k < node_xyzs_rt.size(); k++)
						{
							write_to_log_2(fem, ::std::to_string(node_xyzs_r0[k][0] + node_xyzs_rt[k][0]) + " " + ::std::to_string(node_xyzs_r0[k][1] + node_xyzs_rt[k][1]) + " " + ::std::to_string(node_xyzs_r0[k][2] + node_xyzs_rt[k][2]) + "\n", outFile);
						}

						// DEBUG, output and continue...
						// shouldn't break
						// break;

#else
							// Excessive deformation, mesh distortion!!!!
						throw e; // don't continue execution!!!!
#endif
					}
				}

				// internal virtual work
				write_to_log_2(fem, "calculate internal virtual work.\n", outFile);



				virtualstrainArrayV[index_timestep][index_vf] = ::std::vector<::std::vector<mat3ds>>(task->solution_elementsID.size());

/*            */#pragma omp parallel for
				for (int index_seId = 0; index_seId < task->solution_elementsID.size(); index_seId++)
				{
					int j = task->solution_elementsID[index_seId];
					FEElement& element = *(mesh.Element(j));
					FESolidElement& solidElement = static_cast<FESolidElement&>(element);
					virtualstrainArrayV[index_timestep][index_vf][index_seId] = ::std::vector<mat3ds>(solidElement.GaussPoints());
				}
				const bool CONST_calvvw_ma = true;
				if constexpr (CONST_calvvw_ma)
				{
					::std::vector<::std::vector<double>> vvw_elements(task->solution_elementsID.size());

					pybind11::gil_scoped_release release;

/*                *///#pragma omp parallel for
					for (int index_seId = 0; index_seId < task->solution_elementsID.size(); index_seId++)
					{
						int j = task->solution_elementsID[index_seId];

						FEElement& element = *(mesh.Element(j));
						int nodes_number = element.Nodes();

						const int NELN = FEElement::MAX_NODES;
						vec3d a0[NELN], u0[NELN];

						::std::vector<::std::reference_wrapper<FENode>> nodes;
						for (size_t i = 0; i < element.Nodes(); i++)
						{
							nodes.push_back(mesh.Node(element.m_node[i]));
							auto av = accelerationArray[index_timestep][element.m_node[i]];
							a0[i] = { av[0],av[1],av[2] };
							u0[i] = nodes[i].get().m_rt;
						}

						vvw_elements[index_seId] = ::std::vector<double>(element.GaussPoints());

						// find material by element
						FEMaterial* pmat = fem->GetMaterial(element.GetMatID());

						// convert element to solid element
						FESolidElement& solidElement = static_cast<FESolidElement&>(element);

						// get domain
						FESolidDomain& domain = static_cast<FESolidDomain&>(mesh.Domain(solution_elementsDomainID[index_seId]));


						for (int n = 0; n < solidElement.GaussPoints(); n++)
						{
							// get element's stress
							FEMaterialPoint& mp = *(solidElement.GetMaterialPoint(n));
							FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

							vec3d r0[FEElement::MAX_NODES];
							domain.GetReferenceNodalCoordinates(solidElement, r0);
							vec3d r[FEElement::MAX_NODES];
							domain.GetCurrentNodalCoordinates(solidElement, r);

							double Jc;
							mat3d Ft, Fp;

							// need recalculate
							// calculate inverse jacobian
							double currentJc = trueJArray[index_timestep][index_seId][n];

							try
							{
								Jc = domain_defgrad_GJ(domain, solidElement, Ft, n);
								//domain.defgradp(solidElement, Fp, n);
							}
							catch (NegativeJacobian& e)
							{
								double Jt = Ft.det();
								if (::std::abs(Jt) < 1e-10)
								{

								}
								else
								{
									// convert r0 to vector<vec3d>
									::std::vector<::std::vector<double>> copy_r0(nodes.size());
									::std::vector<::std::vector<double>> copy_r(nodes.size());
									for (int k = 0; k < nodes.size(); k++)
									{
										copy_r0[k] = { r0[k].x, r0[k].y, r0[k].z };
										copy_r[k] = { r[k].x, r[k].y, r[k].z };
									}

									::std::vector<::std::vector<::std::vector<double>>> copy_r_all{ copy_r0, copy_r };
									pybind11::gil_scoped_acquire acquire;
									task->pyfile_module.attr("Display_3D_Element_animation")(copy_r_all, "./temp/debug/virtual strain negative/Timestamp " + ::std::to_string(index_timestep) + " Virtual Strain Exception Jacobian Element " + ::std::to_string(j) + ".gif");

									throw e; // don't continue execution!!!!
								}

							}



							if (j == 35072)
							{
								int awer23dsfcsafwerf = 0;
							}



							double J = pt.m_J;
							double J0 = mp.m_J0;
							
							pt.m_F = Ft;
							pt.m_J = Ft.det();


							// strainCompute infinitesimal strain
							mat3ds infstrain = pt.m_F.sym();
							// strainCompute green-lagrange strain
							mat3ds FEstrain = pt.Strain(); // virtual strain error m_r no real position
							mat3d I = mat3d(1,0,0,0,1,0,0,0,1);
							mat3d lstrain = (pt.m_F.transpose()*pt.m_F - I)*0.5;
							// convert lstrain to mat3ds
							mat3ds strain = mat3ds(lstrain[0][0], lstrain[1][1], lstrain[2][2], lstrain[0][1], lstrain[1][2], lstrain[0][2]);
							
							mat3ds& s = FEstrain;


							if (::std::abs(pt.m_J) > 1e-10 && pt.m_J<0)
							{
								int pweoido2ke233eo = 0;

							}

							if (j == 0)
							{
								int psdoifwejmfkjwe3opewk = 0;
							}

							virtualstrainArrayV[index_timestep][index_vf][index_seId][n] = s;

							mat3ds stress_true(stressArray[index_timestep][j][0] - stressArray[0][j][0], stressArray[index_timestep][j][1] - stressArray[0][j][1], stressArray[index_timestep][j][2] - stressArray[0][j][2], stressArray[index_timestep][j][3] - stressArray[0][j][3], stressArray[index_timestep][j][4] - stressArray[0][j][4], stressArray[index_timestep][j][5] - stressArray[0][j][5]);

							true_internalVirtualWork[index_timestep][index_vf] += stress_true.dotdot(s) * currentJc;

							double density = dynamic_cast<FESolidMaterial*>(pmat)->Density(mp);
							density = 1;
							// acceleration
							vec3d a = solidElement.Evaluate(a0, n);
							vec3d u = solidElement.Evaluate(u0, n);

							double vvw = density * (a * u) * Jc;

							vvw_elements[index_seId][n] = vvw;
						}


					}

					pybind11::gil_scoped_acquire acquire;

					// sum vvw_elements
					for (int i = 0; i < vvw_elements.size(); i++)
					{
						for (int j = 0; j < vvw_elements[i].size(); j++)
						{
							volumeVirtualWork[index_timestep][index_vf] += vvw_elements[i][j];
						}
					}
				}
				// external virtual work
				write_to_log_2(fem, "calculate external virtual work.\n", outFile);

				externalVirtualWork[index_timestep][index_vf] = 0;
				// surface load
				write_to_log_2(fem, "surface load and evw:\n", outFile);

				for (int j = 0; j < task->configure.pressure_load.size(); j++)
				{
					for (int i = 0; i < fem->ModelLoads(); i++)
					{
						FEModelLoad& load = *(fem->ModelLoad(i));
						::std::string loadclassname = load.GetFactoryClass()->GetClassName();
						::std::string loadname = load.GetName();

						if (::std::get<0>(task->configure.pressure_load[j]) == loadclassname && ::std::get<1>(task->configure.pressure_load[j]) == loadname)
						{
							if (loadclassname == "FEPressureLoad")
							{
								FEPressureLoad& pressureLoad = static_cast<FEPressureLoad&>(load);

								auto& surface = pressureLoad.GetSurface();

								::std::vector<::std::vector<double>> evw_surface(surface.Elements());
								::std::vector<bool> evw_surface_flag(surface.Elements(), false);

/*							  *///#pragma omp parallel for
								for (int index_surface_element = 0; index_surface_element < surface.Elements(); index_surface_element++)
								{
									auto& element = surface.Element(index_surface_element);

									auto temp_true_element = element.m_elem[0];
									FEElement* true_element = get_FEElement_p_version(temp_true_element);

									int element_index = mesh.FindElementIndexFromID(true_element->GetID());
									::std::vector<::std::reference_wrapper<FENode>> nodes;
									for (size_t i = 0; i < element.Nodes(); i++)
									{
										nodes.push_back(mesh.Node(element.m_node[i]));
									}

									evw_surface[index_surface_element] = ::std::vector<double>(element.GaussPoints(), 0);

									// judge surfaceelement's element is in solution_elementsID
									int index_seId = -1;
									if (::std::find(task->solution_elementsID.begin(), task->solution_elementsID.end(), element_index) != task->solution_elementsID.end())
									{
										vec3d re[FEElement::MAX_NODES];
										surface.GetReferenceNodalCoordinates(element, re);
										vec3d rv_[FEElement::MAX_NODES];
										surface.NodalCoordinates(element, rv_);

										vec3d rv[FEElement::MAX_NODES];
										for (int i = 0; i < element.Nodes(); ++i)
										{
											rv[i] = rv_[i] - re[i];										
										}



										double Hr[FEElement::MAX_NODES], Hs[FEElement::MAX_NODES];
										element.shape_deriv(Hr, Hs, 0, 0);
										vec3d xr(0, 0, 0), xs(0, 0, 0);
										for (int i = 0; i < element.Nodes(); ++i)
										{
											xr += re[i] * Hr[i];
											xs += re[i] * Hs[i];
										}
										vec3d np = xr ^ xs;
										np.unit();

										Eigen::Vector3d np_eigen(np.x, np.y, np.z);

										mat3ds e_stress(currentStress[element_index][0]-initialStress[element_index][0], currentStress[element_index][1] - initialStress[element_index][1], currentStress[element_index][2] - initialStress[element_index][2], currentStress[element_index][3] - initialStress[element_index][3], currentStress[element_index][4] - initialStress[element_index][4], currentStress[element_index][5] - initialStress[element_index][5]);
										Eigen::Matrix3d e_stress_eigen=convert_mat3ds_EigenMatrix(e_stress);

										double p_e_stress = np_eigen.transpose() * e_stress_eigen * np_eigen;

										//double P = (currentStress[elementId][0] + currentStress[elementId][1] + currentStress[elementId][2]) / 3;
										//double iP = (initialStress[elementId][0] + initialStress[elementId][1] + initialStress[elementId][2]) / 3;
										double P = currentStress[element_index][0] + currentStress[element_index][1] + currentStress[element_index][2];
										double iP = initialStress[element_index][0] + initialStress[element_index][1] + initialStress[element_index][2];

										P -= iP;

										P = p_e_stress;

										// if elementId in surfaceelement_output_element_id
										if (surfaceelement_output_element_id == element_index)
										{
											surfaceelement_node_disp_array[index_timestep][index_vf] = ::std::vector<::std::vector<double>>(element.Nodes(), ::std::vector<double>(3, 0.0));
											for (int i = 0; i < element.Nodes(); ++i)
											{
												int nodeid = element.m_node[i];
												auto& node = mesh.Node(nodeid);

												surfaceelement_node_disp_array[index_timestep][index_vf][i] = ::std::vector<double>(3, 0.0);

												auto m_r0 = node.m_r0;
												auto m_rt = node.m_rt;

												auto disp = m_rt - m_r0;
												surfaceelement_node_disp_array[index_timestep][index_vf][i] = vec3d_to_vector(disp);
												
											}
											if (index_vf == 0)
											{
												surfaceelement_pressure[index_timestep] = vec3d_to_vector(np * P);
												
											}
											// 计算面积
											//if (index_vf == 0)
											{
												switch (element.Nodes()) {
													case 3:
													{
														// 三角形
														auto& nodeA = mesh.Node(element.m_node[0]);
														auto& nodeB = mesh.Node(element.m_node[1]);
														auto& nodeC = mesh.Node(element.m_node[2]);
														surfaceelement_area[index_timestep][index_vf] = CalculateTriangleArea(nodeA.m_r0, nodeB.m_r0, nodeC.m_r0);
														break;
													}

													case 4: 
													{
														// 四边形
														auto& nodeA = mesh.Node(element.m_node[0]);
														auto& nodeB = mesh.Node(element.m_node[1]);
														auto& nodeC = mesh.Node(element.m_node[2]);
														auto& nodeD = mesh.Node(element.m_node[3]);
														surfaceelement_area[index_timestep][index_vf] = CalculateQuadArea(nodeA.m_r0, nodeB.m_r0, nodeC.m_r0, nodeD.m_r0);
														break;
													}

													default: 
														// 处理不支持的单元类型或错误情况
														surfaceelement_area[index_timestep][index_vf] = 0.0; // 或抛出异常
														break;
												}
											}
										}

										for (int n = 0; n < element.GaussPoints(); n++)
										{
											FESurfaceMaterialPoint& pt = *(dynamic_cast<FESurfaceMaterialPoint*>(element.GetMaterialPoint(n)));


											// initialize some material point data
											double* H = element.H(n);
											vec3d rn(0, 0, 0);
											for (int j = 0; j < element.Nodes(); ++j)
											{
												rn += rv[j] * H[j];
											}
											pt.m_rt = rn; // initialized to zero
											//pt.m_rt = rv;

											// calculate initial surface tangents
											double* Gr = element.Gr(n);
											double* Gs = element.Gs(n);

											vec3d dxr(0, 0, 0), dxs(0, 0, 0);
											for (int i = 0; i < element.Nodes(); ++i)
											{
												dxr += re[i] * Gr[i];
												dxs += re[i] * Gs[i];
											}
											pt.dxr = dxr;
											pt.dxs = dxs;

											// initialize the other material point data
											pt.Init();


											vec3d pos = pt.m_rt;



											double J = (pt.dxr ^ pt.dxs).norm(); // real J ?
											//pt.m_Jt = J;

											//vec3d N = (pt.dxr ^ pt.dxs); N.unit();

											//vec3d Nv(0,0,0);
											//if (element.Nodes() == 4)
											//{
											//	vec3d v1 = re[1] - re[0];
											//	vec3d v2 = re[3] - re[0];

											//	vec3d Nv = (v1 ^ v2) * (-1.0 / ((v1 ^ v2).norm()));

											//	if (Nv == N)
											//	{
											//		int xx = 0;
											//	}
											//}

											vec3d t = np * P;

											int temp_debug_omp_index_vf = index_vf;
											int temp_debug_omp_index_timestep = index_timestep;


											double temp_debug_omp_evw_surface = t * pos * J;

											// correct with stress
											double weight_correct = 1;// 1.0 / (10000 * ::std::abs(P) + 0.01);

											evw_surface[index_surface_element][n] = weight_correct * temp_debug_omp_evw_surface;

											if (element_index == 35072 && n == 0)
											{
												int awer23dsfcsafwerf = 0;
											}

											if (index_timestep == 2)
											{
												int sapweoewkkememem = 0;
											}

											double ddd = 222.1;

											if (evw_surface[index_surface_element][n] > 1e-5)
											{
												auto temp_t(t);
												auto temp(P);
												auto temp_pos(pos);
												auto temp_J(J);
											}

											if (pos.x > 100)
											{
												auto temp_t(t);
												auto temp(P);
												auto temp_pos(pos);
												auto temp_J(J);
											}
										}
									}
									else
									{
										evw_surface_flag[index_surface_element] = true;
									}

								}

								// delete evw_surface_flag == true
								{
									::std::vector<::std::vector<double>> evw_surface_temp;
									for (int i = 0; i < evw_surface.size(); i++)
									{
										if (evw_surface_flag[i] == false)
										{
											evw_surface_temp.push_back(evw_surface[i]);
										}
									}
									evw_surface = evw_surface_temp;
								}

								for (int i = 0; i < evw_surface.size(); i++)
								{
									for (int j = 0; j < evw_surface[i].size(); j++)
									{
										externalVirtualWork[index_timestep][index_vf] += evw_surface[i][j];
									}
								}

								// write all evw_surface to file
								{
									::std::ofstream evw_surface_file("./temp/debug/evw_surface_" + ::std::to_string(index_timestep) + "_" + ::std::to_string(index_vf) + ".txt");
									for (int i = 0; i < evw_surface.size(); i++)
									{
										for (int j = 0; j < evw_surface[i].size(); j++)
										{
											evw_surface_file << evw_surface[i][j] << "\n";
										}
									}
									evw_surface_file.close();
								}

								auto evw = externalVirtualWork[index_timestep][index_vf];
							}
							else if (loadclassname == "FENodalForce")
							{
								FENodalLoad& fnl = dynamic_cast<FENodalLoad&>(load);

								int dofs = fnl.GetDOFList().Size();
								::std::vector<double> val(dofs, 0.0);

								auto& nset = *(fnl.GetNodeSet());

								::std::vector<double> evw_nf(nset.Size());

								node_disp_array[index_timestep][index_vf] = ::std::vector<::std::vector<double>>(nset.Size(), ::std::vector<double>(dofs, 0.0));

								for (int i = 0; i < nset.Size(); ++i)
								{
									auto& node = *nset.Node(i);
									int nodeid = node.GetID()-1;
									
									auto m_r0 = node.m_r0;
									auto m_rt = node.m_rt;

									auto disp = m_rt - m_r0;
									node_disp_array[index_timestep][index_vf][i] = vec3d_to_vector(disp);
									// if nodes in solution_node_id_list
									if (::std::find(solution_node_id_list.begin(), solution_node_id_list.end(), nodeid) != solution_node_id_list.end())
									{
										// get the nodal values
										
										vec3d nF = { nodalforceArray[index_timestep][i][0],nodalforceArray[index_timestep][i][1],nodalforceArray[index_timestep][i][2] };
										vec3d inF = { initialNodalForce[i][0],initialNodalForce[i][1],initialNodalForce[i][2] };

										nF -= inF;

										double nFW = nF * disp;

										evw_nf[i] = nFW;
									}
									else
									{
										evw_nf[i] = 0;
									}
								}

								double evw = 0;

								for (int i = 0; i < nset.Size(); i++)
								{
									evw += evw_nf[i];
								}

								externalVirtualWork[index_timestep][index_vf] += evw;
							}
						}
					}

				}

				for (int j = 0; j < task->configure.constraint_load.size(); j++)
				{
					int index_i = 0;
					for (int i = 0; i < fem->NonlinearConstraints(); i++)
					{
						FENLConstraint& fnc = *(fem->NonlinearConstraint(i));
						::std::string fncclassname = fnc.GetFactoryClass()->GetClassName();
						::std::string fncname = fnc.GetName();

						write_to_log_2(fem, ::std::string("fncclassname: ") + fncclassname + " fncname: " + fncname + "\n", outFile);

						if (::std::get<0>(task->configure.constraint_load[j]) == fncclassname && ::std::get<1>(task->configure.constraint_load[j]) == fncname)
						{
							write_to_log_2(fem, ::std::string("activate: ") + ::std::to_string(constraintActivateArray[index_timestep][index_i]) + "\n", outFile);
							if (constraintActivateArray[index_timestep][index_i] == true)
							{
								if (fncclassname == "FEVolumeConstraint")
								{
									FEVolumeConstraint& vc = dynamic_cast<FEVolumeConstraint&>(fnc);

									FESurface& surface = *(vc.GetSurface());

									::std::vector<::std::vector<double>> evw_surface(surface.Elements());
									::std::vector<bool> evw_surface_flag(surface.Elements(), false);

									write_to_log_2(fem, ::std::string("index_surface_elements number: ") + ::std::to_string(surface.Elements()) + "\n", outFile);

									for (int index_surface_element = 0; index_surface_element < surface.Elements(); index_surface_element++)
									{
										auto& element = surface.Element(index_surface_element);

										auto temp_true_element = element.m_elem[0];
										FEElement* true_element = get_FEElement_p_version(temp_true_element);
										

										int element_index = mesh.FindElementIndexFromID(true_element->GetID());
										::std::vector<::std::reference_wrapper<FENode>> nodes;
										for (size_t i = 0; i < element.Nodes(); i++)
										{
											nodes.push_back(mesh.Node(element.m_node[i]));
										}

										evw_surface[index_surface_element] = ::std::vector<double>(element.GaussPoints(), 0);

										write_to_log_2(fem, ::std::string("index_surface_element: ") + ::std::to_string(index_surface_element) + "\n", outFile);

										// judge surfaceelement's element is in solution_elementsID
										int index_seId = -1;
										if (::std::find(task->solution_elementsID.begin(), task->solution_elementsID.end(), element_index) != task->solution_elementsID.end())
										{
											vec3d re[FEElement::MAX_NODES];
											surface.GetReferenceNodalCoordinates(element, re);
											vec3d rv_[FEElement::MAX_NODES];
											surface.NodalCoordinates(element, rv_);

											vec3d rv[FEElement::MAX_NODES];
											for (int i = 0; i < element.Nodes(); ++i)
											{
												rv[i] = rv_[i] - re[i];
											}



											double Hr[FEElement::MAX_NODES], Hs[FEElement::MAX_NODES];
											element.shape_deriv(Hr, Hs, 0, 0);
											vec3d xr(0, 0, 0), xs(0, 0, 0);
											for (int i = 0; i < element.Nodes(); ++i)
											{
												xr += re[i] * Hr[i];
												xs += re[i] * Hs[i];
											}
											vec3d np = xr ^ xs;
											np.unit();

											write_to_log_2(fem, ::std::string("np: ") + ::std::to_string(np.x) + " " + ::std::to_string(np.y) + " " + ::std::to_string(np.z) + "\n", outFile);

											write_to_log_2(fem, ::std::string("index_timestep: ") + ::std::to_string(index_timestep)+"\n", outFile);

											write_to_log_2(fem, ::std::string("constraintPressureArray size: ") + ::std::to_string(constraintPressureArray.size()) + "\n", outFile);

											double P = constraintPressureArray[index_timestep][index_i];
											double iP = initialConstraintPressure[index_i];

											P -= iP;


											write_to_log_2(fem, ::std::string("P: ") + ::std::to_string(P) + "\n", outFile);

											// if elementId in surfaceelement_output_element_id
											if (surfaceelement_output_element_id == element_index)
											{
												surfaceelement_node_disp_array[index_timestep][index_vf] = ::std::vector<::std::vector<double>>(element.Nodes(), ::std::vector<double>(3, 0.0));
												for (int i = 0; i < element.Nodes(); ++i)
												{
													int nodeid = element.m_node[i];
													auto& node = mesh.Node(nodeid);

													surfaceelement_node_disp_array[index_timestep][index_vf][i] = ::std::vector<double>(3, 0.0);

													auto m_r0 = node.m_r0;
													auto m_rt = node.m_rt;

													auto disp = m_rt - m_r0;
													surfaceelement_node_disp_array[index_timestep][index_vf][i] = vec3d_to_vector(disp);

												}
												if (index_vf == 0)
												{
													surfaceelement_pressure[index_timestep] = vec3d_to_vector(np * P);

												}
												// 计算面积
												//if (index_vf == 0)
												{
													switch (element.Nodes()) {
													case 3:
													{
														// 三角形
														auto& nodeA = mesh.Node(element.m_node[0]);
														auto& nodeB = mesh.Node(element.m_node[1]);
														auto& nodeC = mesh.Node(element.m_node[2]);
														surfaceelement_area[index_timestep][index_vf] = CalculateTriangleArea(nodeA.m_r0, nodeB.m_r0, nodeC.m_r0);
														break;
													}

													case 4:
													{
														// 四边形
														auto& nodeA = mesh.Node(element.m_node[0]);
														auto& nodeB = mesh.Node(element.m_node[1]);
														auto& nodeC = mesh.Node(element.m_node[2]);
														auto& nodeD = mesh.Node(element.m_node[3]);
														surfaceelement_area[index_timestep][index_vf] = CalculateQuadArea(nodeA.m_r0, nodeB.m_r0, nodeC.m_r0, nodeD.m_r0);
														break;
													}

													default:
														// 处理不支持的单元类型或错误情况
														surfaceelement_area[index_timestep][index_vf] = 0.0; // 或抛出异常
														break;
													}
												}
											}

											for (int n = 0; n < element.GaussPoints(); n++)
											{
												FESurfaceMaterialPoint& pt = *(dynamic_cast<FESurfaceMaterialPoint*>(element.GetMaterialPoint(n)));


												// initialize some material point data
												double* H = element.H(n);
												vec3d rn(0, 0, 0);
												for (int j = 0; j < element.Nodes(); ++j)
												{
													rn += rv[j] * H[j];
												}
												pt.m_rt = rn; // initialized to zero
												//pt.m_rt = rv;

												// calculate initial surface tangents
												double* Gr = element.Gr(n);
												double* Gs = element.Gs(n);

												vec3d dxr(0, 0, 0), dxs(0, 0, 0);
												for (int i = 0; i < element.Nodes(); ++i)
												{
													dxr += re[i] * Gr[i];
													dxs += re[i] * Gs[i];
												}
												pt.dxr = dxr;
												pt.dxs = dxs;

												// initialize the other material point data
												pt.Init();


												vec3d pos = pt.m_rt;



												double J = (pt.dxr ^ pt.dxs).norm(); // real J ?
												//pt.m_Jt = J;

												//vec3d N = (pt.dxr ^ pt.dxs); N.unit();

												//vec3d Nv(0,0,0);
												//if (element.Nodes() == 4)
												//{
												//	vec3d v1 = re[1] - re[0];
												//	vec3d v2 = re[3] - re[0];

												//	vec3d Nv = (v1 ^ v2) * (-1.0 / ((v1 ^ v2).norm()));

												//	if (Nv == N)
												//	{
												//		int xx = 0;
												//	}
												//}

												vec3d t = np * P;

												int temp_debug_omp_index_vf = index_vf;
												int temp_debug_omp_index_timestep = index_timestep;


												double temp_debug_omp_evw_surface = t * pos * J;

												// correct with stress
												double weight_correct = 1;// 1.0 / (10000 * ::std::abs(P) + 0.01);

												evw_surface[index_surface_element][n] = weight_correct * temp_debug_omp_evw_surface;

												if (element_index == 35072 && n == 0)
												{
													int awer23dsfcsafwerf = 0;
												}

												if (index_timestep == 2)
												{
													int sapweoewkkememem = 0;
												}

												double ddd = 222.1;

												if (evw_surface[index_surface_element][n] > 1e-5)
												{
													auto temp_t(t);
													auto temp(P);
													auto temp_pos(pos);
													auto temp_J(J);
												}

												if (pos.x > 100)
												{
													auto temp_t(t);
													auto temp(P);
													auto temp_pos(pos);
													auto temp_J(J);
												}
											}
										}
										else
										{
											evw_surface_flag[index_surface_element] = true;
										}

									}



									// delete evw_surface_flag == true
									{
										::std::vector<::std::vector<double>> evw_surface_temp;
										for (int i = 0; i < evw_surface.size(); i++)
										{
											if (evw_surface_flag[i] == false)
											{
												evw_surface_temp.push_back(evw_surface[i]);
											}
										}
										evw_surface = evw_surface_temp;
									}

									double evw = 0.0;

									for (int i = 0; i < evw_surface.size(); i++)
									{
										for (int j = 0; j < evw_surface[i].size(); j++)
										{
											evw += evw_surface[i][j];
										}
									}


									externalVirtualWork[index_timestep][index_vf] += evw;
								}
							}
							index_i++;
						}
					}
				}
			}
		}

#pragma region Laplace
		if (task->configure.isLaplaceVFM == true)
		{

			L_externalVirtualWork_LaplcaeTransform(externalVirtualWork, externalVirtualWork_laplace, timeArray, task);

			L_StressPK2_withJc_LaplaceTransform(timeArray, task, fem, solution_elementsDomainID, truedeformationGradientArray, trueJArray, Sepk2_dotdot_vStrain, virtualstrainArrayV);

			fs = L_fs(externalVirtualWork, timeArray, task, fem, solution_elementsDomainID, truedeformationGradientArray, trueJArray, virtualstrainArrayV);
		
			if (task->configure.isTalbotLaplaceVFM_s == true)
			{
				ft = inv_ft(fs, timeArray);
			}
			
		}
#pragma endregion

		// output nodalforceArray[index_timestep][i][0] to csv file
		if(task->configure.optim_function_output_debug_info)
		{
			::std::ofstream nodalforceArray_file;
			nodalforceArray_file.open("./temp/debug/result/nodalforceArray.csv", ::std::ios::ate | ::std::ios::out);
			for (int index_timestep = 0; index_timestep < nodalforceArray.size(); index_timestep++)
			{
				for (int i = 0; i < nodalforceArray[index_timestep].size(); i++)
				{
					nodalforceArray_file << ::std::setprecision(12) << nodalforceArray[index_timestep][i][0] << ",";
				}
				if (index_timestep != nodalforceArray.size() - 1)
				{
					nodalforceArray_file << "\n";
				}
			}
			nodalforceArray_file.close();

			// output node_disp_array[index_timestep][index_vf] to csv file
			::std::ofstream node_disp_array_file;
			node_disp_array_file.open("./temp/debug/result/node_disp_array.csv", ::std::ios::ate | ::std::ios::out);
			// each row is a timestep, only output index_vf 0, each column is a node x,y,z displacement
			for (int index_timestep = 0; index_timestep < node_disp_array.size(); index_timestep++)
			{
				for (int i = 0; i < node_disp_array[index_timestep][0].size(); i++)
				{
					node_disp_array_file << ::std::setprecision(12) << node_disp_array[index_timestep][0][i][0] << "," << node_disp_array[index_timestep][0][i][1] << "," << node_disp_array[index_timestep][0][i][2];
					if (i != node_disp_array[index_timestep][0].size() - 1)
					{
						node_disp_array_file << ",";
					}
				}
				if (index_timestep != node_disp_array.size() - 1)
				{
					node_disp_array_file << "\n";
				}
			}
			node_disp_array_file.close();

			// output virtualstrainArrayV[index_timestep][index_vf][index_seId][n] to csv file
			::std::ofstream virtualstrainArrayV_file;
			virtualstrainArrayV_file.open("./temp/debug/result/virtualstrainArrayV.csv", ::std::ios::ate | ::std::ios::out);
			for (int index_timestep = 0; index_timestep < virtualstrainArrayV.size(); index_timestep++)
			{
				for (int i = 0; i < virtualstrainArrayV[index_timestep][0].size(); i++)
				{
					// sum_strain_gauss
					mat3d sum_strain_gauss(0,0,0,0,0,0,0,0,0);
					for (int n = 0; n < virtualstrainArrayV[index_timestep][0][i].size(); n++)
					{
						sum_strain_gauss += virtualstrainArrayV[index_timestep][0][i][n] / virtualstrainArrayV[index_timestep][0][i].size();
					}
					virtualstrainArrayV_file << ::std::setprecision(12) << sum_strain_gauss[0][0] << "," << sum_strain_gauss[1][1] << "," << sum_strain_gauss[2][2] << "," << sum_strain_gauss[0][1] << "," << sum_strain_gauss[1][2] << "," << sum_strain_gauss[2][0];
					if (i != virtualstrainArrayV[index_timestep][0].size() - 1)
					{
						virtualstrainArrayV_file << ",";
					}
				}
				if (index_timestep != virtualstrainArrayV.size() - 1)
				{
					virtualstrainArrayV_file << "\n";
				}
			}
			virtualstrainArrayV_file.close();

			// output trueJArray[index_timestep][index_seId] to csv file
			::std::ofstream trueJArray_file;
			trueJArray_file.open("./temp/debug/result/trueJArray.csv", ::std::ios::ate | ::std::ios::out);
			for (int index_timestep = 0; index_timestep < trueJArray.size(); index_timestep++)
			{
				for (int i = 0; i < trueJArray[index_timestep].size(); i++)
				{
					for (int j = 0; j < trueJArray[index_timestep][i].size(); j++)
					{
						trueJArray_file << ::std::setprecision(12) << trueJArray[index_timestep][i][j];
						
						if (i != trueJArray[index_timestep].size() - 1 || j != trueJArray[index_timestep][i].size() - 1)
						{
							trueJArray_file << ",";
						}
					}
				}
				if (index_timestep != trueJArray.size() - 1)
				{
					trueJArray_file << "\n";
				}
			}
			trueJArray_file.close();

			// output surfaceelement_node_disp_array to csv file
			::std::ofstream surfaceelement_node_disp_array_file;
			surfaceelement_node_disp_array_file.open("./temp/debug/result/surface_elementIDselected/surfaceelement_node_disp_array.csv", ::std::ios::ate | ::std::ios::out);
			for (int index_timestep = 0; index_timestep < surfaceelement_node_disp_array.size(); index_timestep++)
			{
				for (int i = 0; i < surfaceelement_node_disp_array[index_timestep].size(); i++)
				{
					for (int j = 0; j < surfaceelement_node_disp_array[index_timestep][i].size(); j++)
					{
						for (int k = 0; k < surfaceelement_node_disp_array[index_timestep][i][j].size(); k++)
						{
							surfaceelement_node_disp_array_file << ::std::setprecision(12) << surfaceelement_node_disp_array[index_timestep][i][j][k];
							if (k != surfaceelement_node_disp_array[index_timestep][i][j].size() - 1 || j != surfaceelement_node_disp_array[index_timestep][i].size() - 1 || i != surfaceelement_node_disp_array[index_timestep].size() - 1)
							{
								surfaceelement_node_disp_array_file << ",";
							}
						}
					}
				}
				if (index_timestep != surfaceelement_node_disp_array.size() - 1)
				{
					surfaceelement_node_disp_array_file << "\n";
				}
			}
			surfaceelement_node_disp_array_file.close();

			// output surfaceelement_pressure to csv file
			::std::ofstream surfaceelement_pressure_file;
			surfaceelement_pressure_file.open("./temp/debug/result/surface_elementIDselected/surfaceelement_pressure.csv", ::std::ios::ate | ::std::ios::out);
			for (int index_timestep = 0; index_timestep < surfaceelement_pressure.size(); index_timestep++)
			{
				for (int i = 0; i < surfaceelement_pressure[index_timestep].size(); i++)
				{
					surfaceelement_pressure_file << ::std::setprecision(12) << surfaceelement_pressure[index_timestep][i];
					if (i != surfaceelement_pressure[index_timestep].size() - 1)
					{
						surfaceelement_pressure_file << ",";
					}
					
				}
				if (index_timestep != surfaceelement_pressure.size() - 1)
				{
					surfaceelement_pressure_file << "\n";
				}
			}
			surfaceelement_pressure_file.close();

			// output surfaceelement_area to csv file
			::std::ofstream surfaceelement_area_file;
			surfaceelement_area_file.open("./temp/debug/result/surface_elementIDselected/surfaceelement_area.csv", ::std::ios::ate | ::std::ios::out);
			for (int index_timestep = 0; index_timestep < surfaceelement_area.size(); index_timestep++)
			{
				for (int i = 0; i < surfaceelement_area[index_timestep].size(); i++)
				{
					surfaceelement_area_file << ::std::setprecision(12) << surfaceelement_area[index_timestep][i];
					if (i != surfaceelement_area[index_timestep].size() - 1)
					{
						surfaceelement_area_file << ",";
					}

				}
				if (index_timestep != surfaceelement_area.size() - 1)
				{
					surfaceelement_area_file << "\n";
				}
			}
			surfaceelement_area_file.close();

			// from vritualstrainArrayV find surfaceelement_output_element_id to surfaceelement_virtualstrain

			// find surfaceelement_output_element_id's index in task->solution_elementsID
			int surfaceelement_output_element_id_solution_elementsID_index = -1;
			for (int i = 0; i < task->solution_elementsID.size(); i++)
			{
				if (task->solution_elementsID[i] == surfaceelement_output_element_id)
				{
					surfaceelement_output_element_id_solution_elementsID_index = i;
					break;
				}
			}

			// 用于存储提取的应变数据的向量
			::std::vector<::std::vector<::std::vector<mat3ds>>> extractedVirtualstrainArrayV;

			for (size_t index_timestep = 0; index_timestep < virtualstrainArrayV.size(); ++index_timestep) {
				const auto& dimension1 = virtualstrainArrayV[index_timestep];
					
				std::vector<std::vector<mat3ds>> extractedDim1;

				for (size_t index_vf = 0; index_vf < dimension1.size(); ++index_vf) {
					const auto& dimension2 = dimension1[index_vf];
					std::vector<mat3ds> extractedDim2;

					if (dimension2.size() > surfaceelement_output_element_id_solution_elementsID_index) {
						extractedDim2 = dimension2[surfaceelement_output_element_id_solution_elementsID_index];

						::std::vector<double> sumVector(6, 0.0);
						for (int i = 0; i < extractedDim2.size(); i++)
						{
							::std::vector<double> tempVector{ extractedDim2[i].xx(), extractedDim2[i].yy(),extractedDim2[i].zz(),extractedDim2[i].xy(),extractedDim2[i].yz(),extractedDim2[i].xz() };
							for (int j = 0; j < 6; j++)
							{
								sumVector[j] += tempVector[j]/6;
							}
						}
						surfaceelement_virtualstrain[index_timestep][index_vf] = sumVector;
					}

					if (!extractedDim2.empty()) {
						extractedDim1.push_back(extractedDim2);
					}
				}

				if (!extractedDim1.empty()) {
					extractedVirtualstrainArrayV.push_back(extractedDim1);
				}
			}


			for (int index_timestep = 0; index_timestep < virtualstrainArrayV.size(); index_timestep++)
			{
			}

			// output surfaceelement_virtualstrain to csv file
			::std::ofstream surfaceelement_virtualstrain_file;
			surfaceelement_virtualstrain_file.open("./temp/debug/result/surface_elementIDselected/surfaceelement_virtualstrain.csv", ::std::ios::ate | ::std::ios::out);
			for (int index_timestep = 0; index_timestep < surfaceelement_virtualstrain.size(); index_timestep++)
			{
				for (int i = 0; i < surfaceelement_virtualstrain[index_timestep].size(); i++)
				{
					for (int j = 0; j < surfaceelement_virtualstrain[index_timestep][i].size(); j++)
					{
						surfaceelement_virtualstrain_file << ::std::setprecision(12) << surfaceelement_virtualstrain[index_timestep][i][j];
						if (j != surfaceelement_virtualstrain[index_timestep][i].size() - 1 || i != surfaceelement_virtualstrain[index_timestep].size() - 1)
						{
							surfaceelement_virtualstrain_file << ",";
						}
					}
				}
				if (index_timestep != surfaceelement_virtualstrain.size() - 1)
				{
					surfaceelement_virtualstrain_file << "\n";
				}
			}
			surfaceelement_virtualstrain_file.close();

			// output stressArray to csv file
			::std::ofstream stressArray_file;
			stressArray_file.open("./temp/debug/result/FEBio_stressArray.csv", ::std::ios::ate | ::std::ios::out);
			for (int index_seId = 0; index_seId < task->solution_elementsID.size(); index_seId++)
			{
				int index_elementId = task->solution_elementsID[index_seId];
				for (int index_timestep = 0; index_timestep < stressArray.size(); index_timestep++)
				{
					auto& stress = stressArray[index_timestep][index_elementId];
					for (int i = 0; i < stress.size(); i++)
					{
						stressArray_file << ::std::setprecision(12) << stress[i];
						if (i != stress.size() - 1 || index_timestep != stressArray.size() - 1)
						{
							stressArray_file << ",";
						}
					}
				}
				if (index_seId != task->solution_elementsID.size() - 1)
				{
					stressArray_file << "\n";
				}
			}
			stressArray_file.close();

			// output externalVirtualWork to csv file
			::std::ofstream externalVirtualWork_file;
			externalVirtualWork_file.open("./temp/debug/result/true_externalVirtualWork.csv", ::std::ios::ate | ::std::ios::out);
			for (int index_timestep = 0; index_timestep < externalVirtualWork.size(); index_timestep++)
			{
				for (int index_vf = 0; index_vf < externalVirtualWork[index_timestep].size(); index_vf++)
				{
					auto& externalVirtualWork_value = externalVirtualWork[index_timestep][index_vf];
					externalVirtualWork_file << ::std::setprecision(12) << externalVirtualWork_value;
					if (index_vf != externalVirtualWork[index_timestep].size() - 1)
					{
						externalVirtualWork_file << ",";
					}
					
				}
				if (index_timestep != externalVirtualWork.size() - 1)
				{
					externalVirtualWork_file << "\n";
				}
			}
			externalVirtualWork_file.close();

			// output true_internalVirtualWork to csv file
			::std::ofstream true_internalVirtualWork_file;
			true_internalVirtualWork_file.open("./temp/debug/result/true_internalVirtualWork.csv", ::std::ios::ate | ::std::ios::out);
			for (int index_timestep = 0; index_timestep < true_internalVirtualWork.size(); index_timestep++)
			{
				for (int index_vf = 0; index_vf < true_internalVirtualWork[index_timestep].size(); index_vf++)
				{
					auto& internalVirtualWork_value = true_internalVirtualWork[index_timestep][index_vf];
					true_internalVirtualWork_file << ::std::setprecision(12) << internalVirtualWork_value;
					if (index_vf != true_internalVirtualWork[index_timestep].size() - 1)
					{
						true_internalVirtualWork_file << ",";
					}
					
				}
				if (index_timestep != true_internalVirtualWork.size() - 1)
				{
					true_internalVirtualWork_file << "\n";
				}
			}
			true_internalVirtualWork_file.close();
		}

		write_to_log_2(fem, "Debug Output csv file... \n", outFile);
#pragma endregion


		double forceMPa = 0.0004; // 3mmHg

#pragma endregion


		DumpFile dumpfile(*fem);
		::std::string fundump_path = task->dumpfile + ".fundump";
		if (task->configure.Optim_dump_path.empty() != false && task->configure.Optim_dump_path != "")
		{
			fundump_path = task->configure.Optim_dump_path;
		}
		dumpfile.Create(fundump_path.c_str());

		dumpfile& timeArray& (task->solution_elementsID)& solution_elementsDomainID& trueJArray& truedeformationGradientArray& virtualstrainArrayV& externalVirtualWork& volumeVirtualWork& externalVirtualWork_laplace& Sepk2_dotdot_vStrain& fs;

		dumpfile.Close();

		params.reset(new FunOptimParams{ timeArray, task->solution_elementsID, solution_elementsDomainID, trueJArray, truedeformationGradientArray, virtualstrainArrayV, externalVirtualWork, volumeVirtualWork, task->configure.isLaplaceVFM, task->configure.LaplaceVFM_s, externalVirtualWork_laplace, Sepk2_dotdot_vStrain, fs, task->configure.optim_function_output_debug_info });

		::std::function<double(const ::std::vector<double>&)> fun = [&fem, &outFile, &params](const ::std::vector<double>& ps) {
			double p_E = ps[0];
			double p_g = ps[1];
			double p_t = ps[2];

			double value = fun_for_optim(fem, p_g, p_t, p_E, outFile, *params);


			::std::stringstream ss;
			ss << "E: " << ::std::setprecision(12) << p_E << ", g: " << ::std::setprecision(12) << p_g << ", t: " << ::std::setprecision(12) << p_t << "\nloss: " << ::std::setprecision(12) << value << "\n";
			write_to_log_2(fem, ss.str(), outFile);

			return value;
			};

		function = fun;
	}
	else
	{
		if (task->readfromsaveOptimfunc_index.empty())
		{
			DumpFile dumpfile(*fem);
			::std::string fundump_path = task->dumpfile + ".fundump";
			if (task->configure.Optim_dump_path.empty() != false && task->configure.Optim_dump_path != "")
			{
				fundump_path = task->configure.Optim_dump_path;
			}
			dumpfile.Open(fundump_path.c_str());

			dumpfile& timeArray& (task->solution_elementsID)& solution_elementsDomainID& trueJArray& truedeformationGradientArray& virtualstrainArrayV& externalVirtualWork& volumeVirtualWork& externalVirtualWork_laplace& Sepk2_dotdot_vStrain &fs;

			dumpfile.Close();

			if (task->configure.isSetNewLaplaceVFMs == true && task->configure.isLaplaceVFM == true)
			{
				L_externalVirtualWork_LaplcaeTransform(externalVirtualWork, externalVirtualWork_laplace, timeArray, task);
				L_StressPK2_withJc_LaplaceTransform(timeArray, task, fem, solution_elementsDomainID, truedeformationGradientArray, trueJArray, Sepk2_dotdot_vStrain, virtualstrainArrayV);
				fs = L_fs(externalVirtualWork, timeArray, task, fem, solution_elementsDomainID, truedeformationGradientArray, trueJArray, virtualstrainArrayV);
				
				if (task->configure.isTalbotLaplaceVFM_s == true)
				{
					ft = inv_ft(fs, timeArray);
				}
			}

			params.reset(new FunOptimParams{ timeArray, task->solution_elementsID, solution_elementsDomainID, trueJArray, truedeformationGradientArray, virtualstrainArrayV, externalVirtualWork, volumeVirtualWork, task->configure.isLaplaceVFM, task->configure.LaplaceVFM_s, externalVirtualWork_laplace, Sepk2_dotdot_vStrain, fs, task->configure.optim_function_output_debug_info });

			::std::function<double(const ::std::vector<double>&)> fun = [&fem, &outFile, &params](const ::std::vector<double>& ps) {
				double p_E = ps[0];
				double p_g = ps[1];
				double p_t = ps[2];

				double value = fun_for_optim(fem, p_g, p_t, p_E, outFile, *params);

				::std::stringstream ss;
				ss << "E: " << ::std::setprecision(12) << p_E << ", g: " << ::std::setprecision(12) << p_g << ", t: " << ::std::setprecision(12) << p_t << "\nloss: " << ::std::setprecision(12) << value << "\n";
				write_to_log_2(fem, ss.str(), outFile);

				return value;
				};

			function = fun;
		}
		else
		{
			::std::vector<::std::function<double(const ::std::vector<double>&)>> funs;
			write_to_log_2(fem, (::std::string("Create opt fun ") + ::std::to_string(task->readfromsaveOptimfunc_index.size()) + "\n").c_str(), outFile);
			for (auto& index : task->readfromsaveOptimfunc_index)
			{
				::std::string fundumpfile = task->szfile + "_BMP" + ::std::to_string(index) + "_dump.febdump" + ".fundump";
				write_to_log_2(fem, ("Reading... " + fundumpfile + "\n").c_str(), outFile);


				DumpFile dumpfile(*fem);
				dumpfile.Open(fundumpfile.c_str());

				dumpfile& timeArray& (task->solution_elementsID)& solution_elementsDomainID& trueJArray& truedeformationGradientArray& virtualstrainArrayV& externalVirtualWork& volumeVirtualWork& externalVirtualWork_laplace& Sepk2_dotdot_vStrain &fs;

				dumpfile.Close();

				if (task->configure.isSetNewLaplaceVFMs == true && task->configure.isLaplaceVFM == true)
				{
					L_externalVirtualWork_LaplcaeTransform(externalVirtualWork, externalVirtualWork_laplace, timeArray, task);
					L_StressPK2_withJc_LaplaceTransform(timeArray, task, fem, solution_elementsDomainID, truedeformationGradientArray, trueJArray, Sepk2_dotdot_vStrain, virtualstrainArrayV);
					fs = L_fs(externalVirtualWork, timeArray, task, fem, solution_elementsDomainID, truedeformationGradientArray, trueJArray, virtualstrainArrayV);

					if (task->configure.isTalbotLaplaceVFM_s == true)
					{
						ft = inv_ft(fs, timeArray);
					}
				}

				params.reset(new FunOptimParams{ timeArray, task->solution_elementsID, solution_elementsDomainID, trueJArray, truedeformationGradientArray, virtualstrainArrayV, externalVirtualWork, volumeVirtualWork, task->configure.isLaplaceVFM, task->configure.LaplaceVFM_s, externalVirtualWork_laplace, Sepk2_dotdot_vStrain, fs, task->configure.optim_function_output_debug_info });


				funs.push_back([&fem, &outFile, &params](const ::std::vector<double>& ps) {
					double p_E = ps[0];
					double p_g = ps[1];
					double p_t = ps[2];

					double value = fun_for_optim(fem, p_g, p_t, p_E, outFile, *params);

					return value;
					});
			}
			function = [&fem, &outFile, funs](const ::std::vector<double>& ps) {
				double p_E = ps[0];
				double p_g = ps[1];
				double p_t = ps[2];

				double value = 0;
				for (auto& fun : funs)
				{
					value += fun(ps);
				}

				::std::stringstream ss;
				ss << "E: " << ::std::setprecision(12) << p_E << ", g: " << ::std::setprecision(12) << p_g << ", t: " << ::std::setprecision(12) << p_t << "\nloss: " << ::std::setprecision(12) << value << "\n";
				write_to_log_2(fem, ss.str(), outFile);
				return value;
				};
		}
	}

#pragma region Optim

	//auto xxv = function({ 0.3,8,1 });

	//auto xxxv = function({ 0.2,9,1 });

	//auto xxxxv = function({ 0.5,6,1 });

   	task->configure.optim_function = function;

	task->configure.optim_function_E__gamma_tau = [&fem, &params](double p_g, double p_t) {
		//double p_E = ps[0];
		//double p_g = ps[1];
		//double p_t = ps[2];

		auto Es = fun_optim_E__gamma_tau_Laplace(fem, p_g, p_t, *params);

		return Es;
		};

	task->configure.optim_function_gamma__E_tau = [&fem, &params](double p_E, double p_t) {
		//double p_E = ps[0];
		//double p_g = ps[1];
		//double p_t = ps[2];

		auto gs = fun_optim_gamma__E_tau_Laplace(fem, p_E, p_t, *params);

		return gs;
		};
 
    task->configure.optim_function_tau__E_gamma = [&fem, &params](double p_E, double p_g) {
		//double p_E = ps[0];
		//double p_g = ps[1];
		//double p_t = ps[2];

		auto ts = fun_optim_tau__E_gamma_Laplace(fem, p_E, p_g, *params);

		return ts;
		};

	task->configure.optim_function_E_gamma_tau = [&fem, &params](const ::std::vector<double>& ps) {
		
		double p_E = ps[0];
		double p_g = ps[1];
		double p_t = ps[2];

		double value = fun_for_optim_E_gamma_tau(fem, p_E, p_g, p_t, *params);

		return value;
		};

	task->configure.optim_function_E_gamma_tau_invF = [&fem, &params](const ::std::vector<double>& ps) {
		double p_E = ps[0];
		double p_g = ps[1];
		double p_t = ps[2];

		double value = fun_for_optim_E_gamma_tau_invF(fem, p_E, p_g, p_t, *params);

		return value;
		};

	task->configure.optim_function_elastic_E = [&fem, &params, &outFile](const ::std::vector<double>& ps) {
		double p_E = ps[0];

		double value = fun_for_optim_elastic_E(fem, p_E, outFile, *params);

		return value;
		};

	auto inter_nEvw = cal_internal_normal_linerE_vw(timeArray, task, fem, solution_elementsDomainID, truedeformationGradientArray, trueJArray, virtualstrainArrayV);

	::std::vector<::std::vector<double>> exter_nEvw(params->externalVirtualWork);
	for (int i = 0; i < exter_nEvw.size(); ++i)
	{
		for (int j = 0; j < exter_nEvw[i].size(); ++j)
		{
			exter_nEvw[i][j] = params->externalVirtualWork[i][j] - ::std::get<1>(inter_nEvw)[i][j];
		}
	}

	auto& visco_mask = ::std::get<2>(inter_nEvw);
	::std::ofstream visco_mask_file;
	visco_mask_file.open("./temp/debug/result/ElementInSolution_visco_mask.csv", ::std::ios::ate | ::std::ios::out);
	for (int index_vm=0; index_vm<visco_mask.size(); ++index_vm)
	{
		visco_mask_file << visco_mask[index_vm];
		if (index_vm != visco_mask.size() - 1)
		{
			visco_mask_file << ",";
		}
	}
	visco_mask_file.close();

	auto internal_normal_visco_dse0_strain_Jc = ::std::get<4>(inter_nEvw);

	auto& S_e_0 = ::std::get<3>(inter_nEvw);

	::std::vector<::std::vector<double>> temp_S_e_0(S_e_0[0].size(), ::std::vector<double>(S_e_0.size() * 6,0.0));

	for (int index_seId = 0; index_seId < S_e_0[0].size(); index_seId++)
	{
		for (int index_timestep = 0; index_timestep < S_e_0.size(); index_timestep++)
		{
			for (int index_n = 0; index_n < S_e_0[index_timestep][index_seId].size(); index_n++)
			{
				temp_S_e_0[index_seId][index_timestep * 6] += 1.0 / S_e_0[index_timestep][index_seId].size() * S_e_0[index_timestep][index_seId][index_n](0, 0);
				temp_S_e_0[index_seId][index_timestep * 6 + 1] += 1.0 / S_e_0[index_timestep][index_seId].size() * S_e_0[index_timestep][index_seId][index_n](1, 1);
				temp_S_e_0[index_seId][index_timestep * 6 + 2] += 1.0 / S_e_0[index_timestep][index_seId].size() * S_e_0[index_timestep][index_seId][index_n](2, 2);
				temp_S_e_0[index_seId][index_timestep * 6 + 3] += 1.0 / S_e_0[index_timestep][index_seId].size() * S_e_0[index_timestep][index_seId][index_n](0, 1);
				temp_S_e_0[index_seId][index_timestep * 6 + 4] += 1.0 / S_e_0[index_timestep][index_seId].size() * S_e_0[index_timestep][index_seId][index_n](1, 2);
				temp_S_e_0[index_seId][index_timestep * 6 + 5] += 1.0 / S_e_0[index_timestep][index_seId].size() * S_e_0[index_timestep][index_seId][index_n](2, 0);
			}
		}
	}
	write_vector2D_to_csv(temp_S_e_0, "./temp/debug/result/S_e_0.csv");

	::std::ofstream exter_nEvw_file;
	exter_nEvw_file.open("./temp/debug/result/T_exter_nEvw.csv", ::std::ios::ate | ::std::ios::out);
	for (int index_timestep = 0; index_timestep < exter_nEvw.size(); index_timestep++)
	{
		for (int index_vf = 0; index_vf < exter_nEvw[index_timestep].size(); index_vf++)
		{
			exter_nEvw_file << ::std::setprecision(12) << exter_nEvw[index_timestep][index_vf];
			if (index_vf != exter_nEvw[index_timestep].size() - 1)
			{
				exter_nEvw_file << ",";
			}
		}
		if (index_timestep != exter_nEvw.size() - 1)
		{
			exter_nEvw_file << "\n";
		}
	}
	exter_nEvw_file.close();

	task->configure.optim_function_T = [&fem, &outFile, &timeArray, &exter_nEvw, &inter_nEvw, &params](const ::std::vector<double>& ps) {
		double p_E = ps[0];
		double p_g = ps[1];
		double p_t = ps[2];

		double value = fun_for_optim_T(fem, p_g, p_t, p_E, outFile, timeArray, exter_nEvw, ::std::get<0>(inter_nEvw), ::std::get<2>(inter_nEvw), ::std::get<3>(inter_nEvw), params->virtualstrainArrayV, params->trueJArray, params->truedeformationGradientArray);

		return value;
		};

	const bool CONST_fun_compare_W = true;
	if constexpr (CONST_fun_compare_W == true)
	{
		try
		{
			task->pyfile_module.attr("fun_compare_W")(task->configure, exter_nEvw, ::std::get<0>(inter_nEvw), ::std::get<2>(inter_nEvw), ::std::get<3>(inter_nEvw), timeArray);
		}
		catch (const std::exception& e)
		{
			std::cout << e.what() << '\n';
		}
	}

	const bool CONST_optim_before = true;
	if constexpr (CONST_optim_before == true)
	{
		// py call before_optim
		try
		{
			task->pyfile_module.attr("OptimByPython")(params->LaplaceVFM_s, params->fs);
		}
		catch (const std::exception& e)
		{
			std::cout << e.what() << '\n';
		}

		try
		{
			::std::vector<double> t(timeArray);
			for (int i = 0; i < t.size(); ++i)
			{
				t[i] -= timeArray[0];
			}

			task->pyfile_module.attr("OptimByPython_invF")(t, ft);
		}
		catch (const std::exception& e)
		{
			std::cout << e.what() << '\n';
		}
	}

	// py call before_optim
	try
	{
		task->pyfile_module.attr("BeforeOptim")(task->configure).cast<VFMTask_configure>();
	}
	catch (const std::exception& e)
	{
		std::cout << e.what() << '\n';

		return true;
	}
#pragma region Opt10x10x10

	const bool CONST_Opt10x10x10 = false;
	if constexpr (CONST_Opt10x10x10 == true)
	{
		auto dop = [&](const ::std::vector<::std::vector<double>>& psss, ::std::function<double(const ::std::vector<double>&)> fun, bool ifoutlog = false, bool ifoutTacplot = false, ::std::vector<int> contour_num = { 9,9,1 }, ::std::string outTacplot_filepath = "./temp/debug/Tacplot.dat") {
			double min_value = std::numeric_limits<double>::infinity();
			::std::vector<double> min_ps;

			::std::multimap<double, ::std::vector<double>> map;

			// create a Tacplot .dat txt file
			::std::ofstream TacplotFile;
			if (ifoutTacplot)
			{
				TacplotFile.open(outTacplot_filepath);
			}

			// write Tacplot head
			if (ifoutTacplot)
			{
				TacplotFile << "TITLE = \"Optim\"\n";
				TacplotFile << "VARIABLES = \"E\", \"g\", \"t\", \"loss\"\n";
				TacplotFile << "Zone T=\"Optim\", I=" << contour_num[0] << ", J=" << contour_num[1] << ", K=" << contour_num[2] << ", F = POINT\n";
			}

			for (auto& pss : psss)
			{
				double value = fun(pss);

				// write E, g, t, loss to Tacplot
				if (ifoutTacplot)
				{
					TacplotFile << ::std::setprecision(12) << pss[0] << " " << pss[1] << " " << pss[2] << " " << value << "\n";
				}

				map.insert({ value,pss });

				if (value < min_value)
				{
					min_value = value;
					min_ps = pss;
				}

				if (ifoutlog)
				{
					// output the optimized parameters
					::std::stringstream ss;
					ss << "loss value: \n" << ::std::setprecision(12) << value << "\n";
					write_to_log_2(fem, ss.str(), outFile);
					for (auto& p : pss)
					{
						::std::stringstream ss;
						ss << ::std::setprecision(12) << p << ",";
						write_to_log_2(fem, ss.str(), outFile);
					}
					write_to_log_2(fem, "\n", outFile);
				}

			}

			// close Tacplot file
			if (ifoutTacplot)
			{
				TacplotFile.close();
			}

			if (ifoutlog)
			{
				// output the optimized parameters
				::std::stringstream ss;
				ss << "loss value: \n" << ::std::setprecision(12) << min_value << "\n";
				write_to_log_2(fem, ss.str(), outFile);
				for (auto& p : min_ps)
				{
					::std::stringstream ss;
					ss << ::std::setprecision(12) << p << ",";
					write_to_log_2(fem, ss.str(), outFile);
				}
				write_to_log_2(fem, "\n", outFile);
			}

			return ::std::make_tuple(min_value, min_ps, map);
			};

		auto xxv = function({ 0.3,8,1 });

		auto xxvv = function({ 0.1,2,1 });

		//auto grid = createGrid({ 0.1,1,1 }, { 0.9,9,9 }, { 0.1,1,0.1 });
		auto grid = createGrid({ 0.1,1,1 }, { 0.9,9,9 }, { 0.1,1,1 });

		auto [min_value, min_ps, dop_map] = dop(grid, function, true, true);
		double min_E = min_ps[0];
		double min_g = min_ps[1];
		double min_t = min_ps[2];

		// Convert dop_map to Eigen::MatrixXd [[dop_map.key, dop_map.value[0],dop_map.value[1],dop_map.value[2]],...]
		Eigen::MatrixXd dop_map_mat(dop_map.size(), 4);
		int index_dop_map = 0;
		for (auto& [value, ps] : dop_map)
		{
			dop_map_mat(index_dop_map, 0) = value;
			for (int j = 1; j < 4; ++j)
			{
				dop_map_mat(index_dop_map, j) = ps[j - 1];
			}
			index_dop_map++;
		}

		try
		{
			task->pyfile_module.attr("ShowValues")(dop_map_mat, "./temp/debug/result/dop_map.png");
		}
		catch (const std::exception& e)
		{
			std::cout << e.what() << '\n';

			return true;
		}


		// write to log 2, the first 20 in dop_map
		int i = 0;
		for (auto& [value, ps] : dop_map)
		{
			::std::stringstream ss;
			ss << "loss value: \n" << ::std::setprecision(12) << value << "\n";
			write_to_log_2(fem, ss.str(), outFile);
			for (auto& p : ps)
			{
				::std::stringstream ss;
				ss << ::std::setprecision(12) << p << ",";
				write_to_log_2(fem, ss.str(), outFile);
			}
			write_to_log_2(fem, "\n", outFile);

			i++;
			if (i > 20)
			{
				break;
			}
		}

		write_to_log_2(fem, "min_E: " + ::std::to_string(min_E) + " min_g: " + ::std::to_string(min_g) + " min_t: " + ::std::to_string(min_t) + "\n", outFile);

		//auto grid1 = createGrid({ ::std::max(min_E - 0.1,0.1),::std::max(min_g - 1,1.0),::std::max(min_t - 1,1.0) }, { ::std::min(min_E + 0.1,0.9),::std::min(min_g + 1,9.0),::std::min(min_t + 1,9.0) }, { 0.02,0.2,0.2 });
		auto grid1 = createGrid({ ::std::max(min_E - 0.1,0.1),::std::max(min_g - 1,1.0),1.0 }, { ::std::min(min_E + 0.1,0.9),::std::min(min_g + 1,9.0),1.0 }, { 0.02,0.2,0 });

		auto [min_value1, min_ps1, dop_map1] = dop(grid1, function, true);
		double min_E1 = min_ps1[0];
		double min_g1 = min_ps1[1];
		double min_t1 = min_ps1[2];

		// write to log 2, the first 20 in dop_map
		int i1 = 0;
		for (auto& [value, ps] : dop_map1)
		{
			::std::stringstream ss;
			ss << "loss value: \n" << ::std::setprecision(12) << value << "\n";
			write_to_log_2(fem, ss.str(), outFile);
			for (auto& p : ps)
			{
				::std::stringstream ss;
				ss << ::std::setprecision(12) << p << ",";
				write_to_log_2(fem, ss.str(), outFile);
			}
			write_to_log_2(fem, "\n", outFile);

			i1++;
			if (i1 > 20)
			{
				break;
			}
		}

		write_to_log_2(fem, "min_E1: " + ::std::to_string(min_E1) + " min_g1: " + ::std::to_string(min_g1) + " min_t1: " + ::std::to_string(min_t1) + "\n", outFile);

		//auto grid2 = createGrid({ ::std::max(min_E1 - 0.02,0.1),::std::max(min_g1 - 0.2,1.0),::std::max(min_t1 - 0.2,1.0) }, { ::std::min(min_E1 + 0.02,0.9),::std::min(min_g1 + 0.2,9.0),::std::min(min_t1 + 0.2,9.0) }, { 0.004,0.04,0.04 });
		auto grid2 = createGrid({ ::std::max(min_E1 - 0.02,0.1),::std::max(min_g1 - 0.2,1.0),1.0 }, { ::std::min(min_E1 + 0.02,0.9),::std::min(min_g1 + 0.2,9.0),1.0 }, { 0.004,0.04,0 });

		auto [min_value2, min_ps2, dop_map2] = dop(grid2, function, true);
		double min_E2 = min_ps2[0];
		double min_g2 = min_ps2[1];
		double min_t2 = min_ps2[2];

		// write to log 2, the first 20 in dop_map
		int i2 = 0;
		for (auto& [value, ps] : dop_map2)
		{
			::std::stringstream ss;
			ss << "loss value: \n" << ::std::setprecision(12) << value << "\n";
			write_to_log_2(fem, ss.str(), outFile);
			for (auto& p : ps)
			{
				::std::stringstream ss;
				ss << ::std::setprecision(12) << p << ",";
				write_to_log_2(fem, ss.str(), outFile);
			}
			write_to_log_2(fem, "\n", outFile);

			i2++;
			if (i2 > 20)
			{
				break;
			}
		}

		write_to_log_2(fem, "min_E2: " + ::std::to_string(min_E2) + " min_g2: " + ::std::to_string(min_g2) + " min_t2: " + ::std::to_string(min_t2) + "\n", outFile);

	}
#pragma endregion

#pragma region NLpot_0
	const bool CONST_NLpot_0 = true;
	if constexpr (CONST_NLpot_0 == true)
	{
		write_to_log_2(fem, "Start NLpot_0 :\n", outFile);

		nlopt::opt opt(nlopt::LN_COBYLA, 3);
		//nlopt::opt opt(nlopt::GN_DIRECT, 3);
		//nlopt::opt opt(nlopt::LN_NELDERMEAD, 3);
		//nlopt::opt opt(nlopt::GN_ISRES, 3);
		::std::vector<double> lb(3);
		lb[0] = 0.001; lb[1] = 0.01; lb[2] = 0.01;
		opt.set_lower_bounds(lb);
		::std::vector<double> ub(3);
		ub[0] = 0.9; ub[1] = 9.0; ub[2] = 9.0;
		opt.set_upper_bounds(ub);

		opt.set_min_objective(fun_nlpot, &function);
		//opt.set_min_objective(fun_nlpot00005, &fun);

		opt.set_ftol_rel(1e-6);

		opt.set_xtol_rel(0.00005);

		// Set the initial point
		::std::vector<double> x(3);
		//x[0] = 0.48862; x[1] = 4.44033; x[2] = 1.00015;
		x[0] = 0.1; x[1] = 1; x[2] = 1;
		//x[0] = 0.265891; x[1] = 8.9923; x[2] = 1.00006;
		//x[0] = 0.5; x[1] = 5; x[2] = 2;

		// Perform the optimization
		double minf;
		try {
			nlopt::result result = opt.optimize(x, minf);

			::std::stringstream ss;
			ss << "Optimization succeeded, found minimum at f("
				<< x[0] << "," << x[1] << "," << x[2] << ") = "
				<< minf << std::endl;
			write_to_log_2(fem, ss.str(), outFile);
		}
		catch (std::exception& e) {
			::std::stringstream ss;
			ss << "Optimization failed: " << e.what() << std::endl;
			write_to_log_2(fem, ss.str(), outFile);
		}
	}

#pragma endregion

#pragma region NLpot_more
	const bool CONST_NLpot_more = false;
	if constexpr (CONST_NLpot_more == true)
	{
		write_to_log_2(fem, "Start NLpot_0 :\n", outFile);

		auto ms = { nlopt::GN_DIRECT,nlopt::GN_DIRECT_L,nlopt::GN_CRS2_LM,nlopt::GN_ISRES };
		::std::vector<::std::vector<double>> min_f_xs;
		for (auto m : ms)
		{
			//nlopt::opt opt(nlopt::LN_COBYLA, 3);
			nlopt::opt opt(m, 3);
			//nlopt::opt opt(nlopt::LN_NELDERMEAD, 3);
			//nlopt::opt opt(nlopt::GN_ISRES, 3);
			::std::vector<double> lb(3);
			lb[0] = 0.1; lb[1] = 1.0; lb[2] = 0.1;
			opt.set_lower_bounds(lb);
			::std::vector<double> ub(3);
			ub[0] = 0.9; ub[1] = 9.0; ub[2] = 3.0;
			opt.set_upper_bounds(ub);

			opt.set_min_objective(fun_nlpot, &function);
			//opt.set_min_objective(fun_nlpot00005, &fun);

			opt.set_ftol_rel(1e-6);

			opt.set_xtol_rel(0.00005);

			// Set the initial point
			::std::vector<double> x(3);
			//x[0] = 0.48862; x[1] = 4.44033; x[2] = 1.00015;
			x[0] = 0.3; x[1] = 8; x[2] = 1;
			//x[0] = 0.265891; x[1] = 8.9923; x[2] = 1.00006;
			//x[0] = 0.5; x[1] = 5; x[2] = 2;

			// Perform the optimization
			double minf;
			try {
				nlopt::result result = opt.optimize(x, minf);

				::std::stringstream ss;
				ss << ::std::setprecision(12) << "Optimization succeeded, found minimum at f("
					<< x[0] << "," << x[1] << "," << x[2] << ") = "
					<< minf << std::endl;
				write_to_log_2(fem, ss.str(), outFile);

				auto f_xs = ::std::vector<double>(x);
				f_xs.insert(f_xs.begin(), minf);
				min_f_xs.push_back(f_xs);
			}
			catch (std::exception& e) {
				::std::stringstream ss;
				ss << "Optimization failed: " << e.what() << std::endl;
				write_to_log_2(fem, ss.str(), outFile);
			}


		}

		int index_m = 0;
		for (auto f_xs : min_f_xs)
		{
			::std::stringstream ss;
			ss << ::std::setprecision(12) << index_m << " _ Optimization succeeded, found minimum at f("
				<< f_xs[1] << "," << f_xs[2] << "," << f_xs[3] << ") = "
				<< f_xs[0] << std::endl;
			write_to_log_2(fem, ss.str(), outFile);
			index_m++;
		}
	}

#pragma endregion



#pragma endregion

	outFile.close();

	return true;
}


void L_externalVirtualWork_LaplcaeTransform(std::vector<std::vector<double>>& externalVirtualWork, std::vector<std::vector<Eigen::dcomplex>>& externalVirtualWork_laplace, std::vector<double>& timeArray, VFMTask* task)
{
	int vf_u_size = externalVirtualWork[0].size();

	externalVirtualWork_laplace = ::std::vector<::std::vector<::std::complex<double>>>(vf_u_size);

	for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
	{
		::std::vector<double> externalVirtualWork_VF;
		//从externalVirtualWork中切片[:][index_vf]作为externalVirtualWork_VF
		for (int index_timestep = 0; index_timestep < externalVirtualWork.size(); index_timestep++)
		{
			externalVirtualWork_VF.push_back(externalVirtualWork[index_timestep][index_vf]);
		}
		externalVirtualWork_laplace[index_vf] = laplace_transform_periodic(externalVirtualWork_VF, timeArray, task->configure.LaplaceVFM_s);
	}
}

void L_StressPK2_withJc_LaplaceTransform(std::vector<double>& timeArray, VFMTask* task, FEModel* fem, std::vector<int>& solution_elementsDomainID, std::vector<std::vector<std::vector<mat3d>>>& truedeformationGradientArray, std::vector<std::vector<std::vector<double>>>& trueJArray, std::vector<std::vector<std::vector<std::vector<Eigen::dcomplex>>>>& Sepk2_dotdot_vStrain, std::vector<std::vector<std::vector<std::vector<mat3ds>>>>& virtualstrainArrayV)
{
	FEMesh& mesh=fem->GetMesh();
	auto vf_u_size = virtualstrainArrayV[0].size();
	// Stress Elastic PK2 with Jc in gauss points
	// [index_timestep, index_ElementID, index_gauss]
	::std::vector<::std::vector<::std::vector<mat3ds>>> StressPK2_withJc(timeArray.size());
	for (int index_timestep = 0; index_timestep < timeArray.size(); index_timestep++)
	{
		StressPK2_withJc[index_timestep] = ::std::vector<::std::vector<mat3ds>>(task->solution_elementsID.size());
		for (int index_seId = 0; index_seId < task->solution_elementsID.size(); index_seId++)
		{
			int j = task->solution_elementsID[index_seId];

			FEElement& element = *(mesh.Element(j));
			int nodes_number = element.Nodes();

			// find material by element
			FEMaterial* pmat = fem->GetMaterial(element.GetMatID());

			// convert element to solid element
			FESolidElement& solidElement = static_cast<FESolidElement&>(element);

			// get domain
			FESolidDomain& domain = static_cast<FESolidDomain&>(mesh.Domain(solution_elementsDomainID[index_seId]));

			StressPK2_withJc[index_timestep][index_seId] = ::std::vector<mat3ds>(solidElement.GaussPoints());
			for (int n = 0; n < solidElement.GaussPoints(); n++)
			{
				// get element's stress
				FEMaterialPoint& mp = *solidElement.GetMaterialPoint(n);

				FEElasticMaterialPoint& pt_e = *mp.ExtractData<FEElasticMaterialPoint>();

				FEViscoElasticMaterialPoint* p_pt = mp.ExtractData<FEViscoElasticMaterialPoint>();

				FEViscoElasticMaterialPoint& pt = *p_pt;

				pt_e.m_F = truedeformationGradientArray[index_timestep][index_seId][n];
				pt_e.m_J = truedeformationGradientArray[index_timestep][index_seId][n].det();// trueJArray[index_timestep][index_seId][n];
				double Jc = trueJArray[index_timestep][index_seId][n];

				fem->GetTime().timeIncrement = index_timestep == 0 ? (timeArray[index_timestep + 1] - timeArray[index_timestep]) : (timeArray[index_timestep] - timeArray[index_timestep - 1]);
				mat3ds pts;
				if (p_pt != nullptr)
				{
					FEViscoElasticMaterial* material = static_cast<FEViscoElasticMaterial*>(pmat);
					FECoreBase* elasticb = material->FindProperty("elastic")->get(0);

					// FENeoHookean
					FENeoHookean& elasticm = static_cast<FENeoHookean&>(*elasticb);
					FEParamDouble* material_E = elasticm.GetParameter("E")->pvalue<FEParamDouble>();
					double linerElactical_e = material_E->constValue();

					pts = dynamic_cast<FESolidMaterial*>(pmat)->Stress(mp);
					pts = pt_e.pull_back(pts);

					pts /= linerElactical_e;
				}
				else
				{
					pts = dynamic_cast<FESolidMaterial*>(pmat)->Stress(mp);
					pts = pt_e.pull_back(pts);
				}

				StressPK2_withJc[index_timestep][index_seId][n] = pts * Jc;
			}
		}
	}

	::std::vector<::std::vector<::std::vector<Eigen::Matrix<std::complex<double>, 3, 3>>>> StressPK2_withJc_unitviscoE_laplace(task->solution_elementsID.size());
	for (int index_seId = 0; index_seId < task->solution_elementsID.size(); index_seId++)
	{
		StressPK2_withJc_unitviscoE_laplace[index_seId] = ::std::vector<::std::vector<Eigen::Matrix<std::complex<double>, 3, 3>>>(StressPK2_withJc[0][index_seId].size(), ::std::vector<Eigen::Matrix<std::complex<double>, 3, 3>>(timeArray.size()));
		for (int n = 0; n < StressPK2_withJc[0][index_seId].size(); n++)
		{
			::std::vector<mat3ds> stresspk2_withJc_timeorder(timeArray.size());
			for (int index_timestep = 0; index_timestep < timeArray.size(); index_timestep++)
			{
				stresspk2_withJc_timeorder[index_timestep] = StressPK2_withJc[index_timestep][index_seId][n];
			}
			StressPK2_withJc_unitviscoE_laplace[index_seId][n] = laplace_transform_periodic(stresspk2_withJc_timeorder, timeArray, task->configure.LaplaceVFM_s);
		}
	}

	Sepk2_dotdot_vStrain.resize(vf_u_size);
	for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
	{
		Sepk2_dotdot_vStrain[index_vf] = ::std::vector<::std::vector<::std::vector<::std::complex<double>>>>(task->solution_elementsID.size());
		for (int index_seId = 0; index_seId < task->solution_elementsID.size(); index_seId++)
		{
			int j = task->solution_elementsID[index_seId];

			FEElement& element = *(mesh.Element(j));
			int nodes_number = element.Nodes();

			// find material by element
			FEMaterial* pmat = fem->GetMaterial(element.GetMatID());

			// convert element to solid element
			FESolidElement& solidElement = static_cast<FESolidElement&>(element);

			Sepk2_dotdot_vStrain[index_vf][index_seId] = ::std::vector<::std::vector<::std::complex<double>>>(solidElement.GaussPoints());
			for (int n = 0; n < solidElement.GaussPoints(); n++)
			{
				FEMaterialPoint& mp = *solidElement.GetMaterialPoint(n);

				FEViscoElasticMaterialPoint* p_pt = mp.ExtractData<FEViscoElasticMaterialPoint>();

				auto& vStrain = virtualstrainArrayV[0][index_vf][index_seId][n]; // laplace set vf=vf(0);
				auto vStrain_d = convert_mat3ds_EigenMatrix(vStrain);

				Sepk2_dotdot_vStrain[index_vf][index_seId][n] = ::std::vector<::std::complex<double>>(task->configure.LaplaceVFM_s.size());
				for (int index_s = 0; index_s < task->configure.LaplaceVFM_s.size(); index_s++)
				{
					Sepk2_dotdot_vStrain[index_vf][index_seId][n][index_s] = (StressPK2_withJc_unitviscoE_laplace[index_seId][n][index_s].array() * vStrain_d.array()).sum();
				}
			}
		}
	}
}


::std::vector< ::std::vector<::std::complex<double>>> L_fs(std::vector<std::vector<double>>& externalVirtualWork, std::vector<double>& timeArray, VFMTask* task, FEModel* fem, std::vector<int>& solution_elementsDomainID, std::vector<std::vector<std::vector<mat3d>>>& truedeformationGradientArray, std::vector<std::vector<std::vector<double>>>& trueJArray, std::vector<std::vector<std::vector<std::vector<mat3ds>>>>& virtualstrainArrayV)
{
	FEMesh& mesh = fem->GetMesh();
	auto vf_u_size = virtualstrainArrayV[0].size();

	std::vector<std::vector<Eigen::dcomplex>> externalVirtualWork_laplace_inT = ::std::vector<::std::vector<::std::complex<double>>>(vf_u_size);

	for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
	{
		::std::vector<double> externalVirtualWork_VF;
		//从externalVirtualWork中切片[:][index_vf]作为externalVirtualWork_VF
		for (int index_timestep = 0; index_timestep < externalVirtualWork.size(); index_timestep++)
		{
			externalVirtualWork_VF.push_back(externalVirtualWork[index_timestep][index_vf]);
		}
		externalVirtualWork_laplace_inT[index_vf] = laplace_transform_inT(externalVirtualWork_VF, timeArray, task->configure.LaplaceVFM_s);
	}

	::std::vector<::std::vector<::std::vector<mat3ds>>> StressPK2_withJc(timeArray.size());
	for (int index_timestep = 0; index_timestep < timeArray.size(); index_timestep++)
	{
		StressPK2_withJc[index_timestep] = ::std::vector<::std::vector<mat3ds>>(task->solution_elementsID.size());
		for (int index_seId = 0; index_seId < task->solution_elementsID.size(); index_seId++)
		{
			int j = task->solution_elementsID[index_seId];

			FEElement& element = *(mesh.Element(j));
			int nodes_number = element.Nodes();

			// find material by element
			FEMaterial* pmat = fem->GetMaterial(element.GetMatID());

			// convert element to solid element
			FESolidElement& solidElement = static_cast<FESolidElement&>(element);

			// get domain
			FESolidDomain& domain = static_cast<FESolidDomain&>(mesh.Domain(solution_elementsDomainID[index_seId]));

			StressPK2_withJc[index_timestep][index_seId] = ::std::vector<mat3ds>(solidElement.GaussPoints());
			for (int n = 0; n < solidElement.GaussPoints(); n++)
			{
				// get element's stress
				FEMaterialPoint& mp = *solidElement.GetMaterialPoint(n);

				FEElasticMaterialPoint& pt_e = *mp.ExtractData<FEElasticMaterialPoint>();

				FEViscoElasticMaterialPoint* p_pt = mp.ExtractData<FEViscoElasticMaterialPoint>();

				FEViscoElasticMaterialPoint& pt = *p_pt;

				pt_e.m_F = truedeformationGradientArray[index_timestep][index_seId][n];
				pt_e.m_J = truedeformationGradientArray[index_timestep][index_seId][n].det();// trueJArray[index_timestep][index_seId][n];
				double Jc = trueJArray[index_timestep][index_seId][n];

				fem->GetTime().timeIncrement = index_timestep == 0 ? (timeArray[index_timestep + 1] - timeArray[index_timestep]) : (timeArray[index_timestep] - timeArray[index_timestep - 1]);
				mat3ds pts;
				if (p_pt != nullptr)
				{
					FEViscoElasticMaterial* material = static_cast<FEViscoElasticMaterial*>(pmat);
					FECoreBase* elasticb = material->FindProperty("elastic")->get(0);

					// FENeoHookean
					FENeoHookean& elasticm = static_cast<FENeoHookean&>(*elasticb);
					FEParamDouble* material_E = elasticm.GetParameter("E")->pvalue<FEParamDouble>();
					double linerElactical_e = material_E->constValue();

					auto pts_neohookean = elasticm.Stress(mp);
					auto pts_visco = material->Stress(mp);

					//pts = elasticm.Stress(mp);
					//FESolidMaterial* pmat_so = dynamic_cast<FESolidMaterial*>(pmat);
					
					//pts = pmat_so->Stress(mp);

					pts = pts_neohookean;
					pts = pt_e.pull_back(pts);

					pts /= linerElactical_e;
				}
				else
				{
					pts = dynamic_cast<FESolidMaterial*>(pmat)->Stress(mp);
					pts = pt_e.pull_back(pts);
				}

				StressPK2_withJc[index_timestep][index_seId][n] = pts * Jc;
			}
		}
	}

	::std::vector<::std::vector<::std::vector<Eigen::Matrix<std::complex<double>, 3, 3>>>> StressPK2_withJc_unitviscoE_laplace_inT(task->solution_elementsID.size());
	for (int index_seId = 0; index_seId < task->solution_elementsID.size(); index_seId++)
	{
		StressPK2_withJc_unitviscoE_laplace_inT[index_seId] = ::std::vector<::std::vector<Eigen::Matrix<std::complex<double>, 3, 3>>>(StressPK2_withJc[0][index_seId].size(), ::std::vector<Eigen::Matrix<std::complex<double>, 3, 3>>(timeArray.size()));
		for (int n = 0; n < StressPK2_withJc[0][index_seId].size(); n++)
		{
			::std::vector<mat3ds> stresspk2_withJc_timeorder(timeArray.size());
			for (int index_timestep = 0; index_timestep < timeArray.size(); index_timestep++)
			{
				stresspk2_withJc_timeorder[index_timestep] = StressPK2_withJc[index_timestep][index_seId][n];
			}
			StressPK2_withJc_unitviscoE_laplace_inT[index_seId][n] = laplace_transform_inT(stresspk2_withJc_timeorder, timeArray, task->configure.LaplaceVFM_s);
		}
	}

	std::vector<std::vector<std::vector<std::vector<Eigen::dcomplex>>>> Sepk2_dotdot_vStrain(vf_u_size);

	for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
	{
		Sepk2_dotdot_vStrain[index_vf] = ::std::vector<::std::vector<::std::vector<::std::complex<double>>>>(task->solution_elementsID.size());
		for (int index_seId = 0; index_seId < task->solution_elementsID.size(); index_seId++)
		{
			int j = task->solution_elementsID[index_seId];

			FEElement& element = *(mesh.Element(j));
			int nodes_number = element.Nodes();

			// find material by element
			FEMaterial* pmat = fem->GetMaterial(element.GetMatID());

			// convert element to solid element
			FESolidElement& solidElement = static_cast<FESolidElement&>(element);

			Sepk2_dotdot_vStrain[index_vf][index_seId] = ::std::vector<::std::vector<::std::complex<double>>>(solidElement.GaussPoints());
			for (int n = 0; n < solidElement.GaussPoints(); n++)
			{
				FEMaterialPoint& mp = *solidElement.GetMaterialPoint(n);

				FEViscoElasticMaterialPoint* p_pt = mp.ExtractData<FEViscoElasticMaterialPoint>();

				auto& vStrain = virtualstrainArrayV[0][index_vf][index_seId][n]; // laplace set vf=vf(0);
				auto vStrain_d = convert_mat3ds_EigenMatrix(vStrain);

				Sepk2_dotdot_vStrain[index_vf][index_seId][n] = ::std::vector<::std::complex<double>>(task->configure.LaplaceVFM_s.size());
				for (int index_s = 0; index_s < task->configure.LaplaceVFM_s.size(); index_s++)
				{
					Sepk2_dotdot_vStrain[index_vf][index_seId][n][index_s] = (StressPK2_withJc_unitviscoE_laplace_inT[index_seId][n][index_s].array() * vStrain_d.array()).sum();
				}
			}
		}
	}


	::std::vector< ::std::vector<::std::complex<double>>> fs(vf_u_size, ::std::vector<::std::complex<double>>(task->configure.LaplaceVFM_s.size(), ::std::complex<double>(0, 0)));
	for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
	{
		::std::vector<::std::vector<::std::vector<Eigen::dcomplex>>> stresspk2_withJc_VF_alle_laplace_visco(task->solution_elementsID.size());
		::std::vector<::std::vector<::std::vector<Eigen::dcomplex>>> stresspk2_withJc_VF_alle_laplace_elastic(task->solution_elementsID.size());
		for (int index_seId = 0; index_seId < task->solution_elementsID.size(); index_seId++)
		{
			int j = task->solution_elementsID[index_seId];

			FEElement& element = *(mesh.Element(j));
			int nodes_number = element.Nodes();

			// find material by element
			FEMaterial* pmat = fem->GetMaterial(element.GetMatID());


			// convert element to solid element
			FESolidElement& solidElement = static_cast<FESolidElement&>(element);

			stresspk2_withJc_VF_alle_laplace_visco[index_seId].resize(solidElement.GaussPoints(), ::std::vector<Eigen::dcomplex>(task->configure.LaplaceVFM_s.size(), Eigen::dcomplex(0, 0)));
			stresspk2_withJc_VF_alle_laplace_elastic[index_seId].resize(solidElement.GaussPoints(), ::std::vector<Eigen::dcomplex>(task->configure.LaplaceVFM_s.size(), Eigen::dcomplex(0, 0)));
			for (int index_gauss = 0; index_gauss < solidElement.GaussPoints(); index_gauss++)
			{
				// get element's stress
				FEMaterialPoint& mp = *solidElement.GetMaterialPoint(index_gauss);

				FEElasticMaterialPoint& pt_e = *mp.ExtractData<FEElasticMaterialPoint>();

				FEViscoElasticMaterialPoint* p_pt = mp.ExtractData<FEViscoElasticMaterialPoint>();

				for (int index_s = 0; index_s < task->configure.LaplaceVFM_s.size(); index_s++)
				{
					auto SdS = Sepk2_dotdot_vStrain[index_vf][index_seId][index_gauss][index_s];
					if (p_pt != nullptr)
					{
						stresspk2_withJc_VF_alle_laplace_visco[index_seId][index_gauss][index_s] = SdS;
					}
					else
					{
						stresspk2_withJc_VF_alle_laplace_elastic[index_seId][index_gauss][index_s] = SdS;
					}

				}
			}
		}

		::std::vector<::std::complex<double>> internalVirtualWork_laplace_visco(task->configure.LaplaceVFM_s.size(), ::std::complex<double>(0, 0));
		::std::vector<::std::complex<double>> internalVirtualWork_laplace_elastic(task->configure.LaplaceVFM_s.size(), ::std::complex<double>(0, 0));


		for (int index_s = 0; index_s < task->configure.LaplaceVFM_s.size(); index_s++)
		{
			for (int index_seId = 0; index_seId < task->solution_elementsID.size(); index_seId++)
			{
				for (int index_gauss = 0; index_gauss < stresspk2_withJc_VF_alle_laplace_visco[index_seId].size(); index_gauss++)
				{
					if (stresspk2_withJc_VF_alle_laplace_visco[index_seId][index_gauss].size() > 0)
					{
						internalVirtualWork_laplace_visco[index_s] += stresspk2_withJc_VF_alle_laplace_visco[index_seId][index_gauss][index_s];
					}

					if (stresspk2_withJc_VF_alle_laplace_elastic[index_seId][index_gauss].size() > 0)
					{
						internalVirtualWork_laplace_elastic[index_s] += stresspk2_withJc_VF_alle_laplace_elastic[index_seId][index_gauss][index_s];
					}

				}
			}
			fs[index_vf][index_s] = (externalVirtualWork_laplace_inT[index_vf][index_s] - internalVirtualWork_laplace_elastic[index_s]) / (internalVirtualWork_laplace_visco[index_s]);
		}
	}

	return fs;
}

::std::vector<::std::vector<double>> inv_ft(::std::vector<::std::vector<::std::complex<double>>> fs, ::std::vector<double> timeArray)
{
	int vf_u_size = fs.size();

	::std::vector<::std::vector<double>> ft(vf_u_size, ::std::vector<double>(timeArray.size(), 0));
	double loss_sum = 0;

	int N = 32;
	double shift = CONST_SHIFT;

	for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
	{
		// convert params.fs[index_vf] to Fs[index_time][index_s] each time has N s
		::std::vector<::std::vector<::std::complex<double>>> Fs(timeArray.size(), ::std::vector<::std::complex<double>>(N, ::std::complex<double>(0, 0)));
		for (int index_time = 0; index_time < timeArray.size(); index_time++)
		{
			for (int index_s = 0; index_s < N; index_s++)
			{
				Fs[index_time][index_s] = fs[index_vf][N * index_time + index_s];
			}
		}

		auto invF = talbotInverseLaplaceTransform(Fs, timeArray, N, shift);
		ft[index_vf] = invF;
	}

	return ft;
}


::std::tuple<::std::vector<::std::vector<double>>, ::std::vector<::std::vector<double>>, ::std::vector<int>, ::std::vector<::std::vector<::std::vector<mat3ds>>>, ::std::vector<::std::vector<::std::vector<::std::vector<double>>>>> cal_internal_normal_linerE_vw(std::vector<double>& timeArray, VFMTask* task, FEModel* fem, std::vector<int>& solution_elementsDomainID, std::vector<std::vector<std::vector<mat3d>>>& truedeformationGradientArray, std::vector<std::vector<std::vector<double>>>& trueJArray, std::vector<std::vector<std::vector<std::vector<mat3ds>>>>& virtualstrainArrayV)
{
	FEMesh& mesh = fem->GetMesh();
	auto vf_u_size = virtualstrainArrayV[0].size();

	::std::vector<::std::vector<double>> internal_normal_visco(timeArray.size(), ::std::vector<double>(vf_u_size, 0));
	::std::vector<::std::vector<double>> internal_elastic(timeArray.size(), ::std::vector<double>(vf_u_size, 0));

	::std::vector<::std::vector<::std::vector<::std::vector<double>>>> internal_normal_visco_dse0_strain_Jc(timeArray.size(), ::std::vector<::std::vector<::std::vector<double>>>(vf_u_size, ::std::vector<::std::vector<double>>(task->solution_elementsID.size())));

	::std::vector<::std::vector<::std::vector<mat3ds>>> S_e_0(timeArray.size(), ::std::vector<::std::vector<mat3ds>>(task->solution_elementsID.size()));
	::std::vector<int> visco_mask(task->solution_elementsID.size(), 0); // 0: elastic, 1: visco

	::std::vector<::std::vector<::std::vector<mat3ds>>> StressPK2(timeArray.size());
	for (int index_timestep = 0; index_timestep < timeArray.size(); index_timestep++)
	{
		StressPK2[index_timestep] = ::std::vector<::std::vector<mat3ds>>(task->solution_elementsID.size());
		for (int index_seId = 0; index_seId < task->solution_elementsID.size(); index_seId++)
		{
			int j = task->solution_elementsID[index_seId];

			FEElement& element = *(mesh.Element(j));
			int nodes_number = element.Nodes();

			// find material by element
			FEMaterial* pmat = fem->GetMaterial(element.GetMatID());

			// convert element to solid element
			FESolidElement& solidElement = static_cast<FESolidElement&>(element);

			// get domain
			FESolidDomain& domain = static_cast<FESolidDomain&>(mesh.Domain(solution_elementsDomainID[index_seId]));

			StressPK2[index_timestep][index_seId] = ::std::vector<mat3ds>(solidElement.GaussPoints());
			S_e_0[index_timestep][index_seId] = ::std::vector<mat3ds>(solidElement.GaussPoints());
			
			for (int n = 0; n < solidElement.GaussPoints(); n++)
			{
				// get element's stress
				FEMaterialPoint& mp = *solidElement.GetMaterialPoint(n);

				FEElasticMaterialPoint& pt_e = *mp.ExtractData<FEElasticMaterialPoint>();

				FEViscoElasticMaterialPoint* p_pt = mp.ExtractData<FEViscoElasticMaterialPoint>();

				FEViscoElasticMaterialPoint& pt = *p_pt;

				pt_e.m_F = truedeformationGradientArray[index_timestep][index_seId][n];
				pt_e.m_J = truedeformationGradientArray[index_timestep][index_seId][n].det();// trueJArray[index_timestep][index_seId][n];
				double Jc = trueJArray[index_timestep][index_seId][n];

				fem->GetTime().timeIncrement = index_timestep == 0 ? (timeArray[index_timestep + 1] - timeArray[index_timestep]) : (timeArray[index_timestep] - timeArray[index_timestep - 1]);
				mat3ds pts;
				if (p_pt != nullptr)
				{
					FEViscoElasticMaterial* material = static_cast<FEViscoElasticMaterial*>(pmat);
					FECoreBase* elasticb = material->FindProperty("elastic")->get(0);

					// FENeoHookean
					FENeoHookean& elasticm = static_cast<FENeoHookean&>(*elasticb);
					FEParamDouble* material_E = elasticm.GetParameter("E")->pvalue<FEParamDouble>();
					double linerElactical_e = material_E->constValue();

					auto pts_neohookean = elasticm.Stress(mp);
					auto pts_visco = material->Stress(mp);

					//pts = elasticm.Stress(mp);
					//FESolidMaterial* pmat_so = dynamic_cast<FESolidMaterial*>(pmat);

					//pts = pmat_so->Stress(mp);

					// cal lagrange strain
					// 
					// strainCompute infinitesimal strain
					mat3ds infstrain = pt_e.m_F.sym();
					// strainCompute green-lagrange strain
					mat3ds FEstrain = pt_e.Strain(); // virtual strain error m_r no real position
					mat3d I = mat3d(1, 0, 0, 0, 1, 0, 0, 0, 1);
					mat3d lstrain = (pt_e.m_F.transpose() * pt_e.m_F - I) * 0.5;
					// convert lstrain to mat3ds
					mat3ds strain = mat3ds(lstrain[0][0], lstrain[1][1], lstrain[2][2], lstrain[0][1], lstrain[1][2], lstrain[0][2]);

					mat3ds& s = FEstrain;

					pts = pts_neohookean;
					pts = pt_e.pull_back(pts);

					pts /= linerElactical_e;

					S_e_0[index_timestep][index_seId][n] = pts;
					if (index_timestep == 0)
					{
						visco_mask[index_seId]=1;
					}

					for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
					{
						internal_normal_visco[index_timestep][index_vf] += pts.dotdot(virtualstrainArrayV[index_timestep][index_vf][index_seId][n])*Jc;
					}

					if (index_seId == 1918)
					{
						auto oejdofeospkdsokf = 1;
					}

				}
				else
				{
					mat3ds pts_c = static_cast<FESolidMaterial*>(pmat)->Stress(mp);
					pts = pt_e.pull_back(pts_c);

					//mat3ds pts_iso_e_c = static_cast<FEIsotropicElastic*>(pmat)->Stress(mp);
					//mat3ds pts_iso_e = pt_e.pull_back(pts_iso_e_c);

					//mat3ds pts_iso_e = PK2stress(*dynamic_cast<FEIsotropicElastic*>(pmat), mp);

					//mat3ds pts_error = pts - pts_iso_e;
					//pts = pts_iso_e;

					S_e_0[index_timestep][index_seId][n] = pts;

					for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
					{
						internal_elastic[index_timestep][index_vf] += pts.dotdot(virtualstrainArrayV[index_timestep][index_vf][index_seId][n])*Jc;
					}
				}

				StressPK2[index_timestep][index_seId][n] = pts;
			}

			if (index_seId == 1918)
			{
				auto oejdofeospkdsokf = 1;
			}
		}
	}

	::std::vector<::std::vector<::std::vector<mat3ds>>> dStressPK2_dt(StressPK2);
	for (int index_timestep = 0; index_timestep < timeArray.size(); index_timestep++)
	{
		for (int index_seId = 0; index_seId < task->solution_elementsID.size(); index_seId++)
		{

			for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
			{
				internal_normal_visco_dse0_strain_Jc[index_timestep][index_vf][index_seId] = ::std::vector<double>(StressPK2[index_timestep][index_seId].size());
			}

			for (int n = 0; n < StressPK2[index_timestep][index_seId].size(); n++)
			{
				if (index_timestep == 0)
				{
					dStressPK2_dt[index_timestep][index_seId][n] = (StressPK2[index_timestep][index_seId][n] - StressPK2[StressPK2.size() - 1][index_seId][n]) / (timeArray[index_timestep + 1] - timeArray[index_timestep]);
				}
				else
				{
					dStressPK2_dt[index_timestep][index_seId][n] = (StressPK2[index_timestep][index_seId][n] - StressPK2[index_timestep - 1][index_seId][n]) / (timeArray[index_timestep] - timeArray[index_timestep - 1]);
				}

				for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
				{
					internal_normal_visco_dse0_strain_Jc[index_timestep][index_vf][index_seId][n] = dStressPK2_dt[index_timestep][index_seId][n].dotdot(virtualstrainArrayV[index_timestep][index_vf][index_seId][n]) * trueJArray[index_timestep][index_seId][n];
				}
			}
		}
	}

	return ::std::make_tuple(internal_normal_visco, internal_elastic, visco_mask, S_e_0, internal_normal_visco_dse0_strain_Jc);
}

// 输入：时间序列 ts（严格递增或部分有零尾巴），心率 bpm
// 输出：一对索引 (start_index, end_index)，表示倒数第二个完整周期的区间 [start, end)
std::pair<int, int> find_last_two_cycles(const std::span<double>& ts, double bpm)
{
	if (ts.empty()) {
		throw std::runtime_error("timestep is empty");
	}

	// 1) 找到最后一个有效点（过滤掉末尾的 0 或非递增）
	int last_valid = static_cast<int>(ts.size()) - 1;
	while (last_valid > 0) {
		bool invalid = (ts[last_valid] == 0.0) ||
			!(ts[last_valid] > ts[last_valid - 1]) ||
			std::isnan(ts[last_valid]);
		if (!invalid) break;
		--last_valid;
	}
	if (last_valid <= 0) {
		throw std::runtime_error("no valid timesteps found");
	}

	double period = 60.0 / bpm;

	// 辅助函数：给定一个 end 索引，往前找一个周期的起点
	auto find_cycle_start = [&](int end_index) {
		double target = ts[end_index] - period;
		// 在 [0, end_index] 上找第一个 >= target 的点
		auto first = ts.begin();
		auto it = std::lower_bound(first, first + (end_index + 1), target);
		int idx = static_cast<int>(it - first);
		if (idx > end_index) idx = end_index;
		return idx;
		};

	// 2) 找最后一个周期的起点
	int last_start_index = find_cycle_start(last_valid);
	int last_end_exclusive = last_valid + 1; // 半开区间，包含 last_valid

	// 3) 找倒数第二个周期的起点
	int second_last_end = last_start_index - 1;
	int second_last_start_index = (second_last_end >= 0)
		? find_cycle_start(second_last_end)
		: 0;
	int second_last_end_exclusive = last_start_index;

	// 4) 返回区间 [start, end)
	int start_index = second_last_start_index;
	int end_index = second_last_end_exclusive;

	// 保护：如果没能找到两个完整周期，就退回最后一个周期
	if (start_index >= end_index) {
		start_index = last_start_index;
		end_index = last_end_exclusive;
	}

	return { start_index, end_index };
}

static bool write_tecplot_nodal_position_and_displacement(
	const std::string& filepath,
	const std::span<double>& time,
	const std::vector<std::vector<double>>& initialCoordinate,
	const span3d<double>& timeDisplacement,
	const std::vector<int>& selected_node_list,
    FEMesh& mesh,
	const std::vector<int>& selected_element_list,
	int start_index,
	int end_index_inclusive)
{
	std::ofstream out(filepath, std::ios::out | std::ios::trunc);
	if (!out.is_open()) return false;

    const int nNodes = static_cast<int>(initialCoordinate.size());

	out << "TITLE = \"NodalPositionAndDisplacement\"\n";
	out << "VARIABLES = \"X\",\"Y\",\"Z\",\"U\",\"V\",\"W\",\"Time\",\"NodeId\"\n";

   out << std::setprecision(16);

	// determine which nodes to output (global node indices)
	std::vector<int> output_node_ids;
	output_node_ids.reserve(!selected_node_list.empty() ? selected_node_list.size() : static_cast<size_t>(nNodes));
	if (!selected_node_list.empty())
	{
		output_node_ids.assign(selected_node_list.begin(), selected_node_list.end());
	}
	else
	{
		for (int nid = 0; nid < nNodes; ++nid) output_node_ids.push_back(nid);
	}

	const int output_nodes = static_cast<int>(output_node_ids.size());

	// build global->local node index map (local index is 0-based)
	std::vector<int> global_to_local(nNodes, -1);
	for (int i = 0; i < output_nodes; ++i)
	{
		const int gid = output_node_ids[i];
		if (gid >= 0 && gid < nNodes) global_to_local[gid] = i;
	}

	// build element connectivity using local node indices (Tecplot uses 1-based indices)
	std::vector<std::vector<int>> conn_hex8;
	std::vector<std::vector<int>> conn_penta6;
	std::vector<std::vector<int>> conn_tet4;
	std::vector<std::vector<int>> conn_tri3;
	std::vector<std::vector<int>> conn_seg2;

	const bool use_element_selection = !selected_element_list.empty();
	if (use_element_selection)
	{
		for (int ei = 0; ei < (int)selected_element_list.size(); ++ei)
		{
			const int element_index = selected_element_list[ei];
			FEElement* pe = mesh.Element(element_index);
			if (pe == nullptr) continue;

			const int nn = pe->Nodes();
			std::vector<int> conn(nn);
			bool ok = true;
			for (int k = 0; k < nn; ++k)
			{
				const int gid = pe->m_node[k];
				if (gid < 0 || gid >= nNodes) { ok = false; break; }
				const int lid = global_to_local[gid];
				if (lid < 0) { ok = false; break; }
				conn[k] = lid + 1;
			}
			if (!ok) continue;

			switch (nn)
			{
			case 8: conn_hex8.push_back(std::move(conn)); break;
			case 6: conn_penta6.push_back(std::move(conn)); break;
			case 4: conn_tet4.push_back(std::move(conn)); break;
			case 3: conn_tri3.push_back(std::move(conn)); break;
			case 2: conn_seg2.push_back(std::move(conn)); break;
			default:
				// unsupported element type for Tecplot FE zone in this exporter
				break;
			}
		}
	}

    auto write_zone = [&](int it_step, const double t, const char* zone_suffix, const char* zone_type, const std::vector<std::vector<int>>& conn, int nodes_per_elem)
	{
		if (conn.empty()) return;

		out << "ZONE T=\"t=" << t << " " << zone_suffix << "\""
			<< ", N=" << output_nodes
			<< ", E=" << static_cast<int>(conn.size())
			<< ", ZONETYPE=" << zone_type
			<< ", DATAPACKING=POINT\n";

		for (int local_i = 0; local_i < output_nodes; ++local_i)
		{
			const int gid = output_node_ids[local_i];
            const double u = timeDisplacement[it_step][gid][0];
			const double v = timeDisplacement[it_step][gid][1];
			const double w = timeDisplacement[it_step][gid][2];

			const double x = initialCoordinate[gid][0] + u;
			const double y = initialCoordinate[gid][1] + v;
			const double z = initialCoordinate[gid][2] + w;

			out << x << " " << y << " " << z << " "
				<< u << " " << v << " " << w << " "
				<< t << " " << (gid + 1) << "\n";
		}

		for (const auto& econn : conn)
		{
			if ((int)econn.size() != nodes_per_elem) continue;
			for (int k = 0; k < nodes_per_elem; ++k)
			{
				out << econn[k];
				out << ((k == nodes_per_elem - 1) ? '\n' : ' ');
			}
		}
	};

	for (int it = start_index; it <= end_index_inclusive; ++it)
	{
		const double t = time[it];

		// If no element selection provided (or none matched), still output a node-only zone
		if (!use_element_selection || (conn_hex8.empty() && conn_penta6.empty() && conn_tet4.empty() && conn_tri3.empty() && conn_seg2.empty()))
		{
			out << "ZONE T=\"t=" << t << " scatter\""
				<< ", I=" << output_nodes
				<< ", J=1, K=1"
				<< ", DATAPACKING=POINT\n";

			for (int local_i = 0; local_i < output_nodes; ++local_i)
			{
				const int gid = output_node_ids[local_i];
				const double u = timeDisplacement[it][gid][0];
				const double v = timeDisplacement[it][gid][1];
				const double w = timeDisplacement[it][gid][2];

				const double x = initialCoordinate[gid][0] + u;
				const double y = initialCoordinate[gid][1] + v;
				const double z = initialCoordinate[gid][2] + w;

				out << x << " " << y << " " << z << " "
					<< u << " " << v << " " << w << " "
					<< t << " " << (gid + 1) << "\n";
			}
			continue;
		}

     write_zone(it, t, "HEX8", "FEBRICK", conn_hex8, 8);
		write_zone(it, t, "PENTA6", "FEPRISM", conn_penta6, 6);
		write_zone(it, t, "TET4", "FETETRAHEDRON", conn_tet4, 4);
		write_zone(it, t, "TRI3", "FETRIANGLE", conn_tri3, 3);
		write_zone(it, t, "SEG2", "FELINESEG", conn_seg2, 2);
	}

	return true;
}
