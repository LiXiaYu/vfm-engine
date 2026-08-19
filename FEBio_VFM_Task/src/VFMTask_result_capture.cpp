#include "VFMTask_result_capture.h"
#include "VFMTask.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <tuple>
#include <vector>
bool initialize_result_storage(FEModel* fem, VFMTask* task)
{
	task->recorded_steps = 0;

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

		// Guard against division by zero.
		if (dtmin <= 0) {
			// Fall back to m_ntime when dtmin is not configured.
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
	size_t constraintactivate_size = task->total_steps * task->nconstraint * sizeof(uint8_t); // task->timeconstraintactivate;

	size_t filesize = time_size + displacement_size + stress_size + nodalforce_size + constraintpressure_size + constraintactivate_size;

	if (!task->configure.reuse_saved_result_buffer)
	{
		write_to_log_2(fem, "create new result buffer...\n", outFile);
		// Start with a new buffer. Solve callbacks or the configured setters can
		// populate it later without changing its storage format.
		if (std::filesystem::exists(task->dumpfile)) {
			write_to_log_2(fem, "result buffer exists, replace it: " + task->dumpfile + "\n", outFile);
			std::filesystem::remove(task->dumpfile);
		}
		std::ofstream ofs(task->dumpfile, std::ios::binary | std::ios::trunc);
		ofs.seekp(filesize - 1);
		ofs.write("", 1);
		ofs.close();
	}
	else
	{
		write_to_log_2(fem, "reuse saved result buffer...\n", outFile);
		// Reuse means that an existing buffer is required. Whether Solve runs is
		// controlled independently by run_febio_solve.
		if (!std::filesystem::exists(task->dumpfile)) {
			write_to_log_2(fem, "saved result buffer does not exist: " + task->dumpfile + "\n", outFile);
			throw std::runtime_error("saved result buffer does not exist: " + task->dumpfile);
		}

		//size_t existing_filesize = std::filesystem::file_size(task->dumpfile);
		//if (existing_filesize != filesize) {
		//	std::stringstream ss;
		//	ss << "dump file size mismatch: expected " << filesize << ", got " << existing_filesize;
		//	write_to_log_2(fem, ss.str() + "\n", outFile);
		//	throw std::runtime_error(ss.str());
		//}
	}

	// Map the dump through mio.
	std::error_code ec;
	task->mmap.map(task->dumpfile, 0, mio::map_entire_file, ec);
	if (ec) {
		throw std::runtime_error("mmap failed: " + ec.message());
	}

	char* base = task->mmap.data();

	// Create typed views over the mapped storage.
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

	return everytimestep_withinited_savedata(fem, CB_INIT, task);
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

	if (!task->configure.reuse_saved_result_buffer)
	{
		if (task->recorded_steps >= task->total_steps)
		{
			write_to_log_2(fem, "Recorded timestep count exceeds the allocated result buffer.\n", outFile);
			return false;
		}

		const size_t index_timestep = task->recorded_steps;

		// get timestep
		double currenttime = fem->GetCurrentTime();
		::std::string logs_ct = "currenttime " + to_string(currenttime) + "\n";
		write_log(fem, 0, logs_ct.c_str());
		outFile << logs_ct;

		task->configure.timestep[index_timestep] = currenttime;

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

		task->recorded_steps++;
	}
	else
	{

	}


	outFile.close();

	return true;
}
