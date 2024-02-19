#include "VFMTask.h"

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <filesystem>

#include <signal.h>

#include <omp.h>

#include <nlopt.hpp>

#include "optim.h"
#include "VectorNum.h"

#include "FEBio_refunction.h"

#include "unilang_import.h"

VFMTask::VFMTask(FEModel* pfem):FECoreTask(pfem)
{
}

VFMTask::~VFMTask()
{

}

bool read_inited_information(FEModel* fem, unsigned int when, void* pd)
{
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

	if (task->isRead_FEMresult_fromsavefile == false)
	{
		// get timestep
		double currenttime = fem->GetCurrentTime();
		::std::string logs_ct = "currenttime " + to_string(currenttime) + "\n";
		write_log(fem, 0, logs_ct.c_str());
		outFile << logs_ct;

		task->timestep.push_back(currenttime);
		task->timedisplacement.push_back(::std::vector<::std::vector<double>>()); //ux, uy, uz
		task->timestress.push_back(::std::vector<::std::vector<double>>()); //sx, sy, sz, sxy, syz, szx

		// get displacement
		DataStore& datastore = fem->GetDataStore();
		int data_number = datastore.Size();
		::std::string logs_dn = "data_number " + to_string(data_number) + "\n";
		write_log(fem, 0, logs_dn.c_str());
		outFile << logs_dn;

		//data_number = 1; // only get the first data
		for (int j = 0; j < data_number; j++)
		{
			if (j == task->displacementdataNumber)
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
					::std::vector<double> u(nd);
					for (int k = 0; k < nd; ++k)
					{
						double val = datarecord.Evaluate(datarecord.m_item[i], k);
						//std::stringstream ss;
						//ss << std::setprecision(12) << val;
						//::std::string logs = "item " + to_string(i) + " xyz " + to_string(k) + ":" + ss.str() + "\n";
						//write_log(fem, 0, logs.c_str());
						//outFile << logs;

						u[k] = val;
					}
					(task->timedisplacement.end() - 1)->push_back(u);
				}
			}
			else if (j == task->stressdataNumber)
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
					::std::vector<double> u(nd);
					for (int k = 0; k < nd; ++k)
					{
						double val = datarecord.Evaluate(datarecord.m_item[i], k);
						//std::stringstream ss;
						//ss << std::setprecision(12) << val;
						//::std::string logs = "item " + to_string(i) + " xyz " + to_string(k) + ":" + ss.str() + "\n";
						//write_log(fem, 0, logs.c_str());
						//outFile << logs;

						u[k] = val;
					}
					(task->timestress.end() - 1)->push_back(u);
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

	std::vector<double> timeArray;
	::std::vector<int> solution_elementsDomainID;
	::std::vector<::std::vector<::std::vector<double>>> trueJArray;
	::std::vector<::std::vector<::std::vector<mat3d>>> truedeformationGradientArray;
	::std::vector<std::vector<std::vector<std::vector<mat3ds>>>> virtualstrainArrayV;
	::std::vector<::std::vector<double>> externalVirtualWork;
	::std::vector<::std::vector<double>> volumeVirtualWork;

	::std::function<double(const ::std::vector<double>&)> function;

	if (task->isReadfromsaveOptimfunc == false)
	{

#pragma region getInput
	// get mesh's object points
	FEMesh& mesh = fem->GetMesh();

	write_to_log_2(fem,"read initial Coordinate...\n", outFile);

	int nodes_number = mesh.Nodes();

	task->nodes.resize(nodes_number);
	task->initialCoordinate.resize(nodes_number);

	// This parallel for loop is not necessary, because it is not time-consuming
	//#pragma omp parallel for
	for (int j = 0; j < nodes_number; j++)
	{
		vec3d& node = mesh.Node(j).m_r0;
		//node.x;
		//node.y;
		//node.z;
		task->nodes[j] = node;
		task->initialCoordinate[j] = ::std::vector<double>{ node.x, node.y, node.z };
	}

	write_to_log_2(fem, "get solutions element...\n", outFile);
	// get solutions' element
	for (int i = 0; i < task->solution.size(); i++)
	{
		::std::string solution_type = ::std::get<0>(task->solution[i]);
		::std::string solution_name = ::std::get<1>(task->solution[i]);

		if (solution_type == "elements")
		{
			for (int j = 0; j < mesh.ElementSets(); j++)
			{
				FEElementSet& elementset = mesh.ElementSet(j);

				::std::string elementset_name = elementset.GetName();
				if (elementset_name == solution_name)
				{
					auto toadd = elementset.GetElementIDList();
					// foreach toadd element -=1
					for (int k = 0; k < toadd.size(); k++)
					{
						toadd[k] -= 1;
					}

					task->solution_elementsID.insert(task->solution_elementsID.end(), toadd.begin(), toadd.end());
				}
			}
		}

	}

	write_to_log_2(fem, "get fixed nodes...\n", outFile);
	// get fixed nodes
	for (int i = 0; i < task->fixed.size(); i++)
	{
		::std::string fixed_type = ::std::get<0>(task->fixed[i]);
		::std::string fixed_name = ::std::get<1>(task->fixed[i]);

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

						task->fixednode.push_back(node_id);
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

					task->fixednode.push_back(node_id);
				}
				
			}
		}
	}

#pragma endregion

#pragma region getAllInput
	//number of cycles considered at the end
	int numCycle = 2;
	int nNode = nodes_number;


	//::std::vector<int> sidenodes = { 12, 23, 56, 67, 66, 77, 22, 33 };
	//// sidenodes -= 1
	//for (int i = 0; i < sidenodes.size(); i++)
	//{
	//	sidenodes[i] -= 1;
	//}

	//// joint task->fixednode and sidenodes
	//task->fixednode.insert(task->fixednode.end(), sidenodes.begin(), sidenodes.end());



	//%% displacement at the initial timestep of the last two cycles
	#pragma region readsavefile_dump

	write_to_log_2(fem, "read from savefile dump...\n", outFile);
	if (task->isRead_FEMresult_fromsavefile == true)
	{
		DumpFile dumpfile(*fem);
		dumpfile.Open(task->dumpfile.c_str());

		dumpfile& task->timestep& task->timedisplacement& task->timestress;

		dumpfile.Close();

	}
	else
	{
		DumpFile dumpfile(*fem);
		dumpfile.Create(task->dumpfile.c_str());

		dumpfile& task->timestep& task->timedisplacement& task->timestress;

		dumpfile.Close();
	}

	// log: read save file end
	write_to_log_2(fem, "read save file end\n", outFile);

	#pragma endregion


	FEAnalysis& laststep = *(fem->GetStep(fem->Steps() - 1));

	// log: calculate timestep for VFM begin...
	write_to_log_2(fem, "calculate timestep for VFM begin...\n", outFile);

	double total_cycle_time = 60 / task->BMP;
	double sum_cycle_time = 0.0;
	int start_index = 0;
	for (int i = task->timestep.size() - 1; i >= 0; i--)
	{
		sum_cycle_time += task->timestep[i] - task->timestep[i - 1];

		if (sum_cycle_time >= (total_cycle_time* numCycle))
		{
			start_index= i;
			break;
		}
	}

	// log: calculate timestep for VFM end...
	write_to_log_2(fem, "calculate timestep for VFM end...\n", outFile);
	// log: start_index, sum_sycle_time
	write_to_log_2(fem, "start_index: " + to_string(start_index) + " sum_cycle_time: " + to_string(sum_cycle_time) + "\n", outFile);

	//int steps_per_cycle = 10; // 
	//int nStep_total = task->timestep.size();

	//int start_index = nStep_total - steps_per_cycle * numCycle;



	// copy ::std::vector task->timedisplacement[nNode*(start_index-1)+1:nNode*start_index] as initialDisp
	::std::vector<::std::vector<double>> initialDisp(task->timedisplacement[start_index-1]);
	::std::vector<::std::vector<::std::vector<double>>> displacementArray(task->timedisplacement.begin() + start_index, task->timedisplacement.end());
	timeArray=::std::vector<double>(task->timestep.begin() + start_index, task->timestep.end());

	::std::vector<::std::vector<double>> initialStress(task->timestress[start_index-1]);
	::std::vector<::std::vector<::std::vector<double>>> stressArray(task->timestress.begin() + start_index, task->timestress.end());

	::std::vector<::std::vector<::std::vector<double>>> velocityArray(displacementArray);//(task->timedisplacement.begin() + start_index, task->timedisplacement.end());
	::std::vector<::std::vector<::std::vector<double>>> accelerationArray(displacementArray);//(task->timedisplacement.begin() + start_index, task->timedisplacement.end());

	task->timestep.clear();
	task->timedisplacement.clear();
	task->timestress.clear();

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
	for (int index_initialCoordinate = 0; index_initialCoordinate < task->initialCoordinate.size(); index_initialCoordinate++)
	{
		for (int index_displacement = 0; index_displacement < task->initialCoordinate[index_initialCoordinate].size(); index_displacement++)
		{
			task->initialCoordinate[index_initialCoordinate][index_displacement] += initialDisp[index_initialCoordinate][index_displacement];
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


	task->nVirtualFields = 4;

	#pragma region getAllVirtualFields

	::std::string logs_vf = "Start get all Virtual Fields.\n";
	write_log(fem, 0, logs_vf.c_str());
	outFile << logs_vf;

	::std::vector<::std::function<::std::vector<double>(const ::std::vector<double>&, int )>> vf_u_functions;

	auto radius = 0.8;
	auto height = 0.4;

	double x_r_0 = 0.0;
	double y_r_0 = 0.0;
	double z_r_0 = 0.0;

	double x_r_t = 1.0;
	double y_r_t = 1.0;
	double z_r_t = 1.0;

	// BPM_base_v4_BPM
		#pragma region BPM_base_v4_BPM
	radius = 1.6;
	height = 1;

	x_r_0 = 142.185;
	y_r_0 = 141.782;
	z_r_0 = 0;
	x_r_t = 1.0;
	y_r_t = 1.0;
	z_r_t = 1.0;

	double x_min = 140.185;
	double x_max = 144.185;
	double y_min = 141.289;
	double y_max = 143.289;
	double z_min = -2.0;
	double z_max = 2.0;

	double x_p = 1;
	double y_p = 1;
	double z_p = 1;

	auto it = ::std::remove_if(task->solution_elementsID.begin(), task->solution_elementsID.end(),
		[=, &mesh, &task](int j)->bool {
			FEElement& element = *(mesh.Element(j));

			::std::vector<int> no_inrange_point_number;
			for (int k = 0; k < element.Nodes(); k++)
			{
				int node_id = element.m_node[k];
				auto& node = mesh.Node(node_id);

				double x = task->initialCoordinate[node_id][0];
				double y = task->initialCoordinate[node_id][1];
				double z = task->initialCoordinate[node_id][2];

				double x_r = (x - x_r_0) / x_r_t;
				double y_r = (y - y_r_0) / y_r_t;
				double z_r = (z - z_r_0) / z_r_t;

				double rho = ::std::sqrt(x_r * x_r + z_r * z_r);
				double theta = ::std::atan2(z_r, x_r);


				// if x,y,z not in range, remove this solution element
				if (rho > 0.5 || y < y_min || y > y_max)
				{
					no_inrange_point_number.push_back(k);
				}
			}

			if (no_inrange_point_number.size() == element.Nodes())
			{
				return true;
			}
			else
			{
				// add node id to fixednode
				for (int k = 0; k < no_inrange_point_number.size(); k++)
				{
					task->fixednode.push_back(element.m_node[no_inrange_point_number[k]]);
				}
				return false;
			}
		});
	task->solution_elementsID.erase(it, task->solution_elementsID.end());

	solution_elementsDomainID = ::std::vector<int>(task->solution_elementsID.size());
	#pragma omp parallel for
	for (int index_seId = 0; index_seId < task->solution_elementsID.size(); index_seId++)
	{
		int elementid = task->solution_elementsID[index_seId];

		int domain_id = -1;
		for (int i = 0; i < mesh.Domains(); i++)
		{
			FEDomain& d = mesh.Domain(i);
			auto* elementfound = d.FindElementFromID(elementid);
			if (elementfound != nullptr)
			{
				domain_id = i;
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
		elementset_ss << (elementid+1);
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
	for (int index_fnId = 0; index_fnId < task->fixednode.size(); index_fnId++)
	{
		int nodeid = task->fixednode[index_fnId];
		fixednodeset_ss << (nodeid+1);
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

		#pragma region initialDeformationGradient
	// set initial displacement to nodes
	write_to_log_2(fem, "Set initial displacement to nodes.\n", outFile);

	#pragma omp parallel for
	for (int index_seId = 0; index_seId < task->solution_elementsID.size(); index_seId++)
	{
		int j = task->solution_elementsID[index_seId];

		FEElement& element = *(mesh.Element(j));

		for (int k = 0; k < element.Nodes(); k++)
		{
			int node_id = element.m_node[k];
			auto& node = mesh.Node(node_id);
			node.m_r0 = { task->initialCoordinate[node_id][0], task->initialCoordinate[node_id][1], task->initialCoordinate[node_id][2] };
			node.m_rt = { task->initialCoordinate[node_id][0], task->initialCoordinate[node_id][1], task->initialCoordinate[node_id][2] };
			node.m_rp = { 0,0,0 };
			node.m_d0 = { 0,0,0 };
			node.m_dt = { 0,0,0 };
			node.m_dp = { 0,0,0 };

		}

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
		::std::vector<::std::vector<double>> currentCoordinate(task->initialCoordinate);
		#pragma omp parallel for
		for (int index_coordinate_i = 0; index_coordinate_i < currentCoordinate.size(); index_coordinate_i++)
		{
			for (int index_coordinate_j = 0; index_coordinate_j < currentCoordinate[index_coordinate_i].size(); index_coordinate_j++)
			{
				currentCoordinate[index_coordinate_i][index_coordinate_j] += displacementArray[index_timestep][index_coordinate_i][index_coordinate_j];
			}
		}
		#pragma omp parallel for
		for (int index_seId = 0; index_seId < task->solution_elementsID.size(); index_seId++)
		{
			int j = task->solution_elementsID[index_seId];

			FEElement& element = *(mesh.Element(j));

			for (int k = 0; k < element.Nodes(); k++)
			{
				int node_id = element.m_node[k];
				auto& node = mesh.Node(node_id);
				node.m_r0 = { task->initialCoordinate[node_id][0], task->initialCoordinate[node_id][1], task->initialCoordinate[node_id][2] }; // true displacement as initial displacement.
				node.m_d0 = { 0,0,0 };
				node.m_rt = { currentCoordinate[node_id][0], currentCoordinate[node_id][1], currentCoordinate[node_id][2] }; // move to virtual displacement
				node.m_dt = { 0,0,0 };
				node.m_rp = { 0,0,0 };
				node.m_dp = { 0,0,0 };
			}

			FESolidDomain& domain = static_cast<FESolidDomain&>(mesh.Domain(solution_elementsDomainID[index_seId]));

			domain_init(mesh, static_cast<FESolidElement&>(element), domain);
		}

		trueJArray[index_timestep] = ::std::vector<::std::vector<double>>(task->solution_elementsID.size());
		truedeformationGradientArray[index_timestep] = ::std::vector<::std::vector<mat3d>>(task->solution_elementsID.size());

		#pragma omp parallel for
		for (int index_seId = 0; index_seId < task->solution_elementsID.size(); index_seId++)
		{
			int j = task->solution_elementsID[index_seId];

			FEElement& element = *(mesh.Element(j));
			int nodes_number = element.Nodes();

			// convert element to solid element
			FESolidElement& solidElement = static_cast<FESolidElement&>(element);

			// get domain
			FESolidDomain& domain = static_cast<FESolidDomain&>(mesh.Domain(solution_elementsDomainID[index_seId]));

			trueJArray[index_timestep][index_seId]=::std::vector<double>(solidElement.GaussPoints());
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
						throw e; // don't continue execution!!!!
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

	auto& task_fixednode = task->fixednode;

		#pragma region ONH_sinPulse100s_120_v4
	//vf_u_functions[0] = [&task_fixednode, radius, height, x_r_0, y_r_0, z_r_0, x_r_t, y_r_t, z_r_t](const ::std::vector<double>& coordinate, int index_coordinate) -> ::std::vector<double>
	//{
	//	auto x = coordinate[0];
	//	auto y = coordinate[1];
	//	auto z = coordinate[2];

	//	x = (x - x_r_0) * x_r_t;
	//	y = (y - y_r_0) * y_r_t;
	//	z = (z - z_r_0) * z_r_t;

	//	::std::vector<double> vf_u(3);
	//	vf_u[0] = (radius - ::std::abs(x)) * x / height; //2d
	//	//vf_u[0] = (radius - ::std::abs(x)) * x / height / ::std::sqrt(2); // BPM_base_v4_BPM
	//	vf_u[1] = -1 * (radius - ::std::abs(x)) * ::std::abs(y) / height;
	//	vf_u[2] = 0; //2d
	//	//vf_u[2] = (radius - ::std::abs(z)) * z / height / ::std::sqrt(2); // BPM_base_v4_BPM


	//	// find fixed node and set virtual field to zero
	//	if (::std::find(task_fixednode.begin(), task_fixednode.end(), index_coordinate) != task_fixednode.end())
	//	{
	//		vf_u[0] = 0;
	//		vf_u[1] = 0;
	//		vf_u[2] = 0;
	//	}

	//	return vf_u;
	//};
	//vf_u_functions[1] = [&task_fixednode, radius, height, x_r_0, y_r_0, z_r_0, x_r_t, y_r_t, z_r_t](const ::std::vector<double>& coordinate, int index_coordinate) -> ::std::vector<double>
	//{
	//	auto x = coordinate[0];
	//	auto y = coordinate[1];
	//	auto z = coordinate[2];

	//	x = (x - x_r_0) * x_r_t;
	//	y = (y - y_r_0) * y_r_t;
	//	z = (z - z_r_0) * z_r_t;

	//	::std::vector<double> vf_u(3);

	//	vf_u[0] = 0;
	//	vf_u[1] = (radius - ::std::abs(x)) * ::std::abs(y) / height;
	//	vf_u[2] = 0;

	//	// find fixed node and set virtual field to zero
	//	if (::std::find(task_fixednode.begin(), task_fixednode.end(), index_coordinate) != task_fixednode.end())
	//	{
	//		vf_u[0] = 0;
	//		vf_u[1] = 0;
	//		vf_u[2] = 0;
	//	}

	//	return vf_u;
	//};
	//vf_u_functions[2] = [&task_fixednode, radius, height, x_r_0, y_r_0, z_r_0, x_r_t, y_r_t, z_r_t](const ::std::vector<double>& coordinate, int index_coordinate) -> ::std::vector<double>
	//{
	//	auto x = coordinate[0];
	//	auto y = coordinate[1];
	//	auto z = coordinate[2];

	//	x = (x - x_r_0) * x_r_t;
	//	y = (y - y_r_0) * y_r_t;
	//	z = (z - z_r_0) * z_r_t;

	//	::std::vector<double> vf_u(3);

	//	vf_u[0] = ::std::cos(0.5 * PI * ::std::abs(x) / radius) * x / height; //2d
	//	//vf_u[0] = ::std::cos(0.5 * PI * ::std::abs(x) / radius) * x / height / ::std::sqrt(2); // BPM_base_v4_BPM
	//	vf_u[1] = -1 * ::std::cos(0.5 * PI * ::std::abs(x) / radius) * ::std::abs(y) / height;
	//	vf_u[2] = 0; //2d
	//	//vf_u[2] = ::std::cos(0.5 * PI * ::std::abs(z) / radius) * z / height / ::std::sqrt(2); // BPM_base_v4_BPM

	//	// find fixed node and set virtual field to zero
	//	if (::std::find(task_fixednode.begin(), task_fixednode.end(), index_coordinate) != task_fixednode.end())
	//	{
	//		vf_u[0] = 0;
	//		vf_u[1] = 0;
	//		vf_u[2] = 0;
	//	}

	//	return vf_u;
	//};
	//vf_u_functions[3] = [&task_fixednode, radius, height, x_r_0, y_r_0, z_r_0, x_r_t, y_r_t, z_r_t](const ::std::vector<double>& coordinate, int index_coordinate) -> ::std::vector<double>
	//{
	//	auto x = coordinate[0];
	//	auto y = coordinate[1];
	//	auto z = coordinate[2];

	//	x = (x - x_r_0) * x_r_t;
	//	y = (y - y_r_0) * y_r_t;
	//	z = (z - z_r_0) * z_r_t;

	//	::std::vector<double> vf_u(3);

	//	vf_u[0] = 0;
	//	vf_u[1] = ::std::cos(0.5 * PI * abs(x) / radius) * ::std::abs(y) / height;
	//	vf_u[2] = 0;

	//	// find fixed node and set virtual field to zero
	//	if (::std::find(task_fixednode.begin(), task_fixednode.end(), index_coordinate) != task_fixednode.end())
	//	{
	//		vf_u[0] = 0;
	//		vf_u[1] = 0;
	//		vf_u[2] = 0;
	//	}

	//	return vf_u;
	//};

		#pragma endregion



	auto& task_initialCoordinate = task->initialCoordinate;

	const int CONST_VF_select = 1;// up

	if constexpr (CONST_VF_select == 0) // up
	{
		vf_u_functions.push_back([&task_fixednode, &task_initialCoordinate](const ::std::vector<double>& coordinate, int index_coordinate) -> ::std::vector<double>
			{
				::std::vector<double> vf_u(3);

				vf_u[0] = 0.0001;
				vf_u[1] = 0;
				vf_u[2] = 0;

				// find fixed node and set virtual field to zero
				if (::std::find(task_fixednode.begin(), task_fixednode.end(), index_coordinate) != task_fixednode.end())
				{
					vf_u[0] = 0;
					vf_u[1] = 0;
					vf_u[2] = 0;
				}

				return vf_u;
			});
		vf_u_functions.push_back([&task_fixednode, &task_initialCoordinate](const ::std::vector<double>& coordinate, int index_coordinate) -> ::std::vector<double>
			{
				::std::vector<double> vf_u(3);

				vf_u[0] = 0;
				vf_u[1] = 0.0001;
				vf_u[2] = 0;

				// find fixed node and set virtual field to zero
				if (::std::find(task_fixednode.begin(), task_fixednode.end(), index_coordinate) != task_fixednode.end())
				{
					vf_u[0] = 0;
					vf_u[1] = 0;
					vf_u[2] = 0;
				}

				return vf_u;
			});
		vf_u_functions.push_back([&task_fixednode, &task_initialCoordinate](const ::std::vector<double>& coordinate, int index_coordinate) -> ::std::vector<double>
			{
				::std::vector<double> vf_u(3);

				vf_u[0] = 0;
				vf_u[1] = 0;
				vf_u[2] = 0.0001;

				// find fixed node and set virtual field to zero
				if (::std::find(task_fixednode.begin(), task_fixednode.end(), index_coordinate) != task_fixednode.end())
				{
					vf_u[0] = 0;
					vf_u[1] = 0;
					vf_u[2] = 0;
				}

				return vf_u;
			});
	}
	else
	{
		const bool CONST_bpm_base_v4_true_disp_vf_u = false;
		if constexpr (CONST_bpm_base_v4_true_disp_vf_u == true)
		{

			vf_u_functions.push_back([&task_fixednode, &task_initialCoordinate](const ::std::vector<double>& coordinate, int index_coordinate) -> ::std::vector<double>
				{
					::std::vector<double> vf_u(3);

					vf_u[0] = coordinate[0] - task_initialCoordinate[index_coordinate][0];
					vf_u[1] = coordinate[1] - task_initialCoordinate[index_coordinate][1];
					vf_u[2] = coordinate[2] - task_initialCoordinate[index_coordinate][2];

					// find fixed node and set virtual field to zero
					if (::std::find(task_fixednode.begin(), task_fixednode.end(), index_coordinate) != task_fixednode.end())
					{
						vf_u[0] = 0;
						vf_u[1] = 0;
						vf_u[2] = 0;
					}

					return vf_u;
				});
		}
		else
		{
#pragma region BPM_base_v4_BPM
			const bool CONST_bpm_base_v4_only_one_vf_u = false;
			if constexpr (CONST_bpm_base_v4_only_one_vf_u == false)
			{
				vf_u_functions.push_back([&task_fixednode, radius, height, x_r_0, y_r_0, z_r_0, x_r_t, y_r_t, z_r_t, x_p, y_p, z_p](const ::std::vector<double>& coordinate, int index_coordinate) -> ::std::vector<double>
					{
						auto x = coordinate[0];
						auto y = coordinate[1];
						auto z = coordinate[2];

						x = (x - x_r_0) * x_r_t;
						y = (y - y_r_0) * y_r_t;
						z = (z - z_r_0) * z_r_t;

						double r = ::std::sqrt(x * x + z * z);
						double w = ::std::atan2(z, x);
						double h = y;

						double ur = ::std::sin(2 * PI * r / radius) * h / height;
						double uw = 0;
						double uh = 0;

						::std::vector<double> vf_u(3);
						vf_u[0] = (r + ur) * ::std::cos(w + uw) - x;
						vf_u[1] = uh;
						vf_u[2] = (r + ur) * ::std::sin(w + uw) - z;

						vf_u[0] *= x_p;
						vf_u[1] *= y_p;
						vf_u[2] *= z_p;

						// find fixed node and set virtual field to zero
						if (::std::find(task_fixednode.begin(), task_fixednode.end(), index_coordinate) != task_fixednode.end())
						{
							vf_u[0] = 0;
							vf_u[1] = 0;
							vf_u[2] = 0;
						}

						return vf_u;
					});
				vf_u_functions.push_back([&task_fixednode, radius, height, x_r_0, y_r_0, z_r_0, x_r_t, y_r_t, z_r_t, x_p, y_p, z_p](const ::std::vector<double>& coordinate, int index_coordinate) -> ::std::vector<double>
					{
						auto x = coordinate[0];
						auto y = coordinate[1];
						auto z = coordinate[2];

						x = (x - x_r_0) * x_r_t;
						y = (y - y_r_0) * y_r_t;
						z = (z - z_r_0) * z_r_t;

						double r = ::std::sqrt(x * x + z * z);
						double w = ::std::atan2(z, x);
						double h = y;

						double ur = r * (radius - r) * (radius / 2 - r) * h / height;
						double uw = 0;
						double uh = 0;

						::std::vector<double> vf_u(3);
						vf_u[0] = (r + ur) * ::std::cos(w + uw) - x;
						vf_u[1] = uh;
						vf_u[2] = (r + ur) * ::std::sin(w + uw) - z;

						vf_u[0] *= x_p;
						vf_u[1] *= y_p;
						vf_u[2] *= z_p;

						// find fixed node and set virtual field to zero
						if (::std::find(task_fixednode.begin(), task_fixednode.end(), index_coordinate) != task_fixednode.end())
						{
							vf_u[0] = 0;
							vf_u[1] = 0;
							vf_u[2] = 0;
						}

						return vf_u;
					});
				vf_u_functions.push_back([&task_fixednode, radius, height, x_r_0, y_r_0, z_r_0, x_r_t, y_r_t, z_r_t, x_p, y_p, z_p](const ::std::vector<double>& coordinate, int index_coordinate) -> ::std::vector<double>
					{
						auto x = coordinate[0];
						auto y = coordinate[1];
						auto z = coordinate[2];

						x = (x - x_r_0) * x_r_t;
						y = (y - y_r_0) * y_r_t;
						z = (z - z_r_0) * z_r_t;

						double r = ::std::sqrt(x * x + z * z);
						double w = ::std::atan2(z, x);
						double h = y;

						double ur = 0;
						double uw = ::std::sin(2 * PI * r / radius) * h / height / r;
						double uh = 0;

						::std::vector<double> vf_u(3);
						vf_u[0] = (r + ur) * ::std::cos(w + uw) - x;
						vf_u[1] = uh;
						vf_u[2] = (r + ur) * ::std::sin(w + uw) - z;

						vf_u[0] *= x_p;
						vf_u[1] *= y_p;
						vf_u[2] *= z_p;

						// find fixed node and set virtual field to zero
						if (::std::find(task_fixednode.begin(), task_fixednode.end(), index_coordinate) != task_fixednode.end())
						{
							vf_u[0] = 0;
							vf_u[1] = 0;
							vf_u[2] = 0;
						}

						return vf_u;
					});
				vf_u_functions.push_back([&task_fixednode, radius, height, x_r_0, y_r_0, z_r_0, x_r_t, y_r_t, z_r_t, x_p, y_p, z_p](const ::std::vector<double>& coordinate, int index_coordinate) -> ::std::vector<double>
					{
						auto x = coordinate[0];
						auto y = coordinate[1];
						auto z = coordinate[2];

						x = (x - x_r_0) * x_r_t;
						y = (y - y_r_0) * y_r_t;
						z = (z - z_r_0) * z_r_t;

						double r = ::std::sqrt(x * x + z * z);
						double w = ::std::atan2(z, x);
						double h = y;

						double ur = 0;
						double uw = (radius - r) * (radius / 2 - r) * h / height;
						double uh = 0;

						::std::vector<double> vf_u(3);
						vf_u[0] = (r + ur) * ::std::cos(w + uw) - x;
						vf_u[1] = uh;
						vf_u[2] = (r + ur) * ::std::sin(w + uw) - z;

						vf_u[0] *= x_p;
						vf_u[1] *= y_p;
						vf_u[2] *= z_p;

						// find fixed node and set virtual field to zero
						if (::std::find(task_fixednode.begin(), task_fixednode.end(), index_coordinate) != task_fixednode.end())
						{
							vf_u[0] = 0;
							vf_u[1] = 0;
							vf_u[2] = 0;
						}

						return vf_u;
					});
			}
			vf_u_functions.push_back([&task_fixednode, radius, height, x_r_0, y_r_0, z_r_0, x_r_t, y_r_t, z_r_t, x_p, y_p, z_p](const ::std::vector<double>& coordinate, int index_coordinate) -> ::std::vector<double>
				{
					auto x = coordinate[0];
					auto y = coordinate[1];
					auto z = coordinate[2];

					x = (x - x_r_0) * x_r_t;
					y = (y - y_r_0) * y_r_t;
					z = (z - z_r_0) * z_r_t;

					double r = ::std::sqrt(x * x + z * z);
					double w = ::std::atan2(z, x);
					double h = y;

					double ur = 0;
					double uw = 0;
					double uh = ::std::sin(PI * (radius - r) / 2 / radius) * h / height;

					::std::vector<double> vf_u(3);
					vf_u[0] = (r + ur) * ::std::cos(w + uw) - x;
					vf_u[1] = uh;
					vf_u[2] = (r + ur) * ::std::sin(w + uw) - z;

					vf_u[0] *= x_p;
					vf_u[1] *= y_p;
					vf_u[2] *= z_p;

					// find fixed node and set virtual field to zero
					if (::std::find(task_fixednode.begin(), task_fixednode.end(), index_coordinate) != task_fixednode.end())
					{
						vf_u[0] = 0;
						vf_u[1] = 0;
						vf_u[2] = 0;
					}

					return vf_u;
				});
#pragma endregion
		}

	}

	#pragma endregion


	#pragma region deformationGradient

	// VirtualWork
	externalVirtualWork = ::std::vector<::std::vector<double>>(timeArray.size());

	volumeVirtualWork = ::std::vector<::std::vector<double>>(timeArray.size());

	write_to_log_2(fem, "Attemp to calculate virtual work.\n", outFile);

	virtualstrainArrayV = ::std::vector<::std::vector<::std::vector<::std::vector<mat3ds>>>>(timeArray.size()); // virtual strain : timestep, vf, element

	// foreach timestep
	for (int index_timestep = 0; index_timestep < timeArray.size(); index_timestep++)
	{
		write_to_log_2(fem, "timestep:" + ::std::to_string(index_timestep) + "\n", outFile);

		::std::vector<::std::vector<double>> currentCoordinate(task->initialCoordinate);
		#pragma omp parallel for
		for (int index_coordinate_i = 0; index_coordinate_i < currentCoordinate.size(); index_coordinate_i++)
		{
			for (int index_coordinate_j = 0; index_coordinate_j < currentCoordinate[index_coordinate_i].size(); index_coordinate_j++)
			{
				currentCoordinate[index_coordinate_i][index_coordinate_j] += displacementArray[index_timestep][index_coordinate_i][index_coordinate_j];
			}
		}

		// true stress for parameters
		write_to_log_2(fem, "get true stress.\n", outFile);

		::std::vector<::std::vector<double>>& currentStress=stressArray[index_timestep];

		::std::vector<double> ivw(vf_u_functions.size()); // internalVirtualWork
		::std::vector<double> evw(vf_u_functions.size()); // externalVirtualWork



		// virtual field
		write_to_log_2(fem, "Cal virtual field.\n", outFile);

		virtualstrainArrayV[index_timestep] = ::std::vector<::std::vector<::std::vector<mat3ds>>>(vf_u_functions.size());
		externalVirtualWork[index_timestep] = ::std::vector<double>(vf_u_functions.size());
		volumeVirtualWork[index_timestep] = ::std::vector<double>(vf_u_functions.size());
		for (int index_vf = 0; index_vf < vf_u_functions.size(); index_vf++)
		{
			write_to_log_2(fem, "virtual field " + ::std::to_string(index_vf) + "\n", outFile);

			write_to_log_2(fem, "set virtual displacement.\n", outFile);

			#pragma omp parallel for
			for (int index_seId = 0; index_seId < task->solution_elementsID.size(); index_seId++)
			{
				int j = task->solution_elementsID[index_seId];

				FEElement& element = *(mesh.Element(j));

				for (int k = 0; k < element.Nodes(); k++)
				{
					int node_id = element.m_node[k];
					auto& node = mesh.Node(node_id);

					auto vf_u = vf_u_functions[index_vf](currentCoordinate[node_id], node_id);
					//for (int i = 0; i < 3; i++)
					//{
					//	vf_u[i] -= currentCoordinate[node_id][i];
					//}

					//node.m_r0 = { task->initialCoordinate[node_id][0], task->initialCoordinate[node_id][1], task->initialCoordinate[node_id][2] };

					node.m_r0 = { currentCoordinate[node_id][0], currentCoordinate[node_id][1], currentCoordinate[node_id][2] }; // true displacement as initial displacement.
					node.m_d0 = { 0,0,0 };
					node.m_rt = { vf_u[0],vf_u[1], vf_u[2] }; // move to virtual displacement
					node.m_dt = { 0,0,0 };
					node.m_rp = { 0,0,0 };
					node.m_dp = { 0,0,0 };
				}

				FESolidDomain& domain = static_cast<FESolidDomain&>(mesh.Domain(solution_elementsDomainID[index_seId]));

				domain_init(mesh, static_cast<FESolidElement&>(element), domain);
			}
			
			// internal virtual work
			write_to_log_2(fem, "calculate internal virtual work.\n", outFile);



			virtualstrainArrayV[index_timestep][index_vf] = ::std::vector<::std::vector<mat3ds>>(task->solution_elementsID.size());

			#pragma omp parallel for
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
				#pragma omp parallel for
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

						double Jc;
						mat3d Ft, Fp;

						// need recalculate
						// calculate inverse jacobian
						Jc = trueJArray[index_timestep][index_seId][n];


						double J = pt.m_J;
						double J0 = mp.m_J0;
						pt.m_F = Ft;
						// pt.m_J = Jt;

						// strainCompute infinitesimal strain
						mat3ds infstrain = pt.m_F.sym();
						//mat3ds strain = pt.Strain(); // virtual strain

						mat3ds& s = infstrain;

						virtualstrainArrayV[index_timestep][index_vf][index_seId][n] = s;

						double density = dynamic_cast<FESolidMaterial*>(pmat)->Density(mp);
						density = 1;
						// acceleration
						vec3d a = solidElement.Evaluate(a0, n);
						vec3d u = solidElement.Evaluate(u0, n);

						double vvw = density * (a * u) * Jc;

						vvw_elements[index_seId][n] = vvw;
					}


				}

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

			for (int j = 0; j < task->pressure_load.size(); j++)
			{
				for (int i = 0; i < fem->ModelLoads(); i++)
				{
					FEModelLoad& load = *(fem->ModelLoad(i));
					::std::string loadclassname = load.GetFactoryClass()->GetClassName();
					::std::string loadname = load.GetName();
					
					if (::std::get<0>(task->pressure_load[j]) == loadclassname && ::std::get<1>(task->pressure_load[j]) == loadname)
					{
						if (loadclassname == "FEPressureLoad")
						{
							FEPressureLoad& pressureLoad = static_cast<FEPressureLoad&>(load);

							auto& surface = pressureLoad.GetSurface();

							::std::vector<::std::vector<double>> evw_surface(surface.Elements());
							::std::vector<bool> evw_surface_flag(surface.Elements(), false);

							#pragma omp parallel for
							for (int index_surface_element = 0; index_surface_element < surface.Elements(); index_surface_element++)
							{
								auto& element = surface.Element(index_surface_element);

								auto* true_element = element.m_elem[0];

								int elementId = true_element->GetID()-1;
								::std::vector<::std::reference_wrapper<FENode>> nodes;
								for (size_t i = 0; i < element.Nodes(); i++)
								{
									nodes.push_back(mesh.Node(element.m_node[i]));
								}

								evw_surface[index_surface_element] = ::std::vector<double>(element.GaussPoints(), 0);

								// judge surfaceelement's element is in solution_elementsID
								int index_seId = -1;
								if (::std::find(task->solution_elementsID.begin(), task->solution_elementsID.end(), elementId) != task->solution_elementsID.end())
								{
									vec3d re[FEElement::MAX_NODES];
									surface.GetReferenceNodalCoordinates(element, re);
									vec3d rv[FEElement::MAX_NODES];
									surface.NodalCoordinates(element, rv);


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


										double P = (currentStress[elementId][0] + currentStress[elementId][1] + currentStress[elementId][2]) / 3;
										double iP = (initialStress[elementId][0] + initialStress[elementId][1] + initialStress[elementId][2]) / 3;

										P -= iP;


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

										if (elementId == 35072 && n == 0)
										{
											int awer23dsfcsafwerf = 0;
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

							auto evw=externalVirtualWork[index_timestep][index_vf];
						}

					}
				}

			}

			

		}
	}



	#pragma endregion


	double forceMPa = 0.0004; // 3mmHg

#pragma endregion


		DumpFile dumpfile(*fem);
		dumpfile.Create((task->dumpfile + ".fundump").c_str());

		dumpfile& timeArray& (task->solution_elementsID) & solution_elementsDomainID& trueJArray& truedeformationGradientArray& virtualstrainArrayV& externalVirtualWork& volumeVirtualWork;

		dumpfile.Close();

		auto solution_elementsID(task->solution_elementsID);

		::std::function<double(const ::std::vector<double>&)> fun = [&fem, timeArray, &outFile, solution_elementsID, solution_elementsDomainID, trueJArray, truedeformationGradientArray, virtualstrainArrayV, externalVirtualWork, volumeVirtualWork](const ::std::vector<double>& ps) {
			double p_E = ps[0];
			double p_g = ps[1];
			double p_t = ps[2];

			double value = fun_for_optim(fem, p_g, p_t, p_E, timeArray, outFile, fem->GetMesh(), solution_elementsID, solution_elementsDomainID, trueJArray, truedeformationGradientArray, virtualstrainArrayV, externalVirtualWork, volumeVirtualWork);


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
			dumpfile.Open((task->dumpfile + ".fundump").c_str());

			dumpfile& timeArray& (task->solution_elementsID)& solution_elementsDomainID& trueJArray& truedeformationGradientArray& virtualstrainArrayV& externalVirtualWork& volumeVirtualWork;

			dumpfile.Close();

			auto solution_elementsID(task->solution_elementsID);

			::std::function<double(const ::std::vector<double>&)> fun = [&fem, timeArray, &outFile, solution_elementsID, solution_elementsDomainID, trueJArray, truedeformationGradientArray, virtualstrainArrayV, externalVirtualWork, volumeVirtualWork](const ::std::vector<double>& ps) {
				double p_E = ps[0];
				double p_g = ps[1];
				double p_t = ps[2];

				double value = fun_for_optim(fem, p_g, p_t, p_E, timeArray, outFile, fem->GetMesh(), solution_elementsID, solution_elementsDomainID, trueJArray, truedeformationGradientArray, virtualstrainArrayV, externalVirtualWork, volumeVirtualWork);


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
				write_to_log_2(fem, ("Reading... "+fundumpfile + "\n").c_str(), outFile);


				DumpFile dumpfile(*fem);
				dumpfile.Open(fundumpfile.c_str());

				dumpfile& timeArray& (task->solution_elementsID)& solution_elementsDomainID& trueJArray& truedeformationGradientArray& virtualstrainArrayV& externalVirtualWork& volumeVirtualWork;

				dumpfile.Close();

				auto solution_elementsID(task->solution_elementsID);

				funs.push_back([&fem,timeArray,&outFile, solution_elementsID, solution_elementsDomainID, trueJArray, truedeformationGradientArray, virtualstrainArrayV, externalVirtualWork, volumeVirtualWork](const ::std::vector<double>& ps) {
					double p_E = ps[0];
					double p_g = ps[1];
					double p_t = ps[2];

					double value = fun_for_optim(fem, p_g, p_t, p_E, timeArray, outFile, fem->GetMesh(), solution_elementsID, solution_elementsDomainID, trueJArray, truedeformationGradientArray, virtualstrainArrayV, externalVirtualWork, volumeVirtualWork);

					return value;
					});
			}
			function = [&fem,&outFile,funs](const ::std::vector<double>& ps) {
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

	#pragma region Opt10x10x10

	const bool CONST_Opt10x10x10 = true;
	if constexpr (CONST_Opt10x10x10 == true)
	{
		auto dop = [&](const ::std::vector<::std::vector<double>>& psss, ::std::function<double(const ::std::vector<double>&)> fun, bool ifoutlog = false, bool ifoutTacplot = false,::std::vector<int> contour_num={9,9,1}, ::std::string outTacplot_filepath = "./temp/debug/Tacplot.dat") {
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
					TacplotFile << ::std::setprecision(12)  << pss[0] << " " << pss[1] << " " << pss[2] << " " << value << "\n";
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
		auto grid = createGrid({ 0.1,1,1 }, { 0.9,9,1 }, { 0.1,1,0 });

		auto [min_value, min_ps, dop_map] = dop(grid, function, true, true);
		double min_E = min_ps[0];
		double min_g = min_ps[1];
		double min_t = min_ps[2];

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

		write_to_log_2(fem,"min_E: "+::std::to_string(min_E)+" min_g: "+::std::to_string(min_g)+" min_t: "+::std::to_string(min_t)+"\n", outFile);

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
		auto grid2 = createGrid({ ::std::max(min_E1 - 0.02,0.1),::std::max(min_g1 - 0.2,1.0),1.0 }, { ::std::min(min_E1 + 0.02,0.9),::std::min(min_g1 + 0.2,9.0),1.0}, { 0.004,0.04,0 });

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
	const bool CONST_NLpot_0 = false;
	if constexpr (CONST_NLpot_0 == true)
	{
		write_to_log_2(fem, "Start NLpot_0 :\n", outFile);

		//nlopt::opt opt(nlopt::LN_COBYLA, 3);
		nlopt::opt opt(nlopt::GN_DIRECT, 3);
		//nlopt::opt opt(nlopt::LN_NELDERMEAD, 3);
		//nlopt::opt opt(nlopt::GN_ISRES, 3);
		::std::vector<double> lb(3);
		lb[0] = 0.2; lb[1] = 6.0; lb[2] = 1.0;
		opt.set_lower_bounds(lb);
		::std::vector<double> ub(3);
		ub[0] = 0.4; ub[1] = 8.0; ub[2] = 2.0;
		opt.set_upper_bounds(ub);

		opt.set_min_objective(fun_nlpot, &function);
		//opt.set_min_objective(fun_nlpot00005, &fun);

		opt.set_ftol_rel(1e-6);

		opt.set_xtol_rel(0.00005);

		// Set the initial point
		::std::vector<double> x(3);
		//x[0] = 0.48862; x[1] = 4.44033; x[2] = 1.00015;
		x[0] = 0.3; x[1] = 7; x[2] = 1;
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
			x[0] = 0.3; x[1] = 7; x[2] = 1;
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
			ss << ::std::setprecision(12) << index_m <<" _ Optimization succeeded, found minimum at f("
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


int add(int a, int b)
{
	return a + b;
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
	std::vector<std::vector<double>> subGrid = createGrid(std::vector<double>(xl.begin(),xl.end()-2), std::vector<double>(xu.begin(), xu.end() - 2), std::vector<double>(xd.begin(), xd.end() - 2));
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

	::std::ifstream inFile;
	inFile.open(szfile, ::std::ios::in);
	if (!inFile.is_open())
	{
		feLogError("Open file \"");
		feLogError(szfile);
		feLogError("\" failed.\n");
		return false;
	}
	feLog("successfully open file \"");
	feLog(szfile);
	feLog("\".\n");

	// Parse
		// could import a script language>?
		// normal text file> 
	::std::string unilang_filename = (::std::string(szfile) + ".unilang");

	//judge if the unilang file exists
	if (::std::filesystem::exists(unilang_filename))
	{
		int u_vfm_pm = 0;
		//double u_vfm_bpm = 120;
		unilang_regist_int_value("u_vfm_pm", &u_vfm_pm);
		unilang_regist_double_value("u_vfm_bpm", &(this->BMP));

		unilang_regist_bool_value("u_vfm_isRead_FEMresult_fromsavefile", &(this->isRead_FEMresult_fromsavefile));
		unilang_regist_bool_value("u_vfm_isReadfromsaveOptimfunc", &(this->isReadfromsaveOptimfunc));

		unilang_regist_int_value("u_vfm_displacementdataNumber", &(this->displacementdataNumber));
		unilang_regist_int_value("u_vfm_stressdataNumber", &(this->stressdataNumber));
		unilang_regist_int_value("u_vfm_pressurevalueNumber", &(this->pressurevalueNumber));
		
		std::uintptr_t temp_p_solution = (std::uintptr_t)&(this->solution);
		unilang_regist_uintptr_t_value("u_vfm_solution", (&temp_p_solution));
		std::uintptr_t temp_p_pressure_load = (std::uintptr_t)&(this->pressure_load);
		unilang_regist_uintptr_t_value("u_vfm_pressure_load", (&temp_p_pressure_load));
		std::uintptr_t temp_p_fixed = (std::uintptr_t)&(this->fixed);
		unilang_regist_uintptr_t_value("u_vfm_fixed", (&temp_p_fixed));
		std::uintptr_t temp_p_fixednode = (std::uintptr_t)&(this->fixednode);
		unilang_regist_uintptr_t_value("u_vfm_fixednode", (&temp_p_fixednode));

		std::uintptr_t temp_p_readfromsaveOptimfunc_index= (std::uintptr_t)&(this->readfromsaveOptimfunc_index);
		unilang_regist_uintptr_t_value("u_vfm_readfromsaveOptimfunc_index", (&temp_p_readfromsaveOptimfunc_index));

		feLog("temp_unilang_p_solution:");
		stringstream ssssss_solution;
		ssssss_solution << (std::uintptr_t)temp_p_solution << "\n" << (std::uintptr_t)temp_p_pressure_load << "\n" << (std::uintptr_t)temp_p_fixed << "\n" << (std::uintptr_t)temp_p_fixednode << "\n";

		feLog(ssssss_solution.str().c_str());

		//std::function<int(int,int)> f_add = [](int x, int y) {return x + y; };
		//void* p_f_add = f_add.target<int(*)(int, int)>();
		const char * parem_p_f_add[2] = { "int", "int" };
		
		// print p_f_add address to stringstream

		static void* p_f_add = reinterpret_cast<void*>(add);
		unilang_regist_function_var("u_vfm_add", &p_f_add, "int", 2, parem_p_f_add);
		
		const char* parem_p_f_vt_pushback[3] = { "uintptr_t", "string", "string" };
		static void* p_f_vector__tuple__string_string = reinterpret_cast<void*>(push_back_vector__tuple__string_string);
		unilang_regist_function_var("u_vfm_vt_pushback", &p_f_vector__tuple__string_string, "int", 3, parem_p_f_vt_pushback);

		const char* parem_p_f_vi_pushback[2] = { "uintptr_t", "int" };
		static void* p_f_vector__int = reinterpret_cast<void*>(push_back_vector__int);
		unilang_regist_function_var("u_vfm_vi_pushback", &p_f_vector__int, "int", 2, parem_p_f_vi_pushback);

		unilang_load_file(unilang_filename.data());

	}
	else
	{

		// Text file
		::std::string line1;
		::std::istringstream iss;

		// Get parameters
		::std::getline(inFile, line1);
		iss.str(line1);
		size_t parameters_number = 0;
		iss >> parameters_number;
		iss.clear();

		// set initialize object parameters
		::std::string line2;
		::std::getline(inFile, line2);
		iss.str(line2);
		::std::vector<double> objects(parameters_number);
		for (size_t i = 0; i < parameters_number; ++i)
		{
			iss >> objects[i];
		}
		iss.clear();

		// get BMP
		::std::string line3;
		::std::getline(inFile, line3);
		iss.str(line3);
		this->BMP = 0.0;
		iss >> BMP;
		iss.clear();

		// get is read from save file
		::std::string line4;
		::std::getline(inFile, line4);
		iss.str(line4);
		this->isRead_FEMresult_fromsavefile = false;
		::std::string line4_temp_readFEMresult;
		iss >> line4_temp_readFEMresult >> ::std::boolalpha >> this->isRead_FEMresult_fromsavefile;
		iss.clear();

		::std::string line4x;
		::std::getline(inFile, line4x);
		iss.str(line4x);
		this->isReadfromsaveOptimfunc = false;
		::std::string line4x_temp_readOptimfunc;
		iss >> line4x_temp_readOptimfunc >> ::std::boolalpha >> this->isReadfromsaveOptimfunc;
		iss.clear();

		::std::string line5;
		::std::getline(inFile, line5);
		iss.str(line5);
		this->displacementdataNumber = 0;
		::std::string line5_0, line5_1;
		iss >> line5_0 >> line5_1 >> this->displacementdataNumber;
		iss.clear();

		::std::string line6;
		::std::getline(inFile, line6);
		iss.str(line6);
		this->pressurevalueNumber = 0;
		::std::string line6_0, line6_1;
		iss >> line6_0 >> line6_1 >> this->pressurevalueNumber;
		iss.clear();

		::std::string line;
		while (::std::getline(inFile, line))
		{
			// if line first is solution
			if (line.find("solution") != ::std::string::npos)
			{
				// get solution
				::std::istringstream iss;
				iss.str(line);
				::std::string solution_0, solution_1, solution_2;
				iss >> solution_0 >> solution_1 >> solution_2;
				this->solution.push_back({ solution_1,solution_2 });
			}
			else if (line.find("pressure_load") != ::std::string::npos)
			{
				// get pressure
				::std::istringstream iss;
				iss.str(line);
				::std::string pressure_0, pressure_1, pressure_2;
				iss >> pressure_0 >> pressure_1 >> pressure_2;
				this->pressure_load.push_back({ pressure_1,pressure_2 });
			}
			else if (line.find("fixed") != ::std::string::npos)
			{
				// get fixed
				::std::istringstream iss;
				iss.str(line);
				::std::string fixed_0, fixed_1, fixed_2;
				iss >> fixed_0 >> fixed_1 >> fixed_2;

				if (fixed_2 == "list")
				{
					// read this line as number and number, until the end of line
					int nodeID;
					while (iss >> nodeID)
					{
						this->fixednode.push_back(nodeID - 1);
					}
				}
				else
				{
					this->fixed.push_back({ fixed_1,fixed_2 });
				}


			}
		}

	}

	// number of cycles considered at the end
	int numCycle = 2;

	// create a save file
	this->outSavefile = string(szfile) + "_BMP"+::std::to_string(static_cast<int>(this->BMP)) + "_save.txt";
	if (this->isRead_FEMresult_fromsavefile == false)
	{
		// create a save file
		::std::ofstream outFile;
		outFile.open(this->outSavefile, ::std::ios::out);
		if (!outFile.is_open())
		{
			feLogError("Open file \"");
			feLogError("save.txt");
			feLogError("\" failed.\n");
			return false;
		}

	}

	// create a log file
	this->outlogfile=string(szfile) + "_BMP" + ::std::to_string(static_cast<int>(this->BMP)) +"_log.txt";
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
	this->dumpfile=string(szfile) + "_BMP" + ::std::to_string(static_cast<int>(this->BMP)) +"_dump.febdump";

	// write log
	outFile << endl;
	outFile << "BMP:" << this->BMP << endl;
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
			if (lcname == "BPM" + ::std::to_string(static_cast<int>(this->BMP)))
			{
				// set the load controller
				lc = temp_lc;
				break;
			}
		}

		// find surface_load name="PressureLoad1" and set new lc by BMP
		for (int j = 0; j < this->pressure_load.size(); j++)
		{
			for (int i = 0; i < fem.ModelLoads(); i++)
			{
				FEModelLoad& load0 = *(fem.ModelLoad(i));
				::std::string loadclassname = load0.GetFactoryClass()->GetClassName();
				::std::string loadname = load0.GetName();

				if (::std::get<0>(this->pressure_load[j]) == loadclassname && ::std::get<1>(this->pressure_load[j]) == loadname)
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
	//if (static_cast<int>(this->BMP) < 120)
	//{
	//	step.m_ntime *= (120 / static_cast<int>(this->BMP));

	//	step.m_dt0 = step.m_dt0 * (120 / static_cast<int>(this->BMP)); // add time step
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
	if (this->isRead_FEMresult_fromsavefile == false)
	{
		bool femsolveresult = fem.Solve();
	}

	read_solved_information(&fem, CB_SOLVED, (void*)this); // no need callback
	
	return femsolveresult;
}
