#include "optim.h"

double fun_for_optim(FEModel* fem, double p_g, double p_t, double p_E, const std::vector<double>& timeArray, std::ofstream& outFile, FEMesh& mesh, const ::std::vector<int>& solution_elementsID, const ::std::vector<int>& solution_elementsDomainID, const ::std::vector<::std::vector<::std::vector<double>>>& trueJArray, const ::std::vector<::std::vector<::std::vector<mat3d>>>& truedeformationGradientArray, const ::std::vector<::std::vector<std::vector<std::vector<mat3ds>>>>& virtualstrainArrayV, const ::std::vector<::std::vector<double>>& externalVirtualWork, const ::std::vector<::std::vector<double>>& volumeVirtualWork)
{
	// set new material parameters

	int vf_u_size = externalVirtualWork[0].size();

	FEViscoElasticMaterial* material = nullptr;

	//FEIsotropicElastic* pre_postLC = nullptr;

	for (int i = 0; i < fem->Materials(); i++)
	{
		::std::string materialname = fem->GetMaterial(i)->GetName();
		::std::string materialclassname = fem->GetMaterial(i)->GetFactoryClass()->GetClassName();
		if (materialclassname == "FEViscoElasticMaterial")
		{
			material = static_cast<FEViscoElasticMaterial*>(fem->GetMaterial(i));
			//material.m_g[0] = p_g;
			//material.m_t[0] = p_t;
			material->SetParameter("g1", p_g);
			material->SetParameter("t1", p_t);


			FECoreBase* elasticb = material->FindProperty("elastic")->get(0);
			::std::string elastic_classname = elasticb->GetFactoryClass()->GetClassName();
			if (elastic_classname == "FENeoHookean")
			{
				FENeoHookean& elasticm = static_cast<FENeoHookean&>(*elasticb);
				FEParamDouble* material_E = elasticm.GetParameter("E")->pvalue<FEParamDouble>();
				(*material_E) = p_E;
				//elasticm.SetParameter("E", p_E);
			}

		}
	}

	// calculate true displacement -> strain -> stress

	const int iter_max = 50;

	::std::vector<::std::vector<mat3ds>> last_psts(trueJArray[0].size());
	for (int index_iter = 0; index_iter < iter_max; index_iter++)
	{
		double error_iter = 0;

		::std::vector<::std::vector<double>> error_iters(solution_elementsID.size());

		for (int index_timestep = 0; index_timestep < timeArray.size(); index_timestep++)
		{
			if (index_iter == 0 && index_timestep == timeArray.size() - 1)
			{
				last_psts = ::std::vector<::std::vector<mat3ds>>(trueJArray[index_timestep].size());
			}


#pragma omp parallel for
			for (int index_seId = 0; index_seId < solution_elementsID.size(); index_seId++)
			{
				int j = solution_elementsID[index_seId];

				FEElement& element = *(mesh.Element(j));
				int nodes_number = element.Nodes();

				// find material by element
				FEMaterial* pmat = fem->GetMaterial(element.GetMatID());


				// convert element to solid element
				FESolidElement& solidElement = static_cast<FESolidElement&>(element);

				// get domain
				FESolidDomain& domain = static_cast<FESolidDomain&>(mesh.Domain(solution_elementsDomainID[index_seId]));

				if (index_iter == 0 && index_timestep == timeArray.size() - 1)
				{
					last_psts[index_seId] = ::std::vector<mat3ds>(solidElement.GaussPoints());
				}
				// initilize error_iters as 0
				if (index_timestep == timeArray.size() - 1)
				{
					error_iters[index_seId] = ::std::vector<double>(solidElement.GaussPoints());
					for (int i = 0; i < solidElement.GaussPoints(); i++)
					{
						error_iters[index_seId][i] = 0.0;
					}
				}


				for (int n = 0; n < solidElement.GaussPoints(); n++)
				{
					// get element's stress
					FEMaterialPoint& mp = *solidElement.GetMaterialPoint(n);

					FEElasticMaterialPoint& pt_e = *mp.ExtractData<FEElasticMaterialPoint>();

					FEViscoElasticMaterialPoint* p_pt = mp.ExtractData<FEViscoElasticMaterialPoint>();

					FEViscoElasticMaterialPoint& pt = *p_pt;

					pt_e.m_F = truedeformationGradientArray[index_timestep][index_seId][n];
					pt_e.m_J = truedeformationGradientArray[index_timestep][index_seId][n].det();// trueJArray[index_timestep][index_seId][n];

					fem->GetTime().timeIncrement = index_timestep == 0 ? (timeArray[index_timestep + 1] - timeArray[index_timestep]) : (timeArray[index_timestep] - timeArray[index_timestep - 1]);
					mat3ds pts;
					if (p_pt != nullptr)
					{
						if (index_iter == 0 && index_timestep == 0)
						{
							pt.m_Sep.zero();
							for (int i = 0; i < FEViscoElasticMaterialPoint::MAX_TERMS; i++)
							{
								pt.m_Hp[i].zero();
								pt.m_alphap[i] = 0;//don't need?
							}
							pt.m_sedp = 0;//don't need?

						}
						else
						{
							pt.m_Sep = pt.m_Se;
							for (int i = 0; i < FEViscoElasticMaterialPoint::MAX_TERMS; i++)
							{
								pt.m_Hp[i] = pt.m_H[i];
								pt.m_alphap[i] = pt.m_alpha[i];//don't need?
							}
							pt.m_sedp = pt.m_sed;//don't need?
						}


						pts = material->Stress(mp);

						if (index_timestep == timeArray.size() - 1)
						{
							if (index_iter > 0)
							{
								auto error_pts = pts - last_psts[index_seId][n];
								double error = 0.0;
								for (int i = 0; i < 3; i++)
								{
									for (int j = 0; j < 3; j++)
									{
										error += ::std::abs(error_pts(i,j));
									}

								}
								error_iters[index_seId][n] = error;
							}

						}

					}
					else
					{
						pts = dynamic_cast<FESolidMaterial*>(pmat)->Stress(mp);
					}

					//mat3ds stress = pt.m_s; // error stress for parameters
					if (index_timestep == timeArray.size() - 1)
					{
						last_psts[index_seId][n] = pts;
					}
				}
			}
		}

		// sum error_iters
		for (int index_seId = 0; index_seId < solution_elementsID.size(); index_seId++)
		{
			for (int i = 0; i < error_iters[index_seId].size(); i++)
			{
				error_iter += error_iters[index_seId][i];
			}
		}

		if (index_iter > 0)
		{
			if (error_iter < 1e-15)
			{
				break;
			}
		}
		if (index_iter == iter_max - 1)
		{
			write_to_log_2(fem, "initial stress too long...\n", outFile);
			// return 9999999;
		}
	}

	// internal virtual work
	::std::vector<::std::vector<double>> internalVirtualWork(timeArray.size());
	for (int index_timestep = 0; index_timestep < timeArray.size(); index_timestep++)
	{
		internalVirtualWork[index_timestep] = ::std::vector<double>(vf_u_size); // internalVirtualWork
		for (int i = 0; i < vf_u_size; i++)
		{
			internalVirtualWork[index_timestep][i] = 0;
		}
	}

	::std::vector<::std::vector<mat3ds>> iter_stress(timeArray.size(), ::std::vector<mat3ds>(solution_elementsID.size(), mat3ds(0, 0, 0, 0, 0, 0)));

	::std::vector<::std::vector<mat3ds>> init_stress(solution_elementsID.size());

	for (int index_timestep = 0; index_timestep < timeArray.size(); index_timestep++)
	{
		// internal virtual work
		//logs_string = "calculate internal virtual work.\n";
		//write_log(fem, 0, logs_string.c_str());
		//outFile << logs_string;

		::std::vector<::std::vector<double>> ivw_vf_element(solution_elementsID.size());
		for (int index_seId = 0; index_seId < solution_elementsID.size(); index_seId++)
		{
			ivw_vf_element[index_seId] = ::std::vector<double>(vf_u_size);
			for (int i = 0; i < vf_u_size; i++)
			{
				ivw_vf_element[index_seId][i] = 0;
			}
		}

#pragma omp parallel for
		for (int index_seId = 0; index_seId < solution_elementsID.size(); index_seId++)
		{
			int j = solution_elementsID[index_seId];

			FEElement& element = *(mesh.Element(j));
			int nodes_number = element.Nodes();

			init_stress[index_seId] = ::std::vector<mat3ds>(element.GaussPoints());

			// find material by element
			FEMaterial* pmat = fem->GetMaterial(element.GetMatID());

			// convert element to solid element
			FESolidElement& solidElement = static_cast<FESolidElement&>(element);

			// get domain
			FESolidDomain& domain = static_cast<FESolidDomain&>(mesh.Domain(solution_elementsDomainID[index_seId]));

			for (int n = 0; n < solidElement.GaussPoints(); n++)
			{
				// get element's stress
				FEMaterialPoint& mp = *solidElement.GetMaterialPoint(n);

				FEElasticMaterialPoint& pt_e = *mp.ExtractData<FEElasticMaterialPoint>();

				FEViscoElasticMaterialPoint* p_pt = mp.ExtractData<FEViscoElasticMaterialPoint>();

				FEViscoElasticMaterialPoint& pt = *p_pt;

				pt_e.m_F = truedeformationGradientArray[index_timestep][index_seId][n];
				pt_e.m_J = truedeformationGradientArray[index_timestep][index_seId][n].det();
				double Jc = trueJArray[index_timestep][index_seId][n];

				fem->GetTime().timeIncrement = index_timestep == 0 ? (timeArray[index_timestep + 1] - timeArray[index_timestep]) : (timeArray[index_timestep] - timeArray[index_timestep - 1]);

				mat3ds s_pts;

				mat3ds pts;
				if (p_pt != nullptr)
				{
					pt.m_Sep = pt.m_Se;
					for (int i = 0; i < FEViscoElasticMaterialPoint::MAX_TERMS; i++)
					{
						pt.m_Hp[i] = pt.m_H[i];
						pt.m_alphap[i] = pt.m_alpha[i];//don't need?
					}
					pt.m_sedp = pt.m_sed;//don't need?


					pts = material->Stress(mp);
				}
				else
				{
					::std::string pmname = pmat->GetFactoryClass()->GetClassName();
					pts = dynamic_cast<FESolidMaterial*>(pmat)->Stress(mp);

					s_pts = FEIsotropicElastic_Stress(dynamic_cast<FEIsotropicElastic&>(*pmat), mp);
				}

				mat3ds pts_i = pts - last_psts[index_seId][n]; // stress increment ?

				auto temp_debug_m_F = pt_e.m_F;
				auto temp_debug_m_J = pt_e.m_J;
				auto temp_debug_s_pts = s_pts;

				if (j == 35072 && n==1)
				{
					int awer23dsfcsafwerf = 0;
				}

				auto temp_debug_mpJ0(mp.m_J0);
				auto temp_debug_pts_i(pts_i);

				iter_stress[index_timestep][index_seId].xx() += pts.xx();
				iter_stress[index_timestep][index_seId].yy() += pts.yy();
				iter_stress[index_timestep][index_seId].zz() += pts.zz();
				iter_stress[index_timestep][index_seId].xy() += pts.xy();
				iter_stress[index_timestep][index_seId].yz() += pts.yz();
				iter_stress[index_timestep][index_seId].xz() += pts.xz();
				//mat3ds stress = pt.m_s; // error stress for parameters



				for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
				{
					auto vStrain = virtualstrainArrayV[index_timestep][index_vf][index_seId][n];
					double temp_dd = pts_i.dotdot(vStrain) * Jc;
					double temp_td =(pts_i* vStrain).trace() * Jc;

					if (j == 35072 && n == 0)
					{
						int awer23dsfcsafwerf = 0;
					}

					// correct with stress
					double weight_correct = 1;// 1.0 / (10000 * ::std::abs((pts_i.xx() + pts_i.yy() + pts_i.zz()) / 3) + 0.01);

					ivw_vf_element[index_seId][index_vf] = weight_correct * temp_dd;
				}

			}

			iter_stress[index_timestep][index_seId]/= solidElement.GaussPoints();
		}

		for (int index_seId = 0; index_seId < solution_elementsID.size(); index_seId++)
		{
			for (int i = 0; i < vf_u_size; i++)
			{
				internalVirtualWork[index_timestep][i] += ivw_vf_element[index_seId][i];
			}
		}
	}


	//DumpFile dumpfile(*fem);
	//dumpfile.Create("./stress_BPM60_solution_elementsID.dumpfile");
	//dumpfile& iter_stress;
	//dumpfile.Close();

	//::std::ofstream stress_file;
	//stress_file.open("./stress_BPM60_solution_elementsID.txt", ::std::ios::out);

	//for (int index_seId = 0; index_seId < solution_elementsID.size(); index_seId++)
	//{
	//	int j = solution_elementsID[index_seId];
	//	j++;
	//	stress_file << j;
	//	for (int index_timestep = 0; index_timestep < timeArray.size()-1; index_timestep++)
	//	{
	//		stress_file << "," << ::std::setprecision(12) << iter_stress[index_timestep][index_seId].xx() << "," << iter_stress[index_timestep][index_seId].yy() << "," << iter_stress[index_timestep][index_seId].zz() << "," << iter_stress[index_timestep][index_seId].xy() << "," << iter_stress[index_timestep][index_seId].yz() << "," << iter_stress[index_timestep][index_seId].xz();
	//	}
	//	if (index_seId != solution_elementsID.size() - 1)
	//	{
	//		stress_file << "\n";
	//	}
	//}
	//stress_file.close();

	double loss = 0;
	for (int index_timestep = 0; index_timestep < timeArray.size(); index_timestep++)
	{
		for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
		{
			loss += ::std::abs(internalVirtualWork[index_timestep][index_vf] - externalVirtualWork[index_timestep][index_vf] /*+ volumeVirtualWork[index_timestep][index_vf]*/); // externalVirtualWork opposite of matlab code
		}
	}

	return loss;

}

double fun_nlpot(unsigned n, const double* x, double* grad, void* data)
{
	return (*(static_cast<::std::function<double(const ::std::vector<double>&)>*>(data)))({ x[0],x[1],x[2] });
};

double fun_nlpot00005(unsigned n, const double* x, double* grad, void* data)
{
	// foreach x, 取小数点后五位
	double x00005[3];
	for (int i = 0; i < 3; i++)
	{
		x00005[i] = ::std::floor(x[i] * 100000 + 0.5) / 100000;
	}

	return (*(static_cast<::std::function<double(const ::std::vector<double>&)>*>(data)))({ x[0],x[1],x[2] });
};