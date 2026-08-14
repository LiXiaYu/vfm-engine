#include "optim.h"

double fun_for_optim(FEModel* fem, double p_g, double p_t, double p_E, std::ofstream& outFile, const FunOptimParams& params)
{
	// set new material parameters
	FEMesh& mesh = fem->GetMesh();

	int vf_u_size = params.externalVirtualWork[0].size();

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

	double loss = 0;

	// Laplace the external virtual work
	if (params.isLaplaceVFM == true)
	{
		//define internalVirtualWork_laplace
		std::vector<std::vector<Eigen::dcomplex>> internalVirtualWork_laplace(vf_u_size, std::vector<Eigen::dcomplex>(params.LaplaceVFM_s.size()));

		for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
		{
			::std::vector<::std::vector<::std::vector<Eigen::dcomplex>>> stresspk2_withJc_VF_alle_laplace(params.solution_elementsID.size());
			for (int index_seId = 0; index_seId < params.solution_elementsID.size(); index_seId++)
			{
				int j = params.solution_elementsID[index_seId];

				FEElement& element = *(mesh.Element(j));
				int nodes_number = element.Nodes();

				// find material by element
				FEMaterial* pmat = fem->GetMaterial(element.GetMatID());


				// convert element to solid element
				FESolidElement& solidElement = static_cast<FESolidElement&>(element);

				stresspk2_withJc_VF_alle_laplace[index_seId].resize(solidElement.GaussPoints());
				for (int index_gauss = 0; index_gauss < solidElement.GaussPoints(); index_gauss++)
				{
					// get element's stress
					FEMaterialPoint& mp = *solidElement.GetMaterialPoint(index_gauss);

					FEElasticMaterialPoint& pt_e = *mp.ExtractData<FEElasticMaterialPoint>();

					FEViscoElasticMaterialPoint* p_pt = mp.ExtractData<FEViscoElasticMaterialPoint>();
					
					for(int index_s = 0; index_s < params.LaplaceVFM_s.size(); index_s++)
					{
						auto Sepk2_dotdot_vStrain = params.Sepk2_dotdot_vStrain[index_vf][index_seId][index_gauss][index_s];
						if(p_pt != nullptr)
						{
							::std::complex<double> visco_pe = material->m_g0.constValue();
							for (int i = 0; i < 6; i++)
							{
								visco_pe += material->m_g[i] * params.LaplaceVFM_s[index_s] / (params.LaplaceVFM_s[index_s] + 1.0 / material->m_t[i]);
							}
							
							FEViscoElasticMaterial* material = static_cast<FEViscoElasticMaterial*>(pmat);
							FECoreBase* elasticb = material->FindProperty("elastic")->get(0);

							// FENeoHookean
							FENeoHookean& elasticm = static_cast<FENeoHookean&>(*elasticb);
							FEParamDouble* material_E = elasticm.GetParameter("E")->pvalue<FEParamDouble>();
							double linerElactical_e = material_E->constValue();

							visco_pe*= linerElactical_e * visco_pe;
							stresspk2_withJc_VF_alle_laplace[index_seId][index_gauss].push_back(visco_pe * Sepk2_dotdot_vStrain);
						}
						else
						{
							stresspk2_withJc_VF_alle_laplace[index_seId][index_gauss].push_back(Sepk2_dotdot_vStrain);
						}

					}
				}
			}
		
			for (int index_s = 0; index_s < params.LaplaceVFM_s.size(); index_s++)
			{
				for (int index_seId = 0; index_seId < params.solution_elementsID.size(); index_seId++)
				{
					for (int index_gauss = 0; index_gauss < stresspk2_withJc_VF_alle_laplace[index_seId].size(); index_gauss++)
					{
						internalVirtualWork_laplace[index_vf][index_s] += stresspk2_withJc_VF_alle_laplace[index_seId][index_gauss][index_s];
					}
				}
			}
		}

		

		// claculate loss
		::std::complex<double> loss_d = 0;
		for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
		{
			::std::complex<double> vf_loos = 0;
			for (int index_s = 0; index_s < params.LaplaceVFM_s.size(); index_s++)
			{
				vf_loos += internalVirtualWork_laplace[index_vf][index_s] - params.externalVirtualWork_laplace[index_vf][index_s];
			}
			loss_d += ::std::abs(vf_loos);
		}
		//loss = loss_d.imag();
		loss = ::std::sqrt(::std::pow(loss_d.real(), 2) + ::std::pow(loss_d.imag(), 2));
	}
	else
	{
		// calculate true displacement -> strain -> stress

		const int iter_max = 5000;//50;

		::std::vector<::std::vector<mat3ds>> last_psts(params.trueJArray[0].size());
		for (int index_iter = 0; index_iter < iter_max; index_iter++)
		{
			double error_iter = 0;

			::std::vector<::std::vector<double>> error_iters(params.solution_elementsID.size());

			for (int index_timestep = 0; index_timestep < params.timeArray.size(); index_timestep++)
			{
				if (index_iter == 0 && index_timestep == params.timeArray.size() - 1)
				{
					last_psts = ::std::vector<::std::vector<mat3ds>>(params.trueJArray[index_timestep].size());
				}


#pragma omp parallel for
				for (int index_seId = 0; index_seId < params.solution_elementsID.size(); index_seId++)
				{
					int j = params.solution_elementsID[index_seId];

					FEElement& element = *(mesh.Element(j));
					int nodes_number = element.Nodes();

					// find material by element
					FEMaterial* pmat = fem->GetMaterial(element.GetMatID());


					// convert element to solid element
					FESolidElement& solidElement = static_cast<FESolidElement&>(element);

					// get domain
					FESolidDomain& domain = static_cast<FESolidDomain&>(mesh.Domain(params.solution_elementsDomainID[index_seId]));

					if (index_iter == 0 && index_timestep == params.timeArray.size() - 1)
					{
						last_psts[index_seId] = ::std::vector<mat3ds>(solidElement.GaussPoints());
					}
					// initilize error_iters as 0
					if (index_timestep == params.timeArray.size() - 1)
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

						pt_e.m_F = params.truedeformationGradientArray[index_timestep][index_seId][n];
						pt_e.m_J = params.truedeformationGradientArray[index_timestep][index_seId][n].det();// trueJArray[index_timestep][index_seId][n];

						fem->GetTime().timeIncrement = index_timestep == 0 ? (params.timeArray[index_timestep + 1] - params.timeArray[index_timestep]) : (params.timeArray[index_timestep] - params.timeArray[index_timestep - 1]);
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

							if (index_timestep == params.timeArray.size() - 1)
							{
								if (index_iter > 0)
								{
									auto error_pts = pts - last_psts[index_seId][n];
									double error = 0.0;
									for (int i = 0; i < 3; i++)
									{
										for (int j = 0; j < 3; j++)
										{
											error += ::std::abs(error_pts(i, j));
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
						if (index_timestep == params.timeArray.size() - 1)
						{
							last_psts[index_seId][n] = pts;
						}
					}
				}
			}

			// sum error_iters
			for (int index_seId = 0; index_seId < params.solution_elementsID.size(); index_seId++)
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
		::std::vector<::std::vector<double>> internalVirtualWork(params.timeArray.size());
		for (int index_timestep = 0; index_timestep < params.timeArray.size(); index_timestep++)
		{
			internalVirtualWork[index_timestep] = ::std::vector<double>(vf_u_size); // internalVirtualWork
			for (int i = 0; i < vf_u_size; i++)
			{
				internalVirtualWork[index_timestep][i] = 0;
			}
		}

		::std::vector<::std::vector<mat3ds>> iter_stress(params.timeArray.size(), ::std::vector<mat3ds>(params.solution_elementsID.size(), mat3ds(0, 0, 0, 0, 0, 0)));

		::std::vector<::std::vector<mat3ds>> init_stress(params.solution_elementsID.size());

		for (int index_timestep = 0; index_timestep < params.timeArray.size(); index_timestep++)
		{
			// internal virtual work
			//logs_string = "calculate internal virtual work.\n";
			//write_log(fem, 0, logs_string.c_str());
			//outFile << logs_string;

			::std::vector<::std::vector<double>> ivw_vf_element(params.solution_elementsID.size());
			for (int index_seId = 0; index_seId < params.solution_elementsID.size(); index_seId++)
			{
				ivw_vf_element[index_seId] = ::std::vector<double>(vf_u_size);
				for (int i = 0; i < vf_u_size; i++)
				{
					ivw_vf_element[index_seId][i] = 0;
				}
			}

#pragma omp parallel for
			for (int index_seId = 0; index_seId < params.solution_elementsID.size(); index_seId++)
			{
				int j = params.solution_elementsID[index_seId];

				FEElement& element = *(mesh.Element(j));
				int nodes_number = element.Nodes();

				init_stress[index_seId] = ::std::vector<mat3ds>(element.GaussPoints());

				// find material by element
				FEMaterial* pmat = fem->GetMaterial(element.GetMatID());

				// convert element to solid element
				FESolidElement& solidElement = static_cast<FESolidElement&>(element);

				// get domain
				FESolidDomain& domain = static_cast<FESolidDomain&>(mesh.Domain(params.solution_elementsDomainID[index_seId]));

				for (int n = 0; n < solidElement.GaussPoints(); n++)
				{
					// get element's stress
					FEMaterialPoint& mp = *solidElement.GetMaterialPoint(n);

					FEElasticMaterialPoint& pt_e = *mp.ExtractData<FEElasticMaterialPoint>();

					FEViscoElasticMaterialPoint* p_pt = mp.ExtractData<FEViscoElasticMaterialPoint>();

					FEViscoElasticMaterialPoint& pt = *p_pt;

					pt_e.m_F = params.truedeformationGradientArray[index_timestep][index_seId][n];
					pt_e.m_J = params.truedeformationGradientArray[index_timestep][index_seId][n].det();
					double Jc = params.trueJArray[index_timestep][index_seId][n];

					fem->GetTime().timeIncrement = index_timestep == 0 ? (params.timeArray[index_timestep + 1] - params.timeArray[index_timestep]) : (params.timeArray[index_timestep] - params.timeArray[index_timestep - 1]);

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

					if (j == 35072 && n == 1)
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

					mat3ds pts_i_PK2 = pt_e.pull_back(pts_i);

					for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
					{
						auto& vStrain = params.virtualstrainArrayV[index_timestep][index_vf][index_seId][n];
						double temp_dd = pts_i.dotdot(vStrain) * Jc;
						double temp_td = (pts_i * vStrain).trace() * Jc;

						double temp_dd_PK2 = pts_i_PK2.dotdot(vStrain) * Jc;

						if (j == 35072 && n == 0)
						{
							int awer23dsfcsafwerf = 0;
						}

						if (::std::isnan(temp_dd))
						{
							int awer23dsfcsafwerf = 0;
						}

						// correct with stress
						double weight_correct = 1;// 1.0 / (10000 * ::std::abs((pts_i.xx() + pts_i.yy() + pts_i.zz()) / 3) + 0.01);

						ivw_vf_element[index_seId][index_vf] += weight_correct * temp_dd_PK2;
					}

				}

				iter_stress[index_timestep][index_seId] /= solidElement.GaussPoints();
			}

			for (int index_seId = 0; index_seId < params.solution_elementsID.size(); index_seId++)
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

		//double loss = 0;
		//for (int index_timestep = 0; index_timestep < timeArray.size(); index_timestep++)
		//{
		//	for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
		//	{
		//		loss += ::std::abs(internalVirtualWork[index_timestep][index_vf] - externalVirtualWork[index_timestep][index_vf] /*+ volumeVirtualWork[index_timestep][index_vf]*/); // externalVirtualWork opposite of matlab code
		//		//loss += ::std::abs((internalVirtualWork[index_timestep][index_vf] - externalVirtualWork[index_timestep][index_vf]) / externalVirtualWork[index_timestep][index_vf]);
		//	}
		//}

		for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
		{
			double vf_loss = 0;
			double vf_evw = 0;
			//double vf_evw = 0;
			for (int index_timestep = 0; index_timestep < params.timeArray.size(); index_timestep++)
			{
				vf_loss += ::std::abs(internalVirtualWork[index_timestep][index_vf] - params.externalVirtualWork[index_timestep][index_vf]);
				//vf_evw += externalVirtualWork[index_timestep][index_vf];
				// vf_loss 除以params.externalVirtualWork[index_timestep][index_vf];如果vf_evw为0，则loss为无穷大
				vf_evw += params.externalVirtualWork[index_timestep][index_vf];
			}
			//loss += ::std::abs(vf_loss/vf_evw);
			loss += ::std::abs(vf_loss);
		}

#ifdef DEBUG
		if (params.output_internalVirtualWork_to_file)
		{
			// output internalVirtualWork[i][j] to csv file
			::std::ofstream internalVirtualWork_file;
			internalVirtualWork_file.open("./temp/debug/result/internalVirtualWork_" + ::std::to_string(p_E) + "_" + ::std::to_string(p_g) + "_" + ::std::to_string(p_t) + ".csv", ::std::ios::ate | ::std::ios::out);
			for (int index_timestep = 0; index_timestep < params.timeArray.size(); index_timestep++)
			{
				for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
				{
					internalVirtualWork_file << ::std::setprecision(12) << internalVirtualWork[index_timestep][index_vf] << ",";
				}
				if (index_timestep != params.timeArray.size() - 1)
				{
					internalVirtualWork_file << "\n";
				}
			}
			internalVirtualWork_file.close();
			::std::ofstream externalVirtualWork_file;
			externalVirtualWork_file.open("./temp/debug/result/externalVirtualWork.csv", ::std::ios::ate | ::std::ios::out);
			for (int index_timestep = 0; index_timestep < params.timeArray.size(); index_timestep++)
			{
				for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
				{
					externalVirtualWork_file << ::std::setprecision(12) << params.externalVirtualWork[index_timestep][index_vf] << ",";
				}
				if (index_timestep != params.timeArray.size() - 1)
				{
					externalVirtualWork_file << "\n";
				}
			}
			externalVirtualWork_file.close();

			// output iter_stress[i][j] to csv file
			::std::ofstream iter_stress_file;
			iter_stress_file.open("./temp/debug/result/iter_stress_" + ::std::to_string(p_E) + "_" + ::std::to_string(p_g) + "_" + ::std::to_string(p_t) + ".csv", ::std::ios::ate | ::std::ios::out);
			for (int index_seId = 0; index_seId < params.solution_elementsID.size(); index_seId++)
			{
				for (int index_timestep = 0; index_timestep < params.timeArray.size(); index_timestep++)
				{
					iter_stress_file << ::std::setprecision(12) << iter_stress[index_timestep][index_seId].xx() << "," << iter_stress[index_timestep][index_seId].yy() << "," << iter_stress[index_timestep][index_seId].zz() << "," << iter_stress[index_timestep][index_seId].xy() << "," << iter_stress[index_timestep][index_seId].yz() << "," << iter_stress[index_timestep][index_seId].xz();
					if (index_timestep != params.timeArray.size() - 1)
					{
						iter_stress_file << ",";
					}
				}
				if (index_seId != params.solution_elementsID.size() - 1)
				{
					iter_stress_file << "\n";
				}
			}
			iter_stress_file.close();

			// output iter_stersss[i][j] to csv file
			::std::ofstream iter_stress_elementIDselected_file;
			iter_stress_elementIDselected_file.open("./temp/debug/result/surface_elementIDselected/iter_stress_" + ::std::to_string(p_E) + "_" + ::std::to_string(p_g) + "_" + ::std::to_string(p_t) + ".csv", ::std::ios::ate | ::std::ios::out);
			for (int index_seId = 0; index_seId < params.solution_elementsID.size(); index_seId++)
			{
				if (params.solution_elementsID[index_seId] != 35072)
				{
					continue;
				}

				for (int index_timestep = 0; index_timestep < params.timeArray.size(); index_timestep++)
				{
					iter_stress_elementIDselected_file << ::std::setprecision(12) << iter_stress[index_timestep][index_seId].xx() << "," << iter_stress[index_timestep][index_seId].yy() << "," << iter_stress[index_timestep][index_seId].zz() << "," << iter_stress[index_timestep][index_seId].xy() << "," << iter_stress[index_timestep][index_seId].yz() << "," << iter_stress[index_timestep][index_seId].xz() << ",";
					if (index_timestep != params.timeArray.size() - 1)
					{
						iter_stress_elementIDselected_file << "\n";
					}
				}
				//if (index_seId != solution_elementsID.size() - 1)
				//{
				//	iter_stress_elementIDselected_file << "\n";
				//}
			}
			iter_stress_elementIDselected_file.close();
		}
#endif // DEBUG

	}


	return loss;

}

double fun_for_optim_T(FEModel* fem, double p_g, double p_t, double p_E, std::ofstream& outFile, const ::std::vector<double>& timeArray, const ::std::vector<::std::vector<double>>& exter_nEvw, const ::std::vector<::std::vector<double>>& internal_normal_visco, const ::std::vector<int>& visco_mask, const ::std::vector<::std::vector<::std::vector<mat3ds>>>& S_e_0, const std::vector<std::vector<std::vector<std::vector<mat3ds>>>>& virtualstrainArrayV, const std::vector<std::vector<std::vector<double>>>& trueJArray, const std::vector<std::vector<std::vector<mat3d>>>& truedeformationGradientArray)
{
	int vf_u_size = exter_nEvw[0].size();

	// t = timeArray - timeArray[0]
	auto t(timeArray);
	for (int index_timestep = 0; index_timestep < timeArray.size(); index_timestep++)
	{
		t[index_timestep] = t[index_timestep] - timeArray[0];
	}
	double T = t[t.size() - 1] - t[0];


	::std::vector<::std::vector<double>> internalVirtualWork_visco(t.size(), ::std::vector<double>(vf_u_size, 0));

	::std::vector<::std::vector<::std::vector<mat3ds>>> S(S_e_0);
	::std::vector<::std::vector<::std::vector<mat3ds>>> stress(S_e_0);

	#pragma omp parallel for
	for (int index_timestep = 0; index_timestep < t.size(); index_timestep++)
	{
		for (int index_seId = 0; index_seId < S[index_timestep].size(); index_seId++)
		{
			if(visco_mask[index_seId] == 1)
			{
				for (int index_n = 0; index_n < S[index_timestep][index_seId].size(); index_n++)
				{
					::std::vector<mat3ds> ForIntergate(t.size(), mat3ds(0));
					for(int index_t = 0; index_t < t.size(); index_t++)
					{
						int index_I_t = index_timestep + index_t;
						if (index_timestep + index_t >= t.size())
						{
							index_I_t = index_timestep + index_t - t.size() + 1;
						}
						ForIntergate[index_t] = ::std::exp(t[index_t] / p_t) * (S_e_0[index_I_t][index_seId][index_n]-S_e_0[index_t][index_seId][index_n]);
					}
					
					auto Intergate_eSe0 = simpson_integration(ForIntergate, t);


					S[index_timestep][index_seId][index_n] = p_E * S_e_0[index_timestep][index_seId][index_n] + p_E * p_g * (S_e_0[index_timestep][index_seId][index_n] - 1/ p_t * (::std::exp(-T/ p_t) / (1 - ::std::exp(-T / p_t))) * Intergate_eSe0);

					//auto s_xx = S[index_timestep][index_seId][index_n](0, 0);
					//::std::ostringstream ise0_str;
					//ise0_str << "S:" << ::std::setprecision(12) << s_xx << ::std::endl;
					//write_to_log_2(fem, ise0_str.str(), outFile);

					double Jc = trueJArray[index_timestep][index_seId][index_n];

					for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
					{
						auto& vStrain = virtualstrainArrayV[index_timestep][index_vf][index_seId][index_n];

						double inter_we=S[index_timestep][index_seId][index_n].dotdot(vStrain) * Jc;

						internalVirtualWork_visco[index_timestep][index_vf] += inter_we;
					}

				}
			}
		
			for (int index_n = 0; index_n < S[index_timestep][index_seId].size(); index_n++)
			{
				auto& F = truedeformationGradientArray[index_timestep][index_seId][index_n];
				auto detF = F.det();

				// Cauchy stress
				mat3d stress_cauchy = ((F * S[index_timestep][index_seId][index_n]) * F.transpose()) * (1.0 / detF);
				// save as symmetric matrix
				mat3ds stress_cauchy_sym(stress_cauchy(0, 0), stress_cauchy(1, 1), stress_cauchy(2, 2), stress_cauchy(0, 1), stress_cauchy(1, 2), stress_cauchy(0, 2));
				stress[index_timestep][index_seId][index_n] = stress_cauchy_sym;
			}
		}
	}

	double loss = 0;
	#pragma omp parallel for
	for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
	{
		for (int index_timestep = 0; index_timestep < t.size(); index_timestep++)
		{
			double loss_vw = ::std::abs(exter_nEvw[index_timestep][index_vf] - internalVirtualWork_visco[index_timestep][index_vf]);
			loss += loss_vw;
		}
	}

	// output iter_stress[i][j] to csv file
	::std::vector<::std::vector<double>> temp_iter_stress(S[0].size(), ::std::vector<double>(S.size() * 6, 0.0));
	#pragma omp parallel for
	for (int index_seId = 0; index_seId < S[0].size(); index_seId++)
	{
		for (int index_timestep = 0; index_timestep < S.size(); index_timestep++)
		{
			auto PK2 = S[index_timestep][index_seId][0];
			for (int index_n = 1; index_n < S[index_timestep][index_seId].size(); index_n++)
			{
				PK2+=S[index_timestep][index_seId][index_n];
			}
			PK2 /= S[index_timestep][index_seId].size();
			temp_iter_stress[index_seId][index_timestep * 6] = PK2(0, 0);
			temp_iter_stress[index_seId][index_timestep * 6 + 1] = PK2(1, 1);
			temp_iter_stress[index_seId][index_timestep * 6 + 2] = PK2(2, 2);
			temp_iter_stress[index_seId][index_timestep * 6 + 3] = PK2(0, 1);
			temp_iter_stress[index_seId][index_timestep * 6 + 4] = PK2(1, 2);
			temp_iter_stress[index_seId][index_timestep * 6 + 5] = PK2(0, 2);
		}
	}
	write_vector2D_to_csv(temp_iter_stress, "./temp/debug/result/T_PK2stress_" + ::std::to_string(p_E) + "_" + ::std::to_string(p_g) + "_" + ::std::to_string(p_t) + ".csv");

	// output iter_stress[i][j] to csv file

	::std::vector<::std::vector<double>> temp_T_stress(stress[0].size(), ::std::vector<double>(stress.size() * 6, 0.0));
	#pragma omp parallel for
	for (int index_seId = 0; index_seId < stress[0].size(); index_seId++)
	{
		for (int index_timestep = 0; index_timestep < stress.size(); index_timestep++)
		{
			auto s = stress[index_timestep][index_seId][0];
			for (int index_n = 1; index_n < stress[index_timestep][index_seId].size(); index_n++)
			{
				s += stress[index_timestep][index_seId][index_n];
			}
			s /= stress[index_timestep][index_seId].size();
			temp_T_stress[index_seId][index_timestep * 6] = s(0, 0);
			temp_T_stress[index_seId][index_timestep * 6 + 1] = s(1, 1);
			temp_T_stress[index_seId][index_timestep * 6 + 2] = s(2, 2);
			temp_T_stress[index_seId][index_timestep * 6 + 3] = s(0, 1);
			temp_T_stress[index_seId][index_timestep * 6 + 4] = s(1, 2);
			temp_T_stress[index_seId][index_timestep * 6 + 5] = s(0, 2);
		}
	}
	write_vector2D_to_csv(temp_T_stress, "./temp/debug/result/T_stress_" + ::std::to_string(p_E) + "_" + ::std::to_string(p_g) + "_" + ::std::to_string(p_t) + ".csv");

	// output internalVirtualWork[i][j] to csv file
	write_vector2D_to_csv(internalVirtualWork_visco, "./temp/debug/result/T_internalVirtualWork_visco_" + ::std::to_string(p_E) + "_" + ::std::to_string(p_g) + "_" + ::std::to_string(p_t) + ".csv");


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

double fun_for_optim_E_gamma_tau(FEModel* fem, double E, double gamma, double tau, const FunOptimParams& params)
{
	// set new material parameters
	FEMesh& mesh = fem->GetMesh();

	int vf_u_size = params.externalVirtualWork[0].size();

	::std::vector<::std::vector<::std::complex<double>>> loss_c(vf_u_size, ::std::vector<::std::complex<double>>(params.LaplaceVFM_s.size(), ::std::complex<double>(0, 0)));

	double loss = 0;

	for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
	{
		for (int index_s = 0; index_s < params.LaplaceVFM_s.size(); index_s++)
		{
			auto loss_c_complex = (1.0 + gamma * params.LaplaceVFM_s[index_s] / (params.LaplaceVFM_s[index_s] + 1.0 / tau)) * E - params.fs[index_vf][index_s];
			loss_c[index_vf][index_s] = loss_c_complex;
			if (!std::isnan(loss_c[index_vf][index_s].real()) && !std::isinf(loss_c[index_vf][index_s].real()))
			{
				loss += ::std::abs(loss_c[index_vf][index_s].real());
			}
			else
			{
				auto temp_debug_jwjwjdkji32 = 12;
			}

		}
	}

	return loss;
}

// isTalbotLaplaceVFM_s
double fun_for_optim_E_gamma_tau_invF(FEModel* fem, double E, double gamma, double tau, const FunOptimParams& params)
{
	// set new material parameters
	FEMesh& mesh = fem->GetMesh();

	int vf_u_size = params.externalVirtualWork[0].size();

	::std::vector<::std::vector<double>> loss(vf_u_size, ::std::vector<double>(params.timeArray.size(), 0));
	double loss_sum = 0;

	int N = 32;
	double shift = CONST_SHIFT;

	for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
	{
		// convert params.fs[index_vf] to Fs[index_time][index_s] each time has N s
		::std::vector<::std::vector<::std::complex<double>>> Fs(params.timeArray.size(), ::std::vector<::std::complex<double>>(N, ::std::complex<double>(0, 0)));
		for (int index_time = 0; index_time < params.timeArray.size(); index_time++)
		{
			for (int index_s = 0; index_s < N; index_s++)
			{
				Fs[index_time][index_s] = params.fs[index_vf][N * index_time + index_s];
			}
		}

		auto invF = talbotInverseLaplaceTransform(Fs, params.timeArray, N, shift);

		for(int index_time = 2; index_time < params.timeArray.size(); index_time++)
		{
			double time_r = params.timeArray[index_time] - params.timeArray[0];
			loss[index_vf][index_time] = - gamma / tau * E * ::std::exp(-time_r / tau) + invF[index_time];
			loss_sum += loss[index_vf][index_time];
		}
	}

	return loss_sum;
}

// return E [index_vf][index_s]
::std::vector<::std::vector<double>> fun_optim_E__gamma_tau_Laplace(FEModel* fem, double gamma, double tau, const FunOptimParams& params)
{
	// set new material parameters
	FEMesh& mesh = fem->GetMesh();

	int vf_u_size = params.externalVirtualWork[0].size();

	::std::vector<::std::vector<::std::complex<double>>> Es_c(vf_u_size, ::std::vector<::std::complex<double>>(params.LaplaceVFM_s.size(), ::std::complex<double>(0, 0)));
	// claculate E
	::std::vector<::std::vector<double>> Es(vf_u_size, ::std::vector<double>(params.LaplaceVFM_s.size(), 0));
	for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
	{
		::std::complex<double> vf_f = 0;
		for (int index_s = 0; index_s < params.LaplaceVFM_s.size(); index_s++)
		{
			auto Es_complex = params.fs[index_vf][index_s] / (1.0 + gamma * params.LaplaceVFM_s[index_s] / (params.LaplaceVFM_s[index_s] + 1.0 / tau ));
			//Es[index_vf][index_s] = ::std::sqrt(::std::pow(Es_complex.real(), 2) + ::std::pow(Es_complex.imag(), 2));
			Es[index_vf][index_s] = Es_complex.real();
			Es_c[index_vf][index_s] = Es_complex;
		}
	}

	//// sum all Es_c for each vf
	//::std::vector<::std::complex<double>> sum_Es_c(vf_u_size, ::std::complex<double>(0, 0));
	//for (int index_s = 0; index_s < params.LaplaceVFM_s.size(); index_s++)
	//{
	//	for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
	//	{
	//		sum_Es_c[index_vf] += Es_c[index_vf][index_s];
	//	}
	//}

	//for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
	//{
	//	::std::vector<::std::complex<double>> Es_c_beta;
	//	::std::vector<::std::complex<double>> s_beta;
	//	double beta = 0.0;
	//	for (int index_s = 0; index_s < params.LaplaceVFM_s.size(); index_s++)
	//	{
	//		if (::std::abs(params.LaplaceVFM_s[index_s].real() - beta) < 1e-5)
	//		{
	//			Es_c_beta.push_back(Es_c[index_vf][index_s]);
	//			s_beta.push_back(params.LaplaceVFM_s[index_s]);
	//		}
	//	}

	//	auto inv_Et = inverse_laplace_transform(Es_c_beta, params.timeArray, s_beta);

	//	auto temp_wait_debug_23iopwejkdik = 1221;
	//}

	return Es;
}

::std::vector<::std::vector<double>> fun_optim_gamma__E_tau_Laplace(FEModel* fem, double E, double tau, const FunOptimParams& params)
{
	// set new material parameters
	FEMesh& mesh = fem->GetMesh();

	int vf_u_size = params.externalVirtualWork[0].size();

	::std::vector<::std::vector<::std::complex<double>>> gammas_c(vf_u_size, ::std::vector<::std::complex<double>>(params.LaplaceVFM_s.size(), ::std::complex<double>(0, 0)));
	// claculate E
	::std::vector<::std::vector<double>> gammas(vf_u_size, ::std::vector<double>(params.LaplaceVFM_s.size(), 0));
	for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
	{
		::std::complex<double> vf_f = 0;
		for (int index_s = 0; index_s < params.LaplaceVFM_s.size(); index_s++)
		{
			auto gammas_complex = (params.fs[index_vf][index_s] / E - 1.0) * (1.0 + 1.0 / tau / params.LaplaceVFM_s[index_s]);
			//gammas[index_vf][index_s] = ::std::sqrt(::std::pow(gammas_complex.real(),2)+::std::pow(gammas_complex.imag(),2));
			gammas[index_vf][index_s] = gammas_complex.real();
			gammas_c[index_vf][index_s] = gammas_complex;
		}
	}

	//// sum all gammas_c for each vf
	//::std::vector<::std::complex<double>> sum_gammas_c(vf_u_size, ::std::complex<double>(0, 0));
	//for (int index_s = 0; index_s < params.LaplaceVFM_s.size(); index_s++)
	//{
	//	for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
	//	{
	//		sum_gammas_c[index_vf] += gammas[index_vf][index_s];
	//	}
	//}

	//for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
	//{
	//	::std::vector<::std::complex<double>> gammas_c_beta;
	//	::std::vector<::std::complex<double>> s_beta;
	//	double beta = 0.0;
	//	for (int index_s = 0; index_s < params.LaplaceVFM_s.size(); index_s++)
	//	{
	//		if (::std::abs(params.LaplaceVFM_s[index_s].real() - beta) < 1e-5)
	//		{
	//			gammas_c_beta.push_back(gammas_c[index_vf][index_s]);
	//			s_beta.push_back(params.LaplaceVFM_s[index_s]);
	//		}
	//	}

	//	auto inv_gammat = inverse_laplace_transform(gammas_c_beta, params.timeArray, s_beta);

	//	auto temp_wait_debug_23iopwejkdik = 1221;
	//}

	return gammas;
}

::std::vector<::std::vector<double>> fun_optim_tau__E_gamma_Laplace(FEModel* fem, double E, double gamma, const FunOptimParams& params)
{
	// set new material parameters
	FEMesh& mesh = fem->GetMesh();

	int vf_u_size = params.externalVirtualWork[0].size();

	::std::vector<::std::vector<::std::complex<double>>> taus_c_withoutfracs(vf_u_size, ::std::vector<::std::complex<double>>(params.LaplaceVFM_s.size(), ::std::complex<double>(0, 0)));

	::std::vector<::std::vector<::std::complex<double>>> taus_c(vf_u_size, ::std::vector<::std::complex<double>>(params.LaplaceVFM_s.size(), ::std::complex<double>(0, 0)));

	// claculate E
	::std::vector<::std::vector<double>> taus(vf_u_size, ::std::vector<double>(params.LaplaceVFM_s.size(), 0));
	for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
	{
		::std::complex<double> vf_f = 0;
		for (int index_s = 0; index_s < params.LaplaceVFM_s.size(); index_s++)
		{
			auto taus_complex = 1.0 / (gamma / (params.fs[index_vf][index_s] / E - 1.0) - 1.0) / params.LaplaceVFM_s[index_s];
			
			auto taus_withoutfracs = 1.0 / (gamma / (params.fs[index_vf][index_s] / E - 1.0) - 1.0);

			taus_c_withoutfracs[index_vf][index_s] = taus_withoutfracs;

			taus[index_vf][index_s] = taus_complex.real();
			taus_c[index_vf][index_s] = taus_complex;
		}
	}


	//// sum all taus_c for each vf
	//::std::vector<::std::complex<double>> sum_taus_c(vf_u_size, ::std::complex<double>(0, 0));
	//for (int index_s = 0; index_s < params.LaplaceVFM_s.size(); index_s++)
	//{
	//	for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
	//	{
	//		sum_taus_c[index_vf] += taus_c[index_vf][index_s];
	//	}
	//}

	for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
	{
		::std::vector<::std::complex<double>> taus_c_beta;
		::std::vector<::std::complex<double>> s_beta;
		double beta = 0.0;
		for (int index_s = 0; index_s < params.LaplaceVFM_s.size(); index_s++)
		{
			if(::std::abs(params.LaplaceVFM_s[index_s].real() - beta) < 1e-5)
			{
				taus_c_beta.push_back(taus_c[index_vf][index_s]);
				s_beta.push_back(params.LaplaceVFM_s[index_s]);
			}
		}
                                                                                        
		auto inv_taut = inverse_laplace_transform(taus_c_beta, params.timeArray, s_beta);

		auto temp_wait_debug_23iopwejkdik = 1221;
	}


	return taus;
}

double fun_for_optim_elastic_E(FEModel* fem, double p_E, std::ofstream& outFile, const FunOptimParams& params)
{
	FEMesh& mesh = fem->GetMesh();

	int vf_u_size = params.externalVirtualWork[0].size();

	// t = timeArray - timeArray[0]
	auto t(params.timeArray);
	for (int index_timestep = 0; index_timestep < params.timeArray.size(); index_timestep++)
	{
		t[index_timestep] = t[index_timestep] - params.timeArray[0];
	}
	double T = t[t.size() - 1] - t[0]; // 无用

	FENeoHookean* material = nullptr;

	for (int i = 0; i < fem->Materials(); i++)
	{
		::std::string materialname = fem->GetMaterial(i)->GetName();
		::std::string materialclassname = fem->GetMaterial(i)->GetFactoryClass()->GetClassName();

		if (materialclassname == "FENeoHookean")
		{
			material = static_cast<FENeoHookean*>(fem->GetMaterial(i));

			material->m_E = p_E;
		}
	}

	double loss = 0;

	// internal virtual work
	::std::vector<::std::vector<double>> internalVirtualWork(params.timeArray.size());
	for (int index_timestep = 0; index_timestep < params.timeArray.size(); index_timestep++)
	{
		internalVirtualWork[index_timestep] = ::std::vector<double>(vf_u_size); // internalVirtualWork
		for (int i = 0; i < vf_u_size; i++)
		{
			internalVirtualWork[index_timestep][i] = 0;
		}
	}

	::std::vector<::std::vector<::std::vector<mat3ds>>> S(params.timeArray.size());
	for (int index_timestep = 0; index_timestep < params.timeArray.size(); index_timestep++)
	{
		S[index_timestep] = ::std::vector<::std::vector<mat3ds>>(params.solution_elementsID.size());
		for (int index_seId = 0; index_seId < params.solution_elementsID.size(); index_seId++)
		{
			int j = params.solution_elementsID[index_seId];
			FEElement& element = *(mesh.Element(j));
			FESolidElement& solidElement = static_cast<FESolidElement&>(element);

			S[index_timestep][index_seId] = ::std::vector<mat3ds>(solidElement.GaussPoints());
		}
	}

	for (int index_timestep = 0; index_timestep < params.timeArray.size(); index_timestep++)
	{
		::std::vector<::std::vector<double>> ivw_vf_element(params.solution_elementsID.size());
		for (int index_seId = 0; index_seId < params.solution_elementsID.size(); index_seId++)
		{
			ivw_vf_element[index_seId] = ::std::vector<double>(vf_u_size);
			for (int i = 0; i < vf_u_size; i++)
			{
				ivw_vf_element[index_seId][i] = 0;
			}
		}

//#pragma omp parallel for
		for (int index_seId = 0; index_seId < params.solution_elementsID.size(); index_seId++)
		{
			int j = params.solution_elementsID[index_seId];

			FEElement& element = *(mesh.Element(j));
			int nodes_number = element.Nodes();

			// find material by element
			FEMaterial* pmat = fem->GetMaterial(element.GetMatID());


			// convert element to solid element
			FESolidElement& solidElement = static_cast<FESolidElement&>(element);

			// get domain
			FESolidDomain& domain = static_cast<FESolidDomain&>(mesh.Domain(params.solution_elementsDomainID[index_seId]));

			for (int n = 0; n < solidElement.GaussPoints(); n++)
			{
				// get element's stress
				FEMaterialPoint& mp = *solidElement.GetMaterialPoint(n);

				FEElasticMaterialPoint& pt_e = *mp.ExtractData<FEElasticMaterialPoint>();

				pt_e.m_F = params.truedeformationGradientArray[index_timestep][index_seId][n];
				pt_e.m_J = params.truedeformationGradientArray[index_timestep][index_seId][n].det();
				double Jc = params.trueJArray[index_timestep][index_seId][n];

				mat3ds pts;

				pts = material->Stress(mp);

				double material_E = material->m_E.constValue();
				//write to log 2 material E
				//write_to_log_2(fem, "p_E: " + ::std::to_string(p_E) + " material_E: " + ::std::to_string(material_E), outFile);

				mat3ds pts_PK2 = pt_e.pull_back(pts);

				S[index_timestep][index_seId][n] = pts_PK2;

				for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
				{
					auto& vStrain = params.virtualstrainArrayV[index_timestep][index_vf][index_seId][n];

					double temp_dd_PK2 = pts_PK2.dotdot(vStrain) * Jc;


					ivw_vf_element[index_seId][index_vf] += temp_dd_PK2;
				}
			}
		}

		for (int index_seId = 0; index_seId < params.solution_elementsID.size(); index_seId++)
		{
			for (int i = 0; i < vf_u_size; i++)
			{
				internalVirtualWork[index_timestep][i] += ivw_vf_element[index_seId][i];
			}
		}
	}

	for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
	{
		double vf_loss = 0;
		double vf_evw = 0;
		//double vf_evw = 0;
		for (int index_timestep = 0; index_timestep < params.timeArray.size(); index_timestep++)
		{
			vf_loss += ::std::abs(params.externalVirtualWork[index_timestep][index_vf] - internalVirtualWork[index_timestep][index_vf]);
			//vf_evw += externalVirtualWork[index_timestep][index_vf];
			// vf_loss 除以params.externalVirtualWork[index_timestep][index_vf];如果vf_evw为0，则loss为无穷大
			vf_evw += params.externalVirtualWork[index_timestep][index_vf];
		}
		//loss += ::std::abs(vf_loss/vf_evw);
		loss += ::std::abs(vf_loss);
	}


	// output iter_stress[i][j] to csv file
	::std::vector<::std::vector<double>> temp_iter_stress(S[0].size(), ::std::vector<double>(S.size() * 6, 0.0));
#pragma omp parallel for
	for (int index_seId = 0; index_seId < S[0].size(); index_seId++)
	{
		for (int index_timestep = 0; index_timestep < S.size(); index_timestep++)
		{
			auto PK2 = S[index_timestep][index_seId][0];
			for (int index_n = 1; index_n < S[index_timestep][index_seId].size(); index_n++)
			{
				PK2 += S[index_timestep][index_seId][index_n];
			}
			PK2 /= S[index_timestep][index_seId].size();
			temp_iter_stress[index_seId][index_timestep * 6] = PK2(0, 0);
			temp_iter_stress[index_seId][index_timestep * 6 + 1] = PK2(1, 1);
			temp_iter_stress[index_seId][index_timestep * 6 + 2] = PK2(2, 2);
			temp_iter_stress[index_seId][index_timestep * 6 + 3] = PK2(0, 1);
			temp_iter_stress[index_seId][index_timestep * 6 + 4] = PK2(1, 2);
			temp_iter_stress[index_seId][index_timestep * 6 + 5] = PK2(0, 2);
		}
	}
	write_vector2D_to_csv(temp_iter_stress, "./temp/debug/result/T_PK2stress_" + ::std::to_string(p_E) + ".csv");

	// output internalVirtualWork[i][j] to csv file
	write_vector2D_to_csv(internalVirtualWork, "./temp/debug/result/internalVirtualWork_" + ::std::to_string(p_E) + ".csv");


	return loss;
}
