#include "VFMTask_laplace.h"
#include "optim.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <span>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>
void L_externalVirtualWork_LaplcaeTransform(std::vector<std::vector<double>>& externalVirtualWork, std::vector<std::vector<Eigen::dcomplex>>& externalVirtualWork_laplace, std::vector<double>& timeArray, VFMTask* task)
{
	int vf_u_size = externalVirtualWork[0].size();

	externalVirtualWork_laplace = ::std::vector<::std::vector<::std::complex<double>>>(vf_u_size);

	for (int index_vf = 0; index_vf < vf_u_size; index_vf++)
	{
		::std::vector<double> externalVirtualWork_VF;
		// Slice externalVirtualWork[:, index_vf].
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
		// Slice externalVirtualWork[:, index_vf].
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

// Input: timestep series ts (increasing, possibly with a zero tail) and heart rate bpm.
// Output: the half-open interval [start, end) for the penultimate complete cycle.
std::pair<int, int> find_last_two_cycles(const std::span<double>& ts, double bpm)
{
	if (ts.empty()) {
		throw std::runtime_error("timestep is empty");
	}

	// Find the final valid sample, ignoring zero or non-increasing tail values.
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

	// Find the start of one cycle ending at end_index.
	auto find_cycle_start = [&](int end_index) {
		double target = ts[end_index] - period;
		// Find the first sample at or above target in [0, end_index].
		auto first = ts.begin();
		auto it = std::lower_bound(first, first + (end_index + 1), target);
		int idx = static_cast<int>(it - first);
		if (idx > end_index) idx = end_index;
		return idx;
		};

	// Find the start of the final cycle.
	int last_start_index = find_cycle_start(last_valid);
	int last_end_exclusive = last_valid + 1; // Half-open interval including last_valid.

	// Find the start of the penultimate cycle.
	int second_last_end = last_start_index - 1;
	int second_last_start_index = (second_last_end >= 0)
		? find_cycle_start(second_last_end)
		: 0;
	int second_last_end_exclusive = last_start_index;

	// Return the half-open interval [start, end).
	int start_index = second_last_start_index;
	int end_index = second_last_end_exclusive;

	// Fall back to the final cycle when two complete cycles are unavailable.
	if (start_index >= end_index) {
		start_index = last_start_index;
		end_index = last_end_exclusive;
	}

	return { start_index, end_index };
}
