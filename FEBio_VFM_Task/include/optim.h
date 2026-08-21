#pragma once
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "common_FEBio.h"
#include "FEBio_refunction.h"

struct OptimExecutionOptions {
    bool output_iteration_stress = false;
    bool output_iteration_virtual_work = false;
    bool parallel = true;
    int num_threads = 0;
};

struct FunOptimParams {
    const std::vector<double>& timeArray;
    const std::vector<int>& solution_elementsID;
    const std::vector<int>& solution_elementsDomainID;
    const std::vector<std::vector<std::vector<double>>>& trueJArray;
    const std::vector<std::vector<std::vector<mat3d>>>& truedeformationGradientArray;
    const std::vector<std::vector<std::vector<std::vector<mat3ds>>>>& virtualstrainArrayV;
    const std::vector<std::vector<double>>& externalVirtualWork;
    const std::vector<std::vector<double>>& volumeVirtualWork;
    const bool isLaplaceVFM;
    const std::vector<std::complex<double>>& LaplaceVFM_s;
    const ::std::vector<::std::vector< std::complex<double>>>& externalVirtualWork_laplace;
    const ::std::vector<::std::vector<::std::vector<::std::vector<::std::complex<double>>>>>& Sepk2_dotdot_vStrain;
    const ::std::vector<::std::vector< std::complex<double>>>& fs;
    const bool output_internalVirtualWork_to_file;
    const OptimExecutionOptions execution;
};


double fun_for_optim(FEModel* fem, double p_g, double p_t, double p_E, std::ofstream& outFile, const FunOptimParams& params);

double fun_for_optim_T(FEModel* fem, double p_g, double p_t, double p_E, std::ofstream& outFile, const::std::vector<double>& timeArray, const::std::vector<::std::vector<double>>& exter_nEvw, const::std::vector<::std::vector<double>>& internal_normal_visco, const::std::vector<int>& visco_mask, const::std::vector<::std::vector<::std::vector<mat3ds>>>& S_e_0, const std::vector<std::vector<std::vector<std::vector<mat3ds>>>>& virtualstrainArrayV, const std::vector<std::vector<std::vector<double>>>& trueJArray, const std::vector<std::vector<std::vector<mat3d>>>& truedeformationGradientArray, const OptimExecutionOptions& execution);

template<size_t N>
double fun_call_nlpot(unsigned n, const double* x, double* grad, void* data)
{
	::std::vector<double> x_array(x, x + N);

	return (*(static_cast<::std::function<double(const ::std::vector<double>&)>*>(data)))(x_array);
};

double fun_nlpot(unsigned n, const double* x, double* grad, void* data);

double fun_nlpot00005(unsigned n, const double* x, double* grad, void* data);

double fun_for_optim_E_gamma_tau(FEModel* fem, double E, double gamma, double tau, const FunOptimParams& params);

double fun_for_optim_E_gamma_tau_invF(FEModel* fem, double E, double gamma, double tau, const FunOptimParams& params);

::std::vector<::std::vector<double>> fun_optim_E__gamma_tau_Laplace(FEModel* fem, double gamma, double tau, const FunOptimParams& params);

::std::vector<::std::vector<double>> fun_optim_gamma__E_tau_Laplace(FEModel* fem, double E, double tau, const FunOptimParams& params);

::std::vector<::std::vector<double>> fun_optim_tau__E_gamma_Laplace(FEModel* fem, double E, double gamma, const FunOptimParams& params);

double fun_for_optim_elastic_E(FEModel* fem, double p_E, std::ofstream& outFile, const FunOptimParams& params);

#define CONST_SHIFT 50
