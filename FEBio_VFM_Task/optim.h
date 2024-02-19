#pragma once
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "common_FEBio.h"
#include "FEBio_refunction.h"

double fun_for_optim(FEModel* fem, double p_g, double p_t, double p_E, const std::vector<double>& timeArray, std::ofstream& outFile, FEMesh& mesh, const::std::vector<int>& solution_elementsID, const::std::vector<int>& solution_elementsDomainID, const::std::vector<::std::vector<::std::vector<double>>>& trueJArray, const::std::vector<::std::vector<::std::vector<mat3d>>>& truedeformationGradientArray, const::std::vector<::std::vector<std::vector<std::vector<mat3ds>>>>& virtualstrainArrayV, const::std::vector<::std::vector<double>>& externalVirtualWork, const::std::vector<::std::vector<double>>& volumeVirtualWork);


template<size_t N>
double fun_call_nlpot(unsigned n, const double* x, double* grad, void* data)
{
	::std::vector<double> x_array(x, x + N);

	return (*(static_cast<::std::function<double(const ::std::vector<double>&)>*>(data)))(x_array);
};

double fun_nlpot(unsigned n, const double* x, double* grad, void* data);

double fun_nlpot00005(unsigned n, const double* x, double* grad, void* data);
