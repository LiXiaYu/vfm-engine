#include "common_FEBio.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath> 
#include <iomanip>


void write_to_log_2(FEModel* fem, const std::string& logs_string, std::ofstream& outFile)
{
	write_log(fem, 0, logs_string.c_str());
	outFile << logs_string;
}

// convert vec3d to ::std::vector<double>
::std::vector<double> vec3d_to_vector(const vec3d& v) { return ::std::vector<double>{v.x, v.y, v.z}; }


// 计算三角形的面积
double CalculateTriangleArea(const vec3d& A, const vec3d& B, const vec3d& C) {
	vec3d a = B - A;
	vec3d b = C - A;
	return 0.5 * (a ^ b).norm();
}

// 计算四边形的面积（假设它是平面且非自交的）
// 注意：这个实现假设四边形可以通过对角线划分为两个三角形
double CalculateQuadArea(const vec3d& A, const vec3d& B, const vec3d& C, const vec3d& D) {
	return CalculateTriangleArea(A, B, C) + CalculateTriangleArea(A, C, D);
}

// 计算拉普拉斯变换（周期函数）
std::vector<std::complex<double>> laplace_transform_periodic(const std::vector<double>& data, const std::vector<double>& t, const std::vector<std::complex<double>>& s) {
	std::vector<std::complex<double>> Ls(s.size(), std::complex<double>(0, 0));
	double total_time = t.back() - t.front();
	::std::vector<double> t_begin0(t);
	for (size_t i = 0; i < t.size(); ++i)
	{
		t_begin0[i] -= t.front();
	}
	for (size_t i = 0; i < s.size(); ++i)
	{
		std::vector<std::complex<double>> integrand(data.size());
		for (size_t j = 0; j < data.size(); ++j)
		{
			integrand[j] = data[j] * std::exp(-s[i] * t_begin0[j]);
		}
		
		Ls[i] = simpson_integration(integrand, t_begin0);
		Ls[i] *= 1.0 / (1.0 - std::exp(-s[i] * total_time));

	}
	return Ls;
}

std::complex<double> laplace_transform_periodic(const std::vector<double>& data, const std::vector<double>& t, const std::complex<double>& s)
{
	return laplace_transform_periodic(data, t, std::vector<std::complex<double>>(1, s))[0];
}

// eigne mat3ds with complex<double>
std::vector<Eigen::Matrix<std::complex<double>, 3, 3>> laplace_transform_periodic(const std::vector<mat3ds>& data, const std::vector<double>& t, const std::vector<std::complex<double>>& s)
{
	// convert data ot std::vector<Eigen::Matrix<double, 3, 3>>
	std::vector<Eigen::Matrix<double, 3, 3>> data_d(data.size(), Eigen::Matrix<double, 3, 3>::Zero());
	for (size_t i = 0; i < data.size(); ++i)
	{
		data_d[i] = convert_mat3ds_EigenMatrix(data[i]);
	}

	std::vector<Eigen::Matrix<std::complex<double>, 3, 3>> Ls(s.size(), Eigen::Matrix<std::complex<double>, 3, 3>::Zero());
	double total_time = t.back() - t.front();
	::std::vector<double> t_begin0(t);
	for (size_t i = 0; i < t.size(); ++i)
	{
		t_begin0[i] -= t.front();
	}

	for (size_t i = 0; i < s.size(); ++i )
	{
		std::vector<Eigen::Matrix<std::complex<double>, 3, 3>> integrand(data.size());
		for (size_t j = 0; j < data.size(); ++j)
		{
			integrand[j] = data_d[j] * std::exp(-s[i] * t_begin0[j]);
		}
		Ls[i] = simpson_integration(integrand, t_begin0);
		Ls[i] *= 1.0 / (1.0 - std::exp(-s[i] * total_time));
	}

	return Ls;
}

::std::vector<::std::complex<double>> laplace_transform_inT(const ::std::vector<double>& data, const std::vector<double>& t, const std::vector<std::complex<double>>& s)
{
	std::vector<std::complex<double>> Ls(s.size(), std::complex<double>(0, 0));
	double total_time = t.back() - t.front();
	::std::vector<double> t_begin0(t);
	for (size_t i = 0; i < t.size(); ++i)
	{
		t_begin0[i] -= t.front();
	}
	for (size_t i = 0; i < s.size(); ++i)
	{
		std::vector<std::complex<double>> integrand(data.size());
		for (size_t j = 0; j < data.size(); ++j)
		{
			integrand[j] = data[j] * std::exp(-s[i] * t_begin0[j]);
		}

		Ls[i] = simpson_integration(integrand, t_begin0);
	}
	return Ls;
}

std::vector<Eigen::Matrix<std::complex<double>, 3, 3>> laplace_transform_inT(const std::vector<mat3ds>& data, const std::vector<double>& t, const std::vector<std::complex<double>>& s)
{
	// convert data ot std::vector<Eigen::Matrix<double, 3, 3>>
	std::vector<Eigen::Matrix<double, 3, 3>> data_d(data.size(), Eigen::Matrix<double, 3, 3>::Zero());
	for (size_t i = 0; i < data.size(); ++i)
	{
		data_d[i] = convert_mat3ds_EigenMatrix(data[i]);
	}

	std::vector<Eigen::Matrix<std::complex<double>, 3, 3>> Ls(s.size(), Eigen::Matrix<std::complex<double>, 3, 3>::Zero());
	double total_time = t.back() - t.front();
	::std::vector<double> t_begin0(t);
	for (size_t i = 0; i < t.size(); ++i)
	{
		t_begin0[i] -= t.front();
	}

	for (size_t i = 0; i < s.size(); ++i)
	{
		std::vector<Eigen::Matrix<std::complex<double>, 3, 3>> integrand(data.size());
		for (size_t j = 0; j < data.size(); ++j)
		{
			integrand[j] = data_d[j] * std::exp(-s[i] * t_begin0[j]);
		}
		Ls[i] = simpson_integration(integrand, t_begin0);
	}

	return Ls;
}


// convert mat3ds to Eigen::Matrix<double, 3, 3>
Eigen::Matrix<double, 3, 3> convert_mat3ds_EigenMatrix(const mat3ds& data)
{
	// convert data ot std::vector<Eigen::Matrix<double, 3, 3>>
	Eigen::Matrix<double, 3, 3> data_d;
	for (int j = 0; j < 3; ++j)
	{
		for (int k = 0; k < 3; ++k)
		{
			data_d(j, k) = data(j, k);
		}
	}

	return data_d;
}
Eigen::Matrix<double, 3, 3> convert_mat3d_EigenMatrix(const mat3d& data)
{
	// convert data ot std::vector<Eigen::Matrix<double, 3, 3>>
	Eigen::Matrix<double, 3, 3> data_d;
	for (int j = 0; j < 3; ++j)
	{
		for (int k = 0; k < 3; ++k)
		{
			data_d(j, k) = data(j, k);
		}
	}

	return data_d;
}

// beta_s need in beta + j -Inf → beta + j Inf
::std::vector<double> inverse_laplace_transform(const std::vector<std::complex<double>>& beta_data, const std::vector<double>& t, const std::vector<std::complex<double>>& beta_s) {
	std::vector<double> y(t.size());
	double total_time = t.back() - t.front();
	::std::vector<double> t_begin0(t);
	for (size_t i = 0; i < t.size(); ++i)
	{
		t_begin0[i] -= t.front();
	}

	for (size_t i = 0; i < t_begin0.size(); ++i)
	{
		::std::vector<std::complex<double>> integrand(beta_data.size());
		for (size_t j = 0; j < beta_data.size(); ++j)
		{
			integrand[j] = beta_data[j] * std::exp(beta_s[j] * t_begin0[i]);
		}

		::std::complex<double> yi_c = simpson_integration(integrand, beta_s);
		auto yi_c____ = (yi_c / ::std::complex<double>(0, 2 * PI));
		y[i] = yi_c____.real();
	}

	return y;
}

// 创建Talbot的积分路径
// 对每个t都有一个beta_s
::std::vector<::std::vector<std::complex<double>>> createTalbotPath(const ::std::vector<double>& t, const int N = 32, const double shift = 0)
{
	double total_time = t.back() - t.front();
	::std::vector<double> t_begin0(t);
	for (size_t i = 0; i < t.size(); ++i)
	{
		t_begin0[i] -= t.front();
	}

	::std::vector<::std::vector<std::complex<double>>> beta_s(t.size(), ::std::vector<std::complex<double>>(N, std::complex<double>(0, 0)));

	double h = 2 * PI / N;

	for (size_t index_time = 0; index_time < t.size(); ++index_time)
	{
		for (int k = 0; k < N; ++k)
		{
			double theta = -PI + (k + 0.5) * h;
			std::complex<double> z = shift + static_cast<double>(N) / t_begin0[index_time] * (0.5017 * theta / std::tan(0.6407 * theta) + std::complex<double>(-0.6122, 0.2645 * theta));
			std::complex<double> dz = static_cast<double>(N) / t_begin0[index_time] * (-0.5017 * 0.6407 * theta * (std::pow(1.0 / std::sin(0.6407 * theta), 2.0)) + 0.5017 / std::tan(0.6407 * theta) + std::complex<double>(0.0, 0.2645));

			beta_s[index_time][k] = z;
		}
	}

	return beta_s;
}

// Talbot方法实现
// s is Talbot path, F is F(s)
// N is s.size(), default is 32
// don't need s, only need F(s)
::std::vector<double> talbotInverseLaplaceTransform(const ::std::vector<::std::vector<::std::complex<double>>>& F, const ::std::vector<double>& t, const int N = 32, const double shift = 0)
{
	double total_time = t.back() - t.front();
	::std::vector<double> t_begin0(t);	
	for (size_t i = 0; i < t.size(); ++i)
	{
		t_begin0[i] -= t.front();
	}

	::std::vector<double> f(t.size(), 0);

	for (size_t index_time = 0; index_time < t.size(); ++index_time)
	{
		std::complex<double> ans(0.0, 0.0);

		double h = 2 * PI / N;

		for (int k = 0; k < N; ++k)
		{
			double theta = -PI + (k + 0.5) * h;
			std::complex<double> z = shift + static_cast<double>(N) / t_begin0[index_time] * (0.5017 * theta / std::tan(0.6407 * theta) + std::complex<double>(-0.6122, 0.2645 * theta));
			std::complex<double> dz = static_cast<double>(N) / t_begin0[index_time] * (-0.5017 * 0.6407 * theta * (std::pow(1.0 / std::sin(0.6407 * theta), 2.0)) + 0.5017 / std::tan(0.6407 * theta) + std::complex<double>(0.0, 0.2645));
			ans = ans + std::exp(z * t_begin0[index_time]) * F[index_time][k] * dz;
		}

		f[index_time] = ((h / (std::complex<double>(0.0, 2.0 * PI))) * ans).real();
	}


	return f;
}

void write_vector2D_to_csv(const ::std::vector<::std::vector<double>>& data, const std::string& filename)
{
	std::ofstream outfile(filename, ::std::ios::ate | ::std::ios::out);
	if (outfile.is_open())
	{
		for (size_t i = 0; i < data.size(); ++i)
		{
			for (size_t j = 0; j < data[i].size(); ++j)
			{
				outfile << ::std::setprecision(12) << data[i][j];
				if (j < data[i].size() - 1)
				{
					outfile << ",";
				}
			}
			if (i < data.size() - 1)
			{
				outfile << "\n";
			}
		}
	}
	else
	{
		std::cout << "Unable to open file" << std::endl;
		exit(1);
	}
	outfile.close();
}


mat3ds PK2stress(FEIsotropicElastic& material, const FEMaterialPoint& mp)
{
	const FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	const mat3d& F = pt.m_F;
	double Ji = 1.0 / pt.m_J;

	double E = material.m_E(mp);
	double v = material.m_v(mp);

	// lame parameters
	double lam = Ji * (v * E / ((1 + v) * (1 - 2 * v)));
	double mu = Ji * (0.5 * E / (1 + v));

	mat3ds I(1,1,1,0,0,0);

	// calculate Euler-lagrange strain tensor (ie. b-matrix)
	mat3ds e = 0.5 * (pt.RightCauchyGreen() - I);

	mat3ds pk2 = lam * e.tr() * I + 2.0 * mu * e;

	return pk2;
}


void exit_febio()
{
	// 结束程序
	exit(0);
}