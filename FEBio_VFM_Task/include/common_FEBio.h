#pragma once

#include <FECore/FECoreTask.h>
#include <FECore/log.h>
#include <FECore/Callback.h>
#include <FECore/FEInitialCondition.h>
#include <FECore/FEModelLoad.h>
#include <FECore/FESurface.h>
#include <FECore/FEFacetSet.h>
#include <FECore/DataStore.h>
#include <FECore/FEPlotData.h>
#include <FECore/FEPlotDataStore.h>
#include <FECore/DumpFile.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FETimeStepController.h>
#include <FECore/FEModel.h>
#include <FECore/FENLConstraint.h>
#include <FEBioMech/FEElasticMaterialPoint.h>
#include <FECore/FESolidDomain.h>
#include <FEBioMech/FEElasticMaterial.h>
#include <FEBioMech/FEIsotropicElastic.h>
#include <FEBioMech/FEPressureLoad.h>
#include <FEBioMech/FENodalForce.h>
#include <FECore/FELoadController.h>
#include <FECore/FELoadCurve.h>
#include <FEBioMech/FEViscoElasticMaterial.h>
#include <FEBioMech/FENeoHookean.h>
#include <FEBioMech/FEElasticSolidDomain.h>
#include <FEBioMech/FEVolumeConstraint.h>
#include <FEBioMech/FEIsotropicElastic.h>


#ifndef PI
#define PI 3.141592653589793
#endif

#include <vector>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <span>
#include <complex>
#include <type_traits>

#include <eigen3/Eigen/Dense>
#include <numeric>

#include <mio/mmap.hpp>
#ifdef GetClassName
	#undef GetClassName
#endif
#ifdef GetCurrentTime
	#undef GetCurrentTime
#endif


#ifdef CTEST_ENABLE
	#ifdef IMPORTING_DLL 
		#define DLL_PUBLIC __declspec(dllimport) 
	#else
		#define DLL_PUBLIC __declspec(dllexport) 
	#endif
#else
	#define DLL_PUBLIC
#endif

void write_to_log_2(FEModel* fem, const std::string& logs_string, std::ofstream& outFile);

::std::vector<double> vec3d_to_vector(const vec3d& v);

double CalculateTriangleArea(const vec3d& A, const vec3d& B, const vec3d& C);

double CalculateQuadArea(const vec3d& A, const vec3d& B, const vec3d& C, const vec3d& D);

//// define a concept which type must have +-*/=
//template <typename T>
//concept Arithmetic = requires(T a, T b) {
//	{a + b};
//	{a - b};
//	{a * b};
//	{a / b};
//};

template<typename T>
concept HasZeros = requires(T t)
{
	t.Zero();
};

template <typename T>
std::vector<T> diff(const std::vector<T>& x) {
	std::vector<T> dx(x.size() - 1);
	for (std::size_t i = 1; i < x.size(); ++i) {
		dx[i - 1] = x[i] - x[i - 1];
	}
	return dx;
}

template <typename T, typename Tx>
T _basic_simpson(const std::vector<T>& y, std::size_t start, std::size_t stop,
	const std::vector<Tx>* x, Tx dx, int axis) {
	std::size_t N = y.size();
	std::size_t step = 2;
	T result;
	if constexpr (::std::is_arithmetic_v<T>)
	{
		result = T(0);
	}
	else if constexpr (HasZeros<T>)
	{
		// 如果T.Zero()存在，则使用该函数 
		result = T::Zero();
	}
	else
	{
		result = T{ 0 };
	}
	Tx Tx_ZERO;
	if constexpr (::std::is_arithmetic_v<Tx>)
	{
		Tx_ZERO = Tx(0);
	}
	else if constexpr (HasZeros<Tx>)
	{
		// 如果Tx.Zero()存在，则使用该函数 
		Tx_ZERO = Tx::Zero();
	}
	else
	{
		Tx_ZERO = Tx{ 0 };
	}

	if (!x) {
		for (std::size_t i = start; i < stop; i += step) {
			result += y[i] + 4.0 * y[i + 1] + y[i + 2];
		}
		result *= dx / 3.0;
	}
	else {
		std::vector<Tx> h = diff(*x);
		std::vector<Tx> h0(h.size() / 2);
		std::vector<Tx> h1(h.size() / 2);
		for (std::size_t i = 0; i < h0.size(); ++i) {
			h0[i] = h[start + i * step];
			h1[i] = h[start + i * step + 1];
		}

		std::vector<Tx> hsum(h0.size()), hprod(h0.size()), h0divh1(h0.size());
		for (std::size_t i = 0; i < h0.size(); ++i) {
			hsum[i] = h0[i] + h1[i];
			hprod[i] = h0[i] * h1[i];
			h0divh1[i] = h1[i] != Tx_ZERO ? h0[i] / h1[i] : Tx_ZERO;
		}

		for (std::size_t i = start; i < stop; i += step) {
			result += (hsum[i / 2] / 6.0) *
				(y[i] * (2.0 - (h0divh1[i / 2] != Tx_ZERO ? 1.0 / h0divh1[i / 2] : Tx_ZERO)) +
					y[i + 1] * (hsum[i / 2] * hsum[i / 2] / hprod[i / 2]) +
					y[i + 2] * (2.0 - h0divh1[i / 2]));
		}
	}

	return result;
}

template <typename T, typename Tx>
T simpson(const std::vector<T>& y, const std::vector<Tx>* x = nullptr,
	Tx dx = 1.0, int axis = -1) {
	std::size_t N = y.size();
	if (N < 2) {
		throw std::invalid_argument("There must be at least two data points.");
	}
	T result, val;
	if constexpr (::std::is_arithmetic_v<T>)
	{
		result = T(0);
		val = T(0);
	}
	else if constexpr (HasZeros<T>)
	{
		// 如果T.Zero()存在，则使用该函数 
		result = T::Zero();
		val = T::Zero();
	}
	else
	{
		result = T{ 0 };
		val = T{ 0 };
	}

	if (N % 2 == 0) {
		if (N == 2) {
			val = 0.5 * dx * (y[0] + y[1]);
		}
		else {
			result = _basic_simpson(y, 0, N - 3, x, dx, axis);

			std::vector<Tx> h = { dx, dx };
			if (x) {
				h = diff(*x);
			}

			Tx h1 = h[h.size() - 1];
			Tx h0 = h[h.size() - 2];

			Tx alpha = (2.0 * h1 * h1 + 3.0 * h0 * h1) / (6.0 * (h0 + h1));
			Tx beta = (h1 * h1 + 3.0 * h0 * h1) / (6.0 * h0);
			Tx eta = (h1 * h1 * h1) / (6.0 * h0 * (h0 + h1));

			result += alpha * y[N - 1] + beta * y[N - 2] - eta * y[N - 3];
		}
		result += val;
	}
	else {
		result = _basic_simpson(y, 0, N - 2, x, dx, axis);
	}

	return result;
}


template <typename T, typename Tx>
T simpson_integration(const ::std::vector<T>& y, const ::std::vector<Tx>& x)
{
	return simpson(y, &x);
}

std::vector<std::complex<double>> laplace_transform_periodic(const std::vector<double>& data, const std::vector<double>& t, const std::vector<std::complex<double>>& s);
std::complex<double> laplace_transform_periodic(const std::vector<double>& data, const std::vector<double>& t, const std::complex<double>& s);
std::vector<Eigen::Matrix<std::complex<double>, 3, 3>> laplace_transform_periodic(const std::vector<mat3ds>& data, const std::vector<double>& t, const std::vector<std::complex<double>>& s);

Eigen::Matrix<double, 3, 3> convert_mat3ds_EigenMatrix(const mat3ds& data);

Eigen::Matrix<double, 3, 3> convert_mat3d_EigenMatrix(const mat3d& data);

::std::vector<double> inverse_laplace_transform(const std::vector<std::complex<double>>& beta_data, const std::vector<double>& t, const std::vector<std::complex<double>>& beta_s);

::std::vector<::std::vector<std::complex<double>>> createTalbotPath(const ::std::vector<double>& t, const int N, const double shift);

::std::vector<double> talbotInverseLaplaceTransform(const ::std::vector<::std::vector<::std::complex<double>>>& F, const ::std::vector<double>& t, const int N, const double shift);

void write_vector2D_to_csv(const::std::vector<::std::vector<double>>& data, const std::string& filename);

mat3ds PK2stress(FEIsotropicElastic& material, const FEMaterialPoint& mp);

void exit_febio();

// 序列化 std::complex<double>
template <>
inline DumpStream& DumpStream::operator<<(std::complex<double>& o) {
	if (m_btypeInfo) writeType(TypeID::TYPE_UNKNOWN); // 假设有一个对应于复数的类型ID

	// 序列化复数的实部和虚部
	double real = o.real();
	double imag = o.imag();
	m_bytes_serialized += write(&real, sizeof(double), 1);
	m_bytes_serialized += write(&imag, sizeof(double), 1);

	return *this;
}

// 反序列化 std::complex<double>
template <>
inline DumpStream& DumpStream::operator>>(std::complex<double>& o) {
	if (m_btypeInfo) readType(TypeID::TYPE_UNKNOWN); // 假设有一个对应于复数的类型ID

	// 反序列化复数的实部和虚部
	double real = 0.0, imag = 0.0;
	m_bytes_serialized += read(&real, sizeof(double), 1);
	m_bytes_serialized += read(&imag, sizeof(double), 1);
	o = std::complex<double>(real, imag);

	return *this;
}

::std::vector<::std::complex<double>> laplace_transform_inT(const ::std::vector<double>& data, const std::vector<double>& t, const std::vector<std::complex<double>>& s);

std::vector<Eigen::Matrix<std::complex<double>, 3, 3>> laplace_transform_inT(const std::vector<mat3ds>& data, const std::vector<double>& t, const std::vector<std::complex<double>>& s);


template<typename T>
FEElement* get_FEElement_p_version(T&& element_ref_p)
{
	FEElement* true_element = nullptr;
	// 如果element的类型是FEElement*
	if constexpr (::std::is_same_v<T, FEElement*&>)
	{
		true_element = element_ref_p;
	}
	else
	{
		true_element = element_ref_p.pe;
	}

	return true_element;
}

using ::std::span;

// -------------------- span2d --------------------
template<typename T>
struct span2d {
	T* base;
	size_t dim1, dim2; // dim1 = 行数, dim2 = 列数

	span2d() = default;
	span2d(T* b, size_t d1, size_t d2) : base(b), dim1(d1), dim2(d2) {}

	std::span<T> operator[](size_t i) const {
		return std::span<T>(base + i * dim2, dim2);
	}

	size_t rows() const { return dim1; }
	size_t cols() const { return dim2; }

	using element_type = T;

	// 迭代器 = 指针风格
	struct iterator {
		T* ptr;
		size_t stride;

		using iterator_category = std::random_access_iterator_tag;
		using value_type = std::span<T>;
		using difference_type = std::ptrdiff_t;
		using reference = value_type;

		reference operator*() const { return std::span<T>(ptr, stride); }
		reference operator[](difference_type n) const { return std::span<T>(ptr + n * stride, stride); }

		iterator& operator++() { ptr += stride; return *this; }
		iterator& operator--() { ptr -= stride; return *this; }
		iterator& operator+=(difference_type n) { ptr += n * stride; return *this; }
		iterator& operator-=(difference_type n) { ptr -= n * stride; return *this; }

		friend iterator operator+(iterator it, difference_type n) { it += n; return it; }
		friend iterator operator-(iterator it, difference_type n) { it -= n; return it; }
		friend difference_type operator-(const iterator& a, const iterator& b) {
			if (a.stride == 0) return 0;
			return (a.ptr - b.ptr) / static_cast<difference_type>(a.stride);
		}

		friend bool operator==(const iterator& a, const iterator& b) { return a.ptr == b.ptr; }
		friend bool operator!=(const iterator& a, const iterator& b) { return a.ptr != b.ptr; }
		friend bool operator<(const iterator& a, const iterator& b) { return a.ptr < b.ptr; }
		friend bool operator>(const iterator& a, const iterator& b) { return a.ptr > b.ptr; }
		friend bool operator<=(const iterator& a, const iterator& b) { return a.ptr <= b.ptr; }
		friend bool operator>=(const iterator& a, const iterator& b) { return a.ptr >= b.ptr; }
	};

	iterator begin() const {
		if (dim1 == 0 || dim2 == 0) return iterator{ base, 0 };
		return iterator{ base, dim2 };
	}
	iterator end() const {
		if (dim1 == 0 || dim2 == 0) return iterator{ base, 0 };
		return iterator{ base + dim1 * dim2, dim2 };
	}

	size_t size() const { return dim1; }
	bool empty() const { return dim1 == 0 || dim2 == 0; }
};

// -------------------- span3d --------------------
template<typename T>
struct span3d {
	T* base;
	size_t dim1, dim2, dim3; // dim1 = 第一维长度, dim2 = 第二维长度, dim3 = 第三维长度
	
	span3d() = default;
	span3d(T* b, size_t d1, size_t d2, size_t d3) : base(b), dim1(d1), dim2(d2), dim3(d3) {}

	span2d<T> operator[](size_t i) const {
		return span2d<T>{ base + i * dim2 * dim3, dim2, dim3 };
	}

	size_t dim_first()  const { return dim1; }
	size_t dim_second() const { return dim2; }
	size_t dim_third()  const { return dim3; }

	using element_type = T;

	// 迭代器 = 指针风格
	struct iterator {
		T* ptr;
		size_t stride;
		size_t dim2, dim3;

		using iterator_category = std::random_access_iterator_tag;
		using value_type = span2d<T>;
		using difference_type = std::ptrdiff_t;
		using reference = value_type;

		reference operator*() const { return span2d<T>{ ptr, dim2, dim3 }; }
		reference operator[](difference_type n) const { return span2d<T>{ ptr + n * stride, dim2, dim3 }; }

		iterator& operator++() { ptr += stride; return *this; }
		iterator& operator--() { ptr -= stride; return *this; }
		iterator& operator+=(difference_type n) { ptr += n * stride; return *this; }
		iterator& operator-=(difference_type n) { ptr -= n * stride; return *this; }

		friend iterator operator+(iterator it, difference_type n) { it += n; return it; }
		friend iterator operator-(iterator it, difference_type n) { it -= n; return it; }
		friend difference_type operator-(const iterator& a, const iterator& b) {
			if (a.stride == 0) return 0;
			return (a.ptr - b.ptr) / static_cast<difference_type>(a.stride);
		}

		friend bool operator==(const iterator& a, const iterator& b) { return a.ptr == b.ptr; }
		friend bool operator!=(const iterator& a, const iterator& b) { return a.ptr != b.ptr; }
		friend bool operator<(const iterator& a, const iterator& b) { return a.ptr < b.ptr; }
		friend bool operator>(const iterator& a, const iterator& b) { return a.ptr > b.ptr; }
		friend bool operator<=(const iterator& a, const iterator& b) { return a.ptr <= b.ptr; }
		friend bool operator>=(const iterator& a, const iterator& b) { return a.ptr >= b.ptr; }
	};

	iterator begin() const { 
		if (dim1 == 0 || dim2 == 0 || dim3 == 0) return iterator{ base, 0 , 0 };
		return iterator{ base, dim2 * dim3, dim2, dim3 }; 
	}
	iterator end()   const { 
		if (dim1 == 0 || dim2 == 0 || dim3 == 0) return iterator{ base, 0 , 0 };
		return iterator{ base + dim1 * dim2 * dim3, dim2 * dim3, dim2, dim3 }; 
	}

	size_t size() const { return dim1; }
	bool empty() const { return dim1 == 0 || dim2 == 0 || dim3 == 0; }
};

template<typename Iter>
auto span2d_to_vector2d(Iter begin, Iter end) {
	using T = Iter::value_type::element_type;

	std::vector<std::vector<T>> result;
	result.reserve(std::distance(begin, end)); // 提前分配行数

	for (auto it = begin; it != end; ++it) {
		auto row = *it;
		result.emplace_back(row.begin(), row.end());
	}

	return result;
}
template<typename T>
std::vector<std::vector<T>> span2d_to_vector2d(const span2d<T>& s) {
	return span2d_to_vector2d(s.begin(), s.end());
}

template<typename Iter>
auto span3d_to_vector3d(Iter begin, Iter end) {
	using T = Iter::value_type::element_type;

	std::vector<std::vector<std::vector<T>>> result;
	result.reserve(std::distance(begin, end));

	for (auto it = begin; it != end; ++it) {
		const auto plane = *it;
		std::vector<std::vector<T>> plane_vec;
		plane_vec.reserve(plane.dim1);

		for (size_t j = 0; j < plane.dim1; ++j) {
			auto row = plane[j]; // std::span<T>
			plane_vec.emplace_back(row.begin(), row.end()); // 拷贝一行
		}

		result.emplace_back(std::move(plane_vec));
	}
	return result;
}
template<typename T>
std::vector<std::vector<std::vector<T>>> span3d_to_vector3d(const span3d<T>& s) {
	return span3d_to_vector3d(s.begin(), s.end());
}

template<typename Iter>
auto span_to_vector(Iter begin, Iter end)
{
	using T = Iter::value_type;
	return std::vector<T>(begin, end);
}
template<typename T>
::std::vector<T> span_to_vector(const std::span<T>& s) {
	return span_to_vector(s.begin(), s.end());
}

