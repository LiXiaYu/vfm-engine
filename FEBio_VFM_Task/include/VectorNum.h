#pragma once
#include <vector>

template<typename T>
class VectorNum : public ::std::vector<T>
{
public:
	using ::std::vector<T>::vector;

	VectorNum operator + (const VectorNum& v) const
	{
		VectorNum r;
		r.resize(this->size());
		for (int i = 0; i < this->size(); ++i)
		{
			r[i] = (*this)[i] + v[i];
		}
		return r;
	}
};