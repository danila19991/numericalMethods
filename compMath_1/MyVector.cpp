#include "MyVector.h"
#include <cassert>
#include <random>
#include <chrono>
#include <cmath>

MyVector::MyVector(const size_t size):
				   _data(size)
{
}

MyVector::MyVector(std::vector<double> data):_data(std::move(data))
{
}

double& MyVector::operator[](size_t pos) noexcept
{
	assert(pos < _data.size());
	return _data[pos];
}

const double& MyVector::operator[](size_t pos) const noexcept
{
	assert(pos < _data.size());
	return _data[pos];
}

double MyVector::getNorm() const
{
	double ans = 0.;

	for(auto& elem:_data)
	{
		if (std::isnan(elem))
			return std::numeric_limits<double>::infinity();
		if (ans < abs(elem))
			ans = abs(elem);
	}

	return ans;
}

size_t MyVector::size() const noexcept
{
	return _data.size();
}

bool MyVector::isZero(double eps) const
{
	assert(eps >= 0.);

	for(auto& elem:_data)
	{
		if (abs(elem) < eps)
			return false;
	}

	return true;
}

std::vector<double> MyVector::getData() const
{
	return _data;
}

MyVector MyVector::getRand(const size_t dim,const double minElem,const double maxElem)
{
	assert(maxElem > minElem);
	assert(abs(maxElem) < 100'000.);
	assert(abs(minElem) < 100'000.);

	MyVector ans(dim);
	std::mt19937 gener(std::chrono::system_clock::now().time_since_epoch().count());

	for(auto& elem:ans._data)
	{
		elem = (gener() % (static_cast<int>((maxElem - minElem)*10'000.))) / 10'000. + minElem;
	}

	return ans;
}

MyVector operator*(const Matrix& m, const MyVector& v)
{
	assert(m.numCols() == v.size());
	MyVector ans(m.numRows());
	for(size_t i=0;i<m.numRows();++i)
	{
		for(size_t j=0;j<m.numCols();++j)
		{
			ans[i] += m[i][j] * v[j];
		}
	}
	return ans;
}

Matrix operator*(const MyVector& l, const MyVector& r)
{
	Matrix ans(r.size(), l.size());

	for(size_t i = 0; i<r.size();++i)
	{
		for(size_t j=0;j<l.size();++j)
		{
			ans[i][j] = l[j] * r[i];
		}
	}

	return ans;
}

double operator^(const MyVector& l, const MyVector& r)
{
	double ans = 0.;

	for (size_t i = 0;i < l.size();++i)
		ans += l[i] * r[i];

	return ans;
}

MyVector operator-(const MyVector& l, const MyVector& r)
{
	assert(l.size() == r.size());

	MyVector ans(l.size());

	for (size_t i = 0;i < l.size(); ++i)
		ans[i] = l[i] - r[i];

	return ans;
}

MyVector operator+(const MyVector& l, const MyVector& r)
{
	assert(l.size() == r.size());

	MyVector ans(l.size());

	for (size_t i = 0;i < l.size(); ++i)
		ans[i] = l[i] + r[i];

	return ans;
}

MyVector operator/(const MyVector& vec, double k)
{
	MyVector ans = vec;

	for (size_t i = 0;i < ans.size();++i)
		ans[i] /= k;

	return ans;
}

std::ostream& operator <<(std::ostream& out, const MyVector& vec)
{
	out << vec.size() << '\n';
	for(auto& it:vec._data)
	{
		out << it << ' ';
	}
	return out;
}

std::istream& operator >>(std::istream& in, MyVector& vec)
{
	size_t sz;
	in >> sz;
	vec = MyVector(sz);
	for (auto& elem : vec._data)
	{
		in >> elem;
	}
	return in;
}