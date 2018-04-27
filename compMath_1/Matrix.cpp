#include "Matrix.h"
#include <cassert>
#include <algorithm>
#include "MyVector.h"
#include <random>
#include <chrono>

Matrix::Matrix(const size_t size):
			   _data(size,std::vector<double>(size))
{
}

Matrix::Matrix(const size_t numRows,const size_t numCols):
			   _data(numRows,std::vector<double>(numCols))
{
}

Matrix::Matrix(std::vector<std::vector<double>> data):_data(std::move(data))
{
}

std::vector<double>& Matrix::operator[](const size_t pos) noexcept
{
	return _data[pos];
}

const std::vector<double>& Matrix::operator[](const size_t pos) const noexcept
{
	return _data[pos];
}

double Matrix::getNorm() const
{
	double ans = 0;

	for(auto& row:_data)
	{
		double tmp = 0;
		for (auto& elem : row)
			tmp += abs(elem);
		ans = std::max(ans, tmp);
	}

	return ans;
}

std::vector<double> Matrix::getCol(const size_t pos) const noexcept
{
	assert(pos >= 0 && pos < numCols());
	std::vector<double> ans(numRows());

	for (size_t i = 0;i < numRows();++i)
		ans[i] = _data[i][pos];

	return ans;
}

size_t Matrix::numRows() const noexcept
{
	return _data.size();
}

size_t Matrix::numCols() const noexcept
{
	assert(!_data.empty());
	return _data.begin()->size();
}

void Matrix::swapRows(const size_t l, const size_t r) noexcept
{
	assert(l < numRows() && r < numRows());
	if (l != r) 
	{
		std::swap(_data[l], _data[r]);
	}
}

void Matrix::swapCols(const size_t l, const size_t r) noexcept
{
	assert(l < numCols() && r < numCols());
	if(l!=r)
	{
		for (size_t i = 0;i < numRows();++i) 
		{
			std::swap(_data[i][l], _data[i][r]);
		}
	}
}

void Matrix::print(const size_t pres,std::ostream& out) const noexcept
{
	outputTable(_data, pres, out);
}

Matrix Matrix::trans() const noexcept
{
	Matrix ans(numCols(),numRows());
	for(size_t i=0;i<numRows();++i)
	{
		for(size_t j=0;j<numCols();++j)
		{
			ans[j][i] = _data[i][j];
		}
	}
	return ans;
}

bool Matrix::isZero(const double eps) const
{
	assert(eps >= 0.);

	for(auto& row:_data)
	{
		for(auto& elem:row)
		{
			if (abs(elem) < eps)
				return false;
		}
	}

	return true;
}

Matrix Matrix::onnes(const size_t dim)
{
	Matrix ans(dim);

	for (size_t i = 0;i < dim;++i)
		ans[i][i] = 1;

	return ans;
}

Matrix Matrix::getRandom(const size_t numRows,const size_t numCols,const double minElem,
						 const double maxElem)
{
	assert(numRows >= 0);
	assert(numCols >= 0);
	assert(minElem < maxElem);

	Matrix ans;
	ans._data.reserve(numRows);
	
	for(size_t i=0; i<numRows;++i)
	{
		ans._data.emplace_back(MyVector::getRand(numCols, minElem, maxElem).getData());
	}

	return ans;
}

Matrix Matrix::getRandomWithDiagonalAcpuisition(const size_t size,const double minElem,
	const double maxElem)
{
	assert(size >= 1);
	assert(minElem < maxElem);

	const Matrix tmp = getRandom(size, size - 1, minElem, maxElem);
	Matrix ans(size);

	std::mt19937 gener(std::chrono::system_clock::now().time_since_epoch().count());

	for(size_t i=0;i<size;++i)
	{
		double elemSum = 0.;
		for(size_t j=0; j+1<size;++j)
		{
			elemSum += abs(tmp[i][j]);
			if (j < i)
				ans[i][j] = tmp[i][j];
			else
				ans[i][j + 1] = tmp[i][j];
		}
		elemSum += (gener() % (static_cast<int>((maxElem - std::max(minElem,0.))*10'000.))) /
					10'000. + std::max(minElem,0.);
		if (gener() & 1)
			ans[i][i] = elemSum;
		else
			ans[i][i] = -elemSum;
	}

	return ans;
}

Matrix operator^(const Matrix& l, const Matrix& r)
{
	Matrix ans(l.numRows(), r.numCols()), r1 = r.trans();

	for (size_t i = 0;i < ans.numRows(); ++i)
	{
		for (size_t j = 0; j < ans.numCols(); ++j)
		{
			for (size_t k = 0; k < l.numCols() && k < r1.numCols(); ++k)
			{
				ans[i][j] += l[i][k] * r1[j][k];
			}
		}
	}

	return ans;
}

Matrix operator*(const Matrix& l,const Matrix& r)
{
	assert(l.numCols() == r.numRows());
	Matrix ans(l.numRows(),r.numCols()), r1 = r.trans();

	for(size_t i = 0;i < ans.numRows(); ++i)
	{
		for(size_t j = 0; j < ans.numCols(); ++j)
		{
			for(size_t k = 0; k < l.numCols(); ++k)
			{
				ans[i][j] += l[i][k] * r1[j][k];
			}
		}
	}

	return ans;
}

Matrix operator/(const Matrix& vec, double k)
{
	Matrix ans = vec;

	for (size_t i=0; i<ans.numRows();++i)
	{
		for(size_t j=0;j<ans.numCols();++j)
		{
			ans[i][j] /= k;
		}
	}

	return ans;
}

std::istream& operator>>(std::istream& in, Matrix& m)
{
	size_t szRows, szCols;
	in >> szRows>>szCols;
	m = Matrix(szRows,szCols);
	
	for (auto& row : m._data)
	{
		for (auto& elem : row)
		{
			in >> elem;
		}
	}

	return in;
}

std::ostream& operator<<(std::ostream& out, const Matrix& m)
{
	out << m.numRows() << ' ' << m.numCols() << '\n';
	for(auto& row:m._data)
	{
		for(auto& elem: row)
		{
			out << elem << ' ';
		}
		out << '\n';
	}

	return out;
}
