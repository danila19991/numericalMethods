#include "PLUMatrix.h"
#include <cassert>
#include <algorithm>

const double eps = 1e-3;

PLUMatrix::PLUMatrix(Matrix a):_a(std::move(a))
{
	assert(_a.numRows() != 0 || _a.numCols() != 0);

	const size_t tmpLen = std::min(_a.numCols(), _a.numRows());

	_p.reserve(tmpLen);
	_q.reserve(tmpLen);

	for(size_t i = 0; i < tmpLen;++i)
	{
		size_t posx = i, posy = i;
		double mx = abs(_a[i][i]);
		for(size_t j = i; j < _a.numRows();++j)
		{
			for(size_t k = i; k < _a.numCols(); ++k)
			{
				if(mx < abs(_a[j][k]))
				{
					mx = abs(_a[j][k]);
					posx = k;
					posy = j;
				}
			}
		}

		if (abs(mx) < 1e-9)
			break;

		_a.swapCols(i, posx);
		_a.swapRows(i, posy);

		_p.emplace_back(posx);
		_q.emplace_back(posy);

		for(size_t j = i + 1; j < _a.numRows();++j)
		{
			const double k = _a[j][i] / _a[i][i];
			for(size_t t = i; t < _a.numCols();++t)
			{
				_a[j][t] -= k * _a[i][t];
			}
			_a[j][i] = k;
		}
	}
}

Matrix PLUMatrix::getA() const noexcept
{
	Matrix a(_a.numRows(), _a.numCols()),tmp(_a.trans());
	
	for(int i=0;i<static_cast<int>(a.numRows()); ++i)
	{
		for(int j=0;j<static_cast<int>(a.numCols());++j)
		{
			for(int t=0;t<static_cast<int>(a.numCols()) && t<=i && t<=j;++t)
			{
				if (i == t)
					a[i][j] += tmp[j][t];
				else
					a[i][j] += _a[i][t] * tmp[j][t];
			}
		}
	}

	if (!_p.empty()) {
		for (int i = static_cast<int>(_p.size()) - 1;i >= 0;--i)
		{
			a.swapCols(i, _p[i]);
		}
		for (int i = static_cast<int>(_p.size()) - 1;i >= 0;--i)
		{
			a.swapRows(i, _q[i]);
		}
	}
	
	return a;
}

double PLUMatrix::getDeterminant() const noexcept
{
	assert(_a.numRows() == _a.numCols());
	double ans = 1;
	for(size_t i=0;i<_a.numCols(); ++i)
	{
		ans *= _a[i][i];
	}
	return ans;
}

size_t PLUMatrix::getRank() const noexcept
{
	size_t ans = 0;
	while (ans < _a.numCols() && ans < _a.numRows() && abs(_a[ans][ans]) > eps)
		++ans;
	return ans;
}


Matrix PLUMatrix::obrat() const noexcept
{
	assert(_a.numRows() == _a.numCols());

	assert(abs(getDeterminant()) > eps);

	std::vector<std::vector<double>> ans;
	ans.reserve(_a.numCols());

	for(size_t i=0;i<_a.numCols();++i)
	{
		std::vector<double> res(_a.numRows());
		res[i] = 1.;

		res = solve(res).getData();

		ans.emplace_back(res);
	}

	return Matrix(ans).trans();
}

bool PLUMatrix::hasSolution(MyVector b) const noexcept
{
	for (size_t i = 0;i < _q.size();++i)
		std::swap(b[i], b[_q[i]]);

	const size_t num = getRank();

	for(size_t i = num;i<_a.numRows();++i)
	{
		double res=0;
		for(size_t j = 0;j<num;++j)
		{
			res += _a[i][j] * b[j];
		}
		if (abs(res - b[i]) > eps)
			return false;
	}

	return true;
}

MyVector PLUMatrix::solve(MyVector b) const noexcept
{
	assert(_a.numRows() != 0u);

	assert(hasSolution(b));

	const size_t num = getRank();

	MyVector ans = MyVector(_a.numRows());
	
	for (size_t i = 0;i < _q.size();++i)
		std::swap(b[i], b[_q[i]]);
	
	for (size_t i = 0; i<num;++i)
	{

		double tmp = b[i];
		for (size_t j = 0;j<i;++j)
		{
			tmp -= ans[j] * _a[i][j];
		}
		ans[i] = tmp;
	}

	for (int i = static_cast<int>(num) - 1; i >= 0;--i)
	{
		double tmp = ans[i];
		for (int j = static_cast<int>(num) - 1;j>i;--j)
		{
			tmp -= ans[j] * _a[i][j];
		}
		ans[i] = tmp / _a[i][i];
	}

	if (!_p.empty()) {
		for (int i = static_cast<int>(_p.size()) - 1;i >= 0;--i)
			std::swap(ans[i], ans[_p[i]]);
	}

	return ans;
}

double PLUMatrix::getConditionNumber() const noexcept
{
	return getA().getNorm()*obrat().getNorm();
}

void PLUMatrix::printPLU(std::ostream& out) const noexcept
{
	for (auto& it : _p)
		out << it << ' ';
	out << '\n' << '\n';

	for (auto& it : _q)
		out << it << ' ';
	out << '\n' << '\n';

	out << _a << '\n';
}

double sqrtLess(const double a)
{
	if (a < 0.)
		return -sqrt(-a);
	return sqrt(a);
}

Matrix PLUMatrix::getRandomPositiveDifinite(const size_t size,const double minElem,const double maxElem)
{

	Matrix tmp = Matrix::getRandom(size, size, sqrtLess(minElem), sqrtLess(maxElem));
	
	PLUMatrix plutmp(tmp);

	while(plutmp.getDeterminant() == 0)
	{
		tmp = Matrix::getRandom(size, size, sqrtLess(minElem), sqrtLess(maxElem));
		plutmp = PLUMatrix(tmp);
	}

	return tmp * tmp.trans();
}

std::ostream& operator<<(std::ostream& out, const PLUMatrix& m)
{
	out << m.getA();
	return out;
}

std::istream& operator>>(std::istream& in, PLUMatrix& m)
{
	Matrix tmp;
	in >> tmp;
	m = PLUMatrix(tmp);
	return in;
}
