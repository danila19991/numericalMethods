#include "QRMatrix.h"
#include <cassert>

QRMatrix::QRMatrix(Matrix m):_r(std::move(m))
{
	assert(_r.numRows() == _r.numCols());
	_q = Matrix::onnes(_r.numRows());

	for (size_t i = 0;i + 1 < _r.numRows();++i)
	{
		MyVector tmpU = _r.getCol(i), u(_r.numRows() - i);
		for(size_t j=1; j<=tmpU.size() - i;++j)
		{
			u[u.size() - j] = tmpU[tmpU.size() - j];
		}
		u[0] -= sqrt(u^u);

		Matrix loc = (u * u) / (-1*(u^u) / 2);

		for(size_t j =0 ; j<u.size();++j)
		{
			loc[j][j] += 1;
		}

		Matrix tmp(u.size());

		for(size_t k=0;k<loc.numCols();++k)
		{
			for (size_t t = 0;t < loc.numCols();++t)
			{
				for (size_t j = 0;j < loc.numRows();++j)
				{
					tmp[t][j] += loc[k][j] * _r[k + i][t + i];
				}
			}
		}

		tmp = tmp.trans();

		for(size_t j=0;j<tmp.numRows();++j)
		{
			for(size_t t=0;t<tmp.numCols();++t)
			{
				_r[j + i][t + i] = tmp[j][t];
			}
		}

		tmp = Matrix(_q.numRows(), u.size());

		for (size_t j = 0;j < _r.numRows();++j)
		{
			for (size_t t = 0;t < loc.numCols();++t)
			{
				for (size_t k = 0;k<_r.numCols();++k)
				{
					if (k >= i)
					{
						tmp[j][t] += _q[j][k] * loc[t][k - i];
					}
				}
			}
		}

		for (size_t j = 0;j<tmp.numRows();++j)
		{
			for (size_t t = 0;t<tmp.numCols();++t)
			{
				_q[j][t+i] = tmp[j][t];
			}
		}
	}
}

Matrix QRMatrix::getA() const noexcept
{
	return _q * _r;
}

MyVector QRMatrix::solve(MyVector b) const
{
	assert(b.size() == _q.numRows());

	b = _q.trans()*b;

	MyVector ans(b.size());

	for(int i=static_cast<int>(_r.numRows()) - 1;i>=0;--i)
	{
		double tmp = b[i];
		for(int j= static_cast<int>(_r.numRows()) - 1;j>=i;--j)
		{
			tmp -= ans[j] * _r[i][j];
		}
		ans[i] = tmp / _r[i][i];
	}

	return ans;
}

void QRMatrix::printQR(std::ostream& out) const 
{
	out << "Q Matrix:\n" << _q << "\nR Matrix:\n" << _r << '\n';
}

std::istream& operator>>(std::istream& in, QRMatrix& m)
{
	Matrix tmp;
	in >> tmp;
	m = QRMatrix(tmp);
	return in;
}

std::ostream& operator<<(std::ostream& out, QRMatrix& m)
{
	out << m.getA();
	return out;
}
