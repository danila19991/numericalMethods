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

		Matrix tmp = Matrix::onnes(_r.numRows());

		Matrix loc = (u * u) / ((u^u) / 2);

		for (size_t j = 1;j <= loc.numRows();++j)
		{
			for (size_t t = 1;t <= loc.numCols();++t)
			{
				tmp[tmp.numRows() - j][tmp.numCols() - t] -=
					loc[loc.numRows() - j][loc.numCols() - t];
			}
		}

		_r = tmp * _r;
		_q = _q * tmp;
		//todo:
		//optimise multiplying on matrix with many zeroes.
	}
}

Matrix QRMatrix::getA() const
{
	return _q * _r;
}

MyVector QRMatrix::solve(MyVector b)
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
