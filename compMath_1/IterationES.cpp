#include "IterationES.h"
#include <cassert>
#include <algorithm>

void ReconstructES(Matrix& m, MyVector& vec, double eps)
{
	for(size_t i=0;i<m.numRows();++i)
	{
		if (abs(m[i][i]) > eps) 
		{
			for (size_t j = 0;j < m.numCols();++j)
			{
				if (i != j)
				{
					m[i][j] /= -1.*m[i][i];
				}
			}

			vec[i] /= m[i][i];
		}

		m[i][i] = 0.;
	}
}

MyVector Jakobi(Matrix a, MyVector b,bool& flag, const double eps)
{
	assert(a.numRows() == a.numCols());
	assert(a.numRows() == b.size());

	ReconstructES(a, b, eps);

	MyVector prev(b);
	MyVector nxt = a*prev + b;

	int it=0;

	const double c = a.getNorm();

	while(abs((nxt-prev).getNorm()*c) > abs((1-c)*eps))
	{
		prev = nxt;
		nxt = a * prev + b;
		++it;

		if (it > 100'000)
			break;
	}

	if (it > 100'000)
	{
		std::cout << "Jakobi is not work\n\n";
		flag = false;
		return MyVector(b.size());
	}


	std::cout << "Jakobi in " << it << " steps.\n\n";
	flag = true;
	return nxt;
}

MyVector makeNextInZeidel(const Matrix& c, const MyVector& b, const MyVector& prev)
{
	MyVector ans(b);

	for(size_t i=0;i<b.size();++i)
	{
		for(size_t j = 0;j<i;++j)
		{
			ans[i] += ans[j] * c[i][j];
		}
		for(size_t j=i+1;j<b.size();++j)
		{
			ans[i] += prev[j] * c[i][j];
		}
	}
	return ans;
}

MyVector Zeidel(Matrix a, MyVector b,bool& flag,const double eps)
{
	assert(a.numRows() == a.numCols());
	assert(a.numRows() == b.size());

	Matrix tmpa = a;
	MyVector tmpb = b;

	ReconstructES(a, b, eps);

	MyVector prev(b);
	MyVector nxt = makeNextInZeidel(a, b, prev);

	int it = 0;

	double c = 0.;

	for(size_t i=0;i<a.numRows();++i)
	{
		double tmp = 0.;
		for (size_t j = i + 1;j < a.numRows();++j)
			tmp += abs(a[i][j]);
		c = std::max(tmp, c);
	}

	while (abs((nxt - prev).getNorm()*c) > abs((1-c)*eps))
	{
		prev = nxt;
		nxt = makeNextInZeidel(a, b, prev);
		++it;

		if (it > 100'000)
			break;
	}

	if (it > 100'000) 
	{
		std::cout << "Zeidel is not work\n\nAx-b:\n" << tmpa*nxt-tmpb<<"\n\n";
		flag = false;
		return MyVector(b.size());
	}
	
	std::cout << "Zeidel in " << it << " steps.\n\n";
	flag = true;
	return nxt;
}
