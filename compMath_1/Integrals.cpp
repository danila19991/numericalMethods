#include "Integrals.h"
#include <cassert>

#include "NewtonMethods.h"
#include "QRMatrix.h"

double RectanglesMethod(const std::function<double(double x)>& func, double l,double r, double len)
{
	assert(len > 0.);

	if (l > r)
		std::swap(l, r);

	const long long n = static_cast<long long>((r - l) / len) + 1;

	len = (r - l) / static_cast<double>(n);

	double ans = 0.;

	for(long long i=1;i<n;++i)
	{
		const double x = l + i * len - len / 2.;
		
		ans += len * func(x);
	}

	return ans;
}

/**
 * \brief		Method for calculating a^k(where k is natural)
 * \param[in] a Number of powering.
 * \param[in] d Degree.
 * \return		Result of powering.
 */
double binPow(const double a,const size_t d)
{
	if (d == 0)
		return 1.;
	double b = binPow(a, d / 2);
	b = b * b;
	if (d % 2 == 1)
		b = b*a;
	return b;
}

/**
 * \brief		Calculating intergal of p(x) = (b-x)^(-k) in point x.
 * \param[in] s Parametr of p(x).
 * \param[in] b Parametr og p(x).
 * \param[in] x Point for caluclating.
 * \param[in] k Parametr of p(x).
 * \return		Result of integral in point.
 */
double calcMInPoint(const size_t s, const double b, const double x, const double k)
{
	double ans = 0.;

	double c = 1.;

	double d = pow(b - x, 1 - k);

	double a = binPow(b, s);

	for (size_t i = 0;i <= s;++i)
	{
		if (i % 2 == 0)
			ans += c / (i - k + 1)*d*a;
		else
			ans -= c / (i - k + 1)*d*a;

		c = c / (i + 1) * (s - i);
		d = d * (b - x);
		a = a / b;
	}

	return ans;
}

/**
 * \brief			Calulates vector of moments for \f$ p(x) = (b-x)^(k)\f$
 * \param[in] l		Left side of segment.
 * \param[in] r		Right side of segmant.
 * \param[in] b		Param of weght function.
 * \param[in] k		Power for weght function.
 * \param[in] len	Number of points.
 * \return			Vector of coeficients m.
 */
MyVector makeM(const double l, const double r, const double b, const double k, const size_t len)
{
	assert(k != 1);
	assert(l < r);

	MyVector m(len);

	m[0] = pow(b - l, 1 - k) / (1 - k) - pow(b - r, 1 - k) / (1 - k);

	for (size_t s = 1;s<m.size();++s)
	{
		m[s] = calcMInPoint(s, b, l, k);

		if(r!= b)
		{
			m[s] -= calcMInPoint(s, b, r, k);
		}
	}

	return m;
}

/**
 * \brief			Method for calculating IQAIntegral on segment.
 * \param[in] func	Function for integrating.
 * \param[in] l		Left side of segment.
 * \param[in] r		Rieght side of segment.
 * \param[in] b		Parametr of weight function.
 * \param[in] k		Parametr of weight functron.
 * \return 
 */
double IQAIntagralLoc(const std::function<double(double x)>& func, const double l, const double r,
	const double b, const double k)
{
	assert(k != 1);
	assert(l < r);

	const MyVector x = MyVector({ l, (l + r) / 2., r });

	const MyVector m = makeM(l, r, b, k, x.size());

	Matrix A(x.size());

	for (size_t i = 0; i<A.numRows();++i)
	{
		A[i][0] = 1.;
		for (size_t j = 1;j < A.numCols();++j)
		{
			A[i][j] = A[i][j - 1] * x[i];
		}
	}

	A = A.trans();

	const MyVector a = static_cast<QRMatrix>(A).solve(m);

	MyVector fx = MyVector(x.size());

	for (size_t i = 0; i < fx.size(); i++)
		fx[i] = func(x[i]);

	return a ^ fx;
}

std::pair<double, double> IQAIntagral(const std::function<double(double x)>& func, const double l,
	const double r, const double b, const double k, const double maxDer)
{
	assert(k != 1);
	assert(l < r);

	const MyVector m = makeM(l, r, r, k, 4);

	MyVector m2(4);
	for (size_t i = 0; i<4;++i)
	{
		m2[i] = -calcMInPoint(i, r, (l + r) / 2., 0.25);
	}

	const MyVector c = MyVector({ -l * r*(l + r) / 2. , l*r + (l + r)*(l + r) / 2,
		-3. / 2.*(l + r),1. });

	return std::make_pair(IQAIntagralLoc(func,l,r,b,k), (abs((m2 + m) ^ c) + abs(m2^c))*maxDer / 6);
}

double IntagralNewtonKotsPartial(const std::function<double(double x)>& func, double l, double r,
	const double b, const double k, double len)
{
	assert(len > 0.);

	if (l > r)
		std::swap(l, r);

	const long long n = static_cast<long long>((r - l) / len) + 1;

	len = (r - l) / static_cast<double>(n);

	double ans = 0.;

	for (long long i = 0;i<n;++i)
	{
		const double l1 = l + i * len;

		const double r1 = l + (i + 1)*len;

		const double t = IQAIntagralLoc(func, l1, r1, b, k);

		ans += t;
	}

	return ans;
}

double IntagralNewtonKotsAcurate(const std::function<double(double x)>& func, const double l,
	const double r, const double b, const double k, const double eps)
{
	const double ans = IntagralNewtonKotsPartial(func, l, r, b, k, 0.1);

	const double ans2 = IntagralNewtonKotsPartial(func, l, r, b, k, 0.05);

	const double len = 0.045*pow(eps / abs(ans2 - ans) * 7., 1. / 3.);

	return IntagralNewtonKotsPartial(func, l, r, b, k, len);
}

/**
 * \brief			Method for calcualting integral by Gaus method.
 * \param[in] func	Function for calculating.
 * \param[in] l		Left side of segment.
 * \param[in] r		Rieght siede of segment.
 * \param[in] b		Parametr of weight function.
 * \param[in] k		Parametr of weight function.
 * \return 
 */
std::pair<double, MyVector> GausIntagralLoc(const std::function<double(double x)>& func,
	const double l, const double r, const double b, const double k)
{
	MyVector m = makeM(l, r, b, k, 6);

	MyVector m2(3);

	for (size_t i = 0;i<m2.size();++i)
	{
		m2[i] = -m[i + 3];
	}

	MyVector m1(3);

	for (size_t i = 0;i<m1.size();++i)
	{
		m1[i] = m[i];
	}
	Matrix B(3);

	for (size_t i = 0;i<B.numRows();++i)
	{
		for (size_t j = 0;j<B.numCols();++j)
		{
			B[i][j] = m[i + j];
		}
	}

	const MyVector a2 = static_cast<QRMatrix>(B).solve(m2);


	const auto f = [=](double x)
	{
		double ans = 1;
		for (int i = 2;i >= 0;--i)
		{
			ans = ans * x + a2[i];
		}
		return ans;
	};

	const auto d = [=](double x)
	{
		return 3 * x*x + 2 * x*a2[2] + a2[1];
	};

	const double p = a2[2] / 3, q = a2[1] / 3;

	const double des = p * p - q;

	assert(des > 0);

	const double x1 = -p - sqrt(des), x2 = -p + sqrt(des);

	const MyVector a= MyVector({ Newton(f,d,l,x1,1e-6),Newton(f,d,x1,x2,1e-6) ,
		Newton(f,d,x2,r,1e-6) });


	Matrix A(3);

	for (size_t i = 0; i<A.numRows();++i)
	{
		A[i][0] = 1.;
		for (size_t j = 1;j < A.numCols();++j)
		{
			A[i][j] = A[i][j - 1] * a[i];
		}
	}

	A = A.trans();

	MyVector fx = MyVector(a);

	for (size_t i = 0; i < fx.size(); i++)
		fx[i] = func(a[i]);

	return std::make_pair(static_cast<QRMatrix>(A).solve(m1) ^ fx,a);
}

std::pair<double, double> GausIntagral(const std::function<double(double x)>& func, const double l,
	const double r, const double b, const double k, const double maxDer)
{
	assert(k != 1);
	assert(l < r);
	
	const auto a = GausIntagralLoc(func, l, r, b, k);

	const MyVector m = makeM(l, r, b, k, 7);

	MyVector c1 = MyVector({ -a.second[0] * a.second[1] * a.second[2],
		a.second[0] * a.second[1] + a.second[0] * a.second[2] + a.second[1] * a.second[2],
		-a.second[2] - a.second[1] - a.second[0], 1. });

	MyVector c = MyVector(7);

	for (size_t i = 0;i<c1.size();++i)
	{
		for (size_t j = 0; j<c1.size();++j)
		{
			c[i + j] += c1[i] * c1[j];
		}
	}

	return std::make_pair(a.first, (abs(m^c))*maxDer / 720.);
}

double GausIntagralPartial(const std::function<double(double x)>& func, double l, double r,
	const double b, const double k, double len)
{
	assert(len > 0.);

	if (l > r)
		std::swap(l, r);

	const long long n = static_cast<long long>((r - l) / len) + 1;

	len = (r - l) / static_cast<double>(n);

	double ans = 0.;

	for (long long i = 0;i<n;++i)
	{
		const double l1 = l + i * len;

		const double r1 = l + (i + 1)*len;

		const double t = GausIntagralLoc(func, l1, r1, b, k).first;

		ans += t;
	}

	return ans;
}

double GausIntagralAcurate(const std::function<double(double x)>& func, const double l,
	const double r, const double b, const double k, const double eps)
{
	const double ans = GausIntagralPartial(func, l, r, b, k, 3);

	const double ans2 = GausIntagralPartial(func, l, r, b, k, 1.5);

	const double len = 1.5*pow(eps/ abs(ans2 - ans) * 7., 1. / 3.);

	std::cout << "Length of step in Gaus " << len << '\n';

	return GausIntagralPartial(func, l, r, b, k, len);
}
