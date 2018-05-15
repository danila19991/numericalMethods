#include "Integrals.h"
#include <cassert>

#include "NewtonMethods.h"
#include "QRMatrix.h"

double RectanglesMethod(double(* func)(double x), double l,double r, double len)
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
 * \brief		Calulates vector of weghts for \f$ p(x) = (b-x)^(k)\f$
 * \param[in] l Left side of segment.
 * \param[in] r Right side of segmant.
 * \param[in] b Param of weght function.
 * \param[in] k Power for weght function.
 * \param[in] x Vector of points for calculating.
 * \return		Vector of coeficients A.
 */
MyVector makeWeighsIQA(double l, double r, double b, double k, const MyVector& x)
{
	assert(k != 1);
	assert(l < r);

	const MyVector m = makeM(l, r, b, k, x.size());

	Matrix A(x.size());

	for(size_t i=0; i<A.numRows();++i)
	{ 
		A[i][0] = 1.;
		for (size_t j = 1;j < A.numCols();++j)
		{
			A[i][j] = A[i][j - 1] * x[i];
		}
	}

	A = A.trans();

	return static_cast<QRMatrix>(A).solve(m);
}

/**
* \brief			Calulates vector of error for \f$ p(x) = (b-x)^(k)\f$
* \param[in] maxDer	Maximum absolute value of len's derive of func.
* \param[in] l		Left side of segment.
* \param[in] r		Right side of segmant.
* \param[in] k		Power for weght function.
* \return			Error of integral.
*/
double makeErrorIQA(const double maxDer, const double l, const double r, const double k)
{
	assert(k != 1);
	assert(l < r);

	const MyVector m = makeM(l, r, r, k, 4);

	MyVector m2(4);
	for(size_t i=0; i<4;++i)
	{
		m2[i] = -calcMInPoint(i, r, (l + r) / 2., 0.25);
	}

	const MyVector c = MyVector({ -l * r*(l + r) / 2. , l*r + (l + r)*(l + r) / 2,
		-3. / 2.*(l + r),1. });

	return (abs((m2+m)^c) + abs(m2^c))*maxDer/6;
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
double IQAIntagralLoc(double(*func)(double x), const double l, const double r, const double b,
	const double k)
{
	const MyVector x = MyVector({ l, (l + r) / 2., r });

	const MyVector a = makeWeighsIQA(l, r, b, k, x);

	MyVector fx = MyVector(x.size());

	for (size_t i = 0; i < fx.size(); i++)
		fx[i] = func(x[i]);

	return a ^ fx;
}

std::pair<double, double> IQAIntagral(double (* func)(double x), const double l, const double r,
                                      const double b, const double k)
{
	return std::make_pair(IQAIntagralLoc(func,l,r,b,k), makeErrorIQA(37.146, l, r, k));
}

double IntagralNewtonKotsStep(double (* func)(double x), double l, double r, const double b,
                              const double k, double len)
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

double IntagralNewtonKots(double(* func)(double x), double l, double r, double b, double k,
											 double eps)
{
	const double ans = IntagralNewtonKotsStep(func, l, r, b, k, 0.1);

	const double ans2 = IntagralNewtonKotsStep(func, l, r, b, k, 0.05);

	const double len = 0.045*pow(eps / abs(ans2 - ans) * 7., 1. / 3.);

	return IntagralNewtonKotsStep(func, l, r, b, k, len);
}

/**
 * \brief			Method for founding 3 rational solution of polynom.
 * \param[in] a		Coefficient of polynom with degree 3.(a_n =1)
 * \param[in] l		Left side of segment for searching roots.
 * \param[in] r		Rieght side of segment for searching root.
 * \param[in] eps	Error for foumding root.
 * \return			Vector of 3 roots.
 */
MyVector NewtonPolynomSol(MyVector a, const double l, const double r, const double eps)
{
	assert(a.size() == 3);

	const auto f = [=](double x)
	{
		double ans = 1;
		for(int i=2;i>=0;--i)
		{
			ans = ans * x + a[i];
		}
		return ans;
	};

	const auto d = [=](double x)
	{
		return 3*x*x + 2*x*a[2] + a[1];
	};

	const double p = a[2] / 3, q = a[1] / 3;

	const double des = p * p - q;

	assert(des > 0);

	const double x1 = -p - sqrt(des), x2 = -p + sqrt(des);

	return MyVector({ Newton(f,d,l,x1,eps),Newton(f,d,x1,x2,eps) ,Newton(f,d,x2,r,eps) });
}

/**
 * \brief			Method for founding weights for calculating integral by gaus method.
 * \param[in] l		Left side of segment.
 * \param[in] r		Rieght side of segment.
 * \param[in] b		Parametr of weight function.
 * \param[in] k		Parametr of weight function.
 * \param[in] len	Number of wieghts.
 * \return			Pair of weights and points where function should be found.
 */
std::pair<MyVector, MyVector> makeAGaus(const double l, const double r, const double b,
                                        const double k, const size_t len)
{
	MyVector m = makeM(l, r, b, k, len * 2);

	MyVector m2(len);

	for(size_t i=0;i<m2.size();++i)
	{
		m2[i] = -m[i + len];
	}

	MyVector m1(len);

	for (size_t i = 0;i<m1.size();++i)
	{
		m1[i] = m[i];
	}
	Matrix B(len);

	for(size_t i=0;i<B.numRows();++i)
	{
		for(size_t j=0;j<B.numCols();++j)
		{
			B[i][j] = m[i + j];
		}
	}

	const MyVector a2 = static_cast<QRMatrix>(B).solve(m2);

	
	const MyVector a = NewtonPolynomSol(a2, l, r, 1e-6);


	Matrix A(len);

	for (size_t i = 0; i<A.numRows();++i)
	{
		A[i][0] = 1.;
		for (size_t j = 1;j < A.numCols();++j)
		{
			A[i][j] = A[i][j - 1] * a[i];
		}
	}

	A = A.trans();

	return std::make_pair(static_cast<QRMatrix>(A).solve(m1),a);
}

/**
 * \brief			Method for multipluing two polynoms.
 * \param[in] m1	Coefficients of first polynoms.
 * \param[in] m2	Coefficients of second polynom.
 * \return			Coeeficient of mltiplyed polynoms.
 */
MyVector polyMult(const MyVector& m1, const MyVector& m2)
{
	MyVector ans(m1.size() + m2.size() - 1);

	for(size_t i=0;i<m1.size();++i)
	{
		for(size_t j=0; j<m2.size();++j)
		{
			ans[i+j] += m1[i] * m2[j];
		}
	}

	return ans;
}

/**
 * \brief				Method for calculating method error of Gaus.
 * \param[in] maxDer	Absolutely max value of 2n dervate.
 * \param[in] l			Left side of segmant.
 * \param[in] r			Rieght side of segmant.
 * \param[in] b			Parametr of weight function.
 * \param[in] k			Parametr of weight function.
 * \param[in] x			Points where function was calculated in Gaus method.
 * \return 
 */
double makeGausError(const double maxDer, const double l, const double r, const double b, 
	const double k, const MyVector& x)
{
	assert(k != 1);
	assert(l < r);
	assert(x.size() == 3);

	const MyVector m = makeM(l, r, b, k, 7);

	MyVector c1 = MyVector({ -x[0]*x[1]*x[2] , x[0]*x[1] + x[0]*x[2] + x[1]*x[2],
		-x[2] - x[1] - x[0],1. });

	c1 = polyMult(c1, c1);

	return (abs(m^c1))*maxDer / 720.;
	//690.2
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
std::pair<double, MyVector> GausIntagralLoc(double (*func)(double x), const double l, 
	const double r, const double b, const double k)
{
	const auto A = makeAGaus(l, r, b, k, 3);

	MyVector fx = MyVector(A.second);

	for (size_t i = 0; i < fx.size(); i++)
		fx[i] = func(A.second[i]);

	return std::make_pair(A.first ^ fx,A.second);
}

std::pair<double, double> GausIntagral(double(* func)(double x), const double l, const double r, 
	const double b, const double k)
{
	const auto a = GausIntagralLoc(func, l, r, b, k);

	return std::make_pair(a.first, makeGausError(690.2, l, r, b, k, a.second));
}

double GausIntagralPartialStep(double(*func)(double x), double l, double r, const double b, 
	const double k, double len)
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

double GausIntagral(double(* func)(double x), const double l, const double r, const double b, 
	const double k, const double eps)
{
	const double ans = GausIntagralPartialStep(func, l, r, b, k, 3);

	const double ans2 = GausIntagralPartialStep(func, l, r, b, k, 1.5);

	const double len = 1.5*pow(eps/ abs(ans2 - ans) * 7., 1. / 3.);

	std::cout << "Length of step in Gaus " << len << '\n';

	return GausIntagralPartialStep(func, l, r, b, k, len);
}
