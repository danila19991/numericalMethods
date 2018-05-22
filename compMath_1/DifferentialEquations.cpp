#include "DifferentialEquations.h"
#include <cassert>
#include <algorithm>

/**
 * \brief			Calculates y_n+1 using modified Euler metod.
 * \param[in] func	System of functions for calculating.
 * \param[in] x0	Current x value.
 * \param[in] y0	Current y value.
 * \param[in] len	Length of step.
 * \param[in] f0	Value of function in x0.
 * \return			Value of y_n+1.
 */
MyVector DifferentialEquationsOppStepAfterCalc(
	const std::function<MyVector(double, const MyVector&)>& func, const double x0, 
	const MyVector& y0, const double len, const MyVector& f0)
{
	return y0 + len / 2.*(f0 + func(x0 + len, y0 + len * f0));
}

/**
 * \brief			Calculates y_n+1 using modified Euler metod.
 * \param[in] func	System of functions for calculating.
 * \param[in] x0	Current x value.
 * \param[in] y0	Current y value.
 * \param[in] len	Length of step.
 * \return			Value of y_n+1.
 */
MyVector DifferentialEquationsOppStep(
	const std::function<MyVector(double, const MyVector&)>& func, const double x0,
	const MyVector& y0, const double len)
{
	return DifferentialEquationsOppStepAfterCalc(func, x0, y0, len, func(x0, y0));
}

std::vector<std::pair<double, MyVector>> DifferentialEquationsOpp(
	const std::function<MyVector(double, const MyVector&)>& func,MyVector y0, double len,
	const double l, const double r)
{
	assert(l < r);
	const size_t n = (r - l) / len;
	len = (r - l) / n;

	std::vector<std::pair<double, MyVector>> ans;
	ans.reserve(n+1);

	ans.emplace_back(std::make_pair(l,y0));
	for(size_t i=0;i<n;++i)
	{
		const double x = l + i * len;

		y0 = DifferentialEquationsOppStep(func, x, y0, len);

		ans.emplace_back(std::make_pair(x + len, y0));
	}

	return ans;
}

/**
 * \brief			Calculates y_n+1 using modified Euler metod.
 * \param[in] func	System of functions for calculating.
 * \param[in] x0	Current x value.
 * \param[in] y0	Current y value.
 * \param[in] len	Length of step.
 * \param[in] a		Matrix with constants
 * \param[in] b		Vector with constants.
 * \param[in] c		Vector with constants.
 * \param[in] k0	Value of function in x0.
 * \return			Value of y_n+1.
 */
MyVector RungeKuttStepAgterCalc(const std::function<MyVector(double, const MyVector&)>& func,
	const double x0, MyVector y0, const double len, const Matrix&  a, const MyVector& b,
	const MyVector& c, const MyVector& k0)
{
	assert(b.size() == a.numCols());
	assert(a.numCols() == a.numRows());
	assert(c.size() == b.size());

	std::vector<MyVector> k(c.size());
	k[0] = k0;

	for (size_t i = 1;i<k.size();++i)
	{
		MyVector y1 = y0;
		for (size_t j = 0;j<i;++j)
		{
			y1 = y1 + len * a[i][j] * k[j];
		}
		k[i] = func(x0 + c[i] * len, y1);
	}

	for (size_t i = 0;i<k.size();++i)
	{
		y0 = y0 + len * b[i] * k[i];
	}

	return y0;
}

/**
 * \brief			Calculates y_n+1 using modified Euler metod.
 * \param[in] func	System of functions for calculating.
 * \param[in] x0	Current x value.
 * \param[in] y0	Current y value.
 * \param[in] len	Length of step.
 * \param[in] a		Matrix with constants
 * \param[in] b		Vector with constants.
 * \param[in] c		Vector with constants.
 * \return			Value of y_n+1.
 */
MyVector RungeKuttStep(const std::function<MyVector(double, const MyVector&)>& func,
	const double x0, const MyVector& y0, const double len, const Matrix&  a, const MyVector& b,
	const MyVector& c)
{
	return RungeKuttStepAgterCalc(func,x0,y0,len,a,b,c,func(x0,y0));
}

std::vector<std::pair<double, MyVector>> RungeKutt(
	const std::function<MyVector(double, const MyVector&)>& func, MyVector y0, double len,
	const double l, const double r, const Matrix&  a, const MyVector& b, const MyVector& c)
{
	assert(l < r);
	const size_t n = (r - l) / len;
	len = (r - l) / n;

	std::vector<std::pair<double, MyVector>> ans;
	ans.reserve(n + 1);

	ans.emplace_back(std::make_pair(l, y0));
	for (size_t i = 0;i<n;++i)
	{
		const double x = l + i * len;

		y0 = RungeKuttStep(func, x, y0, len, a, b, c);

		ans.emplace_back(std::make_pair(x + len, y0));
	}

	return ans;
}

/**
 * \brief			Functoin for powerinf in natural degree.
 * \param[in] a		Nember for powering.
 * \param[in] n		Degree for powering.
 * \return			a^n
 */
double binPow1(const double a, const size_t n)
{
	if (n == 0)
		return 1.;
	double b = binPow1(a,n / 2);
	b = b * b;
	if (n % 2 == 1)
		b *= a;
	return b;
}

/**
 * \brief			Finding first stepfor solving.
 * \param[in] y0	Starting condition in left side of segment.
 * \param[in] l		Left side of segmant.
 * \param[in] r		Rieght side of segment.
 * \param[in] rtol	Relative tolerance for one step.
 * \param[in] atol	Absolute tolerance for one step.
 * \param[in] f		Norm of value of func in (x0,y0)
 * \return			Rude value of h1.
 */
double calcFirstStep(const MyVector& y0, const double l, const double r, const double rtol, 
	const double atol, const double f)
{
	const double d = binPow1(1. / std::max(abs(l), abs(r)), 3) + binPow1(f, 3);
	return pow((rtol*y0.getNorm() + atol) / d, 1. / 3.);
}

/**
 * \brief			Finding first stepfor solving.
 * \param[in] func	System of differential equations.
 * \param[in] y0	Starting condition in left side of segment.
 * \param[in] l		Left side of segmant.
 * \param[in] r		Rieght side of segment.
 * \param[in] rtol	Relative tolerance for one step.
 * \param[in] atol	Absolute tolerance for one step.
 * \return			Acurate value of h1.
 */
double firstStep(const std::function<MyVector(double, const MyVector&)>& func,
	const MyVector& y0, const double l, const double r, const double rtol, const double atol)
{
	MyVector tmpVec = func(l, y0);

	const double h1 = calcFirstStep(y0, l, r, rtol, atol, tmpVec.getNorm());

	const MyVector u1 = y0 + h1 * tmpVec;

	return std::min(h1, calcFirstStep(u1, l + h1, r, rtol, atol,func(l+h1,u1).getNorm()));
}

std::tuple<std::vector<std::pair<double, MyVector>>, std::vector<std::pair<double, double>>,
	std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>, size_t> 
	DifferentialEquationsOppAuto( const std::function<MyVector(double, const MyVector&)>& func,
	MyVector y0, double l, const double r, const double rtol, const double atol, const double mlen)
{
	size_t num = 2;

	std::vector<std::pair<double, MyVector>> ans;

	std::vector<std::pair<double, double>> steps1, steps2, errors;
 
	double h = firstStep(func, y0, l, r, rtol, atol);

	while(l<r)
	{

		ans.emplace_back(std::make_pair(l, y0));
		steps1.emplace_back(std::make_pair(l, h));

		num += 5;

		const MyVector f0 = func(l, y0);

		MyVector y1 = DifferentialEquationsOppStepAfterCalc(func, l, y0, h, f0);
		MyVector y05 = DifferentialEquationsOppStepAfterCalc(func, l, y0, h/2, f0);
		MyVector y2 = DifferentialEquationsOppStep(func, l + h/2, y05, h/2);
		double d = (y1 - y2).getNorm()*4. / 3.;
		
		while(d > (rtol*y0.getNorm() + atol)*4)
		{
			num += 3;

			steps2.emplace_back(std::make_pair(l, h));
			h /= 2.;
			y1 = y05;
			y05 = DifferentialEquationsOppStepAfterCalc(func, l, y0, h / 2,f0);
			y2 = DifferentialEquationsOppStep(func, l + h / 2, y05, h / 2);
			d = (y1 - y2).getNorm()*4. / 3.;
		}

		l += h;
		
		if(d< (rtol*y0.getNorm() + atol)/8)
		{
			steps2.emplace_back(std::make_pair(l, h));
			y0 = y1;
			h = std::min(h * 2, mlen);
		}
		else if (d < (rtol*y0.getNorm() + atol))
		{
			y0 = y1;
		}
		else
		{
			steps2.emplace_back(std::make_pair(l, h));
			y0 = y2;
			h /= 2.;
		}

		errors.emplace_back(std::make_pair(l, d));
	}

	ans.emplace_back(std::make_pair(l, y0));

	return std::make_tuple(ans,steps1,steps2, errors, num);
}

std::tuple<std::vector<std::pair<double, MyVector>>, std::vector<std::pair<double, double>>,
	std::vector<std::pair<double, double>>, std::vector<std::pair<double, double>>, size_t> 
	RungeKuttAuto( const std::function<MyVector(double, const MyVector&)>& func, MyVector y0,
	double l, const double r, const double rtol, const double atol, const double mlen,
	const Matrix& a, const MyVector& b, const MyVector& c)
{
	size_t num = 2;

	std::vector<std::pair<double, MyVector>> ans;

	std::vector<std::pair<double, double>> steps1, steps2, errors;

	double h = firstStep(func, y0, l, r, rtol, atol);

	while (l<r)
	{
		ans.emplace_back(std::make_pair(l, y0));
		steps1.emplace_back(std::make_pair(l, h));

		num += 5;

		const MyVector f0 = func(l, y0);

		MyVector y1 = RungeKuttStepAgterCalc(func, l, y0, h, a, b, c,f0);
		MyVector y05 = RungeKuttStepAgterCalc(func, l, y0, h / 2, a, b, c,f0);
		MyVector y2 = RungeKuttStep(func, l + h / 2, y05, h / 2, a, b, c);
		double d = (y1 - y2).getNorm()*4. / 3.;

		while (d >(rtol*y0.getNorm() + atol) * 4)
		{
			num += 3;

			steps2.emplace_back(std::make_pair(l, h));
			h /= 2.;
			y1 = y05;
			y05 = RungeKuttStepAgterCalc(func, l, y0, h / 2, a, b, c,f0);
			y2 = RungeKuttStep(func, l + h / 2, y05, h / 2, a, b, c);
			d = (y1 - y2).getNorm()*4. / 3.;
		}

		l += h;

		if (d< (rtol*y0.getNorm() + atol) / 8)
		{
			steps2.emplace_back(std::make_pair(l, h));
			y0 = y1;
			h = std::min(h * 2, mlen);
		}
		else if (d < (rtol*y0.getNorm() + atol))
		{
			y0 = y1;
		}
		else
		{
			steps2.emplace_back(std::make_pair(l, h));
			y0 = y2;
			h /= 2.;
		}

		errors.emplace_back(std::make_pair(l, d));
	}

	ans.emplace_back(std::make_pair(l, y0));

	return std::make_tuple(ans, steps1, steps2, errors, num);
}
