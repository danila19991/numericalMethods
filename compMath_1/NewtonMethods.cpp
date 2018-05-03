#include "NewtonMethods.h"
#include <cassert>
#include <iostream>
#include "PLUMatrix.h"
#include <tuple>
#include <chrono>

constexpr int operForJakobi = 2267;
constexpr int operForNextVale = 220;
constexpr int calcForJakobi = 2167;
constexpr int calcForNextVale = 210;

double Newton1(double(* func)(double x), double(* driv)(double x), double l, double r, double eps)
{
	assert(func(l)*func(r) <= 0);

	if (l > r)
		std::swap(l, r);

	int steps1 = 0, steps2 = 0;

	while((r-l) > 0.1)
	{
		const double mid = (r + l) / 2;
		if (func(l)*func(mid) > 0)
			l = mid;
		else
			r = mid;

		++steps1;
	}

	double prev = (r + l) / 2;
	double nxt = prev - func(prev) / driv(prev);

	while(abs(prev - nxt) >= eps)
	{
		prev = nxt;
		nxt = prev - func(prev) / driv(prev);

		++steps2;
		if(steps2 > 100000)
		{
			std::cout << "Newton1 is not working";
			return std::numeric_limits<double>::infinity();
		}
	}

	std::cout << "Newton1 work in " << steps1 << " and " << steps2 << '\n';

	return nxt;
}

double AbsoluteLessRoot1(double(* func)(double x), double(* driv)(double x), double eps)
{
	double r=1., l=-1.;

	for (int i = 1; func(eps)*func(i) > 0. && i < 1e6; ++i)
		r = i;

	const double ansr = Newton1(func, driv, eps, r, eps);

	for (int i = -1; func(eps)*func(i) > 0. && i > -1e6; --i)
		l = i;

	const double ansl = Newton1(func, driv, l, -eps, eps);

	return (abs(ansl) < abs(ansr)) ? ansl : ansr;
}

double Newton2(double(* func)(double x), double(* driv)(double x), double l, double r, double eps)
{
	assert(func(l)*func(r) <= 0);
	
	if (l > r)
		std::swap(l, r);
	
	int steps = 0;

	double prev = (r + l) / 2;
	double nxt = nxt = prev - func(prev) / driv(prev);
	while (abs(prev - nxt) >= eps)
	{
		if (func(l)*func(nxt) > 0.)
			l = nxt;
		else
			r = nxt;
		prev = nxt;

		nxt = prev - func(prev) / driv(prev);

		++steps;
		if(steps > 100000)
		{
			std::cout << "Newton2 is not working";
			return std::numeric_limits<double>::infinity();
		}

		if (nxt < l || nxt > r)
			nxt = (l + r) / 2;
	}

	std::cout << "Newton2 work in " << steps << '\n';

	return nxt;
}

double AbsoluteLessRoot2(double(* func)(double x), double(* driv)(double x), double eps)
{
	double r = 1., l = -1.;

	for (int i = 1; func(eps)*func(i) > 0. && i < 1e6; ++i)
		r = i;

	const double ansr = Newton2(func, driv, eps, r, eps);

	for (int i = -1; func(eps)*func(i) > 0. && i > -1e6; --i)
		l = i;

	const double ansl = Newton2(func, driv, l, -eps, eps);

	return (abs(ansl) < abs(ansr)) ? ansl : ansr;
}

std::tuple<MyVector, int, int, int> NewtonES1(MyVector(* func)(const MyVector& x), 
	Matrix(* driv)(const MyVector& x), MyVector fx, const double eps, bool show)
{
	int steps = 1;

	const std::chrono::time_point<std::chrono::steady_clock> start =
		std::chrono::steady_clock::now();

	MyVector nxt = fx - PLUMatrix(driv(fx)).obrat()*func(fx);

	while((fx-nxt).getNorm() > eps)
	{
		++steps;
		
		fx = nxt;
		nxt = fx -PLUMatrix(driv(fx)).obrat()*func(fx);
	}

	const int militime = std::chrono::duration_cast<std::chrono::milliseconds>(
		std::chrono::steady_clock::now() - start).count();

	const int oper = steps * (operForJakobi + operForNextVale);

	const int calc = steps * (calcForJakobi + calcForNextVale);

	if(show)
		std::cout << "Newton in " << steps << " steps\n";

	return std::make_tuple(nxt, militime, oper, calc);
}

std::tuple<MyVector, int, int, int> NewtonES2(MyVector(* func)(const MyVector& x),
	Matrix(* driv)(const MyVector& x), MyVector fx, double eps, bool show)
{
	int steps = 1;

	const std::chrono::time_point<std::chrono::steady_clock> start =
		std::chrono::steady_clock::now();

	//O(n^2 + (2/3*n^3 + n*(n+n^2)) + 4*n^2)=O(6n^2 + 5/3n^3)~operations(2267)
	//O(5n^2 + 5/3n^3)~calculations(2167)
	const Matrix m = PLUMatrix(driv(fx)).obrat();

	//O(n + n + n^2 + 10*n)=O(12n + n^2)~operations(220)
	//O(11n + n^2)~calculations(210)
	MyVector nxt = fx - m*func(fx);

	while ((fx - nxt).getNorm() > eps)//O(2n)
	{
		++steps;

		fx = nxt;//O(n)
		nxt = fx - m*func(fx);
	}

	const int militime = std::chrono::duration_cast<std::chrono::milliseconds>(
		std::chrono::steady_clock::now() - start).count();

	const int oper = steps * operForNextVale + operForJakobi;

	const int calc = steps * calcForNextVale + calcForJakobi;

	if(show)
		std::cout << "Newton in " << steps << " steps\n";

	return std::make_tuple(nxt, militime, oper, calc);
}

std::tuple<MyVector, int, int, int> NewtonES3(MyVector(*func)(const MyVector& x),
	Matrix(*driv)(const MyVector& x), MyVector fx, int k, double eps, bool show)
{
	int steps = 1;

	const std::chrono::time_point<std::chrono::steady_clock> start =
		std::chrono::steady_clock::now();

	Matrix m = PLUMatrix(driv(fx)).obrat();

	MyVector nxt = fx - m * func(fx);

	while ((fx - nxt).getNorm() > eps)
	{
		fx = nxt;
		
		++steps;
		
		if (steps <= k)
			m = PLUMatrix(driv(fx)).obrat();

		nxt = fx - m * func(fx);
	}

	const int militime = std::chrono::duration_cast<std::chrono::milliseconds>(
		std::chrono::steady_clock::now() - start).count();

	int oper;
	if (steps <= k)
		oper = steps * (operForNextVale + operForJakobi);
	else
		oper = steps * operForNextVale + k * operForJakobi;

	int calc;
	if (steps <= k)
		calc = steps * (calcForNextVale + calcForJakobi);
	else
		calc = steps * calcForNextVale + k * calcForJakobi;

	if(show)
		std::cout << "Newton in " << steps << " steps\n";

	return std::make_tuple(nxt, militime, oper, calc);
}

std::tuple<MyVector, int, int, int> NewtonES4(MyVector(* func)(const MyVector& x),
	Matrix(* driv)(const MyVector& x), MyVector fx, int k, double eps, bool show)
{
	assert(k > 0);

	int steps = 1;

	const std::chrono::time_point<std::chrono::steady_clock> start =
		std::chrono::steady_clock::now();

	Matrix m = PLUMatrix(driv(fx)).obrat();

	MyVector nxt = fx - m * func(fx);

	while ((fx - nxt).getNorm() > eps)
	{
		++steps;

		if((steps-1)%k == 0)
			m = PLUMatrix(driv(fx)).obrat();

		fx = nxt;
		nxt = fx - m * func(fx);
	}

	const int militime = std::chrono::duration_cast<std::chrono::milliseconds>(
		std::chrono::steady_clock::now() - start).count();

	const int oper = steps * operForNextVale + ((steps - 1) / k + 1)*operForJakobi;

	const int calc = steps * calcForNextVale + ((steps - 1) / k + 1)*calcForJakobi;

	if(show)
		std::cout << "Newton in " << steps << " steps\n";

	return std::make_tuple(nxt, militime, oper, calc);
}

