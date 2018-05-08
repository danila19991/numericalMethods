#include "NewtonMethods.h"
#include <cassert>
#include <iostream>
#include "PLUMatrix.h"
#include <tuple>
#include <chrono>

constexpr int operForJakobi = 1167;
constexpr int operForNextVale = 220;
constexpr int calcForJakobi = 1067;
constexpr int calcForNextVale = 210;

double Newton(double(* func)(double x), double(* driv)(double x), double l, double r, double eps)
{
	assert(func(l)*func(r) <= 0.);
	
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

double AbsoluteLessRoot(double(*func)(double x), double(*driv)(double x), double eps)
{
	double r = 1., l = -1.;

	int i = 1;

	while (func(eps)*func(i) > 0. && i < 1'000'000)
	{
		++i;
		r = i;
	}

	const double ansr = Newton(func, driv, eps, r, eps);

	i = -1;

	while (func(eps)*func(i) > 0. && i > -1'000'000)
	{
		--i;
		l = i;
	}

	const double ansl = Newton(func, driv, l, -eps, eps);

	return (abs(ansl) < abs(ansr)) ? ansl : ansr;
}

std::tuple<MyVector, int, int, int> NewtonES1(MyVector(* func)(const MyVector& x), 
	Matrix(* driv)(const MyVector& x), MyVector fx, const double eps, bool show)
{
	int steps = 1;

	const std::chrono::time_point<std::chrono::steady_clock> start =
		std::chrono::steady_clock::now();

	MyVector nxt = fx - PLUMatrix(driv(fx)).solve(func(fx));

	while((fx-nxt).getNorm() + func(nxt).getNorm() > eps)
	{
		++steps;

		fx = nxt;
		nxt = fx - PLUMatrix(driv(fx)).solve(func(fx));
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

	//O(n^2 + (2/3*n^3 + 4*n^2)~operations(1167)
	//O(4n^2 + 2/3n^3)~calculations(1067)
	const PLUMatrix m = PLUMatrix(driv(fx));

	//O(n + n + n^2 + 10*n)=O(12n + n^2)~operations(220)
	//O(11n + n^2)~calculations(210)
	MyVector nxt = fx - m.solve(func(fx));

	while ((fx - nxt).getNorm() + func(nxt).getNorm() > eps)//O(2n)
	{
		++steps;

		fx = nxt;//O(n)
		nxt = fx - m.solve(func(fx));
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

	PLUMatrix m = PLUMatrix(driv(fx));

	MyVector nxt = fx - m.solve(func(fx));

	while ((fx - nxt).getNorm() + func(nxt).getNorm() > eps)
	{
		fx = nxt;
		
		++steps;
		
		if (steps <= k)
			m = PLUMatrix(driv(fx));

		nxt = fx - m.solve(func(fx));
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

	PLUMatrix m = PLUMatrix(driv(fx));

	MyVector nxt = fx - m.solve(func(fx));

	while ((fx - nxt).getNorm() + func(nxt).getNorm() > eps)
	{
		++steps;

		if((steps-1)%k == 0)
			m = PLUMatrix(driv(fx));

		fx = nxt;
		nxt = fx - m.solve(func(fx));
	}

	const int militime = std::chrono::duration_cast<std::chrono::milliseconds>(
		std::chrono::steady_clock::now() - start).count();

	const int oper = steps * operForNextVale + ((steps - 1) / k + 1)*operForJakobi;

	const int calc = steps * calcForNextVale + ((steps - 1) / k + 1)*calcForJakobi;

	if(show)
		std::cout << "Newton in " << steps << " steps\n";

	return std::make_tuple(nxt, militime, oper, calc);
}

std::tuple<MyVector, int, int, int> NewtonESh(MyVector(* func)(const MyVector& x),
	Matrix(* driv)(const MyVector& x), MyVector fx, double eps, bool show)
{
	int steps = 1, steps2 = 1;

	const std::chrono::time_point<std::chrono::steady_clock> start =
		std::chrono::steady_clock::now();

	PLUMatrix m = PLUMatrix(driv(fx));

	MyVector nxt = fx - m.solve(func(fx));

	int prev = (log((fx - nxt).getNorm())-log(eps))/log(10.);

	while ((fx - nxt).getNorm() + func(nxt).getNorm() > eps)
	{
		++steps;

		if((log((fx - nxt).getNorm()) - log(eps))/log(10.) + 2 < prev)
		{
			prev = log((log((fx - nxt).getNorm()) -log(eps))/log(10.));
			m = PLUMatrix(driv(fx));

			++steps2;
		}

		fx = nxt;
		nxt = fx - m.solve(func(fx));
	}

	const int militime = std::chrono::duration_cast<std::chrono::milliseconds>(
		std::chrono::steady_clock::now() - start).count();

	const int oper = steps2 * operForJakobi + steps * operForNextVale;

	const int calc = steps2 * calcForJakobi + steps * calcForNextVale;

	if (show)
		std::cout << "Newton in " << steps << " steps and " << steps2 << "\n";

	return std::make_tuple(nxt, militime, oper, calc);
}

