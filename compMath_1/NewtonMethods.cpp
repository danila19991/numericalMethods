#include "NewtonMethods.h"
#include <cassert>
#include <iostream>

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
	double r, l;

	for (r = 1.; func(eps)*func(r) > 0. && r < 1e6; r += 1.);

	const double ansr = Newton1(func, driv, eps, r, eps);

	for (l = -1.; func(-eps)*func(l) > 0. && l > -1e6; r -= 1.);

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
	double r, l;

	for (r = 1.; func(eps)*func(r) > 0. && r < 1e6; r += 1.);

	const double ansr = Newton2(func, driv, eps, r, eps);

	for (l = -1.; func(-eps)*func(l) > 0. && l > -1e6; r -= 1.);

	const double ansl = Newton2(func, driv, l, -eps, eps);

	return (abs(ansl) < abs(ansr)) ? ansl : ansr;
}


