#include <iomanip>
#include <fstream>
#include <conio.h>
#include <iostream>

#include "PLUMatrix.h"
#include "QRMatrix.h"
#include "IterationES.h"
#include <cassert>
#include "NewtonMethods.h"

void first(double(*func)(double x), double(*driv)(double x), double eps)
{
	assert(eps > 0);

	const double ans1 = AbsoluteLessRoot1(func, driv, eps);

	std::cout << ans1 << ' ' << func(ans1) << "\n\n";

	const double ans2 = AbsoluteLessRoot2(func, driv, eps);

	std::cout << ans2 << ' ' << func(ans2) << "\n\n";
}

double f(const double x)
{
	return x * x * x - exp(x) + 1;
}

double d(const double x)
{
	return 3 * x * x - exp(x);
}

int main()
{
	std::cout.precision(15);
	std::cout << std::fixed;

	first(f, d, 1e-4);

	_getch();

	return 0;
}
